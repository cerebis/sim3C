import logging
import numpy as np
from numba import njit

from collections import OrderedDict

# from .empirical_model import EmpiricalDistribution, cdf_geom_unif_ratio, generate_nested_cids, cids_to_blocks
from .empirical_model import EmpiricalDistribution, cdf_geom_unif_ratio
from .exceptions import *
from .site_analysis import AllSites, CutSites
import sim3C.random as random

logger = logging.getLogger(__name__)


@njit("i4(f8[:], f8)")
def random_index(cdf, rv):
    """
    A weighted number draw, functioning like numpy.random.choice. This
    approach is faster than choice, which also cannot be Numba-compiled
    when supplying weights. Additionally, with JIT compilation this
    is a further 3X faster than without. Note: we cannot make the call to
    pcg_uniform() from here as Numba cannot infer its return type.
    """
    return np.searchsorted(cdf, rv)


def cdf_choice(vals, cdf):
    """
    Random selection of an element from an array, biased by the supplied CDF.
    :param vals: the array of potential choices
    :param cdf: the CDF describing each elements probability of selection
    :return: the selected element
    """
    # TODO explore whether moving this method to faster.pyx would be faster
    return vals[random_index(cdf, pcg_uniform())]


class Segment(object):
    """
    A segment of DNA sequence which may represent all or only part of a complete
    replicon (chromosome, plasmid, etc).
    """

    def __init__(self, name, repl, seq, enzyme, anti_rate, create_cids=True, linear=False):
        """
        The definition of a segment of DNA.
        :param name: a unique name for this segment
        :param repl: the parent replicon for this segment
        :param seq: the sequence of this segment
        :param enzyme: the enzyme used to digest DNA in the 3C/HiC library preparation
        :param anti_rate: the rate of anti-diagonal interactions
        :param create_cids: when true, simulate chromosome-interacting-domains
        :param linear: treat replicons as linear
        """
        self.name = name
        self.seq = seq
        self.linear = linear

        self.parent_repl = repl
        repl.register_segment(self)

        # cut-site related properties. These are pre-calculated as a simple
        # means of avoiding performance penalties with repeated calls.
        if not enzyme:
            self.sites = AllSites(len(seq.seq))
        else:
            self.sites = CutSites(enzyme, seq.seq, linear=linear)

        self.length = len(self.seq)
        self.num_sites = self.sites.size
        self.site_density = self.num_sites / self.length

        if create_cids:
            logger.warning('CID model currently disabled')
        self.draw_constrained_site = self._draw_simple_constrained_site

        if self.linear:
            if anti_rate > 0:
                logger.warning('Replicon {} is linear, anti_rate of {} has been set to zero'
                               .format(name, anti_rate))
            self.anti_rate = 0
        else:
            self.anti_rate = anti_rate

        # setup for simple model
        self.empdist = EmpiricalDistribution(self.length,
                                             Replicon.GLOBAL_EMPDIST_BINS, cdf_geom_unif_ratio,
                                             cdf_alpha=Replicon.CDF_ALPHA,
                                             shape=Replicon.GLOBAL_SHAPE_FACTOR)

    def __repr__(self):
        return '{} {}'.format(self.name, self.parent_repl)

    def subseq(self, pos1, length, rev=False):
        """
        Create a subsequence.

        :param pos1: starting genomic position
        :param length: length of subsequence
        :param rev: reverse complement this sequence.
        :return: subseq Seq object
        """

        # handle negative starts as wrapping around or,
        # in linear cases, terminates at first position
        if pos1 < 0:
            pos1 = pos1 % self.length if not self.linear else 0

        pos2 = pos1 + length
        diff = pos2 - self.length
        if diff > 0:
            if self.linear:
                # sequence terminates early
                ss = self.seq[pos1:-1]
                ss.description = Replicon.PART_DESC_FMT.format(rev, pos1 + 1, self.length)
            else:
                # sequence will wrap around
                ss = self.seq[pos1:] + self.seq[:diff]
                ss.description = Replicon.PART_DESC_FMT.format(rev, pos1 + 1, diff)
        else:
            ss = self.seq[pos1:pos2]
            ss.description = Replicon.PART_DESC_FMT.format(rev, pos1 + 1, pos2)

        if rev:
            ss.reverse_complement(id=True, description=True)

        return ss

    def covers_site(self, pos, length):
        """
        Test if the segment defined by a position and length cover any cut-sites
        for this replicon and enzyme digestion.
        :param pos: the beginning position
        :param length: the segment length
        :return: True - this segment covers at least one site.
        """
        return self.sites.covers(pos, length)

    def draw_any_site(self):
        """
        Uniform selection of a random site
        :return: a cut-site
        """
        return self.sites.random_site()

    def draw_any_location(self):
        """
        Uniform selection of any genomic coordinate for this replicon
        :return: a genomic coord (zero based)
        """
        return pcg_integer(self.length)

    @staticmethod
    def get_loc_3c(emp_dist, pos1, length):
        """
        Get a second genomic location (pos2) on this replicon, where the separation |pos2-pos1|
        is constrained by the experimentally determined distribution for 3C/HiC ligation
        products.

        TODO: this method does not explicitly treat linear replicons. To do so, we must
           handle the edge case of drawing beyond first and last positions. Simple approaches
           using redraw are potentially computationally terrible if the first position is
           very close to the end of the replicon. For now, we accept the modulo solution
           for both linear and circular.

        :param emp_dist: empirical distribution of separation
        :param pos1: the first position
        :param length: the length (bp) of the replicon
        :return: pos2: the second position
        """
        # draw a random separation
        delta = int(emp_dist.rand())

        # pve or nve shift, modulo length
        if pcg_integer(2) == 0:
            pos2 = (pos1 - delta) % length
        else:
            pos2 = (pos1 + delta) % length
        return pos2

    def _draw_simple_constrained_site(self, pos1):
        """
        Draw a second site (pos2) relative to the first (x1) which is constrained
        to follow a basic empirical distribution.
        :param pos1: the first location
        :return: pos2: the second location
        """
        pos2 = Segment.get_loc_3c(self.empdist, pos1, self.length)

        # anti-diagonal
        if pcg_uniform() < self.anti_rate:
            pos2 = self.length - pos2

        # return nearest site
        return self.sites.find_nn(pos2)


class Replicon(object):
    """
    A DNA molecule which will be digested. This may be a chromosome, plasmid, etc.
    """

    # Formatting string for sequence descriptions of a part (subsequence) of a replicon.
    # Used in insert/read creation.
    PART_DESC_FMT = '{:d}:{}:{}'

    # Empirical distribution parameters. These might eventually be exposed to users.
    BACKBONE_PROB = 0.25
    GLOBAL_EMPDIST_BINS = 1000
    GLOBAL_SHAPE_FACTOR = 3.5e-6
    CDF_ALPHA = 0.1
    CID_EMPDIST_BINS = 100
    CID_SHAPE_FACTOR = 8.0e-6
    CID_MIN = 3
    CID_MAX = 6
    CID_DEPTH = 2

    def __init__(self, name, cell, cn):
        """
        The definition of a replicon (chromosome, plasmid, etc).
        :param name: a unique name for this replicon
        :param cell: the parent cell for this replicon
        :param cn: the copy-number of this replicon in 'cell'
        """

        self.name = name
        self.copy_number = cn
        # segment registry
        self.segment_registry = OrderedDict()
        self.cdf_extent = None
        self.cdf_sites = None
        self.segment_names = None
        # TODO reimplement CID code for Segments

        # if create_cids:
        #     # setup for more complex simulated CID model
        #     self.draw_constrained_site = self._draw_cid_constrained_site
        #     self.cid_blocks = cids_to_blocks(
        #         generate_nested_cids(self.length, Replicon.BACKBONE_PROB,
        #                              Replicon.GLOBAL_EMPDIST_BINS, Replicon.GLOBAL_SHAPE_FACTOR,
        #                              Replicon.CID_EMPDIST_BINS, Replicon.CID_SHAPE_FACTOR,
        #                              cdf_alpha=Replicon.CDF_ALPHA,
        #                              min_num=Replicon.CID_MIN, max_num=Replicon.CID_MAX,
        #                              recur_depth=Replicon.CID_DEPTH))
        #
        # else:
        #     # setup for simple model
        #     self.draw_constrained_site = self._draw_simple_constrained_site
        #     self.empdist = EmpiricalDistribution(self.length,
        #                                          Replicon.GLOBAL_EMPDIST_BINS, cdf_geom_unif_ratio,
        #                                          cdf_alpha=Replicon.CDF_ALPHA,
        #                                          shape=Replicon.GLOBAL_SHAPE_FACTOR)

        # set bidirection association with containing cell
        self.parent_cell = cell
        cell.register_replicon(self)

    def __repr__(self):
        return '{} {}'.format(self.name, self.parent_cell)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.name == other.name

    def init_prob(self):
        seg_info = self.get_segment_info()
        pdf_extent = np.array(seg_info['length'], dtype='f8')
        pdf_sites = np.array(seg_info['sites'], dtype='f8')
        pdf_extent /= pdf_extent.sum()
        pdf_sites /= pdf_sites.sum()
        self.cdf_extent = np.cumsum(pdf_extent)
        self.cdf_sites = np.cumsum(pdf_sites)
        self.segment_names = np.array(seg_info['name'])

    def get_length(self):
        """
        :return: the length of this replicon
        """
        return sum([seg.length for seg in self.segment_registry.values()])

    def get_num_sites(self):
        """
        :return: the number of cut-sites for this replicon
        """
        return sum([seg.num_sites for seg in self.segment_registry.values()])

    def _max_name_len(self):
        return max([len(seg.name) for seg in self.segment_registry.values()])

    def get_segment_info(self):
        """
        :return: a list of segment site counts
        """
        return np.array([(seg.name, seg.length, seg.num_sites) for seg in self.segment_registry.values()],
                        dtype=np.dtype([('name', f'U{self._max_name_len()}'), ('length', 'i4'), ('sites', 'i4')]))

    def register_segment(self, segment):
        """
        Insert a segment into this replicon's registry
        :param segment: the segment to insert
        :return: the inserted segment
        """
        assert isinstance(segment, Segment),  'Attempted to register invalid class.'
        if segment.name in self.segment_registry:
            raise Sim3CException('duplicate segment names')
        self.segment_registry[segment.name] = segment
        return segment

    def draw_any_segment_by_sites(self):
        """
        Draw any segment from the replicon, based on each segments number of sites
        :return: a Segment
        """
        return self.segment_registry[cdf_choice(self.segment_names, self.cdf_sites)]

    def draw_any_segment_by_extent(self):
        """
        Draw any segment from the replicon, biased by the segment length.
        :return: a Segment
        """
        return self.segment_registry[cdf_choice(self.segment_names, self.cdf_extent)]

    # def _draw_cid_constrained_site(self, x1):
    #     """
    #     Draw a second site (x2) relative to teh first (x1) which is constrained to follow
    #     a nested hierarchy of empirical distributions, intended to simulate chromosomal
    #     interaction domains (CID). Here, CIDs are small regions of a chromosome which interact
    #     (possibly through folding) at a more frequent level than that of the overall backbone.
    #
    #     The range of effect of the CIDs are predetermined at instantiation time and stored as
    #     intervals within an interval-tree. When called, x1 determines which overlapping intervals
    #     are involved and single emp-dist (representing a particular CID) is chosen at random
    #     (randomly weighted).
    #
    #     :param x1: the first location
    #     :return: x2: the second location.
    #     """
    #     block = self.cid_blocks[x1].pop()
    #     ovl_invs = block.data['inv_list']
    #
    #     if len(ovl_invs) == 1:
    #         # only the background distribution governs this block
    #         chosen_inv = ovl_invs[0]
    #
    #         x2 = Replicon.get_loc_3c(chosen_inv.data['empdist'], x1, chosen_inv.length())
    #
    #         # anti-diagonal
    #         if pcg_uniform() < self.anti_rate:
    #             x2 = self.length - x2
    #
    #     else:
    #         # pick a cid or background from those defined for this block
    #         # note, numpy will not accept the list of intervals here
    #         ix = np_choice(len(ovl_invs), p=block.data['prob_list'])
    #         chosen_inv = ovl_invs[ix]
    #
    #         x2 = Replicon.get_loc_3c(chosen_inv.data['empdist'], x1 - chosen_inv.begin, chosen_inv.length())
    #         x2 += chosen_inv.begin
    #
    #         # main backbone gets anti-diagonal treatment
    #         if chosen_inv.data['depth'] == 0:
    #             # antidiagonal
    #             if pcg_uniform() < self.anti_rate:
    #                 x2 = self.length - x2
    #
    #     return self.sites.find_nn(x2)

    def num_segments(self):
        """
        :return: the number of segments for this replicon.
        """
        return len(self.segment_registry)


class Cell(object):
    """
    A cell acts as the container of one or more replicons, where each may have its
    own copy-number in addition to the relative abundance of their containing cell.
    """

    def __init__(self, name, abundance, trans_rate=0.1):
        """
        The definition of a cell.
        :param name: a unique name for this cell.
        :param abundance: the relative abundance of this cell in the community.
        :param trans_rate: the rate of inter-replicon (trans) ligation products
        """

        # if no object supplied, initialise one.
        self.name = name
        self.abundance = abundance
        # replicons are kept in the order they are registered
        self.replicon_registry = OrderedDict()
        self.segment_registry = OrderedDict()
        # inter-replicon (trans) rate
        self.trans_rate = trans_rate

        # probability properties, to be initialised by init_prob()
        self.cdf_extent = None
        self.cdf_sites = None
        # for this cell, indices of segments within or outside replicons
        # indices used with pdf and cdf arrays
        self.cis_segment_indices = {}
        self.trans_segment_indices = {}
        self.replicon_names = None
        self.segment_names = None
        self.cdf_sites_inter = None
        self.cdf_extents_inter = None

    def __repr__(self):
        return repr((self.name, self.abundance, self.num_replicons()))

    def num_replicons(self):
        """
        :return: the number of replicons for this cell.
        """
        return len(self.replicon_registry)

    def num_segments(self):
        """
        :return: the number of segments for this cell.
        """
        return sum([rep.num_segments() for rep in self.replicon_registry.values()])

    def init_prob(self):
        """
        Initialise the selection probabilities. This method should be called after all replicons are
        registered or when a new replicon is added.
        """
        if self.num_replicons() <= 0:
            raise NoRepliconsException(self.name)

        # begin with some empty PDFs
        pdf_extent = np.zeros(self.num_segments())
        pdf_sites = np.zeros(self.num_segments())

        # for each replicon, the PDF for the various modes of selection.
        i = 0
        self.cis_segment_indices = OrderedDict()
        for repl in self.replicon_registry.values():
            repl.init_prob()
            cn = float(repl.copy_number)
            seg_info = repl.get_segment_info()
            self.cis_segment_indices[repl.name] = []
            for si in seg_info:
                pdf_extent[i] = cn * si['length']
                pdf_sites[i] = cn * si['sites']
                self.cis_segment_indices[repl.name].append(i)
                i += 1

        self.trans_segment_indices = OrderedDict()
        for ri in self.cis_segment_indices:
            self.trans_segment_indices[ri] = []
            for rj, seg_idx in self.cis_segment_indices.items():
                if ri != rj:
                    self.trans_segment_indices[ri].extend(seg_idx)

        # normalise the PDFs
        pdf_extent /= pdf_extent.sum()
        pdf_sites /= pdf_sites.sum()

        # CDFs from the PDFs. Selection is accomplished by drawing a
        # unif([0..1]) and mapping that through the relevant CDF.
        self.cdf_extent = np.cumsum(pdf_extent)
        self.cdf_sites = np.cumsum(pdf_sites)

        # array of names for extracting random sequences
        self.replicon_names = np.array([*self.replicon_registry])

        # prepare a registry of all defined segments for this replicon
        for repl in self.replicon_registry.values():
            for seg in repl.segment_registry.values():
                self.segment_registry[seg.name] = seg
        self.segment_names = np.array([*self.segment_registry])

        # One last set of CDFs for "select other" which exclude
        # each replicon in turn.
        if self.num_replicons() > 1:

            self.cdf_sites_inter = {}
            self.cdf_extents_inter = {}

            for rn in self.replicon_names:

                # indices without rn
                xi = self.trans_segment_indices[rn]

                # site probs without rn
                pi = pdf_sites[xi]
                pi /= pi.sum()
                self.cdf_sites_inter[rn] = {'names': self.segment_names[xi], 'prob': np.cumsum(pi)}

                # extent probs without rn
                pi = pdf_extent[xi]
                pi /= pi.sum()
                self.cdf_extents_inter[rn] = {'names': self.segment_names[xi], 'prob': np.cumsum(pi)}

    def register_replicon(self, repl):
        """
        Insert a replicon into this cells registry
        :param repl: the replicon to insert
        :return: the inserted replicon
        """
        assert isinstance(repl, Replicon),  'Attempted to register invalid class.'
        if repl.name not in self.replicon_registry:
            self.replicon_registry[repl.name] = repl
        return repl

    def draw_any_segment_by_extent(self):
        """
        Draw any segment, replicon tuple from this cell. The probability is biased by per-replicon
        genomic extent and copy number.
        :return: a segment,replicon tuple from this cell
        """
        return self.segment_registry[cdf_choice(self.segment_names, self.cdf_extent)]

    def draw_any_segment_by_sites(self):
        """
        Draw any segment, replicon tuple from this cell. The probability is biased by per-replicon
        number of sites and copy number.
        :return: a segment, replicon tuple from this cell
        """
        return self.segment_registry[cdf_choice(self.segment_names, self.cdf_sites)]

    def draw_other_segment_by_sites(self, skip_repl):
        """
        Draw a segment from a different replicon in this cell. The probability is biased by per-replicon
        number of sites and copy number. Note: single replicon cell definitions will
        raise an exception.
        :param skip_repl: the replicon to exclude
        :return: segment of a different replicon in this cell
        """
        if self.num_replicons() <= 1:
            raise MonochromosomalException('inter-replicon events are not possible for monochromosomal cells')

        return self.segment_registry[cdf_choice(self.cdf_sites_inter[skip_repl]['names'],
                                                self.cdf_sites_inter[skip_repl]['prob'])]

    def draw_other_segment_by_extent(self, skip_repl):
        """
        Draw a segment from adifferent replicon in this cell. The probability is biased by per-replicon
        extent and copy number. Note: single replicon cell definitions will
        raise an exception.
        :param skip_repl: the replicon to exclude
        :return: another replicon from this cell
        """
        if self.num_replicons() <= 1:
            raise MonochromosomalException('inter-replicon events are not possible for monochromosomal cells')

        return self.segment_registry[cdf_choice(self.cdf_extents_inter[skip_repl]['names'],
                                                self.cdf_extents_inter[skip_repl]['prob'])]

    def draw_any_site(self):
        """
        Draw a cut-site from any replicon within this cell. The probability
        of drawing a replicon is biased by the per-replicon number of sites
        and copy number, while genomic location is uniform.
        :return: a tuple of (segment, location)
        """
        segment = self.draw_any_segment_by_sites()
        return segment, segment.draw_any_site()

    def print_report(self):
        """
        Print a simple report about this cell.
        """
        print(f'names {self.replicon_names}')
        print(f'cdf_extent {self.cdf_extent}')
        print(f'cdf_site {self.cdf_sites}')

    def is_trans(self):
        """
        Coin toss test for whether an inter-replicon (trans) ligation product was formed. This
        is dictated by the rate supplied at instantiation time (trans_rate).
        :return: True -- treat this as a trans event
        """
        return self.num_replicons() > 1 and pcg_uniform() < self.trans_rate


class Community(object):
    """
    A community represents the entire collection and organisation of DNA molecules (replicons) in a simulation.
    This may be the approximation of an environmental sample, a multi-chromosomal or even monochromosomal
    organism.

    It is the entry point for the supporting reference data in a simulation, including such things as
    the relative abundance profile, DNA sequences and selected restriction enzyme. Additionally, are number
    of simulation parameters are exposed.
    """

    def __init__(self, seq_index, profile, enzyme, anti_rate=0.2, spurious_rate=0.02,
                 trans_rate=0.1, create_cids=True, linear=False):
        """
        Initialise a community.

        :param seq_index: an open sequence index of type _IndexedSeqFileDict as returned from Bio.SeqIO.index
        :param profile: the accompanying abundance profile of all replicons in the community
        :param enzyme: the enzyme used to digest DNA in the 3C/HiC library preparation
        :param anti_rate: the rate of anti-diagonal interactions
        :param spurious_rate: the rate of spurious ligation products
        :param trans_rate: the rate of inter-replicon (trans) ligation products within a cell
        :param create_cids: when true, simulate chromosome-interacting-domains
        :param linear: treat replicons as linear
        """
        global pcg_integer, pcg_integer, pcg_integers, pcg_uniform, pcg_nucleotide, pcg_knockout, pcg_parse_error

        pcg_uniform = random.pcg_random.uniform
        pcg_integer = random.pcg_random.integer

        # global segment, replicon, and cell registries
        self.segm_registry = OrderedDict()
        self.repl_registry = OrderedDict()
        self.cell_registry = OrderedDict()

        # initialise the registries using the community profile
        for row_i in profile.values():
            # register the cell
            cell = self._register_cell(Cell(row_i.cell, row_i.abundance, trans_rate))
            repl = self._register_replicon(Replicon(row_i.molecule, cell, row_i.copy_number))
            try:
                # fetch the sequence from file
                rseq = seq_index[row_i.name].upper()
                # community-wide replicon registry
                self._register_segment(Segment(row_i.name, repl, rseq, enzyme, anti_rate, create_cids, linear))
            except NoCutSitesException:
                logger.warning('Sequence "{}" had no cut-sites for the enzyme {} and will be ignored'.format(
                    row_i.name, str(enzyme)))

        # now we're finished reading the profile, initialise the probabilities for each cell
        for cell in self.cell_registry.values():
            try:
                cell.init_prob()
            except NoRepliconsException:
                logger.warning('Cell "{}" had no usable replicons and will be ignored'.format(cell.name))
                del self.cell_registry[cell.name]

        # now initialise the probs for the whole community
        pdf_extent = np.zeros(len(self.segm_registry))
        pdf_sites = np.zeros(len(self.segm_registry))

        # whether site (prox-lig) or extent (wgs) based, probs are
        # weighted by cellular abundance and copy number
        i = 0
        for repl in self.repl_registry.values():
            cn = float(repl.copy_number)
            seg_info = repl.get_segment_info()
            for si in seg_info:
                pdf_extent[i] = cn * si['length']
                pdf_sites[i] = cn * si['sites']
                i += 1

        # now normalise pdfs
        pdf_extent /= pdf_extent.sum()
        pdf_sites /= pdf_sites.sum()

        # derive cdfs from pdfs, these are what's used for drawing values
        self.cdf_extent = np.cumsum(pdf_extent)
        self.cdf_sites = np.cumsum(pdf_sites)

        # keep a list of names in numpy format
        self.repl_names = np.array([*self.repl_registry])
        self.segm_names = np.array([*self.segm_registry])

        # inter-cellular rate is scaled by the product of all chrom site probs
        self.spurious_rate = spurious_rate

        # keep the number of cells handy
        self.num_cells = len(self.cell_registry)

        # we must have at least one cell defined
        if self.num_cells <= 0:
            logger.error('Community did not contain any useful cell definitions. '
                         'Check that your reference sequence(s) were valid and '
                         'that they each contain at least one cut-site.')
            raise Sim3CException('Community contains no cells')

    def _register_cell(self, cell):
        """
        Add an instance of Cell to the cell registry
        :param cell: a Cell object to add
        :return: the added cell instance
        """
        assert isinstance(cell, Cell), 'Attempted to register invalid class.'
        if cell.name not in self.cell_registry:
            self.cell_registry[cell.name] = cell
        return self.cell_registry[cell.name]

    def _register_replicon(self, repl):
        """
        Add an instance of Replicon to the replicon registry
        :param repl: a Replicon object to add
        :return: the added replicon instance
        """
        assert isinstance(repl, Replicon),  'Attempted to register invalid class.'
        if repl.name not in self.repl_registry:
            self.repl_registry[repl.name] = repl
        return self.repl_registry[repl.name]

    def _register_segment(self, segment):
        """
        Add an instance of Segment to the segment registry
        :param segment:
        :return: the added segment instance
        """
        assert isinstance(segment, Segment),  'Attempted to register invalid class.'
        if segment.name in self.segm_registry:
            raise Sim3CException('duplicate segment names in community')
        self.segm_registry[segment.name] = segment
        return segment

    def get_segment(self, name):
        """
        Return a segment from the registry
        :param name: the name of the segment
        :return: the Segment instance
        """
        return self.segm_registry[name]

    def get_repl(self, name):
        """
        Return a replicon from the registry
        :param name: the name of the replicon
        :return: the Replicon instance
        """
        return self.repl_registry[name]

    def get_cell(self, name):
        """
        Return a cell from the registry
        :param name: the name of the cell
        :return: the Cell instance
        """
        return self.cell_registry[name]

    def draw_any_segment_by_extent(self):
        """
        Draw any replicon from this community. The probability is biased by cellular abundance,
        and per-replicon genomic extent and copy number.
        :return: any replicon from this community
        """
        return self.segm_registry[cdf_choice(self.segm_names, self.cdf_extent)]

    def draw_any_segment_by_sites(self):
        """
        Draw any replicon from this community. The probability is biased by cellular abundance,
        and per-replicon number of cut-sites and copy number.
        :return: any replicon from this community
        """
        return self.segm_registry[cdf_choice(self.segm_names, self.cdf_sites)]

    def draw_any_by_site(self):
        """
        Draw any site from any replicon, biased by abundance, number of sites and copy number.
        :return:
        """
        segment = self.draw_any_segment_by_sites()
        return segment, segment.draw_any_site()

    def draw_any_by_extent(self):
        """
        Draw any site from any replicon, biased by abundance, number of sites and copy number.
        :return:
        """
        segment = self.draw_any_segment_by_extent()
        return segment, segment.draw_any_location()

    def print_report(self):
        print(f'names {self.repl_names}')
        print(f'cdf_ext {self.cdf_extent}')
        print(f'cdf_sit {self.cdf_sites}')

    def is_spurious(self):
        """
        Coin toss test for whether a spurious ligation product was formed. This is
        dictated by the rate supplied at instantiation time (spurious_rate).
        :return: True -- treat this as a spurious event
        """
        return pcg_uniform() < self.spurious_rate
