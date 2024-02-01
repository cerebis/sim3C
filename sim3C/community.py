import logging
import numpy as np
from numba import njit

from collections import OrderedDict

from .abundance import read_toml, read_profile
from .empirical_model import EmpiricalDistribution, cdf_geom_unif_ratio
from .exceptions import *
from .site_analysis import AllSites, CutSites
import sim3C.random as random

logger = logging.getLogger(__name__)

# as the base RNG of PCG returns 32 bit integers, the uniform
# function only possesses 32bit float precision. This _may_
# cause unintended problems if mixed in expressions with
# 64bit floats.
PROB_TYPE = 'f4'

# this constant is used pedantically, to ensure that
# a CDF array never ends in a value less than 1.0
SLIGHTLY_LARGER_THAN_ONE = 1.0 + 1.0e-9

@njit(f'i4({PROB_TYPE}[:], {PROB_TYPE})')
def random_index(cdf, rv):
    """
    A weighted number draw, functioning like numpy.random.choice. This
    approach is faster than choice, which also cannot be Numba-compiled
    when supplying weights. Additionally, with JIT compilation this
    is a further 3X faster than without. Note: we cannot make the call to
    pcg_uniform() from here as Numba cannot infer its return type.
    """
    return np.searchsorted(cdf, rv)


@njit(f'{PROB_TYPE}[:]({PROB_TYPE}[:])')
def pdf2cdf(pdf):
    """
    Convert a PDF to a CDF, making sure that the last element is 1. The PDF need
    not be normalised, as we only require that the CDF's last element 1.
    :return: the CDF array
    """
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    cdf[-1] = SLIGHTLY_LARGER_THAN_ONE
    return cdf


def cdf_choice(vals, cdf):
    """
    Random selection of an element from an array, biased by the supplied CDF.
    :param vals: the array of potential choices
    :param cdf: the CDF describing each elements probability of selection
    :return: the selected element
    """
    # TODO explore whether moving this method to faster.pyx would be faster
    return vals[random_index(cdf, random.pcg_random.uniform())]


class Segment(object):
    """
    A segment of DNA sequence which may represent all or only part of a complete
    replicon (chromosome, plasmid, etc).
    """

    def __init__(self, name, repl, seq, enzyme, create_cids=False):
        """
        The definition of a segment of DNA.
        :param name: a unique name for this segment
        :param repl: the parent replicon for this segment
        :param seq: the sequence of this segment
        :param enzyme: the enzyme used to digest DNA in the 3C/HiC library preparation
        :param create_cids: when true, simulate chromosome-interacting-domains
        """
        self.name = name
        self.seq = seq
        self.length = len(self.seq)
        self.parent_repl = None

        # these are inherited from the containing replicon
        self.linear = repl.linear
        self.anti_rate = repl.anti_rate
        if self.linear:
            if self.anti_rate > 0:
                logger.warning(f'Replicon {name} is linear, anti_rate of {self.anti_rate} has been set to zero')
            self.anti_rate = 0

        # cut-site related properties. These are pre-calculated as a simple
        # means of avoiding performance penalties with repeated calls.
        if enzyme is None:
            self.sites = AllSites(len(seq.seq))
        else:
            try:
                self.sites = CutSites(enzyme, seq.seq, linear=self.linear)
            except NoCutSitesException as e:
                self.num_sites = 0
                self.site_density = 0
                raise e

        self.num_sites = self.sites.size
        self.site_density = self.num_sites / self.length

        if not isinstance(repl, Replicon):
            raise TypeError('repl attribute must be of type Replicon')
        self.parent_repl = repl
        repl.register_segment(self)

        if create_cids:
            logger.warning('CID model currently disabled')
        self.draw_constrained_site = self._draw_simple_constrained_site

        # setup for simple model
        self.empdist = EmpiricalDistribution(self.length,
                                             Replicon.GLOBAL_EMPDIST_BINS, cdf_geom_unif_ratio,
                                             cdf_alpha=Replicon.CDF_ALPHA,
                                             shape=Replicon.GLOBAL_SHAPE_FACTOR)

    def __repr__(self):
        return '{} {}'.format(self.name, self.parent_repl)

    def __eq__(self, other):
        if not isinstance(other, Segment):
            return False
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

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
        return random.pcg_random.integer(self.length)

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
        if random.pcg_random.integer(2) == 0:
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
        if random.pcg_random.uniform() < self.anti_rate:
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

    def __init__(self, name, cell, copy_number, anti_rate, linear):
        """
        The definition of a replicon (chromosome, plasmid, etc).
        :param name: a unique name for this replicon
        :param cell: the parent cell for this replicon
        :param copy_number: the copy-number of this replicon in 'cell'
        :param anti_rate: the rate of anti-diagonal interactions
        :param linear: treat replicons as linear
        """

        self.name = name
        self.copy_number = copy_number
        # segment registry
        self.segment_registry = OrderedDict()
        self.cdf_extent = None
        self.cdf_sites = None
        self.segment_list = None
        self.anti_rate = anti_rate
        self.linear = linear
        self.parent_cell = None

        # set bidirection association with containing cell
        if not isinstance(cell, Cell):
            raise TypeError('cell attribute must be of type Cell')
        self.parent_cell = cell
        cell.register_replicon(self)

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

    def __repr__(self):
        return '{} {}'.format(self.name, self.parent_cell)

    def __eq__(self, other):
        if not isinstance(other, Replicon):
            return False
        return self.name == other.name and self.parent_cell == other.parent_cell

    def __hash__(self):
        return hash(self.name) ^ hash(self.parent_cell)

    def init_prob(self):
        seg_info = self.get_segment_info()
        pdf_extent = np.array(seg_info['length'], dtype=f'{PROB_TYPE}')
        pdf_sites = np.array(seg_info['sites'], dtype=f'{PROB_TYPE}')
        self.cdf_extent = pdf2cdf(pdf_extent)
        self.cdf_sites = pdf2cdf(pdf_sites)
        self.segment_list = list(self.segment_registry)

    def num_segments(self):
        """
        :return: the number of segments for this replicon.
        """
        return len(self.segment_registry)

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
        if not isinstance(segment, Segment):
            raise TypeError('segment parameter must be of type Segment')
        if segment in self.segment_registry:
            raise Sim3CException('Duplicate segment instance. Segment names must be unique across the community')
        self.segment_registry[segment] = segment
        return segment

    def draw_any_segment_by_sites(self):
        """
        Draw any segment from the replicon, based on each segments number of sites
        :return: a Segment
        """
        return cdf_choice(self.segment_list, self.cdf_sites)

    def draw_any_segment_by_extent(self):
        """
        Draw any segment from the replicon, biased by the segment length.
        :return: a Segment
        """
        return cdf_choice(self.segment_list, self.cdf_extent)

    def print_report(self, show_cdf=False):
        """
        Print a simple report about this replicon.
        """
        print(f'segments {[si.name for si in self.segment_list]}')
        if show_cdf:
            print(f'cdf_extent: {self.cdf_extent}')
            print(f'cdf_sites: {self.cdf_sites}')

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
        self.cis_segment_indices = OrderedDict()
        self.trans_segment_indices = OrderedDict()
        self.replicon_list = None
        self.segment_list = None
        self.cdf_sites_inter = None
        self.cdf_extents_inter = None

    def __repr__(self):
        return repr((self.name, self.abundance, self.num_replicons()))

    def __eq__(self, other):
        if not isinstance(other, Cell):
            return False
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

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

        # array of names for extracting random sequences
        # self.replicon_names = np.array([*self.replicon_registry])
        self.replicon_list = list(self.replicon_registry)

        # begin with some empty PDFs
        pdf_extent = np.zeros(self.num_segments(), dtype=f'{PROB_TYPE}')
        pdf_sites = np.zeros(self.num_segments(), dtype=f'{PROB_TYPE}')

        # for each replicon, the PDF for the various modes of selection.
        i = 0
        self.cis_segment_indices = OrderedDict()
        for repl in self.replicon_list:
            repl.init_prob()
            cn = float(repl.copy_number)
            seg_info = repl.get_segment_info()
            self.cis_segment_indices[repl] = []
            for si in seg_info:
                pdf_extent[i] = cn * si['length']
                pdf_sites[i] = cn * si['sites']
                self.cis_segment_indices[repl].append(i)
                i += 1

        self.trans_segment_indices = OrderedDict()
        for repl_i in self.cis_segment_indices:
            self.trans_segment_indices[repl_i] = []
            for repl_j, seg_idx in self.cis_segment_indices.items():
                if repl_i != repl_j:
                    self.trans_segment_indices[repl_i].extend(seg_idx)

        # CDFs from the PDFs. Selection is accomplished by drawing a
        # unif([0..1]) and mapping that through the relevant CDF.
        self.cdf_extent = pdf2cdf(pdf_extent)
        self.cdf_sites = pdf2cdf(pdf_sites)

        # prepare a registry of all defined segments for this replicon
        for repl in self.replicon_registry.values():
            for seg in repl.segment_registry.values():
                self.segment_registry[seg] = seg
        self.segment_list = list(self.segment_registry)

        # One last set of CDFs for "select other" which exclude
        # each replicon in turn.
        if self.num_replicons() > 1:

            self.cdf_sites_inter = OrderedDict()
            self.cdf_extents_inter = OrderedDict()

            for repl_i in self.replicon_list:

                # indices without rn
                xi = self.trans_segment_indices[repl_i]

                # site probs without rn
                pi = pdf_sites[xi]
                self.cdf_sites_inter[repl_i] = {'replicons': [self.segment_list[i] for i in xi],
                                                'prob': pdf2cdf(pi)}

                # extent probs without rn
                pi = pdf_extent[xi]
                self.cdf_extents_inter[repl_i] = {'replicons': [self.segment_list[i] for i in xi],
                                                  'prob': pdf2cdf(pi)}

    def register_replicon(self, repl):
        """
        Insert a replicon into this cells registry
        :param repl: the replicon to insert
        :return: the inserted replicon
        """
        if not isinstance(repl, Replicon):
            raise TypeError('repl parameter must be of type Replicon')
        if repl not in self.replicon_registry:
            self.replicon_registry[repl] = repl
        return repl

    def draw_any_segment_by_extent(self):
        """
        Draw any segment, replicon tuple from this cell. The probability is biased by per-replicon
        genomic extent and copy number.
        :return: a segment,replicon tuple from this cell
        """
        return cdf_choice(self.segment_list, self.cdf_extent)

    def draw_any_segment_by_sites(self):
        """
        Draw any segment, replicon tuple from this cell. The probability is biased by per-replicon
        number of sites and copy number.
        :return: a segment, replicon tuple from this cell
        """
        return cdf_choice(self.segment_list, self.cdf_sites)

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

        return cdf_choice(self.cdf_sites_inter[skip_repl]['replicons'], self.cdf_sites_inter[skip_repl]['prob'])

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

        return cdf_choice(self.cdf_extents_inter[skip_repl]['replicons'], self.cdf_extents_inter[skip_repl]['prob'])

    def draw_any_site(self):
        """
        Draw a cut-site from any replicon within this cell. The probability
        of drawing a replicon is biased by the per-replicon number of sites
        and copy number, while genomic location is uniform.
        :return: a tuple of (segment, location)
        """
        segment = self.draw_any_segment_by_sites()
        return segment, segment.draw_any_site()

    def print_report(self, show_cdf=False):
        """
        Print a simple report about this cell.
        """
        print(f'replicons {[ri.name for ri in self.replicon_list]}')
        if show_cdf:
            print(f'cdf_extent: {self.cdf_extent}')
            print(f'cdf_sites: {self.cdf_sites}')

    def is_trans(self):
        """
        Coin toss test for whether an inter-replicon (trans) ligation product was formed. This
        is dictated by the rate supplied at instantiation time (trans_rate).
        :return: True -- treat this as a trans event
        """
        return self.num_replicons() > 1 and random.pcg_random.uniform() < self.trans_rate


class Community(object):
    """
    A community represents the entire collection and organisation of DNA molecules (replicons) in a simulation.
    This may be the approximation of an environmental sample, a multi-chromosomal or even monochromosomal
    organism.

    It is the entry point for the supporting reference data in a simulation, including such things as
    the relative abundance profile, DNA sequences and selected restriction enzyme. Additionally, are number
    of simulation parameters are exposed.
    """

    # defaults for communities
    ANTI_RATE = 0.2
    SPURIOUS_RATE = 0.02
    TRANS_RATE = 0.1
    LINEAR = False

    def __init__(self, enzyme, spurious_rate=0.02):
        """
        Initialise a community.

        :param enzyme: the enzyme used to digest DNA in the 3C/HiC library preparation
        :param spurious_rate: the rate of spurious ligation products
        """
        # global segment, replicon, and cell registries
        self.segm_registry = OrderedDict()
        self.repl_registry = OrderedDict()
        self.cell_registry = OrderedDict()

        # 3C protocol enzymatic digest
        self.enzyme = enzyme
        # intercellular rate is scaled by the product of all chrom site probs
        self.spurious_rate = spurious_rate

    @staticmethod
    def from_toml(seq_index, profile_toml, enzyme,
                  anti_rate=ANTI_RATE, spurious_rate=SPURIOUS_RATE,
                  trans_rate=TRANS_RATE, linear=LINEAR):
        """
        Construct a community from a TOML file resource.
        :param seq_index:
        :param profile_toml:
        :param enzyme:
        :param anti_rate:
        :param spurious_rate:
        :param trans_rate:
        :param linear:
        :return: Community object
        """

        comm_dict = read_toml(profile_toml, True)

        # override spurious rate if not defined in the TOML
        if 'spurious_rate' not in comm_dict['community']:
            logger.warning(f'No spurious_rate set for community, falling back to default {spurious_rate}')
            comm_dict['community']['spurious_rate'] = spurious_rate
        community = Community(enzyme, comm_dict['community']['spurious_rate'])
        for ci in comm_dict['community']['cells']:
            # override trans_rate if not defined in the TOML
            if 'trans_rate' not in ci:
                logger.warning(f'No trans_rate set for cell {ci["name"]}, falling back to default {trans_rate}')
                ci['trans_rate'] = trans_rate

            cell = community._register_cell(Cell(ci['name'], ci['abundance'], ci['trans_rate']))

            for ri in ci['replicons']:
                # override anti_rate and linear if not defined in the TOML
                if 'anti_rate' not in ri:
                    logger.warning(f'No anti_rate set for replicon {ri["name"]}, falling back to default {anti_rate}')
                    ri['anti_rate'] = anti_rate
                if 'linear' not in ri:
                    logger.warning(f'linear status set for replicon {ri["name"]}, falling back to default {linear}')
                    ri['linear'] = linear

                repl = Replicon(ri['name'], cell, ri['copy_number'], ri['anti_rate'], ri['linear'])

                for si in ri['segments']:
                    try:
                        seg = Segment(si, repl=repl, seq=seq_index[si].upper(), enzyme=enzyme)
                        community._register_segment(seg)
                    except NoCutSitesException:
                        logger.warning(f'Sequence {si} had no cut-sites '
                                       f'for the enzyme {str(enzyme)} and will be ignored')
                if repl.num_segments() > 0:
                    community._register_replicon(repl)
                    logger.debug(f'Replicon {repl.name} was added to community')
                else:
                    logger.warning(f'Replicon {ri["name"]} had no usable segments and will be ignored')
                    del cell.replicon_registry[repl]

        community._init_community()
        return community

    def to_toml(self, output_filename):
        """
        Serialize the community to a TOML file.
        :param output_filename:
        """
        import toml
        comm_dict = OrderedDict()
        comm_dict['community'] = {'spurious_rate': self.spurious_rate, 'cells': []}
        for cell_i in self.cell_registry.values():
            cd = {'name': cell_i.name, 'abundance': cell_i.abundance, 'trans_rate': cell_i.trans_rate, 'replicons': []}
            for repl_i in cell_i.replicon_registry.values():
                rd = {'name': repl_i.name, 'copy_number': repl_i.copy_number, 'anti_rate': repl_i.anti_rate,
                      'linear': repl_i.linear, 'segments': []}
                for seg in repl_i.segment_registry.values():
                    rd['segments'].append(seg.name)
                cd['replicons'].append(rd)
            comm_dict['community']['cells'].append(cd)
        with open(output_filename, 'wt') as output_h:
            encoder = toml.TomlArraySeparatorEncoder(separator='\n')
            toml.dump(comm_dict, output_h, encoder=encoder)

    @staticmethod
    def from_profile(seq_index, profile_table, enzyme,
                     anti_rate=ANTI_RATE, spurious_rate=SPURIOUS_RATE,
                     trans_rate=TRANS_RATE, linear=LINEAR):
        """
        Construct a community from a profile table file resource.
        :param seq_index:
        :param profile_table:
        :param enzyme:
        :param anti_rate:
        :param spurious_rate:
        :param trans_rate:
        :param linear:
        :return: Community object
        """

        comm_dict = read_profile(profile_table, True)
        community = Community(enzyme, spurious_rate)
        for c_name, c_details in comm_dict.items():
            # override trans_rate if not defined in the TOML
            cell = community._register_cell(Cell(c_name, c_details['abundance'], trans_rate))

            for r_name, r_details in c_details['replicons'].items():
                repl = Replicon(r_name, cell, r_details['copy_number'], anti_rate, linear)

                for s_name in r_details['segments']:
                    try:
                        seg = Segment(s_name, repl=repl, seq=seq_index[s_name].upper(), enzyme=enzyme)
                        community._register_segment(seg)
                    except NoCutSitesException:
                        logger.warning(f'Sequence {s_name} had no cut-sites '
                                       f'for the enzyme {str(enzyme)} and will be ignored')
                if repl.num_segments() > 0:
                    community._register_replicon(repl)
                    logger.debug(f'Replicon {repl.name} was added to community')
                else:
                    logger.warning(f'Replicon {r_name} had no usable segments and will be ignored')
                    del cell.replicon_registry[repl]

        community._init_community()
        return community

    def _init_community(self):
        """
        After the community is constructed from a file resource, initialise the probabilities
        """
        # now we're finished reading the profile, initialise the probabilities for each cell
        for cell in self.cell_registry.values():
            try:
                cell.init_prob()
            except NoRepliconsException:
                logger.warning('Cell "{}" had no usable replicons and will be ignored'.format(cell.name))
                del self.cell_registry[cell]

        # now initialise the probs for the whole community
        pdf_extent = np.zeros(len(self.segm_registry), dtype=f'{PROB_TYPE}')
        pdf_sites = np.zeros(len(self.segm_registry), dtype=f'{PROB_TYPE}')

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

        # derive cdfs from pdfs, these are what's used for drawing values
        self.cdf_extent = pdf2cdf(pdf_extent)
        self.cdf_sites = pdf2cdf(pdf_sites)

        # keep a list of names in numpy format
        self.repl_list = list(self.repl_registry)
        self.segm_list = list(self.segm_registry)

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
        if cell not in self.cell_registry:
            self.cell_registry[cell] = cell
        return self.cell_registry[cell]

    def _register_replicon(self, repl):
        """
        Add an instance of Replicon to the replicon registry
        :param repl: a Replicon object to add
        :return: the added replicon instance
        """
        assert isinstance(repl, Replicon), 'Attempted to register invalid class.'
        if repl not in self.repl_registry:
            self.repl_registry[repl] = repl
        return self.repl_registry[repl]

    def _register_segment(self, segment):
        """
        Add an instance of Segment to the segment registry
        :param segment:
        :return: the added segment instance
        """
        assert isinstance(segment, Segment), 'Attempted to register invalid class.'
        if segment in self.segm_registry:
            raise Sim3CException('duplicate segment names in community')
        self.segm_registry[segment] = segment
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
        return cdf_choice(self.segm_list, self.cdf_extent)

    def draw_any_segment_by_sites(self):
        """
        Draw any replicon from this community. The probability is biased by cellular abundance,
        and per-replicon number of cut-sites and copy number.
        :return: any replicon from this community
        """
        return cdf_choice(self.segm_list, self.cdf_sites)

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

    def print_report(self, show_cdf=False):
        print(f'segments {[si.name for si in self.segm_list]}')
        if show_cdf:
            print(f'cdf_extent: {self.cdf_extent}')
            print(f'cdf_sites: {self.cdf_sites}')

    def is_spurious(self):
        """
        Coin toss test for whether a spurious ligation product was formed. This is
        dictated by the rate supplied at instantiation time (spurious_rate).
        :return: True -- treat this as a spurious event
        """
        return random.pcg_random.uniform() < self.spurious_rate
