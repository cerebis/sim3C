import logging
import numpy as np

from collections import OrderedDict

from .empirical_model import EmpiricalDistribution, generate_nested_cids, cids_to_blocks, cdf_geom_unif_ratio
from .exceptions import *
from .site_analysis import AllSites, CutSites

logger = logging.getLogger(__name__)


def choice(rs, vals, cdf):
    """
    Random selection of an element from an array, biased by the supplied CDF.
    :param rs: a numpy RandomState object for random number generation
    :param vals: the array of potential choices
    :param cdf: the CDF describing each elements probability of selection
    :return: the selected element
    """
    return vals[np.searchsorted(cdf, rs.uniform())]


class Replicon:
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

    def __init__(self, name, cell, cn, seq, enzyme, anti_rate, random_state, create_cids=True, linear=False):
        """
        The definition of a replicon (chromosome, plasmid, etc).
        :param name: a unique name for this replicon
        :param cell: the parent cell for this replicon
        :param cn: the copy-number of this replicon in 'cell'
        :param seq: the genomic sequence of this replicon as a Bio.Seq object
        :param enzyme: the enzyme used to digest this replicon as a Bio.Restriction RestrictionType
        :param anti_rate: the rate of anti-diagonal interactions
        :param random_state: a numpy RandomState object for random number generation
        :param create_cids: when true, simulate chromosome-interacting-domains
        :param linear: treat replicon as linear
        """

        # if random state not supplied, initialise one
        self.random_state = random_state
        self.choice = random_state.choice
        self.uniform = random_state.uniform
        self.randint = random_state.randint
        self.linear = linear

        self.name = name
        self.copy_number = cn
        self.seq = seq

        if self.linear:
            if anti_rate > 0:
                logger.warning('Replicon {} is linear, anti_rate of {} has been set to zero'
                               .format(name, anti_rate))
            self.anti_rate = 0
        else:
            self.anti_rate = anti_rate

        # cut-site related properties. These are pre-calculated as a simple
        # means of avoiding performance penalties with repeated calls.
        if not enzyme:
            self.sites = AllSites(len(seq.seq), self.random_state)
        else:
            self.sites = CutSites(enzyme, seq.seq, self.random_state, linear=linear)

        self.length = len(self.seq)
        self.num_sites = self.sites.size
        self.site_density = self.num_sites / float(self.length)

        if create_cids:
            # setup for more complex simulated CID model
            self.draw_constrained_site = self._draw_cid_constrained_site
            self.cid_blocks = cids_to_blocks(
                generate_nested_cids(self.random_state, self.length, Replicon.BACKBONE_PROB,
                                     Replicon.GLOBAL_EMPDIST_BINS, Replicon.GLOBAL_SHAPE_FACTOR,
                                     Replicon.CID_EMPDIST_BINS, Replicon.CID_SHAPE_FACTOR,
                                     cdf_alpha=Replicon.CDF_ALPHA,
                                     min_num=Replicon.CID_MIN, max_num=Replicon.CID_MAX,
                                     recur_depth=Replicon.CID_DEPTH))

        else:
            # setup for simple model
            self.draw_constrained_site = self._draw_simple_constrained_site
            self.empdist = EmpiricalDistribution(self.random_state, self.length,
                                                 Replicon.GLOBAL_EMPDIST_BINS, cdf_geom_unif_ratio,
                                                 cdf_alpha=Replicon.CDF_ALPHA,
                                                 shape=Replicon.GLOBAL_SHAPE_FACTOR)

        # set bidirection association with containing cell
        self.parent_cell = cell
        cell.register_replicon(self)

    def __repr__(self):
        return '{} {} {} {}'.format(self.name, self.parent_cell, self.sites.size, self.length)

    def covers_site(self, x, length):
        """
        Test if the segment defined by a position and length cover any cut-sites
        for this replicon and enzyme digestion.
        :param x: the beginning position
        :param length: the segment length
        :return: True - this segment covers at least one site.
        """
        return self.sites.covers(x, length)

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
        return self.randint(0, self.length)

    @staticmethod
    def get_loc_3c(emp_dist, x1, length):
        """
        Get a second genomic location (x2) on this replicon, where the separation |x2-x1|
        is constrained by the experimentally determined distribution for 3C/HiC ligation
        products.
        Note TODO: this method does not explicitly treat linear replicons. To do so, we must
        handle the edge case of drawing beyond first and last positions. Simple approaches
        using redraw are potentially computationally terrible if the first position is
        very close to the end of the replicon. For now, we accept the modulo solution
        for both linear and circular.

        :param emp_dist: empirical distribution of separation
        :param x1: the first position
        :param length: the length (bp) of the replicon
        :return: x2: the second position
        """
        # draw a random separation
        delta = int(emp_dist.rand())

        # pve or nve shift, modulo length
        if emp_dist.uniform() < 0.5:
            x2 = (x1 - delta) % length
        else:
            x2 = (x1 + delta) % length
        return x2

    def _draw_simple_constrained_site(self, x1):
        """
        Draw a second site (x2) relative to the first (x1) which is constrained
        to follow a basic empirical distribution.
        :param x1: the first location
        :return: x2: the second location
        """
        x2 = Replicon.get_loc_3c(self.empdist, x1, self.length)

        # anti-diagonal
        if self.uniform() < self.anti_rate:
            x2 = self.length - x2

        # return nearest site
        return self.sites.find_nn(x2)

    def _draw_cid_constrained_site(self, x1):
        """
        Draw a second site (x2) relative to teh first (x1) which is constrained to follow
        a nested hierarchy of empirical distributions, intended to simulate chromosomal
        interaction domains (CID). Here, CIDs are small regions of a chromosome which interact
        (possibly through folding) at a more frequent level than that of the overall backbone.

        The range of effect of the CIDs are predetermined at instantiation time and stored as
        intervals within an interval-tree. When called, x1 determines which overlapping intervals
        are involved and single emp-dist (representing a particular CID) is chosen at random
        (randomly weighted).

        :param x1: the first location
        :return: x2: the second location.
        """
        block = self.cid_blocks[x1].pop()
        ovl_invs = block.data['inv_list']

        if len(ovl_invs) == 1:
            # only the background distribution governs this block
            chosen_inv = ovl_invs[0]

            x2 = Replicon.get_loc_3c(chosen_inv.data['empdist'], x1, chosen_inv.length())

            # anti-diagonal
            if self.uniform() < self.anti_rate:
                x2 = self.length - x2

        else:
            # pick a cid or background from those defined for this block
            # note, numpy will not accept the list of intervals here
            ix = self.choice(len(ovl_invs), p=block.data['prob_list'])
            chosen_inv = ovl_invs[ix]

            x2 = Replicon.get_loc_3c(chosen_inv.data['empdist'], x1 - chosen_inv.begin, chosen_inv.length())
            x2 += chosen_inv.begin

            # main backbone gets anti-diagonal treatment
            if chosen_inv.data['depth'] == 0:
                # antidiagonal
                if self.uniform() < self.anti_rate:
                    x2 = self.length - x2

        return self.sites.find_nn(x2)

    def subseq(self, x1, length, rev=False):
        """
        Create a subsequence.

        :param x1: starting genomic position
        :param length: length of subsequence
        :param rev: reverse complement this sequence.
        :return: subseq Seq object
        """

        # handle negative starts as wrapping around or,
        # in linear cases, terminates at first position
        if x1 < 0:
            x1 = x1 % self.length if not self.linear else 0

        x2 = x1 + length
        diff = x2 - self.length
        if diff > 0:
            if self.linear:
                # sequence terminates early
                ss = self.seq[x1:-1]
                ss.description = Replicon.PART_DESC_FMT.format(rev, x1+1, self.length)
            else:
                # sequence will wrap around
                ss = self.seq[x1:] + self.seq[:diff]
                ss.description = Replicon.PART_DESC_FMT.format(rev, x1+1, diff)
        else:
            ss = self.seq[x1:x2]
            ss.description = Replicon.PART_DESC_FMT.format(rev, x1+1, x2)

        if rev:
            ss.reverse_complement(id=True, description=True)

        return ss


class Cell:
    """
    A cell acts as the container of one or more replicons, where each may have its
    own copy-number in addition to the relative abundance of their containing cell.
    """

    def __init__(self, name, abundance, random_state, trans_rate=0.1):
        """
        The definition of a cell.
        :param name: a unique name for this cell.
        :param abundance: the relative abundance of this cell in the community.
        :param random_state: a numpy RandomState object for random number generation
        :param trans_rate: the rate of inter-replicon (trans) ligation products
        """

        # if no object supplied, initialise one.
        self.random_state = random_state
        self.uniform = random_state.uniform
        self.name = name
        self.abundance = abundance
        # replicons are kept in the order they are registered
        self.replicon_registry = OrderedDict()
        # inter-replicon (trans) rate
        self.trans_rate = trans_rate

        # probability properties, to be initialised by init_prob()
        self.cdf = None
        self.pdf = None
        self.pdf_cn = None
        self.pdf_extent = None
        self.pdf_sites = None
        self.cdf_cn = None
        self.cdf_extent = None
        self.cdf_sites = None
        self.replicon_names = None
        self.cdf_sites_inter = None
        self.cdf_extents_inter = None

    def __repr__(self):
        return repr((self.name, self.abundance, self.num_replicons()))

    def num_replicons(self):
        """
        :return: the number of replicons for this cell.
        """
        return len(self.replicon_registry)

    def init_prob(self):
        """
        Initialise the selection probabilities. This method should be called after all replicons are
        registered or when a new replicon is added.
        """
        if self.num_replicons() <= 0:
            raise NoRepliconsException(self.name)

        # begin with some empty PDFs
        self.pdf_cn = np.zeros(self.num_replicons())
        self.pdf_extent = np.zeros(self.num_replicons())
        self.pdf_sites = np.zeros(self.num_replicons())

        # for each replicon, the PDF for the various modes of selection.
        for i, k in enumerate(self.replicon_registry):
            repl = self.replicon_registry[k]
            cn = float(repl.copy_number)
            self.pdf_cn[i] = cn
            self.pdf_extent[i] = cn * repl.length
            self.pdf_sites[i] = cn * repl.num_sites

        # normalise the PDFs
        self.replicon_names = np.array(self.replicon_registry.keys())
        self.pdf_cn /= self.pdf_cn.sum()
        self.pdf_extent /= self.pdf_extent.sum()
        self.pdf_sites /= self.pdf_sites.sum()

        # CDFs from the PDFs. Selection is accomplished by drawing a
        # unif([0..1]) and mapping that through the relevant CDF.
        self.cdf_cn = np.cumsum(self.pdf_cn)
        self.cdf_extent = np.cumsum(self.pdf_extent)
        self.cdf_sites = np.cumsum(self.pdf_sites)

        # One last set of CDFs for "select other" which exclude
        # each replicon in turn.
        if self.num_replicons() > 1:

            self.cdf_sites_inter = {}
            self.cdf_extents_inter = {}

            for ci in self.replicon_names:
                # indices without ci
                xi = np.where(self.replicon_names != ci)

                # site probs without ci
                pi = self.pdf_sites[xi]
                pi /= pi.sum()
                self.cdf_sites_inter[ci] = {'names': self.replicon_names[xi], 'prob': np.cumsum(pi)}

                # extent probs without ci
                pi = self.pdf_extent[xi]
                pi /= pi.sum()
                self.cdf_extents_inter[ci] = {'names': self.replicon_names[xi], 'prob': np.cumsum(pi)}

    def register_replicon(self, repl):
        """
        Insert a replicon into this cells registry
        :param repl: the replicon to insert
        :return: the inserted replicon
        """
        assert isinstance(repl, Replicon),  'Attempted to register invalid class.'
        if repl.name in self.replicon_registry:
            raise Sim3CException('duplicate replicon names')
        self.replicon_registry[repl.name] = repl
        return repl

    def get_replicon(self, name):
        """
        Return a replicon from the registry by name.
        :param name: the name of the replicon
        :return: the replicon
        """
        return self.replicon_registry[name]

    def draw_replicon(self):
        """
        Draw any replicon from this cell, biased only by copy number.
        :return: any replicon from this cell
        """
        return self.get_replicon(choice(self.random_state, self.replicon_names, self.cdf_cn))

    def draw_any_replicon_by_extents(self):
        """
        Draw any replicon from this cell. The probability is biased by per-replicon
        genomic extent and copy number.
        :return: any replicon from this cell
        """
        return self.get_replicon(choice(self.random_state, self.replicon_names, self.cdf_extent))

    def draw_any_replicon_by_sites(self):
        """
        Draw any replicon from this cell. The probability is biased by per-replicon
        number of sites and copy number.
        :return: any replicon from this cell
        """
        return self.get_replicon(choice(self.random_state, self.replicon_names, self.cdf_sites))

    def draw_other_replicon_by_sites(self, skip_repl):
        """
        Draw a different replicon from this cell. The probability is biased by per-replicon
        number of sites and copy number. Note: single replicon cell definitions will
        raise an exception.
        :param skip_repl: the replicon to exclude
        :return: another replicon from this cell
        """
        if self.num_replicons() <= 1:
            raise MonochromosomalException('inter-replicon events are not possible for monochromosomal cells')

        return self.get_replicon(choice(self.random_state,
                                        self.cdf_sites_inter[skip_repl]['names'],
                                        self.cdf_sites_inter[skip_repl]['prob']))

    def draw_other_replicon_by_extents(self, skip_repl):
        """
        Draw a different replicon from this cell. The probability is biased by per-replicon
        extent and copy number. Note: single replicon cell definitions will
        raise an exception.
        :param skip_repl: the replicon to exclude
        :return: another replicon from this cell
        """
        if self.num_replicons() <= 1:
            raise MonochromosomalException('inter-replicon events are not possible for monochromosomal cells')

        return self.get_replicon(choice(self.random_state,
                                        self.cdf_extents_inter[skip_repl]['names'],
                                        self.cdf_extents_inter[skip_repl]['prob']))

    def draw_any_site(self):
        """
        Draw a cut-site from any replicon within this cell. The probability
        of drawing a replicon is biased by the per-replicon number of sites
        and copy number, while genomic location is uniform.
        :return: a tuple of (replicon, location)
        """
        repl = self.draw_any_replicon_by_sites()
        return repl, repl.draw_any_site()

    def print_report(self):
        """
        Print a simple report about this cell.
        """
        print 'names', self.replicon_names
        print 'p_rep', self.pdf_cn
        print 'p_ext', self.pdf_extent
        print 'p_sit', self.pdf_sites

    def is_trans(self):
        """
        Coin toss test for whether an inter-replicon (trans) ligation product was formed. This
        is dictated by the rate supplied at instantiation time (trans_rate).
        :return: True -- treat this as a trans event
        """
        return self.num_replicons() > 1 and self.uniform() < self.trans_rate


class Community:
    """
    A community represents the entire collection and organisation of DNA molecules (replicons) in a simulation.
    This may be the approximation of an environmental sample, a multi-chromosomal or even monochromosomal
    organism.

    It is the entry point for the supporting reference data in a simulation, including such things as
    the relative abundance profile, DNA sequences and selected restriction enzyme. Additionally, are number
    of simulation parameters are exposed.
    """

    def __init__(self, seq_index, profile, enzyme, random_state, anti_rate=0.2, spurious_rate=0.02,
                 trans_rate=0.1, create_cids=True, linear=False):
        """
        Initialise a community.

        :param seq_index: an open sequence index of type _IndexedSeqFileDict as returned from Bio.SeqIO.index
        :param profile: the accompanying abundance profile of all replicons in the community
        :param enzyme: the enzyme used to digest DNA in the 3C/HiC library preparation
        :param anti_rate: the rate of anti-diagonal interactions
        :param random_state: a numpy random state used for random number generation
        :param spurious_rate: the rate of spurious ligation products
        :param trans_rate: the rate of inter-replicon (trans) ligation products within a cell
        :param create_cids: when true, simulate chromosome-interacting-domains
        :param linear: treat replicons as linear
        """

        # init a random state if one was not supplied.
        # keep a reference to uniform handy
        self.random_state = random_state
        self.uniform = random_state.uniform

        # global replicon and cell registries
        self.repl_registry = OrderedDict()
        self.cell_registry = OrderedDict()

        # initialise the registries using the community profile
        for ri in profile.values():
            # register the cell
            cell = self._register_cell(Cell(ri.cell, ri.abundance, self.random_state, trans_rate))
            # fetch the sequence from file
            rseq = seq_index[ri.name].upper()
            try:
                # community-wide replicon registry
                self._register_replicon(Replicon(ri.name, cell, ri.copy_number, rseq, enzyme,
                                                 anti_rate, random_state, create_cids, linear))
            except NoCutSitesException as ex:
                logger.warning('Sequence "{}" had no cut-sites for the enzyme {} and will be ignored'.format(
                    ri.name, str(enzyme)))

        # now we're finished reading replicons, initialise the probs for each cell
        for cell in self.cell_registry.values():
            try:
                cell.init_prob()
            except NoRepliconsException as ex:
                logger.warning('Cell "{}" had no usable replicons and will be ignored'.format(cell.name))
                del self.cell_registry[cell.name]

        # now initialise the probs for the whole community
        self.pdf_repl = np.zeros(len(self.repl_registry))
        self.pdf_extent = np.zeros(len(self.repl_registry))
        self.pdf_sites = np.zeros(len(self.repl_registry))

        # whether site (prox-lig) or extent (wgs) based, probs are
        # weighted by cellular abundance and copy number
        for i, k in enumerate(self.repl_registry):
            repl = self.repl_registry[k]
            abn = repl.parent_cell.abundance
            cn = repl.copy_number
            self.pdf_repl[i] = abn * cn
            self.pdf_extent[i] = abn * cn * repl.length
            self.pdf_sites[i] = abn * cn * repl.num_sites

        # now normalise pdfs
        self.pdf_repl /= self.pdf_repl.sum()
        self.pdf_extent /= self.pdf_extent.sum()
        self.pdf_sites /= self.pdf_sites.sum()

        # derive cdfs from pdfs, these are what's used for drawing values
        self.cdf_repl = np.cumsum(self.pdf_repl)
        self.cdf_extent = np.cumsum(self.pdf_extent)
        self.cdf_sites = np.cumsum(self.pdf_sites)

        # keep a list of chr names in numpy format
        self.repl_names = np.array(self.repl_registry.keys())

        # setup pdfs and cdfs for inter-cellular events
        # each represents a deletion of one chr and renormalisation
        if len(self.repl_registry) > 1:
            self.cdf_sites_inter = {}
            for rni in self.repl_names:
                ix = np.where(self.repl_names != rni)
                pi = self.pdf_sites[ix]
                pi /= pi.sum()
                self.cdf_sites_inter[rni] = {'names': self.repl_names[ix], 'prob': np.cumsum(pi)}

        # inter-cellular rate is scaled by the product of all chrom site probs
        self.spurious_rate = spurious_rate

        # keep the number of cells handy
        self.num_cells = len(self.cell_registry)

        # we must have at least one cell defined
        if self.num_cells <= 0:
            logger.error('Community did not contain any useful cells. '
                         'Check that your reference sequence(s) were '
                         'valid and that they each contain at least one cut-site.')
            raise Sim3CException('Community is empty')

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
        if repl.name in self.repl_registry:
            raise Sim3CException('duplicate replicon names in community')
        self.repl_registry[repl.name] = repl
        return repl

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

    def draw_repl(self):
        """
        Draw any replicon from this community, biased relative abundance and copy number.
        :return: any replicon from this community
        """
        return self.get_repl(choice(self.random_state, self.repl_names, self.cdf_repl))

    def draw_any_repl_by_extent(self):
        """
        Draw any replicon from this community. The probability is biased by cellular abundance,
        and per-replicon genomic extent and copy number.
        :return: any replicon from this community
        """
        return self.get_repl(choice(self.random_state, self.repl_names, self.cdf_extent))

    def draw_any_repl_by_sites(self):
        """
        Draw any replicon from this community. The probability is biased by cellular abundance,
        and per-replicon number of cut-sites and copy number.
        :return: any replicon from this community
        """
        return self.get_repl(choice(self.random_state, self.repl_names, self.cdf_sites))

    def draw_other_repl_by_sites(self, skip_repl):
        """
        Draw a different replicon from this community. The probability is biased by cellular abundance,
        per-replicon number of sites and copy number. Note: single replicon cell definitions will
        raise an exception.
        :param skip_repl: the replicon to exclude
        :return: another replicon from this cell
        """
        return self.get_repl(choice(self.random_state,
                                    self.cdf_sites_inter[skip_repl]['names'],
                                    self.cdf_sites_inter[skip_repl]['prob']))

    def draw_any_by_site(self):
        """
        Draw any site from any replicon, biased by abundance, number of sites and copy number.
        :return:
        """
        repl = self.draw_any_repl_by_sites()
        return repl, repl.draw_any_site()

    def draw_any_by_extent(self):
        """
        Draw any site from any replicon, biased by abundance, number of sites and copy number.
        :return:
        """
        repl = self.draw_any_repl_by_extent()
        return repl, repl.draw_any_location()

    def print_report(self):
        print 'names', self.repl_names
        print 'p_rep', self.pdf_repl
        print 'p_ext', self.pdf_extent
        print 'p_sit', self.pdf_sites

    def is_spurious(self):
        """
        Coin toss test for whether a spurious ligation product was formed. This is
        dictated by the rate supplied at instantiation time (spurious_rate).
        :return: True -- treat this as a spurious event
        """
        return self.uniform() < self.spurious_rate
