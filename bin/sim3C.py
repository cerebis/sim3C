#!/usr/bin/env python
"""
meta-sweeper - for performing parametric sweeps of simulated
metagenomic sequencing experiments.
Copyright (C) 2016 "Matthew Z DeMaere"

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from collections import OrderedDict, namedtuple

import numpy as np
import tqdm
from Bio import Alphabet
from Bio import SeqIO
from Bio.Restriction import Restriction
from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict

import Art
import abundance as abn
import empirical_model as em


class Sim3CException(Exception):
    """Module base exception class"""
    def __init__(self, message):
        super(Sim3CException, self).__init__(message)


class NoCutSitesException(Sim3CException):
    """Occurs when a target template contains no cutsites for a specified restriction enzyme"""
    def __init__(self, seq_name, enz_name):
        super(NoCutSitesException, self).__init__(
            'sequence [{0}] had no cutsites for enzyme [{1}]'.format(seq_name, enz_name))


class OutOfBoundsException(Sim3CException):
    """Raised when coordinates lie out of range of replicon"""
    def __init__(self, pos, maxpos):
        super(OutOfBoundsException, self).__init__(
            "exceeded maximum template length {0} > {1}".format(pos, maxpos))


class EmptyRegistryException(Sim3CException):
    """No registry was empty when attempting to act upon its contents"""
    pass


class MonochromosomalException(Sim3CException):
    """A method require more than one chromosome was invoked on a monochromosomal cell"""
    pass


def choice(rs, vals, cdf):
    """
    Random selection of an element from an array, biased by the supplied CDF.
    :param rs: a numpy RandomState object for random number generation
    :param vals: the array of potential choices
    :param cdf: the CDF describing each elements probability of selection
    :return: the selected element
    """
    return vals[np.searchsorted(cdf, rs.uniform())]


def get_enzyme_instance_ipython(enz_name):
    """
    An alternative method to fetch an instance of a given restriction enzyme by its
    name using a work-around which avoids exception with getattr() in iPython.

    Ackn: StackOverflow: user xbello.
    See: http://stackoverflow.com/questions/20381912/type-object-restrictiontype-has-no-attribute-size

    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    r_type_names = [rt for tid, (rt, rn) in typedict.iteritems() if enz_name in rn][0]
    r_clz = tuple(getattr(Restriction, rt) for rt in r_type_names)
    return Restriction.AbstractCut(enz_name, r_clz, rest_dict[enz_name])


def get_enzyme_instance(enz_name):
    """
    Fetch an instance of a given restriction enzyme by its name.
    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    return getattr(Restriction, enz_name)


class CutSites:
    """
    The cut-sites for a given enzyme on a given template.
    Note: locations are 0-based
    """

    def __init__(self, enzyme, template_seq, random_state, linear=False):
        """
        Initialize the cut-sites of an enzyme on a template sequence.

        Linearity affects both the initial search for recognition sites and
        when requesting the nearest site to a given genomic location.

        :param enzyme: the restriction enzyme (Bio.Restriction RestrictionType object)
        :param template_seq: the template sequence to digest (Bio.Seq object)
        :param random_state: the random state used for draws
        :param linear: treat sequence as linear
        """
        self.random_state = random_state
        self.choice = random_state.choice
        self.randint = random_state.randint

        self.enzyme = enzyme
        self.max_length = len(template_seq) - 1

        # find sites, converting from 1-based.
        self.sites = np.array(enzyme.search(template_seq, linear)) - 1
        self.size = self.sites.shape[0]
        if self.size == 0:
            raise NoCutSitesException(template_seq.id, str(enzyme))

        # method setup
        if linear:
            self.find_nn = self._find_nn_linear
            self.covers = self._covers_site_linear
        else:
            self.find_nn = self._find_nn_circular
            self.covers = self._covers_site_circular
            self.shouldered = self._add_shoulders()

    def _add_shoulders(self):
        """
        Add shoulders to the already determined list of circular chr sites. This allows
        finding positions without logic for edge cases (ix=0 or -1)
        """
        before_first = self.sites[-1] - self.max_length
        after_last = self.max_length + self.sites[0]
        return np.hstack(([before_first], self.sites, [after_last]))

    def random_site(self):
        """
        Select a uniformly random site
        :return: a random site
        """
        return self.sites[self.randint(self.size)]

    def _covers_site_linear(self, x1, length):
        """
        Test whether a specified position and length contains any cut-sites, assuming
        a linear molecule.
        :param x1: starting position
        :param length: length
        :return: True the coordinates contain at least one cut-site
        """
        return ((self.sites > x1) & (self.sites <= x1+length)).any()

    def _covers_site_circular(self, x1, length):
        """
        Test whether a specified position and length contains any cut-sites, assuming
        a circular molecule.
        :param x1: starting position
        :param length: length
        :return: True the coordinates contain at least one cut-site
        """
        x2 = x1 + length
        if x2 > self.max_length:
            ret = ((x2 % self.max_length) >= self.sites) | (x1 <= self.sites)
        else:
            ret = (x1 <= self.sites) & (x2 >= self.sites)
        return ret.any()

    def _find_nn_circular(self, pos):
        """
        Find the nearest cut-site relative to supplied position on a circular
        chromosome.
        :param pos: the position on the chromosome
        :return: nearest site
        """
        if pos > self.max_length:
            OutOfBoundsException(pos, self.max_length)

        cs = self.shouldered
        ix = np.searchsorted(cs, pos)
        x1 = cs[ix - 1]
        x2 = cs[ix]
        # modulo so as we only return values within the range [0..maxlen]
        # this handles the edge case of sites crossing beginning or end.
        if pos - x1 <= x2 - pos:
            return x1 % self.max_length
        else:
            return x2 % self.max_length

    def _find_nn_linear(self, pos):
        """
        Find the nearest cut-site relative to the supplied position on a linear
        chromosome or sequence fragment.
        :param pos: the position on the chromosome
        :return: nearest site
        """
        if pos > self.max_length:
            OutOfBoundsException(pos, self.max_length)

        cs = self.sites
        ix = np.searchsorted(cs, pos)
        # first or last site was closest
        if ix == 0:
            return cs[0]
        elif ix == self.size:
            return cs[1]
        else:
            # pick the closest of nearest neighbours
            x1 = cs[ix - 1]
            x2 = cs[ix]
            if pos - x1 <= x2 - pos:
                return x1
            else:
                return x2

    def find_first(self, pos, origin_site):
        """
        Beginning from pos and searching toward the origin, return the first encountered cut-site.
        This may turn out to be the origin if there is no other intervening site.
        :param pos: a position on the chromosome
        :param origin_site: a valid cut-site used as the relative origin
        :return: the first site encountered between pos and origin_site
        """
        if pos > self.max_length:
            OutOfBoundsException(pos, self.max_length)

        cs = self.sites
        if pos > origin_site:
            # the position begins after the origin
            ix = np.searchsorted(cs, pos, side='right')
            if ix == 0:
                return cs[ix]
            else:
                return cs[ix-1]
        else:
            # position begins before the origin site
            ix = np.searchsorted(cs, pos, side='left')
            return cs[ix]


class AllSites:

    def __init__(self, size, random_state):
        """
        Initialise an instance of AllSites. This represents a completely unconstrained
        model, where any base-position is equally accessible.

        :param size: the maximum length of a sequence
        :param random_state: random state used for draws
        """
        self.random_state = random_state
        self.randint = random_state.randint
        self.size = size

    @staticmethod
    def covers():
        """
        Always true, every position is covered.
        :return: True
        """
        return True

    def random_site(self):
        """
        Draw a random base-position uniformly over the entire extent (0..size-1).
        :return: random base position (0-based)
        """
        return self.randint(self.size)

    @staticmethod
    def find_nn(pos):
        """
        As every site is accessible this is effectively a dummy function returning the input.
        :param pos: the position to find a NN for
        :return: pos, unchanged.
        """
        return pos


class Replicon:
    """
    A DNA molecule which will be digested. This may be a chromosome, plasmid, etc.
    """

    # Formatting string for sequence descriptions of a part (subsequence) of a replicon.
    # Used in insert/read creation.
    PART_DESC_FMT = '{0:d}:{1}:{2}'

    # Empirical distribution parameters. These might eventually be exposed to users.
    BACKBONE_PROB = 0.25
    GLOBAL_EMPDIST_BINS = 1000
    GLOBAL_SHAPE_FACTOR = 4.0e-6
    CID_EMPDIST_BINS = 100
    CID_SHAPE_FACTOR = 8.0e-6
    CID_MIN = 3
    CID_MAX = 6
    CID_DEPTH = 2

    def __init__(self, name, cell, cn, seq, enzyme, anti_rate, random_state, create_cids=True):
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
        """

        # if random state not supplied, initialise one
        self.random_state = random_state
        self.choice = random_state.choice
        self.uniform = random_state.uniform
        self.randint = random_state.randint

        self.name = name
        self.copy_number = cn
        self.seq = seq
        self.anti_rate = anti_rate

        # cut-site related properties. These are pre-calculated as a simple
        # means of avoiding performance penalties with repeated calls.
        if not enzyme:
            self.sites = AllSites(len(seq.seq), self.random_state)
        else:
            self.sites = CutSites(enzyme, seq.seq, self.random_state, linear=False)

        self.length = len(self.seq)
        self.num_sites = self.sites.size
        self.site_density = self.num_sites / float(self.length)

        if create_cids:
            # setup for more complex simulated CID model
            self.draw_constrained_site = self._draw_cid_constrained_site
            self.cid_blocks = em.cids_to_blocks(
                em.generate_nested_cids(self.random_state, self.length, Replicon.BACKBONE_PROB,
                                        Replicon.GLOBAL_EMPDIST_BINS, Replicon.GLOBAL_SHAPE_FACTOR,
                                        Replicon.CID_EMPDIST_BINS, Replicon.CID_SHAPE_FACTOR,
                                        min_num=Replicon.CID_MIN, max_num=Replicon.CID_MAX,
                                        recur_depth=Replicon.CID_DEPTH))

        else:
            # setup for simple model
            self.draw_constrained_site = self._draw_simple_constrained_site
            self.empdist = em.EmpiricalDistribution(self.random_state, self.length,
                                                    Replicon.GLOBAL_EMPDIST_BINS, em.cdf_geom_unif,
                                                    shape=Replicon.GLOBAL_SHAPE_FACTOR)

        # set bidirection association with containing cell
        self.parent_cell = cell
        cell.register_replicon(self)

    def __repr__(self):
        return '{0} {1} {2} {3}'.format(self.name, self.parent_cell, self.sites.size, self.length)

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
        Create a subsequence, where the replicon is always treated as circular.

        :param x1: starting genomic position
        :param length: length of subsequence
        :param rev: reverse complement this sequence.
        :return: subseq Seq object
        """

        # handle negative starts as wrapping around.
        if x1 < 0:
            x1 %= self.length

        x2 = x1 + length
        diff = x2 - self.length
        if diff > 0:
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
        assert self.num_replicons() > 0, 'group contained no registered replicons'

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

    def __init__(self, seq_file, profile, enzyme, random_state, anti_rate=0.2, spurious_rate=0.02,
                 trans_rate=0.1, create_cids=True):
        """
        Initialise a community.

        :param seq_file: the multi-fasta sequences for all replicons in the community
        :param profile: the accompanying abundance profile of all replicons in the community
        :param enzyme: the enzyme used to digest DNA in the 3C/HiC library preparation
        :param anti_rate: the rate of anti-diagonal interactions
        :param random_state: a numpy random state used for random number generation
        :param spurious_rate: the rate of spurious ligation products
        :param trans_rate: the rate of inter-replicon (trans) ligation products within a cell
        :param create_cids: when true, simulate chromosome-interacting-domains
        """

        # init a random state if one was not supplied.
        # keep a reference to uniform handy
        self.random_state = random_state
        self.uniform = random_state.uniform

        # global replicon and cell registeries
        self.repl_registry = OrderedDict()
        self.cell_registry = OrderedDict()

        # reference fasta will be accessed by index.
        seq_index = SeqIO.index(seq_file, 'fasta', alphabet=Alphabet.generic_dna)

        # initialise the registries using the community profile
        for ri in profile.values():
            # register the cell
            cell = self._register_cell(Cell(ri.cell, ri.abundance, self.random_state, trans_rate))
            try:
                rseq = seq_index[ri.name]
            except Exception:
                raise Sim3CException('Error getting sequence {0} from fasta file'.format(ri.name))
            # community-wide replicon registry
            self._register_replicon(Replicon(ri.name, cell, ri.copy_number, rseq, enzyme, anti_rate,
                                             random_state, create_cids))

        # now we're finished reading replicons, initialise the probs for each cell
        for cell in self.cell_registry.values():
            cell.init_prob()

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
        assert self.num_cells > 0, 'Community appears to be empty'

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


class ReadGenerator:
    """
    Generate inserts and subsequent read-pairs for a particular library preparation method. Primarily,
    3C does not generate a duplication of the cut-site, whereas HiC's enriching for ligation products by
    use of biotinylation does produce site duplication during infill of overhangs.

    The desired insert size and its variability are specified here.

    Other sequencing read simulation parameters are supplied here initialise ART.
    """
    def __init__(self, method, enzyme, seed, random_state,
                 prefix='SIM3C', simple=False, machine_profile='EmpMiSeq250',
                 read_length=250, ins_rate=9.e-5, del_rate=1.1e-4,
                 insert_mean=500, insert_sd=100, insert_min=100, insert_max=None):
        """
        Initialise a read generator.
        :param method: The two library preparation methods are: 'meta3c' or 'hic'.
        :param enzyme: The employed restriction enzyme
        :param seed: a random seed - required to initial Art.
        :param random_state: additionally the random state object
        :param prefix: leading string for read names
        :param simple: True: do not simulate sequencing errors (faster), False: fully simulation sequencing
        :param machine_profile: ART Illumina error profile for a particular machine type. Default EmpMiSeq250
        :param read_length: desired read-length
        :param ins_rate: insertion rate (Art default: 9e-5)
        :param del_rate: deletion rate (Art default: 1.1e-4)
        :param insert_mean: mean sequencing insert length (must be > 0 and minimum)
        :param insert_sd: standard deviation of insert length  (must be < mean)
        :param insert_min: minimum allowable insert length (must be > 50)
        :param insert_max: maximum allowable insert length (must be > mean)
        """

        self.method = method
        self.cut_site = enzyme.ovhgseq * 2
        self.prefix = prefix
        self.seq_id_fmt = prefix + ':{seed}:{mode}:1:1:1:{idx} {r1r2}:Y:18:1'
        self.wgs_desc_fmt = 'WGS {repl.name}:{x1}..{x2}:{dir}'
        self._3c_desc_fmt = method.upper() + ' {repl1.name}:{x1} {repl2.name}:{x2}'

        assert insert_min > 50, 'Minimum allowable insert size is 50bp'
        assert insert_min < insert_mean, 'Minimum insert size must be less than expected mean'
        assert insert_mean > 0 and insert_sd > 0, 'Insert mean and stddev must be greater than 0'
        if insert_mean < insert_sd:
            print 'Warning: specified insert mean ({0}) less than stddev ({1})'.format(insert_mean, insert_sd)
        if insert_mean - insert_sd < insert_min:
            print 'Warning: specified insert mean ({0}) and stddev ({1}) will produce many inserts below ' \
                  'the minimum allowable insert length ({2})'.format(insert_mean, insert_sd, insert_min)

        self.insert_mean = insert_mean
        self.insert_sd = insert_sd
        self.insert_min = insert_min
        self.insert_max = insert_max
        self.too_short = 0
        self.too_long = 0

        self.seed = seed
        self.random_state = random_state
        self.uniform = random_state.uniform
        self.normal = random_state.normal
        self.randint = random_state.randint

        # initialise ART read simulator
        self.art = Art.Art(read_length, Art.EmpDist.create(machine_profile), ins_rate, del_rate, seed=self.seed)

        # set the method used to generate reads
        if simple:
            self.next_pair = self.art.next_pair_simple_seq
        else:
            self.next_pair = self.art.next_pair_indel_seq

        try:
            method_switcher = {
                'hic': self._part_joiner_sitedup,
                'meta3c': self._part_joiner_simple,
                'dnase': self._part_joiner_simple
            }
            self._part_joiner = method_switcher[self.method.lower()]
        except Exception:
            raise Sim3CException('unknown library preparation method ({0}) Either: \'3c\' or \'hic\']'.format(method))

    def get_report(self):
        """
        Prepare a small report on the number of inserts which were constrained by either the maximum
        or minimum length settings.
        """
        msg = 'Constrained inserts: '
        msg += '[len < {min}] = {min_count}'.format(min=self.insert_min, min_count=self.too_short)
        if self.insert_max:
            msg += ', [len > {max}] = {max_count}'.format(max=self.insert_max, max_count=self.too_long)
        else:
            msg += ', [no max limit]'
        return msg

    def _part_joiner_simple(self, a, b):
        """
        Join two fragments end to end without any site duplication. The new
        fragment will begin at a_0 and end at b_max. Used when no end fills are
        applied to the library prep (Eg. meta3C)
        Altogether [a_0..a_max] + [b_0..b_max]
        :param a: fragment a
        :param b: fragment b
        :return: a + b
        """
        return a + b

    def _part_joiner_sitedup(self, a, b):
        """
        Join two fragments end to end where the cut-site is duplicated. The new
        fragment will begin at a_0 end at b_max. Used when end-fill is applied
        before further steps in library preparation (Eg. HiC)
        Altogether [a_0..a_max] + [cs_0..cs_max] + [b_0..b_max]
        :param a: fragment a
        :param b: fragment b
        :return: a + b
        """
        return a + self.cut_site + b

    def draw_insert(self):
        """
        Draw an insert length, midpoint and direction. Since this is normally distributed
        even with reasonable mean and stddev, it is possible to return lengths
        less than zero. These cases are handled by redrawing until we get a length
        greater than the specified minimum.
        length.
        :return: length, midpoint and direction tuple. Eg. (100, 56, True)
        """
        length = int(self.normal(self.insert_mean, self.insert_sd))
        if length < self.insert_min:
            self.too_short += 1
            length = self.insert_min
        elif self.insert_max and length > self.insert_max:
            self.too_long += 1
            length = self.insert_max
        midpoint = self.randint(0, length)
        is_fwd = self.uniform() < 0.5
        return length, midpoint, is_fwd

    def make_wgs_readpair(self, repl, x1, ins_len, is_fwd):
        """
        Create a fwd/rev read-pair simulating WGS sequencing.
        :param repl: replicon from which to extract this read-pair
        :param x1: the initial coordinate along the replicon.
        :param ins_len: the insert length
        :param is_fwd: will the insert be off the fwd strand
        :return: a read-pair dict
        """
        frag = repl.subseq(x1, ins_len, is_fwd)
        pair = self.next_pair(str(frag.seq))
        pair['mode'] = 'WGS'
        pair['desc'] = self.wgs_desc_fmt.format(repl=repl, x1=x1, x2=x1+ins_len, dir='F' if is_fwd else 'R')
        return pair

    def make_ligation_readpair(self, repl1, x1, repl2, x2, ins_len, ins_junc):
        """
        Create a fwd/rev read-pair simulating a ligation product (Eg. the outcome of
        HiC or meta3C library prep). As repl1 and repl2 can be the same, these ligation
        products can be inter-rep, intra-rep or spurious products.
        :param repl1: the first replicon
        :param x1:  the location along repl1
        :param repl2: the second replicon
        :param x2: the location along repl2
        :param ins_len: insert length
        :param ins_junc: junction point on insert
        :return: a read-pair dict
        """
        part_a = repl1.subseq(x1 - ins_junc, ins_junc)
        part_b = repl2.subseq(x2, ins_len - ins_junc)

        pair = self.next_pair(str(self._part_joiner(part_a, part_b).seq))
        pair['mode'] = '3C'
        pair['desc'] = self._3c_desc_fmt.format(repl1=repl1, x1=x1, repl2=repl2, x2=x2)
        return pair

    def write_readpair(self, h_out, pair, index, fmt='fastq'):
        """
        Write a read-pair object to a stream.
        :param h_out: the output stream
        :param pair: the read-pair to write
        :param index: a unique identifier for the read-pair. (Eg. an integer)
        :param fmt: the output file format
        """

        # create Bio.Seq objects for read1 (fwd) and read2 (rev)
        read1 = pair['fwd'].read_record(self.seq_id_fmt.format(seed=self.seed, mode=pair['mode'], idx=index, r1r2=1),
                                        desc=pair['desc'])
        read2 = pair['rev'].read_record(self.seq_id_fmt.format(seed=self.seed, mode=pair['mode'], idx=index, r1r2=2),
                                        desc=pair['desc'])
        # write to interleaved file
        SeqIO.write([read1, read2], h_out, fmt)


class SequencingStrategy:
    """
    A SequencingStrategy represents the whole experiment. This includes the reference data from which
    WGS and ligation products are generated, and the simulation of Illumina sequencing reads.

    Experiments are reproducible by supplying the same seed value.
    """

    Strategy = namedtuple('Strategy', 'method run')

    def __init__(self, seed, prof_filename, seq_filename, enz_name, number_pairs,
                 method, read_length, prefix, machine_profile,
                 insert_mean=400, insert_sd=50, insert_min=50, insert_max=None,
                 anti_rate=0.25, spurious_rate=0.02, trans_rate=0.1,
                 efficiency=0.02,
                 ins_rate=9.e-5, del_rate=1.1e-4,
                 create_cids=True, simple_reads=True):
        """
        Initialise a SequencingStrategy.

        :param seed: the random seed for all subsequent calls to numpy.random methods.
        :param prof_filename: the abundance profile for the community
        :param seq_filename: the matching sequence of replicon sequences in Fasta format
        :param enz_name: the restriction enzyme name (case sensitive)
        :param number_pairs: the number of read-pairs to generate
        :param method: the library preparation method (Either: 3c or hic)
        :param read_length: the length of reads
        :param prefix: read-names begin with this string
        :param machine_profile: ART Illumina sequencing machine profile
        :param insert_mean: mean insert length
        :param insert_sd: stddev insert length
        :param insert_min: minimum allowable insert length (must be > 50)
        :param insert_max: maximum allowable insert length (must be > mean)
        :param anti_rate: rate of anti-diagonal interactions
        :param spurious_rate: rate of spurious ligation products
        :param trans_rate: rate of inter-replicon (trans) ligation product
        :param efficiency: for meta3c simulation, the efficiency of ligation product generation.
        :param ins_rate: rate of sequencing insert errors
        :param del_rate: rate of sequencing deletion errors
        :param create_cids: simulate 3D structure, chromosomal interacting domains (CID)
        :param simple_reads: True: sequencing reads do not simulate error (faster), False: full simulation of sequencing
        """
        self.seed = seed
        self.prof_filename = prof_filename
        self.seq_filename = seq_filename
        self.enz_name = enz_name
        self.number_pairs = number_pairs
        self.simple_reads = simple_reads
        self.method = method
        self.read_length = read_length
        self.insert_min = insert_min
        self.insert_max = insert_max
        self.efficiency = efficiency

        # initialise the random state for the simulation
        self.random_state = np.random.RandomState(seed)

        self.enzyme = get_enzyme_instance(enz_name)
        self.profile = abn.read_profile(prof_filename, True)

        # initialise the community for the reference data
        self.community = Community(seq_filename, self.profile, self.enzyme, self.random_state, anti_rate=anti_rate,
                                   spurious_rate=spurious_rate, trans_rate=trans_rate,
                                   create_cids=create_cids)

        # preparate the read simulator for output
        self.read_generator = ReadGenerator(method, self.enzyme, seed, self.random_state,
                                            prefix=prefix, simple=simple_reads, machine_profile=machine_profile,
                                            read_length=read_length, insert_mean=insert_mean,
                                            insert_sd=insert_sd, insert_min=insert_min, insert_max=insert_max,
                                            del_rate=del_rate, ins_rate=ins_rate)

        # the method determines the strategy governing the creation of
        # ligation products and WGS reads.
        try:
            strategy_switcher = {
                'hic': self._simulate_hic,
                'meta3c': self._simulate_meta3c,
                'dnase': self._simulate_dnase
            }
            self._selected_strat = self.Strategy(method, strategy_switcher[method.lower()])
        except Exception:
            raise Sim3CException('unknown library preparation method ({0}) Either: \'3c\' or \'hic\']'.format(method))

    def run(self, ostream):
        """
        Add some pre and post detail to the selected strategy.
        :param ostream: the output stream for reads
        """
        print 'Starting sequencing simulation'
        print 'Library method: {0}'.format(self._selected_strat.method)
        print 'Progress:'
        info = self._selected_strat.run(ostream)
        print 'Finished simulation'
        print 'Run Report:'
        print 'Read counts: WGS reads = {wgs_count}, Ligation products = {lig_count}'.format(**info)
        print '{0}'.format(self.read_generator.get_report())

    def _simulate_meta3c(self, ostream):
        """
        A strategy to simulate the sequencing of a 3C (meta3C) library.

        The most significant differentiator between 3C and HiC is that no biotin pulldown
        is employed for 3C sequencing experiments, thus only a small fraction of reads
        comprise ligation products, with the majority being conventional WGS reads.
        :param ostream: the output stream for reads
        """

        comm = self.community
        uniform = self.random_state.uniform
        read_gen = self.read_generator
        efficiency = self.efficiency

        n_wgs = 0
        n_3c = 0

        for n in tqdm.tqdm(xrange(1, self.number_pairs+1)):

            # pick an replicon, position and insert size
            r1, x1 = comm.draw_any_by_extent()
            ins_len, midpoint, is_fwd = read_gen.draw_insert()

            if uniform() < efficiency and r1.covers_site(x1, midpoint):

                n_3c += 1

                # move x1 to the nearest actual site
                x1 = r1.sites.find_nn(x1)

                # is it spurious ligation
                if comm.is_spurious():

                    r2, x2 = comm.draw_any_by_site()

                # is it an inter-replicon (trans) ligation
                elif r1.parent_cell.is_trans():

                    r2 = r1.parent_cell.draw_other_replicon_by_sites(r1.name)
                    x2 = r2.draw_any_site()

                # otherwise an intra-replicon (cis) ligation
                else:
                    # find the first intervening site between
                    # random constrained position x2 and site x1.
                    x2 = r1.draw_constrained_site(x1)
                    x2 = r1.sites.find_first(x2, x1)
                    r2 = r1

                # randomly permute source/destination
                if uniform() < 0.5:
                    x1, x2 = x2, x1
                    r1, r2 = r2, r1

                pair = read_gen.make_ligation_readpair(r1, x1, r2, x2, ins_len, midpoint)

            # otherwise WGS
            else:
                n_wgs += 1
                # take the already drawn coordinates
                pair = read_gen.make_wgs_readpair(r1, x1, ins_len, is_fwd)

            read_gen.write_readpair(ostream, pair, n)

        assert self.number_pairs - n_wgs == n_3c, 'Error: WGS and 3C pairs did not sum to ' \
                                                  '{0} was did not add'.format(self.number_pairs)

        return {'wgs_count': n_wgs, 'lig_count': self.number_pairs - n_wgs}

    def _simulate_hic(self, ostream):
        """
        A strategy to simulate the sequencing of a HiC library.
        :param ostream: the output stream for reads
        """

        comm = self.community
        uniform = self.random_state.uniform
        read_gen = self.read_generator
        efficiency = self.efficiency

        n_wgs = 0
        n_3c = 0

        for n in tqdm.tqdm(xrange(1, self.number_pairs+1)):

            ins_len, midpoint, is_fwd = read_gen.draw_insert()

            # is HIC pair?
            if uniform() <= efficiency:

                n_3c += 1

                # draw the first replicon and site
                r1, x1 = comm.draw_any_by_site()

                # is it spurious ligation
                if comm.is_spurious():

                    r2, x2 = comm.draw_any_by_site()

                # is it an inter-replicon (trans) ligation
                elif r1.parent_cell.is_trans():

                    r2 = r1.parent_cell.draw_other_replicon_by_sites(r1.name)
                    x2 = r2.draw_any_site()

                # otherwise an intra-replicon (cis) ligation
                else:

                    x2 = r1.draw_constrained_site(x1)
                    r2 = r1

                # randomly permute source/destination
                if uniform() < 0.5:
                    x1, x2 = x2, x1
                    r1, r2 = r2, r1

                # with coordinates, make hic read-pair
                pair = read_gen.make_ligation_readpair(r1, x1, r2, x2, ins_len, midpoint)

            # otherwise WGS
            else:
                n_wgs += 1
                r1, x1 = comm.draw_any_by_extent()
                pair = read_gen.make_wgs_readpair(r1, x1, ins_len, is_fwd)

            read_gen.write_readpair(ostream, pair, n)

        assert self.number_pairs - n_wgs == n_3c, 'Error: WGS and 3C pairs did not sum to ' \
                                                  '{0} was did not add'.format(self.number_pairs)

        return {'wgs_count': n_wgs, 'lig_count': n_3c}

    def _simulate_dnase(self, ostream):
        """
        A strategy to simulate the sequencing of a HiC library without a restriction enzyme.
        :param ostream: the output stream for reads
        """
        comm = self.community
        uniform = self.random_state.uniform
        read_gen = self.read_generator
        efficiency = self.efficiency

        n_wgs = 0
        n_3c = 0

        for n in tqdm.tqdm(xrange(1, self.number_pairs+1)):

            ins_len, midpoint, is_fwd = read_gen.draw_insert()

            # is PLP?
            if uniform() >= efficiency:

                n_3c += 1

                # draw the first replicon and site
                r1, x1 = comm.draw_any_by_extent()

                # is it spurious ligation
                if comm.is_spurious():

                    r2, x2 = comm.draw_any_by_extent()

                # is it an inter-replicon (trans) ligation
                elif r1.parent_cell.is_trans():

                    r2 = r1.parent_cell.draw_other_replicon_by_extents(r1.name)
                    x2 = r2.draw_any_site()

                # otherwise an intra-replicon (cis) ligation
                else:

                    x2 = r1.draw_constrained_site(x1)
                    r2 = r1

                # randomly permute source/destination
                if uniform() < 0.5:
                    x1, x2 = x2, x1
                    r1, r2 = r2, r1

                # with coordinates, make hic read-pair
                pair = read_gen.make_ligation_readpair(r1, x1, r2, x2, ins_len, midpoint)

            # otherwise WGS
            else:
                n_wgs += 1
                r1, x1 = comm.draw_any_by_extent()
                pair = read_gen.make_wgs_readpair(r1, x1, ins_len, is_fwd)

            read_gen.write_readpair(ostream, pair, n)

        assert self.number_pairs - n_wgs == n_3c, 'Error: WGS and 3C pairs did not sum to ' \
                                                  '{0} was did not add'.format(self.number_pairs)

        return {'wgs_count': n_wgs, 'lig_count': n_3c}


if __name__ == '__main__':
    import abundance
    import io_utils

    import argparse
    import sys
    import time
    import os

    #
    # Commandline interface
    #
    parser = argparse.ArgumentParser(description='Simulate HiC read pairs')

    parser.add_argument('-C', '--compress', choices=['gzip', 'bzip2'], default=None,
                        help='Compress output files')

    parser.add_argument('-r', '--seed', metavar='INT', type=int, default=int(time.time()),
                        help="Random seed for initialising number generator")
    parser.add_argument('-m', '--method', default='hic', choices=['hic', 'meta3c', 'dnase'],
                        help='Library preparation method [hic]')
    parser.add_argument('-e', '--enzyme', dest='enzyme_name', default='NlaIII',
                        help='Restriction enzyme (case-sensitive) [NlaIII]')

    parser.add_argument('-n', '--num-pairs', metavar='INT', type=int, required=True,
                        help='Number of read-pairs generate')
    parser.add_argument('-l', '--read-length', metavar='INT', type=int, required=True,
                        help='Length of reads from Hi-C fragments')
    parser.add_argument('--prefix', default='SIM3C', help='Prefix for read names [SIM3C]')
    parser.add_argument('--insert-mean', metavar='INT', type=int, default=400,
                        help='Mean insert size [400]')
    parser.add_argument('--insert-sd', metavar='INT', type=int, default=50,
                        help='Standard deviation of insert sizes [50]')
    parser.add_argument('--insert-min', metavar='INT', type=int, default=100,
                        help='Minimum allowed insert size [100]')
    parser.add_argument('--insert-max', metavar='INT', type=int, default=None,
                        help='Maximum allowed insert size [None]')

    parser.add_argument('--create-cids', default=False, action='store_true',
                        help='Simulate chromosome interacting domains')
    parser.add_argument('--efficiency', metavar='FLOAT', type=float,
                        help='HiC/Meta3C efficiency factor [hic: 0.5 or meta3c: 0.02]')
    parser.add_argument('--anti-rate', metavar='FLOAT', type=float, default=0.2,
                        help='Rate of anti-diagonal fragments [0.2]')
    parser.add_argument('--trans-rate', metavar='FLOAT', type=float, default=0.1,
                        help='Rate of inter-replicon (trans) fragment formation [0.1]')
    parser.add_argument('--spurious-rate', metavar='FLOAT', type=float, default=0.02,
                        help='Rate of spurious fragment formation [0.02]')

    parser.add_argument('-P', '--profile', dest='profile_in', metavar='FILE',
                        help='Community abundance profile')
    parser.add_argument('--profile-name', metavar='FILE', default='profile.tsv',
                        help='Output file name for a procedural community profile', required=False)

    parser.add_argument('--dist', metavar='DISTNAME', choices=['equal', 'uniform', 'lognormal'],
                        help='Abundance profile distribution choices: equal, uniform, lognormal')
    parser.add_argument('--lognorm-mu', metavar='FLOAT', type=float, default='1', required=False,
                        help='Log-normal relative abundance mu parameter')
    parser.add_argument('--lognorm-sigma', metavar='FLOAT', type=float, default='1', required=False,
                        help='Log-normal relative abundance sigma parameter')

    parser.add_argument('--simple-reads', default=False, action='store_true', help='Do not simulate sequencing errors')
    parser.add_argument('--machine-profile', help='An ART sequencing machine profile [EmpMiSeq250]',
                        default='EmpMiSeq250', choices=Art.ILLUMINA_PROFILES.keys())
    parser.add_argument('--ins-rate', type=float, default=9.e-5, help='Insert rate [9e-5]')
    parser.add_argument('--del-rate', type=float, default=1.1e-4, help='Deletion rate [1.1e-4]')

    parser.add_argument(dest='genome_seq', metavar='FASTA',
                        help='Genome sequences for the community')
    parser.add_argument(dest='output_file', metavar='OUTPUT',
                        help='Output Hi-C reads file')
    args = parser.parse_args()

    try:

        if 'community_table' in args and args.dist:
            raise RuntimeError('Cannot define abundance both explicitly as a table (-t) and a distribution (--dist).')

        if args.method == 'dnase' and args.enzyme_name:
            raise RuntimeError('The dnase method does not accept an enyzme specification.')

        #
        # Prepare community abundance profile, either procedurally or from a file
        #
        #   Note: currently, all sequences for a single taxon are
        #   treated equally.
        #
        if not args.profile_in and not args.dist:
            print 'An abundance profile must be supplied either as a file or procedurally'
            sys.exit(1)

        profile = None
        if args.dist:
            # generate a procedural profile.
            # the number of taxa is defined by number of sequences. i.e. monochromosomal organisms

            if os.path.basename(args.profile_name) != args.profile_name:
                print 'Arguments to profile-name should not contain path information'
                sys.exit(1)

            profile_path = os.path.join(os.path.dirname(args.output_file), args.profile_name)
            if os.path.exists(profile_path):
                print 'A previous procedural abundance profile already exists'
                print 'Please delete or move away: {0}'.format(profile_path)
                sys.exit(1)

            seq_names = None
            seq_index = SeqIO.index(args.genome_seq, 'fasta')
            try:
                seq_names = list(seq_index)
            finally:
                seq_index.close()

            profile = abundance.generate_profile(args.seed, seq_names, mode=args.dist,
                                                 lognorm_mu=args.lognorm_mu, lognorm_sigma=args.lognorm_sigma)

            # present result to console
            profile.write_table(sys.stdout)

            # save result to file
            with open(profile_path, 'w') as h_out:
                profile.write_table(h_out)

            # generated profile will be used downstream
            args.profile_in = profile_path

        if not args.efficiency:
            if args.method == 'hic':
                args.efficiency = 0.5
            elif args.method == 'meta3c':
                args.efficiency = 0.02

        # list of CLI arguments to pass as parameters to the simulation
        kw_names = ['prefix', 'machine_profile', 'insert_mean', 'insert_sd', 'insert_min', 'insert_max',
                    'anti_rate', 'spurious_rate', 'trans_rate',
                    'efficiency',
                    'ins_rate', 'del_rate',
                    'create_cids', 'simple_reads']

        # extract these parameters from the parsed arguments
        kw_args = {k: v for k, v in vars(args).items() if k in kw_names}

        # initialise a sequencing strategy for this community
        # and the given experimental parameters
        strategy = SequencingStrategy(args.seed, args.profile_in, args.genome_seq, args.enzyme_name,
                                      args.num_pairs, args.method, args.read_length, **kw_args)

        # Run the simulation
        with io_utils.open_output(args.output_file, mode='w', compress=args.compress) as out_stream:
            strategy.run(out_stream)

    except Exception as ex:
        print 'Error: {0}'.format(ex)
