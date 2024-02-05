import itertools
import logging
import re
import numpy as np
import numba as nb

from collections import namedtuple
from Bio.Restriction import Restriction
from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict

from .exceptions import *
import sim3C.random as random

logger = logging.getLogger(__name__)

# Information pertaining to digestion / private class for readability
digest_info = namedtuple('digest_info', ('enzyme', 'end5', 'end3', 'overhang'))
# Information pertaining to a Hi-C ligation junction
ligation_info = namedtuple('ligation_info',
                           ('enz5p', 'enz3p', 'junction', 'vestigial', 'junc_len',
                            'vest_len', 'pattern', 'cross'))


def get_enzyme_instance_ipython(enz_name):
    """
    An alternative method to fetch an instance of a given restriction enzyme by its
    name using a work-around which avoids exception with getattr() in iPython.

    Ackn: StackOverflow: user xbello.
    See: http://stackoverflow.com/questions/20381912/type-object-restrictiontype-has-no-attribute-size

    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    r_type_names = [rt for tid, (rt, rn) in typedict.items() if enz_name in rn][0]
    r_clz = tuple(getattr(Restriction, rt) for rt in r_type_names)
    return Restriction.AbstractCut(enz_name, r_clz, rest_dict[enz_name])


def get_enzyme_instance(enz_name):
    """
    Fetch an instance of a given restriction enzyme by its name.
    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    return getattr(Restriction, enz_name)


def enzyme_ends(enzyme):
    """
    Determine the 5` and 3` ends of a restriction endonuclease. Here, "ends"
    are the parts of the recognition site not involved in overhang.
    :param enzyme: Biopython ins
    :return: a pair of strings containing the 5' and 3' ends
    """
    end5, end3 = '', ''
    ovhg_size = abs(enzyme.ovhg)
    if ovhg_size > 0 and ovhg_size != enzyme.size:
        a = abs(enzyme.fst5)
        if a > enzyme.size // 2:
            a = enzyme.size - a
        end5, end3 = enzyme.site[:a], enzyme.site[-a:]
    return end5.upper(), end3.upper()


def vestigial_site(enz, junc):
    """
    Determine the part of the 5-prime end cut-site that will remain when a ligation junction is
    created. Depending on the enzymes involved, this may be the entire cut-site or a smaller portion
    and begins from the left.
    """
    i = 0
    while i < enz.size and (enz.site[i] == junc[i] or enz.site[i] == 'N'):
        i += 1
    return str(junc[:i]).upper()


def get_ligation_info(enzyme_a, enzyme_b=None):
    """
    For the enzyme cocktail, generate the set of possible ligation junctions.
    :return: a dictionary of ligation_info objects
    """
    end5, end3 = enzyme_ends(enzyme_a)
    enz_list = [digest_info(enzyme_a, end5, end3, enzyme_a.ovhgseq.upper())]

    if enzyme_b is not None:
        end5, end3 = enzyme_ends(enzyme_b)
        enz_list.append(digest_info(enzyme_b, end5, end3, enzyme_b.ovhgseq.upper()))

    _junctions = {}
    for n, (i, j) in enumerate(itertools.product(range(len(enz_list)), repeat=2), start=1):
        a, b = enz_list[i], enz_list[j]
        # the actual sequence
        _junc_seq = '{}{}{}{}'.format(a.end5, a.overhang, b.overhang, b.end3)
        if _junc_seq not in _junctions:
            _vest_seq = vestigial_site(a.enzyme, _junc_seq)
            _lij = ligation_info(str(a.enzyme), str(b.enzyme),
                                 _junc_seq, _vest_seq,
                                 len(_junc_seq), len(_vest_seq),
                                 re.compile(_junc_seq.replace('N', '[ACGT]')),
                                 str(a.enzyme) == str(b.enzyme))

            _junctions[(i, j)] = _lij

    logger.info('Predicted ligation junctions in the case of Hi-C reactions are: '
                f'{", ".join([_l.junction for _l in _junctions.values()])}')

    return _junctions


@nb.njit('i8[:](i8[:,:], i8, i8)')
def _fast_find_nn_circular(sites, max_length, pos):
    """
    Find the nearest cut-site relative to supplied position on a circular
    chromosome.
    :param pos: the position on the chromosome
    :return: nearest site
    """
    ix = np.searchsorted(sites[:, 0], pos)
    # modulo so as we only return values within the range [0..maxlen]
    # this handles the edge case of sites crossing beginning or end.
    if pos - sites[ix - 1, 0] <= sites[ix, 0] - pos:
        rv = sites[ix - 1].copy()
        rv[0] %= max_length
    else:
        rv = sites[ix].copy()
        rv[0] %= max_length
    return rv


@nb.njit('i8[:](i8[:,:], i8)')
def _fast_find_nn_linear(sites, pos):
    """
    Find the nearest cut-site relative to the supplied position on a linear
    chromosome or sequence fragment.
    :param pos: the position on the chromosome
    :return: nearest site
    """
    ix = np.searchsorted(sites[:, 0], pos)
    # first or last site was closest
    if ix == 0:
        return sites[0]
    elif ix == sites.shape[0]:
        return sites[-1]
    # pick the closest of nearest neighbours
    if pos - sites[ix - 1, 0] <= sites[ix, 0] - pos:
        return sites[ix - 1]
    else:
        return sites[ix]


@nb.njit('i8[:](i8[:,:], i8, i8)')
def _fast_find_first(sites, pos, origin_site):
    """
    Beginning from pos and searching toward the origin, return the first encountered cut-site.
    This may turn out to be the origin if there is no other intervening site.
    :param pos: a position on the chromosome
    :param origin_site: a valid cut-site used as the relative origin
    :return: the first site encountered between pos and origin_site
    """
    if pos > origin_site:
        # the position begins after the origin
        ix = np.searchsorted(sites[:, 0], pos, side='right')
        if ix == 0:
            return sites[ix]
        else:
            return sites[ix-1]
    else:
        # position begins before the origin site
        ix = np.searchsorted(sites[:, 0], pos, side='left')
        return sites[ix]


class CutSites(object):
    """
    The cut-sites for a given enzyme on a given template.
    Note: locations are 0-based
    """

    def __init__(self, template_seq, enzyme_a, enzyme_b=None, linear=False):
        """
        Initialize the cut-sites of an enzyme on a template sequence.

        Linearity affects both the initial search for recognition sites and
        when requesting the nearest site to a given genomic location.

        :param template_seq: the template sequence to digest (Bio.Seq object)
        :param enzyme_a: the first enzyme
        :param enzyme_b: a second enzyme (optional)
        :param linear: treat sequence as linear
        """
        self.enzyme_a = enzyme_a
        self.enzyme_b = enzyme_b
        self.max_length = len(template_seq) - 1

        # find sites, converting from 1-based.
        self.sites = self._find_sites(template_seq, linear)
        if self.sites.shape[0] == 0:
            enz_names = ','.join([str(enz) for enz in [enzyme_a, enzyme_b] if enz is not None])
            raise NoCutSitesException(f'No cut-sites found using: {enz_names}')
        self.size = self.sites.shape[0]

        # method setup
        if linear:
            self.find_nn = self._find_nn_linear
            self.covers = self._covers_site_linear
        else:
            self.find_nn = self._find_nn_circular
            self.covers = self._covers_site_circular
            self.shouldered = self._add_shoulders()

    def _find_sites(self, template_seq, linear):
        sites = []
        for label, enzyme in enumerate([self.enzyme_a, self.enzyme_b]):
            if enzyme is None:
                continue
            cs = np.array(enzyme.search(template_seq, linear), dtype='i8') - 1
            cs = np.vstack([cs, np.zeros_like(cs)]).T
            cs[:, 1] = label
            sites.append(cs)

        # combine all sites into a single array
        sites = np.concatenate(sites)
        # reorder the array by coordinate
        sites = sites[np.argsort(sites[:, 0]),]
        return sites

    def _add_shoulders(self):
        """
        Add shoulders to the already determined list of circular chr sites. This allows
        finding positions without logic for edge cases (ix=0 or -1)
        """
        before_first = self.sites[[-1]]
        before_first[0, 0] -= self.max_length
        after_last = self.sites[[0]]
        after_last[0, 0] += self.max_length
        return np.vstack(([before_first, self.sites, after_last]))

    def random_site(self):
        """
        Select a uniformly random site
        :return: a random site
        """
        return self.sites[random.pcg_random.integer(self.size)]

    def _covers_site_linear(self, x1: np.ndarray, length):
        """
        Test whether a specified position and length contains any cut-sites, assuming
        a linear molecule.
        :param x1: starting position
        :param length: length
        :return: True the coordinates contain at least one cut-site
        """
        return ((self.sites[:, 0] > x1) & (self.sites[:, 0] <= x1 + length)).any()

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
            ret = ((x2 % self.max_length) >= self.sites[:, 0]) | (x1 <= self.sites[:, 0])
        else:
            ret = (x1 <= self.sites[:, 0]) & (x2 >= self.sites[:, 0])
        return ret.any()

    def _find_nn_circular(self, pos):
        """
        Find the nearest cut-site relative to supplied position on a circular
        chromosome.
        :param pos: the position on the chromosome
        :return: nearest site
        """
        return _fast_find_nn_circular(self.shouldered, self.max_length, pos)

    def _find_nn_linear(self, pos):
        """
        Find the nearest cut-site relative to the supplied position on a linear
        chromosome or sequence fragment.
        :param pos: the position on the chromosome
        :return: nearest site
        """
        return _fast_find_nn_linear(self.sites, pos)

    def find_first(self, pos, origin_site):
        """
        Beginning from pos and searching toward the origin, return the first encountered cut-site.
        This may turn out to be the origin if there is no other intervening site.
        :param pos: a position on the chromosome
        :param origin_site: a valid cut-site used as the relative origin
        :return: the first site encountered between pos and origin_site
        """
        return _fast_find_first(self.sites, pos, origin_site)


class AllSites(object):

    def __init__(self, size):
        """
        Initialise an instance of AllSites. This represents a completely unconstrained
        model, where any base-position is equally accessible.

        :param size: the maximum length of a sequence
        """
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
        return random.pcg_random.integer(self.size)

    @staticmethod
    def find_nn(pos):
        """
        As every site is accessible this is effectively a dummy function returning the input.
        :param pos: the position to find a NN for
        :return: pos, unchanged.
        """
        return pos
