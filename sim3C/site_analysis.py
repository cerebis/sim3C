import numpy as np

from collections import namedtuple
from Bio.Restriction import Restriction
from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict

from .exceptions import *
import sim3C.random as random

ligation_info_t = namedtuple('ligation_info', ('enzyme', 'junction', 'vestigial', 'junc_len', 'vest_len'))


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


def get_ligation_info(enz):
    """
    Determine the ligation junction and vestigial site for a given enzyme.
    :param enz: the ezyme to consider
    :return: ligation_info_t namedtuple
    """
    end5, end3 = enzyme_ends(enz)
    overhang = enz.ovhgseq.upper()
    _junc_seq = f'{end5}{overhang}{overhang}{end3}'
    _vest_seq = vestigial_site(enz, _junc_seq)
    return ligation_info_t(str(enz), _junc_seq, _vest_seq, len(_junc_seq), len(_vest_seq))


class CutSites(object):
    """
    The cut-sites for a given enzyme on a given template.
    Note: locations are 0-based
    """

    def __init__(self, enzyme, template_seq, linear=False):
        """
        Initialize the cut-sites of an enzyme on a template sequence.

        Linearity affects both the initial search for recognition sites and
        when requesting the nearest site to a given genomic location.

        :param enzyme: the restriction enzyme (Bio.Restriction RestrictionType object)
        :param template_seq: the template sequence to digest (Bio.Seq object)
        :param linear: treat sequence as linear
        """
        self.enzyme = enzyme
        self.max_length = len(template_seq) - 1

        # find sites, converting from 1-based.
        self.sites = np.array(enzyme.search(template_seq, linear)) - 1
        self.size = self.sites.shape[0]
        if self.size == 0:
            raise NoCutSitesException(str(enzyme))

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
        return self.sites[random.pcg_random.integer(self.size)]

    def _covers_site_linear(self, x1: np.ndarray, length):
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
            return cs[-1]
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
