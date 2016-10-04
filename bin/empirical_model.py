import numpy as np
from intervaltree import Interval, IntervalTree


def cdf_geom_unif(x, length, **kwargs):
    """
    CDF as a linear combination of the geometric and uniform CDFs, covering
    values over the range 0..length. Must supply a geometric shape coeff as
    a keyword arg.
    :param x: 0..1
    :param length: maximum value
    :param kwargs: 'shape' geometric distribution coeff.
    :return: float from 0..length
    """
    return 0.5 * (1.0 - (1.0 - kwargs['shape']) ** x + 1.0/length * x)


class EmpiricalDistribution:
    """
    Defining an empirical distribution, we can then use it to draw random numbers. A user need
    only supply a definition for the CDF for the desired distribution. This will in turn by sampled
    and used as a lookup table when generating random values. For quickly varying distributions,
    a small number of bins may be insufficient. The values are linearly interpolated between the
    nearest two bins.
    """

    def __init__(self, random_state, length, bins, cdf, **coeffs):
        """
        :param random_state: random state from which to draw numbers. If None, then this will be initialized at
        :param shape: distribution shape parameter
        :param length: distribution will be defined over 0..length
        :param bins: number of bins to sample distribution
        :param cdf: cdf from which to generate sampled distribution
        runtime.
        """
        self.random_state = random_state
        self.bins = bins
        self.coeffs = coeffs
        self.cdf = cdf
        self.length = length
        self.xsample = np.linspace(0, length, bins, endpoint=True, dtype=np.float64)
        self.ysample = self.cdf(self.xsample, length, **self.coeffs)
        self.ysample /= self.ysample.max()

    def __add__(self, other):
        """
        Add two instances together, where they must match in terms of CDF, number of bins and coefficients.
        Largely, this involves adding their lengths together and recomputing.
        :param other: an instance of EmpiricalDistribution
        :return:
        """
        assert isinstance(other, EmpiricalDistribution)
        assert self.cdf == other.cdf
        assert self.bins == other.bins
        assert self.coeffs == other.coeffs
        return EmpiricalDistribution(self.length + other.length, self.bins, self.cdf, **self.coeffs)

    def rand(self):
        """
        Using the inverse CDF method, draw a random number for the distribution. This
        method looks up the nearest value of random value x in a sampled representation
        and then interpolates between bin edges. Edge case for the first and last bin
        is to merely use that bin.

        :return: random value following distribution
        """
        return np.interp(self.random_state.uniform(), self.ysample, self.xsample)


def _reducer_cid_data(acc, x):
    """
    IntervalTree merge_overlaps method uses a reducer function to accumulate the Interval data fields. If not
    included, merge_overlaps default behaviour is to throw-away this information.
    :param acc: the initial data
    :param x:  the next data
    :return: combined data
    """
    # probabilities are averaged, empirical distributions are summed.
    return {'prob': 0.5*(acc['prob'] + x['prob']), 'empdist': acc['empdist'] + x['empdist']}


def generate_random_cids(random_state, genome_length, genome_prob=2.0, genome_bins=1000, genome_shape=8.0e-6,
                         min_cid_len=20000, max_cid_len=250000, num_cid=10, cid_bins=100, cid_shape=6.0e-6,
                         merge_overlaps=False):
    """
    Generate a random set of CID intervals for a given genome size. Default values have been set for most
    parameters.

    :param random_state:
    :param genome_length:
    :param genome_prob:
    :param genome_bins:
    :param genome_shape:
    :param min_cid_len:
    :param max_cid_len:
    :param num_cid:
    :param cid_bins:
    :param cid_shape:
    :param merge_overlaps:
    :return:
    """
    # create the list of CID intervals as (cid_begin, cid_length) pairs.
    data = map(lambda _: (random_state.randint(genome_length),
                          random_state.randint(min_cid_len, max_cid_len)), range(num_cid))

    # initialise the interval-tree, where we associate a random probability for each CID interval
    # where we assume their distributions of the same basic type as the full genome.
    cid_tree = IntervalTree(Interval(x, x+y,
                                     {'prob': random_state.uniform(),
                                      'empdist': EmpiricalDistribution(random_state, y, cid_bins, cdf_geom_unif, shape=cid_shape)}
                                     ) for x, y in data if y < genome_length)

    # if requested, simplify the CID tree into non-overlapping regions
    if merge_overlaps:
        cid_tree.merge_overlaps(data_reducer=_reducer_cid_data)

    # add the interval governing the whole genome -- call it the genome-wide CID. mother of all CID? heh
    cid_tree[0:genome_length] = {'prob': genome_prob,
                                 'empdist': EmpiricalDistribution(random_state, genome_length, genome_bins,
                                                                  cdf_geom_unif, shape=genome_shape)}

    return cid_tree


def cids_to_blocks(cid_tree):
    """
    From CID IntervalTree, create a new IntervalTree where now a interval represents a region of constant
    effect. That is, each resulting interval describes range where a constant set of CIDs are involved. Blocks
    do not overlap, but are perfectly adjacent. (zero spacing)
    :param cid_tree: an intervaltree of representing genomic intervals with CIDs
    :return: a blocked CID intervaltree.
    """

    # Get all the begin and end points in ascending order.
    # As they mark where a CID either begins are ends, each therefore
    # marks the end of one block and the beginning of another.
    x = []
    for inv in cid_tree:
        x.append(inv.begin), x.append(inv.end)
    x.sort()

    # interate over the CID coords, making all the block intervals.
    block_tree = IntervalTree()
    for i in xrange(len(x)-1):
        ovl_invs = sorted(cid_tree[x[i]:x[i+1]])  # the CIDs involved in this range
        # normalize probs for the block.
        p = np.fromiter((inv.data['prob'] for inv in ovl_invs), dtype=float)
        p /= p.sum()
        # a block stores the normalized probabilities and originating CID intervals for quick lookup.
        block_tree[x[i]:x[i+1]] = {'prob_list': p, 'inv_list': ovl_invs}

    return block_tree
