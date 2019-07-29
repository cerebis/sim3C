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
import numpy as np

from intervaltree import Interval, IntervalTree

from .random import uniform, randint


def cdf_geom(x, shape):
    return 1. - (1. - shape) ** x


def cdf_geom_unif_ratio(x, length, alpha=0.333, **kwargs):
    """
    CDF as an adjustable linear combination of the geometric and uniform CDFs, covering
    values over the range 0..length. Must supply a geometric shape coeff as
    a keyword arg.
    :param x: 0..1
    :param length: maximum value
    :param alpha: proportional contribution to CDF from uniform distribution term. [0..1]
    :param kwargs: 'shape' geometric distribution coeff
    :return: float from 0..length
    """
    return (1.0 - alpha) * (1.0 - (1.0 - kwargs['shape']) ** x) + alpha / length * x


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
    return 0.5 * (1.0 - (1.0 - kwargs['shape']) ** x + 1.0 / length * x)


def pmf_geom_unif(x, length, **kwargs):
    """
    PMF for the linear combination of the geometric and uniform PMFs.
    :param x: the position at which to evaluate
    :param length: maxium value
    :param kwargs: 'shape' geometric distribution coeff
    :return: pmf value at x
    """
    return 0.5 * (kwargs['shape'] * (1. - kwargs['shape'])**x + 1./length)


class EmpiricalDistribution(object):
    """
    Defining an empirical distribution, we can then use it to draw random numbers. A user need
    only supply a definition for the CDF for the desired distribution. This will in turn by sampled
    and used as a lookup table when generating random values. For quickly varying distributions,
    a small number of bins may be insufficient. The values are linearly interpolated between the
    nearest two bins.
    """

    def __init__(self, length, bins, cdf, **coeffs):
        """
        Initialise an empirical distribution using the supplied CDF and for() the range [0..length].
        The CDF is normalized by 1 / max[CDF(x)].

        :param length: distribution will be defined over 0..length
        :param bins: number of bins to sample distribution
        :param cdf: cdf from which to generate sampled distribution
        :param coeffs: coefficients to cdf function
        runtime.
        """
        self.bins = bins
        self.coeffs = coeffs
        self.cdf = cdf
        self.length = length
        self.xsample = np.linspace(0, length, bins, endpoint=True, dtype=np.float64)
        self.ysample = self.cdf(self.xsample, length, **self.coeffs)
        self.ysample /= self.ysample.max()

    def eval_cdf(self, x):
        """
        Evaluate the CDF at the given position.
        :param x: position at which to evaulate
        :return: CDF value at x
        """
        assert x <= self.length, 'out of bounds {} > {}'.format(x, self.length)
        return self.cdf(x, self.length, **self.coeffs)

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
        return np.interp(uniform(), self.ysample, self.xsample)


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


def generate_random_cids(chr_length, chr_prob=0.5, chr_bins=1000, chr_shape=8.0e-6, cdf_alpha=0.333,
                         min_cid_len=20000, max_cid_len=250000, num_cid=10, cid_bins=100, cid_shape=6.0e-6,
                         merge_overlaps=False):
    """
    Generate a random set of CID intervals for a given genome size. Default values have been set for most
    parameters.

    :param chr_length: length of chromosome
    :param chr_prob: probability of chromosome selection
    :param chr_bins: number of sampling bins used in empirical distribution for chromosome
    :param chr_shape: geometric distribution shape parameter for chromosome
    :param min_cid_len: minimum possible size of CIDs
    :param max_cid_len: maximum possible size of CIDs
    :param num_cid: number of CIDs to generate for chromosome
    :param cid_bins: number of sampling bins used in empirical distribution of CID.
    :param cid_shape: geometric distribution shape parameter for CIDs
    :param cdf_alpha: mixture coefficient
    :param merge_overlaps: if true, overlapping CID are merged. This can result in CID exceeding max_cid_len.
    :return: an intervaltree representing range of effect of each CID, as well as the full chromosome.
    """

    # Although most parameters could have an assert, these two are possibly more useful
    assert chr_prob < 1.0, 'chromosomal probability must be less than 1'
    assert num_cid > 0, 'number of CIDs must be greater than 0'

    # Create the list of CID intervals as (cid_begin, cid_length) pairs.
    data = map(lambda _: (randint(chr_length),
                          randint(min_cid_len, max_cid_len)), range(num_cid))

    # Draw a set of independent probabilities and normalise these along with chr_prob to 1.
    # The closer chr_prob is to 1, the greater its precedence and the less role CIDs will
    # play in determining locations.
    cid_probs = np.array(uniform(size=num_cid))
    cid_probs *= (1.0 - chr_prob) / cid_probs.sum()
    cid_probs_iter = np.nditer(cid_probs)

    # Initialise the interval-tree, and associate a random interaction probability for each CID interval.
    # We assume their distributions are of the same form as the full genome.
    cid_tree = IntervalTree(Interval(x, x+y,
                                     {'prob': next(cid_probs_iter),
                                      'empdist': EmpiricalDistribution(y, cid_bins, cdf_geom_unif_ratio,
                                                                       shape=cid_shape, alpha=cdf_alpha)}
                                     ) for x, y in data if y < chr_length)

    # If requested, simplify the CID tree into non-overlapping regions
    if merge_overlaps:
        cid_tree.merge_overlaps(data_reducer=_reducer_cid_data)

    # Add the interval governing the whole genome -- call it the genome-wide CID. mother of all CID? heh
    cid_tree[0:chr_length] = {'prob': chr_prob,
                              'empdist': EmpiricalDistribution(chr_length, chr_bins, cdf_geom_unif_ratio,
                                                               shape=chr_shape, alpha=cdf_alpha)}

    return cid_tree


def _random_nested_intervals(result, inv, min_len, max_len, min_num, max_num, max_depth, depth=0):
    """
    Recursively divide an interval until we reach 'depth' levels of recursion. An interval is divided
    into a set of smaller intervals, constrained by a proportional min/max length and min/max number.
    :param result: the final list of lists of intervals.
    :param inv: the interval to divide
    :param min_len: the smallest proportional size of an interval [0..1]
    :param max_len: the largest proportional size of an interval [0..1]
    :param min_num: the smallest number of sub-intervals to create for this interval
    :param max_num: the largest number of sub-intervals to create for this interval
    :param max_depth: the maximum recursive depth
    :param depth: current depth (not intended for user)
    """
    if depth < max_depth:
        # draw a set of random points
        x = uniform(min_len * inv.length(), max_len * inv.length(),
                                 size=randint(min_num, max_num + 1))
        # from a sequence from these points and the begin/end of interval
        subseq = np.hstack([[inv.begin],
                            inv.begin + np.cumsum((x / x.sum() * inv.length()).astype(int))[:-1],
                            [inv.end]])
        # adjacent elements become the next level of intervals
        subinvs = [Interval(pi[0], pi[1], data={'depth': depth+1}) for pi in zip(subseq[:-1:1], subseq[1::1])]
        # keep this level
        result.append(subinvs)
        # continue to divide these new intervals
        for ii in subinvs:
            _random_nested_intervals(result, ii, min_len, max_len,
                                     min_num, max_num, max_depth, depth + 1)


def generate_nested_cids(chr_length, chr_prob, chr_bins, chr_shape, cid_bins, cid_shape,
                         cdf_alpha=0.333, min_len=0.05, max_len=0.2, min_num=5, max_num=5, recur_depth=2):
    """cdf
    Generate a set of nested CID intervals for the given genome size. This method better approximates the
    appearance of a bacterial contact map, conceptually that folded regions themselves fold again and so
    create intervals within themselves which interact. Default values to the method appear to produce
    reasonable outcomes.

    :param chr_length: the length of the chromsome
    :param chr_prob: the probability of a regular backbone interaction vs CID interaction
    :param chr_bins: the number of bins over which to sample the backbone CDF
    :param chr_shape: the backbone shape parameter
    :param cdf_alpha: proportion of CDF contributed by the uniform distribution term. [0..1]
    :param cid_bins: the number of bins over which to sample the CID CDFs
    :param cid_shape: the CID shape parameter
    :param min_len: the smallest proportional size of an interval [0..1]
    :param max_len: the largest proportional size of an interval [0..1]
    :param min_num: the smallest number of sub-intervals to create for this interval
    :param max_num: the largest number of sub-intervals to create for this interval
    :param recur_depth: the maximum recursive depth
    :return:
    """

    # recursively create a list of nested intervals
    cid_list = []
    top_inv = Interval(0, chr_length)
    _random_nested_intervals(cid_list, top_inv, min_len, max_len, min_num, max_num, recur_depth)

    # flatten returned list of intervals
    cid_list = [inv for level in cid_list for inv in level]

    # create random probs to assign to each interval, then normalise
    # so that: sum{P_cids} + P_backbone = 1.
    cid_probs = np.array(uniform(size=len(cid_list)))
    cid_probs *= (1.0 - chr_prob) / cid_probs.sum()

    # initialise the tree, where each interval now gets a
    # probability and empirical distribution associated with it.
    cid_probs_iter = np.nditer(cid_probs)
    cid_tree = IntervalTree()
    for inv in cid_list:
        # explicitly cast, avoiding a 1-element np array
        inv.data['prob'] = float(next(cid_probs_iter))
        inv.data['empdist'] = EmpiricalDistribution(inv.length(), cid_bins, cdf_geom_unif_ratio,
                                                    shape=cid_shape, alpha=cdf_alpha)
        cid_tree.add(inv)

    # Add the interval governing the whole genome
    data = {'depth': 0,
            'prob': chr_prob,
            'empdist': EmpiricalDistribution(chr_length, chr_bins, cdf_geom_unif_ratio,
                                             shape=chr_shape, alpha=cdf_alpha)}
    cid_tree.addi(0, chr_length, data=data)

    return cid_tree


def cids_to_blocks(cid_tree):
    """
    Using an IntervalTree as returned by generate_random_cids(), create a new IntervalTree where now the
    intervals represent regions (blocks) of homogeneous effect. That is, each resulting interval defines a
    region where a fixed set of CIDs are involved.

    Blocks, therefore, do not overlap but are instead perfectly adjacent (zero spacing). For a given block
    the independent CID probabilities are normalized to sum to 1, in preparation of selection by random draw.

    :param cid_tree: an IntervalTree representing CIDs and the chromsome.
    :return: an IntervalTree of the homogeneous blocks for this set of CID.
    """

    # Get all the begin and end points in ascending order.
    # As they mark where a CID either begins are ends, each therefore
    # marks the end of one block and the beginning of another.
    x = []
    for inv in cid_tree:
        x.append(inv.begin), x.append(inv.end)
    x = np.unique(x)

    # interate over the CID coords, making all the block intervals.
    block_tree = IntervalTree()
    for i in range(len(x)-1):
        ovl_invs = sorted(cid_tree[x[i]:x[i+1]])  # the CIDs involved in this range

        # normalize probs for the block.
        p = np.fromiter((inv.data['prob'] for inv in ovl_invs), dtype=np.float)
        p /= p.sum()

        # a block stores the normalized probabilities and originating CID intervals for quick lookup.
        block_tree.addi(x[i], x[i+1], {'prob_list': p, 'inv_list': ovl_invs})

    return block_tree
