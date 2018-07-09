import contextlib
import cPickle
import heapq
import itertools
import logging
import os
import subprocess
import sys
import tempfile
import uuid
from collections import OrderedDict, namedtuple, Mapping, Iterable
from functools import partial

import Bio.SeqIO as SeqIO
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import networkx as nx
import numpy as np
import pandas
import pysam
import scipy.sparse as sp
import scipy.stats.mstats as mstats
import seaborn
import tqdm
import yaml

from Bio.Restriction import Restriction
from difflib import SequenceMatcher
from numba import jit, vectorize, int64, int32, float64, void
from numpy import log, pi
from typing import Optional

import louvain_cluster
import io_utils
import ordering
import simple_sparse

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='map_ordering.log',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
console.setFormatter(formatter)
logger = logging.getLogger('map_ordering')
logger.addHandler(console)

MIN_FIELD = 2e-8
P2ALPHA = 0.122123774414444
P2LAMBDA = 13.675170758388262
P2MU = 13.973247315647466


SeqInfo = namedtuple('SeqInfo', ['offset', 'refid', 'name', 'length', 'sites'])


class UnknownEnzymeException(Exception):
    """All sequences were excluded during filtering"""
    def __init__(self, target, similar):
        super(UnknownEnzymeException, self).__init__(
            '{} is undefined, but its similar to: {}'.format(target, ', '.join(similar)))


class UnknownOrientationStateException(Exception):
    """All sequences were excluded during filtering"""
    def __init__(self, ori):
        super(UnknownOrientationStateException, self).__init__('unknown orientation state [{}].'.format(ori))


class NoneAcceptedException(Exception):
    """All sequences were excluded during filtering"""
    def __init__(self):
        super(NoneAcceptedException, self).__init__('all sequences were excluded')


class TooFewException(Exception):
    """LKH requires a minimum of three nodes"""
    def __init__(self):
        super(TooFewException, self).__init__('LKH requires a minimum of three sequences')


class NoReportException(Exception):
    """Clustering does not contain a report"""
    def __init__(self, clid):
        super(NoReportException, self).__init__('Cluster {} did not contain a report'.format(clid))


class ZeroLengthException(Exception):
    """Sequence of zero length"""
    def __init__(self, seq_name):
        super(ZeroLengthException, self).__init__('Sequence [{}] has zero length'.format(seq_name))


class ParsingError(Exception):
    """An error during input parsing"""
    def __init__(self, msg):
        super(ParsingError, self).__init__(msg)


def make_random_seed():
    """
    Provide a random seed value between 1 and 10 million.
    :return: integer random seed
    """
    return np.random.randint(1000000, 10000000)


def make_dir(path):
    """
    Convenience method for making directories with a standard logic.
    An exception is raised when the specified path exists and is not a directory.
    :param path: target path to create
    """
    if not os.path.exists(path):
        os.mkdir(path)
    elif not os.path.isdir(path):
        raise IOError('output path already exists and not a directory')


def ext_path(name):
    """
    Return path to named executable in the external bianries directory

    :param name: name of binary
    :return: absolute path
    """
    return os.path.join(sys.path[0], 'external', name)


@vectorize([float64(float64)])
def piecewise_3c(s):
    pr = MIN_FIELD
    if s < 500e3:
        # Pareto2
        pr = P2ALPHA / P2LAMBDA * (1 + (s - P2MU)/P2LAMBDA)**(- P2ALPHA - 1)
    return pr


@jit(float64(int32[:, :], float64[:, :]), nopython=True)
def poisson_lpmf2(ob, ex):
    """
    Entirely skips terms where no observational counts were recorded.

    :param ob: observed counts
    :param ex: expected counts
    :return: log likelihood
    """
    s = 0.0
    for i in xrange(ob.shape[0]):
        for j in xrange(ob.shape[1]):
            aij = ob[i, j]
            bij = ex[i, j]
            if aij == 0:
                continue
            s += aij * log(aij/bij) + bij - aij + 0.5 * log(2.0 * pi * aij)
    return -s


@jit(float64(int32[:, :], float64[:, :]), nopython=True)
def poisson_lpmf3(ob, ex):
    """
    All terms calculated.

    :param ob: observed counts
    :param ex: expected counts
    :return: log likelihood
    """
    s = 0.0
    for i in xrange(ob.shape[0]):
        for j in xrange(ob.shape[1]):
            aij = ob[i, j]
            bij = ex[i, j]
            if aij == 0:
                s += bij
            else:
                s += aij * log(aij/bij) + bij - aij + 0.5 * log(2.0 * pi * aij)
    return -s


def calc_likelihood(cm):
    """
    For a given order and ContactMap instance, calculate the log likelihood.

    :param cm: instance of ContactMap with matching identifiers
    :return: log likelihood
    """

    borders = cm.grouping.borders
    centers = cm.grouping.centers
    extent_map = cm.extent_map.tocsr().astype(np.int32)
    total_obs = cm.map_weight()

    lengths = cm.order.order['length']
    ori = cm.order.order['ori']

    log_l = 0.0
    for i, j in itertools.combinations(xrange(cm.total_seq), 2):

        # inter-contig separation defined by cumulative
        # intervening contig length.
        gap_length = cm.order.intervening(i, j)

        # contig lengths
        li = lengths[i]
        lj = lengths[j]

        # bin centers for each contig, relative to middle of each contig
        c_ik = centers[i]
        c_jl = centers[j]

        # orientation of sequences
        s_i = ori[i]
        s_j = ori[j]

        # all separations between bins, including the potential intervening distance L
        d_ij = gap_length + 0.5*(li + lj) + s_i * c_jl - s_j * c_ik.T

        # conversion to expected counts
        q_ij = total_obs * piecewise_3c(d_ij)

        # matrix element range which represents cross-terms between contig i and j
        i1, i2 = borders[i]
        j1, j2 = borders[j]

        # observed counts
        # for now this is converted to dense array as we need zeros
        n_ij = extent_map[i1:i2, j1:j2].todense()

        # log likelihood
        log_l += poisson_lpmf3(n_ij, q_ij)

    return log_l


# Cigar codes that are not permitted.
NOT_ALLOWED = {2, 3, 6, 8}


def good_match(cigartuples, min_match=None, match_start=True):
    """
    A confident mapped read.

    :param cigartuples: a read's CIGAR in tuple form
    :param min_match: the minimum number of matched base positions
    :param match_start: alignment must begin at read's first base position
    :return:
    """

    # restrict tuples to a subset of possible conditions
    if len(NOT_ALLOWED & set([t[0] for t in cigartuples])) != 0:
        return False

    # match the first N bases if requested
    elif match_start and cigartuples[0][0] != 0:
        return False

    # impose minimum number of matches
    elif min_match:
        n_matches = sum([t[1] for t in cigartuples if t[0] == 0])
        if n_matches < min_match:
            return False

    return True


def strong_match(mr, min_match=None, match_start=True, min_mapq=None):
    """
    Augment a good_match() by also checking for whether a read is secondary or
    supplementary.

    :param mr: mapped read to test
    :param min_match: the minimum number of matched base positions
    :param match_start: alignment must begin at read's first base position
    :param min_mapq: the minimum acceptable mapping quality. This can be mapper dependent. For BWA MEM
    any read which perfectly aligns in two distinct places, as mapq=0.
    """
    if min_mapq and mr.mapping_quality < min_mapq:
        return False

    elif mr.is_secondary or mr.is_supplementary:
        return False

    # reverse the tuple if this read is reversed, as we're testing from 0-base on read.
    cigtup = mr.cigartuples[::-1] if mr.is_reverse else mr.cigartuples

    return good_match(cigtup, min_match, match_start)


def mean_selector(name):
    """
    Basic Mean functions
    """
    def geometric_mean(x, y):
        return (x*y)**0.5

    def harmonic_mean(x, y):
        return 2*x*y/(x+y)

    def arithmetic_mean(x, y):
        return 0.5*(x+y)

    try:
        mean_switcher = {
            'geometric': geometric_mean,
            'harmonic': harmonic_mean,
            'arithmetic': arithmetic_mean
        }
        return mean_switcher[name]
    except KeyError:
        raise RuntimeError('unsupported mean type [{}]'.format(name))


class SiteCounter(object):

    def __init__(self, enzyme_names, tip_size=None, is_linear=True):
        """
        Simple class to count the total number of enzymatic cut sites for the given
        list if enzymes.

        :param enzyme_names: a list of enzyme names (proper case sensitive spelling a la NEB)
        :param tip_size: when using tip based counting, the size in bp
        :param is_linear: Treat sequence as linear.
        """
        if isinstance(enzyme_names, basestring):
            enzyme_names = [enzyme_names]
        assert isinstance(enzyme_names, Iterable) and \
               not isinstance(enzyme_names, basestring), 'enzyme_names must of a collection of names'
        self.enzymes = [SiteCounter._get_enzyme_instance(en) for en in enzyme_names]
        self.is_linear = is_linear
        self.tip_size = tip_size

    @staticmethod
    def _get_enzyme_instance(enz_name):
        """
        Fetch an instance of a given restriction enzyme by its name.

        :param enz_name: the case-sensitive name of the enzyme
        :return: RestrictionType the enzyme instance
        """
        try:
            # this has to match exactly
            return getattr(Restriction, enz_name)
        except AttributeError as e:
            # since we're being strict about enzyme names, be helpful with errors
            similar = []
            for a in dir(Restriction):
                score = SequenceMatcher(None, enz_name.lower(), a.lower()).ratio()
                if score >= 0.8:
                    similar.append(a)
            raise UnknownEnzymeException(enz_name, similar)

    def _count(self, seq):
        return sum(len(en.search(seq, self.is_linear)) for en in self.enzymes)

    def count_sites(self, seq):
        """
        Count the number of sites found in the given sequence, where sites from
        all specified enzymes are combined

        :param seq: Bio.Seq object
        :return: the total number of sites
        """
        if self.tip_size:
            seq_len = len(seq)
            if seq_len < 2*self.tip_size:
                # small contigs simply divide their extent in half
                l_tip = seq[:seq_len/2]
                r_tip = seq[-seq_len/2:]
            else:
                l_tip = seq[:self.tip_size]
                r_tip = seq[-self.tip_size:]
            # left and right tip counts
            sites = [self._count(l_tip), self._count(r_tip)]
        else:
            # one value for whole sequence
            sites = self._count(seq)

        return sites


class ExtentGrouping:

    def __init__(self, seq_info, bin_size):
        self.bins = []
        self.bin_size = bin_size
        self.map = []
        self.borders = []
        self.centers = []
        self.total_bins = 0

        for n, seq in tqdm.tqdm(enumerate(seq_info), total=len(seq_info), desc='Making bins'):

            if seq.length == 0:
                raise ZeroLengthException(seq.id)

            # integer bin estimation
            num_bins = seq.length / bin_size
            if num_bins == 0:
                num_bins += 1
            # handle non-integer discrepancy by contracting/expanding all bins equally
            # the threshold between contract/expand being half a bin size
            if seq.length % bin_size != 0 and seq.length/float(bin_size) - num_bins >= 0.5:
                num_bins += 1

            edges = np.linspace(0, seq.length, num_bins+1, endpoint=True, dtype=np.int)

            self.bins.append(num_bins)

            # Per reference coordinate pairs (bin_edge, map_index)
            first_bin = self.total_bins
            last_bin = first_bin + num_bins
            self.map.append(np.vstack((edges[1:], np.arange(first_bin, last_bin))).T)
            self.borders.append(np.array([first_bin, last_bin], dtype=np.int))

            self.total_bins += num_bins

            c_nk = edges[:-1] + 0.5*(edges[1] - edges[0]) - 0.5*seq.length
            self.centers.append(c_nk.reshape((1, len(c_nk))))
            # logger.debug('{}: {} bins'.format(n, num_bins))

        self.bins = np.array(self.bins)


class SeqOrder:

    FORWARD = 1
    REVERSE = -1

    ACCEPTED = True
    EXCLUDED = False

    STRUCT_TYPE = np.dtype([('pos', np.int32), ('ori', np.int8), ('mask', np.bool), ('length', np.int32)])
    INDEX_TYPE = np.dtype([('index', np.int32), ('ori', np.int8)])

    def __init__(self, seq_info):
        """
        Initial order is determined by the order of supplied sequence information dictionary. Sequences
        are given surrogate ids using consecutive integers. Member functions expect surrogate ids
        not original names.

        The class also retains orientation and masking state. Orientation defines whether a sequence
        should be in its original direction (as read in) (1) or reverse complemented (-1).

        Masking state defines whether a input sequence shall been excluded from further consideration.
        (accepted=1, excluded=0)

        :param seq_info: sequence information dictionary
        """
        self._positions = None
        _ord = np.arange(len(seq_info), dtype=np.int32)
        self.order = np.array(
            [(_ord[i], SeqOrder.FORWARD, SeqOrder.ACCEPTED, seq_info[i].length) for i in xrange(len(_ord))],
            dtype=SeqOrder.STRUCT_TYPE)

        self._update_positions()

    @staticmethod
    def asindex(_ord):
        """
        Convert a simple list or ndarray of indices, to INDEX_TYPE array with default forward orientation.

        :param _ord: list/ndarray of indices
        :return: INDEX_TYPE array
        """
        assert isinstance(_ord, (list, np.ndarray)), 'input must be a list or ndarray'
        return np.array(zip(_ord, np.ones_like(_ord, dtype=np.bool)), dtype=SeqOrder.INDEX_TYPE)

    def _update_positions(self):
        """
        An optimisation, whenever the positional state changes, this method must be called to
        maintain the current state in a separate array. This avoids unnecessary recalculation
        overhead.
        """
        # Masked sequences last, then by current position.
        sorted_indices = np.lexsort([self.order['pos'], ~self.order['mask']])
        for n, i in enumerate(sorted_indices):
            self.order[i]['pos'] = n
        self._positions = np.argsort(self.order['pos'])

    def remap_gapless(self, gapless_indices):
        """
        Recover the original, potentially sparse (gapped) indices from a dense (gapless) set
        of indices. Gaps originate from sequences being masked in the order. External tools
        often expect and return dense indices. When submitting changes to the current order
        state, it is important to first apply this method and reintroduce any gaps.

        Both a list/array of indices or a INDEX_TYPE array can be passed.

        :param gapless_indices: dense list of indices or an ndarray of type INDEX_TYPE
        :return: remappped indices with gaps (of a similar type to input)
        """
        n = 0
        shift = []
        for oi in self.order:
            if not oi['mask']:
                n += 1
            else:
                shift.append(n)
        shift = np.array(shift)

        if isinstance(gapless_indices, np.ndarray) and gapless_indices.dtype == SeqOrder.INDEX_TYPE:
            remapped = []
            for oi in gapless_indices:
                remapped.append((oi['index'] + shift[oi['index']], oi['ori']))
            return np.array(remapped, dtype=SeqOrder.INDEX_TYPE)

        else:
            # assume an iterable collection
            remapped = []
            for oi in gapless_indices:
                remapped.append(oi + shift[oi])
            return np.array(remapped)

    def accepted_positions(self, copy=True):
        """
        The current positional order of only those sequences which have not been excluded by the mask.
        :param copy: return a copy
        :return: all accepted positons in order of index
        """
        return self.all_positions(copy=copy)[:self.count_accepted()]

    def all_positions(self, copy=True):
        """
        The current positional order of all sequences. Internal logic relegates masked sequences to always come
        last and ascending surrogate id order.

        :param copy: return a copy of the positions
        :return: all positions in order of index, masked or not.
        """
        if copy:
            _p = self._positions.copy()
        else:
            _p = self._positions
        return _p

    @staticmethod
    def double_order(_ord):
        """
        For doublet maps, the stored order must be re-expanded to reference the larger (2x) map.

        :param _ord:
        :return: expanded order
        """
        return np.array([[2*oi, 2*oi+1] for oi in _ord]).ravel()

    def gapless_positions(self):
        """
        A dense index range representing the current positional order without masked sequences. Therefore
        the returned array does not contain surrogate ids, but rather the relative positions of unmasked
        sequences, when all masked sequences have been discarded.

        :return: a dense index range of positional order, once all masked sequences have been discarded.
        """
        # accumulated shift from gaps
        gap_shift = np.cumsum(~self.order['mask'])
        # just unmasked sequences
        # TODO test whether we can reenable this code.
        # _p = self._positions[:self.count_accepted()].copy()
        _p = np.argsort(self.order['pos'])
        _p = _p[:self.count_accepted()]
        # removing gaps leads to a dense range of indices
        _p -= gap_shift[_p]
        return _p

    def set_mask_only(self, _mask):
        """
        Set the mask state of all sequences, where indices in the mask map to
        sequence surrogate ids.

        :param _mask: mask array or list, boolean or 0/1 valued
        """
        _mask = np.asarray(_mask, dtype=np.bool)
        assert len(_mask) == len(self.order), 'supplied mask must be the same length as existing order'
        assert np.all((_mask == SeqOrder.ACCEPTED) | (_mask == SeqOrder.EXCLUDED)), \
            'new mask must be {} or {}'.format(SeqOrder.ACCEPTED, SeqOrder.EXCLUDED)

        # assign mask
        self.order['mask'] = _mask
        self._update_positions()

    def set_order_only(self, _ord, implicit_excl=False):
        """
        Convenience method to set the order using a list or 1D ndarray. Orientations will
        be assumed as all forward (+1).

        :param _ord: a list or ndarray of surrogate ids
        :param implicit_excl: implicitly extend the order to include unmentioned excluded sequences.
        """
        assert isinstance(_ord, (list, np.ndarray)), 'Wrong type supplied, order must be a list or ndarray'
        if isinstance(_ord, np.ndarray):
            _ord = np.ravel(_ord)
            assert np.ndim(_ord) == 1, 'orders as numpy arrays must be 1-dimensional'
        # augment the order to include default orientations
        _ord = SeqOrder.asindex(_ord)
        self.set_order_and_orientation(_ord, implicit_excl=implicit_excl)

    def set_order_and_orientation(self, _ord, implicit_excl=False):
        """
        Set only the order, while ignoring orientation. An ordering is defined
        as a 1D array of the structured type INDEX_TYPE, where elements are the
        position and orientation of each indices.

        NOTE: This definition can be the opposite of what is returned by some
        ordering methods, and np.argsort(_v) should inverse the relation.

        NOTE: If the order includes only active sequences, setting implicit_excl=True
        the method will implicitly assume unmentioned ids are those currently
        masked. An exception is raised if a masked sequence is included in the order.

        :param _ord: 1d ordering
        :param implicit_excl: implicitly extend the order to include unmentioned excluded sequences.
        """
        assert _ord.dtype == SeqOrder.INDEX_TYPE, 'Wrong type supplied, _ord should be of INDEX_TYPE'

        if len(_ord) < len(self.order):
            # some sanity checks
            assert implicit_excl, 'Use implicit_excl=True for automatic handling ' \
                                  'of orders only mentioning accepted sequences'
            assert len(_ord) == self.count_accepted(), 'new order must mention all ' \
                                                       'currently accepted sequences'
            # those surrogate ids mentioned in the order
            mentioned = set(_ord['index'])
            assert len(mentioned & set(self.excluded())) == 0, 'new order and excluded must not ' \
                                                               'overlap when using implicit assignment'
            assert len(mentioned ^ set(self.accepted())) == 0, 'incomplete new order supplied,' \
                                                               'missing accepted ids'
            # assign the new orders
            self.order['pos'][_ord['index']] = np.arange(len(_ord), dtype=np.int32)
            self.order['ori'][_ord['index']] = _ord['ori']
            # mask remaining, unmentioned indices
            _mask = np.zeros_like(self.mask_vector(), dtype=np.bool)
            _mask[_ord['index']] = True
            self.set_mask_only(_mask)
        else:
            # just a simple complete order update
            assert len(_ord) == len(self.order), 'new order was a different length'
            assert len(set(_ord['index']) ^ set(self.accepted())) == 0, 'incomplete new order supplied,' \
                                                                        'missing accepted ids'
            self.order['pos'][_ord['index']] = np.arange(len(_ord), dtype=np.int32)
            self.order['ori'][_ord['index']] = _ord['ori']

        self._update_positions()

    def accepted_order(self):
        """
        :return: an INDEX_TYPE array of the order and orientation of the currently accepted sequences.
        """
        idx = np.where(self.order['mask'])
        ori = np.ones(self.count_accepted(), dtype=np.int)
        return np.array(zip(idx, ori), dtype=SeqOrder.INDEX_TYPE)

    def mask_vector(self):
        """
        :return: the current mask vector
        """
        return self.order['mask']

    def mask(self, _id):
        """
        Mask an individual sequence by its surrogate id

        :param _id: the surrogate id of sequence
        """
        self.order[_id]['mask'] = False
        self._update_positions()

    def count_accepted(self):
        """
        :return: the current number of accepted (unmasked) sequences
        """
        return self.order['mask'].sum()

    def count_excluded(self):
        """
        :return: the current number of excluded (masked) sequences
        """
        return len(self.order) - self.count_accepted()

    def accepted(self):
        """
        :return: the list surrogate ids for currently accepted sequences
        """
        return np.where(self.order['mask'])[0]

    def excluded(self):
        """
        :return: the list surrogate ids for currently excluded sequences
        """
        return np.where(~self.order['mask'])[0]

    def flip(self, _id):
        """
        Flip the orientation of the sequence

        :param _id: the surrogate id of sequence
        """
        self.order[_id]['ori'] *= -1

    def lengths(self, exclude_masked=False):
        # type: (bool) -> np.ndarray
        """
        Sequence lengths

        :param exclude_masked: True include only umasked sequencces
        :return: the lengths of sequences
        """
        if exclude_masked:
            return self.order['length'][self.order['mask']]
        else:
            return self.order['length']

    def shuffle(self):
        """
        Randomize order
        """
        np.random.shuffle(self.order['pos'])
        self._update_positions()

    def before(self, a, b):
        """
        Test if a comes before another sequence in order.

        :param a: surrogate id of sequence a
        :param b: surrogate id of sequence b
        :return: True if a comes before b
        """
        assert a != b, 'Surrogate ids must be different'
        return self.order['pos'][a] < self.order['pos'][b]

    def intervening(self, a, b):
        """
        For the current order, calculate the length of intervening
        sequences between sequence a and sequence b.

        :param a: surrogate id of sequence a
        :param b: surrogate id of sequence b
        :return: total length of sequences between a and b.
        """
        assert a != b, 'Surrogate ids must be different'

        pa = self.order['pos'][a]
        pb = self.order['pos'][b]
        if pa > pb:
            pa, pb = pb, pa
        inter_ix = self._positions[pa+1:pb]
        return np.sum(self.order['length'][inter_ix])


@jit(int64(int64[:, :], int64), nopython=True)
def find_nearest_jit(group_sites, x):
    """
    Find the nearest site from a given position on a contig.

    :param group_sites:
    :param x: query position
    :return: tuple of site and group number
    """
    ix = np.searchsorted(group_sites[:, 0], x)
    if ix == len(group_sites):
        # raise RuntimeError('find_nearest: {} didnt fit in {}'.format(x, group_sites))
        return group_sites[-1, 1]
    return group_sites[ix, 1]


@jit(void(float64[:, :, :, :], float64[:, :], float64[:], float64))
def fast_norm_tipbased_bylength(coords, data, tip_lengths, tip_size):
    """
    In-place normalisation of the sparse 4D matrix used in tip-based maps.

    As tip-based normalisation is slow for large matrices, the inner-loop has been
    moved to a Numba method.

    :param coords: the COO matrix coordinate member variable (4xN array)
    :param data:  the COO matrix data member variable (1xN array)
    :param tip_lengths: per-element min(sequence_length, tip_size)
    :param tip_size: tip size used in map
    """
    for ii in xrange(coords.shape[1]):
        i, j = coords[:2, ii]
        data[ii] *= tip_size**2 / (tip_lengths[i] * tip_lengths[j])


@jit(void(float64[:, :, :, :], float64[:, :], float64[:, :]))
def fast_norm_tipbased_bysite(coords, data, sites):
    """
    In-place normalisation of the sparse 4D matrix used in tip-based maps.

    As tip-based normalisation is slow for large matrices, the inner-loop has been
    moved to a Numba method.

    :param coords: the COO matrix coordinate member variable (4xN array)
    :param data:  the COO matrix data member variable (1xN array)
    :param sites: per-element min(sequence_length, tip_size)
    """
    for ii in xrange(coords.shape[1]):
        i, j, k, l = coords[:, ii]
        data[ii] *= 1.0/(sites[i, k] * sites[j, l])


class IndexedFasta(Mapping):
    """
    Provides indexed access to a sequence file in Fasta format which can be compressed using bgzip.

    The SeqIO.index_db method is used for access, and the temporary index is stored in the system
    tmp directory by default. At closing, the index will be closed and the temporary file removed.

    Wrapped with contextlib.closing(), instances of this object are with() compatible.
    """

    def __init__(self, fasta_file, tmp_path=None):
        """
        :param fasta_file: the input fasta file to access
        :param tmp_path: temporary directory for creating index
        """
        if tmp_path is None:
            tmp_path = tempfile.gettempdir()
        elif not os.path.exists(tmp_path):
            raise IOError('specified temporary path [{}] does not exist'.format(tmp_path))
        self._tmp_file = os.path.join(tmp_path, '{}.seqdb'.format(str(uuid.uuid4())))
        self._fasta_file = fasta_file
        self._index = SeqIO.index_db(self._tmp_file, self._fasta_file, 'fasta')

    def __getitem__(self, _id):
        """
        Access a sequence object by its Fasta identifier.

        :param _id: fasta sequence identifier
        :return: SeqRecord object representing the sequence
        """
        return self._index[_id]

    def __iter__(self):
        """
        :return: iterator over sequence identifiers (keys)
        """
        return iter(self._index)

    def __len__(self):
        """
        :return: the number of sequences in the file
        """
        return len(self._index)

    def close(self):
        """
        Close the index and remove the associated temporary file
        """
        if self._index:
            self._index.close()
        if os.path.exists(self._tmp_file):
            os.unlink(self._tmp_file)


def count_fasta_sequences(file_name):
    """
    Estimate the number of fasta sequences in a file by counting headers. Decompression is automatically attempted
    for files ending in .gz. Counting and decompression is by why of subprocess calls to grep and gzip. Uncompressed
    files are also handled. This is about 8 times faster than parsing a file with BioPython and 6 times faster
    than reading all lines in Python.

    :param file_name: the fasta file to inspect
    :return: the estimated number of records
    """
    if file_name.endswith('.gz'):
        proc_uncomp = subprocess.Popen(['gzip', '-cd', file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_read = subprocess.Popen(['grep', r'^>'], stdin=proc_uncomp.stdout, stdout=subprocess.PIPE)
    else:
        proc_read = subprocess.Popen(['grep', r'^>', file_name], stdout=subprocess.PIPE)

    n = 0
    for _ in proc_read.stdout:
        n += 1
    return n


class ContactMap:

    def __init__(self, bam_file, enzymes, seq_file, min_insert, min_mapq=0, min_len=0, min_sig=1, max_fold=None,
                 random_seed=None, strong=None, bin_size=None, tip_size=None, precount=False, med_alpha=10,
                 verbose=False):

        self.strong = strong
        self.bin_size = bin_size
        self.min_mapq = min_mapq
        self.min_insert = min_insert
        self.min_len = min_len
        self.min_sig = min_sig
        self.max_fold = max_fold
        self.random_state = np.random.RandomState(random_seed)
        self.strong = strong
        self.seq_info = []
        self.seq_map = None
        self.seq_file = seq_file
        self.grouping = None
        self.extent_map = None
        self.order = None
        self.tip_size = tip_size
        self.precount = precount
        self.total_reads = None
        self.cov_info = None
        self.med_alpha = med_alpha
        self.processed_map = None
        self.primary_acceptance_mask = None
        self.bisto_scale = None
        self.seq_analyzer = None
        self.enzymes = enzymes

        # build a dictionary of sites in each reference first
        fasta_info = {}
        with io_utils.open_input(seq_file) as multi_fasta:
            # prepare the site counter for the given experimental conditions
            site_counter = SiteCounter(enzymes, tip_size, is_linear=True)
            # get an estimate of sequences for progress
            fasta_count = count_fasta_sequences(seq_file)
            for seqrec in tqdm.tqdm(SeqIO.parse(multi_fasta, 'fasta'), total=fasta_count, desc='Analyzing sites'):
                if len(seqrec) < min_len:
                    continue
                fasta_info[seqrec.id] = {'sites': site_counter.count_sites(seqrec.seq),
                                         'length': len(seqrec)}

        # now parse the header information from bam file
        with pysam.AlignmentFile(bam_file, 'rb') as bam:

            # test that BAM file is the correct sort order
            if 'SO' not in bam.header['HD'] or bam.header['HD']['SO'] != 'queryname':
                raise IOError('BAM file must be sorted by read name')

            # determine the set of active sequences
            # where the first filtration step is by length
            ref_count = {'seq_missing': 0, 'too_short': 0}
            offset = 0
            logger.info('Reading sequences...')
            for n, (rname, rlen) in enumerate(zip(bam.references, bam.lengths)):

                # minimum length threshold
                if rlen < min_len:
                    ref_count['too_short'] += 1
                    continue

                try:
                    fa = fasta_info[rname]
                except KeyError:
                    logger.info('Sequence: "{}" was not present in reference fasta'.format(rname))
                    ref_count['seq_missing'] += 1
                    continue

                assert fa['length'] == rlen, \
                    'Sequence lengths in {} do not agree: bam {} fasta {}'.format(rname, fa['length'], rlen)

                self.seq_info.append(SeqInfo(offset, n, rname, rlen, fa['sites']))

                offset += rlen

            # total extent covered
            self.total_len = offset
            self.total_seq = len(self.seq_info)

            # all sequences begin as active in mask
            # when filtered, mask[i] = False
            self.current_mask = np.ones(self.total_seq, dtype=np.bool)

            if self.total_seq == 0:
                logger.info('No sequences in BAM found in FASTA')
                raise ParsingError('No sequences in BAM found in FASTA')

            logger.info('Accepted {} sequences covering {} bp'.format(self.total_seq, self.total_len))
            logger.info('References excluded: {}'.format(ref_count))

            if self.bin_size:
                logger.info('Determining bins...')
                self.grouping = ExtentGrouping(self.seq_info, self.bin_size)

            logger.info('Counting reads in bam file...')

            if self.precount:
                self.total_reads = bam.count(until_eof=True)
                logger.info('BAM file contains {0} alignments'.format(self.total_reads))
            else:
                logger.info('Skipping pre-count of BAM file, no ETA will be offered')

            # initialise the order
            self.order = SeqOrder(self.seq_info)

            # accumulate
            self._bin_map(bam)

            # create an initial acceptance mask
            self.set_primary_acceptance_mask(verbose=verbose)

    def _bin_map(self, bam):
        """
        Accumulate read-pair observations from the supplied BAM file.
        Maps are initialized here. Logical control is achieved through initialisation of the
        ContactMap instance, rather than supplying this function arguments.

        :param bam: this instance's open bam file.
        """
        import tqdm

        def _simple_match(r):
            return not r.is_unmapped and r.mapq >= _mapq

        def _strong_match(r):
            return strong_match(r, self.strong, True, _mapq)

        # set-up match call
        _matcher = _strong_match if self.strong else _simple_match

        def _on_tip_withlocs(p1, p2, l1, l2, _tip_size):
            tailhead_mat = np.zeros((2, 2), dtype=np.uint32)
            i = None
            j = None

            # contig1 tips won't overlap
            if l1 > 2 * _tip_size:
                if p1 < _tip_size:
                    i = 0
                elif p1 > l1 - _tip_size:
                    i = 1

            # contig1 tips will overlap
            else:
                # assign to whichever end is closest
                if p1 < l1 - p1:
                    i = 0
                elif l1 - p1 < p1:
                    i = 1

            # only bother with second tip assignment if the first was ok
            if i is not None:

                # contig2 tips won't overlap
                if l2 > 2 * _tip_size:
                    if p2 < _tip_size:
                        j = 0
                    elif p2 > l2 - _tip_size:
                        j = 1

                # contig2 tips will overlap
                else:
                    # assign to whichever end is closest
                    if p2 < l2 - p2:
                        j = 0
                    elif l2 - p2 < p2:
                        j = 1

            tailhead_mat[i, j] = 1
            return i is not None and j is not None, tailhead_mat

        def _always_true(*args):
            return True, 1

        _on_tip = _always_true if not self.is_tipbased() else _on_tip_withlocs

        # initialise a sparse matrix for accumulating the map
        if not self.is_tipbased():
            # just a basic NxN sparse array for normal whole-sequence binning
            _seq_map = simple_sparse.Sparse2DAccumulator(self.total_seq)
        else:
            # each tip is tracked separately resulting in the single count becoming a 2x2 interaction matrix.
            # therefore the tensor has dimension NxNx2x2
            _seq_map = simple_sparse.Sparse4DAccumulator(self.total_seq)

        # if binning also requested, initialise another sparse matrix
        if self.bin_size:
            logger.info('Initialising contact map of {0}x{0} fragment bins, '
                        'representing {1} bp over {2} sequences'.format(self.grouping.total_bins,
                                                                        self.total_len, self.total_seq))
            _extent_map = simple_sparse.Sparse2DAccumulator(self.grouping.total_bins)
            _grouping_map = self.grouping.map
        else:
            _grouping_map = None
            _extent_map = None

        with tqdm.tqdm(total=self.total_reads) as pbar:

            # locals for read filtering
            _min_sep = self.min_insert
            _mapq = self.min_mapq

            _idx = self.make_reverse_index('refid')

            # locals for tip checking
            _len = bam.lengths
            _tip_size = self.tip_size

            counts = OrderedDict({
                'accepted': 0,
                'not_tip': 0,
                'short_insert': 0,
                'ref_excluded': 0,
                'median_excluded': 0,
                'end_buffered': 0,
                'poor_match': 0})

            bam.reset()
            bam_iter = bam.fetch(until_eof=True)
            while True:

                try:
                    r1 = bam_iter.next()
                    pbar.update()
                    while True:
                        # read records until we get a pair
                        r2 = bam_iter.next()
                        pbar.update()
                        if r1.query_name == r2.query_name:
                            break
                        r1 = r2
                except StopIteration:
                    break

                if r1.reference_id not in _idx or r2.reference_id not in _idx:
                    counts['ref_excluded'] += 1
                    continue

                if not _matcher(r1) or not _matcher(r2):
                    counts['poor_match'] += 1
                    continue

                # if no_go[r1.reference_id][r1.reference_start:r1.reference_end] or \
                #         no_go[r2.reference_id][r2.reference_start:r2.reference_end]:
                #     counts['median_excluded'] += 1
                #     continue

                if r1.is_read2:
                    r1, r2 = r2, r1

                # # accept only inward facing read pairs
                # if not ((r1.is_reverse and not r2.is_reverse) or (not r1.is_reverse and r2.is_reverse)):
                #     # counts['skipped'] += 1
                #     continue

                # use 5-prime base depending on orientation
                r1pos = r1.pos if not r1.is_reverse else r1.pos + r1.alen
                r2pos = r2.pos if not r2.is_reverse else r2.pos + r2.alen

                # filter inserts deemed "short" which tend to be heavily WGS signal
                if _min_sep and r1.is_proper_pair:
                    ins_len = r2.pos - r1.pos
                    if ins_len < _min_sep:
                        counts['short_insert'] += 1
                        continue

                # get reference lengths
                l1 = _len[r1.reference_id]
                l2 = _len[r2.reference_id]

                # get internal indices
                ix1 = _idx[r1.reference_id]
                ix2 = _idx[r2.reference_id]

                # maintain just a half-matrix
                if ix2 < ix1:
                    ix1, ix2 = ix2, ix1
                    r1pos, r2pos = r2pos, r1pos
                    l1, l2 = l2, l1

                if _extent_map:
                    b1 = find_nearest_jit(_grouping_map[ix1], r1pos)
                    b2 = find_nearest_jit(_grouping_map[ix2], r2pos)

                    # maintain half-matrix
                    if b1 > b2:
                        b1, b2 = b2, b1

                    # tally all mapped reads for binned map, not just those considered in tips
                    _extent_map[b1, b2] += 1

                # for seq-map, we may reject reads outside of a defined tip region
                tip_info = _on_tip(r1pos, r2pos, l1, l2, _tip_size)
                if not tip_info[0]:
                    counts['not_tip'] += 1
                    continue

                counts['accepted'] += 1

                _seq_map[ix1, ix2] += tip_info[1]

        # default to always making matrices symmetric
        if self.bin_size:
            self.extent_map = _extent_map.get_coo()
            del _extent_map

        self.seq_map = _seq_map.get_coo()
        del _seq_map

        logger.info('Pair accounting: {}'.format(counts))
        logger.info('Total extent map weight {}'.format(self.map_weight()))

    @staticmethod
    def get_fields():
        """
        :return: the list of fields used in seq_info dict.
        """
        return SeqInfo._fields

    def make_reverse_index(self, field_name):
        """
        Make a reverse look-up (dict) from the chosen field in seq_info to the internal index value
        of the given sequence. Non-unique fields will raise an exception.

        :param field_name: the seq_info field to use as the reverse.
        :return: internal array index of the sequence
        """
        rev_idx = {}
        for n, seq in enumerate(self.seq_info):
            fv = getattr(seq, field_name)
            if fv in rev_idx:
                raise RuntimeError('field contains non-unique entries, a 1-1 mapping cannot be made')
            rev_idx[fv] = n
        return rev_idx

    def map_weight(self):
        """
        :return: the total map weight (sum ij)
        """
        return self.seq_map.sum()

    def is_empty(self):
        """
        :return: True if the map has zero weight
        """
        return self.map_weight() == 0

    def is_tipbased(self):
        """
        :return: True if the seq_map is a tip-based 4D tensor
        """
        return self.tip_size is not None

    def find_order(self, _map, inverse_method='inverse', seed=None, runs=5, work_dir='.', verbose=False):
        """
        Using LKH TSP solver, find the best ordering of the sequence map in terms of proximity ligation counts.
        Here, it is assumed that sequence proximity can be inferred from number of observed trans read-pairs, where
        an inverse relationship exists.

        :param _map: the seq map to analyze
        :param inverse_method: the chosen inverse method for converting count (similarity) to distance
        :param seed: specify a random seed, otherwise one is generated at runtime
        :param runs: number of individual runs of lkh to perform
        :param work_dir: working directory
        :param verbose: debug output
        :return: the surrogate ids in optimal order
        """

        # a minimum of three sequences is required to run LKH
        if _map.shape[0] < 3:
            raise TooFewException()

        # we'll supply a partially initialized distance function
        dist_func = partial(ordering.similarity_to_distance, method=inverse_method,
                            alpha=1.2, beta=1, verbose=verbose)

        if seed is None:
            seed = make_random_seed()
            if verbose:
                print 'Using random seed: {}'.format(seed)

        with open(os.path.join(work_dir, 'lkh.log'), 'w+') as stdout:

            control_base_name = os.path.join(work_dir, 'lkh_run')

            if self.is_tipbased():

                lkh_o = ordering.lkh_order(_map, control_base_name, lkh_exe=ext_path('LKH-3.0'), precision=1, seed=seed,
                                           runs=runs, pop_size=50, dist_func=dist_func, special=False, stdout=stdout,
                                           verbose=True, fixed_edges=[(i, i+1) for i in xrange(1, _map.shape[0], 2)])

                # To solve this with TSP, doublet tips use a graph transformation, where each node comes a pair. Pairs
                # possess fixed inter-connecting edges which must be included in any solution tour.
                # Eg. node 0 -> (0,1) or node 1 -> (2,3). The fixed paths are undirected, depending on which direction
                # is traversed, defines the orientation of the sequence.

                # 1.  pair adjacent nodes by reshape the 1D array into a two-column array of half the length
                lkh_o = lkh_o.reshape(lkh_o.shape[0]/2, 2)
                # 2. convert to surrogate ids and infer orientation from paths taken through doublets.
                #   0->1 forward (+1): 1->0 reverse (-1).
                lkh_o = np.fromiter(((oi[0]/2, oi[1]-oi[0]) for oi in lkh_o), dtype=SeqOrder.INDEX_TYPE)

            else:

                lkh_o = ordering.lkh_order(_map, control_base_name, lkh_exe=ext_path('LKH-3.0'), precision=1, seed=seed,
                                           runs=runs, pop_size=50, dist_func=dist_func, special=False, stdout=stdout,
                                           verbose=True)

                # for singlet tours, no orientation can be inferred.
                lkh_o = np.fromiter(((oi, 1) for oi in lkh_o), dtype=SeqOrder.INDEX_TYPE)

        # lkh ordering references the supplied matrix indices, not the surrogate ids.
        # we must map this consecutive set to the contact map indices.
        lkh_o = self.order.remap_gapless(lkh_o)

        return lkh_o

    def get_primary_acceptance_mask(self):
        assert self.primary_acceptance_mask is not None, 'Primary acceptance mask has not be initialized'
        return self.primary_acceptance_mask.copy()

    def set_primary_acceptance_mask(self, min_len=None, min_sig=None, max_fold=None, update=False, verbose=False):
        """
        Determine and set the filter mask using the specified constraints across the entire
        contact map. The mask is True when a sequence is considered acceptable wrt to the
        constraints. The mask is also returned by the function for convenience.

        :param min_len: override instance value for minimum sequence length
        :param min_sig: override instance value for minimum off-diagonal signal (counts)
        :param max_fold: maximum locally-measured fold-coverage to permit
        :param update: replace the current primary mask if it exists
        :param verbose: debug output
        :return: an acceptance mask over the entire contact map
        """
        assert max_fold is None, 'Filtering on max_fold is currently disabled'

        # If parameter based critiera were unset, use instance member values set at instantiation time
        if not min_len:
            min_len = self.min_len
        if not min_sig:
            min_sig = self.min_sig

        assert min_len, 'Filtering criteria min_len is None'
        assert min_sig, 'Filtering criteria min_sig is None'

        if verbose:
            print 'Setting primary acceptance mask'
            print 'Filtering criterion min_len: {} min_sig: {}'.format(min_len, min_sig)

        # simply return the current mask if it has already been determined
        # and an update is not requested
        if not update and self.primary_acceptance_mask is not None:
            if verbose:
                print 'Using existing mask'
            return self.get_primary_acceptance_mask()

        acceptance_mask = np.ones(self.total_seq, dtype=np.bool)

        # mask for sequences shorter than limit
        _mask = self.order.lengths() >= min_len
        if verbose:
            print 'Minimum length threshold removing:', self.total_seq - _mask.sum()
        acceptance_mask &= _mask

        # mask for sequences weaker than limit
        if self.is_tipbased():
            signal = simple_sparse.max_offdiag_4d(self.seq_map)
        else:
            signal = simple_sparse.max_offdiag(self.seq_map)
        _mask = signal >= min_sig
        if verbose:
            print 'Minimum signal threshold removing:', self.total_seq - _mask.sum()
        acceptance_mask &= _mask

        # TODO dropped degen removal for now
        # # mask for sequences considered to have aberrant coverage.
        # if max_fold:
        #     seq_analyzer = SequenceAnalyzer(self.seq_map, self.seq_report, self.seq_info, self.tip_size)
        #     degen_seqs = seq_analyzer.report_degenerates(max_fold, verbose=False)
        #     _mask = ~degen_seqs['status']
        #     if verbose:
        #         print 'Degenerate removal', self.total_seq - _mask.sum()
        #     acceptance_mask &= _mask

        # retain the union of all masks.
        self.primary_acceptance_mask = acceptance_mask

        if verbose:
            print 'Accepted sequences', self.primary_acceptance_mask.sum()

        return self.get_primary_acceptance_mask()

    def prepare_seq_map(self, norm=True, bisto=False, mean_type='geometric', verbose=False):
        """
        Prepare the sequence map (seq_map) by application of various filters and normalisations.

        :param norm: normalisation by sequence lengths
        :param bisto: make the output matrix bistochastic
        :param mean_type: when performing normalisation, use "geometric, harmonic or arithmetic" mean.
        :param verbose: debug output
        """

        if verbose:
            print 'Preparing sequence map'
            print 'Full map size', self.seq_map.shape

        _mask = self.get_primary_acceptance_mask()

        self.order.set_mask_only(_mask)

        if self.order.count_accepted() < 1:
            raise NoneAcceptedException()

        _map = self.seq_map.astype(np.float)

        # apply length normalisation if requested
        if norm:
            _map = self._norm_seq(_map, self.is_tipbased(), mean_type=mean_type, use_sites=True, verbose=verbose)
            if verbose:
                print 'Map normalized'

        # make map bistochastic if requested
        if bisto:
            # TODO balancing may be better done after compression
            _map, scl = self._bisto_seq(_map, self.is_tipbased(), verbose)
            # retain the scale factors
            self.bisto_scale = scl
            if verbose:
                print 'Map balanced'

        # cache the results for optional quick access
        self.processed_map = _map

    def get_subspace(self, external_mask=None, marginalise=False, flatten=True, verbose=False):
        """
        Using an already normalized full seq_map, return a subspace as indicated by an external
        mask or if none is supplied, the full map without filtered elements.

        The supplied external mask must refer to all sequences in the map.

        :param external_mask: an external mask to combine with the existing primary mask
        :param marginalise: Assuming 4D NxNx2x2 tensor, sum 2x2 elements to become a 2D NxN
        :param flatten: convert a NxNx2x2 tensor to a 2Nx2N matrix
        :param verbose: debug output
        :return: subspace map
        """
        assert (not marginalise and not flatten) or np.logical_xor(marginalise, flatten), \
            'marginalise and flatten are mutually exclusive'

        # starting with the normalized map
        _map = self.get_processed_map()

        _mask = self.get_primary_acceptance_mask()
        if verbose:
            print 'Active sequences after primary filtering', _mask.sum()

        # from a union of the sequence filter and external mask
        if external_mask is not None:
            _mask &= external_mask
            if verbose:
                print 'Active sequences after applying external mask', _mask.sum()

        self.order.set_mask_only(_mask)

        # remove masked sequences from the map
        if self.order.count_accepted() < self.total_seq:
            if self.is_tipbased():
                _map = simple_sparse.compress_4d(_map, self.order.mask_vector())
            else:
                _map = simple_sparse.compress(_map.tocoo(), self.order.mask_vector())
            if verbose:
                print 'Subspace map dimensions', _map.shape

        # convert tip-based tensor to other forms
        if self.is_tipbased():
            if marginalise:
                if verbose:
                    print 'Marginalising NxNx2x2 tensor to NxN matrix'
                # sum counts of the 2x2 confusion matrices into 1 value
                _map = _map.sum(axis=(2, 3)).to_scipy_sparse()
            elif flatten:
                if verbose:
                    print 'Flattening NxNx2x2 tensor to 2Nx2N matrix'
                # convert the 4D map into a 2Nx2N 2D map.
                _map = simple_sparse.flatten_tensor_4d(_map)

        return _map

    def get_processed_map(self, permute=False, dtype=np.float, verbose=False):
        """
        Return a copy of the processed map

        :param permute: reorder with current state
        :param dtype: cast to another type
        :param verbose: debug output
        :return: the processed map
        """
        assert self.processed_map is not None, 'processed map has not been prepared!'

        _map = self.processed_map.astype(dtype)
        # reorder if requested
        if permute:
            _map = self._reorder_seq(_map, self.is_tipbased())
            if verbose:
                print 'Map reordered'
        return _map

    def get_extent_map(self, norm=True, bisto=False, permute=False, mean_type='geometric', verbose=False):
        """
        Return the extent map after applying specified processing steps. Masked sequences are always removed.

        :param norm: sequence length normalisation
        :param bisto: make map bistochastic
        :param permute: permute the map using current order
        :param mean_type: length normalisation mean (geometric, harmonic, arithmetic)
        :param verbose: debug information
        :return: processed extent map
        """

        if verbose:
            print 'Prepare extent map'
            print 'Original size', self.extent_map.shape

        _map = self.extent_map.astype(np.float)

        # apply length normalisation if requested
        if norm:
            _map = self._norm_extent(_map, mean_type)
            if verbose:
                print 'Map normalized'

        # if there are sequences to mask, remove them from the map
        if self.order.count_accepted() < self.total_seq:
            _map = self._compress_extent(_map)
            if verbose:
                print 'Post filtration map dimensions', _map.shape

        # make map bistochastic if requested
        if bisto:
            _map, scl = simple_sparse.kr_biostochastic(_map)
            if verbose:
                print 'Map balanced'

        # reorder using current order state
        if permute:
            _map = self._reorder_extent(_map)
            if verbose:
                print 'Map reordered'

        return _map

    def extent_to_seq(self):
        """
        Convert the extent map into a single-pixel per sequence "seq_map". This method
        is useful when only a tip based seq_map has been produced, and an analysis would be
        better done on a full accounting of mapping interactions across each sequences full
        extent.
        :return: a seq_map representing all counts across each sequence
        """
        m = self.extent_map.tocsr()
        m_out = simple_sparse.Sparse2DAccumulator(self.total_seq)
        cbins = np.cumsum(self.grouping.bins)
        a0 = 0
        for i in xrange(len(self.grouping.bins)):
            a1 = cbins[i]
            # sacrifice memory for significant speed up slicing below
            row_i = m[a0:a1, :].todense()
            b0 = 0
            for j in xrange(i, len(self.grouping.bins)):
                b1 = cbins[j]
                mij = row_i[:, b0:b1].sum()
                if mij == 0:
                    continue
                m_out[i, j] = int(mij)
                b0 = b1
            a0 = a1
        return m_out.get_coo()

    def _reorder_seq(self, _map, tip_based):
        """
        Reorder a simple sequence map using the supplied map

        :param _map: the map to reorder
        :return: ordered map
        """
        assert sp.isspmatrix(_map), 'reordering expects a sparse matrix type'

        _order = self.order.gapless_positions()
        if tip_based:
            _order = SeqOrder.double_order(_order)

        assert _map.shape[0] == _order.shape[0], 'supplied map and unmasked order are different sizes'
        p = sp.lil_matrix(_map.shape)
        for i in xrange(len(_order)):
            p[i, _order[i]] = 1.
        p = p.tocsr()
        return p.dot(_map.tocsr()).dot(p.T)

    def _bisto_seq(self, _map, tip_based, verbose=False):
        """
        Make a contact map bistochastic. This is another form of normslisation. Automatically
        handles 2D and 4D maps.
        :param _map: a map to balance (make bistochastic)
        :param tip_based: treat the supplied map as a tip-based tensor
        :param verbose: debug output
        :return: the balanced map
        """
        if verbose:
            print 'Balancing contact map'

        if tip_based:
            print 'Dimension during bisto:', _map.shape
            _map, scl = simple_sparse.kr_biostochastic_4d(_map)
        else:
            _map, scl = simple_sparse.kr_biostochastic(_map)
        return _map, scl

    def _get_sites(self):
        _sites = np.array([si.sites for si in self.seq_info], dtype=np.float)
        # all sequences are assumed to have a minimum of 1 site -- even if not observed
        # TODO test whether it would be more accurate to assume that all sequences are under counted by 1.
        _sites[np.where(_sites == 0)] = 1
        return _sites

    def _norm_seq(self, _map, tip_based, use_sites=True, mean_type='geometric', verbose=False):
        """
        Normalise a simple sequence map in place by the geometric mean of interacting contig pairs lengths.
        The map is assumed to be in starting order.

        :param _map: the target map to apply normalisation
        :param tip_based: treat the supplied map as a tip-based tensor
        :param use_sites: normalise matrix counts using observed sites, otherwise normalise
        using sequence lengths as a proxy
        :param mean_type: for length normalisation, choice of mean (harmonic, geometric, arithmetic)
        :param verbose: debug output
        :return: normalized map
        """
        if use_sites:
            if verbose:
                print 'Doing site based normalisation'

            if tip_based:
                _map = _map.astype(np.float)
                print 'Dimension during norm:', _map.shape
                _sites = self._get_sites()
                fast_norm_tipbased_bysite(_map.coords, _map.data, _sites)

            else:
                # TODO implement non-tip version
                raise RuntimeError('currently unimplemented')

        else:
            if verbose:
                print 'Doing length based normalisation'

            if tip_based:
                _tip_lengths = np.minimum(self.tip_size, self.order.lengths()).astype(np.float)
                fast_norm_tipbased_bylength(_map.coords, _map.data, _tip_lengths, self.tip_size)

            else:
                _mean_func = mean_selector(mean_type)
                _len = self.order.lengths().astype(np.float)
                _map = _map.tolil().astype(np.float)
                for i in xrange(_map.shape[0]):
                    _map[i, :] /= np.fromiter((1e-3 * _mean_func(_len[i],  _len[j])
                                               for j in xrange(_map.shape[0])), dtype=np.float)
                _map = _map.tocsr()

        return _map

    def _norm_extent(self, _map, mean_type='geometric'):
        """
        Normalise a extent map in place by the geometric mean of interacting contig pairs lengths.

        :return: a normalized extent map in lil_matrix format
        """
        assert sp.isspmatrix(_map), 'Extent matrix is not a scipy matrix type'

        if _map.dtype not in {np.float, float}:
            _map = _map.astype(np.float)
        if not sp.isspmatrix_lil(_map):
            _map = _map.tolil()

        _mean_func = mean_selector(mean_type)
        _len = self.order.lengths().astype(np.float)
        _cbins = np.cumsum(self.grouping.bins)
        for row_i, col_dat in enumerate(_map.rows):
            i = np.searchsorted(_cbins, row_i, side='right')
            wi = np.fromiter((1e-3 * _mean_func(_len[i], _len[j])
                              for j in np.searchsorted(_cbins, col_dat, side='right')), dtype=np.float)
            _map.data[row_i] /= wi
        return _map

    def _reorder_extent(self, _map):
        """
        Reorder the extent map using current order.

        :return: sparse CSR format permutation of the given map
        """
        _order = self.order.gapless_positions()
        _bins = self.grouping.bins[self.order.mask_vector()]
        _ori = self.order.order['ori'][np.argsort(self.order.order['pos'])]

        # create a permutation matrix
        p = sp.lil_matrix(_map.shape)
        _shuf_bins = _bins[_order]
        for i, oi in enumerate(_order):
            j_off = _bins[:oi].sum()
            i_off = _shuf_bins[:i].sum()
            if _ori[i] > 0:
                for k in xrange(_bins[oi]):
                    p[i_off+k, j_off+k] = 1
            else:
                # rot90 those with reverse orientation
                _nb = _bins[oi]
                for k in xrange(_nb):
                    p[i_off+_nb-(k+1), j_off+k] = 1

        # permute the extent_map
        p = p.tocsr()
        return p.dot(_map.tocsr()).dot(p.T)

    def _compress_extent(self, _map):
        """
        Compress the extent map for each sequence that is presently masked. This will eliminate
        all bins which pertain to a given masked sequence.

        :return: a scipy.sparse.coo_matrix pertaining to only the unmasked sequences.
        """
        assert sp.isspmatrix(_map), 'Extent matrix is not a scipy sparse matrix type'
        if not sp.isspmatrix_coo(_map):
            _map = _map.tocoo()

        _order = self.order.order
        _bins = self.grouping.bins

        # build a list of every accepted element.
        # TODO this could be done as below, without the memory footprint of realising all elements
        s = 0
        accept_bins = []
        # accept_index = set(np.where(_mask)[0])
        for i in xrange(len(_order)):
            # if i in accept_index:
            if _order[i]['mask']:
                accept_bins.extend([j+s for j in xrange(_bins[i])])
            s += _bins[i]

        # use a hashable container for quicker lookup
        accept_bins = set(accept_bins)

        # collect those values not in the excluded rows/columns
        keep_row = []
        keep_col = []
        keep_data = []
        for i in xrange(_map.nnz):
            if _map.row[i] in accept_bins and _map.col[i] in accept_bins:
                keep_row.append(_map.row[i])
                keep_col.append(_map.col[i])
                keep_data.append(_map.data[i])

        # after removal of those intervening, shift the coordinates of affected bins
        # TODO this could be moved into the loop above
        _shift = np.cumsum((~_order['mask']) * _bins)
        _csbins = np.cumsum(_bins)
        for i in xrange(len(keep_row)):
            # rather than build a complete list of shifts across matrix, we'll
            # sacrifice some CPU and do lookups for the appropriate bin
            ir = np.searchsorted(_csbins, keep_row[i], side='right')
            ic = np.searchsorted(_csbins, keep_col[i], side='right')
            keep_row[i] -= _shift[ir]
            keep_col[i] -= _shift[ic]

        return sp.coo_matrix((keep_data, (keep_row, keep_col)), shape=_map.shape - _shift[-1])

    def plot_seqnames(self, fname, simple=True, permute=False, **kwargs):
        """
        Plot the contact map, annotating the map with sequence names. WARNING: This can often be too dense
        to be legible when there are many (1000s) of sequences.

        :param fname: output file name
        :param simple: True plot seq map, False plot the extent map
        :param permute: permute the map with the present order
        :param kwargs: additional options passed to plot()
        """
        if permute:
            seq_id_iter = self.order.accepted_positions()
        else:
            seq_id_iter = xrange(self.order.count_accepted())

        tick_labs = []
        for i in seq_id_iter:
            if self.order.order[i]['ori'] < 0:
                tick_labs.append('- {}'.format(self.seq_info[i].name))
            else:
                tick_labs.append('+ {}'.format(self.seq_info[i].name))

        if simple:
            step = 2 if self.is_tipbased() else 1
            tick_locs = xrange(2, step*self.order.count_accepted()+step, step)
        else:
            if permute:
                _cbins = np.cumsum(self.grouping.bins[self.order.accepted_positions()])
            else:
                _cbins = np.cumsum(self.grouping.bins[self.order.accepted()])
            tick_locs = _cbins - 0.5

        self.plot(fname, permute=permute, simple=simple, tick_locs=tick_locs, tick_labs=tick_labs, **kwargs)

    def _enable_clusters(self, clustering, cl_list=None, ordered_only=True, min_extent=None, verbose=False):
        """
        Given a clustering and list of cluster ids (or none), enable (unmask) the related sequences in
        the contact map. If a requested cluster has not been ordered, it will be dropped.

        :param clustering: a clustering solution to the contact map
        :param cl_list: a list of cluster ids to enable or None (all ordered clusters)
        :param ordered_only: include only clusters which have been ordered
        :param min_extent: include only clusters whose total extent is greater
        :param verbose: debug output
        :return: the filtered list of cluster ids in ascending numerical order
        """

        # start with all clusters if unspecified
        if cl_list is None:
            cl_list = clustering.keys()

        # drop clusters that are too small
        if min_extent:
            cl_list = [k for k in cl_list if clustering[k]['extent'] >= min_extent]

        if ordered_only:
            # drop any clusters that have not been ordered
            cl_list = [k for k in cl_list if 'order' in clustering[k]]

        # impose a consistent order
        cl_list = sorted(cl_list)

        if ordered_only:
            # use the determined order and orientation
            cmb_ord = np.hstack([clustering[k]['order'] for k in cl_list])
        else:
            # arbitrary order and orientation
            cmb_ord = np.hstack([SeqOrder.asindex(clustering[k]['seq_ids']) for k in cl_list])

        if len(cmb_ord) == 0:
            raise RuntimeError('no requested cluster contained ordering information')

        if verbose:
            print 'The cluster set contains {} sequences'.format(len(cmb_ord))

        # prepare the mask
        _mask = np.zeros_like(self.order.mask_vector(), dtype=np.bool)
        _mask[cmb_ord['index']] = True
        _mask &= self.get_primary_acceptance_mask()
        if verbose:
            print 'After masking ordering contains {} sequences'.format(_mask.sum())
        self.order.set_mask_only(_mask)
        self.order.set_order_and_orientation(cmb_ord, implicit_excl=True)

        return cl_list

    def plot_clusters(self, fname, clustering, cl_list=None, simple=True, permute=False, block_reduction=None,
                      ordered_only=True, min_extent=None, use_taxo=False, verbose=False, **kwargs):
        """
        Plot the contact map, annotating the map with cluster names and boundaries.

        For large contact maps, block reduction can be employed to reduce the size for plotting purposes. Using
        block_reduction=2 will reduce the map dimensions by a factor of 2. Must be integer.

        :param fname: output file name
        :param clustering: the cluster solution
        :param cl_list: the list of cluster ids to include in plot. If none, include all ordered clusters
        :param simple: True plot seq map, False plot the extent map
        :param permute: permute the map with the present order
        :param block_reduction: specify a reduction factor (>1 integer) to reduce map size for plotting
        :param ordered_only: include only clusters which have been ordered
        :param min_extent: include only clusters whose total extent is greater
        :param use_taxo: use taxonomic information within clustering, assuming it exists
        :param verbose: debug output
        :param kwargs: additional options passed to plot()
        """

        if verbose:
            if cl_list is None:
                print 'Plotting heatmap of complete solution'
            else:
                print 'Plotting heatmap of solution to clusters [{}]'.format(cl_list)

        # build a list of relevant clusters and setup the associated mask
        cl_list = self._enable_clusters(clustering,  cl_list=cl_list, ordered_only=ordered_only,
                                        min_extent=min_extent, verbose=verbose)

        if simple:
            # tick spacing simple the number of sequences in the cluster
            tick_locs = np.cumsum([0] + [len(clustering[k]['seq_ids']) for k in cl_list])
            if self.is_tipbased():
                tick_locs *= 2
        else:
            # tick spacing depends on cumulative bins for sequences in cluster
            # cumulative bin count, excluding masked sequences
            csbins = [0]
            for k in cl_list:
                # get the order records for the sequences in cluster k
                _oi = self.order.order[clustering[k]['seq_ids']]
                # count the cumulative bins at each cluster for those sequences which are not masked
                csbins.append(self.grouping.bins[clustering[k]['seq_ids'][_oi['mask']]].sum() + csbins[-1])
            tick_locs = np.array(csbins, dtype=np.int)

        if use_taxo:
            _labels = [clustering[cl_id]['taxon'] for cl_id in cl_list]
        else:
            _labels = [clustering[cl_id]['name'] for cl_id in cl_list]

        self.plot(fname, permute=permute, simple=simple, tick_locs=tick_locs, tick_labs=_labels,
                  block_reduction=block_reduction, verbose=verbose, **kwargs)

    def plot(self, fname, simple=False, tick_locs=None, tick_labs=None, norm=False, permute=False, pattern_only=False,
             dpi=180, width=25, height=22, zero_diag=False, alpha=0.01, robust=False, block_reduction=None,
             verbose=False):
        """
        Plot the contact map. This can either be as a sparse pattern (requiring much less memory but without visual
        cues about intensity), simple sequence or full binned map and normalized or permuted.

        :param fname: output file name
        :param tick_locs: major tick locations (minors take the midpoints)
        :param tick_labs: minor tick labels
        :param simple: if true, sequence only map plotted
        :param norm: normalize intensities by geometric mean of lengths
        :param permute: reorder map to current order
        :param pattern_only: plot only a sparse pattern (much lower memory requirements)
        :param dpi: adjust DPI of output
        :param width: plot width in inches
        :param height: plot height in inches
        :param zero_diag: set bright self-interactions to zero
        :param alpha: log intensities are log (x + alpha)
        :param robust: use seaborn robust dynamic range feature
        :param block_reduction: specify a reduction factor (>1 integer) to reduce map size for plotting
        :param verbose: debug output
        """

        plt.style.use('ggplot')

        fig = plt.figure()
        fig.set_figwidth(width)
        fig.set_figheight(height)
        ax = fig.add_subplot(111)

        if simple:
            # # assume the map must be remade
            if self.processed_map is None:
                self.prepare_seq_map(norm=norm)
            _map = self.get_processed_map(permute=permute, verbose=verbose)
        else:
            _map = self.get_extent_map(norm=norm, permute=permute)

        if pattern_only:
            if zero_diag:
                _map.setdiag(0)
            ax.spy(_map.tocsr(), markersize=5 if simple else 1)

        else:
            # a dense array is necessary here
            if block_reduction is not None:
                if verbose:
                    print 'Down-sampling map for plotting by a factor of: {}'.format(block_reduction)
                full_size = _map.shape[0]
                _map = simple_sparse.downsample(_map, block_reduction)
                tick_locs = np.floor(tick_locs.astype(np.float) / block_reduction)
                if verbose:
                    print 'Map reduced from {} to {}'.format(full_size, _map.shape)
            _map = _map.toarray()
            if zero_diag:
                if verbose:
                    print 'Removing diagonal'
                np.fill_diagonal(_map, 0)
            _map = np.log(_map + alpha)

            if verbose:
                print 'Making raster image'
            seaborn.heatmap(_map, robust=robust, square=True, linewidths=0, ax=ax, cbar=False)

        if tick_locs is not None:

            plt.tick_params(axis='both', which='both',
                            right=False, left=False, bottom=False, top=False,
                            labelright=False, labelleft=False, labelbottom=False, labeltop=False)

            if tick_labs is not None:
                min_labels = ticker.FixedFormatter(tick_labs)
                ax.tick_params(axis='y', which='minor', left=True, labelleft=True, labelsize=10)

                min_ticks = ticker.FixedLocator(tick_locs[:-1] + 0.5 * np.diff(tick_locs))

                ax.yaxis.set_minor_formatter(min_labels)
                ax.yaxis.set_minor_locator(min_ticks)

            # seaborn will not display the grid, so we make our own.
            ax.hlines(tick_locs, *ax.get_xlim(), color='grey', linewidth=0.5, linestyle='-.')
            ax.vlines(tick_locs, *ax.get_ylim(), color='grey', linewidth=0.5, linestyle='-.')

        if verbose:
            print 'Saving plot'
        fig.tight_layout()
        plt.savefig(fname, dpi=dpi)
        plt.close(fig)

    def to_graph(self, norm=True, bisto=False, scale=False, extern_ids=False,
                 min_len=None, min_sig=None, verbose=False):
        """
        Convert the seq_map to a undirected Networkx Graph.

        The contact map is effectively an adjacency matrix, where sequences
        are the nodes and edges weighted by the observed counts. Self-loops
        are not included by default and weights are affected by normalisation
        choices.

        :param norm: normalize weights by length
        :param bisto: normalise using bistochasticity
        :param scale: scale weights (max_w = 1)
        :param extern_ids: use the original external sequence identifiers for node ids
        :param min_len: override minimum sequence length, otherwise use instance's setting)
        :param min_sig: override minimum off-diagonal signal (in raw counts), otherwise use instance's setting)
        :param verbose: debug output
        :return: graph of contigs
        """
        if extern_ids:
            _nn = lambda x: self.seq_info[x].name
        else:
            _nn = lambda x: x

        if not min_len and not min_sig:
            # use the (potentially) existing mask if default criteria
            self.set_primary_acceptance_mask(verbose=verbose)
        else:
            # update the acceptance mask if user has specified new criteria
            self.set_primary_acceptance_mask(min_len, min_sig, update=True, verbose=verbose)

        if self.processed_map is None:
            self.prepare_seq_map(norm=norm, bisto=bisto, verbose=verbose)
        _map = self.get_subspace(marginalise=True, flatten=False, verbose=verbose)

        if verbose:
            print 'Graph will have {} nodes'.format(self.order.count_accepted())

        g = nx.Graph()

        if not sp.isspmatrix_coo(_map):
            _map = _map.tocoo()

        scl = 1.0/_map.max() if scale else 1

        if verbose:
            print 'Building graph from edges'

        import tqdm

        for u, v, w in tqdm.tqdm(itertools.izip(_map.row, _map.col, _map.data), desc='adding edges', total=_map.nnz):
            g.add_edge(_nn(u), _nn(v), weight=w * scl)

        if verbose:
            print nx.info(g)

        return g

    @staticmethod
    def _add_cluster_names(clustering, prefix='CL'):
        """
        Add sequential names beginning from 1 to a clustering in-place.

        A pedantic determination of how many digits are required for the
        largest cluster number is performed, so names will sort conveniently
        in alphanumeric order and align to the eye in output information.

        :param clustering: clustering solution returned from ContactMap.cluster_map
        :param prefix: static prefix of cluster names.
        """
        try:
            num_width = max(1, int(np.ceil(np.log10(max(clustering)+1))))
        except OverflowError:
            num_width = 1

        for cl_id in clustering:
            # names will 1-based
            clustering[cl_id]['name'] = '{0}{1:0{2}d}'.format(prefix, cl_id+1, num_width)

    def cluster_map(self, method='infomap', min_len=None, min_sig=None, work_dir='.', seed=None, verbose=False):
        """
        Cluster a contact map into groups, as an approximate proxy for "species" bins. This is recommended prior to
        applying TSP ordering, minimising the breaking of its assumptions that all nodes should connect and be
        traversed.

        :param method: clustering algorithm to employ
        :param min_len: override minimum sequence length, otherwise use instance's setting)
        :param min_sig: override minimum off-diagonal signal (in raw counts), otherwise use instance's setting)
        :param work_dir: working directory to which files are written during clustering
        :param seed: specify a random seed, otherwise one is generated at runtime.
        :param verbose: debug output
        :return: a dictionary detailing the full clustering of the contact map
        """

        def _read_mcl(pathname):
            """
            Read a MCL solution file converting this to a TruthTable.

            :param pathname: mcl file name
            :return: dict of cluster_id to array of seq_ids
            """

            with open(pathname, 'r') as h_in:
                # read the MCL file, which lists all members of a class on a single line
                # the class ids are implicit, therefore we use line number.
                cl_map = {}
                for cl_id, line in enumerate(h_in):
                    line = line.rstrip()
                    if not line:
                        break
                    cl_map[cl_id] = np.array(sorted([int(tok) for tok in line.split()]))
            return cl_map

        def _read_table(pathname, seq_col=0, cl_col=1):
            # type: (str, Optional[int], int) -> dict
            """
            Read cluster solution from a tabular file, one assignment per line. Implicit sequence
            naming is achieved by setting seq_col=None. The reverse (implicit column naming) is
            not currently supported.

            :param pathname: table file name
            :param seq_col: column number of seq_ids
            :param cl_col: column number of cluster ids
            :return: dict of cluster_id to array of seq_ids
            """
            assert seq_col != cl_col, 'sequence and cluster columns must be different'
            with open(pathname, 'r') as h_in:
                cl_map = {}
                n = 0
                for line in h_in:
                    line = line.strip()
                    if not line:
                        break
                    if seq_col is None:
                        cl_id = int(line)
                        seq_id = n
                        n += 1
                    else:
                        t = line.split()
                        if len(t) != 2:
                            print 'invalid line encountered when reading cluster table: [{}]'.format(line)

                        seq_id, cl_id = int(t[seq_col]), int(t[cl_col])
                    cl_map.setdefault(cl_id, []).append(seq_id)
                for k in cl_map:
                    cl_map[k] = np.array(cl_map[k], dtype=np.int)
                return cl_map

        def _read_tree(pathname):
            """
            Read a tree clustering file as output by Infomap.

            :param pathname: the path to the tree file
            :return: dict of cluster_id to array of seq_ids
            """
            with open(pathname, 'r') as in_h:
                cl_map = {}
                for line in in_h:
                    line = line.strip()
                    if not line:
                        break
                    if line.startswith('#'):
                        continue
                    fields = line.split()
                    hierarchy = fields[0].split(':')
                    # take everything in the cluster assignment string except for the final token,
                    # which is the object's id within the cluster.
                    cl_map.setdefault(tuple(['orig'] + hierarchy[:-1]), []).append(fields[-1])

                # rename clusters and order descending in size
                desc_key = sorted(cl_map, key=lambda x: len(cl_map[x]), reverse=True)
                for n, k in enumerate(desc_key):
                    cl_map[n] = np.array(cl_map.pop(k), dtype=np.int)

            return cl_map

        def _write_edges(g, parent_dir, base_name, sep=' '):
            """
            Prepare an edge-list file from the specified graph. This will be written to the
            specified parent directory, using basename.edges
            :param g: the graph
            :param parent_dir: parent directory of file
            :param base_name: file base name
            :param sep: separator within file
            :return: file name
            """
            edge_file = os.path.join(parent_dir, '{}.edges'.format(base_name))
            nx.write_edgelist(g, edge_file, data=['weight'], delimiter=sep)
            return edge_file

        assert os.path.exists(work_dir), 'supplied output path [{}] does not exist'.format(work_dir)

        if seed is None:
            seed = make_random_seed()

        base_name = 'cm_graph'
        g = self.to_graph(min_len=min_len, min_sig=min_sig, norm=True, bisto=True, scale=True, verbose=verbose)

        method = method.lower()
        if verbose:
            print 'Performing clustering method [{}]'.format(method)

        if method == 'louvain':
            cl_to_ids = louvain_cluster.cluster(g, no_iso=False, ragbag=False)
        elif method == 'mcl':
            with open(os.path.join(work_dir, 'mcl.log'), 'w+') as stdout:
                ofile = os.path.join(work_dir, '{}.mcl'.format(base_name))
                edge_file = _write_edges(g, work_dir, base_name)
                nx.write_edgelist(g, edge_file, data=['weight'])
                subprocess.check_call([ext_path('mcl'),  edge_file, '--abc', '-I', '1.2', '-o', ofile],
                                      stdout=stdout, stderr=subprocess.STDOUT)
                cl_to_ids = _read_mcl(ofile)
        elif method == 'simap':
            with open(os.path.join(work_dir, 'simap.log'), 'w+') as stdout:
                ofile = os.path.join(work_dir, '{}.simap'.format(base_name))
                edge_file = _write_edges(g, work_dir, base_name)
                subprocess.check_call(['java', '-jar', ext_path('simap-1.0.0.jar'), 'mdl', '-s', str(seed),
                                       '-i', '1e-5', '1e-3', '-a', '1e-5', '-g', edge_file, '-o', ofile],
                                      stdout=stdout, stderr=subprocess.STDOUT)
                cl_to_ids = _read_table(ofile)
        elif method == 'infomap':
            with open(os.path.join(work_dir, 'infomap.log'), 'w+') as stdout:
                edge_file = _write_edges(g, work_dir, base_name)
                subprocess.check_call([ext_path('Infomap'), '-u', '-v', '-z', '-i', 'link-list', '-s', str(seed),
                                       '-N', '10', edge_file, work_dir],
                                      stdout=stdout, stderr=subprocess.STDOUT)
                cl_to_ids = _read_tree(os.path.join(work_dir, '{}.tree'.format(base_name)))
        elif method == 'slm':
            with open(os.path.join(work_dir, 'slm.log'), 'w+') as stdout:
                mod_func = '1'
                resolution = '2.0'
                opti_algo = '3'
                n_starts = '10'
                n_iters = '10'
                ofile = os.path.join(work_dir, '{}.slm'.format(base_name))
                verb = '1'
                edge_file = _write_edges(g, work_dir, base_name, sep='\t')
                subprocess.check_call(['java', '-jar', ext_path('ModularityOptimizer.jar'), edge_file, ofile,
                                       mod_func, resolution, opti_algo, n_starts, n_iters, str(seed), verb],
                                      stdout=stdout, stderr=subprocess.STDOUT)
                cl_to_ids = _read_table(ofile, seq_col=None, cl_col=0)
        else:
            raise RuntimeError('unimplemented method [{}]'.format(method))

        if verbose:
            print 'Gathering clustering results'
            print 'There were {} clusters'.format(len(cl_to_ids))

        # standardise the results, where sequences in each cluster
        # are listed in ascending order
        clustering = {}
        for cl_id, _seqs in cl_to_ids.iteritems():
            _ord = SeqOrder.asindex(np.sort(_seqs))
            # IMPORTANT!! sequences are remapped to their gapless indices
            _seqs = self.order.remap_gapless(_ord)['index']

            clustering[cl_id] = {
                'seq_ids': _seqs,
                'extent': self.order.lengths()[_seqs].sum()
                # TODO add other details for clusters here
                # - residual modularity, permanence
            }

        # reestablish clusters in descending order of extent
        sorted_keys = sorted(clustering, key=lambda k: clustering[k]['extent'], reverse=True)
        clustering = {n: clustering[k] for n, k in enumerate(sorted_keys)}

        ContactMap._add_cluster_names(clustering)

        return clustering

    def cluster_report(self, clustering, source_fasta=None, is_spades=True, verbose=False):
        """
        For each cluster, analyze the member sequences and build a report.
        Update the clustering dictionary with this result by adding a "report" for each.

        :param clustering: clustering solution dictionary
        :param source_fasta: source assembly fasta (other than defined at instantiation)
        :param is_spades: if SPAdes output, we can extract coverage information from the sequence names
        :param verbose: debug output
        """
        if verbose:
            print 'Analyzing cluster member sequences'

        seq_info = self.seq_info

        if source_fasta is None:
            source_fasta = self.seq_file

        import Bio.SeqUtils as SeqUtils

        # set up indexed access to the input fasta
        with contextlib.closing(IndexedFasta(source_fasta)) as seq_db:
            # iterate over the cluster set, in the existing order
            for cl_id, cl_info in clustering.iteritems():
                _len = []
                _cov = []
                _gc = []
                for n, _seq_id in enumerate(np.sort(cl_info['seq_ids']), 1):
                    # get the sequence's external name and length
                    _name = seq_info[_seq_id].name
                    _len.append(seq_info[_seq_id].length)
                    # fetch the SeqRecord object from the input fasta
                    _seq = seq_db[_name]
                    _gc.append(SeqUtils.GC(_seq.seq))
                    if is_spades:
                        _cov.append(float(_name.split('_')[-1]))

                if is_spades:
                    report = np.array(zip(_len, _gc, _cov),
                                      dtype=[('length', np.int),
                                             ('gc', np.float),
                                             ('cov', np.float)])
                else:
                    report = np.array(zip(_len, _gc, _cov),
                                      dtype=[('length', np.int),
                                             ('gc', np.float)])
                clustering[cl_id]['report'] = report

    def order_clusters(self, clustering, min_len=None, min_sig=None, max_fold=None, min_extent=None, min_size=1,
                       work_dir='.', seed=None, dist_method='neglog', verbose=False):
        """
        Determine the order of sequences for a given clustering solution, as returned by cluster_map.

        :param clustering: the full clustering solution, derived from the supplied contact map
        :param min_len: within a cluster exclude sequences that are too short (bp)
        :param min_sig: within a cluster exclude sequences with weak signal (counts)
        :param max_fold: within a cluster, exclude sequences that appear to be overly represented
        :param min_size: skip clusters which containt too few sequences
        :param min_extent: skip clusters whose total extent (bp) is too short
        :param work_dir: working directory
        :param seed: random seed
        :param dist_method: method used to transform contact map to a distance matrix
        :param verbose: debug output
        :return: map of cluster orders, by cluster id
        """
        assert os.path.exists(work_dir), 'supplied output path [{}] does not exist'.format(work_dir)

        if verbose:
            print 'Ordering clusters'

        if self.processed_map is None:
            self.set_primary_acceptance_mask(min_len, min_sig, max_fold=max_fold, update=True, verbose=verbose)
            self.prepare_seq_map(norm=True, bisto=True, mean_type='geometric', verbose=verbose)

        for cl_id, cl_info in clustering.iteritems():

            cl_size = len(cl_info['seq_ids'])

            if cl_info['extent'] < min_extent:
                if verbose:
                    print 'Excluding {} too little extent: {} bp'.format(cl_info['name'], cl_info['extent'])
                continue
            elif cl_size < min_size:
                if verbose:
                    print 'Excluding {} too few sequences: {} '.format(cl_info['name'], cl_size)
                continue
            elif verbose:
                print 'Ordering {} extent: {} size: {}'.format(cl_info['name'], cl_info['extent'], cl_size)

            try:
                # we'll consider only sequences in the cluster
                _mask = np.zeros_like(self.order.mask_vector())
                _mask[cl_info['seq_ids']] = True

                _map = self.get_subspace(external_mask=_mask, verbose=verbose)

                print 'Cluster size: {} ordering map size: {}'.format(cl_size, _map.shape)

                _ord = self.find_order(_map, work_dir=work_dir, inverse_method=dist_method, seed=seed, verbose=verbose)

            except NoneAcceptedException as e:
                print '{} : cluster {} will be masked'.format(e.message, cl_info['name'])
                continue
            except TooFewException as e:
                print '{} : ordering not possible for cluster {}'.format(e.message, cl_info['name'])
                _ord = self.order.accepted_order()

            clustering[cl_id]['order'] = _ord

        return clustering

    def write_report(self, fname, clustering):
        """
        Create a tabular report of each cluster from a clustering report. Write the table to CSV.

        :param fname: the CSV output file name
        :param clustering: the input clustering, which contains a report
        """

        def _expect(w, x):
            """
            Weighted expectation of x with weights w. Weights do not need to be
            normalised

            :param w: weights
            :param x: variable
            :return: expectation value of x
            """
            wsum = float(w.sum())
            return np.sum(w * x) / wsum

        df = []
        for k, v in clustering.iteritems():
            try:
                sr = v['report']

                df.append([k,
                           v['name'],
                           len(v['seq_ids']),
                           v['extent'],
                           _expect(sr['length'], sr['gc']), sr['gc'].mean(), np.median(sr['gc']), sr['gc'].std(),
                           _expect(sr['length'], sr['cov']), sr['cov'].mean(), np.median(sr['cov']), sr['cov'].std()])

            except KeyError:
                raise NoReportException(k)

        df = pandas.DataFrame(df, columns=['id', 'name', 'size', 'extent',
                                           'gc_expect', 'gc_mean', 'gc_median', 'gc_std',
                                           'cov_expect', 'cov_mean', 'cov_median', 'cov_std'])
        df.set_index('id', inplace=True)
        df.to_csv(fname, sep=',')

    def write_fasta(self, clustering, output_dir, source_fasta=None, clobber=False, verbose=False):
        """
        Write out multi-fasta for all determined clusters in clustering.

        For each cluster, sequence order and orientation is as follows.
        1. for unordered clusters, sequences will be in descending nucleotide length and
           in original input orientation.
        2. for ordered clusters, sequences will appear in the prescribed order and
           orientation.

        :param clustering: the clustering result, possibly also ordered
        :param output_dir: parent output path
        :param source_fasta: specify a source fasta file, otherwise assume the same path as was used in parsing
        :param clobber: True overwrite files in the output path. Does not remove directories
        :param verbose: debug output
        """

        make_dir(output_dir)

        if verbose:
            print 'Writing output to the path [{}]'.format(output_dir)

        seq_info = self.seq_info

        parent_dir = os.path.join(output_dir, 'fasta')
        make_dir(parent_dir)

        if source_fasta is None:
            source_fasta = self.seq_file

        # set up indexed access to the input fasta
        with contextlib.closing(IndexedFasta(source_fasta)) as seq_db:

            # iterate over the cluster set, in the existing order
            for cl_id, cl_info in clustering.iteritems():

                # Each cluster produces a multi-fasta. Sequences are not joined
                cl_path = os.path.join(parent_dir, '{}.fna'.format(cl_info['name']))

                if not clobber and os.path.exists(cl_path):
                    raise IOError('Output path exists [{}] and overwriting not enabled'.format(cl_path))

                # determine the number of digits required for cluster sequence names
                try:
                    num_width = max(1, int(np.ceil(np.log10(len(cl_info['seq_ids'])+1))))
                except OverflowError:
                    num_width = 1

                with open(cl_path, 'w') as output_h:

                    if verbose:
                        print 'Writing full unordered FASTA for cluster {} to {}'.format(cl_id, cl_path)

                    # iterate simply over sequence ids, while imposing ascending numerical order
                    for n, _seq_id in enumerate(np.sort(cl_info['seq_ids']), 1):

                        # get the sequence's external name and length
                        _name = seq_info[_seq_id].name
                        _length = seq_info[_seq_id].length
                        # fetch the SeqRecord object from the input fasta
                        _seq = seq_db[_name]
                        # orientations are listed as unknown
                        _ori_symb = 'UNKNOWN'

                        # add a new name and description
                        _seq.id = '{0}_{1:0{2}d}'.format(cl_info['name'], n, num_width)
                        _seq.name = _seq.id
                        _seq.description = 'contig:{} ori:{} length:{}'.format(_name, _ori_symb, _length)
                        SeqIO.write(_seq, output_h, 'fasta')

                # write a separate ordered fasta as this is often a subset of all sequences
                if 'order' in cl_info:

                    # Each cluster produces a multi-fasta. Sequences are not joined
                    cl_path = os.path.join(parent_dir, '{}.ordered.fna'.format(cl_info['name']))

                    if not clobber and os.path.exists(cl_path):
                        raise IOError('Output path exists [{}] and overwriting not enabled'.format(cl_path))

                    with open(cl_path, 'w') as output_h:

                        if verbose:
                            print 'Writing ordered FASTA for cluster {} to {}'.format(cl_id, cl_path)

                        # iterate over cluster members, in the determined order
                        for n, _oi in enumerate(cl_info['order'], 1):

                            # get the sequence's external name and length
                            _name = seq_info[_oi['index']].name
                            _length = seq_info[_oi['index']].length
                            # fetch the SeqRecord object from the input fasta
                            _seq = seq_db[_name]

                            # reverse complement as needed
                            if _oi['ori'] == SeqOrder.REVERSE:
                                _seq = _seq.reverse_complement()
                                _ori_symb = '-'
                            elif _oi['ori'] == SeqOrder.FORWARD:
                                _ori_symb = '+'
                            else:
                                raise UnknownOrientationStateException(_oi['ori'])

                            # add a new name and description
                            _seq.id = '{0}_{1:0{2}d}'.format(cl_info['name'], n, num_width)
                            _seq.name = _seq.id
                            _seq.description = 'contig:{} ori:{} length:{}'.format(_name, _ori_symb, _length)
                            SeqIO.write(_seq, output_h, 'fasta')


class SequenceAnalyzer:

    COV_TYPE = np.dtype([('index', np.int16), ('status', np.bool), ('node', np.float),
                         ('local', np.float), ('fold', np.float)])

    @staticmethod
    def read_report(file_name):
        return yaml.load(open(file_name, 'r'))

    def __init__(self, seq_map, seq_report, seq_info, tip_size):
        self.seq_map = seq_map
        self.seq_report = seq_report
        self.seq_info = seq_info
        self.tip_size = tip_size

    def _contact_graph(self):
        g = nx.Graph()
        n_seq = len(self.seq_info)

        for i in xrange(n_seq):
            si = self.seq_info[i]
            d = self.seq_report['seq_info'][si.name]
            if self.tip_size:
                g.add_node(i, _id=si.name,
                           _cov=float(d['coverage']),
                           # this is a 2 element list for tip mapping
                           _sites=d['sites'],
                           _len=int(d['length']))
            else:
                g.add_node(i, _id=si.name,
                           _cov=float(d['coverage']),
                           _sites=int(d['sites']),
                           _len=int(d['length']))

        _m = self.seq_map
        if self.tip_size:
            _m = _m.sum(axis=(2, 3))

        for i in xrange(n_seq):
            for j in xrange(i, n_seq):
                if _m[i, j] > 0:
                    g.add_edge(i, j, weight=float(_m[i, j]))

        return g

    @staticmethod
    def _nlargest(g, u, n, k=0, local_set=None):
        """
        Build a list of nodes of length n within a radius k hops of node u in graph g, which have the
        largest weight. For Hi-C data, after normalisation, high weight can be used as a means of inferring
        proximity.

        :param g: the graph which to analyse
        :param u: the target node
        :param n: the length of the 'nearby nodes' list
        :param k: the maximum number of hops away from u
        :param local_set: used in recursion
        :return: a set of nodes.
        """
        if not local_set:
            local_set = set()

        neighbors = [v[0] for v in heapq.nlargest(n+1, g[u].items(), key=lambda x: x[1]['weight'])]
        local_set.update(neighbors)
        if k > 0:
            for v in neighbors:
                if v == u:
                    continue
                SequenceAnalyzer._nlargest(g, v, n, k-1, local_set)

        return list(local_set)

    def report_degenerates(self, fold_max, min_len=0, verbose=False):
        """
        Making the assumption that degenerate sequences (those sequences which are repeats) have high coverage
        relative to their local region, report those nodes in the graph whose coverage exceeds a threshold.

        :param fold_max: the maximum relative coverage allowed (between a node and its local region)
        :param min_len: the shorest allowable sequence to consider
        :param verbose: debug output
        :return: a report of all sequences degenerate status
        """

        g = self._contact_graph()

        degens = []
        for u in g.nodes_iter():
            if g.node[u]['_len'] < min_len or g.degree(u) == 0:
                continue

            signif_local_nodes = SequenceAnalyzer._nlargest(g, u, 4, 1)
            local_mean_cov = mstats.gmean(np.array([g.node[v]['_cov'] for v in signif_local_nodes]))
            fold_vs_local = g.node[u]['_cov'] / float(local_mean_cov)

            is_degen = True if fold_vs_local > fold_max else False

            degens.append((u, is_degen, g.node[u]['_cov'], local_mean_cov, fold_vs_local))

        degens = np.array(degens, dtype=SequenceAnalyzer.COV_TYPE)

        if verbose:
            if len(degens) == 0:
                print 'No degenerate sequences found'
            else:
                print '\tDegenerate sequence report\n\t--------------------------'
                for di in degens[degens['status']]:
                    print '\t', di

        return degens


def save_object(file_name, obj):
    """
    Serialize an object to a file with gzip compression. .gz will automatically be
    added if missing.

    :param file_name: output file name
    :param obj: object to serialize
    """
    with io_utils.open_output(file_name, compress='gzip') as out_h:
        cPickle.dump(obj, out_h)


def load_object(file_name):
    """
    Deserialize an object from a file with automatic support for compression.

    :param file_name: input file name
    :return: deserialzied object
    """
    with io_utils.open_input(file_name) as in_h:
        return cPickle.load(in_h)


if __name__ == '__main__':
    import argparse


    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    parser = argparse.ArgumentParser(description='Create a 3C fragment map from a BAM file')

    parser.add_argument('-s', '--seed', default=None, help='Random seed')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    parser.add_argument('-f', '--format', choices=['csv', 'h5'], default='csv',
                        help='Input contact map format')
    parser.add_argument('--eta', default=False, action='store_true', help='Precount bam alignments to provide an ETA')
    parser.add_argument('--med-alpha', type=int, default=10, help='Coverage median filter factor.')
    parser.add_argument('--tip-size', type=int, default=None,
                        help='Accept only pairs which map within reference tips (size in bp).')
    parser.add_argument('--strong', type=int, default=None,
                        help='Using strong matching constraint (minimum matches in alignments).')
    parser.add_argument('--bin-size', type=int, required=False, help='Size of bins in bp')
    parser.add_argument('--min-insert', type=int, required=False, help='Minimum pair separation')
    parser.add_argument('--min-mapq', type=int, default=0, help='Minimum acceptable mapping quality [0]')
    parser.add_argument('--min-reflen', type=int, default=1, help='Minimum acceptable reference length [0]')
    parser.add_argument('--min-signal', type=int, default=1, help='Minimum acceptable trans signal [1]')
    parser.add_argument('--max-image', type=int, default=4000, help='Maximum image size for plots [4000]')
    parser.add_argument('--min-order-size', type=int, default=5, help='Minimum cluster size for ordering [5]')
    parser.add_argument('--min-order-extent', type=int, default=50000,
                        help='Minimum cluster extent (kb) for ordering [50000]')
    parser.add_argument('--dist-method', choices=['inverse', 'neglog'], default='inverse',
                        help='Distance method for ordering [inverse]')
    parser.add_argument('--pickle', help='Picked contact map')
    parser.add_argument('-e', '--enzymes', required=True, action='append',
                        help='Case-sensitive enzyme name (NEB), use multiple times for multiple enzymes')
    parser.add_argument('--min-scflen', default=2000,
                        help='Minimum length of sequence to use in scaffolding [2000]')
    parser.add_argument('--skip-scaffolding', default=False, action='store_true',
                        help='Do not attempt to scaffold clusters')
    parser.add_argument('fasta', help='Reference fasta sequence')
    parser.add_argument('bam', help='Input bam file in query order')
    parser.add_argument('out_dir', help='Output directory')

    args = parser.parse_args()

    if not args.seed:
        args.seed = make_random_seed()
        logger.info('Generated random seed: {}'.format(args.seed))
    else:
        logger.info("User set random seed: {}".format(args.seed))

    make_dir(args.out_dir)

    if args.pickle:
        # Load a pre-existing serialized contact map
        logger.info('Loading existing contact map from: {}'.format(args.pickle))
        cm = load_object(args.pickle)
    else:
        # Create a contact map for analysis
        cm = ContactMap(args.bam,
                        args.enzymes,
                        args.fasta,
                        args.min_insert,
                        args.min_mapq,
                        min_len=args.min_reflen,
                        min_sig=args.min_signal,
                        # max_fold=args.max_fold,
                        strong=args.strong,
                        bin_size=args.bin_size,
                        tip_size=args.tip_size,
                        precount=args.eta,
                        med_alpha=args.med_alpha)

        if cm.is_empty():
            logger.info('Stopping as the map is empty')
            sys.exit(1)

        logger.info('Saving contact map instance...')
        save_object(os.path.join(args.out_dir, 'contact_map.p'), cm)

    # cluster the entire map
    clustering = cm.cluster_map(method='infomap', seed=args.seed, work_dir=args.out_dir, verbose=args.verbose)
    # generate report per cluster
    cm.cluster_report(clustering, is_spades=True, verbose=args.verbose)
    # serialize clustering
    save_object(os.path.join(args.out_dir, 'clustering.p'), clustering)
    # write a tabular report
    cm.write_report(os.path.join(args.out_dir, 'cluster_report.csv'), clustering)

    if not args.skip_scaffolding:
        # order each cluster
        cm.order_clusters(clustering, min_reflen=args.min_scflen, min_size=args.min_order_size,
                          min_extent=args.min_order_extent, dist_method=args.dist_method,
                          work_dir=args.out_dir, verbose=args.verbose)
        # save the ordered clustering to another file.
        # TODO remove serialization of highly redundant objects
        save_object(os.path.join(args.out_dir, 'clustering_ordered.p'), clustering)

    # write per-cluster fasta files
    cm.write_fasta(clustering, args.out_dir, verbose=args.verbose)

    # plot a heatmap, while making sure we don't exceed an maximum pixel size.
    image_size = cm.get_primary_acceptance_mask().sum()
    reduce_factor = None
    if image_size > args.max_image:
        from math import ceil
        reduce_factor = int(ceil(image_size / float(args.max_image)))
        logger.info('Reducing image size of {} by {}'.format(image_size, reduce_factor))

    cm.plot_clusters(os.path.join(args.out_dir, 'cluster_plot.png'), clustering,
                     block_reduction=reduce_factor, simple=False, permute=True)
