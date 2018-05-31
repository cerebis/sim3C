import cPickle
import heapq
import itertools
import logging
import os
import subprocess
import tempfile
import uuid
from collections import OrderedDict
from collections import namedtuple
from functools import partial

import Bio.SeqIO as SeqIO
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import networkx as nx
import numpy as np
import pysam
import scipy.sparse as sp
import scipy.stats.mstats as mstats
import seaborn
import yaml
from numba import jit, vectorize, int64, int32, float64, void
from numpy import log, pi

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


SeqInfo = namedtuple('SeqInfo', ['offset', 'refid', 'name', 'length'])


class NoneAcceptedException(Exception):
    """All sequences were excluded during filtering"""
    def __init__(self):
        super(NoneAcceptedException, self).__init__('all sequences were excluded')


class TooFewException(Exception):
    """LKH requires a minimum of three nodes"""
    def __init__(self):
        super(TooFewException, self).__init__('LKH requires a minimum of three sequences')


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

    :param an_order: ordered 1d container of identifiers.
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


class ExtentGrouping:

    def __init__(self, seq_info, bin_size):
        self.bins = []
        self.bin_size = bin_size
        self.map = []
        self.borders = []
        self.centers = []
        self.total_bins = 0

        for n, seq in enumerate(seq_info):

            if seq.length == 0:
                raise RuntimeError('Zero length sequence for {}'.format(seq))

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
            logger.info('{}: {} bins'.format(n, num_bins))

        self.bins = np.array(self.bins)


class SeqOrder:

    FORWARD = 1
    REVERSE = -1

    ACCEPTED = True
    EXCLUDED = False

    STRUCT_TYPE = np.dtype([('pos', np.int32), ('ori', np.int8), ('mask', np.bool), ('length', np.int32)])
    INDEX_TYPE = np.dtype([('index', np.int16), ('ori', np.int8)])

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

    def _update_positions(self):
        """
        An optimisation, whenever the positional state changes, this method must be called to
        maintain the current state in a separate array. This avoids unncessary redetermination of positions
        when calling .intervening() and .before() methods.
        """
        self._positions = np.argsort(self.order['pos'])

    def remap_gapless(self, gapless_indices):
        """
        Recover potentially sparse indices from a dense (gapless) index table. This can happen for
        algorithms which were supplied only filtered data.

        :param gapless_indices: dense index table (index, ori), without gaps from masking
        :return: remappped indices with gaps
        """
        shift = []
        n = 0
        for oi in np.sort(gapless_indices['index']):
            while not self.order[oi + n]['mask']:
                n += 1
            shift.append(n)

        remapped = []
        for oi in gapless_indices:
            remapped.append((oi['index'] + shift[oi['index']], oi['ori']))

        return np.array(remapped, dtype=SeqOrder.INDEX_TYPE)

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
        _p = self._positions[:self.count_accepted()].copy()
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
        # reorder, with masked sequences coming last, otherwise by position
        sorted_indices = np.lexsort([self.order['pos'], ~self.order['mask']])
        for n, i in enumerate(sorted_indices):
            self.order[i]['pos'] = n

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
        # augment the order to incldue default orientations
        _ord = np.array(zip(_ord, np.ones_like(_ord, dtype=np.bool)), dtype=SeqOrder.INDEX_TYPE)
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
        # set up a full mask, reflecting the request
        _mask = np.ones(len(self.order), dtype=np.bool)
        _mask[_id] = False

        # remask the set
        self.set_mask_only(_mask)

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


@jit(int64(int64[:, :], int64[:], int64, int64), nopython=True)
def intervening(_ord, _lengths, a, b):

    ix1 = None
    for i in xrange(len(_ord)):
        if a == _ord[i, 0]:
            ix1 = i
            break

    ix2 = None
    for i in xrange(len(_ord)):
        if b == _ord[i, 0]:
            ix2 = i
            break

    if ix1 == ix2:
        return 0

    if ix1 > ix2:
        ix1, ix2 = ix2, ix1

    return np.sum(_lengths[_ord[ix1 + 1:ix2, 0]])


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


@jit(void(float64[:,:,:,:], float64[:,:], float64[:], float64))
def fast_norm_seq(coords, data, tip_lengths, tip_size):
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


class ContactMap:

    def __init__(self, bam_file, cov_file, seq_file, min_insert, min_mapq=0, min_len=0, min_sig=1, max_fold=None,
                 random_seed=None, strong=None, bin_size=None, tip_size=None, precount=False, med_alpha=10):

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
        self.seq_report = SequenceAnalyzer.read_report(cov_file)
        self.seq_analyzer = None

        with pysam.AlignmentFile(bam_file, 'rb') as bam:

            # test that BAM file is the corret sort order
            if 'SO' not in bam.header['HD'] or bam.header['HD']['SO'] != 'queryname':
                raise IOError('BAM file must be sorted by read name')

            seq_db = None
            tmp_file = '.{}.seqdb'.format(str(uuid.uuid4()))
            try:
                # sqlite3 seems to hate tempfile based temporary files
                # therefore we just make our own
                seq_db = SeqIO.index_db(tmp_file, seq_file, 'fasta')
                # determine the set of active sequences
                # where the first filtration step is by length
                ref_count = {'seq_missing': 0, 'too_short': 0}
                offset = 0
                logger.info('Reading sequences...')
                for n, li in enumerate(bam.lengths):

                    # minimum length threshold
                    if li < min_len:
                        ref_count['too_short'] += 1
                        continue

                    seq_id = bam.references[n]
                    if seq_id not in seq_db:
                        logger.info('Sequence: "{}" not found in reference FASTA'.format(seq_id))
                        ref_count['seq_missing'] += 1
                        continue

                    seq = seq_db[seq_id].seq
                    assert len(seq) == li, 'Sequence lengths in {} do not agree: ' \
                                           'bam {} fasta {}'.format(seq_id, len(seq), li)

                    self.seq_info.append(SeqInfo(offset, n, seq_id, li))

                    logger.debug('Found {} bp sequence for {}'.format(li, seq_id))

                    offset += li

            finally:
                if seq_db:
                    seq_db.close()
                if os.path.exists(tmp_file):
                    os.unlink(tmp_file)

            # total extent covered
            self.total_len = offset
            self.total_seq = len(self.seq_info)

            # all sequences begin as active in mask
            # when filtered, mask[i] = False
            self.current_mask = np.ones(self.total_seq, dtype=np.bool)

            if self.total_seq == 0:
                logger.info('No sequences in BAM found in FASTA')
                raise RuntimeError('No sequences in BAM found in FASTA')

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

            # logger.info('Parsing depth file...')
            # self.cov_info = self.read_covfile(cov_file)
            # for k, ci in self.cov_info.iteritems():
            #     logger.info('RefId: {} depths A/L/R: {}, {}, {}'.format(k, ci['med'], ci['lmed'], ci['rmed']))

            # accumulate
            self._bin_map(bam)

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

        # logger.info('Preparing no-go intervals for median filtering...')
        # no_go = self.depth_intervals(self.med_alpha, 10)
        #
        # for k in no_go:
        #     print k
        #     for i in no_go[k]:
        #         print '\t',i
        #     print

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

    def save(self, fname):
        with io_utils.open_output(fname, compress='gzip') as out_h:
            cPickle.dump(self, out_h)

    @staticmethod
    def load(fname):
        with io_utils.open_input(fname) as in_h:
            return cPickle.load(in_h)

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

    def find_order(self, min_len, min_sig, max_fold, norm=True, bisto=True, mean_type='geometric',
                   inverse_method='inverse', external_mask=None, verbose=False):
        """
        Using LKH TSP solver, find the best ordering of the sequence map in terms of proximity ligation counts.
        Here, it is assumed that sequence proximity can be inferred from number of observed trans read-pairs, where
        an inverse relationship exists.

        :param min_len: the minimum sequence length to include
        :param min_sig: the minimum off-diagonal signal (trans read-pair count)
        :param max_fold: maximum locally-measured fold-coverage to permit
        :param norm: apply length normalisation before analysis
        :param bisto: make the matrix bistochastic before analysis
        :param mean_type: for length normalization, the chosen mean used between sequences
        :param inverse_method: the chosen inverse method for converting count (similarity) to distance
        :param external_mask: include (union with filter) an additional mask on sequences
        :param verbose: debug output
        :return: the surrogate ids in optimal order
        """

        # prepare the input sequence map for analysis
        self.prepare_seq_map(min_len, min_sig, max_fold, norm=norm, bisto=bisto, mean_type=mean_type,
                             external_mask=external_mask, verbose=verbose)

        # we'll supply a partially initialized distance function
        dist_func = partial(ordering.similarity_to_distance, method=inverse_method, alpha=1.2, beta=1,
                            verbose=verbose)

        if self.is_tipbased():
            m = self.get_processed_map()

            # a minimum of three sequences is required to run LKH
            if m.shape[0] < 3:
                raise TooFewException()

            lkh_o = ordering.lkh_order(m, 'lkhdb_run', lkh_exe='../lkh/LKH-3.0/LKH', precision=1, seed=12345, runs=2,
                                       pop_size=50, dist_func=dist_func, special=False, verbose=True,
                                       fixed_edges=[(i, i+1) for i in xrange(1, m.shape[0], 2)])

            # To solve this with TSP, doublet tips use a graph transformation, where each node comes a pair. Pairs
            # possess fixed inter-connecting edges which must be included in any solution tour.
            # Eg. node 0 -> (0,1) or node 1 -> (2,3). The fixed paths are undirected, depending on which direction is
            # traversed, defines the orientation of the sequence.

            # 1.  pair adjacent nodes by reshape the 1D array into a two-column array of half the length
            lkh_o = lkh_o.reshape(lkh_o.shape[0]/2, 2)
            # 2. convert to surrogate ids and infer orientation from paths taken through doublets.
            #   0->1 forward (+1): 1->0 reverse (-1).
            lkh_o = np.fromiter(((oi[0]/2, oi[1]-oi[0]) for oi in lkh_o), dtype=SeqOrder.INDEX_TYPE)

        else:
            m = self.get_processed_map()
            lkh_o = ordering.lkh_order(m, 'lkh_run', lkh_exe='../lkh/LKH-3.0/LKH', precision=1, seed=12345, runs=2,
                                       pop_size=50, dist_func=dist_func, special=False, verbose=True)

            # for singlet tours, no orientation can be inferred.
            lkh_o = np.fromiter(((oi, 1) for oi in lkh_o), dtype=SeqOrder.INDEX_TYPE)

        if verbose:
            print lkh_o

        # lkh ordering references the supplied matrix indices, not the surrogate ids.
        # we must map this consecutive set to the contact map indices.
        lkh_o = self.order.remap_gapless(lkh_o)

        return lkh_o

    def get_primary_acceptance_mask(self):
        return self.primary_acceptance_mask.copy()

    def set_primary_acceptance_mask(self, min_len, min_sig, max_fold=None, update=False, verbose=False):
        """
        Determine and set the filter mask using the specified constraints across the entire
        contact map. The mask is True when a sequence is considered acceptable wrt to the
        constraints. The mask is also returned by the function for convenience.

        :param min_len: minimum sequence length
        :param min_sig: minimum off-diagonal signal (counts)
        :param max_fold: maximum locally-measured fold-coverage to permit
        :param update: replace the current primary mask if it exists
        :param verbose: debug output
        :return: an acceptance mask over the entire contact map
        """
        # simply return the current mask if it has already been determined
        # and an update is not requested
        if not update and self.primary_acceptance_mask is not None:
            if verbose:
                print '\tusing existing mask'
            return self.get_primary_acceptance_mask()

        acceptance_mask = np.ones(self.total_seq, dtype=np.bool)

        # mask for sequences shorter than limit
        _mask = self.order.lengths() >= min_len
        if verbose:
            print '\tmin_len removed', self.total_seq - _mask.sum()
        acceptance_mask &= _mask

        # mask for sequences weaker than limit
        if self.is_tipbased():
            signal = simple_sparse.max_offdiag_4d(self.seq_map)
        else:
            signal = simple_sparse.max_offdiag(self.seq_map)
        _mask = signal >= min_sig
        if verbose:
            print '\tmin_sig removed', self.total_seq - _mask.sum()
        acceptance_mask &= _mask

        # mask for sequences considered to have aberrant coverage.
        if max_fold:
            seq_analyzer = SequenceAnalyzer(self.seq_map, self.seq_report, self.seq_info, self.tip_size)
            degen_seqs = seq_analyzer.report_degenerates(max_fold, verbose=False)
            _mask = ~degen_seqs['status']
            if verbose:
                print '\tdegen removed', self.total_seq - _mask.sum()
            acceptance_mask &= _mask

        # retain the union of all masks.
        self.primary_acceptance_mask = acceptance_mask

        if verbose:
            print '\tAccepted sequences (mask union)', self.order.count_accepted()

        return self.get_primary_acceptance_mask()

    def prepare_seq_map(self, min_len, min_sig, max_fold=None, norm=True, bisto=False, mean_type='geometric',
                        external_mask=None, update_mask=False, verbose=False):
        """
        Prepare the sequence map (seq_map) by application of various filters and normalisations.

        :param min_len: minimum sequence length
        :param min_sig: minimum inter-sequence signal in raw counts.
        :param max_fold: maximum locally-measured fold-coverage to permit
        :param norm: normalisation by sequence lengths
        :param bisto: make the output matrix bistochastic
        :param mean_type: when performing normalisation, use "geometric, harmonic or arithmetic" mean.
        :param external_mask: include (union with filter) an additional mask on sequences
        :param update_mask: replace the existing primary mask, if one alread exists
        :param verbose: debug output
        """
        if verbose:
            print 'Prepare sequence map'
            print '\toriginally', self.seq_map.shape

        m = self.seq_map.astype(np.float)

        # apply length normalisation if requested
        if norm:
            m = self._norm_seq(m, self.is_tipbased(), mean_type=mean_type)
            if verbose:
                print '\tnormed'

        # assign a mask to avoid problematic sequences
        _mask = self.set_primary_acceptance_mask(min_len, min_sig, max_fold=max_fold,
                                                 update=update_mask, verbose=verbose)
        # from a union of the sequence filter and external mask
        if external_mask is not None:
            _mask &= external_mask
            print '\tunion with external', _mask.sum()
        self.order.set_mask_only(_mask)

        if self.order.count_accepted() < 1:
            raise NoneAcceptedException()

        # if there are sequences to mask, remove them from the map
        if self.order.count_accepted() < self.total_seq:
            if self.is_tipbased():
                m = simple_sparse.compress_4d(m, self.order.mask_vector())
            else:
                m = simple_sparse.compress(m.tocoo(), self.order.mask_vector())
            if verbose:
                print '\tfilter reduced', m.shape

        # make map bistochastic if requested
        if bisto:
            m, scl = self._bisto_seq(m, self.is_tipbased())
            # retain the scale factors
            self.bisto_scale = scl
            if verbose:
                print '\tbalanced'

        if self.is_tipbased():
            # lastly, convert the 4D map into a 2Nx2N 2D map.
            m = simple_sparse.flatten_tensor_4d(m)

        self.processed_map = m

    def get_processed_map(self, permute=False, dtype=np.float, verbose=False):
        """
        Return a copy of the processed map

        :param permute: reorder with current state
        :param dtype: cast to another type
        :param verbose: debug output
        :return: the processed map
        """
        _m = self.processed_map.astype(dtype)
        # reorder if requested
        if permute:
            _m = self._reorder_seq(_m, self.is_tipbased())
            if verbose:
                print '\treordered'
        return _m

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
            print '\toriginally', self.extent_map.shape

        m = self.extent_map.astype(np.float)

        # apply length normalisation if requested
        if norm:
            m = self._norm_extent(m, mean_type)
            if verbose:
                print '\tnormed'

        # if there are sequences to mask, remove them from the map
        if self.order.count_accepted() < self.total_seq:
            m = self._compress_extent(m)
            if verbose:
                print '\tfilter reduced', m.shape

        # make map bistochastic if requested
        if bisto:
            m, scl = simple_sparse.kr_biostochastic(m)
            if verbose:
                print '\tbalanced'

        # reorder using current order state
        if permute:
            m = self._reorder_extent(m)
            if verbose:
                print '\treordered'

        return m

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

    def _bisto_seq(self, m, tip_based):
        """
        Make a contact map bistochastic. This is another form of normslisation. Automatically
        handles 2D and 4D maps.
        :param m: a map to balance (make bistochastic)
        :param tip_based: treat the supplied map as a tip-based tensor
        :return: the balanced map
        """
        if tip_based:
            m, scl = simple_sparse.kr_biostochastic_4d(m)
        else:
            m, scl = simple_sparse.kr_biostochastic(m)
        return m, scl

    def _norm_seq(self, _map, tip_based, mean_type='geometric'):
        """
        Normalise a simple sequence map in place by the geometric mean of interacting contig pairs lengths.
        The map is assumed to be in starting order.

        :param _map: the target map to apply normalisation
        :param tip_based: treat the supplied map as a tip-based tensor
        :param mean_type: choice of mean (harmonic, geometric, arithmetic)
        :return: normalized map
        """
        if tip_based:
            _map = _map.astype(np.float)
            _tip_lengths = np.minimum(self.tip_size, self.order.lengths()).astype(np.float)
            fast_norm_seq(_map.coords, _map.data, _tip_lengths, self.tip_size)
            return _map

        else:
            _mean_func = mean_selector(mean_type)
            _len = self.order.lengths().astype(np.float)
            _map = _map.tolil().astype(np.float)
            for i in xrange(_map.shape[0]):
                _map[i, :] /= np.fromiter((1e-3 * _mean_func(_len[i],  _len[j])
                                           for j in xrange(_map.shape[0])), dtype=np.float)
            return _map.tocsr()

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

    def plot_clusters(self, fname, cl_map, simple=True, permute=False, **kwargs):
        """
        Plot the contact map, annotating the map with cluster names and boundaries
        :param fname: output file name
        :param cl_map: the cluster solution
        :param simple: True plot seq map, False plot the extent map
        :param permute: permute the map with the present order
        :param kwargs: additional options passed to plot()
        """
        if simple:
            # tick spacign simple the number of sequences in the cluster
            tick_locs = np.cumsum([0]+[len(cl_map[k]['seq_ids']) for k in cl_map])
            if self.is_tipbased():
                tick_locs *= 2
        else:
            # tick spacing depends on cumulative bins for sequences in cluster
            # cumulative bin count, excluding masked sequences
            csbins = [0]
            for k in cl_map:
                # get the order records for the sequences in cluster k
                _oi = self.order.order[cl_map[k]['seq_ids']]
                # count the cumulative bins at each cluster for those sequences which are not masked
                csbins.append(self.grouping.bins[cl_map[k]['seq_ids'][_oi['mask']]].sum() + csbins[-1])
            tick_locs = np.array(csbins, dtype=np.int)

        self.plot(fname, permute=permute, simple=simple, tick_locs=tick_locs, tick_labs=cluster_names(cl_map), **kwargs)

    def plot(self, fname, simple=False, tick_locs=None, tick_labs=None, norm=False, permute=False, pattern_only=False,
             dpi=180, width=25, height=22, zero_diag=False, alpha=0.01, robust=False):
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
        """

        plt.close()
        plt.style.use('ggplot')
        fig, ax = plt.subplots(1, 1)
        fig.set_figwidth(width)
        fig.set_figheight(height)

        if simple:
            # assume the map must be remade
            self.prepare_seq_map(self.min_len, self.min_sig, norm=norm, update_mask=True)
            m = self.get_processed_map(permute=permute)
        else:
            m = self.get_extent_map(norm=norm, permute=permute)

        if pattern_only:
            if zero_diag:
                m.setdiag(0)
            ax.spy(m.tocsr(), markersize=5 if simple else 1)

        else:
            # a dense array is necessary here
            m = m.toarray()
            if zero_diag:
                np.fill_diagonal(m, 0)
            m = np.log(m + alpha)

            seaborn.heatmap(m, robust=robust, square=True, linewidths=0, ax=ax, cbar=False)

        if tick_locs is not None:

            maj_labels = ticker.FixedFormatter('')
            if tick_labs is not None:
                min_labels = ticker.FixedFormatter(tick_labs)
            else:
                min_labels = maj_labels

            maj_ticks = ticker.FixedLocator(tick_locs[1:-1]-0.5)
            min_ticks = ticker.FixedLocator(tick_locs[:-1] + 0.5*np.diff(tick_locs))

            plt.yticks(fontsize=16)
            plt.xticks(fontsize=16, rotation=45)
            plt.tick_params(axis='both', which='minor', labelsize=20)

            ax.xaxis.set_major_formatter(maj_labels)
            ax.yaxis.set_major_formatter(maj_labels)
            ax.xaxis.set_major_locator(maj_ticks)
            ax.yaxis.set_major_locator(maj_ticks)

            ax.xaxis.set_minor_formatter(min_labels)
            ax.yaxis.set_minor_formatter(min_labels)
            ax.xaxis.set_minor_locator(min_ticks)
            ax.yaxis.set_minor_locator(min_ticks)

            # seaborn will not display the grid, so we make our own.
            ax.hlines(tick_locs, *ax.get_xlim(), color='grey', linewidth=1, linestyle='-.')
            ax.vlines(tick_locs, *ax.get_ylim(), color='grey', linewidth=1, linestyle='-.')

        plt.tight_layout()
        plt.savefig(fname, dpi=dpi)
        plt.close()

    def to_graph(self, norm=True, bisto=False, scale=False, extern_ids=False, convert_extent=False,
                 verbose=False):
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
        :param convert_extent: use the full extent map to calculate a full seq_map (if seq_map is tip based)
        :param verbose: debug output
        :return: graph of contigs
        """
        if extern_ids:
            _nn = lambda x: self.seq_info[x].name
        else:
            _nn = lambda x: x

        g = nx.Graph()
        g.add_nodes_from(_nn(u) for u in xrange(self.total_seq))

        if convert_extent:
            if verbose:
                print 'Converting extent map to seq map'
            m = self.extent_to_seq().astype(np.float)
        else:
            m = self.seq_map.astype(np.float)

        is4d = not convert_extent and self.is_tipbased()
        if norm:
            m = self._norm_seq(m, tip_based=is4d)
        if bisto:
            m, scl = self._bisto_seq(m, tip_based=is4d)

        if is4d:
            # convert 2x2 matrices to a single weight
            m = m.sum(axis=(2, 3)).to_scipy_sparse()

        if not sp.isspmatrix_coo(m):
            m = m.tocoo()

        scl = 1.0/m.max() if scale else 1

        for u, v, w in itertools.izip(m.row, m.col, m.data):
            g.add_edge(_nn(u), _nn(v), weight=w * scl)

        return g


def cluster_names(cl_soln, prefix='CL_'):
    """
    Pedantically work out how many digits are required for the largest cluster number,
    so we can add leading zeros to the remainder.
    :return: standardised cluster names
    """
    try:
        num_width = max(1, int(np.ceil(np.log10(max(cl_soln.keys())))))
    except OverflowError:
        num_width = 1

    return ['{0}{1:{2}d}'.format(prefix, k, num_width) for k in cl_soln]


def cluster_map(cm, method='louvain', convert_extent=False, work_dir=None, infomap_depth=2, verbose=False):
    """
    Cluster a contact map into groups, as an approximate proxy for "species" bins. This is recommended prior to
    applying TSP ordering, minimising the breaking of its assumptions that all nodes should connect and be traversed.
    :param cm: the contact map to cluster
    :param method: clustering algorithm to employ
    :param convert_extent: use the full extent map to calculate a full seq_map (if seq_map is tip based)
    :param work_dir: working directory to which files are written during clustering. Default: use system tmp
    :param infomap_depth: when using Infomap, treat clusters to this depth of the hierarchy.
    :param verbose: debug output
    :return:
    """

    def read_mcl(pathname):
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

    def read_table(pathname, seq_col=0, cl_col=1):
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

    def read_tree(pathname, max_depth=2):
        """
        Read a tree clustering file as output by Infomap.

        :param pathname: the path to the tree file
        :param max_depth: the maximum depth to consider when combining objects into clusters.
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
                cl_map.setdefault(tuple(['orig'] + hierarchy[:max_depth]), []).append(fields[-1])

            # rename clusters and order descending in size
            desc_key = sorted(cl_map, key=lambda x: len(cl_map[x]), reverse=True)
            for n, k in enumerate(desc_key):
                cl_map[n] = np.array(cl_map.pop(k), dtype=np.int)

        return cl_map

    if work_dir is None:
        work_dir = tempfile.gettempdir()
    else:
        assert os.path.exists(work_dir), 'supplied output path [{}] does not exist'.format(work_dir)

    seed = 1234
    base_name = 'cm_graph'
    g = cm.to_graph(norm=True, bisto=True, scale=True, convert_extent=convert_extent, verbose=verbose)

    method = method.lower()
    if method == 'louvain':
        cl_to_ids = louvain_cluster.cluster(g, no_iso=False, ragbag=False)
    elif method == 'mcl':
        ofile = os.path.join(work_dir, '{}.mcl'.format(base_name))
        nx.write_edgelist(g, path='{}.edges'.format(base_name), data=['weight'])
        subprocess.check_call(['mcl-14-137/bin/mcl', '{}.edges'.format(base_name), '--abc',
                               '-I', '1.2', '-o', ofile])
        cl_to_ids = read_mcl(ofile)
    elif method == 'simap':
        ofile = os.path.join(work_dir, '{}.simap'.format(base_name))
        nx.write_edgelist(g, path='{}.edges'.format(base_name), data=['weight'])
        subprocess.check_call(['java', '-jar', 'simap-1.0.0.jar', 'mdl', '-s', str(seed),
                               '-i', '1e-5', '1e-3', '-a', '1e-5',
                               '-g', '{}.edges'.format(base_name),
                               '-o', ofile])
        cl_to_ids = read_table(ofile)
    elif method == 'infomap':
        nx.write_edgelist(g, path='{}.edges'.format(base_name), data=['weight'])
        subprocess.check_call(['../infomap/Infomap', '-u', '-v', '-z', '-i', 'link-list',
                               '{}.edges'.format(base_name), work_dir])
        cl_to_ids = read_tree(os.path.join(work_dir, '{}.tree'.format(base_name)), max_depth=infomap_depth)
    elif method == 'slm':
        mod_func = '1'
        resolution = '1.0'
        opti_algo = '3'
        n_starts = '10'
        n_iters = '10'
        ofile = os.path.join(work_dir, '{}.slm'.format(base_name))
        nx.write_edgelist(g, path='{}.edges'.format(base_name), data=['weight'], delimiter='\t')
        subprocess.check_call(['java', '-jar', 'SLM4J/ModularityOptimizer.jar',
                               '{}.edges'.format(base_name), ofile,
                               mod_func, resolution, opti_algo, n_starts, n_iters, str(seed), '1'])
        cl_to_ids = read_table(ofile, seq_col=None, cl_col=0)
    else:
        raise RuntimeError('unimplemented method [{}]'.format(method))

    # standardise the result
    cl_soln = {cl_id: {'seq_ids': np.array(sorted([seq_id for seq_id in cl_to_ids[cl_id]]))} for cl_id in cl_to_ids}
    # add total extent to each cluster
    for cl_id in cl_soln:
        cl_soln[cl_id]['extent'] = cm.order.lengths()[cl_soln[cl_id]['seq_ids']].sum()

    return cl_soln


def order_clusters(cm, all_cl, min_len=None, min_sig=None, max_fold=None, min_extent=None, min_size=1):
    """
    Determine the order of sequences for a given cluster.

    :param cm: the contact map
    :param all_cl: the sequence clusters derived from the supplied contact map
    :param min_len: within a cluster exclude sequences that are too short (bp)
    :param min_sig: within a cluster exclude sequences with weak signal (counts)
    :param max_fold: within a cluster, exclude sequences that appear to be overly represented
    :param min_size: skip clusters which containt too few sequences
    :param min_extent: skip clusters whose total extent (bp) is too short
    :return: map of cluster orders, by cluster id
    """

    orders = {}
    for cl_id, cl_map in all_cl.iteritems():
        print 'Ordering cluster {}'.format(cl_id)

        # calculate the clusters extent
        print '\t{}bp extent'.format(cl_map['extent']),
        if cl_map['extent'] < min_extent:
            print '-- excluded cluster due to small extent'
            continue
        elif len(cl_map['seq_ids']) < min_size:
            print '-- excluded cluster contains too few sequences'
            continue
        else:
            print

        # enable only sequences in the cluster
        _mask = np.zeros_like(cm.order.mask_vector())
        _mask[cl_map['seq_ids']] = True

        try:
            _ord = cm.find_order(min_len, min_sig, max_fold=max_fold, verbose=True, norm=True, bisto=True,
                                 inverse_method='neglog', external_mask=_mask)
        except NoneAcceptedException as e:
            print '{} : cluster {} will be masked'.format(e.message, cl_id)
            continue
        except TooFewException as e:
            print '{} : ordering not possible for cluster {}'.format(e.message, cl_id)
            _ord = cm.order.accepted_order()

        orders[cl_id] = _ord

    return orders


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


if __name__ == '__main__':
    import argparse


    def out_name(base, suffix):
        return '{}_{}'.format(base, suffix)

    parser = argparse.ArgumentParser(description='Create a 3C fragment map from a BAM file')

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
    parser.add_argument('--max-fold', type=float, default=None, help='Maximum acceptable relative coverage [None]')
    parser.add_argument('--pickle', help='Picked contact map')
    parser.add_argument('fasta', help='Reference fasta sequence')
    parser.add_argument('bam', help='Input bam file in query order')
    parser.add_argument('cov', help='Input depth of coverage file')
    parser.add_argument('outbase', help='Output base file name')

    args = parser.parse_args()

    if args.pickle:
        cm = ContactMap.load(args.pickle)

    else:
        cm = ContactMap(args.bam,
                        args.cov,
                        args.fasta,
                        args.min_insert,
                        args.min_mapq,
                        min_len=args.min_reflen,
                        min_sig=args.min_signal,
                        max_fold=args.max_fold,
                        strong=args.strong,
                        bin_size=args.bin_size,
                        tip_size=args.tip_size,
                        precount=args.eta,
                        med_alpha=args.med_alpha)

        if cm.is_empty():
            import sys
            logger.info('Stopping as the map is empty')
            sys.exit(1)

        logger.info('Saving contact map instance...')
        cm.save(out_name(args.outbase, 'cm.p'))

    # cluster the entire map
    # cl_map = cluster_map(cm, method=args.cluster_method)

    # order each cluster
    # order_clusters(cm, cl_map)
