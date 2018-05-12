from collections import OrderedDict
from numba import jit, vectorize, int64, int32, float64
import simple_sparse
import numpy as np
from numpy import log, pi
import scipy.sparse as sp
import itertools
import logging
import os
import pysam
import Bio.SeqIO as SeqIO
import networkx as nx
import ordering
import matplotlib.pyplot as plt
import cPickle
import seaborn
import uuid
from collections import namedtuple
from functools import partial

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


@vectorize([float64(float64)])
def piecewise_3c(s):
    pr = MIN_FIELD
    if s < 300e3:
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


def calc_likelihood(an_order, cm):
    """
    For a given order and ContactMap instance, calculate the log likelihood.
    :param an_order: ordered 1d container of identifiers.
    :param cm: instance of ContactMap with matching identifiers
    :return: log likelihood
    """

    borders = cm.grouping.borders
    centers = cm.grouping.centers
    lengths = cm.order.lengths
    raw_map = cm.raw_map

    total_obs = cm.map_weight()

    log_l = 0.0
    for i, j in itertools.combinations(xrange(cm.total_seq), 2):

        # inter-contig separation defined by cumulative
        # intervening contig length.
        L = intervening(an_order, lengths, i, j)

        # contig lengths
        li = lengths[i]
        lj = lengths[j]

        # bin centers for each contig, relative to middle of each contig
        c_ik = centers[i]
        c_jl = centers[j]

        # orientation of sequences
        s_i = an_order[i, 1]
        s_j = an_order[j, 1]

        # all separations between bins, including the potential intervening distance L
        d_ij = L + 0.5*(li + lj) + s_i * c_jl - s_j * c_ik.T

        # conversion to expected counts
        q_ij = total_obs * piecewise_3c(d_ij)

        # matrix element range which represents cross-terms between contig i and j
        i1, i2 = borders[i]
        j1, j2 = borders[j]

        # observed counts
        # for now this is converted to dense array as we need zeros
        n_ij = raw_map[i1:i2, j1:j2].toarray()

        # log likelihood
        log_l += poisson_lpmf3(n_ij, q_ij)

    return float(log_l),


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
        self.positions = np.argsort(self.order['pos'])

    def remap_gapless(self, gapless_indices):
        """
        Recover potentially sparse indices from a dense (gapless) index range. This can happen for
        algorithms which were supplied only filtered data.
        :param gapless_indices: dense index range, without gaps from masking
        :return: remappped indices with gaps
        """
        return self.positions[gapless_indices]

    def accepted_indices(self):
        """
        :return: indices (surrogate ids) of accepted sequences.
        """
        return set(np.where(self.order['mask'])[0])

    def excluded_indices(self):
        """
        :return: indices (surrogate ids) of excluded sequences.
        """
        return set(np.where(~self.order['mask'])[0])

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

    def set_only_order(self, _ord, implicit_excl=False):
        """
        Set only the order, while ignoring orientation. An ordering is defined
        as a 1D vectorize whose values are the positions and whose indices are
        the sequence surrogate ids. NOTE: This definition can be the opposite of what
        is returned by some ordering methods, and np.argsort(_v) should inverse
        the relation.

        If the order includes only active sequences, setting implicit_excl=True
        the method will implicitly assume unmentioned ids are those currently
        masked. An exception is raised if a masked sequence is included in the order.

        :param _ord: 1d ordering
        :param implicit_excl: implicitly extend the order to include unmentioned excluded sequences.
        """

        if len(_ord) < len(self.order):
            # some sanity checks
            assert implicit_excl, 'Use implicit_excl=True for automatic handling ' \
                                  'of orders only mentioning accepted sequences'
            assert len(_ord) == self.count_accepted(), 'new order must mention all ' \
                                                       'currently accepted sequences'
            assert len(set(_ord) & self.excluded_indices()) == 0, 'new order and excluded must not ' \
                                                                  'overlap when using implicit assignment'
            assert len(set(_ord) ^ set(self.accepted_indices())) == 0, 'incomplete new order supplied,' \
                                                                       'missing accepted ids'
            # assign the new orders
            self.order['pos'][_ord] = np.arange(len(_ord), dtype=np.int32)
        else:
            # just a simple complete order update
            assert len(_ord) == len(self.order), 'new order was a different length'
            assert len(set(_ord) ^ set(self.accepted_indices())) == 0, 'incomplete new order supplied,' \
                                                                       'missing accepted ids'
            self.order['pos'][_ord] = np.arange(len(_ord), dtype=np.int32)

        self._update_positions()

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

    def lengths(self, masked_only=False):
        if masked_only:
            return self.order['length'][self.order['mask']]
        else:
            return self.order['length']

    def get_length(self, _id):
        """
        :param _id: the surrogate id of sequence
        :return: sequence length
        """
        return self.order[_id]['length']

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
        inter_ix = self.positions[pa+1:pb]
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


class ContactMap:

    def read_covfile(self, fname):
        """
        Read depth of coverge data from samtools depth output. This file should (must) cover the regions of
        interest (tips) but need not contain all positions across each contig. The file can be gzip, bzip2 compressed
        with detection through suffix.
        :param fname: the depth of coverage file from samtools
        :return: dictionary containing per-base coverage over each contig and tip medians.
        """
        import io_utils
        with io_utils.open_input(fname) as in_h:
            rawd = {}
            for line in in_h:
                line = line.strip()
                if not line:
                    break
                _id, _pos, _dp = line.split()
                rawd.setdefault(_id, []).append((int(_pos), int(_dp)))

        cov_info = {}
        for seq in self.seq_info:
            tend = self.tip_size if seq.length >= 2*self.tip_size else seq.length/2
            x = np.array(rawd[seq.name], dtype=np.int)
            left_end = x[:, 0] < tend
            right_end = x[:, 0] > seq.length - tend
            cl = x[np.where(left_end)]
            cr = x[np.where(right_end)]
            x = np.array(x[np.where(left_end | right_end)])
            cov_info[seq.refid] = {'cov': x, 'med': int(np.median(x[:, 1])),
                                   'lmed': int(np.median(cl[:, 1])), 'rmed': int(np.median(cr[:, 1]))}

        return cov_info

    @staticmethod
    def merge_neighbours(x, thres, radius):
        """
        After filtering a 2d array on values larger than 'alpha', merge adjacent or
        nearly adjacent series into a range.
        :param x: the 2d series (0-pos, 1-value)
        :param thres: the threshold value on which to filtering (greater than)
        :param radius: the radius to which adjacent coords which will be merged
        :return: the list of ranges which are deemed contiguous
        """
        x2 = x[(x[:, 1] > thres), 0]
        d = np.ediff1d(x2, to_begin=1)
        ix = np.where(d > radius)[0]
        last = 0
        inv_list = []
        for i in ix:
            inv = x2[last:i]
            inv_list.append([inv.min(), inv.max()+1])
            last = i
        return inv_list

    def depth_intervals(self, alpha, radius):
        from intervaltree import IntervalTree
        no_go = {}
        for k, ci in self.cov_info.iteritems():
            inv_list = ContactMap.merge_neighbours(ci['cov'], alpha * ci['med'], radius)
            no_go[k] = IntervalTree()
            for inv in inv_list:
                no_go[k].addi(*inv)
            no_go[k].merge_overlaps()
            logger.info('RefId: {} excludes {} positions'.format(k, sum(i.end - i.begin for i in no_go[k])))
        return no_go

    def __init__(self, bam_file, cov_file, seq_file, min_insert, min_mapq=0, min_len=0, random_seed=None, strong=None,
                 bin_size=None, tip_size=None, precount=False, med_alpha=10):

        self.strong = strong
        self.bin_size = bin_size
        self.min_mapq = min_mapq
        self.min_insert = min_insert
        self.min_len = min_len
        self.random_state = np.random.RandomState(random_seed)
        self.strong = strong
        self.seq_info = []
        self.seq_map = None
        self.grouping = None
        self.raw_map = None
        self.order = None
        self.tip_size = tip_size
        self.precount = precount
        self.total_reads = None
        self.cov_info = None
        self.med_alpha = med_alpha
        self.processed_map = None
        self.bisto_scale = None

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
        return self.seq_map.sum()

    def is_empty(self):
        return self.map_weight() == 0

    def find_order(self, min_len, min_sig, norm=True, bisto=True, mean_type='geometric', inverse_method='inverse',
                   verbose=False):
        """
        Using LKH TSP solver, find the best ordering of the sequence map in terms of proximity ligation counts.
        Here, it is assumed that sequence proximity can be inferred from number of observed trans read-pairs, where
        an inverse relationship exists.
        :param min_len: the minimum sequence length to include
        :param min_sig: the minimum off-diagonal signal (trans read-pair count)
        :param norm: apply length normalisation before analysis
        :param bisto: make the matrix bistochastic before analysis
        :param mean_type: for length normalization, the chosen mean used between sequences
        :param inverse_method: the chosen inverse method for converting count (similarity) to distance
        :param verbose: debug output
        :return: the surrogate ids in optimal order
        """

        # prepare the input sequence map for analysis
        self.prepare_seq_map(min_len, min_sig, norm=norm, bisto=bisto, mean_type=mean_type, verbose=verbose)

        # we'll supply a partially initialized distance function
        dist_func = partial(ordering.similarity_to_distance,
                            method='inverse', alpha=1.2, beta=1, verbose=verbose)

        lkh_o = ordering.lkh_order(self.get_processed_map(), 'lkh_run', lkh_exe='../lkh/LKH-3.0/LKH', precision=1,
                                   seed=12345, runs=2, pop_size=50, dist_func=dist_func, special=False, verbose=True)

        # lkh ordering references the supplied matrix indices, not the surrogate ids.
        # we must map this consecutive set to the contact map indices.
        lkh_o = self.order.remap_gapless(lkh_o)
        if verbose:
            print lkh_o
        return lkh_o

    def _reorder_seq(self, _map):
        """
        Reorder a simple sequence map using the supplied map
        :param _map: the map to reorder
        :return: ordered map
        """
        assert sp.isspmatrix(_map), 'reordering expects a sparse matrix type'
        _order = self.order.accepted()
        assert _map.shape[0] == _order.shape[0], 'supplied map and unmasked order are different sizes'
        p = sp.lil_matrix(_map.shape)
        for i in xrange(len(_order)):
            p[i, _order[i]] = 1.
        p = p.tocsr()
        return p.dot(_map.tocsr()).dot(p.T)

    def _norm_seq(self, _map, mean_type='geometric'):
        """
        Normalise a simple sequence map in place by the geometric mean of interacting contig pairs lengths.
        The map is assumed to be in starting order.

        :param _map: the target map to apply normalisation
        :param mean_type: choice of mean (harmonic, geometric, arithmetic)
        :return: normalized map
        """

        _mean_func = mean_selector(mean_type)
        _len = self.order.lengths(masked_only=True).astype(np.float)
        _map = _map.tolil().astype(np.float)
        for i in xrange(_map.shape[0]):
            _map[i, :] /= np.fromiter((1e-3 * _mean_func(_len[i],  _len[j])
                                       for j in xrange(_map.shape[0])), dtype=np.float)
        return _map.tocsr()

    def set_filter_mask(self, min_len, min_sig, verbose=False):
        """
        Determine and set the filter mask using the specified constraints.
        :param min_len: minimum sequence length
        :param min_sig: minimum off-diagonal signal (counts)
        :param verbose: debug output
        """

        # boolean array masking sequences shorter than limit
        ix1 = self.order.lengths() >= min_len
        if verbose:
            print '\tmin_len removed', self.total_seq - ix1.sum()

        # boolean array masking sequences weaker than limit
        ix2 = simple_sparse.max_offdiag(self.seq_map) >= min_sig
        if verbose:
            print '\tmin_sig removed', self.total_seq - ix2.sum()

        # ix &= degen_mask
        # print 'without degens', ix.sum()

        #ix1 = np.fromiter((i not in degen_ids for i in xrange(len(m))), dtype=bool)
        #print ' degens', ix1.sum()

        # combine the masks and pass to instance of SeqOrder
        self.order.set_mask_only(ix1 & ix2)
        if verbose:
            print '\tAcceptance union', len(self.order.accepted())

    def prepare_seq_map(self, min_len, min_sig, norm=True, bisto=False, mean_type='geometric', verbose=False):
        """
        Prepare the sequence map (seq_map) by application of various filters and normalisations.

        :param min_len: minimum sequence length
        :param min_sig: minimum inter-sequence signal in raw counts.
        :param norm: normalisation by sequence lengths
        :param bisto: make the output matrix bistochastic
        :param mean_type: when performing normalisation, use "geometric, harmonic or arithmetic" mean.
        :param verbose: debug output
        """
        if verbose:
            print 'Prepare sequence map'
            print '\toriginally', self.seq_map.shape
        m = self.seq_map.astype(np.float)

        # apply mask filter to unwanted sequences
        self.set_filter_mask(min_len, min_sig, verbose)

        # if there are sequences to mask, remove them from the map
        if self.order.count_accepted() < self.total_seq:
            m = simple_sparse.compress(m.tocoo(), self.order.mask_vector())
            if verbose:
                print '\tfilter reduced', m.shape

        # apply length normalisation if requested
        if norm:
            m = self._norm_seq(m, mean_type)
            if verbose:
                print '\tnormed'

        # make map bistochastic if requested
        if bisto:
            m, scl = simple_sparse.kr_biostochastic(m)
            # retain the scale factors
            self.bisto_scale = scl
            if verbose:
                print '\tbalanced'

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
            _m = self._reorder_seq(_m)
            if verbose:
                print '\treordered'
        return _m

    def _bin_map(self, bam):
        import tqdm

        def _simple_match(r):
            return not r.is_unmapped and r.mapq >= _mapq

        def _strong_match(r):
            return strong_match(r, self.strong, True, _mapq)

        # set-up match call
        _matcher = _strong_match if self.strong else _simple_match

        def _on_tip_withlocs(p1, p2, l1, l2, _tip_size):
            # TODO this instantiation is faster/as-fast as taking it outside
            # TODO we would have to fix Sparse4D to just pass back indices i, j
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

        _on_tip = _always_true if self.tip_size is None else _on_tip_withlocs

        # initialise a sparse matrix for accumulating the map
        if not self.tip_size:
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
            _raw_map = simple_sparse.Sparse2DAccumulator(self.grouping.total_bins)
            _grouping_map = self.grouping.map
        else:
            _grouping_map = None
            _raw_map = None

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

                if _raw_map:
                    b1 = find_nearest_jit(_grouping_map[ix1], r1pos)
                    b2 = find_nearest_jit(_grouping_map[ix2], r2pos)

                    # maintain half-matrix
                    if b1 > b2:
                        b1, b2 = b2, b1

                    # tally all mapped reads for binned map, not just those considered in tips
                    _raw_map[b1, b2] += 1

                # for seq-map, we may reject reads outside of a defined tip region
                tip_info = _on_tip(r1pos, r2pos, l1, l2, _tip_size)
                if not tip_info[0]:
                    counts['not_tip'] += 1
                    continue

                counts['accepted'] += 1

                _seq_map[ix1, ix2] += tip_info[1]

        # default to always making matrices symmetric
        if self.bin_size:
            self.raw_map = _raw_map.get_coo()
            del _raw_map

        self.seq_map = _seq_map.get_coo()
        del _seq_map

        logger.info('Pair accounting: {}'.format(counts))
        logger.info('Total raw map weight {}'.format(self.map_weight()))

    def save(self, fname):
        with open(fname, 'wb') as out_h:
            cPickle.dump(self, out_h)

    @staticmethod
    def load(fname):
        with open(fname, 'rb') as in_h:
            return cPickle.load(in_h)

    def plot(self, fname=None, simple=False, norm=False, permute=False, pattern_only=False, with_names=False,
             dpi=180, width=25, height=22, zero_diag=False, alpha=0.01, robust=False):
        """
        Plot the contact map. This can either be as a sparse pattern (requiring much less memory but without visual
        cues about intensity), simple sequence or full binned map and normalized or permuted.
        :param fname: output file name
        :param simple: if true, sequence only map plotted
        :param norm: normalize intensities by geometric mean of lengths
        :param permute: reorder map to current order
        :param pattern_only: plot only a sparse pattern (much lower memory requirements)
        :param with_names: add sequence  names to axes. For large maps, this will be hard to read
        :param dpi: adjust DPI of output
        :param width: plot width in inches
        :param height: plot height in inches
        :param zero_diag: set bright self-interactions to zero
        :param alpha: log intensities are log (x + alpha)
        :param robust: use seaborn robust dynamic range feature
        """

        fig, ax = plt.subplots(1, 1)
        fig.set_figwidth(width)
        fig.set_figheight(height)

        if simple:
            m = self.get_seq_map(norm=norm, permute=permute)
        else:
            m = self.get_raw_map(norm=norm, permute=permute)

        if pattern_only:
            if zero_diag:
                m.setdiag(0)
            ax.spy(m.tocsr(), markersize=5 if simple else 1)

        else:
            # seaborn heatmaps take a little extra work

            _lab = [self.seq_info[i].name for i in self.order.order[:, 0]]

            m = m.toarray()
            if zero_diag:
                np.fill_diagonal(m, 0)
            m = np.log(m + alpha)

            if with_names:
                seaborn.heatmap(m, robust=robust, square=True, xticklabels=_lab, yticklabels=_lab,
                                linewidths=0, ax=ax, cbar=False)
            else:
                seaborn.heatmap(m, robust=robust, square=True, xticklabels=False, yticklabels=False,
                                linewidths=0, ax=ax, cbar=False)

        if with_names:
            if simple:
                ax.set_xticks(xrange(1, self.total_seq+1))
                ax.set_yticks(xrange(1, self.total_seq+1))
            else:
                _cbins = np.cumsum(self.grouping.bins[self.order.order[:, 0]])
                ax.set_xticks(_cbins - 0.5)
                ax.set_yticks(_cbins - 0.5)

                ax.grid(color='grey', linestyle='-.', linewidth=1)

        if fname:
            plt.savefig(fname, dpi=dpi)

    def _norm_raw(self, mean_type='geometric'):
        """
        Normalise a raw map in place by the geometric mean of interacting contig pairs lengths.
        :return:
        """
        _map = self.raw_map
        assert sp.isspmatrix(_map), 'Extent matrix is not a scipy matrix type'
        if not sp.isspmatrix_lil(_map):
            _map = _map.tolil()

        _mean_func = mean_selector(mean_type)
        _len = self.order.lengths
        _cbins = np.cumsum(self.grouping.bins)
        for row_i, col_dat in enumerate(_map.rows):
            i = np.searchsorted(_cbins, row_i, side='right')
            wi = np.fromiter((1e-3 * _mean_func(_len[i], _len[j])
                              for j in np.searchsorted(_cbins, col_dat, side='right')), dtype=np.float)
            _map.data[row_i] /= wi
        return _map

    def compress_extent(self):
        """
        Compress the extent map for each sequence that is presently masked. This will eliminate
        all bins which pertain to a given masked sequence.
        :return: a scipy.sparse.coo_matrix pertaining to only the unmasked sequences.
        """
        _m = self.raw_map
        assert sp.isspmatrix(_m), 'Extent matrix is not a scipy sparse matrix type'
        if not sp.isspmatrix_coo(_m):
            _m = _m.tocoo()

        _order = self.order.initial_order()
        _bins = self.grouping.bins
        _mask = _order[:, 2].astype(np.bool)

        # build a list of every accepted element.
        # TODO this could be done as below, without the memory footprint of realising all elements
        s = 0
        accept_bins = []
        accept_index = set(np.where(_mask)[0])
        for i in _order[:, 0]:
            if i in accept_index:
                accept_bins.extend([j+s for j in xrange(_bins[i])])
            s += _bins[i]
        # use a hashable container for quicker lookup
        accept_bins = set(accept_bins)

        # collect those values not in the excluded rows/columns
        keep_row = []
        keep_col = []
        keep_data = []
        for i in xrange(_m.nnz):
            if _m.row[i] in accept_bins and _m.col[i] in accept_bins:
                keep_row.append(_m.row[i])
                keep_col.append(_m.col[i])
                keep_data.append(_m.data[i])

        # after removal of those intervening, shift the coordinates of affected bins
        # TODO this could be moved into the loop above
        _shift = np.cumsum((~_mask) * _bins)
        _csbins = np.cumsum(_bins)
        for i in xrange(len(keep_row)):
            # rather than build a complete list of shifts across matrix, we'll
            # sacrifice some CPU and do lookups for the appropriate bin
            ir = np.searchsorted(_csbins, keep_row[i])
            ic = np.searchsorted(_csbins, keep_col[i])
            keep_row[i] -= _shift[ir]
            keep_col[i] -= _shift[ic]

        return sp.coo_matrix((keep_data, (keep_row, keep_col)), shape=_m.shape - _shift[-1])

    def _reorder_raw(self, _map):
        """
        Reorder the raw map using current order.
        :return: sparse CSR format permutation of the given map
        """
        _order = self.order.accepted()
        _bins = self.grouping.bins[self.order.mask_vector()]

        # create a permutation matrix
        p = sp.lil_matrix(_map.shape)
        _shuf_bins = _bins[_order]
        for i, oi in enumerate(_order):
            j_off = _bins[:oi].sum()
            i_off = _shuf_bins[:i].sum()
            for k in xrange(_bins[oi]):
                p[i_off+k, j_off+k] = 1

        # permute the raw_map
        p = p.tocsr()
        return p.dot(_map.tocsr()).dot(p.T)

    def get_raw_map(self, norm=False, permute=False):
        """
        Return the binned map, optionally normalised and permuted.
        :param norm: normalise intensities wrt to lengths
        :param permute: permute to the current order.
        :return: sparse matrix
        """
        if norm:
            m = self._norm_raw()
        else:
            m = self.raw_map.copy()

        if permute:
            m = self._reorder_raw(m)

        return m

    def get_extent_map(self, norm=True, bisto=False, permute=False, mean_type='geometric', verbose=False):

        if verbose:
            print 'Prepare extent map'
            print '\toriginally', self.raw_map.shape
        m = self.seq_map.astype(np.float32)

        # if there are sequences to mask, remove them from the map
        if len(self.order.accepted()) < self.total_seq:
            m = self.compress_extent()
            if verbose:
                print '\tfilter reduced', m.shape

        # apply length normalisation if requested
        # if norm:
        #     m = self._norm_seq(m, mean_type)
        #     if verbose:
        #         print '\tnormed'

        # make map bistochastic if requested
        # if bisto:
        #     m = simple_sparse.kr_biostochastic(m)
        #     if verbose:
        #         print '\tbalanced'

        # reorder using current order state
        if permute:
            m = self._reorder_raw(m)
            if verbose:
                print '\treordered'
        return m

    def create_contig_graph(self, norm=True, scale=False, extern_ids=False):
        """
        Create a graph where contigs are nodes and edges linking
        nodes are weighted by the cumulative weight of contacts shared
        between them, normalized by the product of the number of fragments
        involved.

        :param norm: normalize weights by site count, Ni x Nj
        :param scale: scale weights (max_w = 1)
        :param extern_ids: use the original external sequence identifiers for node ids
        :return: graph of contigs
        """
        _order = self.order.order[:, 0]
        _len = self.order.lengths

        if extern_ids:
            _nn = lambda x: self.seq_info[x].name
        else:
            _nn = lambda x: x

        g = nx.Graph()
        for u in _order:
            # as networkx chokes serialising numpy types, explicitly type cast
            g.add_node(_nn(u), length=int(_len[u]))

        max_w = 0
        for i in xrange(len(_order)):

            # number of sites on contig i
            ni = self.grouping.bins[i]
            i1, i2 = self.grouping.borders[i]

            for j in xrange(i + 1, len(_order)):

                j1, j2 = self.grouping.borders[j]

                # as networkx chokes serialising numpy types, explicitly type cast
                obs = float(self.raw_map[i1:i2, j1:j2].sum())
                if obs == 0:
                    continue

                # normalized by the square of the number of sites between i and j
                if norm:
                    # number of sites on contig j
                    nj = self.grouping.bins[j]
                    nrm_w = obs
                    if obs != 0:
                        nrm_w /= ni * nj
                    g.add_edge(_nn(i), _nn(j), weight=nrm_w, rawweight=obs)

                    if nrm_w > max_w:
                        max_w = nrm_w

                else:
                    g.add_edge(_nn(i), _nn(j), weight=obs, rawweight=obs)

                    if obs > max_w:
                        max_w = obs

        if scale:
            for u, v in g.edges_iter():
                if g[u][v] > 0:
                    g[u][v]['weight'] /= max_w

        return g


if __name__ == '__main__':
    import argparse
    import mapio

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
                        strong=args.strong,
                        bin_size=args.bin_size,
                        tip_size=args.tip_size,
                        precount=args.eta,
                        med_alpha=args.med_alpha)

        # if cm.is_empty():
        #     import sys
        #     logger.info('Stopping as the map is empty')
        #     sys.exit(1)

        logger.info('Saving contact map instance...')
        cm.save(out_name(args.outbase, 'cm.p'))

    # logger.info('Saving raw contact map...')
    # mapio.write_map(cm.to_dense(), out_name(args.outbase, 'raw'), args.format)
    #
    # logger.info('Plotting sequence image...')
    # cm.plot(out_name(args.outbase, 'simple.png'), simple=True)
    #
    # logger.info('Plotting starting binned image...')
    # cm.plot(out_name(args.outbase, 'start.png'))
    #
    # logL = calc_likelihood(cm.order.order, cm)
    # logger.info('Initial logL {}'.format(logL[0]))
    #
    # logger.info('Beginning HC ordering...')
    # hc_o = ordering.hc_order(cm.create_contig_graph())
    # cm.order.set_only_order(hc_o)
    # logger.info('HC logL {}'.format(calc_likelihood(cm.order.order, cm)[0]))
    # logger.info('Permuting map...')
    # m = cm.get_ordered_map()
    # logger.info('Plotting HC image...')
    # cm.plot(out_name(args.outbase, 'hc.png'))
    #
    # logger.info('Beginning adhoc ordering...')
    # ah_o = ordering.adhoc_order(cm.create_contig_graph(scale=True))
    # cm.order.set_only_order(ah_o)
    # logL = calc_likelihood(cm.order.order, cm)
    # logger.info('Adhoc logL {}'.format(logL[0]))
    # logger.info('Permuting map...')
    # m = cm.get_ordered_map()
    # logger.info('Plotting adhoc image...')
    # cm.plot(out_name(args.outbase, 'adhoc.png'))

    # logger.info('Beginning LHK ordering...')
    # lkh_o = ordering.lkh_order(cm.get_seq_map(norm=True), args.outbase, precision=1, runs=100)
    # cm.order.set_only_order(lkh_o)
    # logL = calc_likelihood(cm.order.order, cm)
    # logger.info('LKH logL {}'.format(logL[0]))
    # logger.info('Plotting LKH image...')
    # cm.plot(out_name(args.outbase, 'full_lkh.png'), permute=True, norm=True, with_names=False)
    # cm.plot(out_name(args.outbase, 'seq_lkh.png'), simple=True, permute=True, norm=True, with_names=False)
