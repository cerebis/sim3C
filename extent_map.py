from collections import OrderedDict
import numpy as np
from numpy import log, pi
import scipy.sparse as sparse
from numba import jit, vectorize, int64, int32, float64
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

        # all separations between bins, including the potential intervening distance L
        d_ij = L + 0.5*(li + lj) + c_jl - c_ik.T

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


class ExtentGrouping:

    def __init__(self, seq_info, bin_size):
        self.bins = []
        self.bin_size = bin_size
        self.map = []
        self.borders = []
        self.centers = []
        self.total_bins = 0

        for n, seq in enumerate(seq_info):

            if seq['length'] == 0:
                raise RuntimeError('Zeror length sequence ')

            # integer bin estimation
            num_bins = seq['length'] / bin_size
            if num_bins == 0:
                num_bins += 1
            # handle non-integer discrepancy by contracting/expanding all bins equally
            # the threshold between contract/expand being half a bin size
            if seq['length'] % bin_size != 0 and seq['length']/float(bin_size) - num_bins >= 0.5:
                num_bins += 1

            edges = np.linspace(0, seq['length'], num_bins+1, endpoint=True, dtype=np.int)

            self.bins.append(num_bins)

            # Per reference coordinate pairs (bin_edge, map_index)
            first_bin = self.total_bins
            last_bin = first_bin + num_bins
            self.map.append(np.vstack((edges[1:], np.arange(first_bin, last_bin))).T)
            self.borders.append(np.array([first_bin, last_bin], dtype=np.int))

            self.total_bins += num_bins

            c_nk = edges[:-1] + 0.5*(edges[1] - edges[0]) - 0.5*seq['length']
            self.centers.append(c_nk.reshape((1, len(c_nk))))
            logger.info('{}: {} bins'.format(n, num_bins))

        self.bins = np.array(self.bins)


class SeqOrder:

    def __init__(self, seq_info):
        """
        Initial order is determined by the order of supplied sequence information dictionary. Sequenes
        are given new surrogate ids of consecutive integers. Member functions expect surrogate ids
        not original names.

        :param seq_info: sequence information dictionary
        """
        self.index_to_refid = [si['refid'] for si in seq_info]
        self.refid_to_index = {si['refid']: n for n, si in enumerate(seq_info)}
        self.lengths = np.array([si['length'] for si in seq_info])
        _ord = np.arange(len(seq_info))
        _ori = np.zeros(len(seq_info), dtype=np.int)
        self.order = np.vstack((_ord, _ori)).T

    @staticmethod
    def make_ord_and_ori(_ord):
        """
        Create an order/orientation matrix from a simple order list. Here, it is assumed that
        all objects are in their initial orientation.

        :param _ord: a list (or other container) of object identifiers.
        :return: an order and orientation matrix.
        """
        return np.vstack((np.asarray(_ord), np.zeros(len(_ord), dtype=np.int))).T

    def set_only_order(self, _ord):
        """
        Set only the order, assuming orientation is ignored
        :param _ord: 1d ordering
        """
        assert len(_ord) == len(self.order), 'new order was a different length'
        self.order[:, 0] = np.asarray(_ord)

    def set_ord_and_ori(self, _ord):
        """
        Set an order where _ord is a 2d (N x 2) array of numeric id and orientation (0,1)
        :param _ord: 2d ordering and orientation (N x 2)
        """
        assert isinstance(_ord, np.ndarray), 'ord/ori was not a numpy array'
        assert _ord.shape == self.order.shape, 'new ord/ori has different dimensions'
        self.order = _ord

    def get_length(self, i):
        """
        :param i: surrogate id of sequence
        :return: sequence length
        """
        return self.lengths[i]

    def before(self, a, b):
        """
        Test if a comes before another sequence in order.
        :param a: surrogate id of sequence a
        :param b: surrogate id of sequence b
        :return: True if a comes before b
        """
        assert a != b, 'Surrogate ids must be different'

        _o = self.order[:, 0]
        ia = np.argwhere(_o == a)
        ib = np.argwhere(_o == b)
        return ia < ib

    def intervening(self, a, b):
        """
        For the current order, calculate the length of intervening
        sequences between sequence a and sequence b.

        :param a: surrogate id of sequence a
        :param b: surrogate id of sequence b
        :return: total length of sequences between a and b.
        """
        assert a != b, 'Surrogate ids must be different'

        _o = self.order[:, 0]
        ix1, ix2 = np.where((_o == a) | (_o == b))[0]
        if ix1 > ix2:
            ix1, ix2 = ix2, ix1
        return np.sum(self.lengths[_o[ix1 + 1:ix2]])

    def shuffle(self):
        """
        Randomize order
        """
        np.random.shuffle(self.order)


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

    def __init__(self, bam_file, seq_file, bin_size, min_sep, min_mapq, min_len=0, random_seed=None, strong=None):

        self.strong = strong
        self.bin_size = bin_size
        self.min_mapq = min_mapq
        self.min_sep = min_sep
        self.min_len = min_len
        self.random_state = np.random.RandomState(random_seed)
        self.strong = strong
        self.seq_info = []
        self.seq_map = None
        self.grouping = None
        self.raw_map = None
        self.order = None

        with pysam.AlignmentFile(bam_file, 'rb') as bam:

            # test that BAM file is the corret sort order
            if 'SO' not in bam.header['HD'] or bam.header['HD']['SO'] != 'queryname':
                raise IOError('BAM file must be sorted by read name')

            # ':memory:'
            idxfile = '{}.db'.format(seq_file)
            # if os.path.exists(idxfile):
            #     logging.warning('Removed preexisting target path "{}" for sequence database'.format(idxfile))
            #     os.unlink(idxfile)
            if os.path.exists(idxfile):
                logging.warn('Found existing fasta index')
            seq_db = SeqIO.index_db(idxfile, seq_file, 'fasta')

            try:
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
                        ref_count['seq_missing'] += 1
                        continue

                    seq = seq_db[seq_id].seq
                    assert len(seq) == li, 'Reported sequence length in bam and ' \
                                           'fasta do not agree. {}: {} vs {} bp'.format(seq_id, len(seq), li)

                    self.seq_info.append({'offset': offset, 'refid': n, 'name': seq_id, 'length': li})
                    logger.debug('Found {} bp sequence for {}'.format(li, seq_id))

                    offset += li
            finally:
                seq_db.close()

            # total extent covered
            self.total_len = offset
            self.total_seq = len(self.seq_info)

            logger.info('Accepted {} sequences covering {} bp'.format(self.total_seq, self.total_len))
            logger.info('References excluded: {}'.format(ref_count))

            logger.info('Determining binning...')
            self.grouping = ExtentGrouping(self.seq_info, self.bin_size)

            logger.info('Counting reads in bam file...')
            self.total_reads = bam.count(until_eof=True)

            logger.info('BAM file contains {0} alignments'.format(self.total_reads))

            # As sequences from the map file may have been excluded above,
            # we use a dict of dicts to represent the potentially sparse sequence map.
            # memory use could be reduced if symmetry was exploited.
            self.seq_map = sparse.lil_matrix((self.total_seq, self.total_seq), dtype=np.int32)

            # initialise the order
            self.order = SeqOrder(self.seq_info)

            # calculate the mapping
            self._bin_map(bam)

    def _init_map(self):
        n_bins = self.grouping.total_bins

        logger.info('Initialising contact map of {0}x{0} fragment bins, '
                    'representing {1} bp over {2} sequences'.format(n_bins, self.total_len, self.total_seq))

        return sparse.lil_matrix((n_bins, n_bins), dtype=np.int32)

    def map_weight(self):
        return self.raw_map.sum()

    def is_empty(self):
        return self.map_weight() == 0

    def to_dense(self):
        """
        Convert the raw binned map into a contiguous numpy array. For large matrices, this
        can be memory consumptive.
        :return: a contiguous numpy array
        """
        return self.raw_map.toarray()

    def _reorder_seq(self, _map):
        """
        Reorder a simple sequence map using the supplied map
        :param _map: the map to reorder
        :return: ordered map
        """
        assert sparse.isspmatrix(_map), 'reordering expects a sparse matrix type'
        _order = self.order.order[:, 0]
        n = _map.shape[0]
        p = sparse.lil_matrix((n, n))
        for i in xrange(len(_order)):
            p[i, _order[i]] = 1.
        p = p.tocsr()
        return p.dot(_map.tocsr()).dot(p.T)

    def _norm_seq(self):
        """
        Normalise a simple sequence map in place by the geometric mean of interacting contig pairs lengths.
        :param _map: the map to normalise
        """
        _len = self.order.lengths
        _map = self.seq_map.astype(np.float)
        for i, ri in enumerate(_map.rows):
            _map.data[i] /= np.fromiter(((1e-6 * _len[i] * _len[j])**0.5 for j in ri), dtype=np.float)
        return _map

    def get_seq_map(self, norm=False, permute=False):
        """
        Return the simple per-sequence map, optionally normalised and permuted.
        :param norm: normalise intensities wrt to lengths
        :param permute: permute to the current order.
        :return: sparse matrix
        """
        if norm:
            m = self._norm_seq()
        else:
            m = self.seq_map.copy()

        if permute:
            m = self._reorder_seq(m)

        return m

    def _bin_map(self, bam):
        import tqdm

        def _simple_match(r):
            return not r.is_unmapped and r.mapq >= _mapq

        def _strong_match(r):
            return strong_match(r, self.strong, True, _mapq)

        # set-up match call
        _matcher = _strong_match if self.strong else _simple_match

        # initialise a map matrix for fine-binning and seq-binning
        self.raw_map = self._init_map()

        with tqdm.tqdm(total=self.total_reads) as pbar:

            _min_sep = self.min_sep
            _mapq = self.min_mapq
            _idx = self.order.refid_to_index

            _grouping_map = self.grouping.map
            _seq_map = self.seq_map
            _raw_map = self.raw_map

            counts = OrderedDict({
                'accepted': 0,
                'too_close': 0,
                'ref_excluded': 0,
                'skipped': 0,
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

                if r1.is_read2:
                    r1, r2 = r2, r1

                # when possible, estimate separation of pair and exclude if unacceptable
                if r1.is_proper_pair:

                    fwd, rev = (r2, r1) if r1.is_reverse else (r1, r2)
                    ins_len = rev.pos + rev.alen - fwd.pos

                    if fwd.pos <= rev.pos and ins_len < _min_sep:
                        counts['too_close'] += 1
                        continue

                counts['accepted'] += 1

                r1pos = r1.pos if not r1.is_reverse else r1.pos + r1.alen
                r2pos = r2.pos if not r2.is_reverse else r2.pos + r2.alen

                ix1 = _idx[r1.reference_id]
                ix2 = _idx[r2.reference_id]

                # maintain just a half-matrix
                if ix2 < ix1:
                    ix1, ix2 = ix2, ix1
                    r1pos, r2pos = r2pos, r1pos

                # sequence-to-sequence map
                _seq_map[ix1, ix2] += 1

                b1 = find_nearest_jit(_grouping_map[ix1], r1pos)
                b2 = find_nearest_jit(_grouping_map[ix2], r2pos)

                # maintain half-matrix
                if b1 > b2:
                    b1, b2 = b2, b1

                _raw_map[b1, b2] += 1

        # default to always making matrices symmetric
        ContactMap._symm_map(_raw_map)
        ContactMap._symm_map(_seq_map)

        logger.info('Pair accounting: {}'.format(counts))
        logger.info('Total raw map weight {}'.format(self.map_weight()))

    def save(self, fname):
        with open(fname, 'wb') as out_h:
            cPickle.dump(self, out_h)

    @staticmethod
    def load(fname):
        with open(fname, 'rb') as in_h:
            return cPickle.load(in_h)

    @staticmethod
    def _symm_map(_map):
        """
        In place conversion of a half-diagonal representation to fully symmetric.
        Handles both sparse and numpy matrices, assuming the sparse subtype
        supports assignment. It is assumed to be upper half diagonal.
        :param _map:
        """
        if sparse.isspmatrix(_map):
            _map += sparse.tril(_map.T, k=-1)
        else:
            _map += np.tril(_map, k=-1)

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

            _lab = [self.seq_info[i]['name'] for i in self.order.order[:, 0]]

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

    def _norm_raw(self):
        """
        Normalise a raw map in place by the geometric mean of interacting contig pairs lengths.
        :param _map: map to normalise
        :return:
        """
        _map = self.raw_map.astype(np.float64)
        _len = self.order.lengths
        _cbins = np.cumsum(self.grouping.bins)
        for row_i, col_dat in enumerate(_map.rows):
            i = np.searchsorted(_cbins, row_i, side='right')
            wi = np.fromiter(((1e-6 * _len[i] * _len[j])**0.5
                              for j in np.searchsorted(_cbins, col_dat, side='right')), dtype=np.float)
            _map.data[row_i] /= wi
        return _map

    def _reorder_raw(self, _map):
        """
        Reorder the raw map using current order.
        :return: sparse permutation matrix (lil)
        """
        _order = self.order.order[:, 0]
        _bins = self.grouping.bins
        n = _map.shape[0]

        # create a permutation matrix
        p = sparse.lil_matrix((n, n))
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


    # def _determine_block_shifts(self):
    #     """
    #     For the present ordering, calculate the block (whole-group) shifts necessary
    #     to permute the contact matrix to match the order.
    #     :return: list of tuples (start, stop, shift)
    #     """
    #     _order = self.order.order[:, 0]
    #     _bins = self.grouping.bins
    #     _shuf_bins = _bins[_order]
    #
    #     shifts = []
    #     # iterate through current order, determine necessary block shifts
    #     curr_bin = 0
    #     for shuff_i, org_i in enumerate(_order):
    #         intervening_bins = _bins[:org_i]
    #         shft = -(curr_bin - np.sum(intervening_bins))
    #         shifts.append((curr_bin, curr_bin + _bins[org_i], shft))
    #         curr_bin += _shuf_bins[shuff_i]
    #     return shifts
    #
    # def _make_permutation_matrix_old(self):
    #     """
    #     Create the permutation matrix required to reorder the contact matrix to
    #     match the present ordering.
    #     :return: permutation matrix
    #     """
    #     # as a permutation matrix, the identity matrix causes no change
    #     perm_mat = np.identity(self.grouping.total_bins)
    #     block_shifts = self._determine_block_shifts()
    #     for si in block_shifts:
    #         if si[2] == 0:
    #             # ignore blocks which are zero shift.
    #             # as we started with identity, these are
    #             # already defined.
    #             continue
    #         # roll along columns, those rows matching each block with
    #         # a shift.
    #         pi = np.roll(perm_mat, si[2], 1)[si[0]:si[1], :]
    #         # put the result back into P, overwriting previous identity diagonal
    #         perm_mat[si[0]:si[1], :] = pi
    #     return perm_mat

    # def get_ordered_map(self):
    #     """
    #     Reorder the contact matrix to reflect the present order.
    #     :return: reordered contact matrix
    #     """
    #     m = self.symm_map()
    #     p = self._make_permutation_matrix()
    #     return p.dot(m).dot(p.T) #np.dot(np.dot(p, m), p.T)

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
            _nn = lambda x: self.seq_info[x]['name']
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
    parser.add_argument('--strong', type=int, default=None,
                        help='Using strong matching constraint (minimum matches in alignments).')
    parser.add_argument('--bin-size', type=int, required=True, help='Size of bins in bp')
    parser.add_argument('--min-sep', type=int, required=True, help='Minimum pair separation')
    parser.add_argument('--min-mapq', type=int, default=0, help='Minimum acceptable mapping quality [0]')
    parser.add_argument('--min-reflen', type=int, default=0, help='Minimum acceptable reference length [0]')
    parser.add_argument('--pickle', help='Picked contact map')
    parser.add_argument('fasta', help='Reference fasta sequence')
    parser.add_argument('bam', help='Input bam file')
    parser.add_argument('outbase', help='Output base file name')

    args = parser.parse_args()

    if args.pickle:
        cm = ContactMap.load(args.pickle)

    else:
        cm = ContactMap(args.bam,
                        args.fasta,
                        args.bin_size,
                        args.min_sep,
                        args.min_mapq,
                        min_len=args.min_reflen,
                        strong=args.strong)

        if cm.is_empty():
            import sys
            logger.info('Stopping as the map is empty')
            sys.exit(1)

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

    logger.info('Beginning LHK ordering...')
    lkh_o = ordering.lkh_order(cm.simple_map(norm=True), args.outbase, precision=10)
    cm.order.set_only_order(lkh_o)
    logL = calc_likelihood(cm.order.order, cm)
    logger.info('LKH logL {}'.format(logL[0]))
    logger.info('Permuting map...')
    m = cm.get_ordered_map()
    logger.info('Plotting LKH image...')
    cm.plot(out_name(args.outbase, 'lkh.png'))
