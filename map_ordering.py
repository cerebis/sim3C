#!/usr/bin/env python

import pysam
import community as com
import numpy as np
import networkx as nx
from Bio import SeqIO
from Bio.Restriction import Restriction
from collections import OrderedDict
from intervaltree import IntervalTree, Interval

from scoop import futures, shared
import scoop

from numba import jit, int64, float64, boolean


def get_enzyme_instance(enz_name):
    """
    Fetch an instance of a given restriction enzyme by its name.
    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    return getattr(Restriction, enz_name)


def upstream_dist(read, cut_sites, ref_length):
    """
    Upstream distance to nearest cut-site from the end of this
    reads alignment. Note: always assumes circular

    :param read: read in question
    :param cut_sites: the cut-sites for the aligning reference
    :param ref_length: the reference length
    :return:
    """
    if read.is_reverse:
        # reverse reads look left
        r_end = read.pos
        xi = np.searchsorted(cut_sites, r_end, side='left')
        if xi > 0:
            d = r_end - cut_sites[xi - 1]
        else:
            # crossing the origin
            d = ref_length - cut_sites[-1] + r_end
    else:
        # forward reads look right
        r_end = read.pos + read.alen
        xi = np.searchsorted(cut_sites, r_end, side='right')
        if xi < len(cut_sites):
            d = cut_sites[xi] - r_end
        else:
            # crossing the origin
            d = cut_sites[0] - (r_end - ref_length)
    return d


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


def inverse_edge_weights(g):
    """
    Invert the weights on a graph's edges
    :param g: the graph
    """
    for u, v in g.edges():
        g.edge[u][v]['weight'] = 1.0 / g[u][v]['weight']


def edgeiter_to_nodelist(edge_iter):
    """
    Create a list of nodes from an edge iterator
    :param edge_iter: edge iterator
    :return: list of node ids
    """
    nlist = []
    for ei in edge_iter:
        for ni in ei:
            if ni not in nlist:
                nlist.append(ni)
    return nlist


def dfs_weighted(g, source=None):
    """
    Depth first search
    :param g: 
    :param source: 
    :return: 
    """
    if source is None:
        # produce edges for all components
        nodes = g
    else:
        # produce edges for components with source
        nodes = [source]
    visited = set()
    for start in nodes:
        if start in visited:
            continue
        visited.add(start)
        stack = [(start, iter(sorted(g[start], key=lambda x: -g[start][x]['weight'])))]
        while stack:
            parent, children = stack[-1]
            try:
                child = next(children)
                if child not in visited:
                    yield parent, child
                    visited.add(child)
                    stack.append((child, iter(sorted(g[child], key=lambda x: -g[child][x]['weight']))))
            except StopIteration:
                stack.pop()


def decompose_graph(g):
    """
    Using the Louvain algorithm for community detection, as
    implemented in the community module, determine the partitioning
    which maximises the modularity. For each individual partition
    create the sub-graph of g

    :param g: the graph to decompose
    :return: the set of sub-graphs which form the best partitioning of g
    """
    decomposed = []
    part = com.best_partition(g)
    part_labels = np.unique(part.values())

    # for each partition, create the sub-graph
    for pi in part_labels:
        # start with a complete copy of the graph
        gi = g.copy()
        # build the list of nodes not in this partition and remove them
        to_remove = [n for n in g.nodes_iter() if part[n] != pi]
        gi.remove_nodes_from(to_remove)
        decomposed.append(gi)

    return decomposed


class Grouping:
    def __init__(self, seq_info, bin_size, min_sites):
        total_sites = 0
        self.bins = []
        self.map = OrderedDict()
        self.index = OrderedDict()
        self.bin_size = bin_size
        self.min_sites = min_sites

        for n, ix in enumerate(seq_info.keys()):

            num_sites = len(seq_info[ix]['sites'])

            if num_sites == 0:
                raise RuntimeError('Sequence has no sites')

            # required bins is the quotient, but the least number
            # of bins is 1.
            required_bins = num_sites / bin_size if num_sites >= bin_size else 1
            if num_sites > bin_size and float(num_sites % bin_size) / bin_size >= 0.5:
                # add an extra bin when there is a sufficient number of left-overs
                # but do not split a single bin sequence.
                required_bins += 1

            if required_bins == 1 and num_sites < min_sites:
                print '\tSequence had only {0} sites, which is less ' \
                      'than the minimum threshold {1}'.format(required_bins, min_sites)
                # remove the sequence from seq_info
                del seq_info[ix]
                continue

            # starting with integer indices, apply a scale factor to adjust
            # them matching the desired binning. This can then be digitised
            # to determine bin membership
            sclidx = np.arange(num_sites, dtype=np.int32) * float(required_bins) / num_sites

            # determine bin membership, but subtract one from each to get 0-based
            grp = np.digitize(sclidx, np.arange(required_bins)) - 1

            self.index[ix] = n
            self.bins.append(required_bins)
            self.map[ix] = np.column_stack((seq_info[ix]['sites'], grp))

            total_sites += num_sites

            print '\t{0} sites in {1} bins. Bins: {2}'.format(num_sites, required_bins, np.bincount(grp))

        self.bins = np.array(self.bins)

        self.total = total_sites

        # calculate centers of each fragment groupings. This will be used
        # when measuring separation.
        self.centers = {}
        for ix, locs in self.map.iteritems():
            self.centers[ix] = np.array(
                [locs[:, 0][np.where(locs[:, 1] == gi)].mean() for gi in np.unique(locs[:, 1])])

    def total_bins(self):
        return np.sum(self.bins)

    def find_nearest(self, ctg_idx, x):
        """
        Find the nearest site from a given position on a contig.
        :param ctg_idx: contig index
        :param x: query position
        :return: tuple of site and group number
        """
        group_sites = self.map[ctg_idx]
        ix = np.searchsorted(group_sites[:, 0], x)
        if ix == 0:
            return group_sites[0, :]
        elif ix == group_sites.shape[0]:
            return group_sites[-1, :]
        else:
            x1 = group_sites[ix - 1, :]
            x2 = group_sites[ix, :]
            if x - x1[0] <= x2[0] - x:
                return x1
            else:
                return x2

@jit(int64[:](int64[:,:], int64))
def find_nearest_jit(group_sites, x):
    """
    Find the nearest site from a given position on a contig.
    :param group_sites:
    :param x: query position
    :return: tuple of site and group number
    """
    ix = np.searchsorted(group_sites[:, 0], x)
    if ix == 0:
        return group_sites[0, :]
    elif ix == group_sites.shape[0]:
        return group_sites[-1, :]
    else:
        x1 = group_sites[ix - 1, :]
        x2 = group_sites[ix, :]
        if x - x1[0] <= x2[0] - x:
            return x1
        else:
            return x2


class SeqOrder:

    def __init__(self, seq_info):
        """
        Initial order is determined by the order of supplied sequence information dictionary. Sequenes
        are given new surrogate ids of consecutive integers. Member functions expect surrogate ids
        not original names.

        :param seq_info: sequence information dictionary
        """
        self.names = [k for k in seq_info]
        self.lengths = np.array([v['length'] for v in seq_info.values()])
        self.order = np.arange(len(self.names))

    def set_order(self, _ord):
        assert len(_ord) == len(self.order), 'New order was a different length'

        if isinstance(_ord, np.ndarray):
            self.order = _ord
        else:
            self.order = np.array(_ord)

    def get_name(self, _id):
        """
        :param _id: surrogate id of sequence
        :return: original name
        """
        return self.names[_id]

    def get_length(self, _id):
        """
        :param _id: surrogate id of sequence
        :return: sequence length
        """
        return self.lengths[_id]

    def before(self, a, b):
        """
        Test if a comes before another sequence in order.
        :param a: surrogate id of sequence a
        :param b: surrogate id of sequence b
        :return: True if a comes before b
        """
        assert a != b, 'Surrogate ids must be different'

        _o = self.order
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

        _o = self.order
        ix1, ix2 = np.where((_o == a) | (_o == b))[0]
        if ix1 > ix2:
            ix1, ix2 = ix2, ix1
        return np.sum(self.lengths[_o[ix1 + 1:ix2]])

    def shuffle(self):
        """
        Randomize order
        """
        np.random.shuffle(self.order)


@jit(boolean(int64[:], int64, int64))
def before(_ord, a, b):

    ix1 = None
    for i in xrange(len(_ord)):
        if a == _ord[i]:
            ix1 = i
            break

    ix2 = None
    for i in xrange(len(_ord)):
        if b == _ord[i]:
            ix2 = i
            break

    return ix1 < ix2


@jit(int64(int64[:], int64[:], int64, int64))
def intervening(_ord, _lengths, a, b):

    ix1 = None
    for i in xrange(len(_ord)):
        if a == _ord[i]:
            ix1 = i
            break

    ix2 = None
    for i in xrange(len(_ord)):
        if b == _ord[i]:
            ix2 = i
            break

    if ix1 == ix2:
        return 0

    if ix1 > ix2:
        ix1, ix2 = ix2, ix1

    return np.sum(_lengths[_ord[ix1 + 1:ix2]])


class ContactMap:
    def __init__(self, bam, enz_name, seq_file, bin_size, ins_mean, ins_sd, min_mapq, ref_min_len=0,
                 subsample=None, random_seed=None, strong=None, max_site_dist=None, spacing_factor=1.0,
                 linear=True, min_sites=1):

        self.strong = strong
        self.map_weight = None
        self.bin_size = bin_size
        self.ins_mean = ins_mean
        self.ins_sd = ins_sd
        self.min_mapq = min_mapq
        self.subsample = subsample
        self.random_state = np.random.RandomState(random_seed)
        self.strong = strong
        self.max_site_dist = max_site_dist
        self.spacing_factor = spacing_factor
        self.linear = linear
        self.min_sites = min_sites
        self.total_seq = len(bam.references)
        self.total_length = 0

        # test that BAM file is the corret sort order
        if 'SO' not in bam.header['HD'] or bam.header['HD']['SO'] != 'queryname':
            raise IOError('BAM file must be sorted by read name')

        if enz_name:
            self.enzyme = get_enzyme_instance(enz_name)
            print 'Cut-site spacing naive expectation: {0}'.format(4 ** self.enzyme.size)
        else:
            # non-RE based ligation-pairs
            self.enzyme = None
            print 'Warning: no enzyme was specified, therefore ligation pairs ' \
                  'are considered unconstrained in position.'
            raise RuntimeError('not implemented')

        self.unrestricted = False if enz_name else True

        # we'll pull sequences from a store as we go through the bam metadata
        seq_db = SeqIO.index_db(':memory:', seq_file, 'fasta')

        # determine the set of active sequences
        # where the first filtration step is by length
        self.seq_info = OrderedDict()
        offset = 0
        print 'Reading sequences...'
        for n, li in enumerate(bam.lengths):

            # minimum length threshold
            if li < ref_min_len:
                continue

            seq_id = bam.references[n]
            seq = seq_db[seq_id].seq
            if self.enzyme:

                sites = np.array(self.enzyme.search(seq, linear=self.linear))

                # must have at least one site
                if len(sites) < min_sites:
                    continue

                # don't apply spacing analysis to single site sequences
                medspc = 0.0
                if len(sites) > 1:
                    sites.sort()  # pedantic check for order
                    medspc = np.median(np.diff(sites))

                self.seq_info[n] = {'sites': sites,
                                    'invtree': IntervalTree(Interval(si, si + 1) for si in sites),
                                    'max_dist': self.spacing_factor * medspc,
                                    'offset': offset,
                                    'name': seq_id,
                                    'length': li}

                print '\tFound {0} cut-sites for \"{1}\" ' \
                      'med_spc: {2:.1f} max_dist: {3:.1f}'.format(len(sites), seq_id, medspc,
                                                                  self.spacing_factor * medspc)
            else:
                self.seq_info[n] = {'offset': self.total_len}

            offset += li

        self.total_len = offset

        print 'Initially: {0} sequences and {1}bp'.format(len(self.seq_info), self.sum_active())

        print 'Determining binning...'
        self.grouping = Grouping(self.seq_info, self.bin_size, self.min_sites)
        print 'After grouping: {0} sequences and {1}bp'.format(len(self.seq_info), self.sum_active())

        self.raw_map = None

        print 'Counting reads in bam file...'
        self.total_reads = bam.count(until_eof=True)

        print 'BAM file contains:\n' \
              '\t{0} sequences\n' \
              '\t{1}bp total length\n' \
              '\t{2} alignments'.format(self.total_seq, self.total_len, self.total_reads)

        # a simple associative map for binning on whole sequences
        id_set = set(self.seq_info)
        self.seq_map = OrderedDict((n, OrderedDict(zip(id_set, [0] * len(id_set)))) for n in id_set)

        # initialise the ordercm.order.order
        self.order = SeqOrder(self.seq_info)

        self._bin_map(bam)

    def sum_active(self):
        sum_len = 0
        for si, info in self.seq_info.iteritems():
            sum_len += info['length']
        return sum_len

    def num_active(self):
        return len(self.seq_info)

    def _init_map(self):
        n_bins = self.grouping.total_bins()

        print 'Initialising contact map of {0}x{0} fragment bins, ' \
              'representing {1} bp over {2} sequences'.format(n_bins, self.sum_active(), self.num_active())

        _m = {}
        for bi, seq_i in enumerate(self.grouping.map):
            _m[seq_i] = {}
            for bj, seq_j in enumerate(self.grouping.map):
                _m[seq_i][seq_j] = np.zeros((self.grouping.bins[bi], self.grouping.bins[bj]), dtype=np.int32)

        return _m

    def _calc_map_weight(self):
        w = 0
        for i in self.raw_map:
            for j in self.raw_map[i]:
                w += np.sum(self.raw_map[i][j])
        self.map_weight = w

    def _bin_map(self, bam):
        import tqdm

        def _simple_match(r):
            return not r.is_unmapped and r.mapq >= _mapq

        def _strong_match(r):
            return strong_match(r, self.strong, True, _mapq)

        # set-up match call
        _matcher = _strong_match if self.strong else _simple_match

        def _is_toofar_absolute(_x, _abs_max, _medspc):
            return _abs_max and _x > _abs_max

        def _is_toofar_medspc(_x, _abs_max, _medspc):
            return _x > _medspc

        @jit
        def _is_toofar_dummy(_x, _abs_max, _medspc):
            return False

        # set-up the test for site proximity
        if self.unrestricted:
            # all positions can participate, therefore restrictions are ignored.
            # TODO this would make more sense to be handled at arg parsing, since user should
            # be informed, rather tha silently ignoring.
            is_toofar = _is_toofar_dummy
        else:
            if self.max_site_dist:
                is_toofar = _is_toofar_absolute
            else:
                # default to spacing factor
                is_toofar = _is_toofar_medspc

        # initialise a map matrix for fine-binning and seq-binning
        self.raw_map = self._init_map()

        with tqdm.tqdm(total=self.total_reads) as pbar:

            # maximum separation being 3 std from mean
            _wgs_max = self.ins_mean + 2.0 * self.ins_sd

            _hic_max = self.max_site_dist

            # used when subsampling
            uniform = self.random_state.uniform
            sub_thres = self.subsample

            _mapq = self.min_mapq
            _seq_info = self.seq_info
            _active_ids = set(self.seq_info)

            wgs_count = 0
            dropped_3c = 0
            kept_3c = 0

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

                if r1.reference_id not in _active_ids or r2.reference_id not in _active_ids:
                    continue

                if sub_thres and sub_thres < uniform():
                    continue

                if not _matcher(r1) or not _matcher(r2):
                    continue

                if r1.is_read2:
                    r1, r2 = r2, r1

                assume_wgs = False
                if r1.is_proper_pair:

                    fwd, rev = (r2, r1) if r1.is_reverse else (r1, r2)
                    ins_len = rev.pos + rev.alen - fwd.pos

                    if fwd.pos <= rev.pos and ins_len < _wgs_max:
                        # assume this read-pair is a WGS read-pair, not HiC
                        wgs_count += 1
                        continue

                if not assume_wgs:

                    if not self.unrestricted:

                        r1_dist = upstream_dist(r1, _seq_info[r1.reference_id]['sites'], r1.reference_length)
                        if is_toofar(r1_dist, _hic_max, _seq_info[r1.reference_id]['max_dist']):
                            dropped_3c += 1
                            continue

                        r2_dist = upstream_dist(r2, _seq_info[r2.reference_id]['sites'], r2.reference_length)
                        if is_toofar(r2_dist, _hic_max, _seq_info[r2.reference_id]['max_dist']):
                            dropped_3c += 1
                            continue

                    kept_3c += 1
                    self.seq_map[r1.reference_id][r2.reference_id] += 1

                r1pos = r1.pos if not r1.is_reverse else r1.pos + r1.alen
                r2pos = r2.pos if not r2.is_reverse else r2.pos + r2.alen

                # r1_bin = self.grouping.find_nearest(r1.reference_id, r1pos)
                # r2_bin = self.grouping.find_nearest(r2.reference_id, r2pos)
                r1_bin = find_nearest_jit(self.grouping.map[r1.reference_id], r1pos)
                r2_bin = find_nearest_jit(self.grouping.map[r2.reference_id], r2pos)

                s1, s2 = r1.reference_id, r2.reference_id
                bi, bj = r1_bin[1], r2_bin[1]
                if s1 == s2:
                    # same sequence, make sure we apply the same order for diagonal matrix
                    if bi > bj:
                        bi, bj = bj, bi
                elif s1 > s2:
                    # different sequences, always use upper triangle
                    s1, s2 = s2, s1
                    bi, bj = bj, bi

                self.raw_map[s1][s2][bi, bj] += 1

        self._calc_map_weight()

        print 'Assumed {0} were WGS pairs'.format(wgs_count)
        print 'Kept {0} and dropped {1} long-range (3C-ish) pairs'.format(kept_3c, dropped_3c)
        print 'Total raw map weight {0}'.format(self.map_weight)

    def to_dense(self, norm=False):
        """
        Create a dense matrix from the internal representation (dict of dict of submatrices)
        :param norm: normalise matrix by bin cardinality 
        :return: dense matrix, indexes of sequence boundaries
        """
        m = np.vstack(np.hstack(self.raw_map[si].values()) for si in self.raw_map)
        if norm:
            bc = np.hstack(np.bincount(bci[:, 1]) for bci in self.grouping.map.values()).astype(np.float)
            m = np.divide(m, bc)
        return m

    def set_order(self, ord_array):
        self.order.set_order(ord_array)

    def get_ordered_bins(self):
        # make sure the two arrays are np arrays for fancy indexing tricks
        _order = np.array(self.order.order)
        _bins = np.array(self.grouping.bins)
        return _bins[_order]

    def _determine_block_shifts(self):
        """
        For the present ordering, calculate the block (whole-group) shifts necessary
        to permute the contact matrix to match the order.
        :return: list of tuples (start, stop, shift)
        """
        _order = self.order.order
        _bins = self.grouping.bins
        _shuf_bins = _bins[_order]

        shifts = []
        # iterate through current order, determine necessary block shifts
        curr_bin = 0
        for shuff_i, org_i in enumerate(_order):
            intervening_bins = _bins[:org_i]
            shft = -(curr_bin - np.sum(intervening_bins))
            shifts.append((curr_bin, curr_bin + _bins[org_i], shft))
            curr_bin += _shuf_bins[shuff_i]
        return shifts

    def _make_permutation_matrix(self):
        """
        Create the permutation matrix required to reorder the contact matrix to
        match the present ordering.
        :return: permutation matrix
        """
        # as a permutation matrix, the identity matrix causes no change
        perm_mat = np.identity(self.grouping.total_bins())
        block_shifts = self._determine_block_shifts()
        for si in block_shifts:
            if si[2] == 0:
                # ignore blocks which are zero shift.
                # as we started with identity, these are
                # already defined.
                continue
            # roll along columns, those rows matching each block with
            # a shift.
            pi = np.roll(perm_mat, si[2], 1)[si[0]:si[1], :]
            # put the result back into P, overwriting previous identity diagonal
            perm_mat[si[0]:si[1], :] = pi
        return perm_mat

    def get_ordered_map(self, norm=False):
        """
        Reorder the contact matrix to reflect the present order.
        :return: reordered contact matrix
        """
        # this may not be necessary, but we construct the full symmetrical map
        # before reordering, just in case something is reflected over the diagonal
        # and ends up copying empty space
        dm = self.to_dense(norm)
        full_map = np.tril(dm.transpose(), -1) + dm
        pm = self._make_permutation_matrix()
        # two applications of P is required to fully reorder the matrix
        # -- both on rows and columns
        # finally, convert the result back into a triangular matrix.
        return np.triu(np.dot(np.dot(pm, full_map), pm.T))

    def create_contig_graph(self):
        """
        Create a graph where contigs are nodes and edges linking
        nodes are weighted by the cumulative weight of contacts shared
        between them, normalized by the product of the number of fragments
        involved.
    
        :return: graph of contigs
        """
        _order = self.order.order
        _n = self.order.names

        g = nx.Graph()
        g.add_nodes_from(_order)

        for i in xrange(len(_order)):
            for j in xrange(i + 1, len(_order)):
                # networkx can choke on numpy floats when writing graphml format
                obs = float(np.sum(self.raw_map[_n[i]][_n[j]]))
                nsq = len(self.grouping.map[_n[i]]) * len(self.grouping.map[_n[j]])
                if obs > 0:
                    g.add_edge(i, j, weight=obs / nsq)
        return g

    def get_hc_order(self):
        from scipy.cluster.hierarchy import complete
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import dendrogram
        g = self.create_contig_graph()
        inverse_edge_weights(g)
        D = squareform(nx.adjacency_matrix(g).todense())
        Z = complete(D)
        return dendrogram(Z)['leaves']

    def get_adhoc_order(self):
        """
        Attempt to determine an initial starting order of contigs based
        only upon the cross terms (linking contacts) between each using
        graphical techniques.
    
        Beginning with a graph of contigs, where edges are weighted by
        contact weight, it is decomposed using Louvain modularity. Taking
        inverse edge weights, the shortest path of the minimum spanning
        tree of each subgraph is used to define an order. The subgraph
        orderings are then concatenated together to define a full
        ordering of the sample.
    
        Those with no edges, are included by appear in an indeterminate
        order.
    
        :return: order of contigs
        """
        g = self.create_contig_graph()
        decomposed_subgraphs = decompose_graph(g)

        isolates = []
        new_order = []
        for gi in decomposed_subgraphs:
            if gi.order() > 1:
                inverse_edge_weights(gi)
                mst = nx.minimum_spanning_tree(gi)
                inverse_edge_weights(gi)
                new_order.extend(edgeiter_to_nodelist(dfs_weighted(mst)))
            else:
                isolates.extend(gi.nodes())

        return new_order + isolates

    def plot(self, fname, norm=False, with_indexes=False):
        import matplotlib.pyplot as plt
        m = self.get_ordered_map(norm)
        m = m.astype(np.float)
        m += 0.01
        fig = plt.figure()
        img_h = 8
        fig.set_size_inches(img_h, img_h)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        if with_indexes:
            indexes = self.get_ordered_bins()
            ax.set_xticks(indexes.cumsum()-0.5)
            ax.set_yticks(indexes.cumsum()-0.5)
            ax.set_yticklabels([info['name'] for info in self.seq_info.values()])
            ax.set_xticklabels([info['name'] for info in self.seq_info.values()],
                               rotation=45, horizontalalignment='right')
            ax.grid(color='black', linestyle='-', linewidth=1.5)
        fig.add_axes(ax)
        ax.imshow(np.log(m), interpolation='none')
        plt.savefig(fname, dpi=180)

    def save(self, fname):
        import cPickle as p
        with open(fname, 'wb') as out_h:
            p.dump(self, out_h)

    @staticmethod
    def load(fname):
        import cPickle as p
        with open(fname, 'rb') as in_h:
            return p.load(in_h)


import random
from deap import creator, base, tools
from deap.algorithms import varOr
import itertools

from numpy import log, pi
from numba import vectorize, int32


@jit(float64(int32[:, :], float64[:, :]))
def poisson_lpmf(ob, ex):
    s = 0.0
    for i in xrange(ob.shape[0]):
        for j in xrange(ob.shape[1]):
            aij = ob[i, j]
            bij = ex[i, j]
            if aij == 0:
                s += -bij
            else:
                s += aij * log(aij/bij) + bij - aij + 0.5 * log(2.0 * pi * aij)
    return -s


@vectorize([float64(float64)])
def piecewise_3c(s):
    pr = 1./3.e6/10.
    if s < 1000e3:
        pr = 0.5 * (3e-6 * (1 - 3e-6)**s + 1/3e6)
    return pr


def calc_likelihood(an_order):
    """
    Calculate the logLikelihood of a given sequence configuration. The model is adapted from
    GRAAL. Counts are Poisson, with lambda parameter dependent on expected observed contact
    rate as a function of inter-fragment separation. This has been shown experimentally to be
    modelled effectively by a power-law (used here).

    :param an_order: an order with which to calculate the likelihood
    :return:
    """

    cm = shared.getConst('cm')

    #cm.order.order = an_order

    Nd = cm.map_weight

    sumL = 0.0

    _id = cm.order.names
    _centers = cm.grouping.centers
    _lengths = cm.order.lengths
    _map = cm.raw_map

    for i, j in itertools.combinations(cm.grouping.index.values(), 2):

        # inter-contig separation defined by cumulative
        # intervening contig length.
        # L = cm.order.intervening(i, j)
        L = intervening(an_order, _lengths, i, j)

        # bin centers
        centers_i = _centers[_id[i]]
        centers_j = _centers[_id[j]]

        # determine relative origin for measuring separation
        # between sequences. If i comes before j, then distances
        # to j will be measured from the end of i -- and visa versa
        if before(an_order, i, j):
            s_i = _lengths[i] - centers_i
            s_j = centers_j
        else:
            s_i = centers_i
            s_j = _lengths[j] - centers_j

        # separations between contigs in this ordering
        d_ij = np.abs(L + s_i[:, np.newaxis] - s_j)

        q_ij = piecewise_3c(d_ij)

        n_ij = _map[_id[i]][_id[j]]
        sumL += poisson_lpmf(n_ij, Nd * q_ij)

    return float(sumL),


def myMutate(individual, psnv, plarge, inv_size, inv_pb=0.5):
    # chx = np.random.randint(0, 3)
    n = len(individual)
    x = individual

    # inversions and transposition
    if np.random.uniform() < plarge:
        if np.random.uniform() < 0.5:
            # invert a range.
            d = np.random.binomial(inv_size, inv_pb)
            i1 = np.random.randint(0, n-d)
            x[i1:i1+d] = x[i1:i1+d][::-1]

        else:
            # transpose
            d = np.random.binomial(inv_size, inv_pb)
            i1 = np.random.randint(0, n-d)
            while True:
                i2 = np.random.randint(0, n-d)
                if i1 != i2:
                    break
            if i2 < i1:
                i2, i1 = i1, i2

            y = np.hstack((x[:i1], x[i1+d:]))
            y = np.hstack((y[:i2], x[i1:i1+d], y[i2:]))
            x[:] = y[:]
    else:
        # point mutations instead
        xset = set(range(n))
        for i1 in xrange(n):
            if np.random.uniform() < psnv:
                if np.random.uniform() < 0.5:
                    # move one
                    i2 = np.random.choice(list(xset - {i1}))
                    x[i1:i2+1] = np.roll(x[i1:i2+1], -1)
                else:
                    # swap two
                    i2 = np.random.choice(list(xset - {i1}))
                    x[i1], x[i2] = x[i2], x[i1]
    return x,


def myEaMuPlusLambda(population, toolbox, genlog, mu, lambda_, cxpb, mutpb, ngen,
                     stats=None, halloffame=None, verbose=__debug__):

    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    if halloffame is not None:
        halloffame.update(population)

    genlog.append(halloffame[0])

    record = stats.compile(population) if stats is not None else {}
    logbook.record(gen=0, nevals=len(invalid_ind), **record)
    if verbose:
        print logbook.stream

    # Begin the generational process
    for gen in range(1, ngen + 1):

        # Vary the population
        offspring = varOr(population, toolbox, lambda_, cxpb, mutpb)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Update the hall of fame with the generated individuals
        if halloffame is not None:
            halloffame.update(offspring)

        # Select the next generation population
        population[:] = toolbox.select(population + offspring, mu)

        # Update the statistics with the new population
        record = stats.compile(population) if stats is not None else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if verbose:
            print logbook.stream

        genlog.append(tools.selBest(population, 1))

    return population, logbook


def numpy_similar(a, b):
    return np.all(a == b)

creator.create('FitnessML', base.Fitness, weights=(1.0,))
creator.create('Individual', np.ndarray, fitness=creator.FitnessML)


if __name__ == '__main__':
    import argparse
    import pickle

    from datetime import timedelta
    from monotonic import monotonic

    class Timer:
        def __init__(self):
            self.start = monotonic()

        def reset(self):
            self.start = monotonic()

        def elapsed(self):
            end = monotonic()
            return timedelta(seconds=end - self.start)


    parser = argparse.ArgumentParser(description='Create a 3C fragment map from a BAM file')
    parser.add_argument('--strong', type=int, default=None,
                        help='Using strong matching constraint (minimum matches in alignments).')
    parser.add_argument('--bin-size', type=int, required=True, help='Size of bins in numbers of enzymatic sites')
    parser.add_argument('--enzyme', required=True, help='Enzyme used (case-sensitive)')
    parser.add_argument('--insert-mean', type=int, required=True, help='Expected mean of sequencing insert')
    parser.add_argument('--insert-sd', type=int, required=True, help='Expected standard deviation of insert')
    parser.add_argument('--min-sites', type=int, required=True, help='Minimum acceptable number of sites in a bin')
    parser.add_argument('--min-mapq', type=int, default=0, help='Minimum acceptable mapping quality [0]')
    parser.add_argument('--min-reflen', type=int, default=0, help='Minimum acceptable reference length [0]')
    parser.add_argument('--spc-factor', type=float, default=3.0,
                        help='Maximum sigma factor separation from read to nearest site [3]')
    parser.add_argument('fasta', help='Reference fasta sequence')
    parser.add_argument('bam', help='Input bam file')

    args = parser.parse_args()

    with pysam.AlignmentFile(args.bam, 'rb') as bam_file:

        cm = ContactMap(bam_file,
                        args.enzyme,
                        args.fasta,
                        args.bin_size,
                        args.insert_mean, args.insert_sd,
                        args.min_mapq,
                        spacing_factor=args.spc_factor,
                        min_sites=args.min_sites,
                        ref_min_len=args.min_reflen,
                        strong=args.strong)

        t = Timer()

        print 'Saving contact map instance...'
        cm.save('cm.p')
        print t.elapsed()

        print 'Saving raw contact maps as csv...'
        # original map
        t.reset()
        np.savetxt('raw.csv', cm.to_dense(), fmt='%d', delimiter=',')

        print 'Ordering by HC...'
        # hc order
        t.reset()
        o = cm.get_hc_order()
        cm.set_order(o)
        m = cm.get_ordered_map()
        print t.elapsed()
        print 'Saving hc contact maps as csv...'
        t.reset()
        np.savetxt('hc.csv', m, fmt='%d', delimiter=',')
        print t.elapsed()

        # adhoc order
        print 'Ordering by adhoc...'
        t.reset()
        o = cm.get_adhoc_order()
        cm.set_order(o)
        m = cm.get_ordered_map()
        print t.elapsed()
        print 'Saving adhoc contact maps as csv...'
        t.reset()
        np.savetxt('adhoc.csv', m, fmt='%d', delimiter=',')
        print t.elapsed()

        print 'Plotting adhoc image...'
        t.reset()
        cm.plot('adhoc.png', norm=True)
        print t.elapsed()

        print 'Beginning ES ordering...'
        t.reset()

        shared.setConst(cm=cm)

        history = tools.History()
        toolbox = base.Toolbox()
        active = list(cm.order.order)

        toolbox.register('indices', random.sample, active, len(active))
        toolbox.register('individual', tools.initIterate, creator.Individual, toolbox.indices)
        toolbox.register('population', tools.initRepeat, list, toolbox.individual)

        toolbox.register('evaluate', calc_likelihood)
        toolbox.register("mate", tools.cxOrdered)

        toolbox.register("mutate", myMutate, psnv=0.05, plarge=0.1, inv_size=10)
        toolbox.register("select", tools.selSPEA2)

        toolbox.decorate("mate", history.decorator)
        toolbox.decorate("mutate", history.decorator)

        toolbox.register('map', futures.map)

        MU, LAMBDA, NGEN = 50, 50, 200

        population = toolbox.population(n=MU)
        hof = tools.ParetoFront(similar=numpy_similar)
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean, axis=0)
        stats.register("std", np.std, axis=0)
        stats.register("min", np.min, axis=0)
        stats.register("max", np.max, axis=0)

        genlog = []

        pop, logbook = myEaMuPlusLambda(population, toolbox, genlog, mu=MU, lambda_=LAMBDA,
                                        cxpb=0.5, mutpb=0.5, ngen=NGEN,
                                        stats=stats, halloffame=hof)

        with open('gen_best.csv', 'w') as out_h:
            for n, order_i in enumerate(genlog):
                out_h.write('{0}: {1}\n'.format(n, ' '.join([str(v) for v in order_i])))

        cm.set_order(np.array(hof[0]))
        m = cm.get_ordered_map()
        print t.elapsed()

        print 'Saving es contact maps as csv...'
        t.reset()
        np.savetxt('es.csv', m, fmt='%d', delimiter=',')
        print t.elapsed()

        print 'Plotting es image...'
        t.reset()
        cm.plot('es.png', norm=True)
        print t.elapsed()
