import numpy as np
import networkx as nx
import scipy.sparse as sparse
import community
import polo
import lap
import os
import re
import subprocess

from scipy.cluster.hierarchy import ward, complete
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist


def hc_order(g, metric='cityblock', method='ward', use_olo=True):
    """
    Basic hierarchical clustering to determine an order of contigs, using optimal leaf ordering (poor time complexity)
    to adjust tips.
    :param g: the graph to order
    :param metric: any
    :param method: ward or complete
    :param use_olo: use optimal leaf ordering
    :return: an ordering
    """

    d = pdist(nx.adjacency_matrix(g).todense(), metric=metric)
    if method == 'ward':
        z = ward(d)
    elif method == 'complete':
        z = complete(d)
    else:
        raise RuntimeError('unsupported method: {}'.format(method))

    if use_olo:
        z = polo.optimal_leaf_ordering(z, d)

    return np.array(dendrogram(z, no_plot=True)['leaves'])


def adhoc_order(g, alpha=1.0):
    """
    Attempt to determine an ordering based only upon cross-terms
    between contigs using graphical techniques.

    1. Begin with a contig graph where edges are weighted by contact frequency.
    2. The graph is then partitioned into subgraphs using Louvain modularity.
    3. Using inverse edge weights, the shortest path of the minimum spanning
    tree of each subgraph is used to define an order.
    4. The subgraph orderings are then concatenated together to define a full
    ordering of the sample.
    5. TODO add optimal leaf ordering is possible.
    6. Unconnected contigs are included by order of appearance.

    :param g: graph to order
    :param alpha: additive constant used in inverse weighting.
    :return: an ordering
    """
    sg_list = decompose_graph(g)

    # calculate inter-subgraph weights for use in ordering
    w = inter_weight_matrix(g, sg_list, norm=True)

    # perform LAP, but first convert weight matrix to a "cost" matrix
    ord_x, ord_y = lap.lapjv(1.0 / (w + alpha), return_cost=False)

    # reorder the subgraphs, just using row order
    sg_list = [sg_list[i] for i in ord_x]

    # now find the order through each subgraph
    isolates = []
    new_order = []
    for gi in sg_list:
        # if there is more than one node
        if gi.order() > 1:
            inverse_edge_weights(gi)
            mst = nx.minimum_spanning_tree(gi)
            inverse_edge_weights(gi)
            new_order.extend(edgeiter_to_nodelist(dfs_weighted(mst)))
        else:
            isolates.extend(gi.nodes())

    return np.array(new_order + isolates)


def decompose_graph(g, reso=1.0):
    """
    Using the Louvain algorithm for community detection, as
    implemented in the community module, determine the partitioning
    which maximises the modularity. For each individual partition
    create the sub-graph of g

    :param g: the graph to decompose
    :param reso: louvain clustering threshold, smaller -> more partitions
    :return: the set of sub-graphs which form the best partitioning of g
    """

    decomposed = []
    part = community.best_partition(g, resolution=reso)
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


def inter_weight_matrix(g, sg, norm=True):
    """
    Calculate the weight of interconnecting edges between subgraphs identified from
    Louvain decomposition.

    :param g: the original graph
    :param sg: the list of subgraphs
    :param norm: normalize the counts by the number of shared edges
    :return: two matrices, 'w'  the weights of shared edges, 'n' the counts of shared edges
    """

    nsub = len(sg)
    w = np.zeros((nsub, nsub))
    if norm:
        n = np.zeros_like(w, dtype=np.int)

    # for each subgraph i
    for i in xrange(nsub):

        # for every node in subgraph i
        for u in sg[i].nodes_iter():

            # for every other subgraph j
            for j in xrange(i+1, nsub):

                # for every node in subgraph j
                for v in sg[j].nodes_iter():

                    # sum weight of edges connecting subgraphs i and j
                    if g.has_edge(u, v):
                        w[i, j] += g[u][v]['rawweight']
                        if norm:
                            n[i, j] += 1

    if norm:
        # only touch non-zero elements
        ix = np.where(n > 0)
        w[ix] /= n[ix]

    return w


def dfs_weighted(g, source=None):
    """
    Depth first search, guided by edge weights
    :param g: the graph to traverse
    :param source: the starting node used during recursion
    :return: list of nodes
    """
    # either produce edges for all components or only those in list
    if source is None:
        nodes = g
    else:
        nodes = [source]

    visited = set()
    for start in nodes:
        if start in visited:
            continue
        visited.add(start)
        # for node 'start' visit neighbours by edge weight
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


def inverse_edge_weights(g, alpha=1.0):
    """
    Invert the weights on a graph's edges
    :param g: the graph
    :param alpha: additive constant in denominator to avoid DBZ
    """
    for u, v in g.edges():
        g.edge[u][v]['weight'] = 1.0 / (g[u][v]['weight'] + alpha)


def reciprocal_counts(m, alpha=0.1):
    """
    Express a measure of similarity as a distance by taking the reciprocal. The array
    is expected to be strictly positive reals. Additional fiddly steps are taken to
    better support LKH, and its restriction to integers when expressing explicit
    distances.
    1. the diagonal is zeroed
    2. similarity is converted to [0,1] prior to taken recip.
    3. effort is taken to ensure the shortest distance (largest similarity) is 1 and not 0.
    :param m: input similarity matrix, which is also written to.
    :param alpha: additive smoothing factor, which avoids zero.
    :return: distance matrix
    """
    assert m.min() >= 0, 'the input matrix should be stricly positive'
    assert np.issubdtype(m.dtype, np.float64) or np.issubdtype(m.dtype, np.float32), 'matrix must be floats'

    # begin with removal of diagonal self-self interactions
    np.fill_diagonal(m, 0)
    m += alpha
    # scale elements [0,1]
    m = m / m.max()
    # take inverse
    m = 1.0 / m
    # rescale all values
    # As this will be integer truncated, we make sure the smallest value is > 1
    # the reason for this is to protect the shortest distance elements
    # from going to zero when the matrix is converted to integers (requirement of LKH)
    m *= 1.01 / m.min()
    np.fill_diagonal(m, 0)
    return m


def scale_mat(M, _min, _max):
    """
    In-place rescaling of matrix elements to be within the range [_min, _max]
    :param M: the target matrix
    :param _min: the smallest allowable output value
    :param _max: the largest allowable output value
    :return: the rescaled matrix (matrix is changed in-place)
    """
    M[:] = (M - M.min()) / (M.max() - M.min())
    M *= _max - _min
    M += _min
    return M


def similarity_to_distance(M, method, alpha=2, beta=1, headroom=1, verbose=False):
    """
    Convert a matrix representing similarity (larger values indicate increasing association)
    to a distance matrix (or dissimilarity matrix) where larger values indicate decreasing
    association (further away).

    As LKH requires distance as 32bit integers, the transformation must control the size of
    the largest value and make good use of the available range. Therefore similarity zeros
    (which translate to be the largest distances) are constrained to be only a factor
    "alpha" worse than the most distant (originally non-zero) value, while also being
    smaller than a fixed limit imposed by LKH. Exceeding this limit causes LKH to simply
    truncate those values to the limit -- a potential undesirable/unpredictable outcome.

    There are three transformation functions from which to choose:

    "inverse": y =  (1/x)^beta
     "neglog": y = -(log x/xmax)^beta
     "linear": y =  (1 - x/xmax)^beta

    All three functions display different treatment of small and large values.

    :param M: the target similarity matrix
    :param method: "inverse", "neglog" or "linear"
    :param alpha: factor to which similarity zeros are set in the distance beyond the largest distance.
    :param beta: an exponent to raise each element (default =1 ie. no effect)
    :param headroom: further factor of constraint to impose on largest integer allowed
    :param verbose: be verbose
    :return: distance matrix
    """

    assert alpha >= 1, 'alpha cannot be less than 1'
    INT_MAX = 2147483647.0
    largest = INT_MAX / 2.0 / len(M) / headroom
    if verbose:
        print 'Largest available integer: {:d}'.format(int(largest))

    # copy input matrix and remove diagonal
    M = M.astype(np.float)
    np.fill_diagonal(M, 0)

    # remember where zeros were
    zeros = (M == 0)
    if verbose:
        print 'Zero count:', zeros.sum()
        print 'Initial non-zero range: {:.3e} {:.3e}'.format(M[np.where(~zeros)].min(), M.max())

    nzix = np.where(~zeros)

    # transform similarity to distance, avoiding div-zero in some cases
    if method == 'inverse':
        M[nzix] = 1.0/M[nzix]
    elif method == 'linear':
        c = 1.0/M.max()
        M = 1.0 - c * M
    elif method == 'neglog':
        c = 1.0/M.max()
        M[nzix] = - np.log(c * M[nzix])
    else:
        raise RuntimeError('unsupported method: {}'.format(method))

    # apply element-wise power if requested
    if beta != 1:
        M = np.power(M, beta)

    # assign zeros (no observations) as a 'worst-case'
    maxM = M.max()
    if verbose:
        print 'Transformed range: {:.3e} {:.3e}'.format(M[np.where(~zeros)].min(), maxM)

    M[np.where(zeros)] = alpha * maxM
    if verbose:
        print 'Zeros assigned worst case of: {:.3e}'.format(alpha * maxM)

    # rescale to use available integer range
    M = scale_mat(M, 1, largest)
    if verbose:
        print 'Rescaled range: {:.3e} {:.3e}'.format(M.min(), M.max())

    return M


def lkh_order(m, base_name, precision=1, lkh_exe=None, runs=None, seed=None, dist_func=reciprocal_counts,
              fixed_edges=None, special=True, pop_size=None, stdout=None, verbose=False):
    """
    Employ LKH TSP solver to find the best order through a distance matrix. By default, it is assumed that
    LKH is on the path. A CalledProcessError is raised if execution fails. The input to LKH is an explicit
    definition of the full connected distance matrix, therefore sparse matrices will be converted to dense
    representations. For large problems, this can be memory demanding.

    :param m: the distance matrix
    :param base_name: base name of LKH control files
    :param precision: LKH internal precision factor, larger values limit maximum representable number
    :param lkh_exe: Path to binary, otherwise assumed on the path
    :param runs: the number of runs to perform LKH
    :param seed: random seed (milli-time if unspecified)
    :param dist_func: a custom distance function with which to convert matrix m.
    :param fixed_edges: list of edge tuples (u,v) that _must_ occur in the tour
    :param special: use LKH "special" meta-setting
    :param pop_size: population size of tours used in special genetic algorithm component (default: runs/4)
    :param stdout: redirection for stdout of lkh
    :param verbose: LKH will produce runtime information to stdout
    :return: 0-based order as a numpy array
    """
    if sparse.isspmatrix(m):
        m = m.toarray()
    m = dist_func(m.astype(np.float))

    try:
        write_lkh(base_name, m, len(m), max_trials=2*len(m), runs=runs, seed=seed, fixed_edges=fixed_edges,
                  pop_size=pop_size, special=special, verbose=verbose, precision=precision)
        if not lkh_exe:
            lkh_exe = 'LKH'
        subprocess.check_call([lkh_exe, '{}.par'.format(base_name)], stdout=stdout, stderr=subprocess.STDOUT)
        tour = read_lkh('{}.tour'.format(base_name))
    except subprocess.CalledProcessError as e:
        print 'Execution of LHK failed'
        raise e

    return tour['path']


def write_lkh(base_name, m, dim, max_trials=None, runs=None, verbose=False, seed=None, mat_fmt='upper',
              fixed_edges=None, pop_size=None, special=True, precision=1):
    """
    Create the control (.par) and data file (.dat) for the LKH executable. The implementation
    has many additional control parameters which could be included. Refer to LKH-3 documentation.
    :param base_name: base lkh file name
    :param m: the data (2D distance matrix or edges)
    :param dim: the number of nodes (cities)
    :param max_trials: maximum number of trials (default = dim)
    :param runs: number of runs (default = 10)
    :param verbose: enable stdout debug trace
    :param seed: random seed for algorithm (default milliseconds)
    :param mat_fmt: matrix format
    :param fixed_edges: list of edge tuples (u,v) that _must_ occur in the tour
    :param special: use LKH "special" meta-setting
    :param precision: scale integer values
    :param pop_size: population size of tours used in special genetic algorithm component (default: runs/4)
    """

    def write_full_matrix(out_h, m):
        out_h.write('EDGE_WEIGHT_TYPE: EXPLICIT\n')
        out_h.write('EDGE_WEIGHT_FORMAT: FULL_MATRIX\n')
        out_h.write('EDGE_WEIGHT_SECTION\n')
        np.savetxt(out_h, m, fmt='%d')

    def write_upper_row(out_h, m):
        out_h.write('EDGE_WEIGHT_TYPE: EXPLICIT\n')
        out_h.write('EDGE_WEIGHT_FORMAT: UPPER_ROW\n')
        out_h.write('EDGE_WEIGHT_SECTION\n')
        for i in xrange(len(m)-1):
            out_h.write(' '.join([str(int(vi)) for vi in m[i, i+1:]]))
            out_h.write('\n')

    assert isinstance(m, np.ndarray), 'the matrix must be a numpy array'
    if not seed:
        import time
        seed = int(round(time.time() * 1000))
    else:
        assert isinstance(seed, int), 'random seed must be an integer'

    nopath = os.path.basename(base_name)
    control_file = '{}.par'.format(base_name)
    data_file = '{}.dat'.format(base_name)
    with open(control_file, 'w') as out_h:
        if special:
            # SPECIAL is a meta setting for the following
            # out_h.write('GAIN23 = NO\n')
            # out_h.write('KICKS = 1\n')
            # out_h.write('KICK_TYPE = 4\n')
            # out_h.write('MAX_SWAPS = 0\n')
            # out_h.write('MOVE_TYPE = 5 SPECIAL\n')
            out_h.write('SPECIAL\n')
        if pop_size:
            out_h.write('POPULATION_SIZE = {}\n'.format(pop_size))
        out_h.write('PROBLEM_FILE = {}\n'.format(data_file))
        if max_trials:
            out_h.write('MAX_TRIALS = {}\n'.format(max_trials))
        if runs:
            out_h.write('RUNS = {}\n'.format(runs))
        out_h.write('SEED = {}\n'.format(seed))
        out_h.write('OUTPUT_TOUR_FILE = {}.tour\n'.format(base_name))
        out_h.write('PRECISION = {}.tour\n'.format(precision))
        out_h.write('TRACE_LEVEL = {}'.format(int(verbose)))

    with open(data_file, 'w') as out_h:
        out_h.write('NAME: {}\n'.format(nopath))
        out_h.write('TYPE: TSP\n')
        out_h.write('DIMENSION: {}\n'.format(dim))
        out_h.write('SALESMEN: 1\n')

        # currently only supporting one explicit data type
        if mat_fmt == 'full':
            write_full_matrix(out_h, m)
        elif mat_fmt == 'upper':
            write_upper_row(out_h, m)

        if fixed_edges:
            out_h.write('FIXED_EDGES_SECTION:\n')
            for u, v in fixed_edges:
                out_h.write('{} {}\n'.format(u, v))
            out_h.write('-1\n')


def read_lkh(fname):
    """
    Read the resulting solution (output tour) file from LKH.
    :param fname: the solution file name
    :return: dict of the tour information
    """
    tour = {'path': []}
    with open(fname, 'r') as in_h:
        in_tour = False
        for line in in_h:
            line = line.strip()
            if not line or line == 'EOF':
                break
            if line.startswith('TOUR_SECTION'):
                in_tour = True
            elif not in_tour:
                m = re.search('(\w+)[\s:=]+(\S+)', line)
                if not m:
                    continue
                if m.group(1) == 'DIMENSION':
                    tour[m.group(1)] = int(m.group(2))
                else:
                    tour[m.group(1)] = m.group(2)
            else:
                tour.setdefault('path', []).append(int(line))

        # convert ids to 0-based and remove end marker
        tour['path'] = np.array(tour['path'][:-1]) - 1

    return tour
