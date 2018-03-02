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


def lkh_order(m, base_name, precision=1, lkh_exe=None, runs=None, seed=None):
    """
    Employ LKH TSP solver to find the best order through a distance matrix. By default, it is assumed that
    LKH is on the path. A CalledProcessError is raised if execution fails. The input to LKH is an explicit
    definition of the full connected distance matrix, therefore sparse matrices will be converted to dense
    representations. For large problems, this can be memory demanding.

    :param m: the distance matrix
    :param base_name: base name of LKH control files
    :param precision: LKH is limited to integer distances, apply this multiplier to distances prior to truncation.
    :param lkh_exe: Path to binary, otherwise assumed on the path
    :return: 0-based order as a numpy array
    """
    if sparse.isspmatrix(m):
        m = m.toarray()

    m = m.astype(np.float)

    # distance will be inversely proportional to counts
    # first, all elements are non-zero
    np.fill_diagonal(m, 0)
    m += 0.1
    # scale elements [0,1]
    m = m / m.max()
    # take inverse
    m = 1 / m
    # rescale, making sure the shortest distance is > 1
    # the reason for this is to protect the shortest distance elements
    # from going to zero when the matrix is converted to integers (requirement of LKH)
    m *= 1.01 / m.min()

    print m.min(), m.max()

    # remove self-self paths
    # TODO this should perhaps come before the above.
    np.fill_diagonal(m, 0)

    try:
        write_lkh(base_name, m*precision, len(m), max_trials=2*len(m), runs=runs, seed=seed)
        if not lkh_exe:
            lkh_exe = 'LKH'
        subprocess.check_call([lkh_exe, '{}.par'.format(base_name)])
        tour = read_lkh('{}.tour'.format(base_name))
    except subprocess.CalledProcessError as e:
        print 'Execution of LHK failed'
        raise e

    return tour['path']


def write_lkh(base_name, m, dim, max_trials=None, runs=None, trace=0, seed=None):
    """
    Create the control (.par) and data file (.dat) for the LKH executable. The implementation
    has many additional control parameters which could be included. Refer to LKH-3 documentation.
    :param base_name: base lkh file name
    :param m: the data (2D distance matrix or edges)
    :param dim: the number of nodes (cities)
    :param max_trials: maximum number of trials (default = dim)
    :param runs: number of runs (default = 10)
    :param trace: debug level (default = 0 quiet-ish)
    :param seed: random seed for algorithm (default milliseconds)
    """

    def write_full_matrix(out_h, fm):
        out_h.write('EDGE_WEIGHT_TYPE: EXPLICIT\n')
        out_h.write('EDGE_WEIGHT_FORMAT: FULL_MATRIX\n')
        out_h.write('EDGE_WEIGHT_SECTION\n')
        np.savetxt(out_h, fm, fmt='%d')

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
        out_h.write('SPECIAL\n')
        # SPECIAL is a meta setting for the following
        # out_h.write('GAIN23 = NO\n')
        # out_h.write('KICKS = 1\n')
        # out_h.write('KICK_TYPE = 4\n')
        # out_h.write('MAX_SWAPS = 0\n')
        # out_h.write('MOVE_TYPE = 5 SPECIAL\n')
        # out_h.write('POPULATION_SIZE = 10\n')
        out_h.write('PROBLEM_FILE = {}\n'.format(data_file))
        if max_trials:
            out_h.write('MAX_TRIALS = {}\n'.format(max_trials))
        if runs:
            out_h.write('RUNS = {}\n'.format(runs))
        out_h.write('SEED = {}\n'.format(seed))
        out_h.write('OUTPUT_TOUR_FILE = {}.tour\n'.format(nopath))
        out_h.write('TRACE_LEVEL = {}'.format(trace))

    with open(data_file, 'w') as out_h:
        out_h.write('NAME: {}\n'.format(nopath))
        out_h.write('TYPE: TSP\n')
        out_h.write('DIMENSION: {}\n'.format(dim))
        out_h.write('SALESMEN: 1\n')

        # currently only supporting one explicit data type
        write_full_matrix(out_h, m)


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
