#!/usr/bin/env python

import numpy as np
import networkx as nx
import argparse
import zlib

def kolmogorov(s):
    """
    Approximate a measure of Kolmogorov complexity by using lossless compression
    to compare input and output of a string. Here, a graph can be represented
    as a linearised adjacency matrix.

    :param s: a string
    :return: measure of complexity as defined by the degree to which a string can be compressed without loss of
    information [0..1].
    """
    if len(s) < 100:
        raise RuntimeError('Strings shorter than 100 bytes are not effectively '
                           'compressed by zlib, therefore we cannot proceed')
    s_compressed = zlib.compress(s, 9)
    return float(len(s_compressed))/float(len(s))


def find_nonzero_eigenvalues_magnitudes(g, edge_weight='weight'):
    print '  Calculating normalized laplacian ...'
    L = nx.normalized_laplacian_matrix(g, weight=edge_weight)
    print '  Calculating eigenvalues of laplacian ...'
    ev_i = np.absolute(np.linalg.eigvals(L.A.astype(np.float64)))
    nz_ix = ev_i != 0
    if np.sum(nz_ix) <= 0:
        raise RuntimeError('All eigenvalues zero for the Laplacian of supplied graph')
    return ev_i[nz_ix]


def eigen_product(g, s=1.0, edge_weight='weight'):
    """
    Calculate the product of eigenvalues from a graph G. This has been used as a discriminating
    function between graphs.
    :param g: the target graph
    :param s: exponential scale factor
    :param edge_weight: data field used for edge weights. If None, weights =1.
    :return: product eigenvalues
    """
    nz_eig = find_nonzero_eigenvalues_magnitudes(g, edge_weight)
    return np.product(nz_eig ** (1.0/s))


def inverse_eigen_product(g, s=1.0, edge_weight='weight'):
    """
    Calculate the inverse product of eigenvalues from a graph G. This has been used as a discriminating
    function between graphs.
    :param g: the target graph
    :param s: exponential scale factor
    :param edge_weight: data field used for edge weights. If None, weights =1.
    :return: product eigenvalues
    """
    return 1.0/eigen_product(g, s, edge_weight)


def eigen_entropy(g, s=1.0, edge_weight='weight'):
    """
    Calcualte eigenvalue based entropy from a graph G. The eigenvalues are those from the
    normalized laplacian. Log base 2.
    :param s: exponential scale factor
    :param g: the target graph
    :param edge_weight: data field used for edge weights. If None, weights =1.
    :return: entropy value
    """
    ev_i = find_nonzero_eigenvalues_magnitudes(g, edge_weight)
    _pi = ev_i**(1.0/s)
    _pi = _pi / np.sum(_pi)
    return -np.sum(_pi * np.log2(_pi))

parser = argparse.ArgumentParser(description='Graph complexity estimation')
parser.add_argument('-s', '--scale', type=float, default=1.0, help='Exponential scale factor [1.0]')
parser.add_argument('-w', '--weight', default='weight', help='Edge weighting field name [weight]')
parser.add_argument('-m', '--method', choices=['kolmo', 'eigh', 'eigp', 'eigip'], required=True,
                    help='Method to apply')
parser.add_argument('input', help='GraphML format graph file to analyse')
args = parser.parse_args()

print '  Reading graph {0}'.format(args.input)
g = nx.read_graphml(args.input)
print '  Graph contained {0} nodes and {1} edges'.format(g.order(), g.size())

# set weight parameters for edges.
# this is necessary if we want the default to be weighted but still provide
# users the option for equal weighting.
weight = None if args.weight.lower() == 'none' else args.weight

if args.method == 'kolmo':
    print '  Determining adjacency matrix ...'
    g_adj = nx.adjacency_matrix(g)
    g_adj = g_adj.todense()

    # convert the adjacency matrix (which might be weighted) into a binary string.
    # where 0 no connection, 1 connection with weight > 0.
    print '  Converting adjacency matrix to binary string representation ...'
    g_str = ''.join(map(str, g_adj.flatten().astype(bool).astype(int).tolist()[0]))

    try:
        k = kolmogorov(g_str)
    except RuntimeError as er:
        print er
        k = None
    print 'Kolmogorov complexity estimate: {0}'.format(k)

elif args.method == 'eigh':
    try:
        H = eigen_entropy(g, args.scale, weight)
    except RuntimeError as er:
        print er
        H = None
    print 'Eigen value entropy: {0}'.format(H)

elif args.method == 'eigp':
    try:
        P = eigen_product(g, args.scale, weight)
    except RuntimeError as er:
        P = None
    print 'Eigen value entropy: {0}'.format(P)

elif args.method == 'eigip':
    try:
        IP = inverse_eigen_product(g, args.scale, weight)
    except RuntimeError as er:
        IP = None
    print 'Eigen value entropy: {0}'.format(IP)
