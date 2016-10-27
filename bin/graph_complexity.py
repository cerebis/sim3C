#!/usr/bin/env python
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
import json
import zlib

import networkx as nx
import numpy as np
import yaml


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
    return float(np.product(nz_eig ** (1.0/s)))


def inverse_eigen_product(g, s=1.0, edge_weight='weight'):
    """
    Calculate the inverse product of eigenvalues from a graph G. This has been used as a discriminating
    function between graphs.
    :param g: the target graph
    :param s: exponential scale factor
    :param edge_weight: data field used for edge weights. If None, weights =1.
    :return: product eigenvalues
    """
    return float(1.0/eigen_product(g, s, edge_weight))


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
    return float(-np.sum(_pi * np.log2(_pi)))


if __name__ == '__main__':
    import argparse

    def write_output(out_h, d, fmt='json'):
        """
        Write result to stdout or file in specified format
        :param out_h: output stream
        :param d: result of computation
        :param fmt: output format yaml | json | plain
        """
        if fmt == 'yaml':
            yaml.dump(d, out_h, default_flow_style=False)
        elif fmt == 'json':
            json.dump(d, out_h, indent=1)
        elif fmt == 'plain':
            args.output.write('{method} {value}\n'.format(**d))

    def formatter(prog): return argparse.HelpFormatter(prog, width=100, max_help_position=100)

    parser = argparse.ArgumentParser(description='Graph complexity estimation', formatter_class=formatter)
    parser.add_argument('--ofmt', choices=['plain', 'json', 'yaml'], default='plain',
                        help='Output format [plain]')
    parser.add_argument('-s', '--scale', type=float, default=1.0,
                        help='Exponential scale factor [1.0]')
    parser.add_argument('-w', '--weight', default='weight',
                        help='Edge weighting field name [weight]')
    parser.add_argument('--method', choices=['kolmo', 'eigh', 'eigp', 'eigip'],
                        default='eigh', help='Method to apply [eigh]')
    parser.add_argument('input', help='GraphML format graph file to analyse')
    parser.add_argument('output', nargs='?', type=argparse.FileType('w'), default='-',
                        help='Output file')
    args = parser.parse_args()

    print '  Reading graph {0}'.format(args.input)
    g = nx.read_graphml(args.input)
    print '  Graph contained {0} nodes and {1} edges'.format(g.order(), g.size())

    # set weight parameters for edges.
    # this is necessary if we want the default to be weighted but still provide
    # users the option for equal weighting.
    weight = None if args.weight.lower() == 'none' else args.weight

    result = {'method': args.method, 'value': None}

    if args.method == 'kolmo':
        print 'Determining adjacency matrix ...'
        g_adj = nx.adjacency_matrix(g)
        g_adj = g_adj.todense()

        # convert the adjacency matrix (which might be weighted) into a binary string.
        # where 0 no connection, 1 connection with weight > 0.
        print 'Converting adjacency matrix to binary string representation ...'
        g_str = ''.join(map(str, g_adj.flatten().astype(bool).astype(int).tolist()[0]))

        try:
            result['value'] = kolmogorov(g_str)
        except RuntimeError as er:
            print er

    elif args.method == 'eigh':
        try:
            result['value'] = eigen_entropy(g, args.scale, weight)
        except RuntimeError as er:
            print er

    elif args.method == 'eigp':
        try:
            result['value'] = eigen_product(g, args.scale, weight)
        except RuntimeError as er:
            print er

    elif args.method == 'eigip':
        try:
            result['value'] = inverse_eigen_product(g, args.scale, weight)
        except RuntimeError as er:
            print er

    write_output(args.output, result, args.ofmt)
