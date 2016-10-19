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
from Bio import SeqIO
from Bio.Restriction import *
import networkx as nx
import numpy as np
import argparse
import sys

RESTRICTION_BATCH = None
ENZYME = None


def count_sites(seq, linear=True):
    """
    Digest a sequence and count the number of restriction sites found.
    :param seq: sequence to digest
    :param linear: treat the sequence as being linear
    :return: number of sites
    """
    return len(ENZYME.search(seq, linear=linear))


def scale_factor_from_length(g, u, v):
    """
    Use the harmonic mean of the lengths of U and V. Estimate the number of sites
    by assuming probability of a cutsite is dependent only on the length of a cutter's
    recognition site.
    :param g: graph containing u and v
    :param u: node u
    :param v: node v
    :return: scale factor
    """
    l_u = g.node[u]['length']
    l_v = g.node[v]['length']
    est_sites = RECIP_CUTFREQ * 2.0 * (l_u * l_v) / (l_u + l_v)
    return est_sites ** 2.0


def scale_factor_from_sites(g, u, v):
    """
    Scale weight linking U and V by the product of "number of recognition sites" found for each sequence.
    :param g: graph containing u and v
    :param u: node u
    :param v: node v
    :return: scale factor
    """
    s_u = g.node[u]['n_sites']
    s_v = g.node[v]['n_sites']
    return float(s_u * s_v)


parser = argparse.ArgumentParser(description='Filter and normalise graphml file from raw counts')
parser.add_argument('--min-length', type=int, help='Minimum sequence length')
parser.add_argument('--no-isolates', default=False, action='store_true', help='Remove isolated nodes')
parser.add_argument('--no-self', default=False, action='store_true', help='Remove self-loops')
parser.add_argument('-w', '--weight', default=0, type=int, help='Threshold raw edge weight to exclude.[0]')
parser.add_argument('--percentile', default=None, type=float,
                    help='Percentile threshold weight below which to exclude [0..100]')
parser.add_argument('--cutter-length', type=int, help='Length of enzyme cut-site')
parser.add_argument('--cutter', help='Restriction enzyme used to product ligation products')
parser.add_argument('--fasta', help='Fasta file for corresponding node sequences')
parser.add_argument('graph', help='GraphML format graph file to analyse')
parser.add_argument('output', help='Output GraphML file')
args = parser.parse_args()


if args.cutter:
    try:
        RESTRICTION_BATCH = RestrictionBatch([args.cutter])
        ENZYME = RestrictionBatch.get(RESTRICTION_BATCH, args.cutter, add=False)
        if args.cutter_length:
            print 'Warning: option "cutter-length" is ignored if option "cutter" is supplied'
        args.cutter_length = ENZYME.size
    except ValueError:
        print 'Error: failed to find an enzyme named [{0}]'.format(args.cutter)
        sys.exit(1)

RECIP_CUTFREQ = 1.0 / 4 ** args.cutter_length

g = nx.read_graphml(args.graph)
print 'Raw graph contained'
print nx.info(g)

if args.no_self:
    g.remove_edges_from(g.selfloop_edges())
    print 'Removed self-loops'
    print nx.info(g)

if args.min_length:
    to_remove = [u for u, d in g.nodes_iter(data=True) if d['length'] > args.min_length]
    g.remove_nodes_from(to_remove)
    print 'After removing sequences short than {0} bp'.format(args.min_length)
    print nx.info(g)

if args.weight > 0:
    to_remove = []
    for u, v, d in g.edges_iter(data=True):
        if d['weight'] <= args.weight:
            to_remove.append((u, v))
    g.remove_edges_from(to_remove)
    print 'After removing edges <= {0} weight'.format(args.weight)
    print nx.info(g)

if args.fasta:
    if not args.cutter:
        print 'Warning: fasta sequences are ignored if no cutter enzyme is supplied'
    else:
        # for each node in the graph, add the number of sites found.
        seqidx = SeqIO.index(args.fasta, 'fasta')
        vacant_nodes = []
        for u, d in g.nodes_iter(data=True):
            n_sites = count_sites(seqidx[u].seq)
            d.update({'n_sites': n_sites})
            if n_sites == 0:
                vacant_nodes.append(u)

        # nodes with no detected cutsites are removed
        print 'Removing {0} nodes which had no cutsite'.format(len(vacant_nodes))
        print nx.info(g)
        g.remove_nodes_from(vacant_nodes)


print 'Rescaling edges'

# set the function of producing scale factors
scale_factor = scale_factor_from_sites if args.cutter else scale_factor_from_length

# Normalise raw weights
for u, v, d in g.edges_iter(data=True):
    d['weight'] /= scale_factor(g, u, v)

if args.percentile:
    weights = np.fromiter((d['weight'] for u, v, d in g.edges_iter(data=True)), dtype=float)
    thres = np.percentile(weights, args.percentile)
    print 'Thresholding graph at {0}%: {1:.3e}'.format(args.percentile, thres)
    to_remove = []
    for u, v, d in g.edges_iter(data=True):
        if d['weight'] < thres:
            to_remove.append((u, v))
    g.remove_edges_from(to_remove)
    print 'After thresholding'
    print nx.info(g)

if args.no_isolates:
    to_remove = []
    for u in g.nodes():
        if g.degree(u) == 0:
            to_remove.append(u)
    g.remove_nodes_from(to_remove)
    print 'After removing isolated nodes'
    print nx.info(g)

nx.write_graphml(g, args.output)
