#!/usr/bin/env python

import numpy as np
import networkx as nx
import mapio
import argparse


parser = argparse.ArgumentParser('H5 contact map to graph')
parser.add_argument('-f', '--fmt', default='edges', choices=['edges','graphml'], help='Output format [edges]')
parser.add_argument('h5_in', help='Input H5 contact map')
parser.add_argument('out_base', help='Output base file name')
args = parser.parse_args()

cm, n = mapio.read_map(args.h5_in, fmt='h5', has_names=True, make_dense=True)

g = nx.from_numpy_matrix(cm)
nx.relabel_nodes(g, {i: v for i, v in enumerate(n)}, copy=False)

print 'The {0}x{0} contact map transformed into the following graph'.format(len(cm))
print nx.info(g)

if args.fmt == 'edges':
    nx.write_edgelist(g, '{}.edgelist'.format(args.out_base), data=['weight']   )
else:
    nx.write_graphml(g, '{}.graphml'.format(args.out_base))

