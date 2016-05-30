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

import networkx as nx
import community as com
import numpy as np
import argparse
import sys
import os


def formatter(prog): return argparse.HelpFormatter(prog, width=100, max_help_position=100)

parser = argparse.ArgumentParser(description='Graph statistics', formatter_class=formatter)
parser.add_argument('--yaml', action='store_true', default=False,
                    help='Write YAML format')
parser.add_argument('--excl-self', action='store_true', default=False,
                    help='Exclude self-loops')
parser.add_argument('input', help='GraphML format graph file to analyse')
parser.add_argument('output', nargs='?', type=argparse.FileType('w'), default='-',
                    help='Output file')
args = parser.parse_args()

print 'Reading graph ...'
g = nx.read_graphml(args.input)

if g.order() == 0:
    print 'The graph contains no nodes'
    sys.exit(1)

if args.excl_self:
    # no self-loops
    print 'Excluding self-loops ...'
    g.remove_edges_from(g.selfloop_edges())

result = {
    'size': None,
    'order': None,
    'density': None,
    'isolates': None,
    'inclusives': None,
    'mean_deg': None,
    'median_deg': None,
    'modularity': None
}

try:

    result['size'] = g.size()
    result['order'] = g.order()

    # determine Louvain modularity score of entire graph
    if g.size() == 0:
        # no edges, modularity is 1
        result['modularity'] = 1.0
        result['density'] = 0
    else:
        part = com.best_partition(g)
        result['modularity'] = com.modularity(part, g)
        result['density'] = nx.density(g) if g.order() > 0 else 'NA'

except Exception as e:
    print 'An exception occurred during modularity analysis'
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print exc_type, fname, exc_tb.tb_lineno
    sys.exit(1)

# count isolates
result['isolates'] = len(nx.isolates(g))
result['inclusives'] = g.order() - result['isolates']

# mean degree
result['mean_deg'] = np.mean(nx.degree(g).values())

# median degree
result['median_deg'] = np.median(nx.degree(g).values())

print '#file #nodes #iso #inclu Modularity MeanDeg MedianDeg'
print '{0} {1} {2} {3} {4} {5} {6}'.format(sys.argv[1], g.order(), n_iso, g.order() - n_iso, mod, mean_deg, med_deg)

if args.yaml:
    import yaml
    yaml.dump(result, args.output, default_flow_style=False)
else:
    args.output.write('{0}\n'.format('\t'.join(result.keys())))
    args.output.write('{0}\n'.format('\t'.join([str(v) for v in result.values()])))
