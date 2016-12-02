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
import argparse

from dendropy import Tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rescale tree height')
    parser.add_argument('--max-height', type=float, metavar='FLOAT', default=0.1,
                        help='Scale longest branch to max height [0.1]')
    parser.add_argument('--if', dest='input_format', default='newick', choices=['newick', 'nexus', 'nexml'],
                        help='Input tree format [newick]')
    parser.add_argument('--of', dest='output_format', default='newick', choices=['newick', 'nexus', 'nexml'],
                        help='output tree format [newick]')
    parser.add_argument('input', type=argparse.FileType('r'), default='-',
                        help='Input tree')
    parser.add_argument('output', type=argparse.FileType('w'), default='-',
                        nargs='?', help='Output tree [stdout]')
    args = parser.parse_args()

    tr = Tree.get(file=args.input, schema=args.input_format)
    tr.scale_edges(args.max_height / max(tr.calc_node_root_distances()))
    tr.write_to_stream(args.output, args.output_format)
