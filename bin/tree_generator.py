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
import dendropy
import argparse
import random
import sys


def simulate_tree(seed, ntax, max_height, birth_rate, death_rate):
    """
    Simulate a phylogenetic tree using a birth death model
    :param seed: random seed
    :param ntax: number of tips/taxa in resulting tree
    :param max_height: maximum root->tip length
    :param birth_rate: species birth rate
    :param death_rate: speices death rate
    :return: dendropy.Tree object
    """
    random.seed(seed)
    tns = dendropy.TaxonNamespace(map(lambda x: 't{0}'.format(x), range(1, ntax + 1)), label="taxa1")
    tr = dendropy.simulate.treesim.birth_death_tree(birth_rate, death_rate, taxon_namespace=tns, rng=random)
    tr.scale_edges(max_height / max(tr.calc_node_root_distances()))
    return tr


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate random phylogenetic trees')
    parser.add_argument('--suppress-rooting', default=False, action='store_true',
                        help='Suppress printing of rooting token [&R]|[&U]')
    parser.add_argument('--seed', type=int, metavar='INT',
                        help='Random seed')
    parser.add_argument('--max_height', type=float, metavar='FLOAT', default=0.1,
                        help='Maximum height of resulting tree [0.1]')
    parser.add_argument('--birth-rate', type=float, metavar='FLOAT', default=1.0,
                        help='Birth rate [1.0]')
    parser.add_argument('--death-rate', type=float, metavar='FLOAT', default=0.5,
                        help='Death rate [0.5]')
    parser.add_argument('--format', default='newick', choices=['newick', 'nexus', 'nexml'],
                        help='Output tree format [newick]')
    parser.add_argument('--num-taxa', metavar='INT', type=int, required=True,
                        help='Number of taxa in tree')
    parser.add_argument('output', type=argparse.FileType('w'), default=sys.stdout,
                        nargs='?', help='Output file name [stdout]')
    args = parser.parse_args()

    tr = simulate_tree(args.seed, args.num_taxa, args.max_height, args.birth_rate, args.death_rate)
    tr.write_to_stream(args.output, args.format, suppress_rooting=args.suppress_rooting)
