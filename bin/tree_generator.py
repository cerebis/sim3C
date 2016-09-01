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
from dendropy.simulate import treesim
import dendropy
import argparse
import random
import sys


def get_taxon_namespace(ntax, label='taxa', prefix='T'):
    """
    Generate a TaxonNamepsace for use in dendroppy method calls. This represents the names
    of taxa and determines the number of tips generated in a tree simulation.

    :param ntax: the number of taxa/tips
    :param label: a label for the namespace generated.
    :param prefix: a prefix prepending each taxon's name
    :return:
    """

    return dendropy.TaxonNamespace(map(lambda x: '{0}_{1}'.format(prefix, x), range(1, ntax + 1)), label=label)


def rescale_tree(tr, max_height):
    """
    Recale a tree by maximum root->tip height
    :param tr: the tree to scale
    :param max_height: maximum height from root to furtest tip
    :return: rescaled dendropy.Tree
    """
    tr.scale_edges(max_height / max(tr.calc_node_root_distances()))
    return tr


def star_tree(ntax, max_height):
    """
    Generate a star tree, this is a non-random process and so requires no seed.
    :param tns: taxon namespace
    :param max_height: maximum root->tip length
    :return: dendropy.Tree object
    """
    tr = treesim.star_tree(tns)
    # star trees have no length attributes, so we explicitly set them
    for e in tr.edges():
        if e.is_leaf():
            e.length = max_height
    return tr


def simulate_tree(seed, tns, max_height, birth_rate, death_rate):
    """
    Simulate a phylogenetic tree using a birth death model
    :param seed: random seed
    :param tns: taxon namespace
    :param max_height: maximum root->tip length
    :param birth_rate: species birth rate
    :param death_rate: speices death rate
    :return: dendropy.Tree object
    """
    if seed:
        random.seed(seed)

    tr = treesim.birth_death_tree(birth_rate, death_rate, taxon_namespace=tns, rng=random)
    return rescale_tree(tr, max_height)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate random phylogenetic trees')
    parser.add_argument('--suppress-rooting', default=False, action='store_true',
                        help='Suppress printing of rooting token [&R]|[&U]')
    parser.add_argument('--seed', type=int, metavar='INT', help='Random seed')
    parser.add_argument('--mode', choices=['star', 'random'], required=True,
                        help='Tree generation mode [star, random]')
    parser.add_argument('--max-height', type=float, metavar='FLOAT', default=0.1,
                        help='Maximum height of resulting tree [0.1]')
    parser.add_argument('--birth-rate', type=float, metavar='FLOAT', default=1.0,
                        help='Birth rate [1.0]')
    parser.add_argument('--death-rate', type=float, metavar='FLOAT', default=0.5,
                        help='Death rate [0.5]')
    parser.add_argument('--format', default='newick', choices=['newick', 'nexus', 'nexml'],
                        help='Output tree format [newick]')
    parser.add_argument('--prefix', metavar='PREFIX', default='T',
                        help='Taxon name prefix [T]')
    parser.add_argument('--num-taxa', metavar='INT', type=int, required=True,
                        help='Number of taxa in tree')
    parser.add_argument('output', type=argparse.FileType('w'), default=sys.stdout,
                        nargs='?', help='Output file name [stdout]')
    args = parser.parse_args()

    if args.mode == 'random' and not args.seed:
        print 'Warning: random tree without specifying a seed makes repeatability difficult!'

    tns = get_taxon_namespace(args.num_taxa, prefix=args.prefix)

    gen_tr = None
    if args.mode == 'random':
        gen_tr = simulate_tree(args.seed, tns, args.max_height, args.birth_rate, args.death_rate)
    elif args.mode == 'star':
        gen_tr = star_tree(tns, args.max_height)
    else:
        raise RuntimeError('unsupported mode')

    gen_tr.write_to_stream(args.output, args.format, suppress_rooting=args.suppress_rooting)
