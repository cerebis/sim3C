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
from Bio import Phylo
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Apply a global scale factor to all branch lengths in a newick tree')
    parser.add_argument('-a', '--alpha', type=float, required=True, help='Scale factor to apply')
    parser.add_argument('input', metavar='NEWICK', help='Input tree to scale')
    parser.add_argument('output', metavar='NEWICK', help='Scaled output tree')
    args = parser.parse_args()

    # read a tree in newick format. This could obviously be
    # expanded but would require considering output format
    tree = Phylo.read(args.input, 'newick')

    # iterate over clades, scale existing branch lengths
    # at least the root will have no length
    for clade_i in tree.find_clades():
        if clade_i.branch_length:
            clade_i.branch_length *= args.alpha

    # output branch length explicitly set as very small branch lengths
    # will be written as zero otherwise. May not be supported by some
    # downstream tools.
    Phylo.write(tree, args.output, 'newick', format_branch_length="%e")
