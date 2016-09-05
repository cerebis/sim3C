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

import numpy
from Bio import SeqIO

import abundance

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate random abundance profiles for a given multifasta file')
    parser.add_argument('--seed', type=int, metavar='INT', help='Random seed')
    parser.add_argument('--dist', metavar='DISTNAME', choices=['equal', 'uniform', 'lognormal'], required=True,
                        help='Abundance profile distribution [equal, uniform, lognormal]')
    parser.add_argument('--lognorm-mu', metavar='FLOAT', type=float, default='1', required=False,
                        help='Log-normal relative abundance mu parameter')
    parser.add_argument('--lognorm-sigma', metavar='FLOAT', type=float, default='1', required=False,
                        help='Log-normal relative abundance sigma parameter')
    parser.add_argument('input', metavar='FASTA', help='Input fasta sequence')
    parser.add_argument('output', metavar='TABLE', help='Output file name')
    args = parser.parse_args()

    RANDOM_STATE = None
    if args.dist != 'equal' and not args.seed:
        print 'Warning: not specifying a seed makes repeatability difficult!'
        RANDOM_STATE = numpy.random.RandomState()
    else:
        RANDOM_STATE = numpy.random.RandomState(args.seed)

    seq_index = None
    try:
        seq_index = SeqIO.index(args.input, 'fasta')
        if len(seq_index) <= 0:
            raise IOError('Input file contained no sequences')

        profile = abundance.generate_profile(RANDOM_STATE, list(seq_index), mode=args.dist,
                                             lognorm_mu=args.lognorm_mu, lognorm_sigma=args.lognorm_sigma)

        with open(args.output, 'w') as out_h:
            profile.write_table(out_h)

    finally:
        if seq_index:
            seq_index.close()
