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
import numpy as np


def calculate_N50_and_L50(lengths):
    """
    Calculate both N50 and L50 from a given collection of sequence lengths, taken from
    an assembly outcome.

    :param lengths: sequence lengths representing a set of assembly contigs
    :return: n50, l60 dict
    """
    l = np.array(lengths)
    l.sort()
    cl = np.cumsum(l)
    xi = np.searchsorted(cl, 0.5*l.sum(), side='right')
    return {'N50': (int)(l[xi]), 'L50': (int)(len(l) - xi)}

if __name__ == '__main__':
    import argparse

    def formatter(prog): return argparse.HelpFormatter(prog, width=100, max_help_position=100)
    parser = argparse.ArgumentParser(description='Calculate N50 and L50 from a set of sequences',
                                     formatter_class=formatter)
    parser.add_argument('--if', dest='seq_fmt', choices=['fasta', 'fastq'], default='fasta',
                        help='Input sequence format [fasta]')
    parser.add_argument('--yaml', action='store_true', default=False,
                        help='Write YAML format')
    parser.add_argument('input', help='Input sequence file')
    parser.add_argument('output', nargs="?", type=argparse.FileType('w'), default='-',
                        help='Output file')
    args = parser.parse_args()

    # determine sequence lengths from input file
    lengths = [len(si) for si in SeqIO.parse(args.input, args.seq_fmt)]

    # report N50 and L50
    stats = calculate_N50_and_L50(lengths)

    if args.yaml:
        import yaml
        yaml.dump(stats, args.output, default_flow_style=False)
    else:
        args.output.write('{N50}\t{L50}\n'.format(**stats))


