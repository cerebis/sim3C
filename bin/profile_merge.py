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
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate random abundance profiles for a given multifasta file')
    parser.add_argument('input', metavar='FILE', nargs='*', help='Input abundance profile tables')
    parser.add_argument('output', metavar='FILE', help='Output merged table')
    args = parser.parse_args()

    HEADER = '\treplicon\tcell\tabundance'

    all_tables = []
    for fn in args.input:
        with open(fn, 'r') as in_h:
            print 'Reading {0}'.format(fn)
            for line in in_h:
                line = line.strip()
                if len(line) == 0 or line.startswith('#'):
                    continue
                tok = line.split('\t')
                all_tables.append(tok)

    # abundance value column -- normalised
    values = np.array([float(row[3]) for row in all_tables])
    values /= values.sum()
    values = values.reshape((8, 1))

    # middle columns
    middle = np.array(all_tables)[:, 1:3]

    # regenerate index
    index = np.arange(1, middle.shape[0]+1).reshape((8, 1))

    # paste the columns back together
    norm_table = np.hstack((index, middle, values))

    print 'Writing merged profile {0}'.format(args.output)
    np.savetxt(args.output, norm_table, fmt='%s', header=HEADER)
