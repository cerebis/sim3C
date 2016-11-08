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

import abundance

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate random abundance profiles for a given multifasta file')
    parser.add_argument('input', metavar='FILE', nargs='*', help='Input abundance profile tables')
    parser.add_argument('output', metavar='FILE', help='Output merged table')
    args = parser.parse_args()

    profile = abundance.Profile()
    for fn in args.input:
        with open(fn, 'r') as in_h:
            print 'Reading {0}'.format(fn)
            profile.update(abundance.read_profile(in_h))

    # re-establish relative weighting -- bring the sum of all abundances back down to 1.
    profile.normalize()

    print 'Writing merged profile {0}'.format(args.output)
    with open(args.output, 'w') as out_h:
        profile.write_table(out_h)
