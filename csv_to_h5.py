#!/usr/bin/env python
import mapio
import argparse

parser = argparse.ArgumentParser(description='Convert CSV maps to H5')
parser.add_argument('-n', '--no-names', default=False, action='store_true', help='No row and column headers')
parser.add_argument('csv_in', help='CSV input file')
parser.add_argument('h5_out', help='H5 output file')
args = parser.parse_args()

print 'Reading CSV map...'
cmap, names = mapio.read_map(args.csv_in, fmt='csv', has_names=not args.no_names)

print 'Writing H5 map...'
mapio.write_map(cmap, args.h5_out, fmt='h5', names=None if args.no_names else names)
