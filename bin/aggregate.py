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
import io_utils
import yaml
import argparse
import json

parser = argparse.ArgumentParser(description='Aggregate serialized dictionaries')
parser.add_argument('--fmt', choices=['yaml', 'json'], default='json',
                    help='Serialized object format [json]')
parser.add_argument('input_files', metavar='FILE', type=argparse.FileType('r'),
                    nargs='+', help='Serialized data file')
args = parser.parse_args()

load_func = io_utils.json_load_byteified if args.fmt == 'json' else yaml.load
dump_func = json.dumps if args.fmt == 'json' else yaml.safe_dump

d = {}
for fi in args.input_files:
    d_frag = load_func(fi)
    ovl = set(d_frag) & set(d)
    if len(ovl) != 0:
        key_list = ', '.join([str(ki) for ki in ovl])
        raise RuntimeError('Overlapping keys exist between dictionaries: {0}'.format(key_list))
    d.update(d_frag)

print dump_func(d)
