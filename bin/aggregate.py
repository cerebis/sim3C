#!/usr/bin/env python
import pipeline_utils as pu
import yaml
import argparse
import json

parser = argparse.ArgumentParser(description='Aggregate serialized dictionaries')
parser.add_argument('--fmt', choices=['yaml', 'json'], default='json',
                    help='Serialized object format [json]')
parser.add_argument('input_files', metavar='FILE', type=argparse.FileType('r'),
                    nargs='+', help='Serialized data file')
args = parser.parse_args()

load_func = pu.json_load_byteified if args.fmt == 'json' else yaml.load
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
