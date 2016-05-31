#!/usr/bin/env python
import yaml
import argparse
import json

parser = argparse.ArgumentParser(description='Aggregate analysis results')
parser.add_argument('gstats', type=argparse.FileType('r'), help='Yaml output from graph_stats.py')
parser.add_argument('complex', type=argparse.FileType('r'), help='Yaml output from graph_complexity.py')
parser.add_argument('bcubed', type=argparse.FileType('r'), help='Yaml output from bcubed.py')
args = parser.parse_args()

d = {}
d.update(yaml.load(args.gstats))
d.update(yaml.load(args.complex))
d.update(yaml.load(args.bcubed))

ksrt = sorted(d.keys())

print json.dumps(d)
