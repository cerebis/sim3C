#!/bin/env python
import yaml
import pandas
import argparse

parser = argparse.ArgumentParser(description='Convert aggregated results from YAML object to flat csv table')
parser.add_argument('input', help='Input YAML aggregated results file')
parser.add_argument('output', nargs='?', default='-', type=argparse.FileType('w'), help='Output file [stdout]')
args = parser.parse_args()

with open(args.input, 'r') as input_h:
    results = yaml.load(input_h)
    header = ['{0}-{1}'.format(_type, _name) for _type in ['params','asmstat','gstat','geigh','bc'] for _name in results[0][_type].keys()]
    table = [[ _val for _type in ['params','asmstat','gstat','geigh','bc'] for _val in _row[_type].values()] for _row in results]
    pandas.DataFrame(table, columns=header).to_csv(args.output)

