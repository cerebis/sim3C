#!/usr/bin/env python
from Bio import SeqIO
import truthtable as ttable
import argparse
import os
import sys


parser = argparse.ArgumentParser(description='Extract clustered sequences')
parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Verbose standard output')
parser.add_argument('-f', '--force', default=False, action='store_true', help='Force overwriting of files')
parser.add_argument('--cid-list', nargs='*', help='Specific cluster IDs')
parser.add_argument('--clustering', required=True, metavar='CLUSTERING', help='MCL format clustering file')
parser.add_argument('--fasta', required=True, metavar='FASTA', help='Fasta sequences supporting clustering')
parser.add_argument('-o', '--out-dir', help='Output directory')
args = parser.parse_args()

if args.out_dir:
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    elif not os.path.isdir(args.out_dir):
        raise IOError('Output path exists and is not a directory')

if args.verbose:
    print 'Getting sequence index'
seq_index = SeqIO.index(args.fasta, 'fasta')

if args.verbose:
    print 'Reading clustering'
tt = ttable.read_mcl(args.clustering)
cl2seq = tt.invert()


if len(args.cid_list) > 0:
    cid_list = args.cid_list
else:
    cid_list = [ci for ci in cl2seq]

for ci in cid_list:

    if args.verbose:
        print 'Collecting sequences for cluster {0}'.format(ci)

    seqs = []
    for si in cl2seq[int(ci)]:
        seqs.append(seq_index[si])

    if args.verbose:
        print 'Writing {0} sequences for cluster {1}'.format(len(seqs), ci)

    opath = os.path.join(args.out_dir, 'cl{0}.fasta'.format(ci))
    if not args.force and os.path.exists(opath):
        raise IOError('Path {0} already exists. Use --force to overwrite destination files'.format(opath))
    SeqIO.write(seqs, opath, 'fasta')
