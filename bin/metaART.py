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

import argparse
import os
import subprocess
import sys
import gzip
import bz2
import numpy
from numpy import random

TMP_INPUT = 'seq.tmp'
TMP_OUTPUT = 'reads.tmp'

def open_output(fname, compress=None):
    if compress == 'bzip2':
        fh = bz2.BZ2File(fname + '.bz2', 'w')
    elif compress == 'gzip':
        # fix compression level to 6 since this is the norm on Unix. The default
        # of 9 is slow and is still often worse than bzip2.
        fh = gzip.GzipFile(fname + '.gz', 'w', compresslevel=6)
    else:
        fh = open(fname, 'w')
    return fh


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Simulate a metagenomic data set from an abundance profile')
    parser.add_argument('-C', '--compress', choices=['gzip', 'bzip2'], default=None, help='Compress output files')
    parser.add_argument('-n', '--output-name', metavar='PATH', help='Output file base name', required=True)
    parser.add_argument('-t', '--community-table', dest='comm_table', required=False,
                        help='Community profile table', metavar='FILE')
    parser.add_argument('-M', '--max-coverage', metavar='INT', type=int, required=True,
                        help='Coverage of must abundant taxon')
    parser.add_argument('-S', '--seed', metavar='INT', type=int, default='1234', required=True, help='Random seed')
    parser.add_argument('-l', '--read-len', metavar='INT', type=int, required=True, help='Read length')
    parser.add_argument('-m', '--insert-len', metavar='INT', type=int, required=True, help='Insert length')
    parser.add_argument('-s', '--insert-sd', metavar='INT', type=int, required=True, help='Insert standard deviation')
    parser.add_argument('--art-path', default='art_illumina', help='Path to ART executable [default: art_illumina]')
    parser.add_argument('--log', default='metaART.log', type=argparse.FileType('w'), help='Log file name')
    parser.add_argument('-z', '--num-samples', metavar='INT', type=int, default='1', required=True, help='Number of transect samples')
    parser.add_argument('-U', '--lognorm-ra-mu', metavar='FLOAT', type=float, default='1', required=False, help='Lognormal relative abundance mu parameter')
    parser.add_argument('-u', '--lognorm-ra-sigma', metavar='FLOAT', type=float, default='1', required=False, help='Lognormal relative abundance sigma parameter')
    parser.add_argument('fasta', metavar='MULTIFASTA',
                        help='Input multi-fasta of all sequences')
    parser.add_argument('output_dir', metavar='DIR',
                        help='Output directory')
    args = parser.parse_args()

    r1_tmp = os.path.join(args.output_dir, '{0}1.fq'.format(TMP_OUTPUT))
    r2_tmp = os.path.join(args.output_dir, '{0}2.fq'.format(TMP_OUTPUT))
    seq_tmp = os.path.join(args.output_dir, TMP_INPUT)

    seq_index = SeqIO.index(args.fasta, 'fasta')
    base_name = os.path.join(args.output_dir, args.output_name)
    all_R1 = open_output('{0}.r1.fq'.format(base_name), args.compress)
    all_R2 = open_output('{0}.r2.fq'.format(base_name), args.compress)

    # generate N simulated communities
    for n in range(0,args.num_samples):
        profile = {}
    	if args.num_samples == 1 and len(args.comm_table)>0:
            with open(args.comm_table, 'r') as h_table:
                for line in h_table:
                    line = line.rstrip().lstrip()
                    if line.startswith('#') or len(line) == 0:
                        continue
                    field = line.split()
                    if len(field) != 3:
                        print 'sequence table has missing fields at [', line, ']'
                        sys.exit(1)
                    profile[field[0]] = float(field[2])
        else:
            ra_sum = 0
            for seq_id in seq_index:
                profile[seq_id] = numpy.random.lognormal(args.lognorm_ra_mu, args.lognorm_ra_sigma)
                ra_sum += profile[seq_id]

            for seq_id in seq_index:
                profile[seq_id] /= ra_sum
            print "Sample " + str(n) + " rel abundances " + ", ".join(map(str, profile))


        r1_final = '{0}.{1}.r1.fq'.format(base_name,n)
        r2_final = '{0}.{1}.r2.fq'.format(base_name,n)
        r1_tmp = os.path.join(args.output_dir, '{0}1.fq'.format(TMP_OUTPUT,n)) 
        r2_tmp = os.path.join(args.output_dir, '{0}2.fq'.format(TMP_OUTPUT,n))

        with open_output(r1_final, args.compress) as output_R1, open_output(r2_final, args.compress) as output_R2:
            try:
                for seq_id in profile:

                    coverage = float(profile[seq_id] * args.max_coverage)
                    print 'Requesting {0} coverage for {1}'.format(coverage, seq_id)

                    ref_seq = seq_index[seq_id]

                    ref_len = len(ref_seq)
                    SeqIO.write([ref_seq], seq_tmp, 'fasta')

                    try:
                        subprocess.call([args.art_path,
                                     '-p',   # paired-end sequencing
                                     '-na',  # no alignment file
                                     '-rs', str(args.seed),
                                     '-m', str(args.insert_len),
                                     '-s', str(args.insert_sd),
                                     '-l', str(args.read_len),
                                     '-f', str(coverage),
                                     '-i', seq_tmp,
                                     '-o', os.path.join(args.output_dir, TMP_OUTPUT)], stdout=args.log, stderr=args.log)
                    except OSError as ex:
                        print "There was an error starting the art_illumina subprocess."
                        print "You may need to add its location to your PATH or specify it at runtime."
                        raise ex

                    # count generated reads
                    r1_n = 0
                    for seq in SeqIO.parse(r1_tmp, 'fastq'):
                        r1_n += 1

                    r2_n = 0
                    for seq in SeqIO.parse(r2_tmp, 'fastq'):
                        r2_n += 1

                    effective_cov = args.read_len * (r1_n + r2_n) / float(ref_len)
                    print 'Generated {0} paired-end reads for {1}, {2:.3f} coverage'.format(r1_n, seq_id, effective_cov)
                    if r1_n != r2_n:
                        print 'Error: paired-end counts do not match {0} vs {1}'.format(r1_n, r2_n)
                        sys.exit(1)

                    with open(r1_tmp, 'r') as tmp_h:
                        all_R1.write(tmp_h.read())

                    with open(r2_tmp, 'r') as tmp_h:
                        all_R2.write(tmp_h.read())

                    with open(r1_tmp, 'r') as tmp_h:
                        output_R1.write(tmp_h.read())
                        os.remove(tmp_h.name)

                    with open(r2_tmp, 'r') as tmp_h:
                        output_R2.write(tmp_h.read())
                        os.remove(tmp_h.name)

                    os.remove(seq_tmp)
            finally:
                print "all done, let's go home."
