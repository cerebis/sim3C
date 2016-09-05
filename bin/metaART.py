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
import atexit
import os
import subprocess
import sys

import numpy
from Bio import SeqIO

import abundance
import io_utils


TMP_INPUT = 'seq.tmp'
TMP_OUTPUT = 'reads.tmp'

# low-high seeds, giving 5M values
LOW_SEED_VALUE = 1000000
HIGH_SEED_VALUE = 6000000


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Simulate a metagenomic data set from an abundance profile')
    parser.add_argument('-C', '--compress', choices=['gzip', 'bzip2'], default=None, help='Compress output files')
    parser.add_argument('-n', '--output-name', metavar='PATH', help='Output file base name', required=True)
    parser.add_argument('-P', '--profile', dest='profile', required=False,
                        help='Community abundance profile', metavar='FILE')
    parser.add_argument('-M', '--max-coverage', metavar='INT', type=int, required=True,
                        help='Coverage of must abundant taxon')
    parser.add_argument('-S', '--seed', metavar='INT', type=int, required=True, help='Random seed')
    parser.add_argument('-l', '--read-len', metavar='INT', type=int, required=True, help='Read length')
    parser.add_argument('-m', '--insert-len', metavar='INT', type=int, required=True, help='Insert length')
    parser.add_argument('-s', '--insert-sd', metavar='INT', type=int, required=True, help='Insert standard deviation')
    parser.add_argument('--art-path', default='art_illumina', help='Path to ART executable [default: art_illumina]')
    parser.add_argument('--log', default='metaART.log', type=argparse.FileType('w'), help='Log file name')
    parser.add_argument('--coverage-out', metavar='FILE', default='coverage.tsv',
                        help='Output file for simulated genome coverage table', required=False)
    parser.add_argument('-z', '--num-samples', metavar='INT', type=int, default='1', required=True,
                        help='Number of transect samples')
    parser.add_argument('--dist', metavar='DISTNAME', choices=['equal', 'uniform', 'lognormal'],
                        help='Abundance profile distribution [equal, uniform, lognormal]')
    parser.add_argument('--lognorm-mu', metavar='FLOAT', type=float, default='1', required=False,
                        help='Log-normal relative abundance mu parameter')
    parser.add_argument('--lognorm-sigma', metavar='FLOAT', type=float, default='1', required=False,
                        help='Log-normal relative abundance sigma parameter')
    parser.add_argument('fasta', metavar='MULTIFASTA',
                        help='Input multi-fasta of all sequences')
    parser.add_argument('output_dir', metavar='DIR',
                        help='Output directory')
    args = parser.parse_args()

    @atexit.register
    def close_cov():
        coverage_file.close()

    seq_index = SeqIO.index(args.fasta, 'fasta')

    base_name = os.path.join(args.output_dir, args.output_name)

    r1_tmp = os.path.join(args.output_dir, '{0}1.fq'.format(TMP_OUTPUT))
    r2_tmp = os.path.join(args.output_dir, '{0}2.fq'.format(TMP_OUTPUT))
    seq_tmp = os.path.join(args.output_dir, TMP_INPUT)

    all_R1 = io_utils.open_output('{0}.r1.fq'.format(base_name), mode='w', compress=args.compress)
    all_R2 = io_utils.open_output('{0}.r2.fq'.format(base_name), mode='w', compress=args.compress)

    coverage_file = open(os.path.join(args.output_dir, args.coverage_out), 'w')

    RANDOM_STATE = numpy.random.RandomState(args.seed)
    child_seeds = RANDOM_STATE.randint(LOW_SEED_VALUE, HIGH_SEED_VALUE, args.num_samples).tolist()

    if args.profile:
        # if specified, read the static profile table from disk rather than calculate at runtime.
        # this will meant the same abundance profile is used in each sample -- in multisample mode.
        profile = abundance.read_profile(args.profile)

    # generate N simulated communities
    for n in xrange(0, args.num_samples):

        # generate abundance profile from global seeded random state -- if not using a static table
        if not args.profile:
            profile = abundance.generate_profile(RANDOM_STATE, seq_index, mode=args.dist,
                                                 lognorm_mu=args.lognorm_mu, lognorm_sigma=args.lognorm_sigma)

        for i, abn in enumerate(profile.values(), start=1):
            coverage_file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
                n+1, i, abn.chrom, abn.cell, abn.val * args.max_coverage))

        print 'Sample {0} Relative Abundances:'.format(n)
        profile.write_table(sys.stdout)

        r1_final = '{0}.{1}.r1.fq'.format(base_name, n+1)
        r2_final = '{0}.{1}.r2.fq'.format(base_name, n+1)

        r1_tmp = os.path.join(args.output_dir, '{0}1.fq'.format(TMP_OUTPUT, n+1))
        r2_tmp = os.path.join(args.output_dir, '{0}2.fq'.format(TMP_OUTPUT, n+1))

        output_R1 = io_utils.open_output(r1_final, mode='w', compress=args.compress)
        output_R2 = io_utils.open_output(r2_final, mode='w', compress=args.compress)

        try:

            # iteratively call ART for each taxon, accumulate the results
            for seq_id in profile:

                coverage = profile[seq_id].val * args.max_coverage
                print '\tRequesting {0:.4f} coverage for {1}'.format(coverage, seq_id)

                # iteration target for ART
                ref_seq = seq_index[seq_id]
                ref_len = len(ref_seq)
                SeqIO.write([ref_seq], seq_tmp, 'fasta')

                try:

                    subprocess.check_call([args.art_path,
                                           '-p',   # paired-end sequencing
                                           '-na',  # no alignment file
                                           '-rs', str(child_seeds[n]),
                                           '-m', str(args.insert_len),
                                           '-s', str(args.insert_sd),
                                           '-l', str(args.read_len),
                                           '-f', str(coverage),
                                           '-i', seq_tmp,
                                           '-o', os.path.join(args.output_dir, TMP_OUTPUT)],
                                          stdout=args.log, stderr=args.log)

                except OSError as e:
                    print "There was an error executing \"art_illumina\"."
                    print "Check that it is either on your PATH or specify it at runtime."
                    raise e
                except subprocess.CalledProcessError as e:
                    print e
                    raise e

                # count generated reads
                r1_n = 0
                for seq in SeqIO.parse(r1_tmp, 'fastq'):
                    r1_n += 1

                r2_n = 0
                for seq in SeqIO.parse(r2_tmp, 'fastq'):
                    r2_n += 1

                effective_cov = args.read_len * (r1_n + r2_n) / float(ref_len)
                print '\tGenerated {0} paired-end reads for {1}, {2:.3f} coverage'.format(r1_n, seq_id, effective_cov)
                if r1_n != r2_n:
                    print 'Error: paired-end counts do not match {0} vs {1}'.format(r1_n, r2_n)
                    sys.exit(1)

                io_utils.multicopy_tostream(r1_tmp, all_R1, output_R1)
                io_utils.multicopy_tostream(r1_tmp, all_R2, output_R2)

                os.remove(r1_tmp)
                os.remove(r2_tmp)
                os.remove(seq_tmp)

        except Exception as e:
            print e
            print 'Warning!! -- non-zero exit'
            sys.exit(1)
