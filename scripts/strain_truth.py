#!/usr/bin/env python
"""
meta-sweeper - for performing parametric sweeps of simulated
metagenomic sequencing experiments.
Copyright (C) 2016 "Aaron E Darling"

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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Record a ground truth SNV matrix relative to an ancestor')
    parser.add_argument('--mauve-path', default='mauve', required=True, help='Path to Mauve directory [default: mauve]')
    parser.add_argument('-o', '--output', metavar='TSV', required=True, help='Output name')
    parser.add_argument('fasta', metavar='MULTIFASTA', help='Input multi-fasta of evolved sequences')
    parser.add_argument('ref', metavar='MULTIFASTA', help='Input fasta of the reference sequence')
    args = parser.parse_args()


    seq_index = SeqIO.index(args.fasta, 'fasta')
    num_seqs = len(seq_index)
    snv_alleles = dict()
    ref_alleles = dict()
    for seq_id in seq_index:
        # split up the input 
        fname = seq_id + ".fa"
        fh = open(fname, 'w')
        seq = seq_index[seq_id]
        SeqIO.write([seq], fh, 'fasta')
        fh.close()
        
        # align this seq with progressiveMauve and export SNVs
        try:
            subprocess.call([args.mauve_path + "/linux-x64/progressiveMauve",
                         '--output=' + seq_id + ".xmfa",   # xmfa output file name
                         args.ref,  # reference genome file name
                         fname])
            subprocess.call(["java", '-cp',
                         args.mauve_path + "/Mauve.jar",
                         'org.gel.mauve.analysis.SnpExporter',
                         '-f', seq_id + ".xmfa",
                         '-o', seq_id + ".snvs"])
        except OSError as ex:
            print "There was an error starting the art_illumina subprocess."
            print "You may need to add its location to your PATH or specify it at runtime."
            raise ex
    
        # collate SNVs to a single table
        snvfile = open(seq_id + ".snvs")
        header = snvfile.readline()
        for line in snvfile:
            d = line.split("\t")
            if len(d) < 4:
                break
            if not d[3] in snv_alleles:
                snv_alleles[d[3]] = dict()
                ref_alleles[d[3]] = dict()
            snv_alleles[d[3]][seq_id] = d[0][1];
            ref_alleles[d[3]] = d[0][0];

    # create ground truth genotype table
    truth_file = open(args.output, 'w')
    header = "#ref_position"
    for seq_id in seq_index:
        header += "\t" + seq_id
    truth_file.write(header + "\n")
    for pos in snv_alleles:
        cur_line = str(pos)
        for seq_id in seq_index:
            cur_line += "\t"
            if seq_id in snv_alleles[pos]:
                cur_line += snv_alleles[pos][seq_id]
            else:
                cur_line += ref_alleles[pos]

        truth_file.write(cur_line + "\n")
        
    truth_file.close()

