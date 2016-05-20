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
from collections import OrderedDict

import Psl
import truthtable as tt

import numpy as np
import argparse
import sys
import re


# CIGAR matcher for aligning regions
matcher = re.compile(r'([0-9]+)M')


def count_aligned(cigar):
    """
    Sum the length of aligned regions listed in CIGAR pattern
    :param cigar: SAMtools CIGAR string
    :return: total number of aligned bases
    """
    return sum([int(seg_len) for seg_len in matcher.findall(cigar)])


class Alignment:
    """
    Object representing an alignment taken from some alignment format file.

    This object's uses value identity, governed by the values of query and reference names
    NOT memory reference.

    Not all files have direct access to all parameters. In particular, unless
    specified the alignment will not have percentage identity. This will be
    avoided if not defined.
    """
    def __init__(self, query_name, ref_name, bases, query_length, perid=None):
        self.query_name = query_name
        self.ref_name = ref_name
        if not isinstance(bases, (int, long)):
            raise TypeError('bases parameter must be an integer')
        self.align_length = bases
        self.query_length = query_length
        self.perid = perid

    def __repr__(self):
        return self.query_name + ' ' + self.ref_name

    def __str__(self):
        if self.perid:
            return '{0} {1} {2:.4} {3} {4} {5:.4}'.format(
                self.query_name, self.ref_name, self.coverage,
                self.query_length, self.align_length, self.perid)
        else:
            return '{0} {1} {2:.4} {3} {4}'.format(
                self.query_name, self.ref_name, self.coverage,
                self.query_length, self.align_length)


    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.query_name == other.query_name and self.ref_name == other.ref_name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(repr(self))

    def add_bases(self, bases):
        self.align_length += bases

    @property
    def coverage(self):
        return float(self.align_length) / float(self.query_length)


#def ordereddict_rep(dumper, data):
#    return dumper.represent_mapping(u'tag:yaml.org,2002:map', data.items(), flow_style=False)


def write_output(simple, fmt, minlen, mincov, minid, alignment_list, file_name):
    """
    Write result to file in a particular format.

    :param simple: use simple hard tt
    :param fmt: format type
    :param minlen: minimum query length
    :param mincov: minimum alignment coverage
    :param minid: minimum percentage identity
    :param alignment_list: alignment list
    :param file_name: output file name
    """
    if fmt == 'flat':
        write_alignments(minlen, mincov, minid, alignment_list, file_name)
    elif fmt == 'yaml':
        write_truth(simple, minlen, mincov, alignment_list, file_name)


def write_truth(simple, minlen, mincov, alignment_list, file_name):
    """
    Write YAML format truth table, where all relations between queries and subjects
    are included on a single line.

    Eg.

    seq1: [ctgA, ctgB, ctgC]
    seq2: [ctg111, ctg1]
    seq3: [ctg99]

    :param minlen: minimum query length
    :param mincov: minimum alignment coverage
    :param minid: minimum percentage identity
    :param alignment_list: alignment list
    :param file_name: output file name
    """
    truth = tt.TruthTable()
    for aln in alignment_list:
        if aln.align_length > minlen and aln.coverage > mincov:
            t = truth.get(aln.query_name)
            if t is None:
                t = {}
                truth.put(aln.query_name, t)
            t[aln.ref_name] = aln.align_length

    if simple:
        truth.write_hard(file_name)
    else:
        truth.write(file_name)


def write_alignments(minlen, mincov, minid, alignment_list, file_name):
    """
    Write a simple list of assignments, one relation per line.

    Eg.
    seq1 ctgA
    seq1 ctgB
    seq2 ctg111

    :param minlen: minimum query length
    :param mincov: minimum alignment coverage
    :param minid: minimum percentage identity
    :param alignment_list: alignment list
    :param file_name: output file name
    """
    with open(args.output_file[0], 'w') as h_out:
        for aln in alignment_list:
            if aln.align_length > minlen and aln.coverage > mincov:
                h_out.write(str(aln) + '\n')
            else:
                h_out.write('{0} null  \n'.format(aln.query_name))


def parse_sam(concat_sam_file):
    """
    Parse the SAM file for alignment lengths, where each
    query can accumulate total length there are multiple
    record against the same reference.

    :param concat_sam_file:
    :return: dictionary of Alignments objects
    """
    align_repo = OrderedDict()
    with open(concat_sam_file, 'r') as h_in:
        num_sec = 0
        for l in h_in:
            # skip header fields
            if l.startswith('@'):
                continue
            fields = l.rsplit()
            # test if secondary bit set (0x256)
            if (int(fields[1]) & (1 << 8)) != 0:
                num_sec += 1
                continue
            # count aligned bases from cigar matches
            aln_len = count_aligned(fields[5])
            aln = Alignment(fields[0], fields[2], aln_len, seq_length[fields[0]])
            if aln in align_repo:
                align_repo[aln].add_bases(aln_len)
            else:
                align_repo[aln] = aln

        print 'Found {0} secondary alignments, which were skipped'.format(num_sec)

    return align_repo


def parse_last(last_file):
    """
    Parse a LASTAL tabular alignment file for lengths, where
    each alignment of the same query/subject pair accumulate length
    in bases.

    :param last_file: last tabular format alignment file
    :return: dictionary of Alignment objects
    """
    align_repo = OrderedDict()
    with open(last_file, 'r') as h_in:
        for l in h_in:
            # skip header fields
            if l.startswith('#'):
                continue

            fields = l.rsplit()

            qname = fields[6]
            rname = fields[1]
            alen = int(fields[8])
            qlen = int(fields[10])

            # ignore alignment records which fall below mincov with
            # respect to the length of the alignment vs query sequence.
            if float(alen)/float(qlen) < args.mincov:
                continue

            aln = Alignment(qname, rname, alen, qlen)
            if aln in align_repo:
                align_repo[aln].add_bases(alen)
            else:
                align_repo[aln] = aln
    return align_repo


def contiguous_mask(mask):
    """
    Calculate the contiguous regions of the mask for where coverage exists, report their start and stop points.
    :param mask: the coverage mask.
    :return: list of stop/start locations, end inclusive
    """
    dmask = np.diff(mask)
    idx, = dmask.nonzero()
    idx += 1

    if mask[0]:
        idx = np.r_[0, idx]
    if mask[-1]:
        idx = np.r_[idx, mask.size]

    idx.shape = (-1, 2)

    # ends are inclusive.
    idx[:, 1] -= 1

    return idx


def interval_overlap(a, b):
    """
    Test if two intervals overlap.
    :param a: interval a
    :param b0: start of interval b
    :param b1: end of interval b
    :return: True if there is overlap between a and b
    """
    return a[0] <= b[1] and b[0] <= a[1]


def interval_overlap_with_gap(a, b, max_ol):
    """
    Test if two intervals overlap, allowing for a maximum overlapping of max_ol.
    :param a: interval a
    :param b: interval b
    :param max_ol: maximum allowed overlap
    :return: True if there is a sufficient overlap
    """
    if not interval_overlap(a, b):
        return False
    if a[0] < b[0] and a[1] - b[0] < max_ol:
        return False
    if b[0] < a[0] and b[1] - a[0] < max_ol:
        return False
    return True


def parse_psl2(psl_file):
    """
    Calculate a truth table by considering the accumulated coverage of query sequences onto
    the references. The coverage is treated as a binary mask and the total extent a contig
    covers each reference determines the degree of ownership.

    Writes out a full truth table to user specified output path.

    Note: ideally the PSL file should be sorted for descending alignment score.

    :param psl_file:
    :return: None -- this method presently breaks the logical flow of this script.
    """
    with open(psl_file, 'r') as h_in:

        all_hits = 0
        rejected = 0

        aln_masks = {}

        for aln in Psl.parse(h_in):

            all_hits += 1

            if aln.percent_id < 90:
                rejected += 1
                continue

            if aln.q_name not in aln_masks:
                aln_masks[aln.q_name] = {}

            if aln.t_name not in aln_masks[aln.q_name]:
                aln_masks[aln.q_name][aln.t_name] = np.zeros(int(aln.q_size))

            per_id = 0.01 * aln.percent_id
            mask_slice = aln_masks[aln.q_name][aln.t_name][aln.q_start:aln.q_end+1] #= aln.percent_id/100.
            mask_slice[np.where(mask_slice < per_id)] = per_id

        truth = {}
        weights = {}
        for n, ti in enumerate(aln_masks):
            masks = np.vstack(aln_masks[ti].values())
            names = np.array(aln_masks[ti].keys())
            covers = np.mean(masks, 1)
            idx = np.where(covers > 0.96)
            if idx[0].shape[0] > 0:
                truth[ti] = {}
                for i in np.nditer(idx):
                    truth[ti][str(names[i])] = float(covers[i])
            weights[ti] = masks.shape[1]
        t = tt.TruthTable()
        t.update(truth, weights)

        t.write(args.output_file[0])


def parse_psl(psl_file):
    """
    Parse a PSL converted from MAF

    :param psl_file: PSL format alignment file
    :return: dictionary of Alignment objects
    """
    all_hits = 0
    rejected = 0
    alignment_repo = OrderedDict()

    with open(psl_file, 'r') as h_in:

        for aln in Psl.parse(h_in):

            all_hits += 1

            # ignore alignment records which fall below mincov or minid
            # wrt the length of the alignment vs query sequence.
            if aln.coverage < args.mincov or aln.percent_id < args.minid:
                rejected += 1
                continue

            ai = Alignment(aln.q_name, aln.t_name, aln.length, aln.q_size, aln.percent_id)
            if ai in alignment_repo:
                alignment_repo[ai].add_bases(aln.length)
            else:
                alignment_repo[ai] = ai

        print 'Rejected {0}/{1} alignments due to constraints on ID {2} and Coverage {3}'.format(
            rejected, all_hits, args.minid, args.mincov)

    return alignment_repo

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Determine the location of query sequences on the reference')
    parser.add_argument('--simple', default=False, action='store_true',
                        help='Record simple truth table')
    parser.add_argument('--max-overlap', type=int, required=False, default=3,
                        help='Maximum bp overlap between alignment on the same query to the same subject [3]')
    parser.add_argument('--minid', type=float, required=False, default=95.0,
                        help='Minimum percentage identity for alignment')
    parser.add_argument('--minlen', type=int, required=False, default=1000,
                        help='Minimum length in bp')
    parser.add_argument('--mincov', type=float, required=False, default=0.5,
                        help='Minimum coverage of query by alignment')
    parser.add_argument('--afmt', choices=['bwa', 'last', 'psl', 'test'], default='last', help='Alignment file format')
    parser.add_argument('--qf', metavar='FASTA', help='Query fasta sequences')
    parser.add_argument('--ofmt', choices=['flat', 'yaml'], default='yaml', help='Output format')
    parser.add_argument('alignment_file', metavar='ALIGNMENT_FILE', nargs=1,
                        help='concat SAM, PSL or LAST table with query sequences aligned to reference')
    parser.add_argument('output_file', metavar='OUTPUT', nargs=1, help='Output file name')
    args = parser.parse_args()


    align_repo = None
    seq_length = None
    # minimum check that the necessary information has been supplied.
    # if so, fetch lengths of all query sequences
    if args.afmt == 'bwa':
        if args.qf is None:
            raise Exception('BWA alignment format also requires query fasta to be supplied')
        # Calculate the length of query sequences
        seq_length = {rec.id: len(rec) for rec in SeqIO.parse(args.qf, 'fasta')}
        align_repo = parse_sam(args.alignment_file[0])

    elif args.afmt == 'last':
        align_repo = parse_last(args.alignment_file[0])

    elif args.afmt == 'psl':
        align_repo = parse_psl(args.alignment_file[0])

    elif args.afmt == 'test':
        parse_psl2(args.alignment_file[0])
        sys.exit(0)

    if args.ofmt == 'flat':
        print 'Soft results always enabled for flat output format'

    print 'Accepted {0} alignments'.format(len(align_repo))

    #
    # We need to decide on assignment.
    #
    # There are query sequences which are assigned to more than one reference.
    # What shall be done with these? Simply taking the largest alignment
    # as the winner might be acceptable for large differences, but I've
    # seen examples where both are significant alignments and picking one
    # is quite arbitrary.
    #

    if not args.simple:
        '''
        Write out all assignments of queries to references.
        '''
        print 'Writing {0} soft (overlapping) assignments of queries to references'.format(len(align_repo.values()))
        write_output(args.simple, args.ofmt, args.minlen, args.mincov,
                     args.minid, align_repo.values(), args.output_file[0])

    else:
        '''
        Pick a winner

        For each scaffold, determine the best alignment by length. The alignment subject
        will then be considered the true source.
        '''
        winners = OrderedDict()
        for aln in align_repo.values():
            if aln.query_name not in winners or winners[aln.query_name].align_length < aln.align_length:
                winners[aln.query_name] = aln

        print 'Reduced to {0} winning hard assignments'.format(len(winners))
        write_output(args.simple, args.ofmt, args.minlen, args.mincov,
                     args.minid, winners.values(), args.output_file[0])

