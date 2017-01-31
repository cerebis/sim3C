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

import re

import pysam
import tqdm
from Bio import SeqIO
from Bio.Restriction import Restriction
from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict
from intervaltree import IntervalTree
from itertools import product


def get_enzyme_instance_ipython(enz_name):
    """
    An alternative method to fetch an instance of a given restriction enzyme by its
    name using a work-around which avoids exception with getattr() in iPython.

    Ackn: StackOverflow: user xbello.
    See: http://stackoverflow.com/questions/20381912/type-object-restrictiontype-has-no-attribute-size

    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    r_type_names = [rt for tid, (rt, rn) in typedict.iteritems() if enz_name in rn][0]
    r_clz = tuple(getattr(Restriction, rt) for rt in r_type_names)
    return Restriction.AbstractCut(enz_name, r_clz, rest_dict[enz_name])


def get_enzyme_instance(enz_name):
    """
    Fetch an instance of a given restriction enzyme by its name.
    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    return getattr(Restriction, enz_name)


def site_intervaltree(seq, enzyme):
    """
    Initialise an intervaltree representation of an enzyme's cutsites across a given sequence.
    Whether a position involves a cutsite can be queried with tree[x] or tree[x1:x2].
    :param seq: the sequence to digest
    :param enzyme: the restriction enzyme used in digestion
    :return: an intervaltree of cutsites.
    """
    tr = IntervalTree()
    size = enzyme.size
    offset = enzyme.fst3
    for si in enzyme.search(seq):
        start = si + offset - 1
        tr.addi(start, start + size)
    return tr


# Mapping of cigar characters to code values
CODE2CIGAR = dict((y, x) for x, y in enumerate("MIDNSHP=X"))

# Pattern that finds each unit of a cigar i.e. 10M or 9H
CIGAR_ANY = re.compile(r"(\d+)([MIDNSHP=X])")

# Pattern that find only matches (M)
CIGAR_MATCHES = re.compile(r'([0-9]+)M')

# Cigar codes that are not permitted.
NOT_ALLOWED = {2, 3, 6, 8}


def cigar_to_tuple(cigar):
    """
    Convert a CIGAR string into code values
    :param cigar:
    :return:
    """
    return [(CODE2CIGAR[t[1]], int(t[0])) for t in CIGAR_ANY.findall(cigar)]


def count_matches(cigar):
    """
    Sum the length of aligned regions listed in CIGAR pattern
    :param cigar: SAMtools CIGAR string
    :return: total number of aligned bases
    """
    return sum([int(seg_len) for seg_len in CIGAR_MATCHES.findall(cigar)])


def next_pair(it):
    """
    Return a pair (R1/R2) of reacs from a name sorted bam file. Raises StopIteration
    when the end of file reached. Order of reads is not guaranteed. i.e. R1,R2 or R2,R1.
    The method will continue to scan a file until it finds a pair or eof.
    :param it: a Pysam sam/bam file iterator
    :return: tuple R1/R2 or R2/R1
    """
    r1 = it.next()
    while True:
        # read records until we get a pair
        r2 = it.next()
        if r1.query_name == r2.query_name:
            return r1, r2
        r1 = r2


def aln_coverage(aln_list):
    """
    Calculate the coverage across the reported alignments for a given read. This will most
    often involve only a single alignment, but also considers non-overlapping alignments
    reported by BWA MEM scavenged from the XP tag. Reports the number of bases covered
    (<=read_len) and the overlap between them (normally 0).
    :param aln_list: the list of alignments for a read
    :return: dict {coverage: xx, overlap: yy}
    """
    # using an intervaltree for this
    tr = IntervalTree()
    tot = 0
    for ti in aln_list:
        if ti['is_reverse']:
            # reversed reads must be tallied from the opposite end
            n = ti['total']
            for op, nb in ti['cigartuple']:
                if op == 0:
                    tr.addi(n - nb, n)
                    tot += nb
                n -= nb
        else:
            # forward mapped reads tally from start position
            n = 0
            for op, nb in ti['cigartuple']:
                if op == 0:
                    tr.addi(n, n + nb)
                    tot += nb
                n += nb
    # lazy means of merging intervals
    tr.merge_overlaps()
    cov = sum([i.end - i.begin for i in tr])
    return {'coverage': cov, 'overlap': tot - cov, 'has_multi': len(aln_list) > 1}


def infer_inslen(r1, r2, ref_len):
    """
    Infer the length of the insert for R1,R2.
    :param r1: read 1
    :param r2: read 2
    :param ref_len: the total length of the reference sequence.
    :return: insert length in bp.
    """
    # switch reads around so that r1 is always the forward strand
    if r1.is_reverse:
        r2, r1 = r1, r2
    sc = 0
    if r1.pos == 0:
        if r1.cigar[0][0] == 4:
            sc += r1.cigar[0][1]
    il = sc + r2.pos + r2.alen - r1.pos + 1
    if il < 0:
        il += ref_len
    return il


def append_xp_alignments(aln_list, xp_record):
    """
    Add any additional alignments reported in BWA MEM's XP record to the list.
    This string record is doubly delimited, with multiple alignments delimited by
    semi-colons and the individual fields of each alignment delimited by commas.
    :param aln_list: the list of alignments
    :param xp_record: the xp_data for this read
    :return: a list of dicts containing each alignment's details, with any additional alignments appended.
    """
    # split the records
    for aln_i in xp_record.split(';'):
        # skip empty records
        if not aln_i:
            continue

        # split the fields for aln i
        ti = aln_i.split(',')
        pos = int(ti[1])
        if pos < 0:
            # position also encodes direction with sign.
            pos = -pos
            # convert to 0-based
            pos -= 1
            is_rev = True
        else:
            # convert to 0-based
            pos -= 1
            is_rev = False

        # parse the cigar, compute some extra details
        cigtup = cigar_to_tuple(ti[2])
        alen = sum([num for op, num in cigtup if op == 0])
        tot = sum([num for op, num in cigtup])

        aln_list.append({'ref': ti[0],
                         'pos': pos,
                         'is_reverse': is_rev,
                         'cigar': ti[2],
                         'cigartuple': cigtup,
                         'nm': int(ti[3]),
                         'mapq': int(ti[4]),
                         'alen': alen,
                         'total': tot})


def parse_all_alignments(read):
    """
    Parse a read record for all alignments. Both the primary and the "other non-overlapping"
    alignments reported by BWA MEM (since v0.7.3) in the XP tag. The extra alignments are
    important in situations where reads are split -- as is the case for HiC/Meta3C experiments
    when reads cross the ligation junction.
    :param read: the read record to parse
    :return: a list of dicts containing alignment information.
    """
    # primary
    aln_list = [{'ref': read.reference_name,
                 'pos': read.pos,
                 'alen': read.alen,
                 'total': sum([num for op, num in read.cigartuples]),
                 'is_reverse': read.is_reverse,
                 'cigar': read.cigarstring,
                 'cigartuple': read.cigartuples,
                 'nm': None,
                 'mapq': read.mapq}]

    # alternative non-overlapping
    if read.has_tag('XP'):
        append_xp_alignments(aln_list, read.get_tag('XP'))

    return aln_list


def alignment_cutsite_status(aln, ref_sites, gap=1):
    """
    Categorise the cutsites contained in an alignment as either 3 or 5 terminated or internal.
    :param aln: the alignment to inspect
    :param ref_sites: the intervaltree of cutsites for the reference
    :param gap: a gap in bp, permitting adjust of "close enough" situations at the termini.
    :return: a dict of booleans {'5p': t/f, '3p': t/f, 'internal': t/f}
    """
    x1 = aln['pos']
    x2 = x1 + aln['alen'] - 1
    is_left_term = len(ref_sites[x1:x1 + gap + 1]) > 0
    is_right_term = len(ref_sites[x2 - gap:x2 + 1]) > 0

    # internal sites are only reported of wholly contained with the specified range.
    has_internal = len(ref_sites.search(x1 + gap, x2 - gap, True)) > 0

    if aln['is_reverse']:
        return {'5p': is_right_term, '3p': is_left_term, 'internal': has_internal}
    else:
        return {'5p': is_left_term, '3p': is_right_term, 'internal': has_internal}


def read_cutsite_status(aln_terms):
    """
    Determine a read-wide assessment of cutsite termination. This involves inspecting
    all alignment termination details to produce a single boolean of whether the read
    was 3p, 5p terminated or internally contains a cutsite.
    :param aln_terms: list of terminations
    :return:
    """
    return {'5p': any(trm_i['5p'] for trm_i in aln_terms),
            '3p': any(trm_i['3p'] for trm_i in aln_terms),
            'internal': any(trm_i['internal'] for trm_i in aln_terms)}


def parse_bam(bam, ref_seq, enzyme):

    print 'Counting reads...',
    total_reads = bam.count(until_eof=True)
    print ' found {0} reads'.format(total_reads)

    print 'Creating cutsite interval tree...',
    site_tree = site_intervaltree(ref_seq, enzyme)
    print ' found {0} sites'.format(len(site_tree))

    counts = {
        'all_pairs': 0,
        'incomplete': 0,
        '3p_term': 0,
        '3p_trunc': 0,
        'readthru_conf': 0,
        'readthru_multi': 0,
        'readthru': 0,
        'nosite': 0,
        'hassite': 0,
        'proper': 0, }

    ins_len = []

    outh = open('tab.csv', 'w')

    print 'Beginning parsing...'
    with tqdm.tqdm(total=total_reads) as pbar:

        bam.reset()
        bam_iter = bam.fetch(until_eof=True)
        while True:

            try:
                pair = next_pair(bam_iter)
                pbar.update(2)
                counts['all_pairs'] += 1

            except StopIteration:
                break

            r1, r2 = pair

            # get the inferred full length of each read for later
            r1_len = r1.infer_query_length()
            r2_len = r2.infer_query_length()

            if r1.is_unmapped or r2.is_unmapped:
                # ignore incompletely mapped pairs
                counts['incomplete'] += 1
                continue

            # get alignments for R1 and see if any involve a cutsite
            r1_algns = parse_all_alignments(r1)
            r1_count = len([inv for aln in r1_algns for inv in site_tree[aln['pos']:aln['pos'] + aln['alen']]])

            # get alignments for R2 and see if any involve a cutsite
            r2_algns = parse_all_alignments(r2)
            r2_count = len([inv for aln in r2_algns for inv in site_tree[aln['pos']:aln['pos'] + aln['alen']]])

            # if either read involves a cutsite, look a bit deeper
            if r1_count > 0 or r2_count > 0:

                counts['hassite'] += 1

                r1_aln_info = aln_coverage(r1_algns)

                r1_aln_status = [alignment_cutsite_status(aln_i, site_tree) for aln_i in r1_algns]
                r1_status = read_cutsite_status(r1_aln_status)

                r2_aln_info = aln_coverage(r2_algns)

                r2_aln_status = [alignment_cutsite_status(aln_i, site_tree) for aln_i in r2_algns]
                r2_status = read_cutsite_status(r2_aln_status)

                # was there a 3p termination of r1 or r2
                if r1_status['3p'] or r2_status['3p']:

                    counts['3p_term'] += 1

                    # was a read also incompletely aligned
                    if r1.alen < r1_len or r2.alen < r2_len:

                        counts['3p_trunc'] += 1

                # was there both 5p and 3p termination, multiple and short alignments
                if (r1_status['5p'] and r1_status['3p']) or (r2_status['5p'] and r2_status['3p']):

                    counts['readthru'] += 1

                    if (r1_aln_info['has_multi'] and r1.alen < r1_len) or \
                            (r2_aln_info['has_multi'] and r2.alen < r2_len):

                        counts['readthru_multi'] += 1

                        # do R1 and R2 possess an alignment in the same direction
                        r1_rev = [aln['is_reverse'] for aln in r1_algns]
                        r2_rev = [aln['is_reverse'] for aln in r2_algns]
                        if any(map(all, [x for x in product(r1_rev, r2_rev)])):

                            counts['readthru_conf'] += 1

            else:
                counts['nosite'] += 1

            if r1.is_proper_pair:
                counts['proper'] += 1
                ins_len.append(infer_inslen(r1, r2, bam.lengths[0]))

    print counts

    return ins_len

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Parse name-sorted BAM for ligation junctions')
    parser.add_argument('--plot', default=False, action='store_true', help='Plot insert lengths')
    parser.add_argument('enzyme', help='Enzyme name')
    parser.add_argument('seq', help='Reference fasta')
    parser.add_argument('bam', help='BAM file sorted by name')
    args = parser.parse_args()

    try:
        ref_seq = SeqIO.read(args.seq, 'fasta').seq
    except ValueError as er:
        print er.message

    with pysam.AlignmentFile(args.bam, 'rb') as bam:
        enzyme = get_enzyme_instance(args.enzyme)
        lengths = parse_bam(bam, ref_seq, enzyme)

    if args.plot:
        import matplotlib.pyplot as plt
        import numpy as np
        plt.hist(lengths, bins=(np.arange(1, 41) * 25).tolist() + [np.inf])
        plt.xlim(0, 1000)
        plt.show()
