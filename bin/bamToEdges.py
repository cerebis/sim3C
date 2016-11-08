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
import re
import sys
import traceback

import networkx as nx
import pysam

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


def good_match(cigartuples, min_match=None, match_start=True):
    """
    A confident mapped read.
    :param cigartuples: a read's CIGAR in tuple form
    :param min_match: the minimum number of matched base positions
    :param match_start: alignment must begin at read's first base position
    :return:
    """

    # restrict tuples to a subset of possible conditions
    if len(NOT_ALLOWED & set([t[0] for t in cigartuples])) != 0:
        return False

    # match the first N bases if requested
    elif match_start and cigartuples[0][0] != 0:
        return False

    # impose minimum number of matches
    elif min_match:
        n_matches = sum([t[1] for t in cigartuples if t[0] == 0])
        if n_matches < min_match:
            return False

    return True


def strong_match(mr, min_match=None, match_start=True, min_mapq=None):
    """
    Augment a good_match() by also checking for whether a read is secondary or
    supplementary.
    :param mr: mapped read to test
    :param min_match: the minimum number of matched base positions
    :param match_start: alignment must begin at read's first base position
    :param min_mapq: the minimum acceptable mapping quality. This can be mapper dependent. For BWA MEM
    any read which perfectly aligns in two distinct places, as mapq=0.
    """
    if min_mapq and mr.mapping_quality < min_mapq:
        return False

    elif mr.is_secondary or mr.is_supplementary:
        return False

    return good_match(mr.cigartuples, min_match, match_start)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create edge and node tables from a HiC bam file')

    parser.add_argument('--sim', default=False, action='store_true', help='Hi-C simulation read names')

    parser.add_argument('--minlen', type=int, required=False, default=0,
                        help='Minimum length in bp [0]')
    parser.add_argument('--mapq', default=0, type=int, help='Minimum mapping quality [0]')
    parser.add_argument('--strong', default=None, type=int,
                        help='Accept only mapped reads with no disagreements only clipping')

    parser.add_argument('--remove-selfloops', action='store_true', default=False,
                        help='Remove self-loops from nodes')
    parser.add_argument('--preserve-zerodeg', action='store_true', default=False,
                        help='Preserve nodes of zero degree by adding self-loops of weight 1.')
    parser.add_argument('--recover-alts', action='store_true', default=False,
                        help='Recover the alternate alignments from BAM')

    parser.add_argument('--ofmt', default='graphml', choices=['graphml', 'gml'],
                        help='Output graph format [graphml]')
    parser.add_argument('--merged', metavar='BAM', required=True, nargs='+',
                        help='Single merged HiC to contigs bam file')
    parser.add_argument('-o', dest='outbase', metavar='OUTPUT_BASE', required=True, help='Basis of output name')
    args = parser.parse_args()

    if args.remove_selfloops and args.preserve_zerodeg:
        print 'Error: remove-selfloops and preserve-zerodeg are mutually exclusive'
        sys.exit(1)

    # the resulting graph
    g = nx.Graph()

    # an intermediate dictionary recording all linkages between contigs (edges to be)
    linkage_map = {}

    if args.recover_alts:
        if not args.strong:
            print 'Recover alts can only be used with strong matching'
            sys.exit(1)
        print 'Recovery of alternate alignments requested. Please note that with BWA MEM\n' \
              'limits the number of alternates reported to <=5 by default.'

    bam_files = args.merged
    if len(set(bam_files)) != len(bam_files):
        print 'Error: some input bam files were listed twice'
        sys.exit(1)

    alt_rej = 0
    alt_ok = 0
    reject = 0
    accept = 0

    for bam_num, bam_fn in enumerate(bam_files, start=1):
        print 'Parsing {0}...'.format(bam_fn)

        # record the library source file. Note, collection attributes are not
        # supported by some formats, so this creates one per source.
        g.graph['lib{0}'.format(bam_num)] = bam_fn

        print 'Creating contig nodes...'
        with pysam.AlignmentFile(bam_fn, 'rb') as bam_hndl:
            g.add_nodes_from(zip(bam_hndl.references, [{'length': li} for li in bam_hndl.lengths]))

        print 'There were {0} referenced contigs in alignment file'.format(g.order())

        print 'Analyzing mapped reads...'
        # Read the sam file and build a linkage map
        with pysam.AlignmentFile(bam_fn, 'rb') as bam_hndl:
            iter_bam = bam_hndl.fetch()
            for mr in iter_bam:

                # used if alt recovery is requested
                contig_set = None

                if mr.reference_id == -1:
                    reject += 1
                    continue

                if not mr.cigarstring:
                    reject += 1
                    continue

                try:
                    # apply the constraints on alignment coverage and percent identity
                    if args.strong and not strong_match(mr, args.strong, True, args.mapq):
                        reject += 1
                        continue

                except Exception as e:
                    exc_info = sys.exc_info()
                    traceback.print_exception(*exc_info)
                    print 'alignment object: {0}'.format(mr)
                    sys.exit(1)

                accept += 1

                # TODO Scaling: also here, indexing the linkage map by read-ids is expensive with long string ids.
                if args.sim:
                    read = (mr.query_name[:-3], bam_num)
                    suffix = mr.query_name[-3:]
                    if suffix != 'fwd' and suffix != 'rev':
                        raise IOError('simulated read names must end in either "fwd" or "rev"')
                else:
                    read = (mr.query_name, bam_num)

                # TODO Scaling: potentially large memory cost to use reference string names rather than integer ids
                # TODO contig that this alignment line refers to directly
                contig = bam_hndl.getrname(mr.reference_id)

                # if requested, add all acceptable alternate alignments to linkage map
                if args.recover_alts:

                    # first, add the main mapped position -- to what will now be a set of positions
                    contig_set = {contig}

                    try:
                        # XA field contains alternate alignments for read, semi-colon delimited
                        alts_field = mr.get_tag('XA')
                        hit_list = alts_field.split(';')

                        for hit in hit_list:
                            hrec = hit.rstrip(';').split(',')
                            if len(hrec) == 1 and not hrec[0]:
                                continue

                            ctg, pos, cig, nm = hrec

                            # test confidence of alignment
                            if int(nm) > 0 or not good_match(cigar_to_tuple(cig), args.strong, True):
                                alt_rej += 1
                                continue

                            contig_set.add(ctg)
                            alt_ok += 1

                    except KeyError:
                        # if there was no XA tag, just continue
                        pass

                read_dir = True if not mr.is_reverse else False

                if contig_set:
                    # append all the mapped locations after recovery
                    # a list comprehension is faster here than a generator and the list
                    # _should_ not get that large to be a concern
                    linkage_map.setdefault(read, []).extend([(ctg_i, read_dir) for ctg_i in contig_set])
                else:
                    linkage_map.setdefault(read, []).append((contig, read_dir))

            print 'For {0} -- rejected {1} accepted {2}, rejection rate={3:.1f}%'.format(
                    bam_fn, reject, accept, float(reject) / (reject + accept) * 100.)

        print 'Overall: rejected {0} accepted {1}, rejection rate={2:.1f}%'.format(
                reject, accept, float(reject) / (reject + accept) * 100.)

        if args.recover_alts:
            print 'Recover alts: rejected {0} accepted {1}, rejection rate={2:.1f}%'.format(
                    alt_rej, alt_ok, float(alt_rej) / (alt_rej + alt_ok) * 100.)

    # From the set of all linkages, convert this information
    # into inter-contig edges, where the nodes are contigs.
    # Count the number of redundant links as a raw edge weight.
    print 'Creating graph from linkage map'
    
    print 'Input paired and unpaired inserts: {0}'.format(len(linkage_map))

    # TODO: Now we make a graph from all the links. This will likely double the memory
    # TODO footprint and we're in effect copying the information. Can we just build the
    # TODO graph directly? Why the linkage_map dictionary intermediate?
    # TODO Answer: easiest, we would need to read-pairs from BAMS. So either interleaved, expecting
    # TODO successive lines contain pairs, or both files in parallel. Logic needs to deal with
    # TODO iterations not representing pairs.
    paired = 0
    total = 0
    for insert, linkages in linkage_map.iteritems():
        n_links = len(linkages)
        for i in xrange(n_links):
            for j in xrange(i+1, n_links):

                u, udir = linkages[i]
                v, vdir = linkages[j]

                if udir == vdir:
                    # linkages only between opposing directions of ligation product
                    continue

                if g.has_edge(u, v):
                    g[u][v]['weight'] += 1
                else:
                    g.add_edge(u, v, weight=1, lib=insert[1])

                paired += 1

    print 'Paired count: {0}'.format(paired)
    print 'Initial graph stats: order {0} size {1}'.format(g.order(), g.size())

    # prune self-loops edges from ligation products if self-loops were not requested
    if args.remove_selfloops:
        print 'Removing self-loops from all nodes.'
        g.remove_edges_from(g.selfloop_edges())
        print 'Post self-loop pruning: order {0} size {1}'.format(g.order(), g.size())

    # if requested, add weight=1 self-loop to each node.
    # can be necessary when representing graphs in some formats for isolated nodes.
    if args.preserve_zerodeg:
        print 'Checking for isolated zero-degree nodes'
        n_sl = 0
        for v in g.nodes():
            if not g.has_edge(v, v):
                g.add_edge(v, v, weight=1)
                n_sl += 1
        print 'Self-loops were added to {0} isolated zero-degree nodes'.format(n_sl)
        print 'Post node preservation: order {0} size {1}'.format(g.order(), g.size())

    # filter nodes less than min length
    if args.minlen > 0:
        print 'Filtering nodes shorter than {0} bp'.format(args.minlen)
        nlist = g.nodes()
        n_rem = 0
        for n in nlist:
            if g.node[n]['length'] < args.minlen:
                g.remove_node(n)
                n_rem += 1
        print 'Length filtering removed {0} nodes'.format(n_rem)
        print 'Post length filtering: order {0} size {1}'.format(g.order(), g.size())

    print 'Writing edges'
    with open('{0}.edges.csv'.format(args.outbase), 'w') as h_out:
        h_out.write("SOURCE TARGET RAWWEIGHT TYPE\n")
        for u, v, dat in g.edges(data=True):
            h_out.write('{0} {1} {2} UNDIRECTED\n'.format(u, v, dat['weight']))

    print 'Writing nodes'
    with open('{0}.nodes.csv'.format(args.outbase), 'w') as h_out:
        h_out.write('ID LENGTH\n')
        for v, dat in g.nodes(data=True):
            h_out.write('{0} {1[length]}\n'.format(v, dat))

    print 'Writing graph'
    if args.ofmt == 'graphml':
        nx.write_graphml(g, '{0}.graphml'.format(args.outbase))
    elif args.ofmt == 'gml':
        nx.write_gml(g, '{0}.gml'.format(args.outbase))
