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

"""
Convert a SAM file to edge CSV file suitable for importing into Gephi

"""

import argparse
import math
import re
import sys
import traceback

import networkx as nx
import pysam


class Edge:
    """Represents an edge in the network of contigs linked
    by Hi-C read pairs.
    """

    def __init__(self, nodes=[]):
        nodes.sort()
        self.nodes = nodes
        self.weight = 1

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.nodes[0] == other.nodes[0] and self.nodes[1] == other.nodes[1]

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return '{source} {target} {weight} {normweight}'.format(
                source=self.nodes[0].id,
                target=self.nodes[1].id,
                weight=str(self.weight),
                normweight=str(self.norm_weight()))

    def inc_weight(self):
        self.weight += 1

    def get_id(self):
        return self.nodes[0].id + self.nodes[1].id

    def norm_weight(self):
        return self.weight / math.sqrt(self.nodes[0].length * self.nodes[1].length)


class Node:
    """Represents a node in the network
    """

    def __init__(self, id, length, reads):
        self.id = id
        self.length = length
        self.reads = reads

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.id == other.id

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id

    def __le__(self, other):
        return self.id <= other.id

    def __ge__(self, other):
        return self.id >= other.id

    def __str__(self):
        return '{id} {length} {reads}'.format(
                id=self.id, length=self.length, reads=self.reads)


CODE2CIGAR = dict([y, x] for x, y in enumerate("MIDNSHP=X"))
CIGAR_REGEX = re.compile(r"(\d+)([MIDNSHP=X])")
def cigar_to_tuple(cigar):
    return [(CODE2CIGAR[t[1]], int(t[0])) for t in CIGAR_REGEX.findall(cigar)]


def update_linkage_map(l):
    """Parse the line for new information about contig linkages. These
    may be self-self linkages or between inter-contig.
    """
    field = l.rstrip('\n').lstrip().split()
    read = field[0][:-3]
    rdir = field[0][-3:]
    contig = field[2]
    linkage = linkage_map.get(read)
    if linkage is None:
        linkage_map[read] = [(contig, rdir)]
    else:
        linkage.append((contig, rdir))


# Filter lines beginning with '@' and any line where the
# subject sequence is listed as '*'.
def filter(line):
    if line.startswith('@'):
        return True
    fields = line.rsplit()
    if fields[2] == '*':
        return True
    return False


def split_name(query_name):
    """
    Following our HiC naming convention, split query names into their
    components: (id, direction)

    :param query_name: query name to split
    :return: (id, direction)
    """
    return query_name[:-3], query_name[-3:]


cig_matcher = re.compile(r'([0-9]+)M')


def count_matches(cigar):
    """
    Sum the length of aligned regions listed in CIGAR pattern
    :param cigar: SAMtools CIGAR string
    :return: total number of aligned bases
    """
    return sum([int(seg_len) for seg_len in cig_matcher.findall(cigar)])


NOT_ALLOWED = {2, 3, 6, 8}


def good_match(cigartuples, min_match=None, match_start=True):

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
    Check that a mapped read does not contain disagreeing sequence. Only
    clipped regions are permitted.
    """
    if mr.is_secondary or mr.is_supplementary:
        return False

    if min_mapq and mr.mapping_quality < min_mapq:
        return False    

    return good_match(mr.cigartuples, min_match, match_start)

    # if len(NOT_ALLOWED & set([t[0] for t in mr.cigartuples])) == 0:
    #     if match_start and mr.cigartuples[0][0] != 0:
    #         return False
    #
    #     if mr.mapping_quality < min_mapq:
    #         return False
    #
    #     if min_match:
    #         n_matches = sum([t[1] for t in mr.cigartuples if t[0] == 0])
    #         if n_matches < min_match:
    #             return False
    #
    #     return True
    #
    # else:
    #     return False


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create edge and node tables from a HiC bam file')
    parser.add_argument('--recover-alts', action='store_true', default=False,
                        help='Recover the alternate alignments from BAM')
    parser.add_argument('--sim', default=False, action='store_true', help='Hi-C simulation read names')
    parser.add_argument('--afmt', choices=['bam', 'psl'], default='bam', help='Alignment file format (bam)')
    parser.add_argument('--minid', type=float, required=False, default=95.0,
                        help='Minimum percentage identity for alignment (95)')
    parser.add_argument('--minlen', type=int, required=False, default=0,
                        help='Minimum length in bp (1000)')
    parser.add_argument('--mincov', type=float, required=False, default=0.5,
                        help='Minimum coverage of query by alignment (0.5)')

    # parser.add_argument('--wgs', dest='wgs2ctg', metavar='WGS_BAM', nargs=1, help='WGS reads to contigs bam file')
    parser.add_argument('-s', '--add-selfloops', action='store_true', default=False,
                        help='Add self-loops to nodes')
    parser.add_argument('--mapq', default=60, type=int, help='Minimum mapping quality [60]')
    parser.add_argument('--strong', default=None, type=int,
                        help='Accept only mapped reads with no disagreements only clipping')
    parser.add_argument('--graphml', nargs=1, help='Write graphml file')
    parser.add_argument('--split', metavar='BAM', nargs=2, help='Split R1/R2 HiC to contigs bam files')
    parser.add_argument('--merged', metavar='BAM', nargs=1, help='Single merged HiC to contigs bam file')
    parser.add_argument('edge_csv', metavar='EDGE_CSV', nargs=1, help='Edges csv output file')
    parser.add_argument('node_csv', metavar='NODE_CSV', nargs=1, help='Nodes csv output file')
    args = parser.parse_args()

    # TODO We dont need to reference this file.
    # TODO If a scaffold does appear in HiC data then it is not correct to introduce it.

    g = nx.Graph()
    linkage_map = {}

    # input is a BAM file
    if args.afmt == 'bam':

        if args.merged and args.split:
            print 'Error: you must provide either merged OR split bam files'
            sys.exit(1)

        elif args.split and args.split[0] == args.split[1]:
            print 'Error: split files must differ'
            sys.exit(1)

        if not args.strong and args.recover_alts:
            print 'Recover alts can only be used with strong matching'
            sys.exit(1)

        bam_files = args.split if args.split else args.merged

        alt_rej = 0
        alt_ok = 0
        reject = 0
        accept = 0
        for rdir, fn in enumerate(bam_files, start=1):
            print 'Parsing {0}...'.format(fn)

            # rdir = 'fwd' if R == 1 else 'rev'

            print 'Creating contig nodes...'
            with pysam.AlignmentFile(fn, 'rb') as bf:
                g.add_nodes_from(zip(bf.references, [{'length': li} for li in bf.lengths]))
                #for i in xrange(len(bf.references)):
                #    rn = bf.references[i]
                #    if rn not in g:
                #        g.add_node(rn, length=bf.lengths[i])

            print 'There were {0} referenced contigs in alignment file'.format(g.order())

            print 'Analyzing mapped reads...'
            # Read the sam file and build a linkage map
            with pysam.AlignmentFile(fn, 'rb') as bf:
                iter_bam = bf.fetch()
                for mr in iter_bam:

                    contig_set = set()

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

                    #else:
                    #    perid = dict(mr.cigartuples)[0] / float(mr.query_length) * 100.0
                    #    cov = mr.query_alignment_length / float(mr.query_length)
                    #    if cov < args.mincov or perid < args.minid:
                    #        reject += 1
                    #        continue

                    accept += 1

                    if args.sim:
                        read = mr.query_name[:-3]
                        suffix = mr.query_name[-3:]
                        if suffix == 'fwd':
                            rdir = 0
                        elif suffix == 'rev':
                            rdir = 1
                        else:
                            raise IOError('simulated read names must end in "fwd" or "rev"')
                    else:
                        read = mr.query_name

                    # contig that this alignment line refers to directly
                    contig_set.add(bf.getrname(mr.reference_id))

                    # if requested, add also alternate alignment contigs for linkage map
                    if args.recover_alts:
                        try:
                            # XA field contains alternate alignments for read, semi-colon delimited
                            alts_field = mr.get_tag('XA')
                            hit_list = alts_field.split(';')

                            for hit in hit_list:
                                try:
                                    hrec = hit.rstrip(';').split(',')
                                    if len(hrec) == 1 and not hrec[0]:
                                        continue

                                    ctg, pos, cig, nm = hrec

                                    if int(nm) > 0 or not good_match(cigar_to_tuple(cig), args.strong, True):
                                        alt_rej += 1
                                        continue

                                    contig_set.add(ctg)
                                    alt_ok += 1

                                except Exception as e:
                                    print e
                                    raise IOError('alternate alignment did not contain four fields. [{0}] {1}'.format(
                                            alts_field, hit))

                        except KeyError:
                            pass

                    #ctg_assocs = [(ctg, rdir) for ctg in contig_set]
                    rdir = 1 if not mr.is_reverse else 2
                    ctg_assocs = [(ctg, rdir) for ctg in contig_set]

                    linkage = linkage_map.get(read)
                    if linkage is None:
                        linkage_map[read] = ctg_assocs
                        #linkage_map[read] = [(ctg, 1) for ctg in contig_set]
                    else:
                        linkage.extend(ctg_assocs)
                        #linkage.extend([(ctg, 2) for ctg in contig_set])

                print 'For {0} -- rejected {1} accepted {2}, rejection rate={3:.1f}%'.format(
                        fn, reject, accept, float(reject) / (reject + accept) * 100.)

        print 'Overall: rejected {0} accepted {1}, rejection rate={2:.1f}%'.format(
                reject, accept, float(reject) / (reject + accept) * 100.)

        if args.recover_alts:
            print 'Recover alts: rejected {0} accepted {1} rate={2:.1f}%'.format(
                    alt_rej, alt_ok, float(alt_rej) / (alt_rej + alt_ok))

    # input is a PSL file
    elif args.afmt == 'psl':

        if args.recover_alts:
            print 'Recovering alternate alignments is only applicable to BAM file parsing'
            sys.exit(1)

        # line marks a non-header line
        psl_dataline = re.compile(r'^[0-9]+\t')

        with open(args.hic2ctg[0], 'r') as h_in:
            for line in h_in:

                # skip header fields
                if not psl_dataline.match(line):
                    continue

                fields = line.rsplit()

                # qname = fields[9]
                # assumed naming convention of HiC simulated reads.
                read, rdir = split_name(fields[9])

                contig = fields[13]
                alen = int(fields[12]) - int(fields[11]) + 1
                qlen = int(fields[10])

                # Taken from BLAT perl script for calculating percentage identity
                matches = int(fields[0])
                mismatches = int(fields[1])
                repmatches = int(fields[2])
                q_num_insert = int(fields[4])
                perid = (1.0 - float(mismatches + q_num_insert) / float(matches + mismatches + repmatches)) * 100.0

                if not g.has_node(contig):
                    g.add_node(contig, length=int(fields[14]))

                # ignore alignment records which fall below mincov or minid
                # wrt the length of the alignment vs query sequence.
                if float(alen) / float(qlen) < args.mincov or perid < args.minid:
                    continue

                linkage = linkage_map.get(read)
                if linkage is None:
                    linkage_map[read] = [(contig, rdir)]
                else:
                    linkage.append((contig, rdir))

    # From the set of all linkages, convert this information
    # into inter-contig edges, where the nodes are contigs.
    # Count the number of redundant links as a raw edge weight.
    print 'Creating graph from linkage map'
    
    print 'Input paired and unpaired inserts: {0}'.format(len(linkage_map))

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
                    g.add_edge(u, v, weight=1)

                paired += 1

    # prune self-loops edges from ligation products
    g.remove_edges_from(g.selfloop_edges())

    print 'paired={0} order={1} size={2}'.format(paired, g.order(), g.size())

    # if requested, add weight=1 self-loop to each node.
    # can be necessary when representing graphs in some formats for isolated nodes.
    if args.add_selfloops:
        for v in g.nodes():
            if not g.has_edge(v, v):
                g.add_edge(v, v, weight=1)

    # filter nodes on length
    nlist = g.nodes()
    for n in nlist:
        if g.node[n]['length'] < args.minlen:
            g.remove_node(n)

    if args.graphml is not None:
        nx.write_graphml(g, args.graphml[0])

    with open(args.edge_csv[0], 'w') as h_out:
        h_out.write("SOURCE TARGET RAWWEIGHT TYPE\n")
        for u, v, dat in g.edges(data=True):
            h_out.write('{0} {1} {2} UNDIRECTED\n'.format(u, v, dat['weight']))

    with open(args.node_csv[0], 'w') as h_out:
        h_out.write('ID LENGTH\n')
        for v, dat in g.nodes(data=True):
            h_out.write('{0} {1[length]}\n'.format(v, dat))


