#!/usr/bin/env python
import argparse
from collections import OrderedDict

import Bio.SeqIO as SeqIO
import numpy as np
import pysam
from Bio.Restriction import Restriction

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def upstream_dist(read, cut_sites, ref_length):
    """
    Upstream distance to nearest cut-site from the end of this
    reads alignment. Note: always assumes circular

    :param read: read in question
    :param cut_sites: the cut-sites for the aligning reference
    :param ref_length: the reference length
    :return:
    """
    if read.is_reverse:
        # reverse reads look left
        r_end = read.pos
        xi = np.searchsorted(cut_sites, r_end, side='left')
        if xi > 0:
            d = r_end - cut_sites[xi - 1]
        else:
            # crossing the origin
            d = ref_length - cut_sites[-1] + r_end
    else:
        # forward reads look right
        r_end = read.pos + read.alen
        xi = np.searchsorted(cut_sites, r_end, side='right')
        if xi < len(cut_sites):
            d = cut_sites[xi] - r_end
        else:
            # crossing the origin
            d = cut_sites[0] - (r_end - ref_length)
    return d


# Cigar codes that are not permitted.
NOT_ALLOWED = {2, 3, 6, 8}


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


class ContactMap:

    def __init__(self, bam, enz_name, seq_file, bin_size, ins_mean, ins_std, min_mapq, ref_min_len,
                 subsample=None, random_seed=None, strong=None):

        self.bam = bam
        self.bin_size = bin_size
        self.ins_mean = ins_mean
        self.ins_std = ins_std
        self.min_mapq = min_mapq
        self.subsample = subsample
        self.random_state = np.random.RandomState(random_seed)
        self.strong = strong

        print 'Counting reads in bam file...'
        self.total_reads = bam.count(until_eof=True)

        if ref_min_len:
            self.active_seq = []
            self.total_len = 0
            for n, li in enumerate(bam.lengths):
                if li > ref_min_len:
                    self.active_seq.append(bam.references[n])
                    self.total_len += li
            self.total_seq = len(self.active_seq)
        else:
            self.active_seq = bam.references
            self.total_len = sum(bam.lengths)
            self.total_seq = len(self.active_seq)

        self._init_sites(enz_name, seq_file)

        print 'Map based upon mapping containing:\n' \
              '\t{0} sequences\n' \
              '\t{1}bp total length\n' \
              '\t{2} mapped reads'.format(self.total_seq, self.total_len, self.total_reads)

        self.bin_count = int(self.total_len / self.bin_size) + 1
        print 'Map details:\n' \
              '\t{0}bp width\n' \
              '\t{1}x{1} dimension'.format(self.bin_size, self.bin_count)

        self.raw_map = None

    def _init_sites(self, enz_name, seq_file):

        self.cut_sites = OrderedDict()
        enzyme = get_enzyme_instance(enz_name)

        print 'Cut-site spacing naive expectation: {0}'.format(4**enzyme.size)
        offset = 0
        for rseq in SeqIO.parse(seq_file, 'fasta'):
            if rseq.id not in self.active_seq:
                continue
            xi = self.bam.references.index(rseq.id)
            locs = np.array(enzyme.search(rseq.seq, linear=False))
            locs.sort()  # pedantic check for order
            medspc = np.median(np.diff(locs))
            self.cut_sites[xi] = {'locs': locs, 'medspc': medspc, 'offset': offset}
            offset += len(rseq)
            print '{0} cut-sites for \"{1}\" has median spacing: {2:.1f}'.format(args.enzyme, rseq.id, medspc)

    def build_sorted(self):
        import tqdm

        def _simple_match(r):
            return not r.is_unmapped and r.mapq >= _mapq

        def _strong_match(r):
            return strong_match(r, self.strong, True, _mapq)

        # set-up match call
        _matcher = _strong_match if self.strong else _simple_match

        self.raw_map = self._init_map()

        with tqdm.tqdm(total=self.total_reads) as pbar:

            # maximum separation being 3 std from mean
            _wgs_max = self.ins_mean + 3*self.ins_std
            _hic_max = 1 * _wgs_max
            _mapq = self.min_mapq
            _map = self.raw_map
            _sites = self.cut_sites
            wgs_count = 0

            sub_thres = self.subsample
            if self.subsample:
                uniform = self.random_state.uniform

            bam.reset()
            bam_iter = bam.fetch(until_eof=True)
            while True:

                try:
                    r1 = bam_iter.next()
                    pbar.update()
                    while True:
                        # read records until we get a pair
                        r2 = bam_iter.next()
                        pbar.update()
                        if r1.query_name == r2.query_name:
                            break
                        r1 = r2
                except StopIteration:
                    break

                if sub_thres and sub_thres < uniform():
                    continue

                if not _matcher(r1) or not _matcher(r2):
                    continue

                if r1.is_read2:
                    r1, r2 = r2, r1

                assume_wgs = False
                if r1.is_proper_pair:
                    a, b = (r2, r1) if r1.is_reverse else (r1, r2)
                    ins_len = b.pos + b.alen - a.pos
                    if a.pos <= b.pos and ins_len < _wgs_max:
                        assume_wgs = True
                        wgs_count += 1
                        continue

                if not assume_wgs:
                    r1_dist = upstream_dist(r1, _sites[r1.reference_id]['locs'], r1.reference_length)
                    if r1_dist > _hic_max:
                        continue

                    r2_dist = upstream_dist(r2, _sites[r2.reference_id]['locs'], r2.reference_length)
                    if r2_dist > _hic_max:
                        continue

                r1pos = r1.pos if not r1.is_reverse else r1.pos + r1.alen
                r2pos = r2.pos if not r2.is_reverse else r2.pos + r2.alen

                ix1 = int((_sites[r1.reference_id]['offset'] + r1pos) / self.bin_size)
                ix2 = int((_sites[r2.reference_id]['offset'] + r2pos) / self.bin_size)

                if ix1 < ix2:
                    _map[ix1][ix2] += 1
                else:
                    _map[ix2][ix1] += 1

        print 'Assumed {0} wgs pairs'.format(wgs_count)
        print 'Total raw map weight {0}'.format(np.sum(_map))


    def _init_map(self, dt=np.int32):
        print 'Initialising contact map of {0}x{0} from total extent of {1}bp over {2} sequences'.format(
            self.bin_count, self.total_len, self.total_seq)
        return np.zeros((self.bin_count, self.bin_count), dtype=dt)


    def plot_map(self, pname, remove_diag=False, interp='none'):
        """
        Generate a plot (png) for the contact map.
        :param pname: the plot file name
        :param remove_diag: True - remove the central diagonal possibly improving the mapping of colour range
        to dynamic range of contact frequencies.
        """

        # being raw counts, many elements are zero
        # a copy is made to avoid side-effect of adding minimum value before taking logs
        _map = self.raw_map.astype(np.float64)

        # being pure counts, we just add 1 to every element. After log transformation
        # the empty elements will be zero.
        _map += 1.0

        if remove_diag:
            np.fill_diagonal(_map, np.min(_map))

        fig = plt.figure(frameon=False)
        img_h = float(self.bin_count) / 100
        fig.set_size_inches(img_h, img_h)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(np.log(_map), interpolation=interp)
        fig.savefig(pname, dpi=100)

    def write_map(self, oname):
        """
        Write a contact map as an ascii table to a file.
        :param oname: the file name to write
        """
        np.savetxt(oname, self.raw_map)


def get_enzyme_instance(enz_name):
    """
    Fetch an instance of a given restriction enzyme by its name.
    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    return getattr(Restriction, enz_name)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create a HiC contact map from mapped reads')
    parser.add_argument('--strong', default=None, type=float, help='Require a strong alignment match [None]')
    parser.add_argument('--sub-sample', default=None, type=float, help='Threshold probability for sub-sampling BAM')
    parser.add_argument('--ins-mean', default=500, type=int, help='Insert length mean [500]')
    parser.add_argument('--ins-std', default=50, type=int, help='Insert length stddev [50]')
    parser.add_argument('--mapq', default=0, type=int, help='Minimum mapping quality [0]')
    parser.add_argument('--min-len', default=None, type=int, help='Minimum subject sequence length [none]')
    parser.add_argument('--per-contig', default=False, action='store_true', help='Bins are per contig')
    parser.add_argument('--bin-size', type=int, default=25000, help='Bin size in bp (25000)')
    parser.add_argument('--remove-diag', default=False, action='store_true', help='Remove the central diagonal from plot')
    parser.add_argument('--enzyme', required=True, help='Restriction enzyme (case sensitive)')
    parser.add_argument('refseq', metavar='FASTA', help='Reference fasta sequence (in same order)')
    parser.add_argument('bamfile', metavar='BAMFILE', help='Name-sorted BAM file')
    parser.add_argument('output', metavar='OUTPUT_BASE', help='Output base name')
    args = parser.parse_args()

    with pysam.AlignmentFile(args.bamfile, 'rb') as bam:

        contacts = ContactMap(bam, args.enzyme, args.refseq, bin_size=args.bin_size, ins_mean=args.ins_mean,
                              ins_std=args.ins_std, ref_min_len=args.min_len, min_mapq=args.mapq,
                              subsample=args.sub_sample, strong=args.strong)

        contacts.build_sorted()

        print 'Writing raw output'
        contacts.write_map('{0}.raw.cm'.format(args.output))
        contacts.plot_map('{0}.raw.png'.format(args.output), remove_diag=args.remove_diag)
