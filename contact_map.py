#!/usr/bin/env python
import argparse
from collections import OrderedDict

import Bio.SeqIO as SeqIO
import numpy as np
import pysam
from Bio.Restriction import Restriction
from intervaltree import IntervalTree, Interval

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

    # reverse the tuple if this read is reversed, as we're testing from 0-base on read.
    cigtup = mr.cigartuples[::-1] if mr.is_reverse else mr.cigartuples

    return good_match(cigtup, min_match, match_start)


class ContactMap:

    def __init__(self, bam, enz_name, seq_file, bin_size, ins_mean, ins_sd, min_mapq, ref_min_len=0,
                 subsample=None, random_seed=None, strong=None, max_site_dist=None, spacing_factor=1.0,
                 count_reads=True):

        self.bam = bam
        self.bin_size = bin_size
        self.ins_mean = ins_mean
        self.ins_sd = ins_sd
        self.min_mapq = min_mapq
        self.subsample = subsample
        self.random_state = np.random.RandomState(random_seed)
        self.strong = strong
        self.max_site_dist = max_site_dist
        self.spacing_factor = spacing_factor
        self.count_reads = count_reads

        self.total_reads = None
        if count_reads:
            print 'Counting reads in bam file...'
            self.total_reads = bam.count(until_eof=True)

        self.active_seq = OrderedDict()
        self.total_len = 0
        for n, li in enumerate(bam.lengths):
            if li > ref_min_len:
                self.active_seq[bam.references[n]] = n
                self.total_len += li
        self.total_seq = len(self.active_seq)

        self.unrestricted = False if enz_name else True
        self._init_seq_info(enz_name, seq_file)

        print 'Map based upon mapping containing:\n' \
              '\t{0} sequences\n' \
              '\t{1}bp total length\n' \
              '\t{2} alignments'.format(self.total_seq, self.total_len,
                                        self.total_reads if self.count_reads else 'unknown')

        self.bin_count = int(self.total_len / self.bin_size) + 1
        print 'Map details:\n' \
              '\t{0}bp width\n' \
              '\t{1}x{1} dimension'.format(self.bin_size, self.bin_count)

        # a simple associative map for binning on whole sequences
        id_set = set(self.active_seq.values())
        self.seq_map = OrderedDict((n, OrderedDict(zip(id_set, [0]*len(id_set)))) for n in id_set)

        # initialized later
        self.raw_map = None

    def _init_seq_info(self, enz_name, seq_file):
        """
        Initialise information to do with the sites of the given template sequences.
        Each sequence is digested, the sites recorded and statistics collected. In parsing
        the references, the list of active sequences may change.
        :param enz_name: the name of the enzyme (case sensitive)
        :param seq_file: the multi-fasta to digest
        """
        self.seq_info = OrderedDict()

        if enz_name:
            self.enzyme = get_enzyme_instance(enz_name)
            print 'Cut-site spacing naive expectation: {0}'.format(4**self.enzyme.size)
        else:
            # non-RE based ligation-pairs
            self.enzyme = None
            print 'Warning: no enzyme was specified, therefore ligation pairs ' \
                  'are considered unconstrained in position.'

        offset = 0
        for rseq in SeqIO.parse(seq_file, 'fasta'):
            if rseq.id not in self.active_seq:
                continue
            xi = self.bam.references.index(rseq.id)

            if self.enzyme:
                sites = np.array(self.enzyme.search(rseq.seq, linear=False))

                sites.sort()  # pedantic check for order
                medspc = np.median(np.diff(sites))
                self.seq_info[xi] = {'sites': sites,
                                     'invtree': IntervalTree(Interval(si, si+1) for si in sites),
                                     'max_dist': self.spacing_factor * medspc,
                                     'offset': offset}
                print 'Found {0} cut-sites for \"{1}\" ' \
                      'med_spc: {2:.1f} max_dist: {3:.1f}'.format(enz_name, rseq.id, medspc,
                                                                  self.spacing_factor*medspc)
            else:
                self.seq_info[xi] = {'offset': offset}

            offset += len(rseq)

    def _init_map(self, dt=np.int32):
        print 'Initialising contact map of {0}x{0} from total extent of {1}bp over {2} sequences'.format(
            self.bin_count, self.total_len, self.total_seq)
        return np.zeros((self.bin_count, self.bin_count), dtype=dt)

    def build_sorted(self):
        import tqdm

        def _simple_match(r):
            return not r.is_unmapped and r.mapq >= _mapq

        def _strong_match(r):
            return strong_match(r, self.strong, True, _mapq)

        # set-up match call
        _matcher = _strong_match if self.strong else _simple_match

        def _is_toofar_absolute(_x, _abs_max, _medspc):
            return _abs_max and _x > _abs_max

        def _is_toofar_medspc(_x, _abs_max, _medspc):
            return _x > _medspc

        def _is_toofar_dummy(_x, _abs_max, _medspc):
            return False

        # set-up the test for site proximity
        if self.unrestricted:
            # all positions can participate, therefore restrictions are ignored.
            # TODO this would make more sense to be handled at arg parsing, since user should
            # be informed, rather tha silently ignoring.
            is_toofar = _is_toofar_dummy
        else:
            if self.max_site_dist:
                is_toofar = _is_toofar_absolute
            else:
                # default to spacing factor
                is_toofar = _is_toofar_medspc

        # initialise a map matrix for fine-binning and seq-binning
        self.raw_map = self._init_map()

        if self.total_reads:
            progress_bar = tqdm.tqdm(total=self.total_reads)
        else:
            progress_bar = tqdm.tqdm()

        # maximum separation being 3 std from mean
        _wgs_max = self.ins_mean + 2.0*self.ins_sd

        _hic_max = self.max_site_dist

        _mapq = self.min_mapq
        _map = self.raw_map
        _seq_info = self.seq_info
        _active_ids = set(self.active_seq.values())
        wgs_count = 0
        dropped_3c = 0
        kept_3c = 0

        sub_thres = self.subsample
        if self.subsample:
            uniform = self.random_state.uniform

        bam.reset()
        bam_iter = bam.fetch(until_eof=True)
        while True:

            try:
                r1 = bam_iter.next()
                progress_bar.update()
                while True:
                    # read records until we get a pair
                    r2 = bam_iter.next()
                    progress_bar.update()
                    if r1.query_name == r2.query_name:
                        break
                    r1 = r2
            except StopIteration:
                break

            if r1.reference_id not in _active_ids or r2.reference_id not in _active_ids:
                continue

            if sub_thres and sub_thres < uniform():
                continue

            if not _matcher(r1) or not _matcher(r2):
                continue

            if r1.is_read2:
                r1, r2 = r2, r1

            assume_wgs = False
            if r1.is_proper_pair:
                # x1 = r1.pos + r1.alen if r1.is_reverse else r1.pos
                # x2 = r2.pos + r2.alen if r2.is_reverse else r2.pos
                # if x2 < x1:
                #     x1, x2 = x2, x1
                #
                # if not _seq_info[r1.reference_id]['invtree'].overlaps(x1, x2+1):
                #     wgs_count += 1
                #     continue

                fwd, rev = (r2, r1) if r1.is_reverse else (r1, r2)
                ins_len = rev.pos + rev.alen - fwd.pos

                if fwd.pos <= rev.pos and ins_len < _wgs_max:
                    assume_wgs = True
                    wgs_count += 1
                    continue

            if not assume_wgs:

                if not self.unrestricted:

                    r1_dist = upstream_dist(r1, _seq_info[r1.reference_id]['sites'], r1.reference_length)
                    if is_toofar(r1_dist, _hic_max, _seq_info[r1.reference_id]['max_dist']):
                        dropped_3c += 1
                        continue

                    r2_dist = upstream_dist(r2, _seq_info[r2.reference_id]['sites'], r2.reference_length)
                    if is_toofar(r2_dist, _hic_max, _seq_info[r2.reference_id]['max_dist']):
                        dropped_3c += 1
                        continue

                kept_3c += 1
                self.seq_map[r1.reference_id][r2.reference_id] += 1

            r1pos = r1.pos if not r1.is_reverse else r1.pos + r1.alen
            r2pos = r2.pos if not r2.is_reverse else r2.pos + r2.alen

            ix1 = int((_seq_info[r1.reference_id]['offset'] + r1pos) / self.bin_size)
            ix2 = int((_seq_info[r2.reference_id]['offset'] + r2pos) / self.bin_size)

            if ix1 < ix2:
                _map[ix1][ix2] += 1
            else:
                _map[ix2][ix1] += 1

        print 'Assumed {0} were WGS pairs'.format(wgs_count)
        print 'Kept {0} and dropped {1} long-range (3C-ish) pairs'.format(kept_3c, dropped_3c)
        print 'Total raw map weight {0}'.format(np.sum(_map))

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

    def write_map(self, oname, _fmt='%d'):
        """
        Write a contact map as an ascii table to a file.
        :param oname: the file name to write
        :param _fmt: value format, defaults to int
        """
        np.savetxt(oname, self.raw_map, fmt=_fmt)

    def write_seqmap(self, oname):
        """
        Write a map based on entire sequences
        :param oname: output file name
        """
        with open(oname, 'w') as oh:
            rev_dict = dict(zip(self.active_seq.values(), self.active_seq.keys()))
            oh.write(',{0}\n'.format(','.join(self.active_seq.keys())))
            for i in self.seq_map:
                oh.write('{0},{1}\n'.format(rev_dict[i], ','.join([str(j) for j in self.seq_map[i].values()])))

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
    parser.add_argument('--spc-factor', default=None, type=float,
                        help='Maximum distance to nearest upstream site, as a factor '
                             'of each chromosome\'s median spacing')
    parser.add_argument('--max-dist', default=None, type=int,
                        help='Maximum distance to nearest upstream site in absolute terms')
    parser.add_argument('--ins-mean', default=500, type=int, help='Insert length mean [500]')
    parser.add_argument('--ins-sd', default=50, type=int, help='Insert length stddev [50]')
    parser.add_argument('--mapq', default=0, type=int, help='Minimum mapping quality [0]')
    parser.add_argument('--min-len', default=0, type=int, help='Minimum subject sequence length [none]')
    parser.add_argument('--per-contig', default=False, action='store_true', help='Bins are per contig')
    parser.add_argument('--bin-size', type=int, default=25000, help='Bin size in bp (25000)')
    parser.add_argument('--remove-diag', default=False, action='store_true',
                        help='Remove the central diagonal from plot')
    parser.add_argument('--skip-count', default=False, action='store_true', help='Skip initial read count for progress')
    parser.add_argument('-e', '--enzyme', default=None, help='Restriction enzyme (case sensitive)')
    parser.add_argument('refseq', metavar='FASTA', help='Reference fasta sequence (in same order)')
    parser.add_argument('bamfile', metavar='BAMFILE', help='Name-sorted BAM file')
    parser.add_argument('output', metavar='OUTPUT_BASE', help='Output base name')
    args = parser.parse_args()

    if args.spc_factor and args.max_dist:
        import sys
        print 'contactmap.py: error: spc-factor and max-dist are mutually exclusive options'
        sys.exit(1)

    with pysam.AlignmentFile(args.bamfile, 'rb') as bam:

        contacts = ContactMap(bam, args.enzyme, args.refseq, bin_size=args.bin_size, ins_mean=args.ins_mean,
                              ins_sd=args.ins_sd, ref_min_len=args.min_len, min_mapq=args.mapq,
                              subsample=args.sub_sample, strong=args.strong, max_site_dist=args.max_dist,
                              spacing_factor=1.0 if not args.spc_factor else args.spc_factor,
                              count_reads=not args.skip_count)

        contacts.build_sorted()

        print 'Writing output...'
        contacts.write_map('{0}.raw.cm'.format(args.output))
        contacts.write_seqmap('{0}.raw.sm'.format(args.output))
        contacts.plot_map('{0}.raw.png'.format(args.output), remove_diag=args.remove_diag)
