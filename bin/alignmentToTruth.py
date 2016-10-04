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

import numpy as np

import Psl
import truthtable as tt


def parse_psl(psl_file, min_id=0.90, cover_thres=0.96):
    """
    Calculate a truth table by considering the accumulated coverage of query sequences onto
    the references. The coverage is treated as a binary mask and the total extent a contig
    covers each reference determines the degree of ownership.

    Writes out a full truth table to user specified output path.

    Note: ideally the PSL file should be sorted for descending alignment score.

    :param psl_file:
    :param min_id: ignore alignments whose identity is less than this threshold
    :param cover_thres: query mean coverage threshold required for assignment to a ref be accepted.
    :return: None -- this method presently breaks the logical flow of this script.
    """
    with open(psl_file, 'r') as h_in:

        all_hits = 0
        rejected = 0

        aln_masks = {}

        # traverse alignment file, build up the masks for each query to reference[s] assocs.
        for aln in Psl.parse(h_in):

            all_hits += 1

            if aln.percent_id < min_id:
                rejected += 1
                continue

            if aln.q_name not in aln_masks:
                aln_masks[aln.q_name] = {}

            if aln.t_name not in aln_masks[aln.q_name]:
                aln_masks[aln.q_name][aln.t_name] = np.zeros(int(aln.q_size))

            per_id = aln.percent_id
            mask_slice = aln_masks[aln.q_name][aln.t_name][aln.q_start:aln.q_end+1]
            mask_slice[np.where(mask_slice < per_id)] = per_id

        # build dictionary of assignments and weights
        truth = {}
        weights = {}
        for n, ti in enumerate(aln_masks):
            masks = np.vstack(aln_masks[ti].values())
            names = np.array(aln_masks[ti].keys())
            covers = np.mean(masks, 1)
            idx = np.where(covers > cover_thres)
            if idx[0].shape[0] > 0:
                truth[ti] = {}
                for i in np.nditer(idx):
                    truth[ti][str(names[i])] = float(covers[i])
            weights[ti] = masks.shape[1]

        # initialize truthtable
        ttable = tt.TruthTable()
        ttable.update(truth, weights)

        return ttable


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Determine the location of query sequences on the reference')
    parser.add_argument('--ofmt', choices=['yaml', 'json'], default='json',
                        help='Output format [json]')
    parser.add_argument('--simple', default=False, action='store_true',
                        help='Generate a simple 1:1 mapping (hard) truth table')
    parser.add_argument('--min-id', type=float, default=0.90,
                        help='Minimum percentage identity for alignment')
    parser.add_argument('--min-cov', type=float, default=0.96,
                        help='Minimum coverage of query by alignment')
    parser.add_argument('alignment', metavar='ALIGNMENT_FILE',
                        help='Last generated PSL format with query sequences aligned to reference')
    parser.add_argument('output', metavar='OUTPUT', help='Output file name')
    args = parser.parse_args()

    ttable = parse_psl(args.alignment, min_id=args.min_id, cover_thres=args.min_cov)

    if args.simple:
        # Reduce to hard 1:1 mapping
        ttable.write_hard(args.output, fmt=args.ofmt)

    else:
        # Keep all assignments of queries to references.
        ttable.write(args.output, fmt=args.ofmt)
