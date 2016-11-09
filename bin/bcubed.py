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
import numpy as np


def filter_shared(d, shared):
    unshared = set(d.keys()) - shared
    for oi in unshared:
        del d[oi]


def inv_dict(d):
    cl = set()
    for v in d.values():
        cl.update(v)
    inv_d = {ci: set() for ci in cl}
    for oi, cl in d.iteritems():
        for ci in cl:
            inv_d[ci].add(oi)
    return inv_d


def find_overlap_set(oi, inv_map):
    ovl_set = set()
    for cl, obj_set in inv_map.iteritems():
        if oi in obj_set:
            ovl_set |= obj_set
    return ovl_set


def multi_pre(t, c, oi, oj):
    c_size = len(c[oi] & c[oj])
    t_size = len(t[oi] & t[oj])
    return min(t_size, c_size) / float(c_size)


def multi_rec(t, c, oi, oj):
    c_size = len(c[oi] & c[oj])
    t_size = len(t[oi] & t[oj])
    return min(t_size, c_size) / float(t_size)


def bcubed_recall(td, cd, w=None):
    i_td = inv_dict(td)
    R = 0.0
    if w:
        for oi in td:
            ovl = find_overlap_set(oi, i_td)
            ri = 0.0
            for oj in ovl:
                ri += w[oj] * multi_rec(td, cd, oi, oj)
            W = np.sum([w[oj] for oj in ovl])
            R += ri / W
    else:
        for oi in td:
            ovl = find_overlap_set(oi, i_td)
            ri = 0.0
            for oj in ovl:
                ri += multi_rec(td, cd, oi, oj)
            R += ri / len(ovl)
    return R/len(td)


def bcubed_precison(td, cd, w=None):
    i_cd = inv_dict(cd)
    P = 0.0
    if w:
        for oi in cd:
            ovl = find_overlap_set(oi, i_cd)
            pi = 0.0
            for oj in ovl:
                pi += w[oj] * multi_pre(td, cd, oi, oj)
            W = np.sum([w[oj] for oj in ovl])
            P += pi / W
    else:
        for oi in cd:
            ovl = find_overlap_set(oi, i_cd)
            pi = 0.0
            for oj in ovl:
                pi += multi_pre(td, cd, oi, oj)
            P += pi / len(ovl)
    return P/len(cd)


def bcubed_F(td, cd, weights=None):
    cd = cd.copy()
    td = td.copy()
    shared = set(td.keys()) & set(cd.keys())

    # ratio of shared to all objects in analysis
    completeness = len(shared) / float(len(set(cd.keys()) | set(td.keys())))

    filter_shared(cd, shared)
    filter_shared(td, shared)

    pre = bcubed_precison(td, cd, weights)
    rec = bcubed_recall(td, cd, weights)

    return {'pre': pre, 'rec': rec, 'f': 2.0*pre*rec / (pre+rec), 'completeness': completeness}


if __name__ == '__main__':

    import truthtable as tt
    import io_utils
    import argparse
    import sys

    def write_msg(stream, msg):
        stream.write(msg + '\n')

    parser = argparse.ArgumentParser(description='Calculate extended bcubed metric')
    parser.add_argument('truth', metavar='TRUTH', help='Truth table (json format)')
    parser.add_argument('pred', metavar='PREDICTION', help='Prediction table (mcl format)')
    parser.add_argument('output', metavar='OUTPUT', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help='Output file (stdout)')
    parser.add_argument('--tfmt', choices=['json', 'yaml'], default='json',
                        help='Data format of truth table [json]')
    parser.add_argument('--ofmt', choices=['plain', 'json', 'yaml'], default='plain',
                        help='Output format [plain]')
    parser.add_argument('-w', '--weighted', default=False, action='store_true', help='Use truth object weights')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Verbose output')
    parser.add_argument('--hard', action='store_true', default=False, help='Extract hard truth prior to analysis')
    args = parser.parse_args()

    try:
        # read truth and convert to basic soft table
        truth = tt.read_truth(args.truth, args.tfmt)
        if len(truth) == 0:
            raise RuntimeError('Truth table contains no assignments: {0}'.format(args.truth))

        # collect object weights if requested
        weights = truth.get_weights() if args.weighted else None

        if args.verbose:
            print 'Truth Statistics'
            truth.print_tally()

        # convert to a plain dict representation, either soft (1:*) or hard (1:1)
        truth = truth.hard(True) if args.hard else truth.soft(True)

        # read clustering and convert to basic soft table
        clustering = tt.read_mcl(args.pred)
        if len(clustering) == 0:
            raise RuntimeError('Clustering contains no assignments: {0}'.format(args.pred))

        if args.verbose:
            print 'Clustering Statistics'
            clustering.print_tally()
        clustering = clustering.soft(True)

    except RuntimeError as er:
        print er.message
        sys.exit(1)

    result = bcubed_F(truth, clustering, weights)
    io_utils.write_to_stream(args.output, result, fmt=args.ofmt)
