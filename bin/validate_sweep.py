#!/usr/bin/env python
from collections import Counter
import fnmatch
import os
import numpy
import json

def find_files(dir, pattern):
    found = []
    for root, dir_list, file_list in os.walk(dir):
        for fn in fnmatch.filter(file_list, pattern):
            found.append(fn)
    return found


def remove_suffix(file_list, suffix):
    return [fn.split(suffix)[0] for fn in file_list]


def split_count(file_list, sep='-+-'):
    counts = {}
    for fn in file_list:
        levels = fn.split(sep)
        for n, li in enumerate(levels):
            if n not in counts:
                counts[n] = []
            counts[n].append(li)
    return {li + 1: Counter(cni) for li, cni in counts.iteritems()}


def validate_combinations(counters):
    # predicted number of combinations, as the product of observed values per level
    # obviously, one can predic the total from the specified sweep beforehand
    # this is just a convenience.
    goal_by_obs = numpy.prod([len(ci.keys()) for ci in counters.values()])
    print '  Goal by product of observed: {0}'.format(goal_by_obs)

    its = [sum(ci.values()) for ci in counters.values()]
    print '  Iterations per level: {0}'.format(its)
    uniq_its = set(its)
    if len(uniq_its) > 1:
        print '  Iteration count was inconsistent: {0}'.format(uniq_its)
        print '  Something is wrong in sweep'
    else:
        print '  Iteration count was consistent: {0}'.format(uniq_its.pop())


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='HiC combinatoric validation of output files')
    parser.add_argument('-s', '--separator', default='-+-', help='Sweep separator [-+-]')
    parser.add_argument('workdir', help='Working directory of sweep')
    args = parser.parse_args()

    """ HiC Sweep suffixes """
    hic_suffixes = [
        '.evo.fa',
        '.hic.fa.gz',
        '.graphml',
        '.hic2ctg.bam',
        '.graphml',
        '.contigs.fasta',
        '.truth',
        '.wgs2ctg.bam',
        '.wgs2ctg.cov',
        '.wgs.r1.fq.gz',
        '.wgs.r2.fq.gz']

    for suf in hic_suffixes:
        print 'For suffix: {0}'.format(suf)
        file_list = find_files(args.workdir, '*' + suf)
        if len(file_list) <= 0:
            print '  No files found'
            print
            continue

        file_list = remove_suffix(file_list, suf)

        counters = split_count(file_list, sep=args.separator)
        print 'Validation:'
        validate_combinations(counters)
        print 'Raw combination counts:'
        for li, ci in counters.iteritems():
            print '  Level', li, json.dumps(ci, sort_keys=True, indent=4)
        print