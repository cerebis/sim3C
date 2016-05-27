#!/usr/bin/env python
from Bio import SeqIO
import numpy as np
import sys

def calc_N50L50(lengths):
    l = np.array(lengths)
    l.sort()
    cl = np.cumsum(l)
    xi = np.searchsorted(cl, l.sum()/2.0, side='right')
    return l[xi],len(l) - xi

if len(sys.argv) != 2:
    print 'Usage: [fasta-file]'
    sys.exit(0)

lengths = [len(si) for si in SeqIO.parse(sys.argv[1], 'fasta')]

n50, l50 = calc_N50L50(lengths)
print n50, l50

