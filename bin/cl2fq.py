#!/usr/bin/env python
from Bio import SeqIO
import truthtable as ttable
import os
import sys

if len(sys.argv) != 5:
    print 'usage <clid> <clustering> <multi-fassta> <out dir>'
    sys.exit(1)

clid = int(sys.argv[1])

print 'Getting sequence index'
seq_index = SeqIO.index(sys.argv[3], 'fasta')

print 'Reading clustering'
tt = ttable.read_mcl(sys.argv[2])
#tt.print_tally(50)
#tt.filter_extent(0.01, lengths)
cl2seq = tt.invert()

print 'Collecting sequences for cluster {0}'.format(clid)
seqs = []
for si in cl2seq[clid]:
    seqs.append(seq_index[si])

print 'Writing {0} sequences for cluster {1}'.format(len(seqs), clid)
if not os.path.exists(sys.argv[4]):
    os.makedirs(sys.argv[4])

SeqIO.write(seqs, os.path.join(sys.argv[4], 'cl{0}.fasta'.format(clid)), 'fasta')
