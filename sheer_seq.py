import Bio.SeqIO as SeqIO
import numpy as np
import argparse
import os

def log_normal_intervals(extent, mean, sd, min_len, rs):
    invs = []
    inv_sum = 0
    mean = np.log(mean)
    sd = np.log(sd)

    # draw random intervals in batches of 10 for efficiency (maybe)
    done = False
    while not done:
        candidates = [v for v in rs.lognormal(mean, sd, size=10) if v > min_len]
        for ci in candidates:
            invs.append(ci)
            inv_sum += ci
            if inv_sum > extent:
                done = True
                break

    # working from smallest intervals, determine how many to remove
    # beyond extent to best agree.
    invs = sorted(invs, reverse=True)
    trim_i = -1
    while sum(invs[:trim_i]) > extent:
        trim_i -= 1

    # trim these small intervals and randomise
    invs = np.array(invs[:trim_i])
    rs.shuffle(invs)

    # convert the intervals to exactly cover the range 0-extent.
    return np.hstack(([0], np.cumsum(invs/invs.sum()*extent)[:-1], [extent])).astype(np.int)


def sheer_sequence(rs, seq_name, fasta_file, out_file, mean, sd, gap_size=1000, min_size=5000):

    index = SeqIO.index(fasta_file, 'fasta')
    rseq = index[seq_name]
    
    # divide the length uniform random and make into a set of interval pairs.
    invs = log_normal_intervals(len(rseq), mean, sd, min_size, rs)
    invs = np.vstack((np.arange(len(invs)-1), invs[:-1], invs[1:])).T
    
    np.savetxt('{}.txt'.format(os.path.splitext(out_file)[0]), invs, fmt='%d')

    frags = []
    for i, a, b in invs:
        if i == 0:
            b -= gap_size/2
        elif i == len(invs)-1:
            a += gap_size/2
        else:
            a += gap_size/2
            b -= gap_size/2
        
        frg = rseq[a:b]
        frg.id = '{}_{}-{}'.format(i, a, b)
        frg.description = rseq.description
        frags.append(frg)
        
    return frags


parser = argparse.ArgumentParser()
parser.add_argument('--seed', default=1234, type=int, help='Random seed')
parser.add_argument('--mean', type=int, default=2000, help='Log-normal mean size in bases (2000)')
parser.add_argument('--sd', type=int, default=5, help='Log-normal stddev in bases (5)')
parser.add_argument('--min', type=int, default=1000, help='Minimum size in bases (1000)')
parser.add_argument('--gap', type=int, default=50, help='Gap size in bases (50)')
parser.add_argument('id', help='Sequence id of target to sheer')
parser.add_argument('fasta', help='Input (multi)fasta sequence')
parser.add_argument('out', help='Output (multi)fasta sequence')
args = parser.parse_args()

rs = np.random.RandomState(args.seed)

# split sequence into pieces
frags = sheer_sequence(rs, args.id, args.fasta, args.out, mean=args.mean, sd=args.sd, min_size=args.min, 
                       gap_size=args.gap)


#shuffle the order
frags = np.array(frags)
rs.shuffle(frags)

# write out shuffled 
SeqIO.write(frags, args.out, 'fasta')

