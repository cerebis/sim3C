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
import scipy.stats as st
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _clear_list(l):
    del l[:]


class EmpDist:

    FIRST = True
    SECOND = False
    HIGHEST_QUAL = 80
    MAX_DIST_NUMBER = 1.0e6
    CMB_SYMB = '.'
    A_SYMB = 'A'
    C_SYMB = 'C'
    G_SYMB = 'G'
    T_SYMB = 'T'
    N_SYMB = 'N'
    PRIMARY_SYMB = {A_SYMB, C_SYMB, G_SYMB, T_SYMB}
    ALL_SYMB = PRIMARY_SYMB | {CMB_SYMB, N_SYMB}

    # lookup table of probability indexed by quality score
    PROB_ERR = np.apply_along_axis(
        lambda xi: 10.**(-xi/10.), 0, np.arange(HIGHEST_QUAL)).tolist()

    def __init__(self, fname_first, fname_second, seq_quals):
        self.sep_quals = seq_quals
        self.qual_dist_first = dict(zip(EmpDist.ALL_SYMB, [[], [], [], [], [], []]))
        self.qual_dist_second = dict(zip(EmpDist.ALL_SYMB, [[], [], [], [], [], []]))
        self.dist_max = {EmpDist.FIRST: 0, EmpDist.SECOND: 0}
        self.init_dist(fname_first, fname_second)

    def init_dist(self, fname_first, fname_second):
        """
        Initialise the per-base_position distribution of quality scores
        :param fname_first: profile for first read
        :param fname_second: profile for second read
        """
        map(_clear_list, self.qual_dist_first.values())
        map(_clear_list, self.qual_dist_second.values())

        with open(fname_first, 'r') as hndl:
            self.read_emp_dist(hndl, True)

        with open(fname_second, 'r') as hndl:
            self.read_emp_dist(hndl, False)

        if not self.sep_quals:
            self.dist_max[EmpDist.FIRST] = len(self.qual_dist_first[self.CMB_SYMB])
            self.dist_max[EmpDist.SECOND] = len(self.qual_dist_second[self.CMB_SYMB])
        else:
            dist_len = np.array([len(self.qual_dist_first[k]) for k in self.PRIMARY_SYMB])
            assert np.all(dist_len == dist_len.max()), \
                'Invalid first profile, not all symbols represented over full range'
            self.dist_max[EmpDist.FIRST] = dist_len.max()

            dist_len = np.array([len(self.qual_dist_second[k]) for k in self.PRIMARY_SYMB])
            assert np.all(dist_len == dist_len.max()), \
                'Invalid second profile, not all symbols represented over full range'
            self.dist_max[EmpDist.SECOND] = dist_len.max()

    def verify_length(self, length, is_first):
        """
        Verify that profile and requested length are agreeable
        :param length: read length
        :param is_first: first or second read
        :return: True -- supported by profile
        """
        assert length < self.dist_max[is_first], 'Requested length exceeds that of profile'

    def get_read_qual(self, read_len, is_first):
        """
        Read qualities for a given read-length
        :param read_len: length of read to simulate
        :param is_first: first or second read
        :return: simulated qualities
        """
        self.verify_length(read_len, is_first)
        if is_first:
            return self._get_from_dist(self.qual_dist_first[self.CMB_SYMB], read_len)
        else:
            return self._get_from_dist(self.qual_dist_second[self.CMB_SYMB], read_len)

    @staticmethod
    def _get_from_dist(qual_dist_for_symb, read_len):
        """
        Generate simulated quality scores for a given length using the initialised
        distributions. Scores are related to the emporically determined CDFs specified
        at initialisation.
        :param qual_dist_for_symb: combined or separate symbols
        :param read_len: read length to simulate
        :return: simulated quality scores
        """

        sim_quals = []
        for i in xrange(read_len):
            # pull a number between 1 and max
            rnd_pick = Art.random_int(1, EmpDist.MAX_DIST_NUMBER + 1)
            # find the lowerbound for the distribution pertaining to this base 'i'
            # this determines the quality score
            qd = qual_dist_for_symb[i]
            lb = np.argwhere(qd[:, 0] >= rnd_pick)
            sim_quals.append(qd[lb[0], 1][0])
        return sim_quals

    """
    def get_read_qual_1st(self, seq, read_qual):
        read_len = len(seq)
        if read_len == 0:
            return False

        for i in xrange(read_len):
            cumCC = Art.random_int(1, EmpDist.MAX_DIST_NUMBER + 1)

            # if(seq[i]=='A'){
            #     it=a_qual_dist_first[i].lower_bound(cumCC);
            #     read_qual.push_back(it->second);
            # }
            # else if(seq[i]=='C'){
            #     it=c_qual_dist_first[i].lower_bound(cumCC);
            #     read_qual.push_back(it->second);
            # }
            # else if(seq[i]=='G'){
            #     it=g_qual_dist_first[i].lower_bound(cumCC);
            #     read_qual.push_back(it->second);
            # }
            # else if(seq[i]=='T'){
            #     it=g_qual_dist_first[i].lower_bound(cumCC);
            #     read_qual.push_back(it->second);
            # }
            # else{
            #     #return random quality less than 10
            #     read_qual.push_back((short) r_prob()*10);
            # }
        # }
        return True
    """

    def read_emp_dist(self, hndl, is_first):
        """
        Read an empirical distribution from a file.
        :param hndl: open file handle
        :param is_first: first or second read profile
        :return: True -- profile was not empty
        """
        n = 0
        while True:
            line = hndl.readline().strip()

            if not line:
                # end of file
                break
            if len(line) <= 0 or line.startswith('#'):
                # skip empty and comment lines
                continue

            tok = line.split('\t')
            symb, read_pos, values = tok[0], int(tok[1]), np.array(tok[2:], dtype=int)

            # skip lines pertaining to unrequested mode
            if self.sep_quals:
                if symb == self.CMB_SYMB or symb == self.N_SYMB:
                    # requested separate quals but this pertains to combined or N
                    continue
            else:  # if combined
                if symb != self.CMB_SYMB:
                    # requested combined but this pertains to separate
                    continue

            if read_pos != n:
                if read_pos != 0:
                    raise IOError('Error: invalid format in profile at [{0}]'.format(line))
                n = 0

            line = hndl.readline().strip()
            tok = line.split('\t')
            symb, read_pos, counts = tok[0], int(tok[1]), np.array(tok[2:], dtype=int)

            if read_pos != n:
                raise IOError('Error: invalid format in profile at [{0}]'.format(line))

            if len(values) != len(counts):
                raise IOError('Error: invalid format in profile at [{0}]'.format(line))

            dist = np.array([(cc, values[i]) for i, cc in
                             enumerate(np.ceil(counts * EmpDist.MAX_DIST_NUMBER/counts[-1]).astype(int))])

            if dist.size > 0:
                n += 1
                try:
                    if is_first:
                        self.qual_dist_first[symb].append(dist)
                    else:
                        self.qual_dist_second[symb].append(dist)
                except:
                    raise IOError('Error: unexpected base symbol [{0}] linked to distribution'.format(symb))

        return n != 0


class SeqRead:

    def __init__(self, read_len, max_num, ins_prob, del_prob):
        self.max_num = max_num
        self.read_len = read_len
        self.del_prob = del_prob
        self.ins_prob = ins_prob
        self.is_plus_strand = None
        self.seq_ref = []
        self.seq_read = []
        self.bpos = None
        self.indel = {}
        self.substitution = {}
        self.del_rate = self._set_rate(del_prob)
        self.ins_rate = self._set_rate(ins_prob)

    def __str__(self):
        return 'from {0}...{1}bp created {2}'.format(self.seq_ref[0:10], len(self.seq_ref), self.seq_read)

    def read_record(self, ref_id, read_num, quals):
        """
        Create a Biopython SeqRecord appropriate for writing to disk and matching the format
        generated by ART_illumina
        :param ref_id: sequence id of the mother sequence
        :param read_num: an index for this read
        :param quals: simulated quality scores for this read
        :return: Bio.SeqRecord
        """
        rec = SeqRecord(
                Seq(self._read_str(), IUPAC.ambiguous_dna),
                id=self._read_id(ref_id, read_num),
                description=self._read_desc())
        # seems the only means of adding quality scores to a SeqRecord
        rec.letter_annotations['phred_quality'] = quals
        return rec

    def _read_desc(self):
        """
        Create a string description for this read, suitable for inclusion if output
        :return: a string description
        """
        return '{0}{1}'.format(self.bpos, 'F' if self.is_plus_strand else 'R')

    @staticmethod
    def _read_id(ref_id, n):
        """
        Create an id for this read, based on the mother sequence and an index. This follows ART_illumina
        practice.
        :param ref_id: mother sequence id
        :param n: an index for the read
        :return: a string id for this read
        """
        return '{0}-{1}'.format(ref_id, n)

    def _read_str(self):
        """
        Create a string representation of this read's sequence. This is necessary as internally
        the sequence is handled as a list -- since strings are immutable in Python.
        :return:
        """
        return ''.join(self.seq_read)

    def _set_rate(self, prob):
        """
        Set the rate of an error type, returning a list of SeqRead.max_num length
        :param prob: probability of an error
        :return: list
        """
        rates = []
        if self.max_num > self.read_len:
            self.max_num = self.read_len
        for i in xrange(1, self.max_num+1):
            rates.append(1 - st.binom.cdf(i, self.read_len, prob))
        return rates

    def add_error(self, quals):
        """
        ART Illumina uses the simulated reads to potentially modify the read sequence.
        :param quals: simulated quality scores
        :return: number of modified bases.
        """
        assert len(quals) == len(self.seq_read), "The number of bases is not equal to the number of quality scores!\n" \
                                                 "qual size: {0},  read len: {1}".format(len(quals), len(self.seq_read))

        num = 0
        for i in xrange(len(quals)):
            if self.seq_read[i] == EmpDist.N_SYMB:
                quals[i] = 1
                continue

            # if we draw a number less than prob_err for a base, substitute it
            if Art.random_unif() < EmpDist.PROB_ERR[quals[i]]:
                achar = self.seq_read[i].upper()
                achar = Art.random_base(achar)
                self.seq_read[i] = achar
                self.substitution[i] = achar
                num += 1

        return num

    def clear(self):
        """
        Clear the working internal collections.
        """
        self.indel.clear()
        self.substitution.clear()
        del self.seq_ref[:]
        del self.seq_read[:]

    def get_indel(self, read_len):
        """
        Generate insertion and deletions
        :param read_len: length of read
        :return: net change in length, i.e. insertion_length - deletion_length
        """
        self.indel.clear()
        ins_len = 0
        del_len = 0

        # deletion
        for i in xrange(len(self.del_rate)-1, -1, -1):
            if self.del_rate[i] >= Art.random_unif():
                del_len = i+1
                j = i
                while j >= 0:
                    # invalid deletion positions: 0 or read_len-1
                    pos = Art.random_int(0, read_len)
                    if pos == 0:
                        continue
                    if pos not in self.indel:
                        self.indel[pos] = '-'
                        j -= 1
                break

        # insertion
        for i in xrange(len(self.ins_rate)-1, -1, -1):
            # ensure that enough unchanged position for mutation
            if read_len - del_len - ins_len < i+1:
                continue
            if self.ins_rate[i] >= Art.random_unif():
                ins_len = i+1
                j = i
                while j >= 0:
                    pos = Art.random_int(0, read_len)
                    if pos not in self.indel:
                        self.indel[pos] = Art.random_base()
                        j -= 1
                break

        return ins_len - del_len

    # number of deletions <= number of insertions
    def get_indel_2(self, read_len):
        """
        Second method for creating indels. Called in some situations when the first method
        as returned an unusable result.
        :param read_len: length of read
        :return: net change in length, i.e. insertion_length - deletion_length
        """

        # start over
        self.indel.clear()
        ins_len = 0
        del_len = 0

        for i in xrange(len(self.ins_rate)-1, -1, -1):
            if self.ins_rate[i] >= Art.random_unif():
                ins_len = i+1
                j = i
                while j >= 0:
                    pos = Art.random_int(0, read_len)
                    if pos not in self.indel:
                        self.indel[pos] = Art.random_base()
                        j -= 1
                break

        # deletion
        for i in xrange(len(self.del_rate)-1, -1, -1):
            if del_len == ins_len:
                break

            # ensure that enough unchanged position for mutation
            if read_len - del_len - ins_len < i+1:
                continue

            if self.del_rate[i] >= Art.random_unif():
                del_len = i+1
                j = i
                while j >= 0:
                    pos = Art.random_int(0, read_len)
                    if pos == 0:
                        continue
                    if pos not in self.indel:
                        self.indel[pos] = '-'
                        j -= 1
                break

        return ins_len - del_len

    def ref2read(self):
        """
        From the reference (mother) sequence, generating the read's sequence along
        with the indels.
        """
        if len(self.indel) == 0:
            # straight to an result if no indels, where here seq_ref
            # has already been chopped to the read length.
            self.seq_read = self.seq_ref

        else:
            # otherwise, we gotta a little more work to do.
            del self.seq_read[:]

            k = 0
            i = 0
            while i < len(self.seq_ref):
                if k not in self.indel:
                    self.seq_read += self.seq_ref[i]
                    i += 1
                    k += 1
                elif self.indel[k] == '-':
                    # deletion
                    i += 1
                    k += 1
                else:
                    # insertion
                    self.seq_read += self.indel[k]
                    k += 1

            while k in self.indel:
                self.seq_read += self.indel[k]
                k += 1


class Art:

    RANDOM_STATE = np.random.RandomState()

    COMPLEMENT_MAP = {EmpDist.A_SYMB: EmpDist.T_SYMB,
                      EmpDist.C_SYMB: EmpDist.G_SYMB,
                      EmpDist.G_SYMB: EmpDist.C_SYMB,
                      EmpDist.T_SYMB: EmpDist.A_SYMB,
                      EmpDist.N_SYMB: EmpDist.N_SYMB}

    def __init__(self, ref_seq, read_len, seq_read, emp_dist, seed=None):

        # check immediately that read lengths are possible for profile
        emp_dist.verify_length(read_len, True)
        emp_dist.verify_length(read_len, False)

        # initialise random state
        if seed:
            Art.RANDOM_STATE = np.random.RandomState(seed)

        # convert immutable string to list
        if isinstance(ref_seq, str):
            self.ref_seq = list(ref_seq.upper())
        else:
            self.ref_seq = ref_seq.upper()

        self.read_len = read_len
        self.valid_region = len(ref_seq) - read_len
        self.ref_seq_cmp = []
        self.emp_dist = emp_dist
        self.a_read = seq_read

        # complement reference
        for c in self.ref_seq:
            try:
                self.ref_seq_cmp += Art.COMPLEMENT_MAP[c]
            except KeyError:
                # unknown characters get assigned N
                self.ref_seq_cmp += EmpDist.N_SYMB

    def next_read_indel(self):
        """
        Create the next SeqRead and its accompanying quality scores
        :return: SeqRead, qual list
        """
        self.a_read.clear()

        # random position anywhere in valid range
        pos = Art.random_int(0, self.valid_region)
        slen = self.a_read.get_indel(self.read_len)

        # ensure the read fits within the extent of the reference
        if pos + self.read_len - slen > len(self.ref_seq):
            slen = self.a_read.get_indel_2(self.read_len)

        # pick a strand
        self.a_read.is_plus_strand = Art.random_coin_toss()

        if a_read.is_plus_strand:
            self.a_read.seq_ref = self.ref_seq[pos: pos+self.read_len-slen]
        else:
            self.a_read.seq_ref = self.ref_seq_cmp[pos: pos+self.read_len-slen]

        self.a_read.bpos = pos
        self.a_read.ref2read()

        # simulated quality scores from profiles
        quals = self.emp_dist.get_read_qual(self.read_len, True)
        # add quality to read, where it can result in sequence changes
        self.a_read.add_error(quals)

        return self.a_read, quals

    @staticmethod
    def random_int(low, high):
        """
        Random integer between low (incl) and high (excl)
        :param low:
        :param high:
        :return:
        """
        return Art.RANDOM_STATE.randint(low, high)

    @staticmethod
    def random_coin_toss():
        """
        :return: True or False
        """
        return Art.RANDOM_STATE.choice([True, False])

    @staticmethod
    def random_base(excl=None):
        """
        Return a random selection of A,C,G or T. If specified, exclude one of the four.
        :param excl: a base to exclude from the draw
        :return: a random base.
        """
        if not excl:
            return Art.RANDOM_STATE.choice(list(EmpDist.PRIMARY_SYMB))
        else:
            return Art.RANDOM_STATE.choice(list(set(EmpDist.PRIMARY_SYMB) - set(excl)))

    @staticmethod
    def random_unif():
        """
        :return: Uniformly random value between 0 and 1.
        """
        return Art.RANDOM_STATE.uniform()


if __name__ == '__main__':

    from Bio import SeqIO
    import argparse
    import math

    parser = argparse.ArgumentParser(description='Generate Illumina reads')
    parser.add_argument('-S', '--seed', type=int, default=None, help='Random seed')
    parser.add_argument('--profile1', help='ART sequencer profile for R1', required=True)
    parser.add_argument('--profile2', help='ART sequencer profile for R2', required=True)
    parser.add_argument('-l', '--read-length', type=int, help='Read length', required=True)
    parser.add_argument('-X', '--xfold', type=float, help='Depth of coverage')
    parser.add_argument('-N', '--num-reads', type=int, help='Number of reads')
    parser.add_argument('--ins-rate', type=float, default=0.00009, help='Insert rate')
    parser.add_argument('--del-rate', type=float, default=0.00011, help='Deletion rate')
    parser.add_argument('fasta', help='Reference fasta')
    parser.add_argument('output', help='Output fastq file')
    args = parser.parse_args()

    if args.xfold and args.num_reads:
        raise RuntimeError('xfold and num-reads are mutually exclusive options')

    with open(args.output, 'w') as out_h:

        for input_record in SeqIO.parse(args.fasta, 'fasta'):

            a_read = SeqRead(args.read_length, 2, args.del_rate, args.ins_rate)
            emp_dist = EmpDist(args.profile1, args.profile2, False)
            art = Art(str(input_record.seq), args.read_length, a_read, emp_dist, args.seed)

            if args.xfold:
                num_seq = int(math.ceil(len(art.ref_seq) / args.read_length * args.xfold))
            else:
                num_seq = args.num_reads

            print 'Generating {0} reads for {1}'.format(num_seq, input_record.id)

            print_rate = int(num_seq*100./10)

            for n in xrange(num_seq):

                # get next read and quals
                a_read, quals = art.next_read_indel()

                # create file record
                out_record = a_read.read_record(input_record.id, n, quals)

                SeqIO.write([out_record], out_h, 'fastq')

                if ((n+1)*100) % print_rate == 0:
                    print '\tWrote {0} reads'.format(n+1)
