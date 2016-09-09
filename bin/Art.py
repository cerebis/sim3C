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
import scipy.stats as st
import numpy as np


RANDOM_STATE = np.random.RandomState(1234)


def _clear_list(l):
    del l[:]


def _random_int(low, high):
    return RANDOM_STATE.randint(low, high)


def _random_coin_toss():
    return RANDOM_STATE.choice([True, False])


def _random_base(excl=None):
    if not excl:
        return RANDOM_STATE.choice(list(EmpDist.PRIMARY_SYMB))
    else:
        return RANDOM_STATE.choice(list(set(EmpDist.PRIMARY_SYMB) - set(excl)))


def _random_unif():
    return RANDOM_STATE.uniform()


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

    def __init__(self):
        self.sep_quals = False
        self.qual_dist_first = dict(zip(EmpDist.ALL_SYMB, [[], [], [], [], [], []]))
        self.qual_dist_second = dict(zip(EmpDist.ALL_SYMB, [[], [], [], [], [], []]))
        self.dist_max = {EmpDist.FIRST: 0, EmpDist.SECOND: 0}

    def set_dist(self, fname_first, fname_second, sep_quals):
        self.sep_quals = sep_quals
        map(_clear_list, self.qual_dist_first.values())
        map(_clear_list, self.qual_dist_second.values())

        with open(fname_first, 'r') as hndl:
            self.read_emp_dist(hndl, True)

        with open(fname_second, 'r') as hndl:
            self.read_emp_dist(hndl, False)

        if sep_quals:
            self.dist_max[EmpDist.FIRST] = len(self.qual_dist_first[self.CMB_SYMB])
            self.dist_max[EmpDist.SECOND] = len(self.qual_dist_second[self.CMB_SYMB])
        else:
            dlen = np.array([len(self.qual_dist_first[k]) for k in self.PRIMARY_SYMB])
            assert np.all(dlen == dlen.max()), 'Invalid first profile, not all symbols represented over full range'
            self.dist_max[EmpDist.FIRST] = dlen.max()

            dlen = np.array([len(self.qual_dist_second[k]) for k in self.PRIMARY_SYMB])
            assert np.all(dlen == dlen.max()), 'Invalid second profile, not all symbols represented over full range'
            self.dist_max[EmpDist.SECOND] = dlen.max()

    def _verify_length(self, length, is_first):
        # verify that profile and requested length are agreeable
        assert self.dist_max[is_first] < length, 'Requested length exceeds that of profile'

    def get_read_qual(self, read_qual, read_len, is_first):
        self._verify_length(read_len, is_first)
        if is_first:
            return self.get_from_dist(self.qual_dist_first[self.CMB_SYMB], read_qual, read_len)
        else:
            return self.get_from_dist(self.qual_dist_second[self.CMB_SYMB], read_qual, read_len)

    @staticmethod
    def get_from_dist(qual_dist_for_symb, read_qual, read_len):
        # empty the list
        del read_qual[:]

        for i in xrange(read_len):
            # pull a number between 1 and max
            rnd_pick = _random_int(1, EmpDist.MAX_DIST_NUMBER+1)
            # find the lowerbound for the distribution pertaining to this base 'i'
            # this determines the quality score
            # Note, this is slow as we're building this array every call!
            # qd = np.array(qual_dist_for_symb[i].items())
            qd = qual_dist_for_symb[i]
            lb = np.argwhere(qd[:, 0] >= rnd_pick)
            read_qual.append(qd[lb[0], 1][0])

        return True

    def get_read_qual_1st(self, seq, read_qual):
        read_len = len(seq)
        if read_len == 0:
            return False

        for i in xrange(read_len):
            cumCC = _random_int(1, EmpDist.MAX_DIST_NUMBER+1)

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

    def read_emp_dist(self, hndl, is_first):
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

    def __init__(self):
        self.is_plus_strand = None
        self.seq_ref = []
        self.seq_read = []
        self.bpos = None
        self.indel = {}
        self.substitution = {}
        self.del_rate = []
        self.ins_rate = []

    def __str__(self):
        return 'from {0}...{1}bp created {2}'.format(self.seq_ref[0:10], len(self.seq_ref), self.seq_read)

    @staticmethod
    def set_rate(read_len, prob, max_num, rates):
        del rates[:]
        if max_num > read_len:
            max_num = read_len
        for i in xrange(1, max_num+1):
            rates.append(1 - st.binom.cdf(i, read_len, prob))

    def add_error(self, quals):

        if len(quals) != len(self.seq_read):
            raise RuntimeError("Error: the number of bases is not equal to the number of quality scores!\n" +
                               "qual size: {0},  read len: {1}".format(len(quals), len(self.seq_read)))

        num = 0
        for i in xrange(len(quals)):
            if self.seq_read[i] == EmpDist.N_SYMB:
                quals[i] = 1
                continue

            # if we draw a number less than prob_err for a base, substitute it
            if _random_unif() < EmpDist.PROB_ERR[quals[i]]:
                achar = self.seq_read[i].upper()
                achar = _random_base(achar)
                self.seq_read[i] = achar
                self.substitution[i] = achar
                num += 1

        return num

    def clear(self):
        self.indel.clear()
        self.substitution.clear()
        del self.seq_ref[:]
        del self.seq_read[:]

    def get_indel(self, read_len):
        self.indel.clear()
        ins_len = 0
        del_len = 0

        # deletion
        for i in xrange(len(self.del_rate)-1, -1, -1):
            if self.del_rate[i] >= _random_unif():
                del_len = i+1
                j = i
                while j >= 0:
                    # invalid deletion positions: 0 or read_len-1
                    pos = _random_int(0, read_len)
                    if pos == 0:
                        continue
                    if pos not in self.indel:
                        self.indel[pos] = '-'
                        j -= 1
                break

        for i in xrange(len(self.ins_rate)-1, -1, -1):

            # ensure that enough unchanged position for mutation
            if read_len - del_len - ins_len < i+1:
                continue

            if self.ins_rate[i] >= _random_unif():
                ins_len = i+1
                j = i
                while j >= 0:
                    pos = _random_int(0, read_len)
                    if pos not in self.indel:
                        self.indel[pos] = _random_base()
                        j -= 1
                break

        return ins_len - del_len

    # number of deletions <= number of insertions
    def get_indel_2(self, read_len):
        self.indel.clear()
        ins_len = 0
        del_len = 0

        for i in xrange(len(self.ins_rate)-1, -1, -1):
            if self.ins_rate[i] >= _random_unif():
                ins_len = i+1
                j = i
                while j >= 0:
                    pos = _random_int(0, read_len)
                    if pos not in self.indel:
                        self.indel[pos] = _random_base()
                        j -= 1
                break

        # deletion
        for i in xrange(len(self.del_rate)-1, -1, -1):
            if del_len == ins_len:
                break

            # ensure that enough unchanged position for mutation
            if read_len - del_len - ins_len < i+1:
                continue

            if self.del_rate[i] >= _random_unif():
                del_len = i+1
                j = i
                while j >= 0:
                    pos = _random_int(0, read_len)
                    if pos == 0:
                        continue
                    if pos not in self.indel:
                        self.indel[pos] = '-'
                        j -= 1
                break

        return ins_len - del_len

    def ref2read(self):

        if len(self.indel) == 0:
            self.seq_read = self.seq_ref
            return

        del self.seq_read[:]

        k = 0
        i = 0
        while i < len(self.seq_ref):
            if k not in self.indel:
                self.seq_read += self.seq_ref[i]
                i += 1
                k += 1
            elif self.indel[k] == '-':
                i += 1
                k += 1
            else:
                self.seq_read += self.indel[k]
                k += 1

        while k in self.indel:
            self.seq_read += self.indel[k]
            k += 1


class Art:

    COMPLEMENT_MAP = {EmpDist.A_SYMB: EmpDist.T_SYMB,
                      EmpDist.C_SYMB: EmpDist.G_SYMB,
                      EmpDist.G_SYMB: EmpDist.C_SYMB,
                      EmpDist.T_SYMB: EmpDist.A_SYMB,
                      EmpDist.N_SYMB: EmpDist.N_SYMB}

    def __init__(self, ref_seq, read_len):
        if isinstance(ref_seq, str):
            self.ref_seq = list(ref_seq.upper())
        else:
            self.ref_seq = ref_seq.upper()
        self.read_len = read_len
        self.valid_region = len(ref_seq) - read_len
        self.ref_seq_cmp = []

        for n, c in enumerate(self.ref_seq):
            try:
                self.ref_seq_cmp += Art.COMPLEMENT_MAP[c]
                if n % 100000 == 0:
                    print n
            except KeyError:
                # unknown characters get assigned N
                self.ref_seq_cmp += EmpDist.N_SYMB

    def next_read_indel(self, a_read):

        # pos in [0 ..len-1]
        pos = _random_int(0, self.valid_region)
        slen = a_read.get_indel(self.read_len)

        # ensure get a fixed read length
        if pos + self.read_len - slen > len(self.ref_seq):
            slen = a_read.get_indel_2(self.read_len)

        a_read.is_plus_strand = _random_coin_toss()

        if a_read.is_plus_strand:
            a_read.seq_ref = self.ref_seq[pos: pos+self.read_len-slen]
        else:
            a_read.seq_ref = self.ref_seq_cmp[pos: pos+self.read_len-slen]

        a_read.bpos = pos
        a_read.ref2read()

        return True

if __name__ == '__main__':

    from Bio import SeqIO
    import argparse
    pass
