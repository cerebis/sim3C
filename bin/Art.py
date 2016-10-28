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
import string
import os

import numpy as np
import scipy.stats as st
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
The following module was transcribed and adapted from the original project's C++ code:

ART -- Artificial Read Transcription, Illumina Q version
Authors: Weichun Huang 2008-2016
License: GPL v3
"""

def _clear_list(l):
    del l[:]


# path of this (Art.py) source file. It is expected that profiles are co-located with Art.
MODULE_PATH = os.path.dirname(os.path.abspath(__file__))

# A catalog of empirical profiles for Illumina machine types.
ILLUMINA_PROFILES = {
    'Emp100': ('Illumina_profiles/Emp100R1.txt',
               'Illumina_profiles/Emp100R2.txt'),
    'Emp36':('Illumina_profiles/Emp36R1.txt',
             'Illumina_profiles/Emp36R2.txt'),
    'Emp44': ('Illumina_profiles/Emp44R1.txt',
              'Illumina_profiles/Emp44R2.txt'),
    'Emp50': ('Illumina_profiles/Emp50R1.txt',
              'Illumina_profiles/Emp50R2.txt'),
    'Emp75': ('Illumina_profiles/Emp75R1.txt',
              'Illumina_profiles/Emp75R2.txt'),
    'EmpMiSeq250': ('Illumina_profiles/EmpMiSeq250R1.txt',
                    'Illumina_profiles/EmpMiSeq250R2.txt'),
    'EmpR36': ('Illumina_profiles/EmpR36R1.txt',
               'Illumina_profiles/EmpR36R2.txt'),
    'EmpR44': ('Illumina_profiles/EmpR44R1.txt',
               'Illumina_profiles/EmpR44R2.txt'),
    'EmpR50': ('Illumina_profiles/EmpR50R1.txt',
               'Illumina_profiles/EmpR50R2.txt'),
    'EmpR75': ('Illumina_profiles/EmpR75R1.txt',
               'Illumina_profiles/EmpR75R2.txt'),
    'HiSeq2500L125': ('Illumina_profiles/HiSeq2500L125R1.txt',
                      'Illumina_profiles/HiSeq2500L125R2.txt'),
    'HiSeq2500L150': ('Illumina_profiles/HiSeq2500L150R1.txt',
                      'Illumina_profiles/HiSeq2500L150R2.txt'),
    'HiSeq2500L150filt': ('Illumina_profiles/HiSeq2500L150R1filter.txt',
                          'Illumina_profiles/HiSeq2500L150R2filter.txt'),
    'HiSeq2kL100': ('Illumina_profiles/HiSeq2kL100R1.txt',
                    'Illumina_profiles/HiSeq2kL100R2.txt'),
    'HiSeqXPCRfreeL150': ('Illumina_profiles/HiSeqXPCRfreeL150R1.txt',
                          'Illumina_profiles/HiSeqXPCRfreeL150R2.txt'),
    'HiSeqXtruSeqL150': ('Illumina_profiles/HiSeqXtruSeqL150R1.txt',
                         'Illumina_profiles/HiSeqXtruSeqL150R2.txt'),
    'MiSeqv3L250': ('Illumina_profiles/MiSeqv3L250R1.txt',
                    'Illumina_profiles/MiSeqv3L250R2.txt'),
    'NextSeq500v2L75': ('Illumina_profiles/NextSeq500v2L75R1.txt',
                        'Illumina_profiles/NextSeq500v2L75R2.txt')
}

def get_profile(name):
    """
    Return the absolute path to a requested Illumina profile.
    :param name: the name of the profile.
    :return: absolute (full) path
    """
    assert name in ILLUMINA_PROFILES, 'Unknown profile name. Try one of: {0}'.format(
        ', '.join(ILLUMINA_PROFILES.keys()))

    return map(lambda pi: os.path.join(MODULE_PATH, pi), ILLUMINA_PROFILES[name])

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

    @staticmethod
    def create(name, sep_quals=False):
        """
        Instantiate a EmpDist with the specified profile name.
        :param name: empirically derived machine profile
        :param sep_quals: independent quality model per base (A,C,G,T)
        :return: instance of EmpDist
        """
        profile_r1, profile_r2 = get_profile(name)
        return EmpDist(profile_r1, profile_r2, sep_quals)

    def __init__(self, fname_first, fname_second, sep_quals=False):
        """
        :param fname_first: file name of first read profile
        :param fname_second: file name of second read profile
        :param sep_quals: independent quality model per base (A,C,G,T)
        """
        self.sep_quals = sep_quals
        assert not sep_quals, 'Separate base qualities are not currently supported by this python implementation'
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
        assert length <= self.dist_max[is_first], 'Requested length exceeds that of profile'

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
        Generate simulated quality scores for a given length using an initialised
        distribution. Scores are related to the emporically determined CDFs specified
        at initialisation.
        :param qual_dist_for_symb: combined or separate symbols
        :param read_len: read length to simulate
        :return: simulated quality scores
        """

        sim_quals = []
        for pos in xrange(read_len):
            # draw a random integer
            rnd_pick = Art.random_int(1, EmpDist.MAX_DIST_NUMBER + 1)
            # for the quality score distribution at this position,
            # find the lower bound for this random number
            qd = qual_dist_for_symb[pos]
            lb = np.argwhere(qd[:, 0] >= rnd_pick)
            sim_quals.append(qd[lb[0], 1][0])
        return sim_quals

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

    def __init__(self, read_len, ins_prob, del_prob, max_num=2, plus_strand=None):
        self.max_num = max_num
        self.read_len = read_len
        self.del_prob = del_prob
        self.ins_prob = ins_prob
        self.is_plus_strand = plus_strand
        self.seq_ref = []
        self.seq_read = []
        self.quals = None
        self.bpos = None
        self.indel = {}
        self.substitution = {}
        self.del_rate = self._set_rate(del_prob)
        self.ins_rate = self._set_rate(ins_prob)

    def __str__(self):
        return 'from {0}...{1}bp created {2}'.format(self.seq_ref[0:10], len(self.seq_ref), self.seq_read)

    def read_record(self, seq_id, desc=''):
        """
        Create a Biopython SeqRecord appropriate for writing to disk and matching the format
        generated by ART_illumina
        :param seq_id: sequence id for read
        :param desc: sequence description
        :return: Bio.SeqRecord
        """
        rec = SeqRecord(
                Seq(self._read_str(), IUPAC.ambiguous_dna),
                id=seq_id,
                description=desc)
        # seems the only means of adding quality scores to a SeqRecord
        rec.letter_annotations['phred_quality'] = self.quals
        return rec

    def _read_desc(self):
        """
        Create a string description for this read, suitable for inclusion if output
        :return: a string description
        """
        return '{0}{1}'.format(self.bpos, 'F' if self.is_plus_strand else 'R')

    @staticmethod
    def read_id(ref_id, n):
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

    def parse_error(self):
        """
        When analyzed, sequences are potentially modified by the simulated quality scores.
        :return: number of modified bases.
        """
        assert self.quals, 'Quality scores have not been initialized for the read'
        assert len(self.quals) == self.length(), \
            "The number of bases is not equal to the number of quality scores!\n" \
            "qual size: {0},  read len: {1}".format(len(self.quals), self.length())

        num = 0
        for i in xrange(len(self.quals)):
            # if we encounter an undefined base, its quality score goes to 1.
            if self.seq_read[i] == EmpDist.N_SYMB:
                self.quals[i] = 1
                continue

            # if we draw a number less than prob_err for a base, a substitution occurs there.
            if Art.random_unif() < EmpDist.PROB_ERR[self.quals[i]]:
                ch = self.seq_read[i]
                ch = Art.random_base(ch)
                self.seq_read[i] = ch
                self.substitution[i] = ch
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

    def get_indel(self):
        """
        Generate insertion and deletions
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
                    pos = Art.random_int(0, self.read_len)
                    if pos == 0:
                        continue
                    if pos not in self.indel:
                        self.indel[pos] = '-'
                        j -= 1
                break

        # insertion
        for i in xrange(len(self.ins_rate)-1, -1, -1):
            # ensure that enough unchanged position for mutation
            if self.read_len - del_len - ins_len < i+1:
                continue
            if self.ins_rate[i] >= Art.random_unif():
                ins_len = i+1
                j = i
                while j >= 0:
                    pos = Art.random_int(0, self.read_len)
                    if pos not in self.indel:
                        self.indel[pos] = Art.random_base()
                        j -= 1
                break

        return ins_len - del_len

    # number of deletions <= number of insertions
    def get_indel_2(self):
        """
        Second method for creating indels. Called in some situations when the first method
        as returned an unusable result.
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
                    pos = Art.random_int(0, self.read_len)
                    if pos not in self.indel:
                        self.indel[pos] = Art.random_base()
                        j -= 1
                break

        # deletion
        for i in xrange(len(self.del_rate)-1, -1, -1):
            if del_len == ins_len:
                break

            # ensure that enough unchanged position for mutation
            if self.read_len - del_len - ins_len < i+1:
                continue

            if self.del_rate[i] >= Art.random_unif():
                del_len = i+1
                j = i
                while j >= 0:
                    pos = Art.random_int(0, self.read_len)
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

    def length(self):
        """
        Return the actual length of the simulation result. This can be shorter than the requested
        length "read_len" due to short templates.
        :return: length of actual simulated sequence
        """
        return len(self.seq_read)


class Art:

    RANDOM_STATE = np.random.RandomState()

    # translation table, non-standard bases become N
    COMPLEMENT_TABLE = string.maketrans('acgtumrwsykvhdbnACGTUMRWSYKVHDBN',
                                        'TGCAAnnnnnnnnnnnTGCAANNNNNNNNNNN')

    def __init__(self, read_len, emp_dist, ins_prob, del_prob, max_num=2, seed=None, ref_seq=None):

        # check immediately that read lengths are possible for profile
        emp_dist.verify_length(read_len, True)
        emp_dist.verify_length(read_len, False)

        # initialise random state
        if seed:
            Art.RANDOM_STATE = np.random.RandomState(seed)

        # convert immutable string to list
        if ref_seq:
            self.ref_seq = Art.make_mutable(ref_seq)
            self.ref_seq_cmp = list(Art.revcomp(ref_seq))
            self.valid_region = len(ref_seq) - read_len
        else:
            print 'Warning: no reference supplied, calls will have to supply a template'
            self.ref_seq = None
            self.ref_seq_cmp = None

        self.read_len = read_len
        self.emp_dist = emp_dist
        self.ins_prob = ins_prob
        self.del_prob = del_prob
        self.max_num = max_num

    @staticmethod
    def make_mutable(input_seq):
        """
        Pedantic method to assure that an input sequence is a mutable uppercase collection.
        :param input_seq: input sequence to make mutable
        :return: list of ch
        """
        assert isinstance(input_seq, str), 'Error: supplied sequences must be plain strings.'
        return list(input_seq.upper())

    @staticmethod
    def revcomp(seq):
        """
        Reverse complement a string representation of a sequence. This uses string.translate.
        :param seq: input sequence as a string
        :return: revcomp sequence as a string
        """
        return seq.translate(Art.COMPLEMENT_TABLE)[::-1]

    def _new_read(self, rlen=None, plus_strand=True):
        """
        Create a new read object ready for simulation.
        :param rlen: a read length other than what was defined when instantiating Art.
        :param plus_strand: True - forward strand, False - reverse strand
        :return: a new read object
        """
        if not rlen:
            return SeqRead(self.read_len, self.ins_prob, self.del_prob, self.max_num, plus_strand=plus_strand)
        else:
            return SeqRead(rlen, self.ins_prob, self.del_prob, self.max_num, plus_strand=plus_strand)

    def next_pair_simple_seq(self, template):
        """
        Get a fwd/rev pair of simple error-free reads for a template, where each read is sequenced off the ends.
        :param template: the target tempalte to sequencing fwd/rev
        :return: a dict {'fwd': SeqRead, 'rev': SeqRead}
        """
        return {'fwd': self.next_read_simple_seq(template, True),
                'rev': self.next_read_simple_seq(template, False)}

    def next_read_simple_seq(self, template, plus_strand, qual_val=40):
        """
        Generate a simple error-free read and constant quality values.
        :param template: the target template to sequence
        :param plus_strand: forward: True, reverse: False
        :param qual_val: value of constant quality scores
        :return: SeqRead
        """
        read = self._new_read(plus_strand=plus_strand)
        if len(template) < read.read_len:
            # for templates shorter than the requested length, we sequence its total extent
            read.read_len = len(template)

        if read.is_plus_strand:
            read.seq_ref = template[0: self.read_len]
        else:
            rc_temp = Art.revcomp(template)
            read.seq_ref = rc_temp[0: self.read_len]

        read.bpos = 0
        read.ref2read()

        # constant quality scores
        read.quals = [qual_val] * read.read_len

        return read

    def next_pair_indel_seq(self, template):
        """
        Get a fwd/rev pair of reads for a template, where each read is sequenced off the ends.
        :param template: the target tempalte to sequencing fwd/rev
        :return: a dict {'fwd': SeqRead, 'rev': SeqRead}
        """
        return {'fwd': self.next_read_indel_seq(template, True),
                'rev': self.next_read_indel_seq(template, False)}

    def next_read_indel_seq(self, template, plus_strand):
        """
        Generate a read off a supplied target template sequence.
        :param template: the target template to sequence
        :param plus_strand: forward: True, reverse: False
        :return: SeqRead
        """
        read = self._new_read(plus_strand=plus_strand)
        mut_temp = Art.make_mutable(template)
        if len(mut_temp) < read.read_len:
            # for templates shorter than the requested length, we sequence its total extent
            read.read_len = len(mut_temp)

        # indels
        slen = read.get_indel()

        # ensure that this read will fit within the extent of the template
        if self.read_len - slen > len(mut_temp):
            slen = read.get_indel_2()

        if read.is_plus_strand:
            read.seq_ref = mut_temp[0: self.read_len - slen]
        else:
            rc_temp = Art.revcomp(template)
            read.seq_ref = list(rc_temp[0: self.read_len - slen])

        read.bpos = 0
        read.ref2read()

        # simulated quality scores from profiles
        read.quals = self.emp_dist.get_read_qual(read.length(), read.is_plus_strand)
        # the returned quality scores can spawn sequencing errors
        read.parse_error()

        return read

    def next_read_indel_at(self, pos, plus_strand):
        """
        Create a read with an already determined position and direction.
        :param pos: position for read
        :param plus_strand: True = forward, False = reverse
        :return: SeqRead
        """
        read = self._new_read()
        read.is_plus_strand = plus_strand

        # indels
        slen = read.get_indel()

        # ensure that this read will fit within the extent of the reference
        if pos + self.read_len - slen > len(self.ref_seq):
            slen = read.get_indel_2()

        if read.is_plus_strand:
            read.seq_ref = self.ref_seq[pos: pos + self.read_len - slen]
        else:
            read.seq_ref = self.ref_seq_cmp[pos: pos + self.read_len - slen]

        read.bpos = pos
        read.ref2read()

        # simulated quality scores from profiles
        read.quals = self.emp_dist.get_read_qual(read.length(), True)
        # the returned quality scores can spawn sequencing errors
        read.parse_error()

        return read

    def next_read_indel(self):
        """
        Create the next SeqRead and its accompanying quality scores. Position and direction are
        determined by uniform random seletion.

        :return: SeqRead
        """
        # random position anywhere in valid range
        pos = Art.random_int(0, self.valid_region)

        # is it the forward strand?
        plus_strand = Art.random_coin_toss()

        return self.next_read_indel_at(pos, plus_strand)

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
    import numpy as np

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
    parser.add_argument('outbase', help='Output base name')
    args = parser.parse_args()

    if args.xfold and args.num_reads:
        raise RuntimeError('xfold and num-reads are mutually exclusive options')

    with open('{0}.r1.fq'.format(args.outbase), 'w', buffering=262144) as r1_h, \
            open('{0}.r2.fq'.format(args.outbase), 'w', buffering=262144) as r2_h:

        for input_record in SeqIO.parse(args.fasta, 'fasta'):

            # ref to string
            ref_seq = str(input_record.seq)

            # empirical distribution from files
            emp_dist = EmpDist(args.profile1, args.profile2)

            # init Art
            art = Art(args.read_length, emp_dist,
                      args.ins_rate, args.del_rate, seed=args.seed)

            if args.xfold:
                num_seq = int(math.ceil(len(art.ref_seq) / args.read_length * args.xfold))
            else:
                num_seq = args.num_reads

            print 'Generating {0} reads for {1}'.format(num_seq, input_record.id)

            print_rate = int(num_seq*100./10)

            for n in xrange(num_seq):

                ins_len = None
                while True:
                    ins_len = int(np.ceil(np.random.normal(500, 50)))
                    if ins_len > 200:
                        break

                # pick a random position on the chromosome, but we're lazy and don't
                # handle the edge case of crossing the origin
                pos = Art.random_int(0, len(ref_seq)-ins_len)

                # get next read and quals
                pair = art.next_pair_indel_seq(list(ref_seq[pos : pos + ins_len]))

                # create file records
                SeqIO.write(pair['fwd'].read_record(SeqRead.read_id(input_record.id, n)), r1_h, 'fastq')
                SeqIO.write(pair['rev'].read_record(SeqRead.read_id(input_record.id, n)), r2_h, 'fastq')

                if ((n+1)*100) % print_rate == 0:
                    print '\tWrote {0} pairs'.format(n+1)
            break