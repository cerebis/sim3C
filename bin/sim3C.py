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
import os
import sys
import time
from collections import OrderedDict
from distutils.version import StrictVersion

import numpy as np
from Bio import Alphabet
from Bio import SeqIO
from Bio.Restriction import *

import Art
import abundance
import empirical_model
import io_utils

if StrictVersion(np.__version__) < StrictVersion("1.9.0"):
    sys.stderr.write("Error: numpy version 1.9.0 or later required\n")
    sys.stderr.write("If numpy is installed in both your home & system directory, "
                     "you may need to run with python -S\n")
    sys.exit(1)

"""

GLOBALS

"""

# Mixed geom/unif model
MIXED_GEOM_PROB = 6.0e-6
CID_GEOM_PROB = 8.0e-6

# set state for random number generation
RANDOM_STATE = None

# Average & SD size that fragments are sheared (or tagmented) to during adapter ligation
SHEARING_MEAN = 400
SHEARING_SD = 50


"""

Exception classes

"""


class Sim3CException(Exception):
    """Module base exception class"""
    pass


class NoCutSitesException(Sim3CException):
    """Occurs when a target template contains no cutsites for a specified restriction enzyme"""
    def __init__(self, seq_name, enz_name):
        self.message = 'sequence [{0}] had no cutsites for enzyme [{1}]'.format(seq_name, enz_name)


class EmptyRegistryException(Sim3CException):
    pass


class MonochromosomalException(Sim3CException):
    pass


def get_enzyme_instance(enzyme_name):
    """
    Using RestrictionBatch class, convert an enzyme name to
    a concrete restriction enzyme class instance.
    :param enzyme_name:
    :return: class instance (AbstractCut -- maybe)
    """
    rb = RestrictionBatch([enzyme_name])
    enz = RestrictionBatch.get(rb, enzyme_name, add=False)
    return enz


def find_restriction_sites(enzyme_name, seq):
    """
    For supplied enzyme, find all restriction sites in a given sequence
    returns list of sites. Note, positions are converted from 1-based to 0-based.

    :param enzyme_name: name of enzyme to use in digestion
    :param seq: sequence to digest
    :return: list of genomic coordinates
    """
    en = get_enzyme_instance(enzyme_name)
    return np.array(en.search(seq, linear=False)) - 1


def find_priming_sites(oligo, seq):
    """
    For supplied priming sequence, find positions of all matches in a given sequence
    returns list of sites.
    :param oligo: Bio.Seq primer sequence
    :param seq: Bio.Seq template to search
    :return: list of sites
    """
    array = []
    for m in re.finditer(oligo, str(seq)):
        array.append(m.end())
    rc_oligo = Seq(oligo)
    rc_oligo.reverse_complement()
    for m in re.finditer(str(rc_oligo), str(seq)):
        array.append(m.end())
    return array


def make_read(seq, fwd_read, length):
    """
    From sequence, make a forward or reverse read. If the read is longer than the total sequence
    return the entire sequence.
    :param seq: Bio.Seq template from which to generate a read
    :param fwd_read: True if this is on the forward strand
    :param length: read length
    :return: Bio.Seq object representing the read
    """

    # edge case - return whole sequence of read longer than sequence
    if length >= len(seq):
        read = seq
        if not fwd_read:
            read = read.reverse_complement(id=True, description=True)
        return read

    # return forward or reverse read
    if fwd_read:
        read = seq[:length]
    else:
        read = seq[-length:]
        read = read.reverse_complement(id=True, description=True)

    return read


def write_reads(handle, sequences, output_format, dummy_q=False):
    """Append sequence to file in fastq format. Currently no effort is spent
    to stream these records"""
    if dummy_q:
        for s in sequences:
            s.letter_annotations['phred_quality'] = [50] * len(s)

    SeqIO.write(sequences, handle, output_format)


class Part:
    """
    Represents one part of the two-piece 3C ligation fragment.
    """
    def __init__(self, seq, pos1, pos2, fwd, replicon):
        """
        :param seq: template sequence
        :param pos1: beginning
        :param pos2: end
        :param fwd: strand True = forward
        :param replicon: from which this part derives
        """
        self.seq = seq
        self.pos1 = pos1
        self.pos2 = pos2
        self.fwd = fwd
        self.replicon = replicon

    def __repr__(self):
        return repr((self.seq, self.pos1, self.pos2, self.fwd, self.replicon))

    def length(self):
        return self.pos2 - self.pos1 + 1


class Replicon:
    """Represents a replicon which holds a reference to its containing cell"""

    PART_DESC_FMT = '{0:d}:{1}:{2}'

    def __init__(self, name, parent_cell, sequence, cutters):
        """
        Instantiate a Replicon instance
        :param name: name of the replicon
        :param parent_cell: cell to which this replicon belongs
        :param sequence: DNA sequence of this replicon
        :param cutters: restriction enzymes used to digest this replicon
        """
        self.name = name
        self.parent_cell = parent_cell
        self.sequence = sequence

        # for each enzyme, pre-digest the replicon sequence
        self.cut_sites = {}
        for cname in cutters:
            self.cut_sites[cname] = np.array(find_restriction_sites(cname, sequence.seq))
            if len(self.cut_sites[cname]) <= 0:
                raise NoCutSitesException(name, cname)

        # initialise CID blocked empirical model
        self.cid_blocks = empirical_model.cids_to_blocks(
            empirical_model.generate_random_cids(RANDOM_STATE, self.length(),
                                                 chr_prob=BACKBONE_PROB,
                                                 chr_shape=MIXED_GEOM_PROB,
                                                 cid_shape=CID_GEOM_PROB))

    def __repr__(self):
        return repr((self.name, self.parent_cell, self.sequence))

    def __str__(self):
        return str(self.parent_cell) + '.' + self.name

    def length(self):
        """
        Length of the replicon's DNA sequence
        :return:
        """
        return len(self.sequence)

    def is_alone(self):
        """
        :return: True if this is the only replicon in a cell.
        """
        return self.parent_cell.number_replicons() == 1

    def subseq(self, start, length, rev=False):
        """
        Create a subsequence, where the replicon is always treated as circular.

        :param start: starting genomic position
        :param length: length of subsequence
        :param rev: reverse complement this sequence.
        :return: subseq Seq object
        """
        end = start + length
        diff = end - self.length()
        if diff > 0:
            # sequence will wrap around
            ss = self.sequence[start:] + self.sequence[:diff]
            ss.description = Replicon.PART_DESC_FMT.format(rev, start+1, diff)
        else:
            ss = self.sequence[start:end]
            ss.description = Replicon.PART_DESC_FMT.format(rev, start+1, end)

        if rev:
            ss.reverse_complement(id=True, description=True)

        return ss

    def random_cut_site(self, cutter_name):
        """
        Draw a cut-site at random, following a uniform distribution.

        :param cutter_name: enzyme name
        :return: genomic coordinates of cut-site
        """
        cs = self.cut_sites[cutter_name]
        return RANDOM_STATE.choice(cs)

    def nearest_cut_site_above(self, cutter_name, pos):
        """
        Return nearest cut-site for enzyme 'cutter_name' right of the supplied genomic position.

        :param cutter_name: enzyme name
        :param pos: genomic position
        :return:
        """
        cs = self.cut_sites[cutter_name]
        idx = np.searchsorted(cs, pos)
        if idx == len(cs):
            return cs[0]
        else:
            return cs[idx]

    def nearest_cut_site_below(self, cutter_name, pos):
        """
        Return nearest cut-site for enzyme 'cutter_name' left of the supplied genomic position.
        :param cutter_name: enzyme
        :param pos: genomic position
        :return:
        """
        cs = self.cut_sites[cutter_name]
        idx = np.searchsorted(cs, pos)
        if idx == 0:
            return cs[-1]
        else:
            return cs[idx - 1]

    def nearest_cut_site_by_distance(self, cutter_name, pos):
        """
        Find the nearest restriction cut site for the specified cutter type [4cut, 6cut]
        returns genomic position of nearest cut site. Assumes circular DNA for edge cases.

        :param cutter_name: name of restriction enzyme to lookup
        :param pos: position from which to find the nearest site for enzyme "cutter_name"
        :return: position of cutsite
        """

        cs = self.cut_sites[cutter_name]
        above = np.searchsorted(cs, pos)

        if above == len(cs):
            # pos beyond the last site, so wrap around
            prev = above - 1
            return cs[prev if pos - cs[prev] <= self.length() - pos + cs[0] else 0]
        elif above == 0:
            # pos before the first site, so wrap around
            return cs[above if cs[above] - pos <= self.length() - cs[-1] + pos else -1]
        else:
            prev = above - 1
            return cs[prev if pos - cs[prev] <= cs[above] - pos else above]


    @staticmethod
    def get_loc(emp_dist, origin):
        """
        Get a random location from 'origin', following the specified empirical distribution.
        :param emp_dist: EmpiricalDistribution to draw from
        :param origin: relative origin from which a location is derived.
        :return: random position from origin, +ve or -ve in direction
        """
        delta = int(emp_dist.rand())
        # either +ve/-ve relative to origin
        if RANDOM_STATE.uniform() < 0.5:
            return origin - delta
        else:
            return origin + delta

    def constrained_location_cids(self, origin):
        """
        Get a random location along this replicon's extent, relative to 'origin'.

        In our model, this is constrained by both the traditional experimentally determined HiC/3C
        relationship with separation and also our modelling of chromosomal interacting domains (CIDs).
        Whether a CID plays a part is dependent on whether the specified origin falls within a region
        defined as having a CID.

        :param origin: relative origin
        :return: another randomly drawn location.
        """
        block = self.cid_blocks[origin].pop()
        ovl_invs = block.data['inv_list']
        if len(ovl_invs) == 1:
            # if only the background distribution governs this block
            chosen_inv = ovl_invs[0]
        else:
            # pick a cid or background from those defined for this block
            rc = RANDOM_STATE.choice(len(ovl_invs), p=block.data['prob_list'])
            chosen_inv = ovl_invs[rc]

        loc = Replicon.get_loc(chosen_inv.data['empdist'], origin)
        loc -= chosen_inv.begin      # remove shift
        loc %= chosen_inv.length()   # modulo segment size
        loc += chosen_inv.begin      # replace shift
        return loc


class Cell:
    """Represents a cell in the community, that is a container of one or more replicons."""

    def __init__(self, name, abundance):
        """
        :param name: a unique name for this cell
        :param abundance: relative abundance of this cell type within the community
        """
        self.name = name
        self.abundance = float(abundance)
        self.replicon_registry = {}
        self.index_to_name = {}
        self.cdf = None
        self.pdf = None

    def __repr__(self):
        return repr((self.name, self.abundance))

    def __str__(self):
        return self.name

    def init_prob(self):
        """
        Initialize the probabilities for replicon selection from within a cell and
        afterwards, produce a CDF which will be used for random sampling.

        Note: this method should be called after all replicons for a cell have been
        registered.
        """
        # Number of replicons in cell
        n_rep = len(self.replicon_registry)
        if n_rep <= 0:
            raise EmptyRegistryException('cell contained no useful replicons')

        # Uniform probability initially
        prob = np.array([1.0 / n_rep] * n_rep)

        # Scale prob by replicon length and initialise index to name table.
        for n, repA in enumerate(self.replicon_registry.values()):
            self.index_to_name[n] = repA.name
            prob[n] = prob[n] * repA.length()

        # Normalised, this becomes the piece-wise PDF
        self.pdf = prob / prob.sum()

        # Initialize the piece-wise cdf from the pdf, but include leading 0 bin.
        self.cdf = np.hstack((-np.inf, np.cumsum(self.pdf)))

    def register_replicon(self, replicon):
        """
        Add a Replicon to the cell's registry
        :param replicon: Replicon to add
        """
        self.replicon_registry[replicon.name] = replicon

    def number_replicons(self):
        """
        Count the number of replicons known to this cell.
        :return: integer count of replicons
        """
        return len(self.replicon_registry)

    def _select_replicon_index(self, x):
        """
        From a cell, return the index of a replicon by sampling CDF at given value x
        :param x: domain of CDF
        :return: index of replicon
        """
        return np.searchsorted(self.cdf, x, side='right') - 1

    def pick_inter_rep(self, skip_this):
        if self.number_replicons() == 1:
            raise MonochromosomalException('cannot pick another in single replicon cell')

        # Keep trying until we pick a different rep.
        rep = skip_this
        while rep is skip_this:
            ri = self._select_replicon_index(RANDOM_STATE.uniform())
            rep = self.replicon_registry.get(self.index_to_name[ri])

        return rep


class Community:
    """Represents the community, maintaining a registry of cells and replicons.

    Both cells and replicons are expected to be uniquely named across the commmunity. This constraint comes
    from the requirement that a multi-fasta file contain unique sequence identifiers.

    A table is used for a basic community definition, with the following three column format. Lines beginning
    with # and empty lines are ignored.

    [replicon name] [cell name] [abundance]
    """
    def __init__(self, inter_rep_prob, spurious_prob, profile, seq_filename, cutters):
        """
        Instantiate a Community.
        :param inter_rep_prob: probability of inter replicon ligation
        :param spurious_prob: probability of a spurious ligation (noise)
        :param profile: abundance profile of the community members
        :param seq_filename: multi-fasta containing the genomic sequences of all community members
        :param cutters: list of restriction enzymed used in digestion.
        """
        self.pdf = None
        self.cdf = None
        self.totalRawAbundance = 0
        self.replicon_registry = OrderedDict()
        self.index_to_name = None
        self.cell_registry = OrderedDict()
        self.inter_rep_prob = inter_rep_prob
        self.spurious_prob = spurious_prob
        self.cutters = cutters

        # Read in the sequences
        sequences = SeqIO.to_dict(SeqIO.parse(open(seq_filename), 'fasta', Alphabet.generic_dna))

        # build the registries from the defined community profile
        for abn in profile.values():
            try:
                parent_cell = self.register_cell(abn.cell, abn.val)
                self.build_register_replicon(abn.name, parent_cell, sequences.get(abn.name))
            except NoCutSitesException as e:
                print 'Warning: {0}'.format(e.message)

        # init community wide probs
        self._init_replicon_prob()

        # init individual cell probs
        cell_names = self.cell_registry.keys()
        for cn in cell_names:
            try:
                self.cell_registry[cn].init_prob()
            except EmptyRegistryException as e:
                print 'Removing [{0}] as {1}'.format(cn, e)
                del self.cell_registry[cn]
                # refresh abundance and distributions
                self._update_total_abundance()
                self._init_replicon_prob()

        if len(self.cell_registry) <= 0:
            raise EmptyRegistryException('Community is empty. No genomes were registered')

    def build_register_replicon(self, name, parent_cell, sequence):
        """Add a new replicon to a cell type in community"""
        if name in self.replicon_registry:
            raise Sim3CException('Replicon name [{0}] is not unique within the community'.format(name))

        replicon = Replicon(name, parent_cell, sequence, self.cutters)
        self.replicon_registry[name] = replicon
        parent_cell.register_replicon(replicon)
        return replicon

    def register_cell(self, cell_name, cell_abundance):
        """Add a new replicon to the community replicon registry"""
        parent_cell = self.cell_registry.get(cell_name)
        if parent_cell is None:
            parent_cell = Cell(cell_name, cell_abundance)
            self.cell_registry[cell_name] = parent_cell
            self._update_total_abundance()
        return parent_cell

    # Total abundance of all cells in registry.
    def _update_total_abundance(self):
        """Recalculate the total relative abundance specified for the community by referring to the registry"""
        ab = 0
        for ca in self.cell_registry.values():
            ab += ca.abundance
        self.totalRawAbundance = ab

    def _init_replicon_prob(self):
        """Initialize the probabilities for replicon selection, given the abundances, etc.
        Normalization is always applied. Afterwards, produce a CDF which will be used for
        random sampling.
        """
        self.index_to_name = {}
        prob = np.zeros(len(self.replicon_registry))
        i = 0
        for repA in self.replicon_registry.values():
            self.index_to_name[i] = repA.name
            prob[i] = repA.parent_cell.abundance / self.totalRawAbundance * repA.length()
            i += 1

        tot_prob = prob.sum()

        if len(self.replicon_registry) <= 0 or tot_prob == 0:
            raise EmptyRegistryException('No replicons passed preprocessing or the set had zero probability density')

        # normalized piece-wise pdf
        self.pdf = prob / tot_prob

        # initialize the piece-wise cdf from the pdf, adding a 0 element to the beginning
        self.cdf = np.hstack((-np.inf, np.cumsum(self.pdf)))

    def _select_replicon_index(self, x):
        """
        From the entire community, return the index of a replicon by sampling CDF at given value x
        :param x: domain of CDF [0..1]
        :return: nearest corresponding y - (replicon index)
        """
        return np.searchsorted(self.cdf, x, side='right') - 1

    def pick_replicon(self, skip_index=None):
        """Random selection of replicon from community. If skipIndex supplied, do not return
        the replicon with this index.

        return the index"""
        if skip_index is None:
            return self._select_replicon_index(RANDOM_STATE.uniform())
        else:
            ri = skip_index
            while ri == skip_index:
                ri = self._select_replicon_index(RANDOM_STATE.uniform())
            return ri

    def is_spurious_event(self):
        """
        Spurious events are those which link unrelated fragments together during a 3C ligation
        protocol. They are effectively noise within the system.
        :return: True - the next fragment should be constructed as spurious
        """
        return RANDOM_STATE.uniform() < self.spurious_prob

    def is_intra_rep_event(self):
        """
        Choose if the mate is intra or inter replicon associated. This is a simple
        binary paritioning with a chosen threshold frequency.
        """
        return RANDOM_STATE.uniform() > self.inter_rep_prob

    def get_replicon_by_index(self, index):
        return self.replicon_registry.get(self.index_to_name[index])

    def get_replicon_by_name(self, name):
        return self.replicon_registry.get(name)

    def unconstrained_read_location(self, replicon_index):
        """Return a location (position and strand) on a replicon where the position is
        uniformly sampled across the replicons sequence.

        Returns tuple (pos=int, strand=bool)
        """
        replicon = self.get_replicon_by_index(replicon_index)
        return int(RANDOM_STATE.uniform() * replicon.length()), True


"""

Fragment creation methods

"""

def get_partial_fragment_length():
    """
    Return lengths for partial fragments which are half the length of the desired ligation product.
    :return: length for this portion of ligation fragment
    """
    return int(0.5 * RANDOM_STATE.normal(SHEARING_MEAN, SHEARING_SD))


def make_unconstrained_part_a():
    """
    Choose a cut site at random across a randomly selected replicon. Further,
    model a normally distributed DNA fragment length.
    :return: Part
    """

    repl = comm.get_replicon_by_index(comm.pick_replicon())
    pos = repl.random_cut_site(CUTTER_NAME)
    frag_len = get_partial_fragment_length()
    seq = repl.subseq(pos, frag_len)
    return Part(seq, pos, pos + frag_len, True, repl)


def make_any_part_b():
    """
    This selection is at random across the community and not reflective of partA or physical
    organisation. Effectively this is the same as our first selection method.
    :return: Part
    """
    return make_unconstrained_part_a()


def make_unconstrained_part_b(first_part):
    """
    Choose a inter-replicon cut-site, meaning not the same replicon as was first picked but
    from the same cell. The genomic coordinate on this new part is thus unconstrained.
    :param first_part: part A, the already chosen piece of the ligation fragment from a particular replicon.
    :return: another part B, which is not from the same replicon as part A.
    """
    diff_repl = first_part.replicon.parent_cell.pick_inter_rep(first_part.replicon)
    pos = diff_repl.random_cut_site(CUTTER_NAME)
    frag_len = get_partial_fragment_length()
    seq = diff_repl.subseq(pos, frag_len)
    return Part(seq, pos, pos + frag_len, True, diff_repl)


def make_constrained_part_b(first_part):
    """
    Choose a coordinate on the same replicon as the first part, part A. As such, the coordinates
    follow a particular distribution dependent on genomic separation.
    :param first_part: part A, the already chosen piece of the ligation fragment from a particular replicon.
    :return: another part B, on the same replicon which follows the distribution of separation.
    """
    loc = first_part.replicon.constrained_location_cids(first_part.pos1)
    if RANDOM_STATE.uniform() < ANTIDIAG_RATE:
        # an anti-diagonal event
        loc = first_part.replicon.length() - loc

    pos = first_part.replicon.nearest_cut_site_by_distance(CUTTER_NAME, loc)
    frag_len = get_partial_fragment_length()
    seq = first_part.replicon.subseq(pos, frag_len)
    return Part(seq, pos, pos+frag_len, True, first_part.replicon)



if __name__ == '__main__':

    #
    # Commandline interface
    #
    parser = argparse.ArgumentParser(description='Simulate HiC read pairs')

    parser.add_argument('-C', '--compress', choices=['gzip', 'bzip2'], default=None,
                        help='Compress output files')
    parser.add_argument('--sample-name', default='SIM3C', help='Sample name used in reads [SIM3C]')
    parser.add_argument('-f', '--ofmt', dest='output_format', default='fastq', choices=['fasta', 'fastq'],
                        help='Output format')

    parser.add_argument('-r', '--seed', metavar='INT', type=int, default=int(time.time()),
                        help="Random seed for initialising number generator")
    parser.add_argument('--site-dup', default=False, action='store_true',
                        help='HiC style ligation junction site duplication')
    parser.add_argument('--enzyme', default='NlaIII', help='Restriction enzyme [NlaIII]')

    parser.add_argument('-n', '--num-frag', metavar='INT', type=int, required=True,
                        help='Number of Hi-C fragments to generate reads')
    parser.add_argument('-l', '--read-length', metavar='INT', type=int, required=True,
                        help='Length of reads from Hi-C fragments')
    parser.add_argument('--frag-mean', metavar='INT', type=int, default=SHEARING_MEAN,
                        help='Mean fragment size [400]')
    parser.add_argument('--frag-sd', metavar='INT', type=int, default=SHEARING_SD,
                        help='Standard deviation of fragment size [50]')
    parser.add_argument('--min-frag', metavar='INT', type=int, default=200,
                        help='Minimum fragment length [200]')
    parser.add_argument('--max-frag', metavar='INT', type=int, default=1000,
                        help='Maximum fragment length [1000]')

    parser.add_argument('--anti-rate', metavar='FLOAT', type=float, default=0.2,
                        help='Rate of anti-diagonal fragments [0.2]')
    parser.add_argument('--backbone-prob', metavar='FLOAT', type=float, default=0.333,
                        help='Probability of regular chromosomal interaction, rather than CID [0.333]')
    parser.add_argument('--inter-prob', dest='inter_prob', metavar='FLOAT', type=float, default=0.9,
                        help='Probability that a fragment spans two replicons within a single genome [0.9]')
    parser.add_argument('--spurious-prob', dest='spur_prob', metavar='FLOAT', type=float, default=0.01,
                        help='Probability that a spurious fragment is formed [0.01]')

    parser.add_argument('-P', '--profile', dest='profile', metavar='FILE',
                        help='Community abundance profile')
    parser.add_argument('--profile-out', metavar='FILE', default='profile.tsv',
                        help='Output file when generating community profile', required=False)

    parser.add_argument('--dist', metavar='DISTNAME', choices=['equal', 'uniform', 'lognormal'],
                        help='Abundance profile distribution choices: equal, uniform, lognormal')
    parser.add_argument('--lognorm-mu', metavar='FLOAT', type=float, default='1', required=False,
                        help='Log-normal relative abundance mu parameter')
    parser.add_argument('--lognorm-sigma', metavar='FLOAT', type=float, default='1', required=False,
                        help='Log-normal relative abundance sigma parameter')

    parser.add_argument('--no-read-errors', default=False, action='store_true', help='Simulate sequencing errors')
    parser.add_argument('--read-profile', help='The name of an ART sequencing machine profile [EmpMiSeq250]',
                        default='EmpMiSeq250', choices=Art.ILLUMINA_PROFILES.keys())
    parser.add_argument('--ins-rate', type=float, default=9.e-5, help='Insert rate [9e-5]')
    parser.add_argument('--del-rate', type=float, default=1.1e-4, help='Deletion rate [1.1e-4]')

    parser.add_argument(dest='genome_seq', metavar='FASTA',
                        help='Genome sequences for the community')
    parser.add_argument(dest='output_file', metavar='OUTPUT',
                        help='Output Hi-C reads file')
    args = parser.parse_args()

    if 'community_table' in args and args.dist:
        raise RuntimeError('Cannot define abundance both explicitly as a table (-t) and a distribution (--dist).')

    #
    # Main routine
    #
    CUTTER_NAME = args.enzyme
    RANDOM_STATE = np.random.RandomState(args.seed)
    SHEARING_MEAN = args.frag_mean
    SHEARING_SD = args.frag_sd
    ANTIDIAG_RATE = args.anti_rate
    BACKBONE_PROB = args.backbone_prob

    SEQ_ID_FMT = args.sample_name + ':{seed}:{origin}:1:1:1:{idx}'

    #
    # Prepare community abundance profile, either procedurally or from a file
    #
    #   Note: currently, all sequences for a single taxon are
    #   treated equally.
    #
    profile = None
    if args.dist:
        # procedural, number of taxa defined by number of sequences (assumes contiguous genomes)
        seq_index = None
        try:
            seq_index = SeqIO.index(args.genome_seq, 'fasta')
            profile = abundance.generate_profile(RANDOM_STATE, list(seq_index), mode=args.dist,
                                                 lognorm_mu=args.lognorm_mu, lognorm_sigma=args.lognorm_sigma)

            # present result to console
            profile.write_table(sys.stdout)

            # save result to file
            with open(os.path.join(os.path.dirname(args.output_file), args.profile_out), 'w') as prf_h:
                profile.write_table(prf_h)

        finally:
            if seq_index:
                seq_index.close()
    else:
        # read pre-existing abundance profile from file
        with open(args.profile, 'r') as h_table:
            profile = abundance.read_profile(h_table)
            profile.normalize()

    # Initialize community object
    print "Initializing community"

    try:
        comm = Community(args.inter_prob, args.spur_prob, profile, args.genome_seq, [CUTTER_NAME])

        # Junction produced in Hi-C prep
        cut_site = get_enzyme_instance(CUTTER_NAME).site
        hic_junction = cut_site + cut_site

        frag_lengths = {'a': [], 'b': []}

        # Open output file for writing reads
        with io_utils.open_output(args.output_file, mode='w', compress=args.compress) as h_output:

            print "Creating reads"
            skip_count = 0
            overlap_count = 0
            frag_count = 0

            # initialise ART
            art = Art.Art(args.read_length,
                          Art.EmpDist.create(args.read_profile),
                          args.ins_rate, args.del_rate, seed=args.seed)

            # set the method used to generate reads
            next_pair = art.next_pair_simple_seq if args.no_read_errors else art.next_pair_indel_seq

            while frag_count < args.num_frag:
                # Fragment creation

                # Create PartA
                # Steps
                # 1) pick a replicon at random
                # 2) pick a cut-site at random on replicon
                # 3) flip a coin for strand
                part_a = make_unconstrained_part_a()

                # Create PartB
                # Steps
                # 1) choose if intra or inter replicon
                # 2) if intER create partB as above
                # 3) if intRA select from geometric

                if comm.is_spurious_event():
                    # fusion of two randomly chosen fragments -- not reflecting physical organisation
                    part_b = make_any_part_b()

                elif part_a.replicon.is_alone() or comm.is_intra_rep_event():
                    # ligation is between two fragments on same replicon
                    part_b = make_constrained_part_b(part_a)

                else:
                    # ligation crosses replicons
                    part_b = make_unconstrained_part_b(part_a)

                frag_lengths['a'].append(part_a.length())
                frag_lengths['b'].append(part_b.length())

                # Join parts A and B
                if args.site_dup:
                    fragment = part_a.seq + hic_junction + part_b.seq
                else:
                    # meta3C does not create duplicated sites
                    fragment = part_a.seq + part_b.seq

                if len(fragment) < args.min_frag or len(fragment) > args.max_frag:
                    # only accept fragments within a size range
                    skip_count += 1
                    continue

                if part_b.pos1 < part_a.pos2 and part_a.pos2 < part_b.pos2:
                    overlap_count += 1
                    continue

                if part_a.pos1 < part_b.pos2 and part_b.pos2 < part_a.pos2:
                    overlap_count += 1
                    continue

                # create sequencing read pair for fragment
                pair = next_pair(str(fragment.seq))

                # emulate Illumina recent naming convention, where we re-purpose fields as follows:
                #
                # SIM3C:{RANDOM SEED}:{ORIGIN SEQ}:1:1:1:{INDEX} {R1/2}:N:1 {RECOMP}:{PART_START}:{PART_STOP}
                #
                read1 = pair['fwd'].read_record(SEQ_ID_FMT.format(seed=args.seed, origin=part_a.seq.id, idx=frag_count),
                                                desc='1:N:18:1 {0}'.format(part_a.seq.description))

                read2 = pair['rev'].read_record(SEQ_ID_FMT.format(seed=args.seed, origin=part_b.seq.id, idx=frag_count),
                                                desc='2:N:18:1 {0}'.format(part_b.seq.description))

                # write to interleaved file
                write_reads(h_output, [read1, read2], args.output_format, dummy_q=False)

                frag_count += 1

        print "Ignored " + str(skip_count) + " fragments due to length restrictions"
        print "Ignored " + str(overlap_count) + " fragments due to overlap"

        print 'Fragment stats:'
        d = np.array(frag_lengths.values())
        print 'Parts:   means: {0} std: {1}'.format(d.mean(1), d.std(1))
        print 'Overall: means: {0} std: {1}'.format(d.sum(0).mean(), d.sum(0).std())
        rt_fwd = np.count_nonzero(d[0, :] < args.read_length)
        rt_rev = np.count_nonzero(d[1, :] < args.read_length)
        if rt_fwd > 0 or rt_rev > 0:
            print 'Warning!! Read-through occurred on: {0} forward and {1} reverse reads'.format(rt_fwd, rt_rev)

    except Sim3CException as e:
        import traceback
        sys.stderr.write(traceback.format_exception_only(type(e), e)[0])
        sys.exit(1)
