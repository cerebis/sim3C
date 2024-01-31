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
import logging
import tqdm

from Bio import SeqIO
from collections import namedtuple

from .art import Art, EmpDist, ambiguous_base_filter, validator
from .community import Community
from .exceptions import *
from .site_analysis import get_enzyme_instance
from .random import np_normal
import sim3C.random as random


logger = logging.getLogger(__name__)


class ReadGenerator(object):
    """
    Generate inserts and subsequent read-pairs for a particular library preparation method. Primarily,
    3C does not generate a duplication of the cut-site, whereas HiC's enriching for ligation products by
    use of biotinylation does produce site duplication during infill of overhangs.

    The desired insert size and its variability are specified here.

    Other sequencing read simulation parameters are supplied here initialise ART.
    """
    def __init__(self, method, enzyme, prefix='SIM3C', simple=False, machine_profile='EmpMiSeq250',
                 read_length=250, ins_rate=9.e-5, del_rate=1.1e-4,
                 insert_mean=500, insert_sd=100, insert_min=100, insert_max=None):
        """
        Initialise a read generator.
        :param method: The two library preparation methods are: 'meta3c' or 'hic'.
        :param enzyme: The employed restriction enzyme
        :param prefix: leading string for read names
        :param simple: True: do not simulate sequencing errors (faster), False: fully simulation sequencing
        :param machine_profile: ART Illumina error profile for a particular machine type. Default EmpMiSeq250
        :param read_length: desired read-length
        :param ins_rate: insertion rate (Art default: 9e-5)
        :param del_rate: deletion rate (Art default: 1.1e-4)
        :param insert_mean: mean sequencing insert length (must be > 0 and minimum)
        :param insert_sd: standard deviation of insert length  (must be < mean)
        :param insert_min: minimum allowable insert length (must be > 50)
        :param insert_max: maximum allowable insert length (must be > mean)
        """

        self.method = method
        self.prefix = prefix
        self.seq_id_fmt = prefix + ':{mode}:1:1:1:{idx} {r1r2}:Y:18:1'
        self.wgs_desc_fmt = 'WGS {segment.name}:{pos1}..{pos2}:{dir}'
        self._3c_desc_fmt = method.upper() + ' {segment1.name}:{pos1} {segment2.name}:{pos2}'

        assert insert_min > 50, 'Minimum allowable insert size is 50bp'
        assert insert_min < insert_mean, 'Minimum insert size must be less than expected mean'
        assert insert_mean > 0 and insert_sd > 0, 'Insert mean and stddev must be greater than 0'
        if insert_mean < insert_sd:
            logger.warning('Specified insert mean ({}) is less than stddev ({})'.format(insert_mean, insert_sd))
        if insert_mean - insert_sd < insert_min:
            logger.warning('Specified insert mean ({}) and stddev ({}) will produce many inserts below '
                           'the minimum allowable insert length ({})'.format(insert_mean, insert_sd, insert_min))

        if enzyme:
            self.cut_site = enzyme.ovhgseq * 2

        self.insert_mean = insert_mean
        self.insert_sd = insert_sd
        self.insert_min = insert_min
        self.insert_max = insert_max
        self.too_short = 0
        self.too_long = 0

        # initialise ART read simulator
        self.art = Art(read_length, EmpDist.create(machine_profile), ins_rate, del_rate)

        # set the method used to generate reads
        if simple:
            self.next_pair = self.art.next_pair_simple_seq
        else:
            self.next_pair = self.art.next_pair_indel_seq

        try:
            method_switcher = {
                'hic': self._part_joiner_sitedup,
                'meta3c': self._part_joiner_simple,
                'dnase': self._part_joiner_simple
            }
            self._part_joiner = method_switcher[self.method.lower()]
        except Exception:
            raise Sim3CException('unknown library preparation method ({}) Either: \'3c\' or \'hic\']'.format(method))

    def get_report(self):
        """
        Prepare a small report on the number of inserts which were constrained by either the maximum
        or minimum length settings.
        """
        msg = 'Constrained inserts: '
        msg += '[len < {min}] = {min_count}'.format(min=self.insert_min, min_count=self.too_short)
        if self.insert_max:
            msg += ', [len > {max}] = {max_count}'.format(max=self.insert_max, max_count=self.too_long)
        else:
            msg += ', [no max limit]'
        return msg

    @staticmethod
    def _part_joiner_simple(a, b):
        """
        Join two fragments end to end without any site duplication. The new
        fragment will begin at a_0 and end at b_max. Used when no end fills are
        applied to the library prep (Eg. meta3C)
        Altogether [a_0..a_max] + [b_0..b_max]
        :param a: fragment a
        :param b: fragment b
        :return: a + b
        """
        return a + b

    def _part_joiner_sitedup(self, a, b):
        """
        Join two fragments end to end where the cut-site is duplicated. The new
        fragment will begin at a_0 end at b_max. Used when end-fill is applied
        before further steps in library preparation (Eg. HiC)
        Altogether [a_0..a_max] + [cs_0..cs_max] + [b_0..b_max]
        :param a: fragment a
        :param b: fragment b
        :return: a + b
        """
        return a + self.cut_site + b

    def draw_insert(self):
        """
        Draw an insert length, midpoint and direction. Since this is normally distributed
        even with reasonable mean and stddev, it is possible to return lengths
        less than zero. These cases are handled by redrawing until we get a length
        greater than the specified minimum.
        length.
        :return: length, midpoint and direction tuple. Eg. (100, 56, True)
        """
        length = int(np_normal(self.insert_mean, self.insert_sd))
        if length < self.insert_min:
            self.too_short += 1
            length = self.insert_min
        elif self.insert_max and length > self.insert_max:
            self.too_long += 1
            length = self.insert_max
        midpoint = random.pcg_random.integer(length)
        is_fwd = random.pcg_random.integer(2) == 0
        return length, midpoint, is_fwd

    def make_wgs_readpair(self, segment, pos1, ins_len, is_fwd):
        """
        Create a fwd/rev read-pair simulating WGS sequencing.
        :param segment: segment from which to extract this read-pair
        :param pos1: the initial coordinate along the replicon.
        :param ins_len: the insert length
        :param is_fwd: will the insert be off the fwd strand
        :return: a read-pair dict
        """
        frag = segment.subseq(pos1, ins_len, is_fwd)
        pair = self.next_pair(bytes(frag.seq))
        pair['mode'] = 'WGS'
        pair['desc'] = self.wgs_desc_fmt.format(segment=segment, pos1=pos1, pos2=pos1 + ins_len, dir='F' if is_fwd else 'R')
        return pair

    def make_ligation_readpair(self, segment1, pos1, segment2, pos2, ins_len, ins_junc):
        """
        Create a fwd/rev read-pair simulating a ligation product (Eg. the outcome of
        HiC or meta3C library prep). As repl1 and repl2 can be the same, these ligation
        products can be inter-rep, intra-rep or spurious products.
        :param segment1: the first segment
        :param pos1:  the location along repl1
        :param segment2: the second segment
        :param pos2: the location along repl2
        :param ins_len: insert length
        :param ins_junc: junction point on insert
        :return: a read-pair dict
        """
        part_a = segment1.subseq(pos1 - ins_junc, ins_junc)
        part_b = segment2.subseq(pos2, ins_len - ins_junc)

        pair = self.next_pair(bytes(self._part_joiner(part_a, part_b).seq))
        pair['mode'] = '3C'
        pair['desc'] = self._3c_desc_fmt.format(segment1=segment1, pos1=pos1, segment2=segment2, pos2=pos2)
        return pair

    def write_readpair_biopython(self, h_out, pair, index, fmt='fastq'):
        """
        Write a read-pair object to a stream.
        :param h_out: the output stream
        :param pair: the read-pair to write
        :param index: a unique identifier for the read-pair. (Eg. an integer)
        :param fmt: the output file format
        """

        # create Bio.Seq objects for read1 (fwd) and read2 (rev)
        read1 = pair['fwd'].read_record(self.seq_id_fmt.format(mode=pair['mode'], idx=index, r1r2=1),
                                        desc=pair['desc'])
        read2 = pair['rev'].read_record(self.seq_id_fmt.format(mode=pair['mode'], idx=index, r1r2=2),
                                        desc=pair['desc'])
        # write to interleaved file
        SeqIO.write([read1, read2], h_out, fmt)

    def write_readpair_dnaio(self, writer, pair, index):
        """
        Write a read-pair object to a stream.
        :param writer: the output writer
        :param pair: the read-pair to write
        :param index: a unique identifier for the read-pair. (Eg. an integer)
        """

        # create Bio.Seq objects for read1 (fwd) and read2 (rev)
        read1 = pair['fwd'].read_record_dnaio(self.seq_id_fmt.format(mode=pair['mode'], idx=index, r1r2=1),
                                              desc=pair['desc'])
        read2 = pair['rev'].read_record_dnaio(self.seq_id_fmt.format(mode=pair['mode'], idx=index, r1r2=2),
                                              desc=pair['desc'])
        # write to interleaved file
        writer.write(read1, read2)


class SequencingStrategy(object):
    """
    A SequencingStrategy represents the whole experiment. This includes the reference data from which
    WGS and ligation products are generated, and the simulation of Illumina sequencing reads.

    Experiments are reproducible by supplying the same seed value.
    """

    Strategy = namedtuple('Strategy', 'method run')

    def __init__(self, profile_filename, seq_filename, enzyme_name, number_pairs,
                 method, read_length, prefix, machine_profile,
                 insert_mean=400, insert_sd=50, insert_min=50, insert_max=None,
                 anti_rate=0.25, spurious_rate=0.02, trans_rate=0.1,
                 efficiency=0.02,
                 ins_rate=9.e-5, del_rate=1.1e-4,
                 create_cids=True, simple_reads=True, linear=False, convert_symbols=False,
                 profile_format='table'):
        """
        Initialise a SequencingStrategy.

        :param profile_filename: the abundance profile for the community
        :param seq_filename: the matching sequence of replicon sequences in Fasta format
        :param enzyme_name: the restriction enzyme name (case sensitive)
        :param number_pairs: the number of read-pairs to generate
        :param method: the library preparation method (Either: 3c or hic)
        :param read_length: the length of reads
        :param prefix: read-names begin with this string
        :param machine_profile: ART Illumina sequencing machine profile
        :param insert_mean: mean insert length
        :param insert_sd: stddev insert length
        :param insert_min: minimum allowable insert length (must be > 50)
        :param insert_max: maximum allowable insert length (must be > mean)
        :param anti_rate: rate of anti-diagonal interactions
        :param spurious_rate: rate of spurious ligation products
        :param trans_rate: rate of inter-replicon (trans) ligation product
        :param efficiency: for meta3c simulation, the efficiency of ligation product generation.
        :param ins_rate: rate of sequencing insert errors
        :param del_rate: rate of sequencing deletion errors
        :param create_cids: simulate 3D structure, chromosomal interacting domains (CID)
        :param simple_reads: True: sequencing reads do not simulate error (faster), False: full simulation of sequencing
        :param linear: treat replicons as linear
        :param convert_symbols: if true, unsupported (by Art) symbols in the input sequences are converted to N
        :param profile_format: the format of the abundance profile (Either: table or toml)
        """
        self.profile_filename = profile_filename
        self.profile_format = profile_format
        self.seq_filename = seq_filename
        self.enzyme_name = enzyme_name
        self.number_pairs = number_pairs
        self.simple_reads = simple_reads
        self.method = method
        self.read_length = read_length
        self.insert_min = insert_min
        self.insert_max = insert_max
        self.efficiency = efficiency

        self.enzyme = None if not enzyme_name else get_enzyme_instance(enzyme_name)

        # reference sequences will be accessed via an SeqIO index. Optionally
        # overriding getter for validation and base filtering
        try:
            seq_index = SeqIO.index(seq_filename, 'fasta')
            if convert_symbols:
                # remove IUPAC ambiguity symbols
                seq_index = ambiguous_base_filter(seq_index)
            elif not simple_reads:
                # getter will now validate but not change sequence symbols
                seq_index = validator(seq_index)
        except (ValueError, TypeError):
            raise FastaException(seq_filename)

        # initialise the community for the reference data
        if profile_format == 'table':
            self.community = Community.from_profile(seq_index, profile_filename, self.enzyme,
                                                    anti_rate=anti_rate, spurious_rate=spurious_rate,
                                                    trans_rate=trans_rate, linear=linear)
        elif profile_format == 'toml':
            self.community = Community.from_toml(seq_index, profile_filename, self.enzyme,
                                                 anti_rate=anti_rate, spurious_rate=spurious_rate,
                                                 trans_rate=trans_rate, linear=linear)
        else:
            raise Sim3CException(f'unknown profile format ({profile_format}) Either: "table" or "toml"]')

        # preparate the read simulator for output
        self.read_generator = ReadGenerator(method, self.enzyme,
                                            prefix=prefix, simple=simple_reads, machine_profile=machine_profile,
                                            read_length=read_length, insert_mean=insert_mean,
                                            insert_sd=insert_sd, insert_min=insert_min, insert_max=insert_max,
                                            del_rate=del_rate, ins_rate=ins_rate)

        # the method determines the strategy governing the creation of
        # ligation products and WGS reads.
        try:
            strategy_switcher = {
                'hic': self._simulate_hic,
                'meta3c': self._simulate_meta3c,
                'dnase': self._simulate_dnase
            }
            self._selected_strat = self.Strategy(method, strategy_switcher[method.lower()])
        except Exception:
            raise Sim3CException('unknown library preparation method ({}) Either: \'3c\' or \'hic\']'.format(method))

    def run(self, ostream):
        """
        Add some pre and post detail to the selected strategy.
        :param ostream: the output stream for reads
        """
        logger.info('Starting sequencing simulation')
        logger.info('Library method: {}'.format(self._selected_strat.method))
        info = self._selected_strat.run(ostream)
        logger.info('Finished simulation')
        logger.info('Read counts: WGS reads = {wgs_count}, ligation products = {lig_count}'.format(**info))
        logger.info(self.read_generator.get_report())

    def _simulate_meta3c(self, ostream):
        """
        A strategy to simulate the sequencing of a 3C (meta3C) library.

        The most significant differentiator between 3C and HiC is that no biotin pulldown
        is employed for 3C sequencing experiments, thus only a small fraction of reads
        comprise ligation products, with the majority being conventional WGS reads.
        :param ostream: the output stream for reads
        """

        comm = self.community
        read_gen = self.read_generator
        efficiency = self.efficiency

        n_wgs = 0
        n_3c = 0

        for n in tqdm.tqdm(range(1, self.number_pairs+1)):

            # pick an replicon, position and insert size
            seg1, pos1 = comm.draw_any_by_extent()
            ins_len, midpoint, is_fwd = read_gen.draw_insert()

            if random.pcg_random.uniform() < efficiency and seg1.covers_site(pos1, midpoint):

                n_3c += 1

                # move pos1 to the nearest actual site
                pos1 = seg1.sites.find_nn(pos1)

                # is it spurious ligation
                if comm.is_spurious():

                    seg2, pos2 = comm.draw_any_by_site()

                # is it an inter-replicon (trans) ligation
                elif seg1.parent_repl.parent_cell.is_trans():

                    seg2 = seg1.parent_repl.parent_cell.draw_other_segment_by_sites(seg1.parent_repl)
                    pos2 = seg2.draw_any_site()

                # otherwise an intra-replicon (cis) ligation
                else:

                    seg2 = seg1.parent_repl.draw_any_segment_by_sites()
                    if seg2 == seg1:
                        # find the first intervening site between
                        # random constrained position pos2 and site pos1.
                        pos2 = seg1.draw_constrained_site(pos1)
                        pos2 = seg1.sites.find_first(pos2, pos1)
                    else:
                        # unconstrained, any site on seg2
                        pos2 = seg2.draw_any_site()

                # randomly permute source/destination
                if random.pcg_random.integer(2) == 0:
                    pos1, pos2 = pos2, pos1
                    seg1, seg2 = seg2, seg1

                pair = read_gen.make_ligation_readpair(seg1, pos1, seg2, pos2, ins_len, midpoint)

            # otherwise WGS
            else:

                n_wgs += 1
                # take the already drawn coordinates
                pair = read_gen.make_wgs_readpair(seg1, pos1, ins_len, is_fwd)

            read_gen.write_readpair_dnaio(ostream, pair, n)

        assert self.number_pairs - n_wgs == n_3c, 'Error: WGS and 3C pairs did not sum to ' \
                                                  '{} was did not add'.format(self.number_pairs)

        return {'wgs_count': n_wgs, 'lig_count': self.number_pairs - n_wgs}

    def _simulate_hic(self, ostream):
        """
        A strategy to simulate the sequencing of a HiC library.
        :param ostream: the output stream for reads
        """

        comm = self.community
        read_gen = self.read_generator
        efficiency = self.efficiency

        n_wgs = 0
        n_3c = 0

        for n in tqdm.tqdm(range(1, self.number_pairs+1)):

            ins_len, midpoint, is_fwd = read_gen.draw_insert()

            # is HIC pair?
            if random.pcg_random.uniform() <= efficiency:

                n_3c += 1

                # draw the first replicon and site
                seg1, pos1 = comm.draw_any_by_site()

                # is it spurious ligation
                if comm.is_spurious():

                    # draw any segment and position
                    seg2, pos2 = comm.draw_any_by_site()

                # is it an inter-replicon (trans) ligation
                elif seg1.parent_repl.parent_cell.is_trans():

                    # draw another segment from any other replicon in the cell
                    seg2 = seg1.parent_repl.parent_cell.draw_other_segment_by_sites(seg1.parent_repl)
                    pos2 = seg2.draw_any_site()

                # otherwise an intra-replicon (cis) ligation
                else:

                    # draw another segment from the same replicon
                    seg2 = seg1.parent_repl.draw_any_segment_by_sites()
                    if seg2 == seg1:
                        # must follow defined distribution of Hi-C pair separation
                        pos2 = seg1.draw_constrained_site(pos1)
                    else:
                        # unconstrained, any site on seg2
                        pos2 = seg2.draw_any_site()

                # randomly permute source/destination
                if random.pcg_random.integer(2) == 0:
                    pos1, pos2 = pos2, pos1
                    seg1, seg2 = seg2, seg1

                # with coordinates, make hic read-pair
                pair = read_gen.make_ligation_readpair(seg1, pos1, seg2, pos2, ins_len, midpoint)

            # otherwise WGS
            else:

                # TODO this process will always draw WGS pairs on the same segment
                #   whereas WGS pairs are capable of spanning segments
                n_wgs += 1
                seg1, pos1 = comm.draw_any_by_extent()
                pair = read_gen.make_wgs_readpair(seg1, pos1, ins_len, is_fwd)

            read_gen.write_readpair_dnaio(ostream, pair, n)

        assert self.number_pairs - n_wgs == n_3c, 'Error: WGS and 3C pairs did not sum to ' \
                                                  '{} was did not add'.format(self.number_pairs)

        return {'wgs_count': n_wgs, 'lig_count': n_3c}

    def _simulate_dnase(self, ostream):
        """
        A strategy to simulate the sequencing of a HiC library without a restriction enzyme.
        :param ostream: the output stream for reads
        """
        comm = self.community
        read_gen = self.read_generator
        efficiency = self.efficiency

        n_wgs = 0
        n_3c = 0

        for n in tqdm.tqdm(range(1, self.number_pairs+1)):

            ins_len, midpoint, is_fwd = read_gen.draw_insert()

            # is PLP?
            if random.pcg_random.uniform() <= efficiency:

                n_3c += 1

                # draw the first replicon and site
                seg1, pos1 = comm.draw_any_by_extent()

                # is it spurious ligation
                if comm.is_spurious():

                    seg2, pos2 = comm.draw_any_by_extent()

                # is it an inter-replicon (trans) ligation
                elif seg1.parent_repl.parent_cell.is_trans():

                    seg2 = seg1.parent_repl.parent_cell.draw_other_segment_by_extent(seg1.parent_repl)
                    pos2 = seg2.draw_any_site()

                # otherwise an intra-replicon (cis) ligation
                else:

                    seg2 = seg1.parent_repl.draw_any_segment_by_extent()
                    if seg2 == seg1:
                        # must follow defined distribution of Hi-C pair separation
                        pos2 = seg1.draw_constrained_site(pos1)
                    else:
                        pos2 = seg2.draw_any_site()

                # randomly permute source/destination
                if random.pcg_random.integer(2) == 0:
                    pos1, pos2 = pos2, pos1
                    seg1, seg2 = seg2, seg1

                # with coordinates, make hic read-pair
                pair = read_gen.make_ligation_readpair(seg1, pos1, seg2, pos2, ins_len, midpoint)

            # otherwise WGS
            else:

                n_wgs += 1
                seg1, pos1 = comm.draw_any_by_extent()
                pair = read_gen.make_wgs_readpair(seg1, pos1, ins_len, is_fwd)

            read_gen.write_readpair_dnaio(ostream, pair, n)

        assert self.number_pairs - n_wgs == n_3c, 'Error: WGS and 3C pairs did not sum to ' \
                                                  '{} was did not add'.format(self.number_pairs)

        return {'wgs_count': n_wgs, 'lig_count': n_3c}
