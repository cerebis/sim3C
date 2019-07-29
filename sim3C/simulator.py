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

from .abundance import read_profile
from .art import Art, EmpDist, ambiguous_base_filter, validator
from .community import Community
from .exceptions import *
from .random import uniform, randint, normal
from .site_analysis import get_enzyme_instance

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
        self.wgs_desc_fmt = 'WGS {repl.name}:{x1}..{x2}:{dir}'
        self._3c_desc_fmt = method.upper() + ' {repl1.name}:{x1} {repl2.name}:{x2}'

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

    def _part_joiner_simple(self, a, b):
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
        length = int(normal(self.insert_mean, self.insert_sd))
        if length < self.insert_min:
            self.too_short += 1
            length = self.insert_min
        elif self.insert_max and length > self.insert_max:
            self.too_long += 1
            length = self.insert_max
        midpoint = randint(0, length)
        is_fwd = uniform() < 0.5
        return length, midpoint, is_fwd

    def make_wgs_readpair(self, repl, x1, ins_len, is_fwd):
        """
        Create a fwd/rev read-pair simulating WGS sequencing.
        :param repl: replicon from which to extract this read-pair
        :param x1: the initial coordinate along the replicon.
        :param ins_len: the insert length
        :param is_fwd: will the insert be off the fwd strand
        :return: a read-pair dict
        """
        frag = repl.subseq(x1, ins_len, is_fwd)
        pair = self.next_pair(str(frag.seq))
        pair['mode'] = 'WGS'
        pair['desc'] = self.wgs_desc_fmt.format(repl=repl, x1=x1, x2=x1+ins_len, dir='F' if is_fwd else 'R')
        return pair

    def make_ligation_readpair(self, repl1, x1, repl2, x2, ins_len, ins_junc):
        """
        Create a fwd/rev read-pair simulating a ligation product (Eg. the outcome of
        HiC or meta3C library prep). As repl1 and repl2 can be the same, these ligation
        products can be inter-rep, intra-rep or spurious products.
        :param repl1: the first replicon
        :param x1:  the location along repl1
        :param repl2: the second replicon
        :param x2: the location along repl2
        :param ins_len: insert length
        :param ins_junc: junction point on insert
        :return: a read-pair dict
        """
        part_a = repl1.subseq(x1 - ins_junc, ins_junc)
        part_b = repl2.subseq(x2, ins_len - ins_junc)

        pair = self.next_pair(str(self._part_joiner(part_a, part_b).seq))
        pair['mode'] = '3C'
        pair['desc'] = self._3c_desc_fmt.format(repl1=repl1, x1=x1, repl2=repl2, x2=x2)
        return pair

    def write_readpair(self, h_out, pair, index, fmt='fastq'):
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


class SequencingStrategy(object):
    """
    A SequencingStrategy represents the whole experiment. This includes the reference data from which
    WGS and ligation products are generated, and the simulation of Illumina sequencing reads.

    Experiments are reproducible by supplying the same seed value.
    """

    Strategy = namedtuple('Strategy', 'method run')

    def __init__(self, prof_filename, seq_filename, enz_name, number_pairs,
                 method, read_length, prefix, machine_profile,
                 insert_mean=400, insert_sd=50, insert_min=50, insert_max=None,
                 anti_rate=0.25, spurious_rate=0.02, trans_rate=0.1,
                 efficiency=0.02,
                 ins_rate=9.e-5, del_rate=1.1e-4,
                 create_cids=True, simple_reads=True, linear=False, convert_symbols=False):
        """
        Initialise a SequencingStrategy.

        :param prof_filename: the abundance profile for the community
        :param seq_filename: the matching sequence of replicon sequences in Fasta format
        :param enz_name: the restriction enzyme name (case sensitive)
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
        """
        self.prof_filename = prof_filename
        self.seq_filename = seq_filename
        self.enz_name = enz_name
        self.number_pairs = number_pairs
        self.simple_reads = simple_reads
        self.method = method
        self.read_length = read_length
        self.insert_min = insert_min
        self.insert_max = insert_max
        self.efficiency = efficiency

        self.enzyme = None if not enz_name else get_enzyme_instance(enz_name)
        self.profile = read_profile(prof_filename, True)

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
        self.community = Community(seq_index, self.profile, self.enzyme, anti_rate=anti_rate,
                                   spurious_rate=spurious_rate, trans_rate=trans_rate,
                                   create_cids=create_cids, linear=linear)

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
            r1, x1 = comm.draw_any_by_extent()
            ins_len, midpoint, is_fwd = read_gen.draw_insert()

            if uniform() < efficiency and r1.covers_site(x1, midpoint):

                n_3c += 1

                # move x1 to the nearest actual site
                x1 = r1.sites.find_nn(x1)

                # is it spurious ligation
                if comm.is_spurious():

                    r2, x2 = comm.draw_any_by_site()

                # is it an inter-replicon (trans) ligation
                elif r1.parent_cell.is_trans():

                    r2 = r1.parent_cell.draw_other_replicon_by_sites(r1.name)
                    x2 = r2.draw_any_site()

                # otherwise an intra-replicon (cis) ligation
                else:
                    # find the first intervening site between
                    # random constrained position x2 and site x1.
                    x2 = r1.draw_constrained_site(x1)
                    x2 = r1.sites.find_first(x2, x1)
                    r2 = r1

                # randomly permute source/destination
                if uniform() < 0.5:
                    x1, x2 = x2, x1
                    r1, r2 = r2, r1

                pair = read_gen.make_ligation_readpair(r1, x1, r2, x2, ins_len, midpoint)

            # otherwise WGS
            else:
                n_wgs += 1
                # take the already drawn coordinates
                pair = read_gen.make_wgs_readpair(r1, x1, ins_len, is_fwd)

            read_gen.write_readpair(ostream, pair, n)

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
            if uniform() <= efficiency:

                n_3c += 1

                # draw the first replicon and site
                r1, x1 = comm.draw_any_by_site()

                # is it spurious ligation
                if comm.is_spurious():

                    r2, x2 = comm.draw_any_by_site()

                # is it an inter-replicon (trans) ligation
                elif r1.parent_cell.is_trans():

                    r2 = r1.parent_cell.draw_other_replicon_by_sites(r1.name)
                    x2 = r2.draw_any_site()

                # otherwise an intra-replicon (cis) ligation
                else:

                    x2 = r1.draw_constrained_site(x1)
                    r2 = r1

                # randomly permute source/destination
                if uniform() < 0.5:
                    x1, x2 = x2, x1
                    r1, r2 = r2, r1

                # with coordinates, make hic read-pair
                pair = read_gen.make_ligation_readpair(r1, x1, r2, x2, ins_len, midpoint)

            # otherwise WGS
            else:
                n_wgs += 1
                r1, x1 = comm.draw_any_by_extent()
                pair = read_gen.make_wgs_readpair(r1, x1, ins_len, is_fwd)

            read_gen.write_readpair(ostream, pair, n)

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
            if uniform() <= efficiency:

                n_3c += 1

                # draw the first replicon and site
                r1, x1 = comm.draw_any_by_extent()

                # is it spurious ligation
                if comm.is_spurious():

                    r2, x2 = comm.draw_any_by_extent()

                # is it an inter-replicon (trans) ligation
                elif r1.parent_cell.is_trans():

                    r2 = r1.parent_cell.draw_other_replicon_by_extents(r1.name)
                    x2 = r2.draw_any_site()

                # otherwise an intra-replicon (cis) ligation
                else:

                    x2 = r1.draw_constrained_site(x1)
                    r2 = r1

                # randomly permute source/destination
                if uniform() < 0.5:
                    x1, x2 = x2, x1
                    r1, r2 = r2, r1

                # with coordinates, make hic read-pair
                pair = read_gen.make_ligation_readpair(r1, x1, r2, x2, ins_len, midpoint)

            # otherwise WGS
            else:
                n_wgs += 1
                r1, x1 = comm.draw_any_by_extent()
                pair = read_gen.make_wgs_readpair(r1, x1, ins_len, is_fwd)

            read_gen.write_readpair(ostream, pair, n)

        assert self.number_pairs - n_wgs == n_3c, 'Error: WGS and 3C pairs did not sum to ' \
                                                  '{} was did not add'.format(self.number_pairs)

        return {'wgs_count': n_wgs, 'lig_count': n_3c}
