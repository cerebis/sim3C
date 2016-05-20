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
import re

class Alignment:
    """
    Basic class representing a line from a PSL format file.
    Where field definitions were taken from: http://asia.ensembl.org/info/website/upload/psl.html
    """
    def __init__(self, line):
        fields = line.strip().split()
        # Number of matching bases that aren't repeats.
        self.matches = int(fields[0])
        # Number of bases that don't match.
        self.mismatches = int(fields[1])
        # Number of matching bases that are part of repeats.
        self.repmatches = int(fields[2])
        # Number of 'N' bases.
        self.nCount = int(fields[3])
        # Number of inserts in query.
        self.q_num_insert = int(fields[4])
        # Number of bases inserted into query.
        self.q_base_insert = int(fields[5])
        # Number of inserts in target.
        self.t_num_insert = int(fields[6])
        # Number of bases inserted into target.
        self.t_base_insert = int(fields[7])
        # defined as + (forward) or - (reverse) for query strand.
        # In mouse, a second '+' or '-' indecates genomic strand.
        self.strand = fields[8]
        # Query sequence name.
        self.q_name = fields[9]
        # Query sequence size.
        self.q_size = int(fields[10])
        # Alignment start position in query.
        self.q_start = int(fields[11])
        # Alignment end position in query.
        self.q_end = int(fields[12])
        # Target sequence name.
        self.t_name = fields[13]
        # Target sequence size.
        self.t_size = int(fields[14])
        # Alignment start position in query.
        self.t_start = int(fields[15])
        # Alignment end position in query.
        self.t_end = int(fields[16])
        # Number of blocks in the alignment.
        self.block_count = int(fields[17])
        # Comma-separated list of sizes of each block.
        self.block_sizes = _delim_string_to_int_list(fields[18])
        # Comma-separated list of start position of each block in query.
        self.q_starts = _delim_string_to_int_list(fields[19])
        # Comma-separated list of start position of each block in target.
        self.t_starts = _delim_string_to_int_list(fields[20])

    @property
    def length(self):
        '''
        Length of the alignment
        '''
        return self.q_end - self.q_start + 1

    @property
    def coverage(self):
        '''
        Coverage fraction of alignment relative to query size
        '''
        return self.length / float(self.q_size)

    @property
    def percent_id(self):
        '''
        Percentage identity of alignment
        As defined by BLAT perl script for calculating percentage identity
        '''
        return (1.0 - float(self.mismatches + self.q_num_insert) /
                float(self.matches + self.mismatches + self.repmatches)) * 100.0


_dataline_pattern = re.compile(r'^[0-9]+\t')


def parse(handle):
    """
    Generator function for parsing a PSL file.
    Any line which does not conform to a basic test of fields will
    be silently skipped.
    :param handle: the open file handle to parse
    :return: Alignment instance
    """
    for line in handle:
        if not _dataline_pattern.match(line):
            continue
        yield Alignment(line)


def _delim_string_to_int_list(str, delim=','):
    """
    Convert a delimited string to a list of ints
    :param str: character delimited string
    :param sep: field delimiter
    :return: list of ints
    """
    return [int(si) for si in str.split(delim) if si]
