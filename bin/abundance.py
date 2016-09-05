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
from collections import OrderedDict

import io

import numpy as np


def generate_profile(random_state, taxa, mode, **kwargs):
    """
    Generate a relative abundance profile.

    :param random_state: numpy.RandomState object from which to draw values
    :param taxa: the number of taxa in the profile or a list of names (chrom) or tuples (chrom, cell)
    :param mode: selected mode [equal, uniform or lognormal]
    :param kwargs: additional options for mode. Log-normal requires lognorm_mu, lognorm_sigma
    :return: array of abundance values
    """

    ntax = None
    is_named = False
    if isinstance(taxa, int):
        ntax = taxa
    elif isinstance(taxa, (list, tuple)):
        ntax = len(taxa)
        print ntax
        # use the first element to determine if we've been passed a list of scalars
        if not isinstance(taxa[0], (list, tuple)):
            # convert single chrom names to (chrom, None) explicit tuples with no cell name.
            tx = []
            for ti in taxa:
                tx.append((ti, None))
            taxa = tx
        is_named = True
    else:
        raise RuntimeError('taxa parameter must be a integer or a list/tuple of names')

    # obtain the set of values from the chosen distribution
    if mode == 'equal':
        abn_val = np.full(ntax, 1.0/ntax, dtype=np.float64)
    elif mode == 'uniform':
        abn_val = random_state.uniform(size=ntax)
        abn_val /= abn_val.sum()
    elif mode == 'lognormal':
        abn_val = random_state.lognormal(kwargs['lognorm_mu'], kwargs['lognorm_sigma'], size=ntax)
        abn_val /= abn_val.sum()
    else:
        raise RuntimeError('unsupported mode [{0}]'.format(mode))

    if is_named:
        # names to be inserted in alphabetical order
        ordered_names = sorted(taxa)
        profile = Profile()
        for n, (chrom, cell) in enumerate(ordered_names):
            profile.add(chrom, abn_val[n], cell)
        return profile
    else:
        # otherwise just return a plain list of values
        return abn_val.tolist()


class Abundance:
    """
    An entry in a profile, where object identity is keyed by both chromosome and cell name. Cell names are explicit
    and important for supporting multi-chromosomal genome definitions in simulations of 3C read-pairs, where inter
    and intra chromsomal sampling behaviour is drastically different. In situations where a community simulation is
    entirely monochromosomal.

    The abundance 'value' is expected to be a number.
    """
    def __init__(self, chrom, val, cell=None):
        """

        :param chrom: chromsome name
        :param val: abundance value [0,1]
        :param cell: cell/species name. If none, takes on the name of the chromosome
        """
        self.chrom = chrom
        self.cell = chrom if not cell else cell
        try:
            self.val = float(val)
        except ValueError:
            ValueError('Error: abundance value was not a number [{0}]'.format(val))

    @property
    def long_name(self):
        return '{0}-{1}'.format(self.cell, self.chrom)
    
    def __hash__(self):
        return hash(self.chrom)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.chrom == other.chrom

    def __ne__(self, other):
        return not self.__eq__(other)

    def __cmp__(self, other):
        return self.cell + self.chrom > other.cell + other.chrom

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return self.chrom


class Profile(OrderedDict):
    """
    The Profile class represents an abundance profile of a single community. Relative abundance
    values are keyed by both cell and chromsome, thereby permitting more than mono-chromosomal
    cell definitions and varying copy number.

    Abundance entries are kept in order of entry but are sorted by default when writing
    tabular output.
    """
    def __init__(self, *args, **kwargs):
        super(Profile, self).__init__(*args, **kwargs)

    def add(self, chrom, val, cell=None):
        """
        Convenience method adding an entry to the profile. The Abundance object
        is created internally.
        :param chrom: chromosome name
        :param val: abundance value
        :param cell: cell/species name, defaults to chrom if not specified
        """
        self.addAbundance(Abundance(chrom, val, cell))

    def addAbundance(self, abdn):
        """
        Add an abundance object to the Profile.
        :param abdn: abundance object to add
        """
        if abdn in self:
            raise RuntimeError('Error: duplicate abundances with identity [{0}] in profile'.format(abdn._id()))
        self[abdn] = abdn

    def to_table(self, sort=True):
        """
        Table form, were rows are in the order: [chrom, cell, value]
        :param sort: sort entries before creating table
        :return: table representation of profile
        """
        t = []
        keys = sorted(self.keys()) if sort else self.keys()
        for n, k in enumerate(keys):
            t.append([self[k].chrom, self[k].cell, self[k].val])
        return t

    def write_table(self, hndl, sort=True):
        """
        Write a profile to an output stream as a tab-delimited table.
        :param hndl: output stream handle
        :param sort: Sort by names prior to writing
        """
        t = self.to_table(sort)
        hndl.write('#chrom\tcell\tabundance\n')
        for row in t:
            hndl.write('{0}\t{1}\t{2}\n'.format(row[0], row[1], row[2]))

    def normalize(self):
        """
        Normalize a profile so that all abundance entries sum to 1.
        """
        val_sum = sum([ai.val for ai in self.values()])
        for ai in self.values():
            ai.val /= val_sum


def read_profile(hndl):
    """
    Read a profile from an input stream.
    :param hndl: the input file name or file object
    :return: Profile object
    """

    close_handle = False
    try:
        if isinstance(hndl, basestring):
            hndl = open(hndl, 'r')
            close_handle = True

        profile = Profile()
        n = 1
        for line in hndl:
            line = line.strip()
            if len(line) <= 0:
                continue
            if line.startswith('#'):
                continue
            try:
                chrom, cell, val = line.split('\t')
                profile.add(chrom, val, cell)
            except:
                raise IOError('Error: invalid table at line {0} [{0}]'.format(n, line))
            n += 1
        return profile
    finally:
        if close_handle:
            hndl.close()
