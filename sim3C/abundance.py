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

import logging
import numpy as np
import re

from .random import np_uniform, np_lognormal


logger = logging.getLogger(__name__)


def generate_profile(taxa, mode, **kwargs):
    """
    Generate a relative abundance profile.

    :param taxa: the number of taxa in the profile or a list of names (chrom) or tuples (chrom, cell)
    :param mode: selected mode [equal, uniform or lognormal]
    :param kwargs: additional options for mode. Log-normal requires lognorm_mu, lognorm_sigma
    :return: array of abundance values
    """

    is_named = False
    if isinstance(taxa, int):
        ntax = taxa
    elif isinstance(taxa, (list, tuple)):
        ntax = len(taxa)
        if ntax == 0:
            logger.error('The supplied list of taxa was empty')
        logger.info('Profile will be over {} taxa'.format(ntax))
        # use the first element to determine if we've been passed a list of scalars
        if not isinstance(taxa[0], (list, tuple)):
            # convert single chrom names to (chrom, None) explicit tuples with no cell name.
            tx = []
            for ti in taxa:
                tx.append((ti, None, None))
            taxa = tx
        is_named = True
    else:
        raise RuntimeError('taxa parameter must be a integer or '
                           'a list/tuple of names. was [{}]'.format(taxa.__class__))

    # obtain the set of values from the chosen distribution
    if mode == 'equal':
        abn_val = np.full(ntax, 1 / ntax, dtype=np.float64)
    elif mode == 'uniform':
        abn_val = np_uniform(size=ntax)
        abn_val /= abn_val.sum()
    elif mode == 'lognormal':
        abn_val = np_lognormal(kwargs['lognorm_mu'], kwargs['lognorm_sigma'], size=ntax)
        abn_val /= abn_val.sum()
    else:
        raise RuntimeError('unsupported mode [{}]'.format(mode))

    if is_named:
        # names to be inserted in alphabetical order
        ordered_names = sorted(taxa)
        profile = Profile()
        for n, (chr_name, cell, molecule) in enumerate(ordered_names):
            profile.add(chr_name, abn_val[n], 1, cell, molecule)
        return profile
    else:
        # otherwise just return a plain list of values
        return abn_val.tolist()


class ChromAbundance(object):
    """
    An entry in a profile, where object identity is keyed by both chromosome and cell name. Cell names are explicit
    and important for supporting multi-chromosomal genome definitions in simulations of 3C read-pairs, where inter
    and intra chromsomal sampling behaviour is drastically different. In situations where a community simulation is
    entirely monochromosomal.
    """
    def __init__(self, name, abundance, copy_number, cell=None, molecule=None):
        """

        :param name: chromsome name
        :param abundance: cellular abundance [0..1]
        :param copy_number: chromsome/replicon copy number
        :param cell: cell/species name. If none, takes on the name of the chromosome
        """
        self.name = name
        self.cell = name if not cell else cell
        self.molecule = name if not molecule else molecule

        try:
            self.abundance = float(abundance)
        except ValueError:
            raise ValueError('abundance expected to be number')

        try:
            self.copy_number = int(copy_number)
        except ValueError:
            raise ValueError('copy_number expected be an integer')

    def effective_abundance(self):
        """
        Product of cellular abundance and this chromosome's copy-number.
        :return: abundance * copy_number
        """
        return self.abundance * self.copy_number

    @property
    def long_name(self):
        return '{}-{}-{}'.format(self.cell, self.molecule, self.name)
    
    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.name == other.name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __cmp__(self, other):
        return self.cell + self.molecule + self.name > other.cell + other.molecule + other.name

    def __lt__(self, other):
        return self.cell + self.molecule + self.name < other.cell + other.molecule + other.name

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return self.name


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

    def add(self, chr_name, abundance, copy_number, cell=None, molecule=None):
        """
        Convenience method adding an entry to the profile. The Abundance object
        is created internally.
        :param chr_name: chromosome name
        :param abundance: cellular abundance float: [0..1]
        :param copy_number: chromosome/replicon copy number int:[>0]
        :param cell: cell/species name, defaults to chrom if not specified
        :param molecule: molecule name, defaults to chrom if not specified
        """
        self.add_abundance(ChromAbundance(chr_name, abundance, copy_number, cell, molecule))

    def add_abundance(self, abn):
        """
        Add an abundance object to the Profile.
        :param abn: abundance object to add
        """
        if abn in self:
            raise RuntimeError('Error: duplicate abundances with identity [{}] in profile'.format(abn))
        self[abn] = abn

    def to_table(self, sort=True):
        """
        Table form, where rows are in the order: [chrom, cell, abundance, copy_number]
        :param sort: sort entries before creating the table
        :return: a table representation of the profile
        """
        t = []
        keys = sorted([*self]) if sort else [*self]
        for n, k in enumerate(keys):
            t.append([self[k].name, self[k].cell, self[k].molecule, self[k].abundance, self[k].copy_number])
        return t

    def write_table(self, hndl, sort=True):
        """
        Write a profile to an output stream as a tab-delimited table.
        :param hndl: output stream handle
        :param sort: Sort by names prior to writing
        """
        t = self.to_table(sort)
        hndl.write('#chrom\tcell\tmolecule\tabundance\tcopy_number\n')
        for row in t:
            hndl.write('{}\t{}\t{}\t{:f}\t{:d}\n'.format(*row))

    def normalize(self):
        """
        Normalize a profile so that all abundance entries sum to 1.
        """
        val_sum = sum([ai.abundance for ai in self.values()])
        for ai in self.values():
            ai.abundance /= val_sum


def read_profile(hndl, normalise=False):
    """
    Read a profile from an input stream.
    :param hndl: the input file name or file object
    :param normalise: when true, normalise after reading
    :return: Profile object
    """

    close_handle = False
    try:
        if isinstance(hndl, str):
            hndl = open(hndl, 'rt')
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
                chrom, cell, molecule, abn, cn = re.split(r'[\s,]+', line)
                profile.add(chrom, abn, cn, cell, molecule)
            except Exception:
                raise IOError('Error: invalid table at line {} [{}]'.format(n, line))
            n += 1

        if normalise:
            profile.normalize()

        return profile
    finally:
        if close_handle:
            hndl.close()
