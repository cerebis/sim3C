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

import numpy as np
import re
import yaml


defaults = {}
# read default values
with open('sim3C_config.yaml', 'r') as hndl:
    defaults = yaml.load(hndl)
    defaults.setdefault('anti_rate', 0.1)
    defaults.setdefault('trans_rate', 0.1)
    defaults.setdefault('copy_number', 1)
    defaults.setdefault('linear', False)
    defaults.setdefault('create_cids', False)


def generate_profile(seed, taxa, mode, **kwargs):
    """
    Generate a relative abundance profile.

    :param seed: random state initialisation seed
    :param taxa: the number of taxa in the profile or a list of names (chrom) or tuples (chrom, cell)
    :param mode: selected mode [equal, uniform or lognormal]
    :param kwargs: additional options for mode. Log-normal requires lognorm_mu, lognorm_sigma
    :return: array of abundance values
    """

    random_state = np.random.RandomState(seed)

    is_named = False
    if isinstance(taxa, int):
        ntax = taxa
    elif isinstance(taxa, (list, tuple)):
        ntax = len(taxa)
        print 'Profile will be over {0} taxa'.format(ntax)
        # use the first element to determine if we've been passed a list of scalars
        if not isinstance(taxa[0], (list, tuple)):
            # convert single chrom names to (chrom, None) explicit tuples with no cell name.
            tx = []
            for ti in taxa:
                tx.append((ti, None))
            taxa = tx
        is_named = True
    else:
        raise RuntimeError('taxa parameter must be a integer or '
                           'a list/tuple of names. was [{0}]'.format(taxa.__class__))

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
        for n, (chr_name, cell) in enumerate(ordered_names):
            profile.add(chr_name, abn_val[n], 1, cell)
        return profile
    else:
        # otherwise just return a plain list of values
        return abn_val.tolist()


def as_type(v, type_func, accept_none, msg='Unexpected type'):
    try:
        if accept_none:
            return v
        else:
            return type_func(v)
    except:
        raise ValueError(msg)


def as_float(v, accept_none=False, msg='Float expected'):
    return as_type(v, float, accept_none, msg)


def as_int(v, accept_none=False, msg='Int expected'):
    return as_type(v, int, accept_none, msg)


def as_bool(v, accept_none=False, msg='Bool expected'):
    return as_type(v, bool, accept_none, msg)


class RepliconDefinition:
    """
    An entry in a profile, where object identity is keyed by both chromosome and cell name. Cell names are explicit
    and important for supporting multi-chromosomal genome definitions in simulations of 3C read-pairs, where inter
    and intra chromsomal sampling behaviour is drastically different. In situations where a community simulation is
    entirely monochromosomal.
    """
    def __init__(self, cell, name, copy_number, anti_rate, create_cids, linear, centromere):
        """
        :param cell: containing CellDefintion
        :param name: replicon name
        :param copy_number: chromsome/replicon copy number
        :param anti_rate: anti-diagonal interaction rate
        :param create_cids: create CIDS/TAD interactions
        :param linear: treat replicon as linear
        :param centromere: set centromere location for simulation
        """
        self.cell = cell
        self.name = name
        self.copy_number = as_int(copy_number, msg='copy_number expected be an integer')
        self.anti_rate = as_float(anti_rate, msg='anti_rate expected to be a number')
        self.linear = as_bool(linear, msg='linear expected to be True or False')
        self.create_cids = as_bool(create_cids, msg='create_cids expected to be True or False')
        self.centromere = as_int(centromere, True, msg='centromere expected to be an integer')

    def effective_abundance(self):
        """
        Product of cellular abundance and this chromosome's copy-number.
        :return: abundance * copy_number
        """
        return self.cell.abundance * self.copy_number

    @property
    def long_name(self):
        return '{0}-{1}'.format(self.cell, self.name)

    def to_dict(self):
        d = {'name': self.name,
             'copy_number': self.copy_number,
             'anti_rate': self.anti_rate,
             'linear': self.linear,
             'create_cids': self.create_cids,
             'centromere': self.centromere}
        return d

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.name == other.name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __cmp__(self, other):
        return self.cell > other.cell and self.name > other.name

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return self.name


class CellDefinition:

    def __init__(self, name, seq_file, abundance, trans_rate, replicons):
        self.name = name
        self.seq_file = seq_file
        self.abundance = as_float(abundance, msg='abundance expeected to be a number')
        self.trans_rate = as_float(trans_rate, msg='trans_rate expected to be a number')
        self.repl_repository = OrderedDict()
        for repl_name, repl_info in replicons.iteritems():
            self.add_replicon(RepliconDefinition(cell=self, name=repl_name, **repl_info))

    def add_replicon(self, repl):
        if not isinstance(repl, RepliconDefinition):
            raise ValueError('repl must be a RepliconDefinition')
        if repl in self.repl_repository:
            raise ValueError('Duplicate replicon name [{}]'.format(repl.name))
        self.repl_repository[repl] = repl

    def to_dict(self):
        d = {'seq_file': self.seq_file,
             'abundance': self.abundance,
             'trans_rate': self.trans_rate,
             'replicons': []}
        for ri in self.repl_repository:
            d['replicons'].append(ri.to_dict())
        return d

    def __repr__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.name == other.name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __cmp__(self, other):
        return self.name > other.name

    def __str__(self):
        return repr(self)


class Profile:
    """
    The Profile class represents an abundance profile of a single community. Relative abundance
    values are keyed by both cell and chromsome, thereby permitting more than mono-chromosomal
    cell definitions and varying copy number.

    Abundance entries are kept in order of entry but are sorted by default when writing
    tabular output.
    """
    def __init__(self, ):
        self.cell_repository = OrderedDict()

    def add_cell(self, cell):
        if not isinstance(cell, CellDefinition):
            raise ValueError('cell must be a CellDefinition')
        if cell in self.cell_repository:
            raise ValueError('Duplicate cell name [{}]'.format(cell.name))
        self.cell_repository[cell] = cell

    def to_table(self):
        """
        Table form, were rows are in the order: [chrom, cell, abundance, copy_number]
        :return: table representation of profile
        """
        t = [['name', 'cell', 'seq_file', 'abundance', 'copy_number', 'anti_rate',
              'trans_rate', 'linear', 'create_cids', 'centromere']]
        for ci in sorted(self.cell_repository):
            for ri in sorted(ci.repl_repository):
                t.append([ri.name, ri.cell.name, ci.seq_file, ci.abundance, ri.copy_number, ri.anti_rate,
                          ci.trans_rate, ri.linear, ri.create_cids, ri.centromere])
        return t

    def write_table(self, hndl):
        """
        Write a profile to an output stream as a tab-delimited table.
        :param hndl: output stream handle
        """
        tbl = self.to_table()
        hndl.write('#{}\n'.format('\t'.join(str(v) for v in tbl[0])))
        if len(tbl) > 1:
            for row in tbl[1:]:
                hndl.write('{}\n'.format('\t'.join(str(v) for v in row)))

    def write_yaml(self, hndl):
        out = {}
        for ci in self.cell_repository:
            out[ci.name] = ci.to_dict()
        yaml.dump(out, hndl, default_flow_style=False, indent=4)

    def normalize(self):
        """
        Normalize a profile so that all abundance entries sum to 1.
        """
        val_sum = sum([cell_i.abundance for cell_i in self.cell_repository.values()])
        for cell_i in self.cell_repository:
            cell_i.abundance /= val_sum


def read_table(hndl, normalise=False):
    """
    Read a profile from an input stream. Column position dictates field and must be in the expected
    order -- as written out by write_table()
    :param hndl: the input file name or file object
    :param normalise: when true, normalise after reading
    :return: Profile object
    """
    close_handle = False
    try:
        if isinstance(hndl, basestring):
            hndl = open(hndl, 'r')
            close_handle = True

        profile = Profile()
        for n, line in enumerate(hndl, start=1):
            line = line.strip()
            print line
            if len(line) <= 0:
                continue
            if line.startswith('#'):
                continue
            try:
                # split line and associate with field names
                fields = re.split('[\s,]+', line)
                labels = ['name', 'cell', 'seq_file', 'abundance', 'copy_number', 'anti_rate',
                          'trans_rate', 'linear', 'create_cids', 'centromere']
                t = dict(zip(labels, fields))

                # instantiate cell->repl
                cell = CellDefinition(t['cell'], t['seq_file'], t['abundance'], t['trans_rate'], {})
                cell = profile.cell_repository.setdefault(cell, cell)
                repl = RepliconDefinition(cell, t['name'], t['copy_number'], t['anti_rate'],
                                          t['create_cids'], t['linear'],
                                          t['centromere'])
                cell.add_replicon(repl)
            except Exception as e:
                import traceback
                traceback.print_exc()
                raise IOError('invalid table at line {0} [{0}]'.format(n, line))

        if normalise:
            profile.normalize()

        return profile
    finally:
        if close_handle:
            hndl.close()


def read_yaml(fname):
    """
    Read a community profile definition file in YAML format. The method returns a ordered
    dict of cells, which each contain an ordered dict of replicons, where all fields
    pertinent to the object are populated. Non-mandatory fields not specified in the file
    are returned with values set to None.
    
    Unspecified mandatory fields will raise an IOError.
    
    :param fname: a yaml format profile
    :return: an ordered dict of dicts
    """

    with open(fname, 'r') as hndl:
        cfg = yaml.load(hndl)

    com_def = OrderedDict()
    for cell_name, cell_info in sorted(cfg.items(), key=lambda x: x[0]):

        if 'abundance' not in cell_info:
            raise IOError('Invalid community definition: cell [{}] has no abundance specified'.format(cell_name))

        if 'seq_file' not in cell_info:
            raise IOError('Invalid community definition: cell [{}] has no sequence file specified'.format(cell_name))

        com_def[cell_name] = {
            'trans_rate': cell_info.get('trans_rate', defaults['trans_rate']),
            'abundance': cell_info['abundance'],
            'seq_file': cell_info['seq_file']
        }

        rd = OrderedDict()
        com_def[cell_name]['replicons'] = rd
        for repl_info in sorted(cell_info['replicons'], key=lambda x: x['name']):
            rd[repl_info['name']] = {
                'anti_rate': repl_info.get('anti_rate', defaults['anti_rate']),
                'create_cids': repl_info.get('create_cids', defaults['create_cids']),
                'linear': repl_info.get('linear', defaults['linear']),
                'copy_number': repl_info.get('copy_number', defaults['copy_number']),
                'centromere': repl_info.get('centromere', None)
            }

    profile = Profile()
    for cell_name, cell_info in com_def.iteritems():
        profile.add_cell(CellDefinition(name=cell_name, **cell_info))

    return profile
