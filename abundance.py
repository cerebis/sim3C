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
from Bio import SeqIO

import numpy as np
import re
import yaml


"""
Initialise default values
"""
defaults = {}
# read default values
with open('sim3C_config.yaml', 'r') as hndl:
    defaults = yaml.load(hndl)
    defaults.setdefault('anti_rate', 0.1)
    defaults.setdefault('trans_rate', 0.1)
    defaults.setdefault('copy_number', 1)
    defaults.setdefault('linear', False)
    defaults.setdefault('create_cids', False)
    defaults.setdefault('centromere', None)


def as_type(v, type_func, accept_none, msg='Unexpected type'):
    """
    Cast to type with a message on error
    :param v: value to cast
    :param type_func: casting method
    :param accept_none: accept None as well
    :param msg: message on error
    :return: type cast value
    """
    try:
        if accept_none:
            return v
        else:
            return type_func(v)
    except:
        raise ValueError(msg)


def as_float(v, accept_none=False, msg='Float expected'):
    """
    Cast to float
    :param v: value
    :param accept_none: accept None 
    :param msg: message on error
    :return: float
    """
    return as_type(v, float, accept_none, msg)


def as_int(v, accept_none=False, msg='Int expected'):
    """
    Cast to int
    :param v: value
    :param accept_none: accept None 
    :param msg: message on error
    :return: int
    """
    return as_type(v, int, accept_none, msg)


def as_bool(v, accept_none=False, msg='Bool expected'):
    """
    Cast to boolean
    :param v: value
    :param accept_none: accept None 
    :param msg: message on error
    :return: boolean
    """
    return as_type(v, bool, accept_none, msg)


def generate_profile(seed, seq_file, mode, **kwargs):
    """
    Generate a relative abundance profile.
    :param seed: random state initialisation seed
    :param seq_file: input sequences from which to generate profile
    :param mode: selected mode [equal, uniform or lognormal]
    :param kwargs: additional options for mode. Log-normal requires lognorm_mu, lognorm_sigma
    :return: Profile objecttrans
    """
    random_state = np.random.RandomState(seed)

    seq_index = SeqIO.index(seq_file, 'fasta')
    try:
        taxa = list(seq_index)
    finally:
        seq_index.close()

    ntax = len(taxa)
    print 'Profile will be over {0} taxa'.format(ntax)

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

    # names to be inserted in alphabetical order
    repl_names = sorted(taxa)

    profile = Profile()
    for n, name in enumerate(repl_names):
        cell = CellDefinition(name, seq_file, abn_val[n], defaults['trans_rate'], {})
        cell = profile.cell_repository.setdefault(cell, cell)
        repl = RepliconDefinition(cell, name, defaults['copy_number'], defaults['anti_rate'],
                                  defaults['create_cids'], defaults['linear'], defaults['centromere'])
        cell.add_replicon(repl)

    return profile


class RepliconDefinition:
    """
    Relevant values required for simulation, pertinent to an individual replicon.
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
    """
    Relevant values required for simulation, pertinent to an individual Cell
    """

    def __init__(self, name, seq_file, abundance, trans_rate, replicons):
        """
        :param name: cell name (must by unique across a profile) 
        :param seq_file: file from which all replicons belong to this cell can be found
        :param abundance: abundance of cell
        :param trans_rate: trans (inter-chr) probability
        :param replicons: dictionary of replicons, indexed by themselves. ie. {replicon: replicon}
        """
        self.name = name
        self.seq_file = seq_file
        self.abundance = as_float(abundance, msg='abundance expeected to be a number')
        self.trans_rate = as_float(trans_rate, msg='trans_rate expected to be a number')
        self.repl_repository = OrderedDict()
        for repl_name, repl_info in replicons.iteritems():
            self.add_replicon(RepliconDefinition(cell=self, name=repl_name, **repl_info))

    def add_replicon(self, repl):
        """
        Add a replicon to this cell. Must by unique by name.
        :param repl: RepliconDefinition instance
        """
        if not isinstance(repl, RepliconDefinition):
            raise ValueError('repl must be a RepliconDefinition')
        if repl in self.repl_repository:
            raise ValueError('Duplicate replicon name [{}]'.format(repl.name))
        self.repl_repository[repl] = repl

    def to_dict(self):
        """
        Convert a cell to a dictionary representation. This is used to get a simple
        serialisation with YAML.
        :return: dict
        """
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
    values are keyed by both cell and replicon, thereby permitting more than mono-chromosomal
    cell definitions and varying copy number.

    Abundance entries are kept in order of entry but are sorted by default when writing
    tabular output.
    """
    def __init__(self, ):
        self.cell_repository = OrderedDict()

    def add_cell(self, cell):
        """
        Add a cell to this profile. Cell name must be unique across profile.
        :param cell: CellDefinition instance
        """
        if not isinstance(cell, CellDefinition):
            raise ValueError('cell must be a CellDefinition')
        if cell in self.cell_repository:
            raise ValueError('Duplicate cell name [{}]'.format(cell.name))
        self.cell_repository[cell] = cell

    def to_table(self):
        """
        Table form of this profile, where column definition is positional.
        :return: table representation of profile
        """
        t = [['name', 'cell', 'seq_file', 'abundance', 'copy_number', 'anti_rate',
              'trans_rate', 'linear', 'create_cids', 'centromere']]
        for ci in sorted(self.cell_repository):
            for ri in sorted(ci.repl_repository):
                t.append([ri.name, ri.cell.name, ci.seq_file, ci.abundance, ri.copy_number, ri.anti_rate,
                          ci.trans_rate, ri.linear, ri.create_cids, ri.centromere])
        return t

    def write(self, output, fmt='yaml'):
        """
        Write a profile to output to either 'yaml' or 'table' format.
        :param output: either a filename or open output stream
        :param fmt: yaml or table format
        """
        close_output = False
        try:
            if isinstance(output, basestring):
                output = open(output, 'w')
                close_output = True
            if fmt == 'yaml':
                self.write_yaml(output)
            elif fmt == 'table':
                self.write_table(output)
            else:
                raise RuntimeError('unknown profile format [{}]'.format(fmt))
        finally:
            if close_output:
                output.close()

    def write_table(self, hndl):
        """
        Serialise this profile to a file in tabular form.
        :param hndl: output stream handle
        """
        tbl = self.to_table()
        hndl.write('#{}\n'.format('\t'.join(str(v) for v in tbl[0])))
        if len(tbl) > 1:
            for row in tbl[1:]:
                hndl.write('{}\n'.format('\t'.join(str(v) for v in row)))

    def write_yaml(self, hndl):
        """
        Serialise this profile to a file in YAML format.
        :param hndl: 
        """
        out = {}
        for ci in self.cell_repository:
            out[ci.name] = ci.to_dict()
        yaml.dump(out, hndl, default_flow_style=False, indent=4)

    def normalize(self):
        """
        Normalize a profile so that cell abundances sum to 1.
        """
        val_sum = sum([cell_i.abundance for cell_i in self.cell_repository.values()])
        for cell_i in self.cell_repository:
            cell_i.abundance /= val_sum


def read(input, fmt='yaml', normalize=False):
    """
    Read a profile from an open input stream, either in yaml (fmt=yaml) or tabular format (fmt=table).
    :param input: 
    :param fmt: 
    :param normalize: 
    :return: 
    """
    close_input = False
    if isinstance(input, basestring):
        input = open(input, 'r')
        close_input = True

    try:
        if fmt == 'yaml':
            return read_yaml(input, normalize)
        elif fmt == 'table':
            return read_table(input, normalize)
        else:
            raise RuntimeError('unknown profile format [{}]'.format(fmt))
    finally:
        if close_input:
            input.close()


def read_table(hndl, normalize=False):
    """
    Read a profile from an input stream. Column position dictates field and must be in the expected
    order -- as written out by write_table()
    :param hndl: the input stream
    :param normalize: when true, normalise after reading
    :return: Profile object
    """
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

    if normalize:
        profile.normalize()

    return profile


def read_yaml(hndl, normalize=False):
    """
    Read a community profile definition file in YAML format. The method returns a ordered
    dict of cells, which each contain an ordered dict of replicons, where all fields
    pertinent to the object are populated. Non-mandatory fields not specified in the file
    are returned with values set to None.
    
    Unspecified mandatory fields will raise an IOError.
    
    :param hndl: input stream
    :param normalize: when true, normalise after reading
    :return: an ordered dict of dicts
    """

    cfg = yaml.load(hndl)

    comm_def = OrderedDict()
    for cell_name, cell_info in sorted(cfg.items(), key=lambda x: x[0]):

        if 'abundance' not in cell_info:
            raise IOError('Invalid community definition: cell [{}] has no abundance specified'.format(cell_name))

        if 'seq_file' not in cell_info:
            raise IOError('Invalid community definition: cell [{}] has no sequence file specified'.format(cell_name))

        comm_def[cell_name] = {
            'trans_rate': cell_info.get('trans_rate', defaults['trans_rate']),
            'abundance': cell_info['abundance'],
            'seq_file': cell_info['seq_file']
        }

        rd = OrderedDict()
        comm_def[cell_name]['replicons'] = rd
        for repl_info in sorted(cell_info['replicons'], key=lambda x: x['name']):
            rd[repl_info['name']] = {
                'anti_rate': repl_info.get('anti_rate', defaults['anti_rate']),
                'create_cids': repl_info.get('create_cids', defaults['create_cids']),
                'linear': repl_info.get('linear', defaults['linear']),
                'copy_number': repl_info.get('copy_number', defaults['copy_number']),
                'centromere': repl_info.get('centromere', None)
            }

    profile = Profile()
    for cell_name, cell_info in comm_def.iteritems():
        profile.add_cell(CellDefinition(name=cell_name, **cell_info))

    if normalize:
        profile.normalize()

    return profile
