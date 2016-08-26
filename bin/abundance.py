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
from collections import Iterable, OrderedDict
import numpy as np


def relative_profile(random_state, taxa, mode, **kwargs):
    """
    Generate a relative abundance profile.

    :param random_state: numpy.RandomState object from which to draw values
    :param taxa: the number of taxa in the profile or a list of names
    :param mode: selected mode [equal, uniform or lognormal]
    :param kwargs: additional options for mode. Log-normal requires lognorm_mu, lognorm_sigma
    :return: array of abundance values
    """

    ntax = None
    is_named = False
    if isinstance(taxa, int):
        ntax = taxa
    elif isinstance(taxa, Iterable):
        ntax = len(taxa)
        is_named = True

    if mode == 'equal':
        prf = np.full(ntax, 1.0/ntax, dtype=np.float64)
    elif mode == 'uniform':
        prf = random_state.uniform(size=ntax)
        prf /= prf.sum()
    elif mode == 'lognormal':
        prf = random_state.lognormal(kwargs['lognorm_mu'], kwargs['lognorm_sigma'], size=ntax)
        prf /= prf.sum()
    else:
        raise RuntimeError('unsupported mode [{0}]'.format(mode))

    # just return a plain Python list
    if is_named:
        named_prf = OrderedDict()
        ordered_names = sorted(list(taxa))
        for n, ti in enumerate(ordered_names):
            named_prf[ti] = prf[n]
        return named_prf
    else:
        return prf.tolist()


def print_abundance(profile, ostream):
    for k, v in profile.iteritems():
        print '{0}\t{1:.3f}'.format(k, v)
