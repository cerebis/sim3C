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
