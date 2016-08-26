from collections import Iterable
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
    isNamed = False
    if isinstance(taxa, int):
        ntax = taxa
    elif isinstance(taxa, Iterable):
        ntax = len(taxa)
        isNamed = True

    if mode == 'flat':
        if kwargs:
            raise RuntimeWarning('mode=flat does not accept additional options [{0}]'.format(kwargs))
        prf = np.full(ntax, 1.0/ntax, dtype=np.float64)
    elif mode == 'uniform':
        if kwargs:
            raise RuntimeWarning('mode=uniform does not accept additional options [{0}]'.format(kwargs))
        prf = random_state.uniform(size=ntax)
        prf /= prf.sum()
    elif mode == 'lognormal':
        prf = random_state.lognormal(kwargs['lognorm_mu'], kwargs['lognorm_sigma'], size=ntax)
        prf /= prf.sum()
    else:
        raise RuntimeError('unsupported mode [{0}]'.format(mode))

    # just return a plain Python list
    if isNamed:
        named_prf = {}
        for n, ti in taxa:
            named_prf[ti] = prf[n]
        return named_prf
    else:
        return prf.tolist()
