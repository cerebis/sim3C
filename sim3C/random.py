import logging
import random as py_random
import numpy.random as np_random

logger = logging.getLogger(__name__)


# make local references to chosen methods.
# TODO profile np.random against random. This would best be done after unit tests created, so as to validate changes.
uniform = np_random.uniform
randint = np_random.randint
choice = np_random.choice
normal = np_random.normal
lognormal = np_random.lognormal


def init_random_state(seed):
    """
    Initialise random generator APIs using a specific seed value.
    Currently initialisation the entire Numpy API.

    :param seed: the integer seed to use for initialisation
    """
    np_random.seed(seed)
    py_random.seed(seed)
    logger.debug('Initialized random state using seed {}'.format(seed))
