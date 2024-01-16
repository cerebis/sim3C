import logging
import random as py_random
import numpy.random as np_random
from .faster import *

logger = logging.getLogger(__name__)

# make local references to chosen methods.
np_uniform = np_random.uniform
np_randint = np_random.randint
np_choice = np_random.choice
np_normal = np_random.normal
np_lognormal = np_random.lognormal

pcg_random = None

LOW_SEED = 1_000_000
HIGH_SEED = 10_000_000


def init_state(seed):
    """
    Initialise random generator APIs using a specific seed value.
    Currently, initialisation the entire Numpy API.

    :param seed: the integer seed to use for initialisation
    """
    global pcg_random

    if seed is None:
        seed = np_random.randint(LOW_SEED, HIGH_SEED)
        logger.debug(f'Drawing random seed for state initialisation')

    np_random.seed(seed)
    py_random.seed(seed)
    # draw two values from out initial state for PCG
    state, inc = np_randint(low=LOW_SEED, high=HIGH_SEED, size=2)
    pcg_random = PCGRandom(state, inc)

    logger.debug(f'Random state initialized with seed: {seed}')
