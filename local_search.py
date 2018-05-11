from random import randrange
import numpy as np
from numba import jit, int64, int32, int16


def moc_crossover(a, b, place_holder=-1):
    """
    Modified Order Crossover
    """
    ri = np.random.randint(2, len(a) - 2)

    a2 = np.zeros_like(a)
    a2[:] = place_holder
    b2 = np.zeros_like(a)
    b2[:] = place_holder

    for i in xrange(len(a)):
        if np.isin(a[i, 0], b[:ri, 0]):
            a2[i] = a[i]
        if np.isin(b[i, 0], a[:ri, 0]):
            b2[i] = b[i]

    n, m = 0, 0
    for i in xrange(len(a)):
        if a2[i, 0] == place_holder:
            a2[i] = b[ri + n]
            n += 1
        if b2[i, 0] == place_holder:
            b2[i] = a[ri + m]
            m += 1

    a[:, :] = a2[:, :]
    b[:, :] = b2[:, :]

    return a, b


def uopx_crossover(a, b, p=0.5, place_holder=-1):
    """
    Uniform Order Preserving Crossover
    """
    # mask of keepers
    mask = np.random.choice([False, True], size=len(a), p=[1 - p, p])

    a2 = np.zeros_like(a)
    a2[:] = place_holder
    b2 = np.zeros_like(a)
    b2[:] = place_holder

    # keepers
    a2[mask, :] = a[mask, :]
    b2[mask, :] = b[mask, :]
    # those not kept
    a_unused = a[~mask, 0]
    b_unused = b[~mask, 0]
    # extract "those not kept" in other order
    b_to_a = []
    a_to_b = []
    for i in xrange(len(a)):
        if np.isin(b[i, 0], a_unused):
            b_to_a.append(b[i])
        if np.isin(a[i, 0], b_unused):
            a_to_b.append(a[i])

    # transfer
    a2[~mask] = b_to_a[:]
    b2[~mask] = a_to_b[:]

    a[:, :] = a2[:, :]
    b[:, :] = b2[:, :]

    return a, b


def ox1_crossover(a, b):
    """
    OX1 crossover operator
    """

    a2 = np.zeros_like(a)
    b2 = np.zeros_like(a)

    i1, i2 = np.random.choice(range(len(a)), replace=False, size=2)
    if i1 > i2:
        i1, i2 = i2, i1

    # preserve state range
    a2[i1:i2 + 1] = a[i1:i2 + 1]
    b2[i1:i2 + 1] = b[i1:i2 + 1]

    # prep fill sources with 0-base
    a_donor = np.roll(a, -(i2 + 1), axis=0)
    b_donor = np.roll(b, -(i2 + 1), axis=0)

    i = (i2 + 1) % len(a)
    for bi in b_donor:
        if not np.isin(bi[0], a2):
            a2[i] = bi
            i += 1
            i %= len(a)

    i = (i2 + 1) % len(a)
    for ai in a_donor:
        if not np.isin(ai[0], b2):
            b2[i] = ai
            i += 1
            i %= len(a)

    return a2, b2


def ox2_crossover(a, b):
    """
    OX2 Crossover operator
    """
    a2 = a.copy()
    b2 = b.copy()

    # selected indices for order transfer
    min_d = max(1, len(a) / 4)
    max_d = max(2, len(a) / 2)
    _size = np.random.randint(min_d, max_d + 1)
    ix = np.random.choice(range(len(a)), replace=False, size=_size)
    ix.sort()

    # values ordered in b
    val_order_b = b2[ix]
    n = 0
    for i in xrange(len(a2)):
        # transfer b order to a
        if np.isin(a2[i, 0], val_order_b[:, 0]):
            print a2[i], val_order_b[n]
            a2[i] = val_order_b[n]
            n += 1

    # values ordered in a
    val_order_a = a2[ix]
    n = 0
    for i in xrange(len(b2)):
        # transfer a order to b
        if np.isin(b2[i, 0], val_order_a[:, 0]):
            b2[i] = val_order_a[n]
            n += 1

    return a2, b2


def mpx_crossover(a, b, place_holder=-1):
    """
    Maximal preservative crossover
    """
    assert len(a) > 4, 'minimum individual length is 4'

    min_d = max(1, len(a) / 4)
    max_d = max(2, len(a) / 2)

    ri = np.random.randint(len(a))
    d = np.random.randint(min_d, max_d + 1)

    a2 = np.empty_like(a)
    a2[:] = place_holder
    b2 = np.empty_like(b)
    b2[:] = place_holder

    ix = np.arange(ri, ri + d) % len(a)
    a2[ix] = a[ix]

    ix = np.arange(ri, ri + d) % len(b)
    b2[ix] = b[ix]
    not_ix = np.setxor1d(np.arange(len(b)), ix)

    donor_a = []
    donor_b = []
    not_a = a[not_ix]
    not_b = b[not_ix]
    for i in xrange(len(a)):
        if np.isin(b[i, 0], not_a[:, 0]):
            donor_a.append(b[i])
        if np.isin(a[i, 0], not_b[:, 0]):
            donor_b.append(a[i])

    a2[not_ix] = donor_a
    b2[not_ix] = donor_b

    a[:, :] = a2[:, :]
    b[:, :] = b2[:, :]

    return a, b


def double_bridge_move(a):
    """
    Double bridge permutation
    """
    seg_len = len(a) / 4
    p1 = 1 + randrange(0, seg_len)
    p2 = p1 + 1 + randrange(0, seg_len)
    p3 = p2 + 1 + randrange(0, seg_len)
    return np.vstack((a[:p1], a[p3:], a[p2:p3], a[p1:p2]))


@jit([int16[:, :](int16[:, :], int16, int16),
      int32[:, :](int32[:, :], int32, int32),
      int64[:, :](int64[:, :], int64, int64)])
def two_opt_move(_a, _i, _j):
    b = _a.copy()
    b[_i:_j+1] = _a[_i:_j+1][::-1]
    return b


@jit([int16[:, :](int16[:, :], int16, int16),
      int32[:, :](int32[:, :], int32, int32),
      int64[:, :](int64[:, :], int64, int64)])
def two_opt_move_flip(_a, _i, _j):
    b = _a.copy()
    b[_i:_j+1] = _a[_i:_j+1][::-1]
    b[_i:_j+1, 1] = (b[_i:_j+1, 1] == 0)
    return b


@jit([int16[:, :](int16[:, :], int16),
      int32[:, :](int32[:, :], int32),
      int64[:, :](int64[:, :], int64)])
def flip(_a, _i):
    b = _a.copy()
    b[_i, 1] = (_a[_i, 1] == 0)
    return b


from tqdm import trange


def two_opt(tour, fitness, temp):
    """
    2-opt local search algorithm
    """
    unif = np.random.uniform
    best_tour = tour
    best_fit = fitness(tour)[0]

    while True:
        improved = False
        t_outer = xrange(len(tour))
        #t_outer.set_description('Current {}'.format(fit_cur))
        for i in t_outer:
            for j in xrange(i + 1, len(tour)):
                # if unif() < 0.0:
                #     a_alt = flip(a, j)
                # else:
                tour_alt = two_opt_move(tour, i, j)
                fit_alt = fitness(tour_alt)[0]
                if fit_alt > best_fit:
                    improved = True
                    best_tour = tour_alt
                    best_fit = fit_alt
                    #t_outer.set_description('Current {}'.format(fit_cur))
                # elif np.log(np.random.uniform()) < (fit_alt - fit_cur)/temp:
                #     tour = tour_alt
                #     fit_cur = fit_alt
                #     print 'RND', fit_cur
                    #t_outer.set_description('Random  {}'.format(fit_cur))
        if improved:
            print 'Improved:', best_fit
            tour = best_tour
        else:
            print 'No improvement found, stop'
            break

    return tour

