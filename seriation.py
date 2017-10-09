#!/usr/bin/env python
from collections import deque
from operator import lt, gt
import numpy as np
import lap


def sk_bistochastic(m, round_to=8, verbose=False, max_iter=1000):
    """
    Normalise a matrix to be bistochastic using Sinkorn-Knopp algorithm.

    Note: normalisation is done in place.

    This implementation has an awkward means of determinining convergence, where
    each matrix element is rounded to a specified precision and compared to unity.
    That is, at convergence each element should actually be equal to 1.

    :param m: the input matrix
    :param round_to: rounding precision in decimal places
    :param verbose: print debug info
    :param max_iter: maximum number of iterations before abandoning.
    :return: the normalised matrix
    """
    rs = np.empty(len(m), dtype=np.float)
    cs = np.empty(len(m), dtype=np.float)
    i = 0
    while i < max_iter:
        np.sum(m, 1, dtype=np.float, out=rs)
        m = (m.T / rs).T
        np.sum(m, 0, dtype=np.float, out=cs)
        m /= cs
        # continue until convergence in both rows and columns
        if np.all(np.round(rs, round_to) == 1.0) and np.all(np.round(cs, round_to) == 1.0):
            break
        i += 1

    if i == max_iter:
        print 'Weight matrix did to converge to doubly stochastic in 1000 iterations'
    if verbose:
        print 'It took {} iterations to achieve bistochasticity'.format(i)

    return m


def kr_bistochastic(m, tol=1e-6, x0=None, delta=0.1, Delta=3, verbose=False, max_iter=1000):
    """
    Normalise a matrix to be bistochastic using Knight-Ruiz algorithm. This method is expected
    to converge more quickly.

    :param m: the input matrix
    :param tol: precision tolerance
    :param x0: an initial guess
    :param delta: how close balancing vector can get
    :param Delta: how far balancing vector can get
    :param verbose: print debug info
    :param max_iter: maximum number of iterations before abandoning.
    """
    n = len(m)
    e = np.ones(n)

    if not x0:
        x0 = e

    g = 0.9
    etamax = 0.1
    eta = etamax
    stop_tol = tol * 0.5

    x = x0
    rt = tol ** 2
    v = x * np.dot(m, x)

    rk = 1 - v
    rho_km1 = np.dot(rk.T, rk)  # transpose possibly implicit
    rout = rho_km1
    rold = rout

    n_iter = 0
    i = 0

    while rout > rt and n_iter < max_iter:

        i += 1
        k = 0
        y = e
        inner_tol = np.maximum(rout * eta ** 2, rt)

        while rho_km1 > inner_tol:

            k += 1
            if k == 1:
                Z = rk / v
                p = Z
                rho_km1 = np.dot(rk.T, Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p

            w = x * np.dot(m, x * p) + v * p
            alpha = rho_km1 / np.dot(p.T, w)
            ap = alpha * p

            ynew = y + ap

            if np.amin(ynew) <= delta:
                if delta == 0:
                    break
                ind = np.where(ap < 0)[0]
                gamma = np.amin((delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            if np.amax(ynew) >= Delta:
                ind = np.where(ynew > Delta)[0]
                gamma = np.amin((Delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            y = ynew
            rk = rk - alpha * w
            rho_km2 = rho_km1

            Z = rk * v
            rho_km1 = np.dot(rk.T, Z)

        x = x * y
        v = x * np.dot(m, x)
        rk = 1 - v
        rho_km1 = np.dot(rk.T, rk)
        rout = rho_km1
        n_iter += k + 1

        rat = rout / rold
        rold = rout
        res_norm = np.sqrt(rout)
        eta_o = eta
        eta = g * rat

        if g * eta_o ** 2 > 0.1:
            eta = np.maximum(eta, g * eta_o ** 2)
        eta = max(min(eta, etamax), stop_tol / res_norm)

    if verbose:
        print 'It took {} iterations to achieve bistochasticity'.format(n_iter)

    if n_iter >= max_iter:
        print 'Warning: maximum number of iterations ({}) reached without convergence'.format(max_iter)

    X = np.diagflat(x)
    return np.dot(X.T, np.dot(m, X))


def make_symmetric(m, inplace=True):
    """
    Assuming a matrix is symmetric but only represented by its upper triangle, copies
    the transpose to lower triangle.
    :param m: the matrix to make symmetric
    :param inplace: do the copy to the supplied matrix
    :return: the full symmetric (possibly copy)
    """
    if not inplace:
        m_out = m.copy()
    else:
        m_out = m
    ix = np.tril_indices_from(m_out, k=-1)
    m_out[ix] = m_out.T[ix]
    return m_out


def additive_recip(x, a):
    """
    Elementwise reciprocal of a matrix, where a scalar value is added to each prior to
    becoming the denominator. Eg. 1/(x+a)
    :param x: the input matrix
    :param a: the scalar additive
    :return: the elementwise reciprocal
    """
    return 1.0/(x+a)


def log_additive_recip(x, a):
    """
    Elementwise natural log reciprocal of a matrix, where a scalar is added to each prior to
    becoming the denomiator. Eg. 1/log(x+a)
    :param x: the matrix
    :param a: the scalar additive
    :return: the elementwise reciprocal
    """
    return 1.0/np.log(x+a)


def create_large_weight(n, sigma):
    """
    Ugly but with much better scaling as N becomes large (>1000). Ultimately
    this method is faster than the more succinct single line used below.
    :param n: size of square matrix
    :param sigma: weight factor controlling width of effect
    :return: unnormalized weight matrix
    """
    # make a 1d vector representing the entire double-sided weight function.
    v = np.fromfunction(lambda i, c0: np.exp(-c0*i**2), (n,), dtype=np.float, c0=1.0/(n*sigma))
    v = np.hstack((v[:0:-1], v))

    # roll the weight vector as we drop down the rows, initialising the entire row
    x = np.empty((n, n))
    for i in xrange(n):
        x[i, :] = np.roll(v, i)[n-1:]
    return x


def create_weight(n, sigma, verbose=False, max_iter=1000, use_kr=True):
    """
    Create a doubly stochastic (balanced) weight matrix using the traditional algorithm.
    :param n: size of matrix
    :param sigma: weight factor
    :param verbose: print info
    :param max_iter: maximum iterations in balancing
    :param use_kr: use Knight-Ruiz to balance matrix, otherwise Sinkhorn-Knopp
    :return:
    """
    # switch methods depending on requested size.
    if n < 500:
        w = np.fromfunction(lambda i, j, c0: np.exp(-c0*((i+1)-(j+1))**2), (n, n), dtype=np.float, c0=1.0/(n*sigma))
    else:
        w = create_large_weight(n, sigma)

    if use_kr:
        w = kr_bistochastic(w, tol=1.0e-12, delta=0.0, verbose=verbose, max_iter=1000)
    else:
        w = sk_bistochastic(w, round_to=8, verbose=verbose, max_iter=1000)

    return w


def linear_assignment_jv(m, **kwargs):
    """
    Determining the permutation matrix requires solving a linear assignment problem, matching
    rows to columns in the matrix M = DW. The following method uses the Jonker-Volgenant
    algorithm, which is an exact solution. Time complexity may prove to be an issue for large
    matrices, but this is offset by the fact that each iteration of seriation should move
    to a smaller value.

    :param m: the matrix to solve LAP
    :param kwargs: additional arguments for generic support
    :return: order of columns matched to rows.
    """
    return lap.lapjv(m, return_cost=False)[1]


def linear_assignment_heuristic(m, cost_func, random_state, **kwargs):
    """
    A heuristic (inexact) solution to the linear assignment problem -- as proposed by the authors
    of SPIN_NH (Tsafrir et al, 2005). The solution may not always find the minimum cost and
    therefore more iterations are often required. Time complexity is better than that for JV.

    heuristic for the linear assignment problem.
    primary sequence (param 2) is column arg_func (min/max) and
    tie-breaking made possible by secondary sequence (param 1) of random numbers.

    :param m: the matrix to solve LAP
    :param cost_func: eg np.argmax or np.argmin.
    :param random_state: used to draw permutations
    :param kwargs: additional arguments for generic support
    :return: order of columns approx matched to rows.
    """
    return np.lexsort((random_state.permutation(len(m)), cost_func(m, axis=1)))


def seriate_spin_nh(x, sigma=None, max_step_iter=20, weight_func=create_weight, verbose=False, maximize=False,
                    sigma_steps=10, max_sigma=20., min_sigma=1., energy_tol=1e-6, sample_space='linear',
                    use_jv=False, seed=None):
    """
    Matrix sorting by SPIN_NH (Sorting points into neighborhoods) algorithm. (Tsafir et al, 2005)

    :param x: the input matrix
    :param sigma: explicit sigma values in descending order.
    :param max_step_iter: maximum number of iterations per sigma
    :param weight_func: weighting function. Default gaussian diagonal (toeplitz)
    :param verbose: verbose output
    :param maximize: maximise instead of minimise the energy
    :param sigma_steps: number of steps between maximum and minimum requested sigma
    :param max_sigma: maximum (starting) sigma
    :param min_sigma: minimum (final) sigma
    :param energy_tol: change in energy tolerance, below which the iterations at sigma s stops
    :param sample_space: choose to sample between max and min sigma linear or log space
    :param use_jv: during linear assignment, use more precise but computationally expensive O(n^3) Jonker-Volgenant
    :return: the lowest energy permutation matrix
    """
    if not seed:
        random_state = np.random.RandomState()
    else:
        random_state = np.random.RandomState(seed)

    max_small = 3
    max_noimp = 5

    # when not sigma supplied, sample between min and max.
    if not sigma:
        if sample_space == 'linear':
            sigma = np.linspace(max_sigma, min_sigma, sigma_steps, endpoint=True)
        elif sample_space == 'log':
            sigma = np.logspace(np.log10(max_sigma), np.log10(min_sigma), sigma_steps, endpoint=True)

    # sigmas as a queue, we'll just pop them as we go
    sigma = deque(sigma)

    d = np.asarray(x)
    d = d.astype(np.float)
    n = len(d)

    # first sigma
    s = sigma.popleft()

    # create the weight matrix
    if verbose:
        print "\nBeginning sigma: {}".format(s)
    w_orig = weight_func(n, s, verbose)
    w = w_orig.copy()

    p_best = np.empty_like(d)
    m = np.empty_like(d, np.float)

    # normally SPIN minimises energy of a dissimilarity matrix, but you
    # could just as well maximise a similarity. This has not be tested
    # extensively.
    if maximize:
        if use_jv:
            raise RuntimeError('Maximisation not currently supported when using JV algorithm')
        energy_best = -np.inf
        arg_func = np.argmax
        best_energy = gt
    else:
        energy_best = np.inf
        arg_func = np.argmin
        best_energy = lt

    # set up the linear assignment solver.
    if use_jv:
        lap_order = linear_assignment_jv
    else:
        lap_order = linear_assignment_heuristic

    # for making animations of sort
    # n_anim = 1
    # plt.imshow(np.log10(cmap+1), cmap='PuRd', interpolation=None)
    # plt.title('start')
    # plt.savefig('anim/00001.png')
    # plt.close()

    # if we get too many small changes in a row, we stop.
    # Here "small" is defined by energy tolerance.
    small_change = 0
    noimp = 0

    iter_all = 0
    iter_step = 0
    while True:

        if verbose:
            print 'Iteration {:3d}:'.format(iter_all+1),

        # original paper had M.t = D x W.t
        # M = D x W
        np.matmul(d, w, m)

        o = lap_order(m, cost_func=arg_func, random_state=random_state)

        # TODO p could just be cleared rather than instantiated each round
        p = np.zeros((n, n))
        for _j in xrange(n):
            p[_j, o[_j]] = 1.

        # original paper E = trace( P.t x M.t )
        energy_new = np.einsum('ii', np.matmul(p, m))
        if verbose:
            print "best energy: {} new energy: {}".format(energy_best, energy_new)

        # stopping threshold test
        energy_delta = np.fabs(energy_new - energy_best)
        if energy_delta < energy_tol:
            small_change += 1
        else:
            small_change = 0

        # was energy improved?
        if best_energy(energy_new, energy_best):
            energy_best = energy_new
            p_best[:] = p[:]
            # for making animations of sort
            # n_anim += 1
            # plt.imshow(np.log10(np.dot(np.dot(p, cmap), p.T)+1), cmap='PuRd', interpolation=None)
            # plt.title('{:.3f} {} {:.2f}'.format(s, i+1, energy_new))
            # plt.savefig('anim/{:05d}.png'.format(n_anim))
            # plt.close()
        else:
            noimp += 1

        # adapt sigma
        if small_change >= max_small or noimp >= max_noimp or iter_step >= max_step_iter-1:

            if len(sigma) == 0:
                if verbose:
                    print "\nFinished in {} iterations.".format(iter_all+1),
                break

            # next sigma
            s = sigma.popleft()

            if verbose:
                print "\nReducing sigma to: {}".format(s)

            # new weight function from sigma
            w_orig = weight_func(n, s, verbose)

            # recalculate best energy
            np.matmul(p.T, w_orig, w)
            np.matmul(d, w, m)
            energy_best = np.einsum('ii', np.matmul(p, m))

            if verbose:
                print "best energy is now: {} \n".format(energy_best)

            small_change = 0
            noimp = 0
            iter_step = 0

        else:
            np.matmul(p.T, w_orig, w)
            iter_step += 1

        iter_all += 1

    if verbose:
        print "Final energy: {}".format(energy_best)

    return p_best


def weak_to_end(m, names, threshold):
    """
    Move columns with weak interaction to end of matrix. Weak column are defined as
    those whose diagonal element is more than "threshold" of the total weight of the
    column. These are sequences with very low interaction with other sequecnes and
    should be ignored during seriation.

    In a metagenomic setting, weak interaction does not indicate an error as a closed
    mono-chromosomal genome will largely only interact with itself.

    This is done in place.

    :param m: matrix to permute
    :param names: sequence names corresponding to matrix m
    :param threshold: relative threshold above which
    :return: 4-tuple -- filtered matrix, filtered names, number of weak entries, number of zero-count entries
    """

    # first remove zero count sequences
    ix = m.sum(axis=0) > 0
    non_empty = (m[ix])[:, ix]
    names = names[ix]
    n_empty = len(ix) - np.sum(ix)
    if n_empty > 0:
        print 'After zero-count removal {}x{}'.format(*non_empty.shape)

    # now find weak sequences and remove
    diags = np.diagonal((non_empty+1) / non_empty.sum(axis=0))
    wk_ix = np.where(diags > threshold)[0]
    n_wk = len(wk_ix)

    if n_wk > 0:
        o = np.concatenate((np.where(diags <= threshold)[0], wk_ix))
        p = np.zeros_like(non_empty)
        for i in xrange(len(p)):
            p[i, o[i]] = 1.

        m = np.dot(np.dot(p, non_empty), p.T)
        names = names[o]

    return m, names, n_wk, n_empty


if __name__ == '__main__':
    import argparse
    import matplotlib.pyplot as plt
    import mapio
    import time

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--seed', metavar='INT', type=int, default=int(time.time()),
                        help="Random seed for initialising number generator")
    parser.add_argument('--headers', action='store_true', default=False,
                        help='Input CSV file has matching row and column headers list sequence names')
    parser.add_argument('--delim', default=',', help='Matrix delimiter [,]')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    parser.add_argument('--exclude-weak', metavar='THRES', default=None, type=float,
                        help='Exclude weak near-singletons during seriation [keep all]')
    parser.add_argument('--reduce', default=None, type=int,
                        help='Reduce matrix size by taking the first N rows and columns')
    parser.add_argument('--maximize', default=False, action='store_true', help='Maximize energy rather than minimize')
    parser.add_argument('--smoothing', default=5, type=float,
                        help='Per-element additive smoothing scalar used in reciprocal [5]')
    parser.add_argument('--log-recip', default=False, action='store_true', help='Use log-scale elementwise reciprocal')
    parser.add_argument('--use-jv', default=False, action='store_true', help='Use more precise JV linear assignment')
    parser.add_argument('--symm', default=False, action='store_true', help='Make half-matrix symmetric')
    parser.add_argument('--sample-space', default='linear', choices=['linear','log'],
                        help='Sample spapce for sigma [linear]')
    parser.add_argument('--steps', type=int, default=10, help='Number of steps in simga reduction [10]')
    parser.add_argument('--etol', type=float, default=1e-6, help='Convergence tolerance [1e-6]')
    parser.add_argument('--max-sigma', type=float, default=20., help='Maximum initial sigma value [20]')
    parser.add_argument('--min-sigma', type=float, default=1., help='Minimum final sigma value [1]')
    parser.add_argument('--sigma-iter', type=int, default=5, help='Number of settling iterations per sigma value [5]')
    parser.add_argument('-f', '--format', choices=['csv', 'h5'], default='csv',
                        help='Input contact map format')
    parser.add_argument('--cmap-name', default='PuRd', help='Colour palette used for heatmap (matplotlib) [PuRd]')
    parser.add_argument('map', help='Contact map')
    parser.add_argument('output', help='Output base')
    args = parser.parse_args()

    if args.steps < 1:
        raise RuntimeError('There must be at least one step.')

    if args.max_sigma <= args.min_sigma:
        raise RuntimeError('Sigma maximum must be greater than minimum')

    # read in contact map
    cmap, seq_names = mapio.read_map(args.map, args.format, delim=args.delim, names=args.headers)

    # set aside seq name information if supplied in input
    if not args.headers:
        seq_names = np.arange(len(cmap))

    print 'Initial matrix size {}x{}'.format(*cmap.shape)

    if args.symm:
        print 'Making symmetric'
        cmap = make_symmetric(cmap)

    if args.exclude_weak:
        cm_sort, seq_names, n_wk, n_empty = weak_to_end(cmap, seq_names, args.exclude_weak)
        if n_wk > 0:
            cmap = cm_sort[:-n_wk, :-n_wk]
            seq_names = seq_names[:-n_wk]
            print 'Post weak entry removal {}x{}'.format(*cmap.shape)

    if args.reduce:
        if args.reduce > len(cmap):
            raise RuntimeError('Reduction size is larger than starting matrix')
        cmap = cmap[:args.reduce, :args.reduce]
        seq_names = seq_names[:args.reduce]
        print 'User requested reduction to {}x{}'.format(*cmap.shape)

    target_map = cmap.copy()
    if not args.maximize:
        if args.log_recip:
            if args.smoothing <= 1:
                raise RuntimeError('Log-scale additive smoothing requires requires a smoothing factor > 1.')
            target_map = log_additive_recip(target_map.astype(np.float), args.smoothing)
        else:
            if args.smoothing < 0:
                raise RuntimeError('Additive smoothing requires non-negative smoothing factor.')
            target_map = additive_recip(target_map.astype(np.float), args.smoothing)

    perm = seriate_spin_nh(target_map, verbose=args.verbose, maximize=args.maximize, sigma_steps=args.steps,
                           max_sigma=args.max_sigma, min_sigma=args.min_sigma, max_step_iter=args.sigma_iter,
                           sample_space=args.sample_space, use_jv=args.use_jv, energy_tol=args.etol, seed=args.seed)

    # extract indices order from permutation matrix
    order_ix = np.argwhere(perm == 1)[:, 1]

    # write a two-column table of original and final orders.
    np.savetxt('{}_order.csv'.format(args.output), np.vstack((seq_names, seq_names[order_ix])).T, fmt='%s')

    # permute the input contact map
    cmap = np.dot(np.dot(perm, cmap), perm.T)

    # write permuted map
    mapio.write_map(cmap, args.output, args.format, names=seq_names[order_ix])

    fig = plt.figure()
    fig.set_size_inches(10, 10)
    plt.imshow(np.log10(cmap+1), cmap=args.cmap_name, interpolation=None)
    plt.savefig('{}.png'.format(args.output), dpi=360)
