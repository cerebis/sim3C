import numpy as np
import random


def make_symmetric(m, inplace=True):
    if not inplace:
        m_out = m.copy()
    else:
        m_out = m
    ix = np.tril_indices_from(m_out, k=-1)
    m_out[ix] = m_out.T[ix]
    return m_out


dissim_reciprocal = np.vectorize(lambda x: 1/x if x > 0 else 4)


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


def create_weight(n, sigma, verbose=False, max_iter=1000):
    """
    Create a doubly stochastic (balanced) weight matrix using the traditional algorithm.
    :param n: size of matrix
    :param sigma: weight factor
    :param verbose: print info
    :param max_iter: maximum iterations in balancing
    :return: 
    """
    # switch methods depending on requested size.
    if n < 500:
        w = np.fromfunction(lambda i, j, c0: np.exp(-c0*((i+1)-(j+1))**2), (n, n), dtype=np.float, c0=1.0/(n*sigma))
    else:
        w = create_large_weight(n, sigma)

    rs = np.empty(n, dtype=np.float)
    cs = np.empty(n, dtype=np.float)
    for i in xrange(max_iter):
        np.sum(w, 1, dtype=np.float, out=rs)
        w = (w.T / rs).T
        np.sum(w, 0, dtype=np.float, out=cs)
        w /= cs
        # continue until convergence in both rows and columns
        if np.all(np.round(rs, 5) == 1.0) and np.all(np.round(cs, 5) == 1.0):
            break
    if i == 999:
        print 'Weight matrix did to converge to doubly stochastic in 1000 iterations'
    if verbose:
        print 'It took {} iterations to make W doubly stochastic'.format(i)
    return w


def seriate_spin_nh(x, sigma=None, step_iters=5, weight_func=create_weight, verbose=False, maximize=False,
                    sigma_steps=10, max_sigma=20., min_sigma=1.):

    if not sigma:
        sigma = np.linspace(max_sigma, min_sigma, sigma_steps)

    d = x.copy()
    d = d.astype(np.float)

    n = len(d)

    # weight matrix
    w_orig = weight_func(n, sigma[0], verbose)
    w = w_orig.copy()

    p_best = np.empty_like(d)
    m = np.empty_like(d, np.float)

    if maximize:
        energy_best = -np.inf
        argfunc = np.argmax
        best_energy = lambda e, best: e > best
    else:
        energy_best = np.inf
        argfunc = np.argmin
        best_energy = lambda e, best: e < best

    for i in xrange(len(sigma)*step_iters):

        if verbose:
            print 'Iteration {}...'.format(i+1),

        np.matmul(d, w, m)

        # heuristic for the linear assignment problem
        # (second argument to order breaks ties randomly)
        o = np.lexsort((random.sample(xrange(n), n), argfunc(m, axis=1)))
        p = np.zeros((n, n))
        for _j in xrange(n):
            p[_j, o[_j]] = 1.

        energy_new = np.einsum('ii', np.matmul(p, m))
        if verbose:
            print "best energy: {} new energy: {}".format(energy_best, energy_new)

        # was energy improved?
        if best_energy(energy_new, energy_best):
            energy_best = energy_new
            p_best[:] = p[:]

        # adapt sigma
        if ((i+1) % step_iters) == 0 and (i+1) != len(sigma)*step_iters:
            s = sigma[i / step_iters + 1]
            if verbose:
                print "\nReducing sigma to: {}".format(s)

            w_orig = weight_func(n, s, verbose)

            # recalculate best energy
            np.matmul(p.T, w_orig, w)
            np.matmul(d, w, m)
            energy_best = np.einsum('ii', np.matmul(p, m))

            if verbose:
                print "best energy is now: {} \n".format(energy_best)

        else:
            np.matmul(p.T, w_orig, w)

    if verbose:
        print "Best Energy: {}".format(energy_best)

    return p_best


if __name__ == '__main__':
    import argparse
    import matplotlib.pyplot as plt
    import mapio

    parser = argparse.ArgumentParser()
    parser.add_argument('--headers', action='store_true', default=False,
                        help='Input CSV file has matching row and column headers list sequence names')
    parser.add_argument('--delim', default=',', help='Matrix delimiter [,]')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    parser.add_argument('--reduce', default=None, type=int,
                        help='Reduce matrix size by taking the first N rows and columns')
    parser.add_argument('--maximize', default=False, action='store_true', help='Maximize energy rather than minimize')
    parser.add_argument('--recip', default=False, action='store_true', help='Apply per-element reciprocal')
    parser.add_argument('--symm', default=False, action='store_true', help='Make half-matrix symmetric')
    parser.add_argument('--steps', type=int, default=10, help='Number of steps in simga reduction [10]')
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
    cmap = mapio.read_map(args.map, args.format, delim=args.delim, names=args.headers)

    # set aside seq name information if supplied in input
    if args.headers:
        cmap, rnames, cnames = cmap
        # just keeping one record of input sequence order, once we establish that our two sources agree.
        if rnames != cnames:
            raise RuntimeError('Row and column headers do not agree')
        seq_names = np.array(rnames)
        del cnames
        del rnames

    print 'Contact matrix size {}x{}'.format(*cmap.shape)

    if args.reduce:
        if args.reduce > len(cmap):
            raise RuntimeError('Reduction size is larger than starting matrix')
        cmap = cmap[:args.reduce, :args.reduce]
        seq_names = seq_names[:args.reduce]

    if args.symm:
        cmap = make_symmetric(cmap)

    target_map = cmap.copy()
    if args.recip:
        target_map = dissim_reciprocal(target_map.astype(np.float))

    perm = seriate_spin_nh(target_map, verbose=args.verbose, maximize=args.maximize, sigma_steps=args.steps,
                           max_sigma=args.max_sigma, min_sigma=args.min_sigma, step_iters=args.sigma_iter)

    # extract indices order from permutation matrix
    order_ix = np.argwhere(perm == 1)[:, 1]
    # write a two-column table of original and final orders.
    np.savetxt('the_order.csv', np.vstack((seq_names, seq_names[order_ix])).T, fmt='%s')

    # permute the input contact map
    cmap = np.dot(np.dot(perm, cmap), perm.T)

    # write permuted map
    mapio.write_map(cmap, args.output, args.format, field_as=np.int)

    fig = plt.figure()
    fig.set_size_inches(10, 10)
    plt.imshow(np.log10(cmap+1), cmap=args.cmap_name, interpolation=None)
    plt.savefig('{}.png'.format(args.output), dpi=360)
