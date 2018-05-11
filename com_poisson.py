from numba import jit
from math import log, exp
 
@jit
def Z(lam, nu, thres):
    """
    Calculate a truncation approximation of the infinite sum Z(lambda, nu) from COM-Poisson distribution.
    An attempt has been made to do this in an efficient way.
    """
    s = 1
    zlast = 0.0
    zsum = 1.0
    loglam = log(lam)
    slogfac = 0.0

    # stop at threshold small change.
    while abs(zsum - zlast) > thres:
        zlast = zsum
        # ca
        zsum += exp(s * loglam - slogfac * nu)
        s += 1
        # log factorial is a sum, since s is incremented, 
        # we can just add each new term to the existing value.
        slogfac += log(s)

    return zsum

def logfac(n):
    if n < 2:
        return 0.0
    return sum(log(n) for n in xrange(2, n+1))

def com_logL(y, lam, nu):
    return log(lam)*sum(y) - nu * sum(logfac(yi) for yi in y) - len(y) * log(Z(lam, nu))

