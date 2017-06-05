import numpy

try:
    import scipy.misc
    def log_choose(n,k):
        return numpy.log10(scipy.misc.comb(n,k))
except ImportError:
    def log_choose(n, k):
        r = 0.0
        # swap for efficiency if k is more than half of n
        if k * 2 > n:
            k = n - k

        for  d in range(1,k+1):
            r += numpy.log10(n)
            r -= numpy.log10(d)
            n -= 1

        return r
