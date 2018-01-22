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


def phred_to_prob(phred, phred_scale):
    if phred < 0:
        phred = 0
    return 10.0 ** (phred / -phred_scale)

def prob_to_phred(phred, scale):
    if phred <= 0:
        return numpy.inf

    return -scale * numpy.log10(phred)


def log_sum_exp(array, base=10.0):
    base = float(base)
    array = numpy.asarray(array)
    m = array.max()
    diff = array - m
    sum_of_exps = (base**(diff)).sum()
    return m + numpy.log(sum_of_exps)/numpy.log(base)