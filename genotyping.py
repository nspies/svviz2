import numpy

import mapq
import utilities

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


def calculate_genotype_likelihoods(ref, alt, priors=[0.05, 0.5, 0.95], max_qual=200):
    """
    calculates the bayesian genotype likelihoods as per Chiang et al (2015)
    """

    ref = int(ref)
    alt = int(alt)

    log_combo = log_choose(ref+alt, alt)

    log_prob_homref = log_combo + alt * numpy.log10(priors[0]) + ref * numpy.log10(1-priors[0])
    log_prob_het    = log_combo + alt * numpy.log10(priors[1]) + ref * numpy.log10(1-priors[1])
    log_prob_homalt = log_combo + alt * numpy.log10(priors[2]) + ref * numpy.log10(1-priors[2])

    # This is the "genotype likelihoods", aka GL
    log_probs = numpy.array([log_prob_homref, log_prob_het, log_prob_homalt])
    
    log_prob_sum = numpy.log10((10**log_probs).sum())
    genotype_qualities = 1-(10**log_probs/10**log_prob_sum)
    # print("::", genotype_qualities)
    with numpy.errstate(divide="ignore"):
        phred_genotype_qualities = numpy.abs(-10 * numpy.log10(genotype_qualities))
    phred_genotype_qualities[phred_genotype_qualities>max_qual] = max_qual
    return log_probs, phred_genotype_qualities


def assign_reads_to_alleles(aln_sets, ref_breakpoint_collection, alt_breakpoint_collection):
    def get_best_score(_aln_set, _allele):
        alignments = _aln_set.get_alns(_allele)
        if len(alignments) > 0:
            return alignments[0].mapq
        return 0

    ref_total = 0
    alt_total = 0

    print(ref_breakpoint_collection)
    print(alt_breakpoint_collection)
    for aln_set in aln_sets:
        ref_score = get_best_score(aln_set, "ref")
        alt_score = get_best_score(aln_set, "alt")

        if aln_set.name == "D00360:99:C8VWFANXX:4:2310:5190:27306":
            print("ALT:", [x.mapq for x in aln_set.get_alns("alt")])
            print("REF:", [x.mapq for x in aln_set.get_alns("ref")])

        aln_set.supports_allele = "amb"
        aln_set.support_prob = 0
        aln_set.supporting_aln = None

        if ref_score - alt_score > 1:
            aln = aln_set.get_alns("ref")[0]
            if utilities.overlaps(aln.locus, ref_breakpoint_collection):
                aln_set.supports_allele = "ref"
                aln_set.support_prob = (1 - mapq.phred_to_prob(ref_score, 10.0))
                aln_set.supporting_aln = aln
                ref_total += aln_set.support_prob
        elif alt_score - ref_score > 1:
            aln = aln_set.get_alns("alt")[0]
            if utilities.overlaps(aln.locus, alt_breakpoint_collection):
                aln_set.supports_allele = "alt"
                aln_set.support_prob = (1 - mapq.phred_to_prob(alt_score, 10.0))
                aln_set.supporting_aln = aln
                alt_total += aln_set.support_prob

    return ref_total, alt_total


def test():
    print(calculate_genotype_likelihoods(10,2))
    print(calculate_genotype_likelihoods(13,2))
    print(calculate_genotype_likelihoods(2,40))
    print(calculate_genotype_likelihoods(25,26))


if __name__ == '__main__':
    test()