import numpy

import json
from svviz2.utility import intervals, statistics



def calculate_genotype_likelihoods(ref, alt, priors=[0.05, 0.5, 0.95], max_qual=200):
    """
    calculates the bayesian genotype likelihoods as per Chiang et al (2015)
    """

    ref = int(ref)
    alt = int(alt)

    log_combo = statistics.log_choose(ref+alt, alt)

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


def get_overlaps(read_locus, unsequenced_insert_locus, breakpoints):
    overlaps = {}

    for breakpoint in breakpoints:
        if not read_locus.overlapsAnysense(breakpoint):
            continue
        if len(breakpoint) > 1:
            raise NotImplementedError("breakpoints with size > 1")

        cur_overlap = min([
            breakpoint.start - read_locus.start,
            read_locus.end - breakpoint.start])

        extension = read_locus.end - breakpoint.start

        if unsequenced_insert_locus:
            overlaps_sequence = not unsequenced_insert_locus.overlapsAnysense(breakpoint)
        else:
            overlaps_sequence = True

        overlaps[str(breakpoint)] = (cur_overlap, overlaps_sequence, extension)
        # best_overlap = max(cur_overlap, best_overlap)
        # if cur_overlap > best_overlap:
        #     best_breakpoint = breakpoint
        #     best_overlap = cur_overlap

    return overlaps#best_overlap, best_breakpoint

def set_read_supports_allele(aln_set, aln, allele, score, read_stats, breakpoint_collection, min_overlap):
    if not aln.concordant(read_stats):
        return 0

    assert len(aln.loci) == 1
    aln_locus = aln.loci[0]
    try:
        chrom = aln.aln1.chrom
        start = max(aln.aln1.reference_start, aln.aln2.reference_start)
        end = min(aln.aln1.reference_end, aln.aln2.reference_end)

        unsequenced_insert_locus = intervals.Locus(chrom, start, end, "+")
    except AttributeError:
        unsequenced_insert_locus = None

    try:
        if aln.insert_size > read_stats.max_reasonable_insert_size():
            return 0
        if aln.insert_size < read_stats.min_reasonable_insert_size():
            return 0
    except (IndexError, AttributeError):
        pass

    overlaps = get_overlaps(aln_locus, unsequenced_insert_locus, breakpoint_collection)
    if len(overlaps) == 0:
        return 0

    best_overlap = max(list(zip(*overlaps.values()))[1])

    aln_set.supports_allele = allele
    aln_set.support_prob = score / 40.0 #(1 - mapq.phred_to_prob(score, 10.0))
    aln_set.supporting_aln = aln

    aln.set_tag("OV", best_overlap)
    aln.set_tag("Ov", json.dumps(overlaps))
    aln.overlap = best_overlap

    return aln_set.support_prob


def assign_reads_to_alleles(aln_sets, ref_breakpoint_collection, alt_breakpoint_collection, read_stats):
    def get_best_score(_aln_set, _allele):
        if _allele == "ref":
            alignments = _aln_set.ref_pairs
        elif _allele == "alt":
            alignments = _aln_set.alt_pairs
        if len(alignments) > 0:
            return alignments[0].mapq
        return 0

    ref_total = 0
    alt_total = 0

    for aln_set in aln_sets:
        ref_score = get_best_score(aln_set, "ref")
        alt_score = get_best_score(aln_set, "alt")

        # if aln_set.name == "D00360:99:C8VWFANXX:4:2310:5190:27306":
        aln_set.supports_allele = "amb"
        aln_set.support_prob = 0
        aln_set.supporting_aln = None

        if ref_score - alt_score > 1:
            # print(aln_set.name)
            # print(">REF<")
            # for aln in aln_set.ref_pairs:
            #     print(" ", aln.aln1.locus, aln.aln1.cigarstring, aln.aln1.score)
            #     print(" ", aln.aln2.locus, aln.aln2.cigarstring, aln.aln2.score)
            #     print(" ", aln.mapq)
            # print(">ALT<")
            # for aln in aln_set.alt_pairs:
            #     print(" ", aln.aln1.locus, aln.aln1.cigarstring, aln.aln1.score)
            #     print(" ", aln.aln2.locus, aln.aln2.cigarstring, aln.aln2.score)
            #     print(" ", aln.mapq)

            aln = aln_set.ref_pairs[0]
            ref_total += set_read_supports_allele(
                aln_set, aln, "ref", ref_score, read_stats, ref_breakpoint_collection, min_overlap=4)


        elif alt_score - ref_score > 1:
            aln = aln_set.alt_pairs[0]
            alt_total += set_read_supports_allele(
                aln_set, aln, "alt", alt_score, read_stats, alt_breakpoint_collection, min_overlap=4)
        elif len(aln_set.ref_pairs) > 0:
            aln = aln_set.ref_pairs[0]
            aln_set.supports_allele = "amb"
            aln_set.supporting_aln = aln


    return ref_total, alt_total


def test():
    print(calculate_genotype_likelihoods(0, 53))
    # print(calculate_genotype_likelihoods(13,2))
    # print(calculate_genotype_likelihoods(2,40))
    # print(calculate_genotype_likelihoods(25,26))


if __name__ == '__main__':
    test()