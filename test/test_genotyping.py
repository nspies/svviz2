import collections
import numpy
import pytest
import random

import conftest
from svviz2.io import readstatistics
from svviz2.remap import genotyping
from svviz2.utility.intervals import Locus

def get_read_stats(isize=400):
    stats = readstatistics.ReadStatistics(None)
    stats.insertSizes = numpy.random.normal(400, 20, 2000).astype(int)
    stats.orientations = ["+-"]
    return stats

def test_gt(genome_source, genome_source_deletion):
    genome_source_deletion, deletion_length = genome_source_deletion

    refseq = genome_source.names_to_contigs["chr2"]
    altseq = genome_source_deletion.names_to_contigs["chr2"]

    print("")

    coverage = 50
    read_length = 150
    ref_reads = conftest.simulate_read_pairs(refseq, int(len(refseq)/(read_length*2)*coverage))
    alt_reads = conftest.simulate_read_pairs(altseq, int(len(altseq)/(read_length*2)*coverage))
    print(len(ref_reads), len(alt_reads))

    combined_reads = []

    for i, _, pair in ref_reads:
        if 4000-500 < i < 4000+500+deletion_length:
            pair._allele = "ref"
            combined_reads.append(pair)
    for i, _, pair in alt_reads:
        if 4000-500 < i < 4500:
            pair._allele = "alt"
            combined_reads.append(pair)

    for pair in combined_reads:
        pair.realign([genome_source], [genome_source_deletion])

    ref_breakpoints = [Locus("chr2", 4000, 4000, "+"),
                       Locus("chr2", 4000+deletion_length, 4000+deletion_length, "+")]
    alt_breakpoints = [Locus("chr2", 4000, 4000, "+")]

    ref_count, alt_count = genotyping.assign_reads_to_alleles(
        combined_reads, ref_breakpoints, alt_breakpoints, get_read_stats())

    print(":::::", ref_count, alt_count)
