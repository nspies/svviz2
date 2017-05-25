import numpy
import pytest
import random

import conftest

def test_realignment(genome_source, genome_source_deletion, heterozygous_readpairs_deletion):
    print("")
    genome_source_deletion, deletion_length = genome_source_deletion

    genome_source.aligner_type = "bwa"
    genome_source_deletion.aligner_type = "bwa"

    for allele, pairs in heterozygous_readpairs_deletion.items():
        # if allele == "amb":
        sample_size = 1000
        if len(pairs) > sample_size:
            pairs = random.sample(pairs, sample_size)

        print(allele, len(pairs))

        for i, pair in enumerate(pairs):
            # if i % 100 == 0:
            #     print(i)
            pair.realign([genome_source], [genome_source_deletion])
            ref_score = max([x.score for x in pair.ref_pairs]) if len(pair.ref_pairs) > 0 else -numpy.inf
            alt_score = max([x.score for x in pair.alt_pairs]) if len(pair.alt_pairs) > 0 else -numpy.inf
            if allele == "alt":
                assert ref_score < alt_score
            elif allele == "ref":
                assert ref_score > alt_score
            elif allele == "amb":
                assert ref_score == pytest.approx(alt_score)
            # scores.append(max(pair.ref_pairs+pair.alt_pairs, key=lambda x: x.score))

        # print(numpy.mean(scores))

def test_slight_overlap(genome_source, genome_source_deletion):
    genome_source_deletion, deletion_length = genome_source_deletion

    refseq = genome_source.names_to_contigs["chr2"]
    altseq = genome_source_deletion.names_to_contigs["chr2"]

    isize = 400
    read_stats = conftest._ReadStats()
    print("")

    print("REF "*10)
    for i in range(min(30, deletion_length)):
        cur_reads = [
            conftest.simulate_read_pair(refseq, 4000+deletion_length-i, isize=isize),
            conftest.simulate_read_pair(refseq, 4000-isize+i, isize=isize)
            ]

        for read1, read2 in cur_reads:
            pair = conftest.ReadPair(conftest.Alignment(read1), conftest.Alignment(read2), read_stats)
            pair.realign([genome_source], [genome_source_deletion])
            ref_score = max([x.score for x in pair.ref_pairs]) if len(pair.ref_pairs) > 0 else -numpy.inf
            alt_score = max([x.score for x in pair.alt_pairs]) if len(pair.alt_pairs) > 0 else -numpy.inf

            print(i, ref_score, alt_score, ref_score-alt_score)

            if i < 5:
                assert 0 <= (ref_score-alt_score) < 2
            elif i < 8:
                assert 0 <= (ref_score-alt_score) < 10
            elif i > 10:
                assert (ref_score-alt_score) > 5

    print("ALT "*10)
    for i in range(15):
        cur_reads = [
            conftest.simulate_read_pair(altseq, 4000-i, isize=isize),
            conftest.simulate_read_pair(altseq, 4000-isize+i, isize=isize)
            ]

        for read1, read2 in cur_reads:
            pair = conftest.ReadPair(conftest.Alignment(read1), conftest.Alignment(read2), read_stats)
            pair.realign([genome_source], [genome_source_deletion])
            ref_score = max([x.score for x in pair.ref_pairs]) if len(pair.ref_pairs) > 0 else -numpy.inf
            alt_score = max([x.score for x in pair.alt_pairs]) if len(pair.alt_pairs) > 0 else -numpy.inf

            print(i, ref_score, alt_score, ref_score-alt_score)

            if i < 5:
                assert 0 <= (alt_score-ref_score) < 2
            elif i < 8:
                assert 0 <= (alt_score-ref_score) < 10
            elif i > 10:
                assert (alt_score-ref_score) > 5

