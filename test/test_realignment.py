import numpy
import pytest
import random

import conftest

def test_realignment(genome_source, genome_source_deletion, heterozygous_readpairs_deletion):
    print("")
    genome_source_deletion, deletion_length = genome_source_deletion

    # in theory we could make a fixture that runs this using both bwa and ssw,
    # but that'd be slow so let's skip it for development purposes
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
    """
    make sure that reads slightly overlapping the breakpoint get only a slightly 
    better score than the other allele
    """
    genome_source_deletion, deletion_length = genome_source_deletion

    refseq = genome_source.names_to_contigs["chr2"]
    altseq = genome_source_deletion.names_to_contigs["chr2"]

    isize = 400
    print("")

    # the reads derived from the reference allele
    print("REF "*10)
    for i in range(min(30, deletion_length)):
        cur_reads = [
            conftest.simulate_read_pair(refseq, 4000+deletion_length-i, isize=isize),
            conftest.simulate_read_pair(refseq, 4000-isize+i, isize=isize)
            ]

        for pair in cur_reads:
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

    # reads derived from the alternate allele
    print("ALT "*10)
    for i in range(15):
        cur_reads = [
            conftest.simulate_read_pair(altseq, 4000-i, isize=isize),
            conftest.simulate_read_pair(altseq, 4000-isize+i, isize=isize)
            ]

        for pair in cur_reads:
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


def test_repetitive(genome_source):
    """
    let's:
    1. create a deletion
    2. replicate the sequence around the deletion nearly identically on another chromosome
    3. make sure the reads mapping to the deletion allele map there with low mapq
    """
    alt_contigs = genome_source.names_to_contigs.copy()
    alt_contigs["chr2"] = alt_contigs["chr2"][:4000] + alt_contigs["chr2"][4100:]
    mutated_junction = list(alt_contigs["chr2"][3600:4400])
    r = random.Random()
    r.seed(5125)
    for pos in [10, 100, 300, 390, 440, 500]:
        # pos = r.randint(0, len(mutated_junction))
        mutated_junction[pos] = r.choice(list(set(list("ACGT")) - set(mutated_junction[pos])))

    alt_contigs["chr1"] = alt_contigs["chr1"] + "".join(mutated_junction + [r.choice("ACGT") for i in range(1000)])
    alt = conftest.GenomeSource(alt_contigs)

    altseq = alt_contigs["chr2"]

    for i in [300, 125, 75, 10]:
        pair = conftest.simulate_read_pair(altseq, 4000-i)

        pair.realign([genome_source], [alt])

        ref_scores = sorted([x.score for x in pair.ref_pairs], reverse=True)
        alt_scores = sorted([x.score for x in pair.alt_pairs], reverse=True)

        print(ref_scores, alt_scores)

        assert (alt_scores[0]-alt_scores[1]) < 4
        assert alt_scores[0] > ref_scores[0]
