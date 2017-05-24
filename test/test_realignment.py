import numpy
import pytest
import random

def test_realignment(genome_source, genome_source_deletion, heterozygous_readpairs_deletion):
    print("")
    genome_source_deletion, deletion_length = genome_source_deletion

    genome_source.aligner_type = "ssw"
    genome_source_deletion.aligner_type = "ssw"

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
