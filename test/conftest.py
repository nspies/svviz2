import pytest

from genosv.remap.new_realignment import GenomeSource

def _seq1():
    seq = "AGCGGCGAGATCTATTTTATCGAGAGACCACACCCCATCTCTCGAGCTTCGATCGGGGCTCTAGGATCTCTTCTAGGGGGC" \
          "TAAGGCGGAGGGAGCGGGGGGGAAAATTCTATATACCCCCCACAAATATATAAAGAAAGAAAAGCCCTCTTTTATCGGAAA"
    return seq

def _seq2():
    import random
    r = random.Random()
    r.seed(135)
    return "".join(r.choice("ACGT") for i in range(10000))

@pytest.fixture
def genome_source():
    gs = GenomeSource({"chr1":_seq1(),
                       "chr2":_seq2()})
    return gs