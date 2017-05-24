import pysam
import pytest

from genosv.remap.new_realignment import Alignment
import conftest


def test_get_seq(genome_source):
    seq2 = conftest._seq2()
    seq1 = conftest._seq1()
    
    # test chroms
    assert genome_source.get_seq("chr1", 6, 38, "+") == seq1[6:39]
    assert genome_source.get_seq("chr2", 13, 40, "+") == seq2[13:41]

    # test reverse complement
    rc = seq2[13:41]
    rc = rc[::-1]
    _convert = {"A":"T", "T":"A", "C":"G", "G":"C"}
    rc = "".join(_convert[x] for x in rc)

    assert genome_source.get_seq("chr2", 13, 40, "-") == rc

def test_get_names(genome_source):
    assert set(genome_source.keys()) == set(["chr1", "chr2"])

def test_get_bwa(genome_source):
    bwa = genome_source.bwa
    genome_source.bwa.align
    assert genome_source.bwa == bwa

def test_pickle(genome_source):
    import pickle
    genome_source.bwa
    pickle.dumps(genome_source)



@pytest.mark.parametrize("chrom, start, end", [
    ("chr1", 10, 29),
    ("chr1", 120, 150),
    ("chr2", 5000, 5150),
    ("chr2", 8000, 8400),
])
def test_realign(genome_source, chrom, start, end):
    read = pysam.AlignedSegment()
    read.query_sequence = genome_source.get_seq(chrom, start, end, "+")

    alns = genome_source.align(Alignment(read))
    assert len(alns) == 1
    assert alns[0].cigarstring == "{}M".format(end-start+1)
    assert alns[0].chrom == chrom
    assert alns[0].reference_start == start
    assert alns[0].reference_end == end+1

def test_realign_insertion(genome_source):
    read = pysam.AlignedSegment()
    read.query_sequence = genome_source.get_seq("chr2", 30, 50, "+") + \
                            "GG" + genome_source.get_seq("chr2", 51, 65, "+")

    alns = genome_source.align(Alignment(read))
    assert len(alns) == 1
    assert alns[0].cigarstring == "21M2I15M"
    assert alns[0].reference_start == 30
    assert alns[0].reference_end == 66 # the past-end position


def test_realign_rc(genome_source):
    read = pysam.AlignedSegment()
    read.query_sequence = genome_source.get_seq("chr1", 30, 50, "-")

    alns = genome_source.align(Alignment(read))
    assert len(alns) == 1
    assert alns[0].cigarstring == "21M"
    assert alns[0].reference_start == 30
    assert alns[0].reference_end == 51
    assert alns[0].is_reverse

    qs = "<<<<<<<:<9/,&,22;;<<<"
    read.query_qualities = pysam.qualitystring_to_array(qs)
    alns = genome_source.align(Alignment(read))

    import warnings
    with warnings.catch_warnings():
        # this is a python 2/3 incompatibility I think, where the warning
        # indicates array.tostring() is deprecated but array.tobytes()
        # only exists in py3
        warnings.simplefilter("ignore")
        assert pysam.qualities_to_qualitystring(alns[0].query_qualities) == qs[::-1]

    # raise Exception(1)

def test_realign_align_params(genome_source):
    read = pysam.AlignedSegment()
    read.query_sequence = genome_source.get_seq("chr1", 101, 114, "+")

    # the default seed length should be too short for this to align
    alns = genome_source.align(Alignment(read))
    assert len(alns) == 0

    genome_source.bwa.SetMinSeedLength(13)

    # but now it should align
    alns = genome_source.align(Alignment(read))
    assert len(alns) == 1
    assert alns[0].cigarstring == "14M"
    assert alns[0].reference_start == 101
