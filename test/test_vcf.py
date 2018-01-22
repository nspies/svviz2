import numpy
import pytest
from unittest import mock

from svviz2.app.datahub import DataHub
from svviz2.io.vcfparser import VCFParser

from svviz2.app.variants import Deletion, SequenceDefinedVariant

VCFHEADER = """
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=13,length=115169878>
##contig=<ID=17,length=81195210>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##SAMPLE=<ID=s1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
""".lstrip()

# This is straight out of the VCF spec
VCFEXAMPLE_BND = """
2	321681	bnd_W	G	G]17:198982]	6	PASS	SVTYPE=BND;MATEID=bnd_Y	GT	0/1
2	321682	bnd_V	T	]13:123456]T	6	PASS	SVTYPE=BND;MATEID=bnd_U	GT	0/1
13	123456	bnd_U	C	C[2:321682[	6	PASS	SVTYPE=BND;MATEID=bnd_V	GT	0/1
13	123457	bnd_X	A	[17:198983[A	6	PASS	SVTYPE=BND;MATEID=bnd_Z	GT	0/1
17	198982	bnd_Y	A	A]2:321681]	6	PASS	SVTYPE=BND;MATEID=bnd_W	GT	0/1
17	198983	bnd_Z	C	[13:123457[C	6	PASS	SVTYPE=BND;MATEID=bnd_X	GT	0/1
""".lstrip()

VCFEXAMPLE_SEQUENCE_DEFINED = """
1	30100011	EVENT1	G	GAACCGCAGAAACCCACACCTCTTTCTCGGGACCGAACGGGCCCC	100	.	SVTYPE=INS	GT	0/1
2	220010533	EVENT2	AATCCCAGCTA	AAAGAAAAGAATTATCCACCCAG	100	.	SVTYPE=INS;END=220010543	GT	0/1
17	45895177	EVENT3	G	GGGTACCACACCCGTGT	100	.	SVTYPE=INS	GT	0/1
2	31506461	EVENT4	TCTGATTACCTGCTCCACCTGACTCATT	A	.	.	SVLEN=27;SVTYPE=DEL;END=31506488	GT	0/1
""".lstrip()

VCFEXAMPLE_DEL = """
2	321682	EVENT5	T	<DEL>	6	PASS	SVTYPE=DEL;END=321887;SVLEN=-205	GT:GQ	0/1:12
""".lstrip()

@pytest.fixture
def breakend_vcf(tmpdir):
    vcf_path = str(tmpdir.join("svs.vcf"))

    with open(vcf_path, "w") as vcf:
        vcf.write(VCFHEADER)
        vcf.write(VCFEXAMPLE_BND)

    return vcf_path

def compare_breakpoints(observed, expected):
    """
    this is a bit complicated because the strand is asymmetric, depending on order
    (that is, if chr1:a+/chr2:b+ and chr2:b-/chr1:a- are the same event, but 
    opposite order and therefore opposite signs)
    """
    assert len(expected) == len(observed) == 2

    if (observed[0].chrom, observed[0].start, observed[0].strand) == expected[0]:
        if (observed[1].chrom, observed[1].start, observed[1].strand) == expected[1]:
            return
        else:
            raise Exception("failed to match: {} vs {}".format(expected, observed))
    else:
        observed = observed[1], observed[0]

        assert observed[0].chrom == expected[0][0]
        assert observed[0].start == expected[0][1]
        assert observed[0].strand != expected[0][2]

        assert observed[1].chrom == expected[1][0]
        assert observed[1].start == expected[1][1]
        assert observed[1].strand != expected[1][2]

        return


def test_breakend_vcf(breakend_vcf):
    """
    just test to make sure the correct events are hooked up together with
    the proper coordinates and orientations
    """
    datahub = DataHub()
    datahub.args = mock.Mock(variants=breakend_vcf)
    datahub.align_distance = 1000

    print("")
    parser = VCFParser(datahub)
    for v in parser.get_variants():
        if v.name in ["bnd_W", "bnd_Y"]:
            compare_breakpoints(v.breakpoints, [("2", 321680, "+"), ("17", 198981, "-")])
        elif v.name in ["bnd_V", "bnd_U"]:
            compare_breakpoints(v.breakpoints, [("2", 321681, "-"), ("13", 123455, "-")])
        elif v.name in ["bnd_X", "bnd_Z"]:
            compare_breakpoints(v.breakpoints, [("13", 123456, "-"), ("17", 198982, "+")])

@pytest.fixture
def sequence_defined_vcf(tmpdir):
    vcf_path = str(tmpdir.join("svs.vcf"))

    with open(vcf_path, "w") as vcf:
        vcf.write(VCFHEADER)
        vcf.write(VCFEXAMPLE_SEQUENCE_DEFINED)

    # print(open(vcf_path).read())
    return vcf_path

def test_sequence_defined_vcf(sequence_defined_vcf):
    datahub = DataHub()
    datahub.args = mock.Mock(variants=sequence_defined_vcf)
    datahub.align_distance = 1000

    print("")
    parser = VCFParser(datahub)
    for v in parser.get_variants():
        chrom, start = v.breakpoints[0].chrom, v.breakpoints[0].start
        end = v.breakpoints[0].end
        if v.name == "EVENT1":
            assert (chrom, start) == ("1", 30100010)
        elif v.name == "EVENT2":
            assert (chrom, start, end) == ("2", 220010532, 220010532+11-1)


@pytest.fixture
def deletion_vcf(tmpdir):
    vcf_path = str(tmpdir.join("svs.vcf"))

    with open(vcf_path, "w") as vcf:
        vcf.write(VCFHEADER)
        vcf.write(VCFEXAMPLE_DEL)

    # print(open(vcf_path).read())
    return vcf_path


def test_sequence_deletion_vcf(deletion_vcf):
    datahub = DataHub()
    datahub.args = mock.Mock(variants=deletion_vcf)
    datahub.align_distance = 1000

    print("")
    parser = VCFParser(datahub)
    variants = list(parser.get_variants())

    assert len(variants) == 1
    for v in variants:
        chrom, start, end = v.breakpoints[0].chrom, v.breakpoints[0].start, v.breakpoints[1].start
        if v.name == "EVENT5":
            assert (chrom, start, end) == ("2", 321681, 321886)



def random_sequence(n):
    return "".join(list(numpy.random.choice(list("ACGT"),n)))


########################### DELETIONS ###########################

def get_deletion_variant(ref, pos, length):
    def validate(variant):
        assert isinstance(variant, Deletion)
        assert variant.deletionLength() == length
        assert variant.segments("ref")[1].start == pos, variant.segments("ref")

    result = ["2",
              pos, # output is one-based, position before deletion starts
              "deletion_1",
              ref[pos-1:pos+length],
              ref[pos-1],
              ".",
              "PASS",
              "SVTYPE=DEL",
              ".",
              "."]

    return "\t".join(map(str, result))+"\n", validate

@pytest.fixture(params=[
    ("chr1", 300, 25),
    ("chr1", 600, 120),
    ("chr2", 3000, 300),
])
def deletion(genome_source, request):
    chrom, pos, length = request.param
    refseq = genome_source.names_to_contigs[chrom]
    deletion, validate = get_deletion_variant(refseq, pos, pos+length)

    return deletion, validate

def test_deletions(tmpdir, genome_source, deletion):
    do_variant_test(tmpdir, genome_source, deletion)


########################### INSERTIONS ###########################

def get_insertion_variant(ref, pos, length):
    def validate(variant):
        assert isinstance(variant, SequenceDefinedVariant)
        assert len(variant.segments("alt")[1]) == length+1

    import random
    r = random.Random()
    r.seed(551)
    insertion_sequence =  "".join(r.choice("ACGT") for i in range(length))

    result = ["2",
              pos, # output is one-based, position before deletion starts
              "deletion_1",
              ref[pos],
              ref[pos]+insertion_sequence,
              ".",
              "PASS",
              "SVTYPE=INS",
              ".",
              "."]

    return "\t".join(map(str, result))+"\n", validate

@pytest.fixture(params=[
    ("chr1", 300, 25),
    ("chr1", 600, 120),
    ("chr2", 3000, 300),
])
def insertion(genome_source, request):
    chrom, pos, length = request.param
    refseq = genome_source.names_to_contigs[chrom]
    insertion, validate = get_insertion_variant(refseq, pos, pos+length)

    return insertion, validate

def test_inserts(tmpdir, genome_source, insertion):
    do_variant_test(tmpdir, genome_source, insertion)    



def do_variant_test(tmpdir, genome_source, variant):
    variant_line, variant_validator = variant
    vcf_path = str(tmpdir.join("svs.vcf"))

    with open(vcf_path, "w") as vcf:
        vcf.write(VCFHEADER)
        vcf.write(variant_line)

    print(variant)
    datahub = DataHub()
    datahub.args = mock.Mock(variants=vcf_path)
    datahub.genome = genome_source
    datahub.align_distance = 1000

    print("")
    parser = VCFParser(datahub)
    for v in parser.get_variants():
        variant_validator(v)


# def test_sim():
#     ref = random_sequence(100)
#     for s, read in simulate_long_reads(ref, 5, 20):
#         print(ref)
#         print(" "*(s-1), read)

# if __name__ == '__main__':
#     test_sim()