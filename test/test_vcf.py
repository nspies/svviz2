import pytest
from unittest import mock

from genosv.app.datahub import DataHub
from genosv.io.vcfparser import VCFParser

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
VCFEXAMPLE="""
2	321681	bnd_W	G	G]17:198982]	6	PASS	SVTYPE=BND;MATEID=bnd_Y	GT	0/1
2	321682	bnd_V	T	]13:123456]T	6	PASS	SVTYPE=BND;MATEID=bnd_U	GT	0/1
13	123456	bnd_U	C	C[2:321682[	6	PASS	SVTYPE=BND;MATEID=bnd_V	GT	0/1
13	123457	bnd_X	A	[17:198983[A	6	PASS	SVTYPE=BND;MATEID=bnd_Z	GT	0/1
17	198982	bnd_Y	A	A]2:321681]	6	PASS	SVTYPE=BND;MATEID=bnd_W	GT	0/1
17	198983	bnd_Z	C	[13:123457[C	6	PASS	SVTYPE=BND;MATEID=bnd_X	GT	0/1
""".lstrip()

@pytest.fixture
def breakend_vcf(tmpdir):
    vcf_path = str(tmpdir.join("svs.vcf"))

    with open(vcf_path, "w") as vcf:
        vcf.write(VCFHEADER)
        vcf.write(VCFEXAMPLE)

    return vcf_path

def test_breakend_vcf(breakend_vcf):
    datahub = DataHub()
    datahub.args = mock.Mock(variants=breakend_vcf)

    parser = VCFParser(datahub)
    for v in parser.get_variants():
        print(v.name, v)