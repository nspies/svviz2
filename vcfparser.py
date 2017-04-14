import logging
import pysam
import re

import utilities
import variants

logger = logging.getLogger(__name__)

class VCFParser(object):
    def __init__(self, datahub):
        self.datahub = datahub
        self.vcf = pysam.VariantFile(datahub.args.variants)

    def get_variants(self):
        breakends = {}
        for variant in self.vcf:
            if not "SVTYPE" in variant.info:
                print("variant does not appear to be a structural variant, skipping:{}".format(variant))
                continue

            sv_type = variant.info["SVTYPE"]

            if sv_type == "BND":
                mateid = variant.info["MATEID"]
                if "," in mateid:
                    NotImplementedError("we currently don't support ambiguous breakends")

                if mateid in breakends:
                    breakend = get_breakend(variant, breakends.pop(mateid), self.datahub)
                    if breakend is not None:
                        yield breakend
                else:
                    assert not variant.id in breakends
                    breakends[variant.id] = variant

        if len(breakends) > 0:
            logger.warn("found {} unpaired breakends".format(len(breakends)))

def get_breakend(first, second, datahub):
    # return "{}\n{}".format(_parse_breakend(first), _parse_breakend(second))
    return parse_breakend(first, second, datahub)



def _parse_breakend(record):
    ref = record.ref
    alt = record.alts[0]
    if not ("[" in alt or "]" in alt):
        return None
    
    orientation = None
    altre1 = "(\[|\])(\w*):(\w*)(\[|\])(.*)"
    match = re.match(altre1, alt)
    if match:
        dir1, other_chrom, other_pos, dir2, alt_seq = match.groups()
        assert dir1 == dir2, alt
        if alt_seq != ref:
            raise Exception("not yet implemented: complex event")
        chrom, pos = record.chrom, record.pos

        if dir1 == "]": # eg ]13:123456]T bnd_V
            orientation = "--"
        else:           # eg [17:198983[A bnd_X 
            orientation = "-+"
    else:
        altre2 = "(.*)(\[|\])(\w*):(\w*)(\[|\])"
        match = re.match(altre2, alt)
        if match:
            alt_seq, dir1, other_chrom, other_pos, dir2 = match.groups()
            assert dir1 == dir2, (dir1, dir2)
            if alt_seq != ref:
                raise Exception("not yet implemented: complex event")
            chrom, pos = record.chrom, record.pos

            if dir1 == "]": # eg G]17:198982 bnd_W
                orientation = "+-"
            else:           # eg C[2:321682[ bnd_U
                orientation = "++"

    if orientation is None:
        return None
    else:
        return chrom, pos, other_chrom, int(other_pos), orientation, alt

def parse_breakend(record1, record2, datahub):
    chrom, pos, other_chrom, other_pos, orientation, alt = _parse_breakend(record1)

    _chrom, _pos, _other_chrom, _other_pos, _orientation, _alt = _parse_breakend(record2)
    if not (_chrom == other_chrom and _pos == other_pos):
        print(chrom, pos, other_chrom, other_pos, orientation)
        print(_chrom, _pos, _other_chrom, _other_pos, _orientation)
        logger.error("Malformed VCF: breakends do not appear to match:\n{}\n{}".format(
            record1, record2))
        return None

    if chrom == other_chrom and abs(pos-other_pos) < datahub.align_distance*5:
        logger.error("Can't yet handle nearby breakends; skipping")
        return None

    # convert from 1-based to 0-based coordinates
    breakpoint1 = utilities.Locus(chrom, pos-1, pos-1, orientation[0])
    breakpoint2 = utilities.Locus(chrom, other_pos-1, other_pos-1, orientation[1])

    print(breakpoint1, breakpoint2)
    return variants.Breakend(breakpoint1, breakpoint2, datahub)



# import logging
# import pyfaidx

# from svviz import genomesource
# from svviz import utilities
# from svviz import variants

# class VCFParserError(Exception):
#     pass

# class VCFRecord(object):
#     def __init__(self, fields, info):
#         self.chrom = fields[0]
#         self.start = int(fields[1])

#         self.svtype = info["SVTYPE"].upper()
#         self.alt = fields[4]

#         self.fields = fields
#         self.info = info

#         self.end = self.start
#         if "END" in info:
#             self.end = int(info["END"])

#         if self.svtype in ["INV", "DEL"]:
#             if self.end - self.start <= 0:
#                 raise VCFParserError("Inversions and deletions need to have start < end")

#     def __str__(self):
#         return "{}::{}:{}-{},{}".format(self.svtype, self.chrom, self.start, self.end, self.end-self.start)


# class VCFFile(object):
#     def __init__(self, vcf_path_or_buffer, data_hub):
#         try:
#             self.vcf = open(vcf_path_or_buffer)
#         except TypeError:
#             self.vcf = vcf_path_or_buffer

#         self.data_hub = data_hub

#     def get_variants(self):
#         for line in self.vcf:
#             if line.startswith("#") or len(line) < 10:
#                 continue
            
#             sv = parseVCFLine(line)            

# def getVariants(dataHub):
#     vcfpath = dataHub.args.breakpoints[0]
#     try:
#         vcfFile = open(vcfpath)
#     except IOError:
#         raise Exception("Could not open vcf file:{}".format(vcfpath))

#     svs = []
#     for line in vcfFile:
#         if line.startswith("#") or len(line) < 10:
#             continue
#         sv = parseVCFLine(line, dataHub)
#         if sv is not None:
#             svs.append(sv)
#     return svs

# def getMobileElementFasta(dataHub):
#     if not "repeats" in dataHub.sources:
#         dataHub.sources["repeats"] = genomesource.GenomeSource(dataHub.args.fasta)
#         # dataHub.sources["repeats"] = pyfaidx.Fasta(dataHub.args.fasta, as_raw=True)
#     return dataHub.sources["repeats"]

# def parseInfo(infoString):
#     info = {}
#     for curInfo in infoString.split(";"):
#         value = True
#         if "=" in curInfo:
#             key, value = curInfo.split("=")
#         else:
#             key = curInfo
#         info[key.upper()] = value
#     return info

# def parseVCFLine(line, dataHub):
#     try:
#         fields = line.strip().split()

#         info = parseInfo(fields[7])

#         record = VCFRecord(fields, info)

#         if record.svtype == "INS":
#             return parseInsertion(record, dataHub)
#         elif record.svtype == "DEL":
#             return parseDeletion(record, dataHub)
#         elif record.svtype == "INV":
#             return parseInversion(record, dataHub)
#         elif record.svtype == "TRA":
#             return parseTranslocation(record, dataHub)
#         raise VCFParserError("Unsupported variant type:{}".format(record.svtype))
#     except Exception as e:
#         logging.error("\n== Failed to load variant: {} ==".format(e))
#         logging.error(str(line.strip()))
#         return None

# def parseDeletion(record, dataHub):
#     deletion = variants.Deletion.from_breakpoints(record.chrom, record.start, record.end, dataHub.alignDistance, dataHub.genome)
#     if dataHub.args.max_deletion_size and deletion.deletionLength() > dataHub.args.max_deletion_size:
#         deletion = variants.LargeDeletion.from_breakpoints(record.chrom, record.start, record.end, dataHub.alignDistance, dataHub.genome)
#         print("*"*1000, deletion)
#     return deletion

# def parseInversion(record, dataHub):
#     region = utilities.Locus(record.chrom, record.start, record.end, "+")
#     return variants.Inversion(region, dataHub.alignDistance, dataHub.genome)
    
# def parseInsertion(record, dataHub):
#     altchars = set(record.alt.upper())
#     breakpoint = utilities.Locus(record.chrom, record.start, record.end, "+")

#     if altchars <= set("ACGTN"):
#         insertSeq = record.alt.upper()

#         variant = variants.Insertion(breakpoint, insertSeq, dataHub.alignDistance, dataHub.genome)
#     elif "MEINFO" in record.info:
#         meinfo = record.info["MEINFO"].split(",")
#         meName = meinfo[0]
#         meStart = utilities.getListDefault(meinfo, 1, 0)
#         meEnd = utilities.getListDefault(meinfo, 2, 1e100)
#         meStrand = utilities.getListDefault(meinfo, 3, "+")

#         meLocus = utilities.Locus(meName, meStart, meEnd, meStrand)
#         variant = variants.MobileElementInsertion(breakpoint, meLocus, getMobileElementFasta(dataHub), 
#             dataHub.alignDistance, dataHub.genome)
#     else:
#         raise VCFParserError("Unknown insertion sequence")

#     return variant


# if __name__ == '__main__':
#     t = """
# #CHROM  POS   ID  REF ALT   QUAL  FILTER  INFO  FORMAT  NA00001
# 1 2827693   . CCGTGGATGCGGGGACCCGCATCCCCTCTCCCTTCACAGCTGAGTGACCCACATCCCCTCTCCCCTCGCA  C . PASS  SVTYPE=DEL;END=2827680;BKPTID=Pindel_LCS_D1099159;HOMLEN=1;HOMSEQ=C;SVLEN=-66 GT:GQ 1/1:13.9
# 2 321682    . T <DEL>   6 PASS    IMPRECISE;SVTYPE=DEL;END=321887;SVLEN=-105;CIPOS=-56,20;CIEND=-10,62  GT:GQ 0/1:12
# 2 14477084  . C <DEL:ME:ALU>  12  PASS  IMPRECISE;SVTYPE=DEL;END=14477381;SVLEN=-297;MEINFO=AluYa5,5,307,+;CIPOS=-22,18;CIEND=-12,32  GT:GQ 0/1:12
# 3 9425916   . C <INS:ME:L1> 23  PASS  IMPRECISE;SVTYPE=INS;END=9425916;SVLEN=6027;CIPOS=-16,22;MEINFO=L1HS,1,6025,- GT:GQ 1/1:15
# 3 9425916   . C GGATGCTGATCGTAGCTG 23  PASS  IMPRECISE;SVTYPE=INS;END=9425916;SVLEN=6027;CIPOS=-16,22 GT:GQ 1/1:15
# 3 9425916   . C <INS:ME:L1> 23  PASS  IMPRECISE;SVTYPE=INS;END=9425916;SVLEN=6027;CIPOS=-16,22 GT:GQ 1/1:15
# 3 12665100  . A <DUP>   14  PASS  IMPRECISE;SVTYPE=DUP;END=12686200;SVLEN=21100;CIPOS=-500,500;CIEND=-500,500   GT:GQ:CN:CNQ  ./.:0:3:16.2
# 4 18665128  . T <DUP:TANDEM>  11  PASS  IMPRECISE;SVTYPE=DUP;END=18665204;SVLEN=76;CIPOS=-10,10;CIEND=-10,10  GT:GQ:CN:CNQ  ./.:0:5:8.3
#     """
#     class MockDataHub(object):
#         pass
#     dh = MockDataHub()
#     dh.alignDistance = 1000
#     dh.genome = "(*genome(*"

#     for line in t.split("\n"):
#         if len(line) <= 10 or line.startswith("#"):
#             continue
#         print(parseVCFLine(line, dh))