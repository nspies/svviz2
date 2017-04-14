import collections
import pysam

import utilities


ATTRIBS = ["query_name",
           "flag",
           "reference_id",
           "reference_start",
           "mapping_quality",
           "cigarstring",
           "next_reference_id",
           "next_reference_start",
           "query_alignment_length",
           "query_sequence",
           "query_qualities",
           "tags"]

# def check_overlaps():
#   if check_overlaps:
#       for ref_mapping in ref_mappings:
#           ref_chrom = self.ref_mapper.ChrIDToName(ref_mapping.reference_id)
#           ref_mapping_locus = utilities.Locus(
#               ref_chrom, ref_mapping.reference_start,
#               ref_mapping.reference_start+ref_mapping.reference_length, "+")
#           overlaps = False
#           for breakpoint in self.breakpoints:
#               if (ref_mapping_locus.overlaps(breakpoint) or ref_mapping_locus.overlapsAntisense(breakpoint)):
#                   overlaps = True
#                   break
#           if overlaps:
#               ref_mapping.set_tag("OV", True)
#           else:
#               ref_mapping.set_tag("OV", False)

class Alignment(object):
    """
    this is a pickle-able wrapper for reads from a bam file
    """
    def __init__(self, read):
        self._read = read
        self._storage = None

    def original_sequence(self):
        if self.is_reverse:
            return utilities.reverse_comp(self.query_sequence)
        return self.query_sequence

    def original_qualities(self):
        if self.is_reverse:
            return self.query_qualities[::-1]
        return self.query_qualities

    @property
    def locus(self):
        chrom = "chrMissing"
        if self.source is not None:
            chrom = self.source.get_chrom_by_id(self._read.reference_id)
        start = self._read.reference_start
        end = self._read.reference_end
        locus = utilities.Locus(chrom, start, end, "+")

        return locus

    def fix_flags(self):
        if self.mapq is not None:
            self._read.mapq = int(self.mapq)

    def __getattr__(self, name):
        # kinda inelegant -- maybe we can fix this in the future?

        if name in ["_read", "_storage"]:
            return
        if self._read is None and self._storage is not None:
            self._unflatten()
        return getattr(self._read, name)

    def _flatten(self):
        self._storage = {}
        for name in ATTRIBS:
            self._storage[name] = getattr(self._read, name)
        self._read = None

    def _unflatten(self):
        if self._read is None:
            self._read = pysam.AlignedSegment()
            for name in ATTRIBS:
                setattr(self._read, name, self._storage[name])

    def __getstate__(self):
        self._flatten()
        state = self.__dict__.copy()

        return state


_orients = {False: "+", True:"-"}


class PairAlignment(object):
    """
    This is a concrete instance of a single possible pairing of read alignments
    """
    def __init__(self, read1, read2):
        assert read1.query_name == read2.query_name

        self.query_name = read1.query_name

        self.read1 = read1
        self.read2 = read2

        self._insert_size = None
        self._orientation = None

        self.score = None
        self.posterior = None
        self.mapq = None

        self.source = None

    @property
    def locus(self):
        assert self.read1.reference_id == self.read2.reference_id

        chrom = "chrMissing"
        if self.source is not None:
            chrom = self.source.get_chrom_by_id(self.read1.reference_id)
        start = min(self.read1.reference_start, self.read2.reference_start)
        end = max(self.read1.reference_end, self.read2.reference_end)
        locus = utilities.Locus(chrom, start, end, "+")

        return locus
        
    @property
    def insert_size(self):
        if self._insert_size is None:
            self._calculate_pairing_stats()
        return self._insert_size

    @property
    def orientation(self):
        if self._orientation is None:
            self._calculate_pairing_stats()
        return self._orientation
        
    def _calculate_pairing_stats(self):
        if self.read1.reference_id != self.read2.reference_id:
            return

        aln1, aln2 = self.read1, self.read2

        if aln1.reference_start > aln2.reference_start:
            aln1, aln2 = aln2, aln1

        orientation = (aln1.is_reverse, aln2.is_reverse)
        self._orientation = "".join(_orients[x] for x in orientation)

        self._insert_size = (aln2.reference_start + aln2.reference_length) - aln1.reference_start

    def concordant(self, read_stats):
        if self.read1.reference_id != self.read2.reference_id:
            return False

        if self.orientation not in read_stats.orientations:
            return False
        if self.insert_size is not None and self.insert_size > read_stats.maxInsertSize() * 1.5:
            return False

        return True

    def fix_flags(self):
    #     pair.read1.query_qualities = pysam.qualitystring_to_array("I"*pair.read1.query_length)#numpy.repeat(40, pair.read1.query_length)
    #     pair.read2.query_qualities = pysam.qualitystring_to_array("I"*pair.read2.query_length)#numpy.repeat(40, pair.read2.query_length)
        self.read1.is_paired = True
        self.read2.is_paired = True

        self.read1.is_read1 = True
        self.read2.is_read1 = False
        self.read1.is_read2 = False
        self.read2.is_read2 = True

        self.read1.mate_is_reverse = self.read2.is_reverse
        self.read2.mate_is_reverse = self.read1.is_reverse

        self.read1.next_reference_id = self.read2.reference_id
        self.read2.next_reference_id = self.read1.reference_id

        self.read1.next_reference_start = self.read2.reference_start
        self.read2.next_reference_start = self.read1.reference_start

        self.read1.template_length = self.insert_size
        self.read2.template_length = self.insert_size

        self.read1.is_secondary = False
        self.read1.is_supplementary = False

        self.read2.is_secondary = False
        self.read2.is_supplementary = False

        if self.mapq is not None:
            self.read1.mapq = int(self.mapq)
            self.read2.mapq = int(self.mapq)


class MultiReferenceAlignmentSet(object):
    """
    represents all the alignments of a read(-pair) against each allele
    """
    def __init__(self, ref_source, alt_source):
        self.name = None
        self.ref_source = ref_source
        self.alt_source = alt_source

        self.alignments = None

    def add_alignments(self, alignments, allele):
        if self.alignments is None:
            self.alignments = {"ref":[], "alt":[]}

        for aln in alignments:
            if not self.name:
                self.name = aln.query_name
            assert self.name == aln.query_name

        source = self.ref_source if allele == "ref" else self.alt_source
        for alignment in alignments:
            alignment.source = source

        self.alignments[allele].extend(alignments)

    def get_alns(self, allele=None):
        alleles = [allele]
        if allele is None:
            alleles = ["ref", "alt"]
        alns = []
        for allele in alleles:
            alns.extend(self.alignments[allele])

        return alns

class MultiReferenceReadEndAlignmentSet(MultiReferenceAlignmentSet):
    """
    represents all the alignments of the read ends against one genome
    this is used before read end alignments have been paired together
    """
    def add_alignments(self, alignments, allele, end=1):
        if self.alignments is None:
            self.alignments = {"ref":{}, "alt":{}}

        if not end in self.alignments[allele]:
            self.alignments[allele][end] = []

        for aln in alignments:
            if not self.name:
                self.name = aln.query_name
            assert self.name == aln.query_name

        self.alignments[allele][end].extend(alignments)

    def get_alns(self, allele=None, end=None):
        alleles = [allele]
        if allele is None:
            alleles = ["ref", "alt"]
        ends = [end]
        if end is None:
            ends = sorted(set(self.alignments["ref"].keys()).union(set(self.alignments["alt"].keys())))

        alns = []
        for allele in alleles:
            for end in ends:
                alns.extend(self.alignments[allele][end])

        return alns

    def info(self):
        return "END 1:: ref:{}  alt:{}    END 2::  ref:{}  alt:{}".format(
            len(self.get_alns("ref", 1)),len(self.get_alns("alt", 1)),
            len(self.get_alns("ref", 2)),len(self.get_alns("alt", 2)))

    # def extend(self, alns):
    #     for aln in alns:
    #         self.append(aln)

    # def append(self, aln):
    #     if self.name is not None:
    #         self.name = aln.query_name
    #     assert self.name == aln.query_name

    #     self.alignments.append(aln)




