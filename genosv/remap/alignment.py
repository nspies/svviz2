import collections
import pysam

from genosv.utility import misc
from genosv.remap import mapq


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

class Alignment(object):
    """
    this is a pickle-able wrapper for reads from a bam file
    """
    def __init__(self, read):
        self._read = read
        self._storage = None
        self.chrom = None

        self.ref_pairs = []
        self.alt_pairs = []


    def original_sequence(self):
        if self.is_reverse:
            return misc.reverse_comp(self.query_sequence)
        return self.query_sequence

    def original_qualities(self):
        if self.query_qualities is None:
            return None
        if self.is_reverse:
            return self.query_qualities[::-1]
        return self.query_qualities

    @property
    def locus(self):
        chrom = self.chrom
        if chrom is None:
            chrom = self._read.reference_name

        start = self._read.reference_start
        end = self._read.reference_end

        # if self.cigartuples[0][0] == 4:
        #     start += self.cigartuples[0][1]
        # if self.cigartuples[-1][0] == 4:
        #     end -= self.cigartuples[-1][1]

        locus = misc.Locus(chrom, start, end, "-" if self.is_reverse else "+")

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



    def realign_against_allele(self, genome_sources):
        """ use this to realign against one allele """
        alns = []
        for genome_source in genome_sources:
            # this calculates the alignment score
            cur_alns = genome_source.align(self)
            for aln in cur_alns:
                aln.concordant = lambda x: True
                aln.loci = [aln.locus]
                aln.name = self._read.query_name
            alns.extend(cur_alns)

        return alns

    def realign(self, ref_genome_sources, alt_genome_sources):
        self.ref_pairs = self.realign_against_allele(ref_genome_sources)
        self.alt_pairs = self.realign_against_allele(alt_genome_sources)

        mapq.set_mapqs(self.ref_pairs+self.alt_pairs)
        self.ref_pairs.sort(key=lambda x: x.score, reverse=True)
        self.alt_pairs.sort(key=lambda x: x.score, reverse=True)



_orients = {False: "+", True:"-"}

class AlignmentPair(object):
    def __init__(self, aln1, aln2, name=None):
        self.aln1 = aln1
        self.aln2 = aln2

        self._insert_size = None
        self._orientation = None

        self.score = None
        self.posterior = None
        self.mapq = None

        self.name = name

    @property
    def loci(self):
        if self.aln1.reference_id == self.aln2.reference_id:
            chrom = self.aln1.chrom
            start = min(self.aln1.reference_start, self.aln2.reference_start)
            end = max(self.aln1.reference_end, self.aln2.reference_end)
            locus = misc.Locus(chrom, start, end, "+")

            return [locus]
        else:
            return [self.aln1.locus, self.aln2.locus]
        
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
        if self.aln1.reference_id != self.aln2.reference_id:
            return

        aln1, aln2 = self.aln1, self.aln2

        if aln1.reference_start > aln2.reference_start:
            aln1, aln2 = aln2, aln1

        orientation = (aln1.is_reverse, aln2.is_reverse)
        self._orientation = "".join(_orients[x] for x in orientation)

        self._insert_size = (aln2.reference_start + aln2.reference_length) - aln1.reference_start

    def concordant(self, read_stats):
        if self.aln1.reference_id != self.aln2.reference_id:
            return False

        if self.orientation not in read_stats.orientations:
            return False
        if self.insert_size is not None and self.insert_size > read_stats.maxInsertSize() * 1.5:
            return False

        return True

    def fix_flags(self):
        self.aln1._read.query_name = self.name
        self.aln2._read.query_name = self.name
        
        self.aln1._read.is_paired = True
        self.aln2._read.is_paired = True

        self.aln1._read.is_read1 = True
        self.aln2._read.is_read1 = False
        self.aln1._read.is_read2 = False
        self.aln2._read.is_read2 = True

        self.aln1._read.mate_is_reverse = self.aln2.is_reverse
        self.aln2._read.mate_is_reverse = self.aln1.is_reverse

        self.aln1._read.next_reference_id = self.aln2.reference_id
        self.aln2._read.next_reference_id = self.aln1.reference_id

        self.aln1._read.next_reference_start = self.aln2.reference_start
        self.aln2._read.next_reference_start = self.aln1.reference_start

        if self.insert_size is not None:
            if self.aln1.reference_start < self.aln2.reference_start:
                self.aln1._read.template_length = self.insert_size
                self.aln2._read.template_length = -self.insert_size
            else:
                self.aln1._read.template_length = -self.insert_size
                self.aln2._read.template_length = self.insert_size

        self.aln1._read.is_secondary = False
        self.aln1._read.is_supplementary = False

        self.aln2._read.is_secondary = False
        self.aln2._read.is_supplementary = False

        if self.mapq is not None:
            self.aln1._read.mapq = int(self.mapq)
            self.aln2._read.mapq = int(self.mapq)




