import collections
import logging
import pyfaidx
import pysam

import seqlib
from genosv.utility import misc
from genosv.remap import mapq

logger = logging.getLogger(__name__)

class ReadPair(object):
    """ Represents a read-pair

    only one ReadPair exists per read-pair, but many Alignments and AlignmentPairs can exist
    per read-pair
    """
    def __init__(self, original_read1, original_read2, read_stats):
        self.original_read_ends = {
            "1": original_read1,
            "2": original_read2,
        }

        assert original_read1.query_name == original_read2.query_name
        self.name = original_read1.query_name

        self.ref_pairs = []
        self.alt_pairs = []

        # self.chosen_alignment_pair = None

        self.read_stats = read_stats

    def realign_against_allele(self, genome_sources):
        """ use this to realign against one allele """
        end_alns = collections.defaultdict(list)
        for end, original_read in self.original_read_ends.items():
            for genome_source in genome_sources:
                # this calculates the alignment score
                cur_alns = genome_source.align(original_read)
                end_alns[end].extend(cur_alns)

        pair_alns = []
        for aln1 in end_alns["1"]:
            for aln2 in end_alns["2"]:
                pair = AlignmentPair(aln1, aln2, name=self.name)
                self.read_stats.score_read_pair(pair)
                pair_alns.append(pair)

        return pair_alns

    def realign(self, ref_genome_sources, alt_genome_sources):
        self.ref_pairs = self.realign_against_allele(ref_genome_sources)
        self.alt_pairs = self.realign_against_allele(alt_genome_sources)

        mapq.set_mapqs(self.ref_pairs+self.alt_pairs)
        self.ref_pairs.sort(key=lambda x: x.score, reverse=True)
        self.alt_pairs.sort(key=lambda x: x.score, reverse=True)

        
        # print(self.name)
        # print("REF")
        # for pair in self.ref_pairs:
        #     print(pair.loci, pair.score, pair.aln1.cigarstring, pair.aln2.cigarstring, 
        #         pair.aln1.score, pair.aln2.score,
        #         pair.concordant(self.read_stats),
        #         pair.mapq)
        # print("ALT")
        # for pair in self.alt_pairs:
        #     print(pair.loci, pair.score, pair.aln1.cigarstring, pair.aln2.cigarstring, 
        #         pair.aln1.score, pair.aln2.score,
        #         pair.concordant(self.read_stats),
        #         pair.mapq)

        # print("///")


        # if len(self.ref_pairs)+len(self.alt_pairs) > 0:
        #     self.chosen_alignment_pair = max(self.ref_pairs+self.alt_pairs, 
        #                                      key=lambda x: x.score)


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




def match_chrom_format(chrom, keys):
    if chrom in keys:
        return chrom
    if "chr" in chrom:
        chrom2 = chrom.replace("chr", "")
    else:
        chrom2 = "chr{}".format(chrom)

    if chrom2 in keys:
        return chrom2
    return chrom


PARAMS = {
    "illumina":{
        "min_seed_length":19,
        "min_chain_weight":0,
        "gap_open":6,
        "gap_extension":1,
        "mismatch_penalty":4,
        "3prime_clipping_penalty":5,
        "5prime_clipping_penalty":5,
        "reseed_trigger":1.5,
        },
    "pacbio": {
        "min_seed_length":17,
        "min_chain_weight":40,
        "gap_open":1,
        "gap_extension":1,
        "mismatch_penalty":1,
        "3prime_clipping_penalty":0,
        "5prime_clipping_penalty":0,
        "reseed_trigger":10,
    },
    "nanopore": {
        "min_seed_length":14,
        "min_chain_weight":20,
        "gap_open":1,
        "gap_extension":1,
        "mismatch_penalty":1,
        "3prime_clipping_penalty":0,
        "5prime_clipping_penalty":0,
        "reseed_trigger":10,
    }
}

class GenomeSource(object):
    def __init__(self, names_to_contigs, blacklist=None):
        self.names_to_contigs = collections.OrderedDict(names_to_contigs)
        self._bwa = None
        self.blacklist = blacklist

    def get_seq(self, chrom, start, end, strand):
        seq = self.names_to_contigs[chrom][start:end+1]
        if strand == "-":
            seq = misc.reverse_comp(seq)
        return seq

    def keys(self):
        return list(self.names_to_contigs.keys())


    def align(self, read):
        alns = []
        qualities = read.original_qualities()

        for aln in self.bwa.align(read.original_sequence(), hardclip=False):
            aln = Alignment(aln)
            if aln.is_reverse and qualities is not None:
                aln._read.query_qualities = qualities[::-1]
            else:
                aln._read.query_qualities = qualities

            aln.chrom = self.keys()[aln.reference_id]

            if self.blacklist is None or not misc.overlaps(aln.locus, self.blacklist):
                aln.source = self
                aln.chrom = self.keys()[aln.reference_id]
                self.score_alignment(aln)

                alns.append(aln)
    
        return alns

    def score_alignment(self, aln):
        # TODO: move the mapqcalculator code to here

        mc = mapq.MAPQCalculator(self)
        aln.score = mc.get_alignment_end_score(aln)

    @property
    def bwa(self):
        """
        pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
        ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
        """
        if self._bwa is None:
            self._bwa = seqlib.BWAWrapper()
            self._bwa.makeIndex(self.names_to_contigs)

        return self._bwa

    def set_aligner_params(self, sequencer):
        params = PARAMS[sequencer]
        if "min_seed_length" in params:
            self.bwa.SetMinSeedLength(params["min_seed_length"])
        if "min_chain_weight" in params:
            self.bwa.SetMinChainWeight(params["min_chain_weight"])

        if "mismatch_penalty" in params:
            self.bwa.SetMismatchPenalty(params["mismatch_penalty"])

        if "gap_open" in params:
            self.bwa.SetGapOpen(params["gap_open"])
        if "gap_extension" in params:
            self.bwa.SetGapExtension(params["gap_extension"])

        if "3prime_clipping_penalty" in params:
            self.bwa.Set3primeClippingPenalty(params["3prime_clipping_penalty"])
        if "5prime_clipping_penalty" in params:
            self.bwa.Set5primeClippingPenalty(params["5prime_clipping_penalty"])

        if "reseed_trigger" in params:
            self.bwa.SetReseedTrigger(params["reseed_trigger"])


    def __getstate__(self):
        state = self.__dict__.copy()
        state["_bwa"] = None
        return state


class FastaGenomeSource(GenomeSource):
    """ pickle-able wrapper for pyfaidx.Fasta """
    def __init__(self, path, blacklist=None):
        self.path = path
        self._fasta = None
        self._bwa = None
        self.blacklist = blacklist
        
    def get_seq(self, chrom, start, end, strand):
        chrom = match_chrom_format(chrom, list(self.fasta.keys()))

        seq = self.fasta[chrom][start:end+1]
        if strand == "-":
            seq = misc.reverse_comp(seq)
        return seq

    def keys(self):
        return list(self.fasta.keys())
        
    @property
    def fasta(self):
        if self._fasta is None:
            self._fasta = pyfaidx.Fasta(self.path, as_raw=True)
        return self._fasta

    @property
    def bwa(self):
        if self._bwa is None:
            logger.info("Loading bwa index from file {}...".format(self.path))
            
            self._bwa = seqlib.BWAWrapper()
            self._bwa.loadIndex(self.path)
            
            logger.info("Loading bwa index done.")

        return self._bwa

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_fasta"] = None
        state["_bwa"] = None
        return state