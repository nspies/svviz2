import collections
import logging
# import numpy
import pyfaidx
import seqlib

from svviz2.utility import intervals, misc
# from svviz2.remap import mapq
from svviz2.remap import ssw_aligner
from svviz2.remap.alignment import Alignment


from svviz2.remap import _mapq

logger = logging.getLogger(__name__)



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

## These are pretty much identical to the versions in genomeview, but each method
## adds some functionality to deal with bwa/ssw alignment

class GenomeSource:
    def __init__(self, names_to_contigs, aligner_type="bwa"):
        self.names_to_contigs = collections.OrderedDict(names_to_contigs)
        self._bwa = None
        self._ssw = None
        self._blacklist = None

        self.aligner_type = aligner_type

    def get_seq(self, chrom, start, end, strand):
        seq = self.names_to_contigs[chrom][start:end+1]
        if strand == "-":
            seq = misc.reverse_comp(seq)
        return seq

    def keys(self):
        return list(self.names_to_contigs.keys())

    @property
    def blacklist(self):
        return self._blacklist

    @blacklist.setter
    def blacklist(self, blacklist_loci):
        self._blacklist = []

        for locus in blacklist_loci:
            cur_chrom = misc.match_chrom_format(locus.chrom, list(self.keys()))
            self._blacklist.append(intervals.Locus(cur_chrom, locus.start, locus.end, locus.strand))

    def align(self, read):
        alns = []
        qualities = read.original_qualities()

        raw_alns = self.aligner.align(read.original_sequence())

        for aln in raw_alns:
            aln = Alignment(aln)
            if aln.is_reverse and qualities is not None:
                aln._read.query_qualities = qualities[::-1]
            else:
                aln._read.query_qualities = qualities

            aln.chrom = self.keys()[aln.reference_id]
            aln._read.query_name = read.query_name

            # if self.blacklist is not None:
                # print("....", aln.locus, self.blacklist, misc.overlaps(aln.locus, self.blacklist))
            if self.blacklist is None or not intervals.overlaps(aln.locus, self.blacklist):
                aln.source = self
                aln.chrom = self.keys()[aln.reference_id]
                self.score_alignment(aln)

                aln.set_tag("mq", read.mapq)
                alns.append(aln)
    
        return alns

    def score_alignment(self, aln):
        # TODO: move the mapqcalculator code to here

        # mc = mapq.MAPQCalculator(self)
        # aln.score = mc.get_alignment_end_score(aln)

        ref_seq = self.get_seq(aln.chrom, aln.reference_start, aln.reference_end, "+").upper()
        aln.score = _mapq.get_alignment_end_score(aln._read, ref_seq)
        # aln.score = s2

        # assert numpy.isclose(aln.score, s2, rtol=1e-5), "{} :: {}".format(aln.score, s2)

    @property
    def aligner(self):
        if self.aligner_type == "bwa":
            return self.bwa
        elif self.aligner_type == "ssw":
            return self.ssw

    @property
    def ssw(self):
        if self._ssw is None:
            self._ssw = ssw_aligner.Aligner(self.names_to_contigs)
        return self._ssw

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
        if self.aligner_type != "bwa":
            print("not bwa... skipping setting aligner settings")
            return

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
    def __init__(self, path, aligner_type="bwa"):
        self.path = path
        self._fasta = None
        self._bwa = None
        self._blacklist = None
        self.aligner_type = aligner_type
        
    def get_seq(self, chrom, start, end, strand):
        chrom = misc.match_chrom_format(chrom, list(self.fasta.keys()))

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
            result = self._bwa.loadIndex(self.path)

            if not result:
                raise IOError("Failed to load bwa index from file {}".format(self.path))
            
            logger.info("Loading bwa index done.")

        return self._bwa

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_fasta"] = None
        state["_bwa"] = None
        return state
