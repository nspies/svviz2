import collections
import logging
import pyfaidx
import seqlib

from genosv.utility import misc
from genosv.remap import mapq
from genosv.remap import ssw_aligner
from genosv.remap.alignment import Alignment

logger = logging.getLogger(__name__)


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
            cur_chrom = match_chrom_format(locus.chrom, list(self.keys()))
            self._blacklist.append(misc.Locus(cur_chrom, locus.start, locus.end, locus.strand))

    # @profile
    def align(self, read):
        alns = []
        qualities = read.original_qualities()

        raw_alns = self.aligner.align(read.original_sequence(), hardclip=False)

        for aln in raw_alns:
            aln = Alignment(aln)
            if aln.is_reverse and qualities is not None:
                aln._read.query_qualities = qualities[::-1]
            else:
                aln._read.query_qualities = qualities

            aln.chrom = self.keys()[aln.reference_id]

            # if self.blacklist is not None:
                # print("....", aln.locus, self.blacklist, misc.overlaps(aln.locus, self.blacklist))
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