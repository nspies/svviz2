#export DYLD_LIBRARY_PATH=/usr/local/lib/python3.6/site-packages/pysam/
import collections
import logging
import pyfaidx

import seqlib
import utilities

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


class GenomeSource(object):
    """
    this is a uniform interface to a reference genome, either on disk eg a FastaGenomeSource,
    or in memory (this GenomeSource class)
    """
    def __init__(self, names_to_contigs):
        self.names_to_contigs = collections.OrderedDict(names_to_contigs)

    def get_chrom_by_id(self, chrom_id):
        return list(self.names_to_contigs.keys())[chrom_id]

    def get_seq(self, chrom, start, end, strand):
        seq = self.names_to_contigs[chrom][start:end+1]
        if strand == "-":
            seq = utilities.reverse_comp(seq)
        return seq

    def get_bwa(self):
        """
        pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
        ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
        """
        bwa = seqlib.BWAWrapper()
        bwa.makeIndex(self.names_to_contigs)

        return bwa


class FastaGenomeSource(GenomeSource):
    """ pickle-able wrapper for pyfaidx.Fasta """
    def __init__(self, path):
        self.path = path
        self._fasta = None

    def get_chrom_by_id(self, chrom_id):
        return list(self.fasta.keys())[chrom_id]
        
    def get_seq(self, chrom, start, end, strand):
        chrom = match_chrom_format(chrom, list(self.fasta.keys()))

        seq = self.fasta[chrom][start:end+1]
        if strand == "-":
            seq = utilities.reverse_comp(seq)
        return seq

    @property
    def fasta(self):
        if self._fasta is None:
            self._fasta = pyfaidx.Fasta(self.path, as_raw=True)
        return self._fasta

    def get_bwa(self):
        logger.info("Loading bwa index...")
        
        bwa = seqlib.BWAWrapper()
        bwa.loadIndex(self.path)
        
        logger.info("Loading bwa index done.")

        return bwa

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_fasta"] = None
        return state

def make_genome_source(ref_path=None, ref_seqs=None):
    if ref_path is not None:
        return FastaGenomeSource(ref_path)
    else:
        return GenomeSource(ref_seqs)