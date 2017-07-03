import logging
import numpy
import pysam
import sys

from genosv.io.readstatistics import ReadStatistics
from genosv.utility.bam import bam_sort_index
from genosv.utility.misc import str_to_bool

logger = logging.getLogger(__name__)



def _get_bam_headers(variant, allele):
    seqs = variant.seqs(allele)
    header = {"HD":{"VN":1.3,"SO":"unsorted"}}
    sq = []
    for name in seqs:
        sq.append({"SN":name.replace("/", "__"), "LN":len(seqs[name])})
    header["SQ"] = sq
    return header

def get_sequencer_from_bam_header(bam):
    sequencer = None
    for rg in bam.header.get("RG", []):
        if "PL" in rg:
            sequencer = rg["PL"].lower()
    return sequencer

class Sample(object):
    def __init__(self, name, bam_path, datahub, extra_args=None):
        self.name = name
        self.datahub = datahub

        self.single_ended = False
        self.orientations = None
        self._search_distance = None

        self.bam_path = bam_path
        self._bam = None

        self.outbams = {}
        self.outbam_paths = {}

        self.read_statistics = ReadStatistics(self.bam)
        if self.read_statistics.orientations == "any":
            self.single_ended = True

        self.sequencer = "illumina"
        if self.single_ended:
            mismatches = numpy.mean(self.read_statistics.number_mismatches)
            lengths = numpy.mean(self.read_statistics.readLengths)
            mismatch_rate = mismatches / lengths

            # these error rates aren't really accurate in terms of describing the 
            # two platforms anymore, but they do correspond to the presets that
            # bwa mem has, which we're mimicking
            if mismatch_rate > 0.10:
                self.sequencer = "nanopore"
            elif mismatch_rate > 0.01:
                self.sequencer = "pacbio"
            elif numpy.isnan(mismatches) and lengths > 1000:
                self.sequencer = "pacbio"

        sequencer = get_sequencer_from_bam_header(self.bam)
        if sequencer in ["illumina", "pacbio", "nanopore"]:
            self.sequencer = sequencer

        if "sequencer" in extra_args:
            self.sequencer = extra_args["sequencer"].lower()
            assert self.sequencer in ["illumina", "pacbio", "nanopore"]
            self.single_ended = True
        if "single_ended" in extra_args:
            self.single_ended = str_to_bool(extra_args["single_ended"])

        print("ALIGNMENT PARAMS:::", self.sequencer)


    @property
    def search_distance(self):
        if self._search_distance is None:
            if self.single_ended:
                longest_reads = numpy.percentile(self.read_statistics.readLengths, 99)
                if longest_reads > 1e5:
                    self._search_distance = 10000
                self._search_distance = 1000
            else:
                search_distance = numpy.percentile(self.read_statistics.insertSizes, 99)
                self._search_distance = int(search_distance)
        return self._search_distance

    @property
    def align_distance(self):
        if self.single_ended:
            longest_reads = numpy.percentile(self.read_statistics.readLengths, 95) * 1.5
            return int(longest_reads)
        else:
            longest_inserts = numpy.percentile(self.read_statistics.insertSizes, 98) * 1.5
            return int(longest_inserts)

    @property
    def bam(self):
        if self._bam is None:
            self._bam = pysam.AlignmentFile(self.bam_path, "rb")
            try:
                self._bam.fetch()
            except ValueError:
                logger.error("ERROR: Need to create index for input bam file: {}".format(self.bam_path))
                sys.exit(0)
        return self._bam

    def add_realignments(self, aln_sets):
        for allele in ["alt", "ref", "amb"]:
            self.outbam(allele, "w")

        for aln_set in aln_sets:
            # if aln_set.supports_allele != "amb":
            if aln_set.supporting_aln is None: continue
            aln_set.supporting_aln.fix_flags()

            outbam = self.outbam(aln_set.supports_allele, "w")
            if self.single_ended:
                outbam.write(aln_set.supporting_aln._read)
            else:
                outbam.write(aln_set.supporting_aln.aln1._read)
                outbam.write(aln_set.supporting_aln.aln2._read)

    def finish_writing_realignments(self):
        for allele in ["alt", "ref", "amb"]:
            self.outbams[allele].close()
            self.outbams.pop(allele)
            try:
                bam_sort_index(self.outbam_paths[allele])
            except:
                print("ERROR!"*30)


    def outbam(self, allele, mode):
        if mode == "w":
            if not allele in self.outbams:
                self.outbams[allele] = pysam.AlignmentFile(self.outbam_paths[allele], "wb",
                    header=_get_bam_headers(self.datahub.variant, allele))
            return self.outbams[allele]

        assert not allele in self.outbams, "forgot to close outbam before re-opening"

        return pysam.AlignmentFile(self.outbam_paths[allele].replace(".bam", ".sorted.bam"))


    def __getstate__(self):
        """ allows pickling of Samples()s """
        state = self.__dict__.copy()
        del state["bam"]
        return state

