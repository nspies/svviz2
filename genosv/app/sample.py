import logging
import numpy
import pysam
import sys

from genosv.io.readstatistics import ReadStatistics

logger = logging.getLogger(__name__)




class Sample(object):
    def __init__(self, name, bam_path):
        self.name = name

        self.single_ended = False
        self.orientations = None
        self._search_distance = None

        self.bam_path = bam_path
        self._bam = None

        self.outbam = None

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

    # def set_bwa_params(self, realigner):
    #     for bwa in [realigner.ref_mapper, realigner.alt_mapper]:
    #         if self._sequencer == "illumina":
    #             set_illumina_params(bwa)
    #         elif self._sequencer == "pacbio":
    #             set_pacbio_params(bwa)
    #         elif self._sequencer == "minion":
    #             set_minion_params(bwa)


    def __getstate__(self):
        """ allows pickling of Samples()s """
        state = self.__dict__.copy()
        del state["bam"]
        return state

