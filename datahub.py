import collections
import logging
import os
import pysam
import sys

# from svviz import annotations
# from svviz import gff
import genomesource
from sample import Sample
import vcfparser

logger = logging.getLogger(__name__)


def name_from_bam_path(bampath):
    return os.path.basename(bampath).replace(".bam", "").replace(".sorted", "").replace(".sort", "").replace(".", "_").replace("+", "_")
# def name_from_bed_path(bampath):
#     return os.path.basename(bampath).replace(".bed", "").replace(".sorted", "").replace(".sort", "").replace(".", "_").replace("+", "_").replace(".gz", "")

def _get_bam_headers(variant, allele):
    seqs = variant.seqs(allele)
    header = {"HD":{"VN":1.3,"SO":"unsorted"}}
    sq = []
    for name in seqs:
        sq.append({"SN":name.replace("/", "__"), "LN":len(seqs[name])})
    header["SQ"] = sq
    return header


class DataHub(object):
    def __init__(self):
        self.args = None
        self.align_distance = 0
        self.samples = collections.OrderedDict()
        self.genome = None
        # self.annotationSets = collections.OrderedDict()


    def __getstate__(self):
        """ allows pickling of DataHub()s """
        state = self.__dict__.copy()
        # del state["args"]
        del state["genome"]
        return state

    def get_variants(self):
        vcf = vcfparser.VCFParser(self)
        for variant in vcf.get_variants():
            self.set_cur_variant(variant)
            yield variant

    def set_cur_variant(self, variant):
        self.variant = variant
        ref_genome_source = genomesource.GenomeSource(self.variant.seqs("ref"))
        self.realigner.set_ref_genome_source(ref_genome_source)

        alt_genome_source = genomesource.GenomeSource(self.variant.seqs("alt"))
        self.realigner.set_alt_genome_source(alt_genome_source)

        # TODO: fix this...
        for sample_name, sample in self.samples.items():
            sample.out_alt_bam = pysam.AlignmentFile("{}.alt.realigned.bam".format(sample_name), "wb",
                header=_get_bam_headers(self.variant, "alt"))
            sample.out_ref_bam = pysam.AlignmentFile("{}.ref.realigned.bam".format(sample_name), "wb",
                header=_get_bam_headers(self.variant, "ref"))

                # template=sample.bam)

        for allele in ["alt", "ref"]:
            with open("{}_genome.{}.fa".format(allele, variant.short_name()), "w") as genome_file:
                for name, seq in self.variant.seqs(allele).items():
                    genome_file.write(">{}\n{}\n".format(name.replace("/", "__"), seq))

        # with open("ref_genome.{}.fa".format(variant), "w") as ref_genome_file:
        #     for name, seq in self.variant.ref_seqs().items():
        #         ref_genome_file.write(">{}\n{}\n".format(name.replace("/", "__"), seq))


    def set_args(self, args):
        self.args = args

        self.genome = genomesource.FastaGenomeSource(args.ref)

        for bamPath in self.args.bam:
            name = name_from_bam_path(bamPath)

            # get unique name by appending _i as needed
            i = 0
            while name in self.samples:
                i += 1
                curname = "{}_{}".format(name, i)
                if curname not in self.samples:
                    name = curname
                    break

            sample = Sample(name, bamPath)
            self.samples[name] = sample

        # if self.args.annotations:
        #     for annoPath in self.args.annotations:
        #         name = nameFromBedPath(annoPath)
        #         if annoPath.endswith(".bed") or annoPath.endswith(".bed.gz"):
        #             self.annotationSets[name] = annotations.AnnotationSet(annoPath)
        #         else:
        #             if not (annoPath.endswith(".gff") or annoPath.endswith(".gff.gz") \
        #                 or annoPath.endswith(".gtf") or annoPath.endswith(".gtf.gz")):
        #                 logging.warn("Unknown annotation file extension; trying to parse as if GTF/GFF format: '{}'".format(annoPath))
        #             self.annotationSets[name] = gff.GeneAnnotationSet(annoPath)


    def __iter__(self):
        return iter(list(self.samples.values()))

