import collections
import logging
import os
import pysam


from genosv.app import genomesource
from genosv.app.sample import Sample
from genosv.app import variants
from genosv.io import getreads
from genosv.io import vcfparser
from genosv.io import saverealignments
from genosv.remap import maprealign
from genosv.remap import genotyping
from genosv.utility import misc


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

        self.alleleTracks = collections.defaultdict(collections.OrderedDict)
        # self.annotationSets = collections.OrderedDict()

        self.aligner_type = "bwa"


    def genotype_cur_variant(self):
        temp_storage = {}

        for sample_name, sample in self.samples.items():
            ref_count = 0
            alt_count = 0

            # sample.set_bwa_params(datahub.realigner)
            temp_storage[sample_name] = []

            for batch in getreads.get_read_batch(sample, self):
                if sample.single_ended:
                    logger.info("Analyzing {} reads".format(len(batch)))
                else:
                    logger.info("Analyzing {} read pairs".format(len(batch)))

                aln_sets = maprealign.map_realign(batch, self, sample)
                # aln_sets = map_realign(batch, datahub.realigner, sample)

                cur_ref_count, cur_alt_count = genotyping.assign_reads_to_alleles(
                    aln_sets,
                    variants.get_breakpoints_on_local_reference(self.variant, "ref"),
                    variants.get_breakpoints_on_local_reference(self.variant, "alt"),
                    sample.read_statistics)
                ref_count += cur_ref_count
                alt_count += cur_alt_count

                if self.args.savereads:
                    saverealignments.save_realignments(aln_sets, sample, self)

                temp_storage[sample_name].extend(aln_sets)

            print("REF:", ref_count, "ALT:", alt_count)

        return temp_storage

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

        local_coords_in_full_genome = self.variant.search_regions(self.align_distance)
        self.genome.blacklist = local_coords_in_full_genome
        
        self.local_ref_genome_source = genomesource.GenomeSource(
            self.variant.seqs("ref"), aligner_type=self.aligner_type)
        self.local_alt_genome_source = genomesource.GenomeSource(
            self.variant.seqs("alt"), aligner_type=self.aligner_type)

        if self.args.savereads:
            for sample_name, sample in self.samples.items():
                sample.alt_bam_path = os.path.join(
                    self.args.outdir, "{}.{}.alt.bam".format(variant.short_name(), sample_name))
                sample.out_alt_bam = pysam.AlignmentFile(sample.alt_bam_path, "wb",
                    header=_get_bam_headers(self.variant, "alt"))

                sample.ref_bam_path = os.path.join(
                    self.args.outdir, "{}.{}.ref.bam".format(variant.short_name(), sample_name))
                sample.out_ref_bam = pysam.AlignmentFile(sample.ref_bam_path, "wb",
                    header=_get_bam_headers(self.variant, "ref"))

            for allele in ["alt", "ref"]:
                outpath = os.path.join(self.args.outdir, "{}.genome.{}.fa".format(variant.short_name(), allele))
                with open(outpath, "w") as genome_file:
                    for name, seq in self.variant.seqs(allele).items():
                        genome_file.write(">{}\n{}\n".format(name.replace("/", "__"), seq))



    def set_args(self, args):
        self.args = args

        self.aligner_type = args.aligner
        assert self.aligner_type in ["bwa", "ssw"]
        
        self.genome = genomesource.FastaGenomeSource(args.ref)

        if self.args.outdir is None:
            self.args.outdir = os.getcwd()
        misc.ensure_dir(self.args.outdir)

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

