import logging
import sys
import time

from genosv.app import commandline
from genosv.app import variants
from genosv.io import getreads
from genosv.app.datahub import DataHub
from genosv.remap import genotyping

logger = logging.getLogger(__name__)


def get_datahub():
    args = commandline.parse_args(sys.argv[1:])

    datahub = DataHub()
    datahub.set_args(args)
    datahub.align_distance = 0
    for sample_name, sample in datahub.samples.items():
        logger.info("Search distance: {:,}bp".format(sample.search_distance))
        datahub.align_distance = max(datahub.align_distance, sample.align_distance)
    logger.info("Align distance: {:,}bp".format(sample.align_distance))

    return datahub


def map_realign(batch, realigner, sample):
    if sample.single_ended:
        return map_realign_unpaired(batch, realigner)
    else:
        return map_realign_pairs(batch, realigner, sample)

def map_realign_pairs(batch, datahub, sample):
    ref_genome_sources = [datahub.genome, datahub.local_ref_genome_source]
    alt_genome_sources = [datahub.local_alt_genome_source]

    for genome_source in ref_genome_sources+alt_genome_sources:
        genome_source.set_aligner_params(sample.sequencer)

    import tqdm
    # for pair in batch:
    for pair in tqdm.tqdm(batch):
        # if pair.name == "ST-E00130:359:HGV3HCCXX:1:1120:26098:65265":
        pair.realign(ref_genome_sources, alt_genome_sources)

    return batch

# def map_realign_pairs(batch, realigner, sample):
#     results = []
#     for pair in batch:
#         realigned_pair = realigner.realign_pair(pair, sample.read_statistics)
#         results.append(realigned_pair)

#     return results

# def map_realign_unpaired(batch, realigner):
#     import tqdm
#     results = []
#     # for read in tqdm.tqdm(batch):
#     for read in batch:
#         realigned = realigner.realign_single_end(read)
#         results.append(realigned)

#     return results


def save_realignments(aln_sets, sample, datahub):
    #TODO: check if we're supposed to save based on arguments
    for aln_set in aln_sets:
        if aln_set.supports_allele != "amb":
            aln_set.supporting_aln.fix_flags()

        if aln_set.supports_allele == "ref":
            if sample.single_ended:
                sample.out_ref_bam.write(aln_set.supporting_aln._read)
            else:
                sample.out_ref_bam.write(aln_set.supporting_aln.aln1._read)
                sample.out_ref_bam.write(aln_set.supporting_aln.aln2._read)

        elif aln_set.supports_allele == "alt":
            if sample.single_ended:
                sample.out_alt_bam.write(aln_set.supporting_aln._read)
            else:
                sample.out_alt_bam.write(aln_set.supporting_aln.aln1._read)
                sample.out_alt_bam.write(aln_set.supporting_aln.aln2._read)

    sample.out_alt_bam.close()
    sample.out_ref_bam.close()


def genotype_variant(datahub):
    temp_storage = {}

    for sample_name, sample in datahub.samples.items():
        ref_count = 0
        alt_count = 0

        # sample.set_bwa_params(datahub.realigner)
        temp_storage[sample_name] = []

        for batch in getreads.get_read_batch(sample, datahub):
            if sample.single_ended:
                logger.info("Analyzing {} reads".format(len(batch)))
            else:
                logger.info("Analyzing {} read pairs".format(len(batch)))

            aln_sets = map_realign(batch, datahub, sample)
            # aln_sets = map_realign(batch, datahub.realigner, sample)

            cur_ref_count, cur_alt_count = genotyping.assign_reads_to_alleles(
                aln_sets,
                variants.get_breakpoints_on_local_reference(datahub.variant, "ref"),
                variants.get_breakpoints_on_local_reference(datahub.variant, "alt"),
                sample.read_statistics)
            ref_count += cur_ref_count
            alt_count += cur_alt_count

            save_realignments(aln_sets, sample, datahub)

            temp_storage[sample_name].extend(aln_sets)

        print("REF:", ref_count, "ALT:", alt_count)

    return temp_storage

def visualize(datahub, temp_storage):
    from genosv.visualize import track
    from genosv.io import export
    for sample_name, sample in datahub.samples.items():
        sample.tracks = {}
        for allele in ["alt", "ref"]:
            cur_alns = [aln.supporting_aln for aln in temp_storage[sample_name] if aln.supports_allele==allele]
            sample.tracks[allele] = track.Track(
                datahub.variant.chrom_parts(allele), cur_alns, 3000, 4000, datahub.variant, allele, False, True)

        for allele in ["alt", "ref"]:
            axis = track.Axis(sample.tracks[allele].scale, datahub.variant, allele)
            datahub.alleleTracks[allele]["axis"] = axis

    e = export.TrackCompositor(datahub)
    
    with open("temp.pdf", "wb") as outf:
        d = export.convertSVG(e.render(), "pdf", "webkittopdf")
        print(d[:100])
        outf.write(d)


def run(datahub):
    """ this runs the app on the provided datahub """
    for variant in datahub.get_variants():
        t0 = time.time()
        temp_storage = genotype_variant(datahub)
        t1 = time.time()
        print("TIME:::", t1-t0)
        visualize(datahub, temp_storage)

def main():
    """ entry point from command line """
    logging.basicConfig(level=logging.DEBUG)
    datahub = get_datahub()
    run(datahub)



if __name__ == '__main__':
    main()