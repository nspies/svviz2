import logging
import sys
import time

import alignment
import commandline
import genotyping
import realignment
import pairedreaditer
import utilities
import variants

from datahub import DataHub


logger = logging.getLogger(__name__)


def get_datahub():
    args = commandline.parse_args(sys.argv[1:])

    datahub = DataHub()
    datahub.set_args(args)
    datahub.align_distance = 1000
    for sample_name, sample in datahub.samples.items():
        sample.search_distance = 2400
    datahub.realigner = realignment.Realigner(datahub.genome)

    return datahub

def get_read_batch(sample, datahub):
    if sample.single_ended:
        yield from get_reads_unpaired(sample, datahub)
    else:
        yield from get_read_pairs(sample, datahub)

def get_reads_unpaired(sample, datahub):
    logger.info("Loading more reads...")
    cur_reads = []
    search_regions = datahub.variant.search_regions(sample.search_distance)

    for region in search_regions:
        chrom, start, end = region.chrom, region.start, region.end
        for read in sample.bam.fetch(chrom, start, end):
            # if read.query_name != "m150105_192231_42177R_c100761782550000001823161607221526_s1_p0/138972/39862_46995":
            #     continue

            cur_reads.append(alignment.Alignment(read))
            if datahub.args.batch_size is not None and len(cur_reads) >= datahub.args.batch_size:
                yield cur_reads
                logger.info("Loading more reads...")
                cur_reads = []

    yield cur_reads



def get_read_pairs(sample, datahub):
    """
    get batches of read-pairs -- this allows us to exhaustively search for mates without
    keeping everyone in memory
    """
    logger.info("Loading more read pairs...")
    cur_read_pairs = []
    search_regions = datahub.variant.search_regions(sample.search_distance)
    paired_read_iter = pairedreaditer.PairedReadIter(sample.bam, search_regions)

    for read_pair in paired_read_iter:
        if read_pair.query_name == "HA2WPADXX:19:5:1605415:0":
            print("!"*30)
            print(read_pair.read1._read)
            print(read_pair.read2._read)
            print("!"*30)
            import genomesource
            g = genomesource.FastaGenomeSource("/Volumes/frida/nspies/data/bwa-hg19/hg19.fasta")
            realigner = realignment.Realigner(g)
            realigner.set_alt_genome_source(g)
            realigner.realign_pair(read_pair, sample.read_statistics)

        cur_read_pairs.append(read_pair)
        if datahub.args.batch_size is not None and len(cur_read_pairs) >= datahub.args.batch_size:
            yield cur_read_pairs
            logger.info("Loading more read pairs...")
            cur_read_pairs = []

    yield cur_read_pairs

    print("Reads with only N:", paired_read_iter.N_count)


def map_realign(batch, realigner, sample):
    if sample.single_ended:
        return map_realign_unpaired(batch, realigner)
    else:
        return map_realign_pairs(batch, realigner, sample)

def map_realign_pairs(batch, realigner, sample):
    results = []
    for pair in batch:
        realigned_pair = realigner.realign_pair(pair, sample.read_statistics)
        results.append(realigned_pair)

    return results

def map_realign_unpaired(batch, realigner):
    import tqdm
    results = []
    # for read in tqdm.tqdm(batch):
    for read in batch:
        realigned = realigner.realign_single_end(read)
        results.append(realigned)

    return results


def save_realignments(aln_sets, sample, datahub):
    #TODO: check if we're supposed to save based on arguments
    for aln_set in aln_sets:
        if aln_set.supports_allele != "amb":
            aln_set.supporting_aln.fix_flags()

        if aln_set.supports_allele == "ref":
            if sample.single_ended:
                sample.out_ref_bam.write(aln_set.supporting_aln._read)
            else:
                sample.out_ref_bam.write(aln_set.supporting_aln.read1)
                sample.out_ref_bam.write(aln_set.supporting_aln.read2)

        elif aln_set.supports_allele == "alt":
            if sample.single_ended:
                sample.out_alt_bam.write(aln_set.supporting_aln._read)
            else:
                sample.out_alt_bam.write(aln_set.supporting_aln.read1)
                sample.out_alt_bam.write(aln_set.supporting_aln.read2)




def genotype_variant(datahub):
    for sample_name, sample in datahub.samples.items():
        ref_count = 0
        alt_count = 0

        sample.set_bwa_params(datahub.realigner)

        for batch in get_read_batch(sample, datahub):
            if sample.single_ended:
                logger.info("Analyzing {} reads".format(len(batch)))
            else:
                logger.info("Analyzing {} read pairs".format(len(batch)))

            aln_sets = map_realign(batch, datahub.realigner, sample)

            cur_ref_count, cur_alt_count = genotyping.assign_reads_to_alleles(
                aln_sets,
                variants.get_breakpoints_on_local_reference(datahub.variant, "ref"),
                variants.get_breakpoints_on_local_reference(datahub.variant, "alt"))
            ref_count += cur_ref_count
            alt_count += cur_alt_count

            save_realignments(aln_sets, sample, datahub)

        print("REF:", ref_count, "ALT:", alt_count) 


def run(datahub):
    """ this runs the app on the provided datahub """
    for variant in datahub.get_variants():
        genotype_variant(datahub)
        afgjhkl()

def main():
    """ entry point from command line """
    logging.basicConfig(level=logging.DEBUG)
    datahub = get_datahub()
    run(datahub)



if __name__ == '__main__':
    main()