import collections
import json
import numpy
import os
import pandas

from genosv.app.variants import non_negative
from genosv.utility import statistics
from genosv.remap import genotyping


def report(datahub):
    results = []

    results.extend(tally_support(datahub))
    results.extend(tally_segments(datahub))
    results.extend(tally_nearby_polymorphisms(datahub))

    result_df = pandas.DataFrame(results, columns=["sample", "allele", "key", "value"])

    print(result_df)

    report_filename = "{}.report.tsv".format(datahub.variant.short_name())
    report_path = os.path.join(datahub.args.outdir, report_filename)

    result_df.to_csv(report_path, sep="\t", index=False)


# def get_genotype_likelihoods(datahub):
#     results = []

#     for sample_name, sample in datahub.samples.items():
#         allele_counts1 = {}
#         allele_counts2 = {}


#         for allele in ["alt", "ref"]:
#             bam = sample.outbam(allele, "r")
#             cur_count1 = 0
#             cur_count2 = 0
#             for read in bam:
#                 if not read.is_paired or read.is_read1:
#                     cur_count1 += read.mapq
#                     cur_count2 += read.mapq/40


#             allele_counts1[allele] = cur_count1
#             allele_counts2[allele] = cur_count2







def tally_support(datahub):
    results = []

    for sample_name, sample in datahub.samples.items():
        allele_count = {}
        allele_mapq_sum = {}
        allele_weighted_count = {}

        for allele in ["alt", "ref"]:
            bam = sample.outbam(allele, "r")
            cur_results, count, mapq_sum, weighted_count = _tally_support(bam)

            allele_count[allele] = count
            allele_mapq_sum[allele] = mapq_sum
            allele_weighted_count[allele] = weighted_count

            for cur_result in cur_results:
                results.append((sample_name, allele) + cur_result)

        for name, value in [("count", allele_count), ("mapq", allele_mapq_sum), ("weighted", allele_weighted_count)]:
            print(name, value)
            cur_genotype = genotyping.calculate_genotype_likelihoods(
                value["ref"], value["alt"])
            results.append((sample_name, "", "GL_{}".format(name), cur_genotype[0]))
            results.append((sample_name, "", "GQ_{}".format(name), cur_genotype[1]))

            gt = max(zip(cur_genotype[0], ["0/0", "0/1", "1/1"]), key=lambda x:x[0])[1]
            results.append((sample_name, "", "GT_{}".format(name), gt))

    return results

def tally_segments(datahub):
    """
    gets region-specific statistics (ie, per-segment)
    """

    results = []

    for sample_name, sample in datahub.samples.items():
        for allele in ["alt", "ref"]:
            bam = sample.outbam(allele, "r")
            # parts = datahub.variant.chrom_parts(allele)
            # if len(parts) > 1: continue
            # part = list(parts)[0]
            # segments = part.segments

            segment_overlaps = collections.defaultdict(list)
            # cur_offset = 0

            # for segment in segments:
            #     start, end = cur_offset, cur_offset+len(segment)

            for part, start, end, segment in iter_segments(datahub, allele):
                for read in bam.fetch(part.id, start, end):
                    overlap = read.get_overlap(start, end)
                    segment_overlaps[segment].append(overlap)
                # cur_offset += len(segment)

            for _,_,_,segment in iter_segments(datahub, allele):
                overlaps = segment_overlaps.get(segment, [0])
                results.append((sample_name, allele, "region_{}".format(segment), numpy.mean(overlaps)))

    return results



def _tally_support(bam):
    count = 0
    weighted_count = 0
    mapq_sum = 0
    breakpoint_overlaps = collections.defaultdict(list)
    breakpoint_counts = collections.Counter()
    extensions = collections.defaultdict(list)

    for read in bam:
        if not read.is_paired or read.is_read1:
            count += 1
            weighted_count += read.mapq/40.
            mapq_sum += (1 - statistics.phred_to_prob(read.mapq, 10))

            cur_breakpoint_overlaps = json.loads(read.get_tag("Ov"))
            for breakpoint, info in cur_breakpoint_overlaps.items():
                overlap, overlaps_sequence, extension = info
                breakpoint_overlaps[breakpoint].append(overlap)
                breakpoint_counts[(breakpoint, overlaps_sequence)] += 1
                extensions[breakpoint].append(extension)

    results = [("count", count),
               ("weighted_count", weighted_count)]

    for breakpoint, overlaps in breakpoint_overlaps.items():
        key = "overlap_{}".format(breakpoint)
        results.append((key, numpy.mean(overlaps)))

    for breakpoint, cur_count in breakpoint_counts.items():
        breakpoint, overlaps_sequence = breakpoint
        key = "count_{}_{}".format(breakpoint, "seq" if overlaps_sequence else "pair")
        results.append((key, cur_count))

    for breakpoint, cur_extensions in extensions.items():
        key = "extension_{}".format(breakpoint)
        results.append((key, numpy.mean(cur_extensions)))


    return results, count, mapq_sum, weighted_count


def iter_segments(datahub, allele):
    parts = datahub.variant.chrom_parts(allele)

    for part in parts:
        segments = part.segments

        segment_overlaps = collections.defaultdict(list)
        cur_offset = 0

        for segment in segments:
            start, end = cur_offset, cur_offset+len(segment)
            cur_offset += len(segment)
            yield part, start, end, segment


def tally_nearby_polymorphisms(datahub):
    cur_part = None

    distance = 100
    results = []

    for sample_name, sample in datahub.samples.items():
        for allele in ["alt", "ref"]:
            bam = sample.outbam(allele, "r")

            parts_to_segments = collections.defaultdict(list)
            for part, start, end, segment in iter_segments(datahub, allele):
                parts_to_segments[part.id].append((part, start, end, segment))

            for segments in parts_to_segments.values():
                for i, (part, start, end, segment) in enumerate(segments):
                    if i == 0:
                        cur_count = _tally_polymorphisms(bam, part, non_negative(end-distance), end)
                    elif i == len(segments)-1:
                        cur_count = _tally_polymorphisms(bam, part, start, start+distance)
                    else:
                        cur_count = _tally_polymorphisms(bam, part, start, end)

                    kinds = ["n_mismatch_low", "n_mismatch_high", "n_indel_low", "n_indel_high", "n_total"]
                    for kind, value in zip(kinds, cur_count):
                        key = "{}_{}".format(segment, kind)
                        results.append((sample_name, allele, key, value))

    return results

def _tally_polymorphisms(bam, part, start, end):
    n_mismatch_low = 0
    n_mismatch_high = 0

    n_indel_low = 0
    n_indel_high = 0

    denominator = 0

    region_seq = part.get_seq()

    for pileupcolumn in bam.pileup(part.id, start, end, truncate=True):
        # pileupcolumn.n
        cur_counts = collections.Counter()
        deletions = 0
        insertions = 0

        denominator += 1

        column_counts = 0
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_refskip:
                continue

            column_counts += 1

            if pileupread.is_del:
                deletions += 1
            else:
                nuc =  pileupread.alignment.query_sequence[pileupread.query_position]
                if nuc != "N":
                    if nuc == region_seq[pileupcolumn.pos]:
                        cur_counts["ref"] += 1
                    else:
                        cur_counts[nuc] += 1
            if pileupread.indel > 0:
                insertions += 1

        if column_counts < 15: continue

        total = sum(cur_counts.values())
        cur_counts.pop("ref", 0)
        if len(cur_counts) > 0:
            best_mismatch = max(cur_counts.values())
            if best_mismatch > total * 0.2:
                n_mismatch_low += 1
            if best_mismatch > total * 0.8:
                n_mismatch_high += 1

        if (insertions > total*0.2) or (deletions > total*0.2):
            n_indel_low += 1

        if (insertions > total*0.8) or (deletions > total*0.8):
            n_indel_high += 1

    return n_mismatch_low, n_mismatch_high, n_indel_low, n_indel_high, denominator