import collections
import json
import numpy
import os
import pandas

def report(datahub):
    """
    sample  allele  key value
    """

    results = []

    results.extend(tally_support(datahub))

    result_df = pandas.DataFrame(results, columns=["sample", "allele", "key", "value"])

    print(result_df)

    report_filename = "{}.report.tsv".format(datahub.variant.short_name())
    report_path = os.path.join(datahub.args.outdir, report_filename)

    result_df.to_csv(report_path, sep="\t", index=False)

def tally_support(datahub):
    results = []

    for sample_name, sample in datahub.samples.items():
        for allele in ["alt", "ref"]:
            bam = sample.outbam(allele, "r")
            cur_results = _tally_support(bam)
            for cur_result in cur_results:
                results.append((sample_name, allele) + cur_result)

    return results

def _tally_support(bam):
    count = 0
    weighted_count = 0
    breakpoint_overlaps = collections.defaultdict(list)
    breakpoint_counts = collections.Counter()
    extensions = collections.defaultdict(list)

    for read in bam:
        if not read.is_paired or read.is_read1:
            count += 1
            weighted_count += read.mapq/40.

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

    for breakpoint, count in breakpoint_counts.items():
        breakpoint, overlaps_sequence = breakpoint
        key = "count_{}_{}".format(breakpoint, "seq" if overlaps_sequence else "pair")
        results.append((key, count))

    for breakpoint, cur_extensions in extensions.items():
        key = "extension_{}".format(breakpoint)
        results.append((key, numpy.mean(cur_extensions)))


    return results