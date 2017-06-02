import argparse
import sys

import genosv


def parse_args(input_args):
    parser = argparse.ArgumentParser(description="genosv version {}".format(genosv.__version__),
        usage="%(prog)s [options] [demo] --ref REF --variants VARIANTS BAM [BAM2 ...]",
        formatter_class=argparse.RawTextHelpFormatter)

    required_args = parser.add_argument_group("Required arguments")

    required_args.add_argument("bam", action="append", help=
        "sorted, indexed bam file containing reads of interest to plot; can be specified multiple \n"
        "times to load multiple samples")

    required_args.add_argument("--ref", "-r", help=
        "reference fasta file (a .faidx index file will be created if it doesn't exist so you need \n"
        "write permissions for this directory)")

    required_args.add_argument("--variants", "-V", help=
        "the variants to analyze, in vcf or bcf format (vcf files may be compressed with gzip)")


    optional_args = parser.add_argument_group("Optional arguments")

    optional_args.add_argument("--batch-size", type=int, default=10000, help=
        "Number of reads to analyze at once; larger batch-size values may run more quickly \n"
        "but will require more memory (default=10000)")

    optional_args.add_argument("--aligner", type=str, default="bwa", help=
        "The aligner to use for realigning reads; either ssw (smith-waterman) or \n"
        "bwa (default=bwa)")

    if len(input_args)<1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(input_args)

    return args