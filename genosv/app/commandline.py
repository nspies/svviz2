import argparse
import sys

import genosv

def visualization_file_format(original_string):
    string = original_string.lower()
    if string not in ["pdf", "png", "svg"]:
        msg = "Invalid output file format: '{}'".format(string)
        raise argparse.ArgumentTypeError(msg)
    return string

def parse_args(input_args):
    parser = argparse.ArgumentParser(description="genosv version {}".format(genosv.__version__),
        usage="%(prog)s [options] [demo] --ref REF --variants VARIANTS BAM [BAM2 ...]",
        formatter_class=argparse.RawTextHelpFormatter)

    required_args = parser.add_argument_group("Required arguments")

    required_args.add_argument("bam", nargs="+", help=
        "sorted, indexed bam file containing reads of interest to plot; can be specified multiple \n"
        "times to load multiple samples")

    required_args.add_argument("--ref", "-r", help=
        "reference fasta file (a .faidx index file will be created if it doesn't exist so you need \n"
        "write permissions for this directory)")

    required_args.add_argument("--variants", "-V", help=
        "the variants to analyze, in vcf or bcf format (vcf files may be compressed with gzip)")


    optional_args = parser.add_argument_group("Optional arguments")

    optional_args.add_argument("--outdir", "-o", type=str, help=
        "output directory for visualizations, summaries, etc (default: current working directory)")

    optional_args.add_argument("--format", type=visualization_file_format, default="pdf", help=
        "format for output visualizations; must be one of pdf, png or svg (default: pdf)")

    optional_args.add_argument("--savereads", action="store_true", help=
        "output the read realignments against the appropriate alt or ref allele (default: false)")

    optional_args.add_argument("--align-distance", type=int, help=
        "sequence upstream and downstream of breakpoints to include when performing re-alignment"
        "(default: infer from data)")
    
    optional_args.add_argument("--batch-size", type=int, default=10000, help=
        "Number of reads to analyze at once; larger batch-size values may run more quickly \n"
        "but will require more memory (default=10000)")

    optional_args.add_argument("--aligner", type=str, default="bwa", help=
        "The aligner to use for realigning reads; either ssw (smith-waterman) or \n"
        "bwa (default=bwa)")

    optional_args.add_argument("--only-realign-locally", action="store_true", help=
        "Only when using bwa as the aligner backend, when this option is enabled,\n"
        "reads will only be aligned locally around the breakpoints and not also against\n"
        "the full reference genome (default: False)")

    optional_args.add_argument("--render-only", action="store_true", help=
        "")

    optional_args.add_argument("--only-plot-context", type=int, help=
        "Only show this many nucleotides before the first breakpoint, and the last breakpoint\n"
        "in each region (default: show as much context as needed to show all reads fully)")

    optional_args.add_argument("--also-plot-context", type=int, help=
        "Generates two plots per event, one using the default settings, and one generated \n"
        "by zooming in on the breakpoints as per the --only-plot-context option")

    if len(input_args)<1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(input_args)
    if args.only_plot_context and args.also_plot_context:
        parser.print_help()
        print("ERROR: you only want to use one of --only-plot-context or --also-plot-context")
        sys.exit(1)

    return args