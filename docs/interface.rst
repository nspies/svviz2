Interface
=========

svviz2 requires three types of files as input:

- one or more sorted, indexed BAM (or CRAM) files containing read data to be assessed 
- a reference genome, including a fasta file that's been indexed using `bwa <https://github.com/lh3/bwa>`_
- a VCF file containing variants of interest. See below for more notes about this file.

**Note about CRAM files** svviz2 uses `pysam <http://pysam.readthedocs.io>`_ under the hood to access read data from BAM files. In theory pysam can read from CRAM files, but I have found its implementations to be somewhat less reliable. You're probably best off ensuring you have the latest version of pysam installed if you're working with CRAM files!


Representing Structural Variants in VCF files
---------------------------------------------

Please note that the input VCF file is parsed using htslib and so the formatting must conform quite closely to the `spec <https://samtools.github.io/hts-specs/>`_.

The following types of events are currently recognized by svviz2. Please submit an `issue <https://github.com/nspies/svviz2/issues>`_ if you would like to request support for additional event types or formats.

``Sequence-defined`` An arbitrary event can be specified by including the exact reference and alternate sequence in the ref and alt fields of the event. These two fields must thus include only the letters A, C, G and T. For example, a deletion of 6 nucleotides, replacing them with an insertion of 8 nucleotides would be specified as "CAGGTCA" (ref) and "CTACGAAGT" (alt) where the pos coordinate of the event points to the position of the C (upstream of the deletion and insertion). Note that any arbitrary length sequence can be placed in these fields, so an insertion of 8kb (for example, a LINE element) could include the full 8,000 base sequence of the inserted LINE in the alt field of the variant.

``Deletion`` Deletions can be specified either by including the exact reference and alternate sequence in the ref and alt columns of the VCF file, or by specifying ``SVTYPE=DEL;END=<end coordinate>`` in the INFO field of the variant record; please note that the end coordinate is the last genomic position that is deleted (ie, **inclusive**)

``Insertion`` Insertions may only be specified by including the exact reference and alternate sequence in the ref and alt fields. Only sequence-resolved insertions are supported!

``Breakend`` Any complex type of event such as a translocation can be represented as a pair of "breakends" which together specify the position and orientation of the two halves of the event. Please see the `VCF spec <https://samtools.github.io/hts-specs/VCFv4.3.pdf>`_ for a detailed description of BND-type complex structural variants.

``Inversion`` Similar to deletions, inversions can be specified by including ``SVTYPE=INV;END=<end coordinate>`` in the INFO field of a variant. Remember that end coordinates are inclusive.


Full command interface
----------------------

A brief summary of all of svviz2's arguments and options can be obtained by running ``svviz2`` without any arguments at the command line:

.. code-block:: none

    ssw library not found
    usage: svviz2 [options] --ref REF --variants VARIANTS BAM [BAM2 ...]

    svviz2 version 2.0a3

    optional arguments:
      -h, --help            show this help message and exit

    Required arguments:
      bam                   sorted, indexed bam file containing reads of interest to plot; can be
                            specified multiple times to load multiple samples
      --ref REF, -r REF     reference fasta file (a .faidx index file will be created if it doesn't
                            exist so you need write permissions for this directory)
      --variants VARIANTS, -V VARIANTS
                            the variants to analyze, in vcf or bcf format (vcf files may be
                            compressed with gzip)

    Optional arguments:
      --outdir OUTDIR, -o OUTDIR
                            output directory for visualizations, summaries, etc (default: current
                            working directory)
      --format FORMAT       format for output visualizations; must be one of pdf, png or svg
                            (default: pdf, or svg if no suitable converter is found)
      --savereads           output the read realignments against the appropriate alt or ref allele
                            (default: false)
      --min-mapq MIN_MAPQ   only reads with mapq>=MIN_MAPQ will be analyzed; when analyzing
                            paired-end data, at least one read end must be near the breakpoints
                            with this mapq (default:0)
      --align-distance ALIGN_DISTANCE
                            sequence upstream and downstream of breakpoints to include when
                            performing re-alignment (default: infer from data)
      --batch-size BATCH_SIZE
                            Number of reads to analyze at once; larger batch-size values may run
                            more quickly but will require more memory (default=10000)
      --downsample DOWNSAMPLE
                            Ensure the total number of reads per event per sample does not exceed
                            this number by downsampling (default: infinity)
      --aligner ALIGNER     The aligner to use for realigning reads; either ssw (smith-waterman) or
                            bwa (default=bwa)
      --only-realign-locally
                            Only when using bwa as the aligner backend, when this option is enabled,
                            reads will only be aligned locally around the breakpoints and not also
                            against the full reference genome (default: False)
      --fast                More aggressively skip reads that are unlikely to overlap
                            the breakpoints (default: false)
      --first-variant FIRST_VARIANT
                            Skip all variants before this variant; counting starts with first
                            variant in input VCF as 0 (default: 0)
      --last-variant LAST_VARIANT
                            Skip all variants after this variant; counting starts with first
                            variant in input VCF as 0 (default: end of vcf)
      --render-only
      --no-render
      --dotplots-only
      --no-dotplots
      --report-only
      --no-report
      --only-plot-context ONLY_PLOT_CONTEXT
                            Only show this many nucleotides before the first breakpoint, and the
                            last breakpoint in each region (default: show as much context as needed
                            to show all reads fully)
      --also-plot-context ALSO_PLOT_CONTEXT
                            Generates two plots per event, one using the default settings, and one
                            generatedby zooming in on the breakpoints as per the
                            --only-plot-context option
