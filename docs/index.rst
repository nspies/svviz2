:tocdepth: 2

================================================================
svviz2 - to evaluate and visualize candidate structural variants
================================================================

.. Table of Contents
.. -----------------

.. toctree::
   :maxdepth: 1
   :hidden:

   Overview <self>
   install
   Visualization <visualization>
   Automated analysis <report>
   details
   changelog
   license

svviz assesses evidence for a structural variant found in high-throughput sequencing data. svviz is free and open source, available at `<https://github.com/nspies/svviz2>`_. Please submit issues, questions or feature requests using the `github issue tracker <https://github.com/nspies/svviz2/issues>`_.

Some typical use cases, with links to relevant documentation:

- :ref:`visually assessing <visualization>` whether tricky candidate SVs are likely to be true events
- :ref:`automated large-scale evaluation of SVs <report>`; for example, taking candidate events found using Illumina and assessing evidence in PacBio or Nanopore data
- :ref:`genotyping structural variants <genotypes>`, potentially in new datasets


Quick start
-----------

svviz2 requires python 3.3 or later. To install the latest version of svviz from github, run:

``pip3 install -U git+git://github.com/nspies/svviz2.git``

Detailed instructions, including how to ensure that all prerequisites are installed, are available on the :ref:`installation` page.

svviz2 takes as input a reference genome (bwa indexed), a list of variants in VCF format, and assesses read evidence for each variant for one or more samples provided in bam format:

``svviz2 [options] --ref <REF> --variants <VARIANTS.vcf> <BAM> [<BAM2> ...]``

Several files will be created for each event, including a pdf (or svg) visualizing the read evidence from each sample, a machine-readable report with information about the strength of evidence (eg number of reads supporting variant or reference allele), and optionally dotplots used to assess repetitive structures near structural variants.

How it works
------------

svviz takes as input one or several bam files with sequencing data, a genome fasta file and information about a putative structural variant (that has been identified using other methods). From these inputs, svviz uses a realignment process to identify reads supporting the reference allele, reads supporting the structural variant (or alternate allele), and reads that are not informative one way or the other (ambiguous).

Reads are then plotted, as in a standard genome browser, but along the sequence of either the alternate allele or the reference allele. The user can thus assess the support for the putative structural variant visually, determine if the breakpoints appear accurate, or estimate a likely genotype for each sample at the given structural variant.

svviz differs from existing genome browsers in being able to display arbitrary types of structural variants, such as large insertions, deletions, inversions and translocations. Genome browsers such as IGV are poorly suited to displaying the read data supporting these types of structural variants, which can differ dramatically from the reference genome sequence. 


New Features in svviz 2.0
-------------------------

- uses `bwa mem <https://github.com/lh3/bwa>`_ under the hood for realignments

  - substantial improvements in reliability and speed when realigning long reads
  - enables realignment against entire genome, identifying potential second-best hits
  - calculates a quantitative mapping quality score taking account of ref and alt hits genome-wide
  - uses weighted mapq scores to calculate evidence for ref and alt alleles, including genotype likelihoods

- substantially improved visualizations

  - "quick consensus" reduces background error rate in pacbio/nanopore and other long-read technologies
  - optionally uses tandem repeat finder (trf) to identify tandem repeats near candidate SV
  - visualization engine has been refactored into a separate `genomeview <https://github.com/nspies/genomeview>`_ module, facilitating future improvements

- integrated dotplots

  - visualizes ref vs alt, allowing for visual identification of tandem repeats and other complex sequence
  - if bwa is being used for realignment, visualizes any second-best hit regions against candidate SV locus
  - if long-reads are provided as input, picks several long reads to plot as dotplots against ref and alt


Citation
--------

svviz has been `published in Bioinformatics <http://dx.doi.org/10.1093/bioinformatics/btv478>`_. If you found this software useful for your research, please cite svviz as follows:

Spies N, Zook JM, Salit M, Sidow A. 2015. svviz: a read viewer for validating structural variants. Bioinformatics doi:bioinformatics/btv478.

svviz was developed by `Noah Spies <http://noahspies.org/>`_, a member of the `Sidow lab <http://www.sidowlab.org>`_ at Stanford and part of the `Joint Initiative for Metrology in Biology (JIMB) <http://jimb.stanford.edu/>`_. Funding was provided by the `National Institute of Standards and Technology (NIST) <http://www.nist.gov>`_.