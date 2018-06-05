:tocdepth: 2

================================================================
svviz2 - to evaluate and visualize candidate structural variants
================================================================


.. toctree::
   :maxdepth: 1

   self
   intro
   install
   report
   details
   changelog
   license



svviz assesses evidence for a structural variant found in high-throughput sequencing data. svviz is free and open source, available at `<https://github.com/nspies/svviz2>`_. Please submit issues, questions or feature requests using the `github issue tracker <https://github.com/nspies/svviz2/issues>`_.

Some typical use cases:

- visually assessing whether tricky candidate SVs are likely to be true events
- large-scale evaluation of SVs; for example, taking candidate events found using Illumina and assessing evidence in PacBio or Nanopore data
- genotyping events


Quick start
-----------

svviz2 requires python 3.3 or later. To install the latest version of svviz from github, run:

``pip3 install -U git+git://github.com/nspies/svviz2.git``

Detailed instructions, including how to ensure that all prerequisites are installed, are available on the :ref:`installation` page.

svviz2 takes as input a reference genome (bwa indexed), a list of variants in VCF format, and assesses read evidence for each variant for one or more samples provided in bam format:

``svviz2 [options] --ref REF --variants VARIANTS.vcf BAM [BAM2 ...]``

Several files will be created for each event, including a pdf (or svg) visualizing the read evidence from each sample, a machine-readable report with information about the strength of evidence (eg number of reads supporting variant or reference allele), and optionally dotplots used to assess repetitive structures near structural variants.

Citation
--------

svviz has been `published in Bioinformatics <http://dx.doi.org/10.1093/bioinformatics/btv478>`_. If you found svviz useful for your research, please cite svviz as follows:

Spies N, Zook JM, Salit M, Sidow A. 2015. svviz: a read viewer for validating structural variants. Bioinformatics doi:bioinformatics/btv478.

svviz was developed by `Noah Spies <http://stanford.edu/~nspies/>`_, a member of the `Sidow lab <http://mendel.stanford.edu/SidowLab/index.html>`_ at Stanford and part of the `Joint Initiative for Metrology in Biology (JIMB) <http://jimb.stanford.edu/>`_. Funding was provided by the `National Institute of Standards and Technology (NIST) <http://www.nist.gov>`_.