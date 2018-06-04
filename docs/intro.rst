Introduction to svviz
=====================

The program takes as input one or several bam files with sequencing data, a genome fasta file, information about a putative structural variant (that has been identified using other methods), and (optionally) a bed file with genomic annotations. From these inputs, ``svviz`` uses a realignment process to identify reads supporting the reference allele, reads supporting the structural variant (or alternate allele), and reads that are not informative one way or the other (ambiguous).

Reads are then plotted, as in a standard genome browser, but along the sequence of either the alternate allele or the reference allele. The user can thus assess the support for the putative structural variant visually, determine if the breakpoints appear accurate, or estimate a likely genotype for each sample at the given structural variant.

``svviz`` differs from existing genome browsers in being able to display arbitrary types of structural variants, such as large insertions, deletions, inversions and translocations. Genome browsers such as IGV are poorly suited to displaying the read data supporting these types of structural variants, which can differ dramatically from the reference genome sequence. 


New Features in svviz 2.0
-------------------------

- uses `bwa mem`<https://github.com/lh3/bwa>_ under the hood for realignments
  - substantial improvements in reliability and speed when realigning long reads
  - enables realignment against entire genome, identifying potential second-best hits
  - calculates a quantitative mapping quality score taking account of ref and alt hits genome-wide
  - uses weighted mapq scores to calculate evidence for ref and alt alleles, including genotype likelihoods
- substantially improved visualizations
  - "quick consensus" reduces background error rate in pacbio/nanopore and other long-read technologies
  - optionally uses tandem repeat finder (trf) to identify tandem repeats near candidate SV
  - visualization engine has been refactored into a separate `genomeview`<https://github.com/nspies/genomeview>_ module, facilitating future improvements
- integrated dotplots
  - visualizes ref vs alt, allowing for visual identification of tandem repeats and other complex sequence
  - if bwa is being used for realignment, visualizes any second-best hit regions against candidate SV locus
  - if long-reads are provided as input, picks several long reads to plot as dotplots against ref and alt