# svviz2

[![Build Status](https://travis-ci.org/nspies/svviz2.svg?branch=master)](https://travis-ci.org/nspies/svviz2)


This is a near complete rewrite of [svviz1](https://github.com/svviz/svviz). New features:

- uses bwa under the hood for realignments
  - substantial improvements in reliability and speed when realigning long reads
  - enables realignment against entire genome, identifying potential second-best hits
  - calculates a quantitative mapping quality score taking account of ref and alt hits genome-wide
  - uses weighted mapq scores to calculate evidence for ref and alt alleles, including genotype likelihoods
- substantially improved visualizations
  - "quick consensus" reduces background error rate in pacbio/nanopore and other long-read technologies
  - optionally uses tandem repeat finder (trf) to identify tandem repeats near candidate SV
- integrated dotplots
  - visualizes ref vs alt, allowing for visual identification of tandem repeats and other complex sequence
  - if bwa is being used for realignment, visualizes any second-best hit regions against candidate SV locus
  - if long-reads are provided as input, picks several long reads to plot as dotplots against ref and alt
