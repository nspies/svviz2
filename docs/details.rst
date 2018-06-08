Low-level details
=================

Mapping quality calculations
----------------------------

Each read analyzed by svviz2 is assigned a mapping quality (mapq) indicating the support for the indicated allele. For example, a read strongly supporting the alt allele might receive a mapq of 35 whereas a read that is ambiguous might receive a mapq of 0.

First, an alignment score is calculated according to the methods described in Zheng and Grice (PLoS Comput Biol, 2016). Briefly, each position that is aligned correctly (ie the bases match) adds probability proportional to 1-baseQ, ie proportional to the probability of a sequencing error at that position. Each position that is in a mismatch, indel or is clipped adds probability proportional to baseQ. A fixed penalty is added for each type of gap-open/gap-extend/clip. For paired-end reads, the resulting probabilities are multiplied together with the empirical insert size probability.

The alignment scores are then used to calculate a mapq, which is equal to the probability of a given alignment divided by the sum of the probabilities of all alignments. A read may align to one or more locations in the reference region (ie close to the breakpoints), the alt region, or elsewhere in the genome. All of these mappings are considered when calculating mapq.

.. _weighted_mapq:

**A note about mapqs**

In our experience, the mapq values are poorly calibrated -- that is, the actual values are not really proportional to the true probabilities. For example, a mapq of 40 doesn't really mean that a read is 99.99% likely to have derived from the given genomic locus. There are many possible reasons for this: the simple mapq model doesn't accurately model novel non-reference sequence; the scores for indels aren't normalized correctly to the mismatch scores; the base sequencing probabilities are incorrect; etc. In any case, we have found empirically that a "weighted mapq" score typically captures these additional uncertainties better than the raw mapq.

The weighted mapq is calculated as a proportion: :math:`Q/40`, where :math:`Q` is the (Phred-scaled) mapq value, and 40 is the maximum mapq value.


.. _genotypes:

Genotype likelihood and quality scores
--------------------------------------

Please note that while svviz2 makes some calculations for genotypes, the underlying model is extremely simplified and doesn't take into account the complexities and heterogeneity of structural variants. In particular, given that the alleles can differ substantially in size, the expected ratios of reads supporting the ref and alt alleles are rarely 1:1 (as is true for SNPs, but which is expected by the model). For example, paired-end reads that span a 10bp deletion (ie read1 ends upstream of event and read2 starts downstream) are difficult to assign with any probability to one allele or another. However, as the length of the deletion begins to approach the average insert size it becomes much more obvious which allele the read derives from. All that said, I'm not sure it's particularly worse than anything else out there!

Genotype likelihoods are calculated using a simple Bayesian formula similar to that employed by the GATK (Genome Analysis Toolkit, from the Broad Institute) for genotyping single nucleotide variants. Briefly, a prior is defined as:

.. math::

    S(g) = 
    \begin{cases}
        0.05    ,& g=\text{hom ref} \\
        0.5     ,& g=\text{het} \\
        0.95    ,& g=\text{hom alt}
    \end{cases}


Then, the probability of observing :math:`A` reads supporting the alternate allele and :math:`R` reads supporting the reference allele, given genotype :math:`g` is:

.. math::

    P(A,R | g) = {{A+R}\choose{A}} S(g)^A \big(1-S(g)\big)^R

and the resulting likelihood of genotype :math:`g` is thus calculated as:

.. math::

    P(g|A,R) = \frac{P(A,R|g) P(g)}{\sum_g P(A,R|g)P(g)}

svviz2 reports two values: the genotype likelihoods (GL), which are the log10-transformed numerator values from the above formula, and the genotype qualities, which are the Phred-scaled (:math:`Q = -10\log_{10} P`) likelihoods.

GL and GQ values are reported as triplets: (homozygous-reference, heterozygous, homozygous-alternate). The largest of these GQ values is the most likely genotype, and the difference in this value and the second-highest value can be used to determine confidence in the given genotype.

**An important note**

Genotype likelihoods and qualities assume the correctness of the event being analyzed. Prior to using the GQ/GL values, it should be determined whether an event is true based on other factors such as the visual features, the amount of read overlap for the breakpoints, the context around the event, such as the number of nearby single nucleotide polymorphisms (SNPs) and insertions/deletions (indels), etc.
