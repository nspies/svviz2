Low-level details
=================

Mapping quality calculations
----------------------------

Each read analyzed by svviz2 is assigned a mapping quality (mapq) indicating the support for the indicated allele. For example, a read strongly supporting the alt allele might receive a mapq of 35 whereas a read that is ambiguous might receive a mapq of 0.

First, an alignment score is calculated according to the methods described in Zheng and Grice (PLoS Comput Biol, 2016). Briefly, each position that is aligned correctly (ie the bases match) adds probability proportional to 1-baseQ, ie proportional to the probability of a sequencing error at that position. Each position that is in a mismatch, indel or is clipped adds probability proportional to baseQ. A fixed penalty is added for each type of gap-open/gap-extend/clip. For paired-end reads, the resulting probabilities are multiplied together with the empirical insert size probability.

The alignment scores are then used to calculate a mapq, which is equal to the probability of a given alignment divided by the sum of the probabilities of all alignments. A read may align to one or more locations in the reference region (ie close to the breakpoints), the alt region, or elsewhere in the genome. All of these mappings are considered when calculating mapq.

Genotype likelihood and quality scores
--------------------------------------

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

