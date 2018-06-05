Reports - quantitative assessment of SVs
----------------------------------------

svviz2 generates a number of useful statistics for each event and sample it analyzes. The resulting reports can easily be automatically parsed and analyzed by downstream applications. A separate report is generated for each event, with filename ``<event_id>.<coordinates>.report.tsv``. The report contains four columns:

- sample - which sample being analyzed, corresponding to an input bam file
- allele - the allele being reported; this is either alt, ref or amb; or empty in the case of a few overall descriptors such as genotype
- key - the name of the statistic
- value - the value of the statistic

A partial description of the different statistics generated follows:

``count``
    The simple count of reads supporting the given allele.

``weighted_count``
    The weighted count, where each read contributes mapq/40.0 to the sum (see :ref:`weighted mapq <weighted_mapq>`). For example, a read supporting the alt allele with mapq 30 would contribute 3/4=0.75 to the weighted count.

``GL_count``
    The simple :ref:`genotype likelihood <genotypes>`.

``GQ_count``
    The simple :ref:`genotype quality <genotypes>`.

``GT_count``
    The genotype inferred from ``GQ_count``

``GL_mapq``
    The genotype likelihood treating the mapq as a simple probability. For example, a read supporting the alt allele with mapq 30 would be presumed to have derived from the alt allele with probability 99.9%

``GQ_mapq``
    The genotype quality treating mapq as a simple probability

``GT_count``
    The genotype inferred from ``GQ_mapq``

``GL_weighted``
    The genotype likelihood using the proportional weighted mapq/40 value treated as a probability (see :ref:`here <weighted_mapq>` for more info)

``GQ_weighted``
    The corresponding genotype quality using the proportional weighted mapq

``GT_weighted``
    The genotype inferred from ``GQ_weighted``

``<locus>_n_mismatch_low``
    The number of putative SNPs in the region of the breakpoint; high mismatch values suggest the variant is in a segmental duplication. This is inferred by looking at each position upstream or downstream 100 bp of the breakpoints, as well as the regions between breakpoints, then identifying positions where at least 20% of reads show a specific mismatch to the reference.

``<locus>_n_mismatch_high``
    The number of putative SNPs in the region of the breakpoint; high mismatch values suggest the variant is in a segmental duplication. This is similar to ``*_n_mismatch_low``, except 80% of reads at a given position must show a specific mismatch to the reference to count as a putative SNP.

``<locus>_n_indel_low``
    The number of putative insertion/deletions near the breakpoints. 20% of reads must support an indel at a given position to be counted as an indel.

``<locus>_n_indel_high``
    The number of putative insertion/deletions near the breakpoints. 80% of reads must support an indel at a given position to be counted as an indel.

``<locus>_n_total``
    The number of positions assessed for SNPs and indels

``overlap_<allele>:<coordinate>``
    The mean number of bases that reads extend across the breakpoint at <coordinate> supporting <allele>. For example, a 50bp read with 46bp to the left of the breakpoint and 4bp to the right of the breakpoint would be considered to extend 4bp (the lesser of 46 and 4). In general, a low number here indicates most reads don't extend particularly far into the SV, suggesting these reads are either mismapping or the sequence or coordinates of the SV are incorrect.

``extension_<allele>:<coord>``
    The mean number of bases that reads extend to the right of the breakpoint at <coord>. This is similar to the overlap_* statistic above, except we're no longer taking the lesser of 46 and 4; we just take the second value.

``count_<allele>.<coord>_seq``
    The number of reads whose sequence spans the breakpoint at <coord>. 

``count_<allele>.<coord>_pair``
    The number of reads whose mate pairs span the breakpoint at <coord>. Typically, true events should have both ``*_seq`` and ``*_pair`` counts in paired-end datasets, although complex sequence patterns surrounding the breakpoints may mean this isn't always true.
