from genosv.utility.bam import bam_sort_index

def save_realignments(aln_sets, sample, datahub):
    for aln_set in aln_sets:
        if aln_set.supports_allele != "amb":
            aln_set.supporting_aln.fix_flags()

        if aln_set.supports_allele == "ref":
            if sample.single_ended:
                sample.out_ref_bam.write(aln_set.supporting_aln._read)
            else:
                sample.out_ref_bam.write(aln_set.supporting_aln.aln1._read)
                sample.out_ref_bam.write(aln_set.supporting_aln.aln2._read)

        elif aln_set.supports_allele == "alt":
            if sample.single_ended:
                sample.out_alt_bam.write(aln_set.supporting_aln._read)
            else:
                sample.out_alt_bam.write(aln_set.supporting_aln.aln1._read)
                sample.out_alt_bam.write(aln_set.supporting_aln.aln2._read)

    sample.out_alt_bam.close()
    sample.out_ref_bam.close()

    bam_sort_index(sample.ref_bam_path)
    bam_sort_index(sample.alt_bam_path)