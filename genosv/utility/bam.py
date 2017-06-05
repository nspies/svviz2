import os
import pysam

def bam_sort_index(unsorted_bam_path):
    sorted_bam_path = unsorted_bam_path.replace(".bam", "") + ".sorted.bam"
    pysam.sort("-o", sorted_bam_path, unsorted_bam_path)
    pysam.index(sorted_bam_path)

    os.remove(unsorted_bam_path)