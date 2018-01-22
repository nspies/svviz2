import collections
import pysam

from svviz2.utility.misc import reverse_comp

try:
    from ssw import ssw_wrap
except ImportError:
    print("ssw library not found")


class Aligner(object):
    def __init__(self, names_to_contigs):
        self.names_to_aligners = collections.OrderedDict()
        for name, contig in names_to_contigs.items():
            self.names_to_aligners[name] = ssw_wrap.Aligner(
                contig, report_cigar=True, report_secondary=True)

    def align(self, seq):
        alns = []
        revseq = reverse_comp(seq)

        for i, name in enumerate(self.names_to_aligners):
            aligner = self.names_to_aligners[name]
            faln = aligner.align(seq)
            raln = aligner.align(revseq)

            cur_aln = pysam.AlignedSegment()
            cur_aln.reference_id = i
            if faln.score > raln.score:
                cur_aln.query_sequence = seq
                cur_aln.reference_start = faln.ref_begin
                cur_aln.set_tag("AS", faln.score)
                cur_aln.cigarstring = faln.cigar_string
            else:
                cur_aln.query_sequence = revseq
                cur_aln.reference_start = raln.ref_begin
                cur_aln.set_tag("AS", raln.score)
                cur_aln.cigarstring = raln.cigar_string

                cur_aln.is_reverse = True

            alns.append(cur_aln)

        return alns