import time

from alignment import MultiReferenceReadEndAlignmentSet, MultiReferenceAlignmentSet
import mapq
import combinedreference

class Realigner(object):
    def __init__(self, ref_genome_source, breakpoints=None):
        self.ref_genome_source = ref_genome_source
        self.local_ref_genome = None
        self.local_alt_genome = None
        self.breakpoints = breakpoints

        self._ref_mapper = None
        self._alt_mapper = None

        self._ref_mapq_calculator = None
        self._alt_mapq_calculator = None


    def _realign(self, seq, name):
        extra_args = {"hardclip":False, "secondary_hit_cutoff":0.6, "max_secondary":50}
        ref_mappings = self.ref_mapper.align(seq, name, **extra_args)
        alt_mappings = self.alt_mapper.align(seq, name, **extra_args)

        return ref_mappings, alt_mappings

    def realign_single_end(self, read):
        seq = read.original_sequence()
        name = read.query_name

        t0 = time.time()
        ref_alns, alt_alns = self._realign(seq, name)
        t1 = time.time()

        tx0 = time.time()
        ref_alns = self.ref_mapq_calculator.score_single_end_alns(ref_alns)
        alt_alns = self.alt_mapq_calculator.score_single_end_alns(alt_alns)

        mapq.set_mapqs(ref_alns+alt_alns)
        # if len(ref_pairs+alt_pairs) > 1:
        # print(["{:.2g}|{}".format(p.score, p.mapq) for p in ref_pairs],
        #       ["{:.2g}|{}".format(p.score, p.mapq) for p in alt_pairs])


        # if read.query_name == "m150105_192231_42177R_c100761782550000001823161607221526_s1_p0/138972/39862_46995":
        # print("*"*200)
        # print("...ref...")
        # for aln in ref_alns:
        #     print(aln.get_tag("AS"), aln.mapq)#, aln._read)
        #     # print(aln.reference_start, aln.cigarstring, aln.get_tag("AS"))
        # print("...alt...")
        # for aln in alt_alns:
        #     print(aln.get_tag("AS"), aln.mapq)#, aln._read)
        #     # print(aln.reference_start, aln.cigarstring, aln.get_tag("AS"))
        # print(".........")

        aln_set = MultiReferenceAlignmentSet(self.ref_genome_source, self.alt_genome_source)
        aln_set.add_alignments(ref_alns, "ref")
        aln_set.add_alignments(alt_alns, "alt")
        tx1 = time.time()

        # print("TIME:", t1-t0, tx1-tx0)
        return aln_set


    def realign_pair(self, pair, read_statistics):
        seq1 = pair.read1.original_sequence()
        seq2 = pair.read2.original_sequence()

        end_aln_set = self._realign_pair(seq1, seq2, pair.read1.query_name)
        if pair.read1.query_name == "HA2WPADXX:19:5:1605415:0":
            print("ENDS SEPARATELY:")
            print("...ref...")
            for end in [1,2]:
                print("{}".format(end)*30)
                for aln in end_aln_set.get_alns("ref", end):
                    print(aln)
            print("...alt...")
            for end in [1,2]:
                for aln in end_aln_set.get_alns("alt", end):
                    print(aln)
            print(".........")

        ref_pairs = self.ref_mapq_calculator.score_pairs(
            end_aln_set.get_alns("ref", 1), end_aln_set.get_alns("ref", 2), read_statistics)
        alt_pairs = self.alt_mapq_calculator.score_pairs(
            end_aln_set.get_alns("alt", 1), end_aln_set.get_alns("alt", 2), read_statistics)

        mapq.set_mapqs(ref_pairs+alt_pairs)
        # if len(ref_pairs+alt_pairs) > 1:
        # print(["{:.2g}|{}".format(p.score, p.mapq) for p in ref_pairs],
        #       ["{:.2g}|{}".format(p.score, p.mapq) for p in alt_pairs])

        paired_set = MultiReferenceAlignmentSet(end_aln_set.ref_source, end_aln_set.alt_source)
        paired_set.add_alignments(ref_pairs, "ref")
        paired_set.add_alignments(alt_pairs, "alt")

        return paired_set

    def _realign_pair(self, seq1, seq2, name):
        aln_set = MultiReferenceReadEndAlignmentSet(self.ref_genome_source, self.alt_genome_source)
        
        seqs = {1:seq1, 2:seq2}
        for end in [1,2]:
            ref_mappings, alt_mappings = self._realign(seqs[end], name)
            if end == 2:
                for mapping in ref_mappings+alt_mappings:
                    mapping.is_read2 = True
            aln_set.add_alignments(ref_mappings, "ref", end)
            aln_set.add_alignments(alt_mappings, "alt", end)

        return aln_set

    @property
    def ref_mapper(self):
        if self._ref_mapper is None:
            # self._ref_mapper = self.ref_genome_source.get_bwa()
            self._ref_mapper = combinedreference.CombinedAligner(
                self.ref_genome_source,
                self.local_ref_genome,
                self.local_coords_in_full_genome
                )
        return self._ref_mapper

    @property
    def alt_mapper(self):
        if self._alt_mapper is None:
            self._alt_mapper = self.local_alt_genome.get_bwa()
        return self._alt_mapper

    @property
    def ref_mapq_calculator(self):
        if self._ref_mapq_calculator is None:
            self._ref_mapq_calculator = mapq.MAPQCalculator(self.ref_genome_source)
        return self._ref_mapq_calculator

    @property
    def alt_mapq_calculator(self):
        if self._alt_mapq_calculator is None:
            self._alt_mapq_calculator = mapq.MAPQCalculator(self.alt_genome_source)
        return self._alt_mapq_calculator

    def set_local_alt_genome(self, local_alt_genome_source):
        self.local_alt_genome = local_alt_genome_source
        self._alt_mapper = None
        self._alt_mapq_calculator = None

    def set_local_ref_genome(self, local_ref_genome_source, local_coords_in_full_genome):
        self.local_ref_genome = local_ref_genome_source
        self.local_coords_in_full_genome = local_coords_in_full_genome
        self._ref_mapper = None
        self._ref_mapq_calculator = None

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["_ref_mapper"]
        del state["_alt_mapper"]
        return state



# class BWAMapper(object):
#     """
#     thin wrapper for the pybind11/C++ BWAWrapper class
#     """
#     def __init__(self, ref_path=None, ref_seqs=None, mapq_calculator=None, as_raw=False):
#         self.reference = genomesource.make_genome_source(ref_path=ref_path, ref_seqs=ref_seqs)
#         self.bwa = self.reference.get_bwa()
#         self.as_raw = as_raw

#     def align(self, seq, read_name="read"):
#         alns = self.bwa.align(
#             seq, read_name=read_name, hardclip=False, secondary_hit_cutoff=0.1, max_secondary=10)

#         if not self.as_raw:
#             alns = [Alignment(aln) for aln in alns]

#         return alns





def test():
    # mapper = BWAMapper(ref_path="/Volumes/frida/nspies/data/bwa-hg19/hg19.fasta")
    mapper = BWAMapper(ref_seqs={"chr1":"TTTGAGCTGAGAGCGGCGTATGATGCTGTATATGCGCGGCGGGA"
                                        "GGAGCCGGAGAGCGGATTAATTCTATATATCGCGGGCGAGAGGC"
                                        "GCGATCTTAGGCGGGCAACCCACGCGATATCGAGG"
                                        "GGAGCGAGAGTTTCTTATATATTCTTATATAAACACACCACCCC"
                                        "GGGGGGAGAGGAAACCCTTTCTATATCTCTTTTATCTTCTCTTC"
                                        "ACCACCCCACATTTCTCTCTTTATCTTAGATGCTGTAGTCGTCG"
                                        })
    alignments = mapper.align("GCGATTGCTTAGGCGGGCAACCCACGCGATATCGAGG")

    for aln in alignments:
        print(aln.get_tags(), aln.cigarstring)#("YS"))

if __name__ == '__main__':
    test()
