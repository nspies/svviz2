import collections
import logging


# from svviz2.remap import mapq
from svviz2.remap import alignment

logger = logging.getLogger(__name__)

class ReadPair(object):
    """ Represents a read-pair

    only one ReadPair exists per read-pair, but many Alignments and AlignmentPairs can exist
    per read-pair
    """
    def __init__(self, original_read1, original_read2, read_stats):
        self.original_read_ends = {
            "1": original_read1,
            "2": original_read2,
        }

        assert original_read1.query_name == original_read2.query_name
        self.name = original_read1.query_name

        self.ref_pairs = []
        self.alt_pairs = []

        # self.chosen_alignment_pair = None

        self.read_stats = read_stats

    def realign_against_allele(self, genome_sources):
        """ use this to realign against one allele """
        end_alns = collections.defaultdict(list)
        for end, original_read in self.original_read_ends.items():
            for genome_source in genome_sources:
                # this calculates the alignment score
                cur_alns = genome_source.align(original_read)
                end_alns[end].extend(cur_alns)

                # if self.original_read_ends["1"].query_name == "HA2WPADXX:21:3:257317:0":    
                #     for aln in cur_alns:
                #         print(aln._read)



        pair_alns = []
        for aln1 in end_alns["1"]:
            for aln2 in end_alns["2"]:
                pair = alignment.AlignmentPair(aln1, aln2, name=self.name)
                self.read_stats.score_read_pair(pair)
                pair_alns.append(pair)

        return pair_alns

    def realign(self, ref_genome_sources, alt_genome_sources):
        self.ref_pairs = self.realign_against_allele(ref_genome_sources)
        self.alt_pairs = self.realign_against_allele(alt_genome_sources)

        alignment.set_mapqs(self.ref_pairs+self.alt_pairs)
        self.ref_pairs.sort(key=lambda x: x.score, reverse=True)
        self.alt_pairs.sort(key=lambda x: x.score, reverse=True)

        # if self.original_read_ends["1"].query_name == "HA2WPADXX:21:3:257317:0":    
        #     print(self.name)
        #     print("REF")
        #     for pair in self.ref_pairs:
        #         print(pair.loci, pair.score, pair.aln1.cigarstring, pair.aln2.cigarstring, 
        #             pair.aln1.score, pair.aln2.score,
        #             pair.concordant(self.read_stats),
        #             pair.mapq)
        #     print("ALT")
        #     for pair in self.alt_pairs:
        #         print(pair.loci, pair.score, pair.aln1.cigarstring, pair.aln2.cigarstring, 
        #             pair.aln1.score, pair.aln2.score,
        #             pair.concordant(self.read_stats),
        #             pair.mapq)

        #     print("///")


        # if len(self.ref_pairs)+len(self.alt_pairs) > 0:
        #     self.chosen_alignment_pair = max(self.ref_pairs+self.alt_pairs, 
        #                                      key=lambda x: x.score)

