import collections
import logging
import numpy
# import time

from svviz2.remap.alignment import Alignment
from svviz2.utility.statistics import phred_to_prob, prob_to_phred

logger = logging.getLogger(__name__)

MIN_PHRED_QUALITY = 3.0
BAM_CHARD_CLIP = 5

TAG_END_SCORE = "Es"



class MAPQCalculator(object):
    def __init__(self, reference):
        self.reference = reference
        self.scoring = {
            "match" : 0.0,       # when set to zero, these probabilities are calculated
            "mismatch": 0.0,     # exclusively from the base qualities

            "gap_open": -1.0,
            "gap_extend": -1.0,
            "clipping_penalty": 0.0,
            "min_clip_length": 10,
            "discordant_penalty": -15,

            "phred_scale": 10.0,

            "min_base_quality": 0,
            "max_base_quality": 10
        }

        self.rescale_base_qualities = True

    def get_qualities(self, aln):
        if aln.query_qualities is None:
            qualities = numpy.repeat(40, aln.query_length)
        else:
            qualities = numpy.asarray(aln.query_qualities)
        qualities[qualities < MIN_PHRED_QUALITY] = MIN_PHRED_QUALITY

        if self.rescale_base_qualities:
            scale = self.scoring["max_base_quality"] - self.scoring["min_base_quality"]
            qualities = (qualities / 40.0) * scale + self.scoring["min_base_quality"]

        return qualities

    # @profile
    def get_alignment_end_score(self, aln, add_tag=True):
        if aln.has_tag(TAG_END_SCORE):
            return float(aln.get_tag(TAG_END_SCORE))

        qualities = self.get_qualities(aln)

        aln_start = aln.reference_start
        aln_end = aln_start + aln.reference_length
        # ref_chrom = self.reference.get_chrom_by_id(aln.reference_id)
        ref_chrom = aln.chrom
        ref_seq = self.reference.get_seq(ref_chrom, aln_start, aln_end, "+").upper()

        # print(aln.get_tag("AS"), aln.cigarstring)
        log10_score = 0
        in_gap = False
        phred_scale = self.scoring["phred_scale"]

        i = 0
        nm = 0

        clip_left_adjust =  min(1.0, aln.cigartuples[0][1] / self.scoring["min_clip_length"])
        clip_right_adjust = min(1.0, aln.cigartuples[-1][1] / self.scoring["min_clip_length"])

        penalties = collections.defaultdict(int)

        query_sequence = aln.query_sequence
        query_end = aln.query_alignment_start+aln.query_alignment_length
        for read_pos, ref_pos in aln.get_aligned_pairs():#matches_only=True):
            read_seq = None
            if read_pos is not None:
                read_seq = query_sequence[read_pos]
                read_quality = qualities[read_pos]

            cur_ref_seq = None
            if ref_pos is not None:
                cur_ref_seq = ref_seq[ref_pos-aln.reference_start]

            # print(read_seq, cur_ref_seq)

            # print(i, read_pos, read_seq, ref_pos, cur_ref_seq, log10_score)

            if i < aln.query_alignment_start:
                penalties["clip_left"] += (clip_left_adjust) * (read_quality / -phred_scale + self.scoring["clipping_penalty"])
            elif i >= query_end:
                penalties["clip_right"] += (clip_right_adjust) * (read_quality / -phred_scale + self.scoring["clipping_penalty"])
            elif read_seq is None:
                # deletion
                if in_gap:
                    penalties["deletion_extend"] += read_quality / -phred_scale + self.scoring["gap_extend"]
                else:
                    penalties["deletion_open"] += read_quality / -phred_scale + self.scoring["gap_open"]
                    in_gap = True
            elif cur_ref_seq is None:
                # insertion
                if in_gap:
                    penalties["insertion_extend"] += read_quality / -phred_scale + self.scoring["gap_extend"]
                else:
                    penalties["insertion_open"] += read_quality / -phred_scale + self.scoring["gap_open"]
                    in_gap = True
            elif read_seq != cur_ref_seq:
                # mismatch
                in_gap = False
                penalties["mismatch"] += read_quality / -phred_scale + self.scoring["mismatch"]
                nm += 1
            else:
                # match
                in_gap = False
                prob = 1 - phred_to_prob(read_quality, self.scoring["phred_scale"])
                penalties["match"] += prob_to_phred(prob, -1) + self.scoring["match"]

            i += 1

            # print(read_pos, read_seq, ref_pos, ref_seq, match)
        if aln.cigartuples[0][0] == BAM_CHARD_CLIP:
            raise Exception("hard clipping not yet supported") # should match soft-clipping
            # log10_score += read_quality / -phred_scale + self.scoring["gap_open"]
            penalties["hard_clip_left"] += (read_quality / -phred_scale + self.scoring["clipping_penalty"]) * (aln.cigartuples[0][1])
        if aln.cigartuples[-1][0] == BAM_CHARD_CLIP:
            raise Exception("hard clipping not yet supported") # should match soft-clipping
            # log10_score += read_quality / -phred_scale + self.scoring["gap_open"]
            log10_score += (read_quality / -phred_scale + self.scoring["clipping_penalty"]) * (aln.cigartuples[-1][1])

        # print("*"*23)
        # print("ADJUST:", aln.cigartuples, clip_left_adjust, clip_right_adjust)
        # print("GOOD:", prob_to_phred(1 - phred_to_prob(max(qualities), self.scoring["phred_scale"]), -1))
        # print("BAD:", max(qualities) / -phred_scale + self.scoring["mismatch"])
        # print(aln.locus)
        # for key in sorted(penalties):
        #     print("{:>20}  {:>6.2f}".format(key, penalties[key]))

        log10_score = sum(penalties.values())

        if add_tag:
            aln.set_tag(TAG_END_SCORE, str(log10_score))
            if not aln.has_tag("NM"):
                aln.set_tag("NM", nm)

        return log10_score


    def score_single_end_alns(self, alns):
        alns = [Alignment(aln) for aln in alns]
        for aln in alns:
            aln.score = self.get_alignment_end_score(aln)
        alns.sort(key=lambda x: x.score, reverse=True)
        alns = [aln for aln in alns if aln.score!=0]
        return alns


    def score_pairs_new(self, alns1, alns2, read_stats):
        pairs = []
        for aln1 in alns1:
            for aln2 in alns2:
                pair = PairAlignment(aln1, aln2)
                pair.score = self._score_pair(pair, read_stats)

                pairs.append(pair)

        pairs.sort(key=lambda x: x.score, reverse=True)
        pairs = [p for p in pairs if p.score!=0]   

        return pair

    def score_pairs(self, alns1, alns2, read_stats):
        pairs = get_concordant_pairs(alns1, alns2, read_stats)
        for pair in pairs:
            pair.score = self._score_pair(pair, read_stats)
        pairs.sort(key=lambda x: x.score, reverse=True)
        pairs = [pair for pair in pairs if pair.score!=0]
        return pairs

    def _score_pair(self, read_pair, read_stats):
        # these are log10 likelihoods
        score1 = self.get_alignment_end_score(read_pair.read1)
        score2 = self.get_alignment_end_score(read_pair.read2)

        if not read_pair.concordant(read_stats):
            return score1 + score2 + self.scoring["discordant_penalty"]

        # this is a probability
        insert_size_prob = read_stats.scoreInsertSize(read_pair.insert_size)

        with numpy.errstate(divide="ignore"):
            log10_pair_prob = numpy.log10(insert_size_prob) + score1 + score2

        # print("PAIR:", score1, score2, pair_prob, prob_to_phred(pair_prob, 10))
        # print(pair_prob, prob_to_phred(pair_prob, 10))

        return log10_pair_prob


def get_concordant_pairs(alns1, alns2, read_stats):
    pairs = []
    for aln1 in alns1:
        for aln2 in alns2:
            pair = PairAlignment(aln1, aln2)
            if pair.concordant(read_stats):
                pairs.append(pair)

            # if pair.query_name == "D00360:99:C8VWFANXX:3:1307:7153:24671":
            #     print(pair.concordant(read_stats))
            #     print(pair.read1)
            #     print(pair.read2)
            #     print()
    return pairs
