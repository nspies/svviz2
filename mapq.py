import logging
import numpy
import time

from alignment import Alignment, PairAlignment

logger = logging.getLogger(__name__)

MIN_PHRED_QUALITY = 3.0
BAM_CHARD_CLIP = 5

TAG_END_SCORE = "Es"


def phred_to_prob(phred, phred_scale):
    if phred < 0:
        phred = 0
    return 10.0 ** (phred / -phred_scale)

def prob_to_phred(phred, scale):
    if phred <= 0:
        return numpy.inf

    return -scale * numpy.log10(phred)


def log_sum_exp(array, base=10.0):
    base = float(base)
    array = numpy.asarray(array)
    m = array.max()
    diff = array - m
    sum_of_exps = (base**(diff)).sum()
    return m + numpy.log(sum_of_exps)/numpy.log(base)


class MAPQCalculator(object):
    def __init__(self, reference):
        self.reference = reference
        self.scoring = {
            "match" : 1.0,
            "mismatch": -2.0,
            "gap_open": -3.0,
            "gap_extend": -1.0,
            "clipping_penalty": 0.0,
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

    def get_alignment_end_score(self, aln, add_tag=True):
        if aln.has_tag(TAG_END_SCORE):
            return float(aln.get_tag(TAG_END_SCORE))

        qualities = self.get_qualities(aln)

        aln_start = aln.reference_start
        aln_end = aln_start + aln.reference_length
        ref_chrom = self.reference.get_chrom_by_id(aln.reference_id)
        ref_seq = self.reference.get_seq(ref_chrom, aln_start, aln_end, "+").upper()

        # print(aln.get_tag("AS"), aln.cigarstring)
        log10_score = 0
        in_gap = False
        phred_scale = self.scoring["phred_scale"]

        i = 0
        nm = 0
        for read_pos, ref_pos in aln.get_aligned_pairs():#matches_only=True):
            read_seq = None
            if read_pos is not None:
                read_seq = aln.query_sequence[read_pos]
                read_quality = qualities[read_pos]

            cur_ref_seq = None
            if ref_pos is not None:
                cur_ref_seq = ref_seq[ref_pos-aln.reference_start]

            # print(read_seq, cur_ref_seq)

            # print(i, read_pos, read_seq, ref_pos, cur_ref_seq, log10_score)

            if i < aln.query_alignment_start or i >= aln.query_alignment_start+aln.query_alignment_length:
                # print("  clipped")
                log10_score += read_quality / -phred_scale + self.scoring["clipping_penalty"]
                pass
            elif read_seq is None:
                # deletion
                if in_gap:
                    log10_score += read_quality / -phred_scale + self.scoring["gap_extend"]
                else:
                    log10_score += read_quality / -phred_scale + self.scoring["gap_open"]
                    in_gap = True
            elif cur_ref_seq is None:
                # insertion
                if in_gap:
                    log10_score += read_quality / -phred_scale + self.scoring["gap_extend"]
                else:
                    log10_score += read_quality / -phred_scale + self.scoring["gap_open"]
                    in_gap = True
            elif read_seq != cur_ref_seq:
                # mismatch
                in_gap = False
                log10_score += read_quality / -phred_scale
                nm += 1
            else:
                # match
                in_gap = False
                prob = 1 - phred_to_prob(read_quality, self.scoring["phred_scale"])
                log10_score += prob_to_phred(prob, -1)

            i += 1

            # print(read_pos, read_seq, ref_pos, ref_seq, match)
        if aln.cigartuples[0][0] == BAM_CHARD_CLIP:
            # log10_score += read_quality / -phred_scale + self.scoring["gap_open"]
            log10_score += (read_quality / -phred_scale + self.scoring["clipping_penalty"]) * (aln.cigartuples[0][1])
        if aln.cigartuples[-1][0] == BAM_CHARD_CLIP:
            # log10_score += read_quality / -phred_scale + self.scoring["gap_open"]
            log10_score += (read_quality / -phred_scale + self.scoring["clipping_penalty"]) * (aln.cigartuples[-1][1])

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



    def score_pairs(self, alns1, alns2, read_stats):
        pairs = get_concordant_pairs(alns1, alns2, read_stats)
        for pair in pairs:
            pair.score = self._score_pair(pair, read_stats)
        pairs.sort(key=lambda x: x.score, reverse=True)
        pairs = [pair for pair in pairs if pair.score!=0]
        return pairs

    def _score_pair(self, read_pair, read_stats):
        # t0 = time.time()

        # these are log10 likelihoods
        score1 = self.get_alignment_end_score(read_pair.read1)
        score2 = self.get_alignment_end_score(read_pair.read2)

        # this is a probability
        insert_size_prob = read_stats.scoreInsertSize(read_pair.insert_size)

        with numpy.errstate(divide="ignore"):
            log10_pair_prob = numpy.log10(insert_size_prob) + score1 + score2

        # print("PAIR:", score1, score2, pair_prob, prob_to_phred(pair_prob, 10))
        # print(pair_prob, prob_to_phred(pair_prob, 10))

        # t1 = time.time()
        # logger.info("time:{:.4f}".format(t1-t0))

        return log10_pair_prob

def set_mapqs(alns):
    if len(alns) == 0:
        return

    # first, let's make sure all the scores are within 300 of one another; scores that are
    # smaller than best_score-300 will be zeroed out
    scores = numpy.array([aln.score for aln in alns])
    best_score = scores.max()

    if best_score < -300:
        correction = best_score + 300
        scores = scores - correction
        scores[scores>0] = 0

    probs = 10 ** (scores)
    total = probs.sum()

    # print("::::", total)#, [aln.score for aln in alns])

    for aln, prob in zip(alns, probs):
        aln.posterior = prob / total
        aln.mapq = int(min(prob_to_phred(1-aln.posterior, 10), 40))

        # print(aln.get_tag("AS"), aln.score, aln.posterior, prob_to_phred(1-aln.posterior, 10), aln.mapq)

        # # if len(pairs) > 1:
        # print(":::::::::::::::::::")
        # for pair in pairs:
        #     print(pair.insert_size)#, read_stats.scoreInsertSize(pair.insert_size))
        #     print(pair.read1)
        #     print(pair.read2)
        #     print(pair.score, pair.posterior, "mapq:", float(pair.mapq))


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