import collections
import numpy
import pysam
import pytest
import random

from svviz2.remap.readpair import ReadPair
from svviz2.remap.alignment import Alignment
from svviz2.app.genomesource import GenomeSource
from svviz2.utility.misc import reverse_comp
from svviz2.app.sample import ReadStatistics

class _ReadStats(ReadStatistics):
    def __init__(self):
        self.insertSizes = [int(x) for x in numpy.random.normal(400, 25, 50000)]
        self.readLengths = []
        self.orientations = ["+-"]
        self.number_mismatches = []

        self._insertSizeKDE = None
        self.singleEnded = False

        self._insertSizeScores = {} # cache
        self._max_insert_size = None
read_stats = _ReadStats()


def simulate_read_pair(sequence, start, length=150, isize=400, flip=False):
    r1 = pysam.AlignedSegment()
    r1.query_sequence = sequence[start:start+length]

    r2 = pysam.AlignedSegment()
    pos2 = start+isize
    r2.query_sequence = reverse_comp(sequence[pos2-length:pos2])

    if flip:
        r1,r2 = r2,r1
    return ReadPair(Alignment(r1), Alignment(r2), read_stats)

def simulate_read_pairs(sequence, n_pairs, length=150, isize=400, r=random.Random()):
    reads = []
    for i in range(n_pairs):
        start = r.randint(0, len(sequence)-isize)
        
        pair = simulate_read_pair(sequence, start, length, isize)

        reads.append((start, start+isize-1, pair))
    return reads

# def get_error_type():
#     r = numpy.random.rand()
#     if r < 0.4:
#         return "insertion"
#     elif r < 0.8:
#         return "deletion"
#     else:
#         return "mismatch"

# def simulate_long_reads(seq, n, mean_length=10000, error_rate=0.0):
#     reads = []

#     for i in range(n):
#         cur_length = min(numpy.random.poisson(mean_length), len(seq)-1)
#         cur_start = numpy.random.randint(0, len(seq)-cur_length-1)

#         cur_seq = list(seq[cur_start:(cur_start+cur_length)])

#         errors = int(cur_length * error_rate)

#         while errors > 0:
#             error_type = get_error_type()
#             position = numpy.random.randint(0, len(cur_seq))

#             if error_type == "insertion":
#                 ins_length = numpy.random.poisson(2.7)
#                 cur_seq.insert(position, random_sequence(ins_length))
#                 errors -= ins_length
#             elif error_type == "deletion":
#                 del_length = numpy.random.poisson(2.7)
#                 cur_seq = cur_seq[:position] + cur_seq[position+del_length:]
#                 errors -= del_length
#             elif error_type == "mismatch":
#                 cur_seq[position] = numpy.random.choice(list(set("ACGT")-set(cur_seq[position])))
#                 errors -= 1

#         reads.append((cur_start, "".join(cur_seq)))

#     return reads


def _seq1():
    import random
    r = random.Random()
    r.seed(551)
    return "".join(r.choice("ACGT") for i in range(2000))


def _seq2():
    import random
    r = random.Random()
    r.seed(135)
    return "".join(r.choice("ACGT") for i in range(10000))

@pytest.fixture
def genome_source():
    gs = GenomeSource({"chr1":_seq1(),
                       "chr2":_seq2()})
    return gs

@pytest.fixture(params=[10, 50, 150, 500, 1000, 2000])
def genome_source_deletion(request):
    seq2 = _seq2()
    seq2 = seq2[:4000] + seq2[4000+request.param:]
    gs = GenomeSource({"chr1":_seq1(),
                       "chr2":seq2})
    return gs, request.param

@pytest.fixture
def heterozygous_readpairs_deletion(genome_source_deletion):
    r = random.Random()
    r.seed(525)

    ref = genome_source()
    alt, deletion_length = genome_source_deletion
    alleles_to_reads = collections.defaultdict(list)

    for allele, gs in [("ref", ref), ("alt", alt)]:
        for chrom in gs.keys():
            seq = gs.names_to_contigs[chrom]
            n = int(len(seq) / 300 * 1000)
            cur_pairs = simulate_read_pairs(seq, n, r=r)
            event_start = 4000
            event_end = 4000+deletion_length
            if allele == "alt":
                event_end = 4000

            for p1, p2, read_pair in cur_pairs:
                # read_pair = ReadPair(read1, read2, read_stats)
                if p1 > event_end or p2 < event_start:
                    alleles_to_reads["amb"].append(read_pair)
                elif p1 > event_end-5 or p2 < event_start+5:
                    ...
                else:
                    alleles_to_reads[allele].append(read_pair)

    return alleles_to_reads

if __name__ == '__main__':
    alleles_to_reads = heterozygous_readpairs_deletion()
    for key in alleles_to_reads:
        print(key, len(alleles_to_reads[key]))
