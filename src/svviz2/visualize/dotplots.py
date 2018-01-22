import collections
import itertools
import logging
import numpy
import os
import re
import subprocess
# import tempfile

from svviz2.utility import intervals
from svviz2.utility import misc
from svviz2.visualize import trf
logger = logging.getLogger(__name__)

try:
    from rpy2 import robjects as ro
    from rpy2.robjects import numpy2ri
    numpy2ri.activate()
except ImportError:
    ro = None

class YassException(Exception):
    pass


def can_generate_dotplots():
    if ro is None:
        logger.warn("rpy2 could not be imported; dotplots will not be generated")
        return False
    # try:
    #     subprocess.check_call("yass --version", stderr=subprocess.PIPE, shell=True)
    # except subprocess.CalledProcessError:
    #     logger.warn("yass helper program could not be run; dotplots will not be generated")
    #     return False

    return True

_CAN_GENERATE_DOTPLOTS = can_generate_dotplots()

nucs = list("ACGT")
dinucs = ["".join(x) for x in itertools.product(nucs, repeat=2) if len(set(x))!=1]
trinucs = ["".join(x) for x in itertools.product(nucs, repeat=3) if len(set(x))!=1]

def detect_simple_repeats(seq):
    repeats = trf.run_trf({"a":seq})
    if repeats is not None:
        repeats = repeats.pop("a", [])
        return [(s,e,r) for (s,e,r) in repeats]

    return None
    
    # patterns = {1: ["({}){{10,}}".format(nuc) for nuc in nucs],
    #             2: ["({}){{5,}}".format(dinuc) for dinuc in dinucs],
    #             3: ["({}){{3,}}".format(trinuc) for trinuc in trinucs]}

    # repeats = []

    # for unit, pattern in patterns.items():
    #     pattern = "|".join(pattern)
    #     for match in re.finditer(pattern, seq):
    #         repeats.append((match.start(), match.end(), unit))

    # return repeats

def plot_simple_repeats(s1, s2):
    repeats1 = detect_simple_repeats(s1)
    repeats2 = detect_simple_repeats(s2)

    from_axis = max([len(s1), len(s2)]) * 0.01

    for repeat in repeats1:
        ro.r.segments(repeat[0], from_axis, repeat[1], from_axis, lwd=5, lend=1, col="red")
    for repeat in repeats2:
        ro.r.segments(from_axis, repeat[0], from_axis, repeat[1], lwd=5, lend=1, col="red")


def generate_dotplots(datahub):
    if not _CAN_GENERATE_DOTPLOTS:
        return

    parts = collections.OrderedDict()

    variant = datahub.variant
    old_align_distance = variant.align_distance
    variant.align_distance = 2000

    for allele in ["alt", "ref"]:
        for part in variant.chrom_parts(allele):
            parts[part.id] = part
            print(part.id, part.get_seq())

    import time
    outpath = os.path.join(
        datahub.args.outdir, "{}.dotplots.pdf".format(datahub.variant.short_name()))

    ro.r.pdf(outpath)
    ro.r.par(pty="s") # makes the output graphs square, irrespective of the dimensions of the device

    t0 = time.time()
    for i, part1 in enumerate(parts.keys()):
        for part2 in list(parts.keys())[i:]:
            draw_chrompart_dotplot(parts[part1], parts[part2])

    variant.align_distance = old_align_distance

    for sample_name, sample in datahub.samples.items():
        if sample.sequencer in ["pacbio", "nanopore"]:
            draw_read_dotplots(datahub, sample)

    if not datahub.args.only_realign_locally:
        do_homology_search(datahub)

    ro.r["dev.off"]()

    t1 = time.time()
    logging.debug("TIME for dotplots: {:.1f}s".format(t1-t0))



def draw_chrompart_dotplot(part1, part2):
    """
    compares part1 to part2 -- eg, ref_part1 vs ref_part2
    """

    # XXX refactor
    breakpoints1 = get_breakpoints(part1)
    breakpoints2 = get_breakpoints(part2)

    # print("::", part1.id, part2.id)

    # yass_dotplot(part1.get_seq(), part2.get_seq(),
    #              breakpoints1, breakpoints2,
    #              part1.id, part2.id)

    # plot_simple_repeats(part1.get_seq(), part2.get_seq())

    dotplot = simple_dotplot(part1.get_seq(), part2.get_seq())

    draw_simple_dotplot(dotplot,
        (0,len(part1)), (0,len(part2)),
        breakpoints1, breakpoints2,
        part1.id, part2.id)

    # dotplot2(part1.get_seq(), part2.get_seq())


def get_breakpoints(part):
    print("::", part.id, [len(segment) for segment in part.segments])
    return numpy.cumsum([len(segment) for segment in part.segments])[:-1]

def get_interesting_region(part):
    breakpoints = get_breakpoints(part)
    return breakpoints[0], breakpoints[-1]

def get_interesting_reads(sample, allele, part):
    start, end = get_interesting_region(part)
    bam = sample.outbam(allele, "r")

    reads = []
    for read in bam.fetch(part.id, start, end+1):
        if read.reference_start < start and read.reference_end > end:
            reads.append(read)

    reads.sort(key=lambda x: (x.mapq, x.reference_length), reverse=True)
    return reads


def draw_read_dotplots(datahub, sample):
    reads_to_plot = []
    for allele in ["alt", "ref"]:
        for part in datahub.variant.chrom_parts(allele):
            reads_to_plot.extend(get_interesting_reads(sample, allele, part)[:2])


    for read in reads_to_plot:
        for allele in ["alt", "ref"]:
            for part in datahub.variant.chrom_parts(allele):
                dotplot = simple_dotplot(
                    part.get_seq(), read.query_sequence)

                draw_simple_dotplot(
                    dotplot,
                    (0,len(part)), (0, len(read.query_sequence)),
                    get_breakpoints(part), [],
                    part.id, read.query_name[:15])


# def yass_dotplot(s1, s2, breakpoints1, breakpoints2, label1, label2):
#     tempDir = tempfile.mkdtemp()
    
#     outpaths = []
    
#     for i, seq in enumerate([s1,s2]):
#         tempFastaPath = os.path.join(tempDir, "seq{}.fa".format(i))
#         outpaths.append(tempFastaPath)
#         with open(tempFastaPath, "w") as tempFastaFile:
#             tempFastaFile.write(">seq\n{}".format(seq))
#             tempFastaFile.close()

#     tempYASSResult = os.path.join(tempDir, "result.txt")
    
#     gapExtend = -int(max([len(s1), len(s2)]) / 2      # half the length the entire sequence
#                                              / 10     # 10 nt insertion
#                                              * 5      # the match bonus
#                      )
#     # print(gapExtend)
    
#     # yassCommand = "yass -d 1 -G -50,{} -E 1 {} {}".format(gapExtend, outpaths[0], outpaths[1])
#     # subprocess.check_call(yassCommand, shell=True)

#     yassCommand = "yass -d 3 -G -50,{} -E 1 -o {} {} {}".format(gapExtend, tempYASSResult, outpaths[0], outpaths[1])
#     proc = subprocess.Popen(yassCommand, shell=True,
#         stderr=subprocess.PIPE)
#     resultCode = proc.wait()
    
#     if resultCode != 0:
#         raise YassException("Check that yass is installed correctly")
#     stderr = proc.stderr.readlines()[0].decode()
#     if "Error" in stderr:
#         print("Error running yass: '{}'".format(stderr))
#         raise YassException("Error running yass")

#     ro.r.plot(ro.IntVector([0]), ro.IntVector([0]), type="n", 
#               main="{} : {}".format(label1, label2),
#               xaxs="i", yaxs="i", 
#               xlab="Position in {}".format(label1),
#               ylab="Position in {}".format(label2),
#               xlim=ro.IntVector([0,len(s1)]),
#               ylim=ro.IntVector([0,len(s2)]))
    
#     # for breakpoint in breakpoints1:
#     ro.r.abline(v=breakpoints1, lty=2, col="gray")
        
#     # for breakpoint in breakpoints2:
#     ro.r.abline(h=breakpoints2, lty=2, col="gray")
        
#     for line in open(tempYASSResult):
#         if line.startswith("#"):continue
            
#         res = line.strip().split()
#         if res[6]=="f":
#             ro.r.segments(int(res[0]), int(res[2]), int(res[1]), int(res[3]), col=ro.r.rgb(102, 0, 198, maxColorValue=255), lwd=1)
#         else:
#             ro.r.segments(int(res[1]), int(res[2]), int(res[0]), int(res[3]), col=ro.r.rgb(0, 198, 46, maxColorValue=255), lwd=1)




def simple_dotplot(s1, s2, wordsize=8, scale=650):
    # scale is the final size of the output matrix for visualization

    l1 = int((len(s1)-wordsize))
    l2 = int((len(s2)-wordsize))

    width = int(numpy.ceil(l1/max([l1,l2]) * scale))
    height = int(numpy.ceil(l2/max([l1,l2]) * scale))
    
    mat = numpy.zeros((height, width))
    binsize = l1/(width-1)
    
    kmertopos1 = collections.defaultdict(list)
    
    # get positions of kmers in s1
    for i in range(l1):
        kmer = s1[i:i+wordsize]
        kmertopos1[kmer].append(i)

    # find all matching kmers from s2
    for i in range(l2):
        kmer = s2[i:i+wordsize]
        positions = kmertopos1[kmer]
        
        positions = (numpy.array(positions)/binsize).astype(int)
        y = int(i/binsize)

        mat[y, positions] += 1

    # find all rev-comp kmer matches from s2
    for i in range(l2):
        kmer = misc.reverse_comp(s2[i:i+wordsize])
        positions = kmertopos1[kmer]

        positions = (numpy.array(positions)/binsize).astype(int)
        y = int(i/binsize)

        mat[y, positions] += 1
    
    mat = mat[::-1,]

    return mat

def adjust_boundaries(x1, x2, y1, y2):
    diff = (x2-x1) - (y2-y1)
    if diff > 2:
        y1 -= int(diff/2)
        y2 += int(diff/2)
    elif diff < -2:
        x1 -= int(diff/2)
        x2 += int(diff/2)

    return numpy.array([x1,x2]), numpy.array([y1,y2])

def draw_simple_dotplot(mat, xlim=None, ylim=None, breakpointsx=None, breakpointsy=None, labelx="", labely=""):
    if xlim is None:
        x1 = 0
        x2 = mat.shape[0]
    else:
        x1, x2 = xlim
    if ylim is None:
        y1 = 0
        y2 = mat.shape[1]
    else:
        y1, y2 = ylim

    xlim, ylim = adjust_boundaries(x1, x2, y1, y2)

    main = ""
    if labelx and labely:
        main = "{} : {}".format(labelx, labely)

    ro.r.plot(
           numpy.array([0]),
           xlim=xlim,
           ylim=ylim,
           type="n", bty="n",
           main=main, xlab=labelx, ylab=labely)

    if mat.max() > mat.min():
        mat = mat.max()-mat

        rasterized = ro.r["as.raster"](mat, max=mat.max())
        ro.r.rasterImage(rasterized, x1, y1, x2, y2)

    if breakpointsx is not None:
        ro.r.abline(v=numpy.array(breakpointsx), lty=2, col="gray")
    if breakpointsy is not None:
        ro.r.abline(h=numpy.array(breakpointsy), lty=2, col="gray")




def cluster_loci(loci):
    loci_by_chrom = collections.defaultdict(list)
    for locus in loci:
        loci_by_chrom[locus.chrom].append(locus)

    clustered = []
    for chrom in loci_by_chrom:
        clustered.extend(
            intervals.unionLoci(loci_by_chrom[chrom], 500)
            )

    return clustered

def find_homologous_regions(seq, genomesource, segments, window_size=500, offset=500):
    homologous_regions = []

    for i in range(0, max(1, len(seq)-window_size), offset):
        # print("---", i, "---")

        curseq = seq[i:i+window_size]

        cur_alns = genomesource.bwa.align(curseq, secondary_hit_cutoff=0.5)

        if len(cur_alns) > 0: print("BEST:", cur_alns[0])

        for i, cur_aln in enumerate(cur_alns):
            chrom = genomesource.bwa.ChrIDToName(cur_aln.reference_id)
            locus = intervals.Locus(chrom, cur_aln.reference_start, cur_aln.reference_end, "+")

            if not intervals.overlaps(locus, segments):
                homologous_regions.append(locus)

            # print(cur_aln)
        # print()

    clustered = cluster_loci(homologous_regions)
    # for l in clustered:
        # print(l)

    return clustered

def do_homology_search(datahub):
    logger.info("Finding homologous genomic regions (segmental duplications)...")

    variant = datahub.variant

    segment_loci = []
    for part in variant.chrom_parts("ref"):
        for segment in part.segments:
            cur_chrom = misc.match_chrom_format(segment.chrom, datahub.genome.keys())
            segment = intervals.Locus(cur_chrom, max(0, segment.start-100), segment.end+100, "+")
            segment_loci.append(segment)

    for part in variant.chrom_parts("ref"):
        breakpoints = numpy.cumsum([len(segment) for segment in part.segments])[:-1]

        seq = part.get_seq()
        homologous_regions = find_homologous_regions(
            seq, datahub.genome, segment_loci)
        plot_homologous_regions(seq, homologous_regions, datahub.genome, part.id, breakpoints)

    # also look for any reasonably-sized segments that are unique to the alt allele
    ref_ids = set(segment.id for part in variant.chrom_parts("ref") for segment in part.segments)
    alt_ids = set(segment.id for part in variant.chrom_parts("alt") for segment in part.segments)
    alt_only_ids = alt_ids - ref_ids
    alt_only_segments = [segment for part in variant.chrom_parts("alt") for segment in part.segments
                         if segment.id in alt_only_ids]

    for segment in alt_only_segments:
        if len(segment) >= 50:
            seq = datahub.variant.sources[segment.chrom].get_seq(
                segment.chrom, segment.start, segment.end, segment.strand)

            homologous_regions = find_homologous_regions(
                seq, datahub.genome, [])
            plot_homologous_regions(seq, homologous_regions, datahub.genome, 
                label="segment_{}".format(segment))


def plot_homologous_regions(seq, homologous_regions, genomesource, label, breakpoints=None):
    # from rpy2 import robjects as ro

    # ro.r.pdf("temp.pdf")
    from svviz2.visualize import dotplots

    for region in cluster_loci(homologous_regions):
        if len(region) < len(seq) * 0.9:
            len_diff = int((len(seq)-len(region))/2)
            region._start = max(0, region.start-len_diff)
            region._end += len_diff
            
        other_seq = genomesource.get_seq(
            region.chrom, region.start, region.end, "+").upper()  
        
        d = dotplots.simple_dotplot(seq, other_seq)
        dotplots.draw_simple_dotplot(d, xlim=(0, len(seq)), ylim=(region.start, region.end),
            labelx=label, labely=region.chrom, breakpointsx=breakpoints)

    # ro.r["dev.off"]()
