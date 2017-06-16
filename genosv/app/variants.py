import collections
# import logging

# from svviz.utilities import Locus, getListDefault
# import genomesource
from genosv.utility.misc import safe_file_name
from genosv.utility.intervals import Locus
from genosv.app import genomesource


def non_negative(x):
    return max(x, 0)
    
def get_breakpoints_on_original_reference(sv):
    chrom_parts = sv.chrom_parts("alt")

    breakpoints = []

    for part in chrom_parts:
        for i in range(len(part.segments)):#segment in part.segments[:-1]:
            segment = part.segments[i]
            print("SEGMENT", segment)
            if i > 0:
                # breakpoint = Locus(segment.chrom, segment.start, segment.start, "+")
                breakpoint = segment.fiveEndLocus()
                breakpoints.append(breakpoint)
            elif i < len(part.segments)-1:
                # breakpoint = Locus(segment.chrom, segment.end, segment.end, "+")
                breakpoint = segment.threeEndLocus()
                breakpoints.append(breakpoint)

    return breakpoints


def get_breakpoints_on_local_reference(sv, allele):
    chrom_parts = sv.chrom_parts(allele)

    breakpoints = []

    for part in chrom_parts:
        cur_pos = 0
        for segment in part.segments[:-1]:
            cur_pos += len(segment)
            breakpoints.append(Locus(part.id, cur_pos, cur_pos, "+"))

            # print("SEGMENT", segment)
            # if i > 0:
            #     # breakpoint = Locus(segment.chrom, segment.start, segment.start, "+")
            #     breakpoint = segment.fiveEndLocus()
            #     breakpoints.append(breakpoint)
            # elif i < len(part.segments)-1:
            #     # breakpoint = Locus(segment.chrom, segment.end, segment.end, "+")
            #     breakpoint = segment.threeEndLocus()
            #     breakpoints.append(breakpoint)

    return breakpoints


class ChromPart(object):
    def __init__(self, regionID, segments, sources):
        self.id = regionID
        self.segments = segments
        self.sources = sources
        self._seq = None

    def get_seq(self, start=0, end=None):
        if self._seq is None:
            seqs = []
            for segment in self.segments:
                seq = self.sources[segment.source].get_seq(segment.chrom, segment.start, segment.end, segment.strand)
                seqs.append(seq)
            self._seq = "".join(seqs).upper()
        if end is None:
            end = len(self._seq)
        return self._seq[start:end]

    def __len__(self):
        return len(self.get_seq())

    def __repr__(self):
        return "{}:{}".format(self.id, self.segments)

class ChromPartsCollection(object):
    def __init__(self, parts=None):
        self.parts = collections.OrderedDict()
        if parts is not None:
            for part in parts:
                self.parts[part.id] = part

    def getPart(self, id_):
        return self.parts[id_]

    def get_seq(self, id_, *args, **kwdargs):
        return self.parts[id_].get_seq(*args, **kwdargs)

    def __iter__(self):
        return iter(self.parts.values())
    def __len__(self):
        return len(self.parts)




class Segment(Locus):
    colors = {0:"red", 1:"blue", 2:"gray", 3:"orange", 4:"brown"}

    def __init__(self, chrom, start, end, strand, id_, source="genome"):
        self._chrom = chrom
        if start > end:
            start, end = end, start

        if start < 0:
            start = 0
        if end < 0:
            raise Exception("Segment end coordinate cannot be negative: {}".format(start))
        self._start = start
        self._end = end
        self._strand = strand
        self.id = id_
        self.source = source

    def __len__(self):
        return abs(self.end - self.start)

    def color(self):
        return self.colors[self.id]

    def antisense(self):
        antisense = {"+":"-", "-":"+"}
        return Segment(self.chrom, self.start, self.end, antisense[self.strand], self.id, self.source)

    def __repr__(self):
        return "<Segment {} {}:{}-{}{} (len={};{})>".format(self.id, self.chrom, self.start, self.end, self.strand, 
            len(self), self.source)
    


class StructuralVariant(object):
    def __init__(self, breakpoints, datahub, name):
        self.breakpoints = sorted(breakpoints, key=lambda x: (x.chrom, x.start))
        self.align_distance = datahub.align_distance

        self.sources = {"genome":datahub.genome}
        self.name = safe_file_name(name)

        self._seqs = {}

    def __getstate__(self):
        """ allows pickling of StructuralVariant()s """
        for allele in ["alt", "ref"]:
            for part in self.chrom_parts(allele):
                part.get_seq()
        state = self.__dict__.copy()
        return state


    def __str__(self):
        return "{}:{}({};{})".format(self.name, self.__class__.__name__, self.breakpoints, self.align_distance)
    def short_name(self):
        return "{}.{}_{}_{}".format(self.name, self.__class__.__name__[:3].lower(), self.breakpoints[0].chrom, self.breakpoints[0].start)

    def search_regions(self):
        pass    

    def seqs(self, allele):
        """ this is basically the 'reference' genome for the alternate allele """
        names_to_references = {}
        chrom_parts_collection = self.chrom_parts(allele)
        for name, chrom_part in chrom_parts_collection.parts.items():
            names_to_references[chrom_part.id] = chrom_part.get_seq()
        return names_to_references

    def chrom_parts(self, allele):
        """ overload this method for multi-part variants """
        segments = self.segments(allele)
        
        name = "{}_part".format(allele)
        if allele == "amb":
            name = "ref_part"

        parts = [ChromPart(name, segments, self.sources)]
        return ChromPartsCollection(parts)   


    def _segments(self, allele):
        segments = []
        for part in self.chrom_parts(allele):
            segments.extend(part.segments)
        return segments

    # def common_segments(self):
    #     """ return the segment IDs of the segments that are identical between 
    #     the ref and alt alleles (eg, flanking regions) """
    #     common = []
    #     refCounter = collections.Counter((segment.id for segment in self._segments("ref")))
    #     altCounter = collections.Counter((segment.id for segment in self._segments("alt")))
    #     if max(refCounter.values()) > 1 or max(altCounter.values()) > 1:
    #         logging.warn(" Same genomic region repeated multiple times within one allele; "
    #             "all flanking reads will be marked as ambiguous")
    #         return []


    #     refSegments = dict((segment.id, segment) for segment in self._segments("ref"))
    #     altSegments = dict((segment.id, segment) for segment in self._segments("alt"))

    #     for segmentID, refSegment in refSegments.items():
    #         if not segmentID in altSegments:
    #             continue
    #         altSegment = altSegments[segmentID]

    #         # Could remove the requirement to have the strand be the same
    #         # allowing the reads within the inversion to be plotted too
    #         if refSegment.chrom==altSegment.chrom and \
    #             refSegment.start == altSegment.start and \
    #             refSegment.end == altSegment.end and \
    #             refSegment.strand == altSegment.strand and \
    #             refSegment.source == altSegment.source:
    #             common.append(segmentID)

    #     return common


class SequenceDefinedVariant(StructuralVariant):
    def __init__(self, chrom, start, end, alt_seq, datahub, name):
        breakpoint = Locus(chrom, start, end, "+")
        super(SequenceDefinedVariant, self).__init__([breakpoint], datahub, name)

        self.sources["insertion"] = genomesource.GenomeSource({"insertion":alt_seq})
        self.insertionLength = len(alt_seq)

    def search_regions(self, searchDistance):
        chrom = self.breakpoints[0].chrom
        return [Locus(chrom, non_negative(self.breakpoints[0].start-searchDistance), 
                      self.breakpoints[-1].end+searchDistance, "+")]

    def segments(self, allele):
        breakpoint = self.breakpoints[0]
        chrom = breakpoint.chrom

        # If breakpoint has no length, we make the insertion before the breakpoint coordinate
        deletionOffset = 0
        if len(breakpoint) > 1:
            # If we're deleting some bases in addition to inserting, we'll make sure to start
            # the last segment after the deleted bases
            deletionOffset = 1

        if allele in ["ref", "amb"]:
            refSegments = []
            refSegments.append(Segment(chrom, breakpoint.start-self.align_distance, breakpoint.start-1, "+", 0))
            if len(breakpoint) > 1:
                refSegments.append(Segment(chrom, breakpoint.start, breakpoint.end, "+", 3))

            refSegments.append(Segment(chrom, breakpoint.end+deletionOffset, breakpoint.end+self.align_distance, "+", 2))
            return refSegments
        elif allele == "alt":
            return [Segment(chrom, breakpoint.start-self.align_distance, breakpoint.start-1, "+", 0),
                    Segment("insertion", 0, self.insertionLength, "+", 1, source="insertion"),
                    Segment(chrom, breakpoint.end+deletionOffset, breakpoint.end+self.align_distance, "+", 2)]

    def __str__(self):
        if len(self.breakpoints[0]) > 1:
            return "{}:{}::{}:{:,}-{:,};len={}".format(
                                                    self.name,
                                                    self.__class__.__name__, 
                                                    self.breakpoints[0].chrom, 
                                                    self.breakpoints[0].start, 
                                                    self.breakpoints[0].end,
                                                    self.insertionLength)

        return "{}:{}::{}:{:,};len={}".format(self.name, self.__class__.__name__,
            self.breakpoints[0].chrom, self.breakpoints[0].start, self.insertionLength)
       
    def short_name(self):
        return "{}.{}.{}_{}-{}".format(self.name, self.__class__.__name__,
            self.breakpoints[0].chrom, 
            self.breakpoints[0].start, 
            self.breakpoints[-1].end)


class Breakend(StructuralVariant):
    def __init__(self, breakpoint1, breakpoint2, datahub, name):
        super(Breakend, self).__init__([breakpoint1, breakpoint2], datahub, name)

        self.breakpoints = [breakpoint1, breakpoint2]
        self.chrom_parts("alt")

    def search_regions(self, searchDistance):
        search_regions = []

        for breakpoint in self.breakpoints:
            cur_locus = Locus(breakpoint.chrom,
                              non_negative(breakpoint.start-searchDistance), 
                              breakpoint.end+searchDistance,
                              breakpoint.strand)
            search_regions.append(cur_locus)

        return search_regions

    def chrom_parts(self, allele):
        b1 = self.breakpoints[0]
        b2 = self.breakpoints[1]

        segments = []
        for i, breakpoint in enumerate(self.breakpoints):
            segments.append(Segment(breakpoint.chrom, breakpoint.start-self.align_distance, 
                                    breakpoint.start, "+", 0+i*2))
            segments.append(Segment(breakpoint.chrom, breakpoint.start, 
                                    breakpoint.start+self.align_distance, "+", 1+i*2))

        parts = []
        if allele in ["ref", "amb"]:
            name = "ref_{}".format(b1.chrom)
            parts.append(ChromPart(name, [segments[0], segments[1]], self.sources))

            name = "ref_{}".format(b2.chrom)
            if b1.chrom == b2.chrom: name += "b"
            parts.append(ChromPart(name, [segments[2], segments[3]], self.sources))

        else:
            if b1.strand == "+": s1 = segments[0]
            else: s1 = segments[1].antisense()

            if b2.strand == "+": s2 = segments[3]
            else: s2 = segments[2].antisense()
            
            l1 = Locus(s1.chrom, s1.start, s1.end, "+")
            l2 = Locus(s2.chrom, s2.start, s2.end, "+")
            if l1.overlaps(l2) or l1.overlapsAntisense(l2):
                raise Exception("Not yet implemented - breakend-breakpoints near one another")

            # loci = [Locus(s.chrom, s.start, s.end, "+") for s in segments]
            # for i in range(len(loci)-1):
            #     for j in range(i+1, len(loci)):
            #         if loci[i].overlaps(loci[j]):
            #             raise Exception("Not yet implemented - breakend-breakpoints near one another")

            name = "alt_{}/{}".format(b1.chrom, b2.chrom)
            parts.append(ChromPart(name, [s1, s2], self.sources))

        return ChromPartsCollection(parts) 

    def __str__(self):
        chrom1 = self.breakpoints[0].chrom
        chrom2 = self.breakpoints[1].chrom
        if not chrom1.startswith("chr"):
            chrom1 = "chr{}".format(chrom1)
            chrom2 = "chr{}".format(chrom2)
        return "{}|{}::{}:{:,}/{}:{:,}".format(
            self.name, self.__class__.__name__, 
            chrom1, self.breakpoints[0].start, chrom2, self.breakpoints[1].start)


class Deletion(StructuralVariant):
    @classmethod
    def from_breakpoints(class_, chrom, first, second, datahub, name):
        breakpointLoci = [Locus(chrom, first, first, "+"), Locus(chrom, second, second, "+")]
        return class_(breakpointLoci, datahub, name)

    def search_regions(self, searchDistance):
        chrom = self.breakpoints[0].chrom
        deletionRegion = Locus(chrom, non_negative(self.breakpoints[0].start-searchDistance), 
            self.breakpoints[-1].end+searchDistance, "+")
        return [deletionRegion]

    def deletionLength(self):
        length = self.breakpoints[1].end - self.breakpoints[0].start
        return length

    def segments(self, allele):
        chrom = self.breakpoints[0].chrom

        if allele in ["ref", "amb"]:
            return [Segment(chrom, self.breakpoints[0].start-self.align_distance, self.breakpoints[0].start-1, "+", 0),
                    Segment(chrom, self.breakpoints[0].start, self.breakpoints[1].end, "+", 1),
                    Segment(chrom, self.breakpoints[1].end+1, self.breakpoints[1].end+self.align_distance, "+", 2)]
        elif allele == "alt":
            return [Segment(chrom, self.breakpoints[0].start-self.align_distance, self.breakpoints[0].start-1, "+", 0),
                    Segment(chrom, self.breakpoints[1].end+1, self.breakpoints[1].end+self.align_distance, "+", 2)]

    def __str__(self):
        return "{}::{}:{:,}-{:,}({})".format(self.__class__.__name__, self.breakpoints[0].chrom, self.breakpoints[0].start, 
            self.breakpoints[1].end, self.deletionLength())


# class Inversion(StructuralVariant):
#     def __init__(self, region, align_distance, fasta):
#         breakpoints = [Locus(region.chrom, region.start, region.start, "+"), Locus(region.chrom, region.end, region.end, "+")]
#         super(Inversion, self).__init__(breakpoints, align_distance, fasta)

#         self.region = region

#     def chrom(self):
#         return self.region.chrom

#     def search_regions(self, searchDistance):
#         chrom = self.chrom()

#         if len(self.region) < 2*searchDistance:
#             # return a single region
#             return [Locus(chrom, nonNegative(self.region.start-searchDistance), self.region.end+searchDistance, "+")]
#         else:
#             # return two regions, each around one of the ends of the inversion
#             search_regions = []
#             search_regions.append(Locus(chrom, nonNegative(self.region.start-searchDistance), 
#                 self.region.start+searchDistance, "+"))
#             search_regions.append(Locus(chrom, nonNegative(self.region.end-searchDistance), 
#                 self.region.end+searchDistance, "+"))
#             return search_regions

#     def segments(self, allele):
#         chrom = self.chrom()

#         if allele in ["ref", "amb"]:
#             return [Segment(chrom, self.region.start-self.align_distance, self.region.start-1, "+", 0),
#                     Segment(chrom, self.region.start, self.region.end, "+", 1),
#                     Segment(chrom, self.region.end+1, self.region.end+self.align_distance, "+", 2)]
#         elif allele == "alt":
#             return [Segment(chrom, self.region.start-self.align_distance, self.region.start-1, "+", 0),
#                     Segment(chrom, self.region.start, self.region.end, "-", 1),
#                     Segment(chrom, self.region.end+1, self.region.end+self.align_distance, "+", 2)]

                
#     def __str__(self):
#         return "{}::{}:{:,}-{:,}".format(self.__class__.__name__, self.region.chrom, self.region.start, self.region.end)


# class Insertion(StructuralVariant):
#     def __init__(self, breakpoint, insertSeq, align_distance, fasta):
#         super(Insertion, self).__init__([breakpoint], align_distance, fasta)
#         self.sources["insertion"] = genomesource.GenomeSource(insertSeq)
#         self.insertionLength = len(insertSeq)

#     def search_regions(self, searchDistance):
#         chrom = self.breakpoints[0].chrom
#         return [Locus(chrom, nonNegative(self.breakpoints[0].start-searchDistance), 
#             self.breakpoints[-1].end+searchDistance, "+")]

#     def segments(self, allele):
#         breakpoint = self.breakpoints[0]
#         chrom = breakpoint.chrom

#         # If breakpoint has no length, we make the insertion before the breakpoint coordinate
#         deletionOffset = 0
#         if len(breakpoint) > 1:
#             # If we're deleting some bases in addition to inserting, we'll make sure to start
#             # the last segment after the deleted bases
#             deletionOffset = 1

#         if allele in ["ref", "amb"]:
#             refSegments = []
#             refSegments.append(Segment(chrom, breakpoint.start-self.align_distance, breakpoint.start-1, "+", 0))
#             if len(breakpoint) > 1:
#                 refSegments.append(Segment(chrom, breakpoint.start, breakpoint.end, "+", 3))

#             refSegments.append(Segment(chrom, breakpoint.end+deletionOffset, breakpoint.end+self.align_distance, "+", 2))
#             return refSegments
#         elif allele == "alt":
#             return [Segment(chrom, breakpoint.start-self.align_distance, breakpoint.start-1, "+", 0),
#                     Segment("insertion", 0, self.insertionLength, "+", 1, source="insertion"),
#                     Segment(chrom, breakpoint.end+deletionOffset, breakpoint.end+self.align_distance, "+", 2)]

#     def __str__(self):
#         if len(self.breakpoints[0]) > 1:
#             return "{}::{}:{:,}-{:,};len={}".format(self.__class__.__name__, 
#                                                   self.breakpoints[0].chrom, 
#                                                   self.breakpoints[0].start, 
#                                                   self.breakpoints[0].end,
#                                                   self.insertionLength)

#         return "{}::{}:{:,};len={}".format(self.__class__.__name__, self.breakpoints[0].chrom, self.breakpoints[0].start, self.insertionLength)
       

# class MobileElementInsertion(StructuralVariant):
#     def __init__(self, breakpoint, insertedSeqLocus, insertionFasta, align_distance, refFasta):
#         super(MobileElementInsertion, self).__init__([breakpoint], align_distance, refFasta)

#         self.sources["repeats"] = insertionFasta
#         self.insertedSeqLocus = insertedSeqLocus


#     def search_regions(self, searchDistance):
#         chrom = self.breakpoints[0].chrom
#         return [Locus(chrom, nonNegative(self.breakpoints[0].start-searchDistance), 
#             self.breakpoints[-1].end+searchDistance, "+")]

#     def segments(self, allele):
#         chrom = self.breakpoints[0].chrom

#         if allele in ["ref", "amb"]:
#             return [Segment(chrom, self.breakpoints[0].start-self.align_distance, self.breakpoints[0].start-1, "+", 0),
#                     Segment(chrom, self.breakpoints[0].end, self.breakpoints[0].end+self.align_distance, "+", 2)]
#         elif allele == "alt":
#             return [Segment(chrom, self.breakpoints[0].start-self.align_distance, self.breakpoints[0].start-1, "+", 0),
#                     Segment(self.insertedSeqLocus.chrom, self.insertedSeqLocus.start, 
#                         self.insertedSeqLocus.end, self.insertedSeqLocus.strand, 1, source="repeats"),
#                     Segment(chrom, self.breakpoints[0].end, self.breakpoints[0].end+self.align_distance, "+", 2)]

#     def __str__(self):
#         return "{}::{}({});{})".format(self.__class__.__name__, self.insertedSeqLocus.chrom, self.breakpoints, len(self.insertedSeqLocus))
#     def short_name(self):
#         return "{}_{}_{}".format("mei", self.breakpoints[0].chrom, self.breakpoints[0].start)


# class Translocation(StructuralVariant):
#     def __init__(self, breakpoint1, breakpoint2, align_distance, refFasta):
#         super(Translocation, self).__init__([breakpoint1, breakpoint2], align_distance, refFasta)

#         self.breakpoints = [breakpoint1, breakpoint2]

#     def search_regions(self, searchDistance):
#         search_regions = []

#         for breakpoint in self.breakpoints:
#             search_regions.append(Locus(breakpoint.chrom, nonNegative(breakpoint.start-searchDistance), 
#                 breakpoint.end+searchDistance, breakpoint.strand))

#         return search_regions

#     def chrom_parts(self, allele):
#         parts = []
#         b1 = self.breakpoints[0]
#         b2 = self.breakpoints[1]

#         segments = []
#         for i, breakpoint in enumerate(self.breakpoints):
#             segments.append(Segment(breakpoint.chrom, breakpoint.start-self.align_distance, 
#                                     breakpoint.start-1, breakpoint.strand, 0+i*2))
#             segments.append(Segment(breakpoint.chrom, breakpoint.start, 
#                                     breakpoint.start+self.align_distance, breakpoint.strand, 1+i*2))
#             # assert breakpoint.strand == "+", breakpoint

#         if b2.strand == "-":
#             if allele in ["ref", "amb"]:
#                 name = "ref_{}".format(b1.chrom)
#                 parts.append(ChromPart(name, [segments[0], segments[1]], self.sources))

#                 name = "ref_{}".format(b2.chrom)
#                 if b1.chrom == b2.chrom: name += "b"
#                 parts.append(ChromPart(name, [segments[3], segments[2]], self.sources))


#             elif allele == "alt":
#                 name = "alt_{}/{}".format(b1.chrom, b2.chrom)
#                 parts.append(ChromPart(name, [segments[0], segments[2]], self.sources))

#                 name = "alt_{}/{}".format(b2.chrom, b1.chrom)
#                 if b1.chrom == b2.chrom: name += "b"

#                 parts.append(ChromPart(name, [segments[3], segments[1]], self.sources))
#         else:
#             if allele in ["ref", "amb"]:
#                 name = "ref_{}".format(b1.chrom)
#                 part = ChromPart(name, segments[:2], self.sources)
#                 parts.append(part)

#                 name = "ref_{}".format(b2.chrom)
#                 if b1.chrom == b2.chrom: name += "b"

#                 part = ChromPart(name, segments[2:], self.sources)
#                 parts.append(part)

#             elif allele == "alt":
#                 name = "alt_{}/{}".format(b1.chrom, b2.chrom)
#                 parts.append(ChromPart(name, [segments[0], segments[3]], self.sources))

#                 name = "alt_{}/{}".format(b2.chrom, b1.chrom)
#                 if b1.chrom == b2.chrom: name += "b"

#                 parts.append(ChromPart(name, [segments[2], segments[1]], self.sources))


#         return ChromPartsCollection(parts) 

#     def __str__(self):
#         chrom1 = self.breakpoints[0].chrom
#         chrom2 = self.breakpoints[1].chrom
#         if not chrom1.startswith("chr"):
#             chrom1 = "chr{}".format(chrom1)
#             chrom2 = "chr{}".format(chrom2)
#         return "{}::{}:{:,}/{}:{:,}".format(self.__class__.__name__, chrom1, self.breakpoints[0].start, chrom2, self.breakpoints[1].start)



# class LargeDeletion(Breakend):
#     @classmethod
#     def from_breakpoints(class_, chrom, first, second, align_distance, fasta):
#         breakpoint1 = Locus(chrom, first, first, "+")
#         breakpoint2 = Locus(chrom, second, second, "+")
#         return class_(breakpoint1, breakpoint2, align_distance, fasta)

#     def deletionLength(self):
#         return self.breakpoints[1].end - self.breakpoints[0].start

#     def __str__(self):
#         return "{}::{}:{:,}-{:,}({})".format(self.__class__.__name__, self.breakpoints[0].chrom, self.breakpoints[0].start, 
#             self.breakpoints[1].end, self.deletionLength())