import collections
import itertools
import math
import numpy
import re
from genosv.visualize.svg import SVG
# from svviz import utilities
# from svviz import variants

class Scale(object):
    def __init__(self, chromPartsCollection, pixelWidth, dividerSize=25):
        # length is in genomic coordinates, starts is in pixels
        self.dividerSize = dividerSize
        self.partsToLengths = collections.OrderedDict()
        self.partsToStartPixels = collections.OrderedDict()
        self.chromPartsCollection = chromPartsCollection

        for part in chromPartsCollection:
            self.partsToLengths[part.id] = len(part)

        self.pixelWidth = pixelWidth

        totalLength = sum(self.partsToLengths.values()) + (len(self.partsToLengths)-1)*dividerSize
        self.basesPerPixel = totalLength / float(pixelWidth)

        curStart = 0
        for regionID in self.partsToLengths:
            self.partsToStartPixels[regionID] = curStart
            curStart += (self.partsToLengths[regionID]+dividerSize) / self.basesPerPixel

    def topixels(self, g, regionID=None):
        pts = 0
        if regionID != None:
            pts = self.partsToStartPixels[regionID]
        else:
            assert len(self.partsToStartPixels) == 1


        pos = g / float(self.basesPerPixel) + pts
        return pos

    def relpixels(self, g):
        dist = g / float(self.basesPerPixel)
        return dist

    def getBreakpointPositions(self, regionID):
        breakpoints = []

        part = self.chromPartsCollection.parts[regionID]

        curpos = 0
        for segment in part.segments[:-1]:
            curpos += len(segment)
            breakpoints.append(curpos)

        return breakpoints


class Axis(object):
    def __init__(self, scale, variant, allele):
        self.scale = scale
        self.allele = allele
        self.variant = variant
        self.chromPartsCollection = variant.chrom_parts(allele)
        self.height = 75


    def baseHeight(self):
        return 75

    def render(self, scaleFactor=1.0, spacing=1.0, height=None, thickerLines=False):
        self.height = height
        if height == None:
            self.height = 75 * scaleFactor

        self.svg = SVG(self.scale.pixelWidth, self.height, yrelto="top", headerExtras="""preserveAspectRatio="none" """)

        # dividers! (multi-part only)
        if len(self.chromPartsCollection) > 1:
            for part in list(self.chromPartsCollection)[1:]:
                divWidth = self.scale.relpixels(self.scale.dividerSize)
                x = self.scale.partsToStartPixels[part.id] - divWidth
                self.svg.rect(x, 0, divWidth, self.height, fill="#B2B2B2")
                self.svg.line(x, 0, x, self.height, stroke="black", **{"stroke-width":1*scaleFactor})
                self.svg.line(x+divWidth, 0, x+divWidth, self.height, stroke="black", **{"stroke-width":1*scaleFactor})

            for i, part in enumerate(self.chromPartsCollection):
                if i == 0:
                    x = self.scale.topixels(self.scale.partsToLengths[part.id], part.id) - 15
                    anchor = "end"
                else:
                    x = self.scale.partsToStartPixels[part.id] + 15
                    anchor = "start"
                self.svg.text(x, 18*scaleFactor, part.id, anchor=anchor, size=18*scaleFactor)



        # tick marks and coordinates
        for regionID in self.scale.partsToStartPixels:
            for tick in self.getTicks(regionID):
                x = self.scale.topixels(tick, regionID)

                self.svg.rect(x, 35*scaleFactor, 1*scaleFactor, 15*scaleFactor, fill="black")
                label = tick
                if tick > 1e6:
                    label = "{:.1f}MB".format(tick/1e6)
                elif tick > 1e3:
                    label = "{:.1f}KB".format(tick/1e3)

                if x < 50:
                    x = 50
                elif x > self.scale.pixelWidth - 50:
                    x = self.scale.pixelWidth - 50
                extras = {}
                if thickerLines:
                    extras["font-weight"] = "bold"
                self.svg.text(x, self.height-4*scaleFactor, label, size=18*scaleFactor, **extras)


        # segment arrows
        for part in self.chromPartsCollection:
            curOffset = self.scale.partsToStartPixels[part.id]
            for segment in part.segments:
                start = curOffset
                end = self.scale.relpixels(len(segment)) + curOffset
                curOffset = end

                arrowDirection = "right"

                if segment.strand == "-":
                    start, end = end, start
                    arrowDirection = "left"

                y = 35*scaleFactor

                self.svg.line(start, y, end, y, stroke=segment.color(), **{"stroke-width":8*scaleFactor})
                self.svg.lineWithInternalArrows(start, y, end, y, stroke=segment.color(), direction=arrowDirection,
                    arrowKwdArgs={"class":"scaleArrow"}, **{"stroke-width":3*scaleFactor})

        # breakpoints
        previousPosition = None
        for part in self.chromPartsCollection:
            for vline in self.scale.getBreakpointPositions(part.id):
                thickness = 1*scaleFactor
                if thickerLines:
                    thickness *= 2
                x = self.scale.topixels(vline, part.id)
                self.svg.line(x, 20*scaleFactor, x, 55*scaleFactor, 
                    stroke="black", **{"stroke-width":thickness})
                
                if previousPosition is None or x-previousPosition > 250*scaleFactor:     
                    self.svg.text(x-(scaleFactor/2.0), 18*scaleFactor, "breakpoint", size=18*scaleFactor, fill="black")
                previousPosition = x


        return str(self.svg)

    def getTicks(self, regionID):
        ticks = []
        start = 0
        width = self.scale.partsToLengths[regionID]
        end = start + width

        res = (10 ** round(math.log10(end - start))) / 10.0
        if width / res > 15:
            res *= 2.5
        elif width / res < 5:
            res /= 2.0

        roundStart = start - (start%res)

        for i in range(int(roundStart), end, int(res)):
            ticks.append(i)

        return ticks

NYIWarned = False
class ReadRenderer(object):
    def __init__(self, rowHeight, scale, chromPartsCollection, thickerLines=False, colorCigar=True, mismatch_counts=None):
        self.rowHeight = rowHeight
        self.svg = None
        self.scale = scale
        self.chromPartsCollection = chromPartsCollection
        self.thickerLines = thickerLines
        self.colorCigar = colorCigar

        self.nucColors = {"A":"blue", "C":"orange", "G":"green", "T":"black", "N":"gray"}
        self.colorsByStrand = {False:"purple", True:"red"}
        self.insertionColor = "cyan"
        self.deletionColor = "gray"
        self.overlapColor = "lime"

        self.mismatch_counts = mismatch_counts

    def render(self, alns, layout):
        # yoffset = alignmentSet.yoffset
        read_name = alns[0].query_name
        regionID = alns[0].reference_name

        yoffset = layout[read_name]

        if len(alns) > 1:
            assert alns[0].reference_id == alns[1].reference_id

        pstart = self.scale.topixels(alns[0].reference_start, regionID)
        pend = self.scale.topixels(alns[-1].reference_end, regionID)

        # isFlanking = (alignmentSet.parentCollection.why == "flanking")

        thinLineWidth = 5
        curColor = "#CCCCCC"
        extras = {}
        # if isFlanking:
        #     extras["class"] = "flanking"
        #     curColor = "#EEEEEE"

        self.svg.rect(pstart, yoffset-(self.rowHeight/2.0)+thinLineWidth/2.0, pend-pstart, thinLineWidth, fill=curColor, **extras)
        _overlap_start_x = pstart

        positionCounts = collections.Counter()

        if self.thickerLines:
            # extra "bold":
            ystart = yoffset+3
            height = self.rowHeight+6
        else:
            ystart = yoffset
            height = self.rowHeight


        # # for alignment in alignmentSet.getAlignments():
        # try:
        #     alns = [alignmentSet.aln1, alignmentSet.aln2]
        # except AttributeError:
        #     alns = [alignmentSet]

        for aln in alns:
            # alignment = end.locus
            # alignment.cigar = end.cigarstring
            # alignment.seq = end.query_sequence

            for position in range(aln.reference_start, aln.reference_end+1):
                positionCounts[position] += 1

            pstart = self.scale.topixels(aln.reference_start, regionID)
            pend = self.scale.topixels(aln.reference_end, regionID)

            curColor = self.colorsByStrand[aln.is_reverse]
            extras = {"class":"read"}#, "data-cigar":alignment.cigar,"data-readid":alignment.name}
            # if isFlanking:
            #     extras["class"] = "read flanking"
            #     curColor = "#AAAAAA"

            self.svg.rect(pstart, ystart, pend-pstart, height, fill=curColor, 
                          **extras)

            if self.colorCigar:
                self._drawCigar(aln, ystart, height)#, isFlanking)

        highlightOverlaps = True
        if highlightOverlaps:# and not isFlanking:
            self._highlightOverlaps(positionCounts, ystart, height, regionID, aln.query_name)#, isFlanking)

        self.svg.text(_overlap_start_x, yoffset, "{},{}".format(aln.mapq, aln.get_tag("OV")),
            anchor="end")


    def _drawCigar(self, alignment, yoffset, height):#, isFlanking):
        eachNuc = False # this gets to be computationally infeasible to display in the browser
        pattern = re.compile('([0-9]*)([MIDNSHP=X])')

        genomePosition = alignment.reference_start
        sequencePosition = 0

        region_id = alignment.reference_name
        chromPartSeq = self.chromPartsCollection.get_seq(region_id)

        extras = {}
        # if isFlanking:
            # extras = {"class":"flanking"}

        # for length, code in pattern.findall(alignment.cigar):
        alnseq = alignment.seq
        for code, length in alignment.cigartuples:
            length = int(length)
            if code == 0: #"M":
                for i in range(length):
                    curstart = self.scale.topixels(genomePosition+i, alignment.reference_name)
                    curend = self.scale.topixels(genomePosition+i+1, alignment.reference_name)

                    color = self.nucColors[alnseq[sequencePosition+i]]

                    alt = alnseq[sequencePosition+i]
                    ref = chromPartSeq[genomePosition+i]
                    
                    if eachNuc or alt!=ref:
                        if not self.mismatch_counts or alt=="N" or self.mismatch_counts.query(region_id, alt, genomePosition+i):
                            width = curend-curstart
                            width = max(width, self.scale.pixelWidth*1e-4)
                            midpoint = (curstart+curend)/2
                            # print(width, self.scale.pixelWidth, width/self.scale.pixelWidth)
                            self.svg.rect(midpoint-width/2, yoffset, width, height, fill=color, **extras)
                            # self.svg.rect(curstart, yoffset, curend-curstart, height, fill=color, **extras)

                sequencePosition += length
                genomePosition += length
            elif code == 2: #in "D":
                if not self.mismatch_counts or self.mismatch_counts.query(region_id, "DEL", genomePosition, genomePosition+length+1):
                    curstart = self.scale.topixels(genomePosition, alignment.reference_name)
                    curend = self.scale.topixels(genomePosition+length+1, alignment.reference_name)
                    self.svg.rect(curstart, yoffset, curend-curstart, height, fill=self.deletionColor, **extras)

                genomePosition += length
            elif code in [1, 4, 5]: #"IHS":
                # TODO: always draw clipping, irrespective of consensus sequence or mode
                if not self.mismatch_counts or self.mismatch_counts.query(region_id, "INS", genomePosition-2, genomePosition+2):
                    curstart = self.scale.topixels(genomePosition-0.5, alignment.reference_name)
                    curend = self.scale.topixels(genomePosition+0.5, alignment.reference_name)

                    width = curend-curstart
                    width = max(width, self.scale.pixelWidth*1e-4)
                    midpoint = (curstart+curend)/2

                    self.svg.rect(midpoint-width/2, yoffset, width, height, fill=self.insertionColor, **extras)
                    # self.svg.rect(curstart, yoffset, curend-curstart, height, fill=self.insertionColor, **extras)

                sequencePosition += length


    def _highlightOverlaps(self, positionCounts, yoffset, height, regionID, readID):#, isFlanking):
        overlapSegments = [list(i[1]) for i in itertools.groupby(sorted(positionCounts), lambda x: positionCounts[x]) if i[0] > 1]

        for segment in overlapSegments:
            start = min(segment)
            end = max(segment)

            curstart = self.scale.topixels(start, regionID)
            curend = self.scale.topixels(end, regionID)

            curColor = self.overlapColor
            # if isFlanking:
            #     curColor = "#88FF88"
            self.svg.rect(curstart, yoffset, curend-curstart, height, fill=curColor, 
                **{"class":"read", "data-readid":readID})

class MismatchCounts(object):
    """
    keeps track of how many of each nucleotide (or insertion/deletion) are present at each position
    -> used for the quick consensus mode
    """
    types_to_id = {"A":0, "C":1, "G":2, "T":3, "DEL":5}

    def __init__(self, chromPartsCollection):
        self.counts_by_region = {}
        self.insertions_by_region = {}
        for part in chromPartsCollection:
            self.counts_by_region[part.id] = numpy.zeros([6, len(part)])#, dtype="uint8")
            self.insertions_by_region[part.id] = numpy.zeros(len(part))#, dtype="uint8")

    def tally_reads(self, bam):
        depths = []
        for region_id in bam.references:
            for pileupcolumn in bam.pileup(region_id):
                depths.append(pileupcolumn.n)
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_refskip:
                        continue
                    elif pileupread.is_del:
                        self.add_count(region_id, pileupcolumn.pos, "DEL")
                    else:
                        nuc =  pileupread.alignment.query_sequence[pileupread.query_position]
                        if nuc != "N":
                            self.add_count(region_id, pileupcolumn.pos, nuc)
                    if pileupread.indel > 0:
                        self.add_count(region_id, pileupcolumn.pos, "INS")

    def add_count(self, region_id, position, type_):
        if type_ == "INS":
            self.insertions_by_region[region_id][position] += 1
        else:
            row = self.types_to_id[type_]
            self.counts_by_region[region_id][row,position] += 1

    def counts(self, region_id, position):
        return self.counts_by_region[region_id][position,:]

    def query(self, region_id, type_, start, end=None):
        if end is None:
            end = start

        total = self.counts_by_region[region_id][:,start:(end+1)].sum(axis=0)
        total = total.astype(float)

        if (type_ == "INS"):
            insertions = self.insertions_by_region[region_id][start:(end+1)]
            if (insertions.sum() / total.sum()) > 0.2:
               # print("INSERTION!", region_id, type_, start, end, insertions, total)
               return True

        else:        
            row = self.types_to_id[type_]
            this_type = self.counts_by_region[region_id][row,start:(end+1)]

            if type_ == "DEL":
                if  ((this_type / total) > 0.3).any():
                # if type_ == "DEL":
                #     print("DEL", region_id, type_, start, end, this_type, total)
                    return True
            elif ((this_type / total) > 0.2).any():
                return True

        return False


class Track(object):
    def __init__(self, chromPartsCollection, bam, height, width, variant, allele, thickerLines, colorCigar, paired,
                 quick_consensus=True):
        self.chromPartsCollection = chromPartsCollection
        self.height = height
        self.width = width

        self.scale = Scale(chromPartsCollection, width)

        self.rowHeight = 5
        self.rowMargin = 1

        self.bam = bam

        self.svg = None
        self.rendered = None

        self.variant = variant
        self.allele = allele

        self.rows = []
        self._axis = None

        self.xmin = None
        self.xmax = None

        self.thickerLines = thickerLines

        self.layout = {}
        self.paired = paired

        self.mismatch_counts = MismatchCounts(chromPartsCollection) if quick_consensus else None
        self.readRenderer = ReadRenderer(self.rowHeight, self.scale, self.chromPartsCollection,
                                         thickerLines, colorCigar, self.mismatch_counts)

    def findRow(self, start, end, regionID):
        for currow in range(len(self.rows)):
            if self.rows[currow] is None or (start - self.rows[currow]) >= 2:
                self.rows[currow] = end
                break
        else:
            self.rows.append(end)
            currow = len(self.rows)-1

        return currow

    def dolayout(self):
        self.rows = [None]

        self.xmin = 1e100
        self.xmax = 0

        for regionID in self.bam.references:
            cur_read_coords = collections.defaultdict(list)

            for read in self.bam.fetch(reference=regionID):
                cur_read_coords[read.query_name].append((read.reference_start, read.reference_end))

            for read_name, coords in cur_read_coords.items():
                start = coords[0][0]
                end = coords[-1][1]

                start = self.scale.topixels(start, regionID)
                end = self.scale.topixels(end, regionID)

                currow = self.findRow(start, end, regionID)
                yoffset = (self.rowHeight+self.rowMargin) * currow
                self.layout[read_name] = yoffset

                self.xmin = min(self.xmin, start)
                self.xmax = max(self.xmax, end)

        self.height = (self.rowHeight+self.rowMargin) * len(self.rows)

        if self.mismatch_counts is not None: # quick consensus
            print("using quick consensus mode; tallying...")
            self.mismatch_counts.tally_reads(self.bam)
            print("tallying done.")

    def render(self):        
        # if len(self.getAlignments()) == 0:
        if self.bam.count() == 0:
            xmiddle = self.scale.pixelWidth / 2.0

            self.height = xmiddle/20.0

            self.svg = SVG(self.width, self.height)
            self.svg.text(xmiddle, self.height*0.05, "No reads found", size=self.height*0.9, fill="#999999")
            self.rendered = self.svg.asString()
            return self.rendered

        self.dolayout()

        self.svg = SVG(self.width, self.height)
        self.readRenderer.svg = self.svg

        lineWidth = 1 if not self.thickerLines else 3
        lineWidth = lineWidth * ((self.xmax-self.xmin)/1200.0)

        for part in self.chromPartsCollection:
            for vline in self.scale.getBreakpointPositions(part.id):
                x = self.scale.topixels(vline, part.id)
                y1 = -20
                y2 = self.height+20
                self.svg.line(x, y1, x, y2, stroke="black", **{"stroke-width":lineWidth})


        # dividers!
        for part in list(self.chromPartsCollection)[1:]:
            divWidth = self.scale.relpixels(self.scale.dividerSize)
            x = self.scale.partsToStartPixels[part.id] - divWidth
            self.svg.rect(x, self.height+40, divWidth, self.height+40, fill="#B2B2B2")
            self.svg.line(x, 0, x, self.height+40, stroke="black", **{"stroke-width":4})
            self.svg.line(x+divWidth, 0, x+divWidth, self.height+40, stroke="black", **{"stroke-width":4})

        self.svg.rect(0, self.svg.height+20, self.scale.pixelWidth, 
            self.height+40, opacity=0.0, zindex=0)
        self.rendered = str(self.svg)

        # alnSets = self.getAlignments()
        # flankingAlignments = [alnSet for alnSet in alnSets if alnSet.parentCollection.why == "flanking"]
        # nonFlankingAlignments = [alnSet for alnSet in alnSets if alnSet.parentCollection.why != "flanking"]

        # for alignmentSet in alnSets:#flankingAlignments+nonFlankingAlignments:
            # self.readRenderer.render(alignmentSet)

        for regionID in self.bam.references:
            read_buffer = {}
            for read in self.bam.fetch(reference=regionID):
                if self.paired and read.query_name in read_buffer:
                    other_read = read_buffer.pop(read.query_name)
                    cur_reads = [other_read, read]
                elif self.paired:
                    read_buffer[read.query_name] = read
                    continue
                else:
                    cur_reads = [read]

                self.readRenderer.render(cur_reads, self.layout)

        return self.rendered


class AnnotationTrack(object):
    def __init__(self, annotationSet, scale, variant, allele):
        self.chromPartsCollection = variant.chromParts(allele)
        self.annotationSet = annotationSet
        self.scale = scale
        self.height = None
        self.variant = variant
        self.allele = allele

        self._annos = None
        self.rows = [None]
        self.svg = None

        self.rowheight = 20

    def _topixels(self, gpos, segment, psegoffset):
        if segment.strand == "+":
            pos = self.scale.relpixels(gpos - segment.start) + psegoffset
        elif segment.strand == "-":
            pos = self.scale.relpixels(segment.end - gpos) + psegoffset
        return pos

    def findRow(self, start, end):
        for currow in range(len(self.rows)):
            if self.rows[currow] is None or (start - self.rows[currow]) >= 2:
                self.rows[currow] = end
                break
        else:
            self.rows.append(end)
            currow = len(self.rows)-1

        return currow

    def baseHeight(self):
        if self._annos is not None and len(self._annos) == 0:
            return 0
        return ((len(self.rows)+2) * self.rowheight) + 20

    def dolayout(self, scaleFactor, spacing):
        # coordinates are in pixels not base pairs
        self.rows = [None]

        self._annos = []
        for part in self.chromPartsCollection:
            segmentStart = self.scale.partsToStartPixels[part.id]

            for segment in variants.mergedSegments(part.segments):
                curWidth = len(segment)
                curAnnos = self.annotationSet.getAnnotations(segment.chrom, segment.start, segment.end, clip=True)
                if segment.strand == "-":
                    curAnnos = sorted(curAnnos, key=lambda x:x.end, reverse=True)

                for anno in curAnnos:
                    start = max(anno.start, segment.start)
                    end = min(anno.end, segment.end)

                    start = self._topixels(start, segment, segmentStart)
                    end = self._topixels(end, segment, segmentStart)
                    if end < start:
                        start, end = end, start
                    textLength = len(anno.label)*self.rowheight/1.0*scaleFactor*spacing
                    rowNum = self.findRow(start, end+textLength)

                    anno.coords = {}
                    anno.coords["row"] = rowNum
                    anno.coords["start"] = start
                    anno.coords["end"] = end
                    anno.coords["strand"] = anno.strand if segment.strand=="+" else utilities.switchStrand(anno.strand)
                    anno.coords["segment"] = segment
                    anno.coords["segmentStart"] = segmentStart

                    self._annos.append(anno)

                segmentStart += self.scale.relpixels(curWidth)

    def drawBox(self, start, end, segment, segmentStart, y, height, scaleFactor, color, anno):
        start = self._topixels(start, segment, segmentStart)
        end = self._topixels(end, segment, segmentStart)

        if end < start:
            start, end = end, start
            
        width = end - start

        ystart = y-((self.rowheight-height)/2.0)*scaleFactor
        self.svg.rect(start, ystart, width, height*scaleFactor, fill=color)



    def _drawGenes(self, scaleFactor):
        for anno in self._annos:
            anno.txExons # check to see if this is a gff or a bed

            color = "blue" if anno.coords["strand"] == "+" else "darkorange"

            y = ((anno.coords["row"]+1) * self.rowheight + 20) * scaleFactor
            width = anno.coords["end"] - anno.coords["start"]

            self.svg.rect(anno.coords["start"], y-((self.rowheight-2.0)/2.0)*scaleFactor, width, 2.0*scaleFactor, fill=color)     

            for txExon in anno.txExons:
                start, end = txExon
                self.drawBox(start, end, anno.coords["segment"], anno.coords["segmentStart"], y, 5, scaleFactor, color, anno)

            for cdExon in anno.cdExons:
                start, end = cdExon
                self.drawBox(start, end, anno.coords["segment"], anno.coords["segmentStart"], y, self.rowheight, scaleFactor, color, anno)

            textSize = (self.rowheight-2)*scaleFactor
            if anno.coords["end"] + len(anno.label)*textSize*0.70 > self.scale.pixelWidth:
                w = len(anno.label)*textSize*0.70
                self.svg.rect(self.scale.pixelWidth - w, y, w, self.rowheight*scaleFactor, fill="white", **{"fill-opacity":0.40})
                self.svg.text(self.scale.pixelWidth, y-((self.rowheight-1)*scaleFactor), anno.label, size=textSize, fill=color, anchor="end", **{"fill-opacity":0.70})
            else:
                self.svg.text(anno.coords["end"]+(self.rowheight/2.0), y-((self.rowheight-1)*scaleFactor), 
                    anno.label, size=textSize, anchor="start", fill=color)   

    def _drawBED(self, scaleFactor):
        for anno in self._annos:
            color = "blue" if anno.coords["strand"] == "+" else "darkorange"
            y = ((anno.coords["row"]+1) * self.rowheight + 20) * scaleFactor
            width = anno.coords["end"] - anno.coords["start"]

            self.svg.rect(anno.coords["start"], y, width, self.rowheight*scaleFactor, fill=color)
            self.svg.text(anno.coords["end"]+(self.rowheight/2.0), y-((self.rowheight-1)*scaleFactor), 
                anno.label, size=(self.rowheight-2)*scaleFactor, anchor="start", fill=color)

    def render(self, scaleFactor=1.0, spacing=1, height=None, thickerLines=False):
        self.dolayout(scaleFactor, spacing)

        self.height = self.baseHeight()*scaleFactor
        self.svg = SVG(self.scale.pixelWidth, self.height)

        # dividers! (multi-part only)
        for part in list(self.chromPartsCollection)[1:]:
            divWidth = self.scale.relpixels(self.scale.dividerSize)
            x = self.scale.partsToStartPixels[part.id] - divWidth
            self.svg.rect(x, self.height, divWidth, self.height, fill="#B2B2B2")
            self.svg.line(x, 0, x, self.height, stroke="black", **{"stroke-width":1*scaleFactor})
            self.svg.line(x+divWidth, 0, x+divWidth, self.height, stroke="black", **{"stroke-width":1*scaleFactor})

        try:
            self._drawGenes(scaleFactor)
        except Exception as e:
            self._drawBED(scaleFactor)


        for part in self.chromPartsCollection:
            for vline in self.scale.getBreakpointPositions(part.id):
                x = self.scale.topixels(vline, part.id)-scaleFactor/2.0
                y1 = 0
                y2 = self.height
                thickness = 1*scaleFactor
                if thickerLines:
                    thickness *= 2
                self.svg.line(x, y1, x, y2, stroke="black", **{"stroke-width":thickness})


