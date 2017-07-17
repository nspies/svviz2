import logging
import os

from genosv.visualize import track
from genosv.io import export
from genosv.app import genomesource

logger = logging.getLogger(__name__)


# TODO: this is genosv.visualize.visualize.visualize -- awkward!
def visualize(datahub):#, temp_storage):
    if datahub.args.also_plot_context:
        logger.info("Now plotting zoomed-in")
        _visualize(datahub, context=datahub.args.also_plot_context)

    if datahub.args.only_plot_context:
        _visualize(datahub, context=datahub.args.only_plot_context)
    else:
        _visualize(datahub)


def _visualize(datahub, context=None):
    from genomeview import Document, ViewRow, GenomeView
    from genomeview.axis import Axis

    doc = Document(1500)
    for allele in ["alt", "ref"]:
        row = ViewRow(allele)

        for part in datahub.variant.chrom_parts(allele):
            source = genomesource.GenomeSource({part.id:part.get_seq()})
            genome_view = GenomeView(part.id, part.id, 0, len(part), "+", source)

            for sample_name, sample in datahub.samples.items():
                bam_path = sample.outbam_paths[allele].replace(".bam", ".sorted.bam") ## XXX todo: refactor
                bam_track = SVPairedEndBAMTrack(bam_path, part.segments)
                genome_view.add_track(bam_track)

            # axis = Axis(part.id)
            axis = ChromSegmentAxis(part.id, part.segments)
            genome_view.add_track(axis)
            row.add_view(genome_view)

        doc.elements.append(row)

    with open("temp.svg", "w") as f:
        for l in doc.render():
            f.write(l+"\n")
        #     for sample_name, sample in datahub.samples.items():
        # sample.tracks = {}




# def _visualize(datahub, context=None):
    for sample_name, sample in datahub.samples.items():
        sample.tracks = {}
        for allele in ["alt", "ref"]:
            bam = sample.outbam(allele, "r")
            sample.tracks[allele] = track.Track(
                datahub.variant.chrom_parts(allele), bam, 3000, 4000, datahub.variant, allele, False, True,
                (not sample.single_ended), (not sample.sequencer=="illumina"), zoomed=bool(context))

        for allele in ["alt", "ref"]:
            axis = track.Axis(sample.tracks[allele].scale, datahub.variant, allele, zoomed=bool(context))
            datahub.alleleTracks[allele]["axis"] = axis


            simple_repeats_track = track.SimpleRepeatsTrack(
                sample.tracks[allele].scale, datahub.variant, allele, zoomed=bool(context))
            datahub.alleleTracks[allele]["Simple Repeats"] = simple_repeats_track


    e = export.TrackCompositor(datahub, zoomed=context)
    
    file_format = datahub.args.format
    converter = export.getExportConverter(file_format)

    context_label = ".zoomed{}".format(context) if context else ""

    outpath = os.path.join(datahub.args.outdir,
                           "{}{}.{}".format(datahub.variant.short_name(), context_label, file_format))
    mode = "w" if file_format=="svg" else "wb"

    with open(outpath, mode) as outf:
        d = export.convertSVG(e.render(), file_format, converter)
        outf.write(d)

from genomeview.axis import Axis, get_ticks
from genomeview.bamtrack import PairedEndBAMTrack

class SVPairedEndBAMTrack(PairedEndBAMTrack):
    def __init__(self, bam_path, segments):
        super().__init__(bam_path)
        self.segments = segments

    def render(self, renderer):
        cur_coord = 0
        for segment in self.segments[:-1]:
            cur_coord += len(segment)
            cur_pos = self.scale.topixels(cur_coord)
            yield from renderer.line(cur_pos, 0, cur_pos, self.height)

        yield from super().render(renderer)

class ChromSegmentAxis(Axis):
    def __init__(self, name, chrom_segments):
        super().__init__(name)
        self.chrom_segments = chrom_segments
        self.height

    def render(self, renderer):
        start, end = self.scale.start, self.scale.end
        # yield from renderer.line(self.scale.topixels(start), 5, self.scale.topixels(end), 5)
        
        cur_start = 0
        arrow_size = 20
        arrow_spacing = 30
        y = arrow_size/2
        bar_height = 5

        ticks = get_ticks(start, end, max(1,self.scale.pixel_width/120))
        for i, (tick, label) in enumerate(ticks):
            x = self.scale.topixels(tick)
            if x < 0 or x > self.scale.pixel_width: continue

            yield from renderer.line(x, y, x, 20)
            anchor = "middle"
            if x < 50 and i == 0:
                anchor = "start"
                x = min(x, 5)
            elif x > self.scale.pixel_width-50 and i == len(ticks)-1:
                anchor = "end"
                x = max(x, self.scale.pixel_width-5)


            yield from renderer.text(x, y+30, label, anchor=anchor, size=16)



        for segment in self.chrom_segments:
            cur_end = cur_start + len(segment)
            left = self.scale.topixels(cur_start)
            width = self.scale.topixels(cur_end) - left
            yield from renderer.rect(left, y-bar_height/2, width, bar_height,
                                     fill=segment.color(), stroke="none")

            direction = "right" if segment.strand=="+" else "left"
            x = 0
            while x < width:
                yield from renderer.arrow(left+x+arrow_size, y, direction,
                    color=segment.color(), scale=arrow_size/10)
                x += arrow_size+arrow_spacing

            cur_start += len(segment)

