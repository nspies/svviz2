import logging
import os

from genosv.visualize import track
from genosv.io import export

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
