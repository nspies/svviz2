import os

from genosv.visualize import track
from genosv.io import export

# TODO: this is genosv.visualize.visualize.visualize -- awkward!
def visualize(datahub):#, temp_storage):
    for sample_name, sample in datahub.samples.items():
        sample.tracks = {}
        for allele in ["alt", "ref"]:
            # cur_alns = [aln.supporting_aln for aln in temp_storage[sample_name] if aln.supports_allele==allele]
            # bam = sample.out_alt_bam if allele == "alt" else sample.out_ref_bam
            bam = sample.outbam(allele, "r")
            sample.tracks[allele] = track.Track(
                datahub.variant.chrom_parts(allele), bam, 3000, 4000, datahub.variant, allele, False, True,
                (not sample.single_ended))

        for allele in ["alt", "ref"]:
            axis = track.Axis(sample.tracks[allele].scale, datahub.variant, allele)
            datahub.alleleTracks[allele]["axis"] = axis

    e = export.TrackCompositor(datahub)
    
    file_format = datahub.args.format
    converter = export.getExportConverter(file_format)

    outpath = os.path.join(datahub.args.outdir,
                           "{}.{}".format(datahub.variant.short_name(), file_format))
    mode = "w" if file_format=="svg" else "wb"

    with open(outpath, mode) as outf:
        d = export.convertSVG(e.render(), file_format, converter)
        outf.write(d)
