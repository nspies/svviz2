from genosv.visualize import track
from genosv.io import export

# TODO: this is genosv.visualize.visualize.visualize -- awkward!
def visualize(datahub, temp_storage):
    for sample_name, sample in datahub.samples.items():
        sample.tracks = {}
        for allele in ["alt", "ref"]:
            cur_alns = [aln.supporting_aln for aln in temp_storage[sample_name] if aln.supports_allele==allele]
            sample.tracks[allele] = track.Track(
                datahub.variant.chrom_parts(allele), cur_alns, 3000, 4000, datahub.variant, allele, False, True)

        for allele in ["alt", "ref"]:
            axis = track.Axis(sample.tracks[allele].scale, datahub.variant, allele)
            datahub.alleleTracks[allele]["axis"] = axis

    e = export.TrackCompositor(datahub)
    
    converter = export.getExportConverter("pdf")

    with open("{}.pdf".format(datahub.variant.short_name()), "wb") as outf:
        d = export.convertSVG(e.render(), "pdf", converter)
        outf.write(d)
