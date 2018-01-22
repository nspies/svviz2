import logging
import os
import subprocess
import sys
import tempfile

from svviz2.utility import misc


def export(sv_document, datahub, context=None):
    file_format = datahub.args.format
    converter = getExportConverter(file_format)

    context_label = ".zoomed{}".format(context) if context else ""

    file_name = "{}{}".format(datahub.variant.short_name(), context_label)

    if file_format == "svg":
        outpath = os.path.join(datahub.args.outdir, file_name+".svg")
    else:
        temp_outdir = os.path.join(datahub.args.outdir, "svviz2-temp")
        misc.ensure_dir(temp_outdir)
        outpath = os.path.join(temp_outdir, file_name+".svg")

        final_outpath = os.path.join(datahub.args.outdir, file_name+"." + file_format)

    with open(outpath, "w") as f:
        for l in sv_document.render():
            f.write(l+"\n")

    if file_format != "svg":
        convertSVG(outpath, final_outpath, file_format, converter)



def getExportFormat(args):
    formats = [None, "png", "pdf", "svg"]
    if args.type == "batch" or args.format is not None:
        exportFormat = args.format
        if exportFormat is None:
            exportFormat = "pdf"
    else:
        exportFormat = args.export.partition(".")
        if len(exportFormat[2]) > 0:
            exportFormat = exportFormat[2]
            if exportFormat not in formats:
                logging.warn("= File suffix {} not recognized; exporting as .svg =".format(exportFormat))
                exportFormat = "svg"
        else:
            exportFormat = "svg"

    exportFormat = exportFormat.lower()
    return exportFormat

def getExportConverter(exportFormat, requested_converter=None):
    if requested_converter == "webkittopdf" and exportFormat=="png":
        logging.error("webkitToPDF does not support export to PNG; use librsvg or inkscape instead, or "
            "export to PDF")
        sys.exit(1)

    if exportFormat == "png" and requested_converter is None:
        return "librsvg"

    if requested_converter == "rsvg-convert":
        return "librsvg"

    if requested_converter in [None, "webkittopdf"]:
        if checkWebkitToPDF():
            return "webkittopdf"

    if requested_converter in [None, "librsvg"]:
        if checkRSVGConvert():
            return "librsvg"

    if requested_converter in [None, "inkscape"]:
        if checkInkscape():
            return "inkscape"

    return None



def checkWebkitToPDF():
    try:
        subprocess.check_call("webkitToPDF", stderr=subprocess.PIPE, shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

def checkRSVGConvert():
    try:
        subprocess.check_call("rsvg-convert -v", stdout=subprocess.PIPE, shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

def checkInkscape():
    try:
        subprocess.check_call("inkscape --version", stdout=subprocess.PIPE, shell=True)
        return True
    except subprocess.CalledProcessError:
        return False


def convertSVG(inpath, outpath, outformat, converter):
    if converter == "webkittopdf":
        exportData = _convertSVG_webkitToPDF(inpath, outpath, outformat)
    elif converter == "librsvg":
        exportData = _convertSVG_rsvg_convert(inpath, outpath, outformat)
    elif converter == "inkscape":
        exportData = _convertSVG_inkscape(inpath, outpath, outformat)

    return exportData

def _convertSVG_webkitToPDF(inpath, outpath, outformat):
    if outformat.lower() != "pdf":
        return None

    try:
        cmd = "webkitToPDF {} {}".format(inpath, outpath)
        subprocess.check_call(cmd, shell=True)#, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        return None

    return open(outpath, "rb").read()

def _convertSVG_inkscape(inpath, outpath, outformat):
    options = ""
    outformat = outformat.lower()
    if outformat == "png":
        options = "--export-dpi 150 --export-background white"

    try:
        subprocess.check_call("inkscape {} {} --export-{}={}".format(options, inpath, outformat, outpath), 
            shell=True)
    except subprocess.CalledProcessError as e:
        print("EXPORT ERROR:", str(e))

    return open(outpath, "rb").read()


def _convertSVG_rsvg_convert(inpath, outpath, outformat):
    options = ""
    outformat = outformat.lower()
    if outformat == "png":
        options = "-a --background-color white"

    try:
        subprocess.check_call("rsvg-convert -f {} {} -o {} {}".format(outformat, options, outpath, inpath), shell=True)
    except subprocess.CalledProcessError as e:
        print("EXPORT ERROR:", str(e))

    return open(outpath, "rb").read()


# def test():
#     base = """  <svg><rect x="10" y="10" height="100" width="100" style="stroke:#ffff00; stroke-width:3; fill: #0000ff"/><text x="25" y="25" fill="blue">{}</text></svg>"""
#     svgs = [base.format("track {}".format(i)) for i in range(5)]

#     tc = TrackCompositor(200, 600)
#     for i, svg in enumerate(svgs):
#         tc.addTrack(svg, i, viewbox="0 0 110 110")

#     outf = open("temp.svg", "w")
#     outf.write(tc.render())
#     outf.flush()
#     outf.close()

#     pdfPath = convertSVGToPDF("temp.svg")
#     subprocess.check_call("open {}".format(pdfPath), shell=True)

# if __name__ == '__main__':
#     test()

#     import sys
#     print(canConvertSVGToPDF(), file=sys.stderr)