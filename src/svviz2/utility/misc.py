import errno
import os
import string


def ensure_dir(directory):
    try:
        os.makedirs(directory)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise


def safe_file_name(text):
    text = text.replace(":", "_")
    whitelist = "-_.()" + string.ascii_letters + string.digits

    return "".join(t for t in text if t in whitelist)



try:
    comp = str.maketrans('ATCGNatcgn','TAGCNtagcn')
except AttributeError:
    comp = string.maketrans('ATCGNatcgn','TAGCNtagcn')
    
def reverse_comp(st):
    """ Returns the reverse complement of a DNA sequence; non ACGT bases will be ignored. """
    return str(st)[::-1].translate(comp)


def str_to_bool(s):
    return s.lower() in ["1", "true", "t", "yes","y"]



def match_chrom_format(chrom, keys):
    if chrom in keys:
        return chrom
    if "chr" in chrom:
        chrom2 = chrom.replace("chr", "")
    else:
        chrom2 = "chr{}".format(chrom)

    if chrom2 in keys:
        return chrom2
    return chrom