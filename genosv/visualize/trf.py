"""
A tandem-repeats finder (trf) wrapper
"""

import collections
import logging
import os
import subprocess
import tempfile

logger = logging.getLogger(__name__)

ALREADY_WARNED_TRF_MISSING = False

def run_trf(seqs):
    temp_dir = tempfile.mkdtemp()

    fasta_path = os.path.join(temp_dir, "seqs.fa")
    with open(fasta_path, "w") as fasta_file:
        for name, seq in seqs.items():
            fasta_file.write(">"+name+"\n"+seq+"\n")

    # the important parameter appears to be the "minscore" parameter,
    # which we'll set to 30 so that we can find repeats of length 15 and up
    # (with a match score of 2)
    trf_command = "trf {} 2 5 7 80 10 30 200 -l 6 -h -ngs".format(fasta_path)

    trf_proc = subprocess.Popen(trf_command, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result_code = trf_proc.wait()

    if result_code != 0:
        logger.warn("Failed to run trf (tandem repeat finder); please install from "
                    "http://tandem.bu.edu/trf/trf.download.html, rename binary to 'trf' and "
                    "move into your $PATH")
        return

    seqs_to_trfs = collections.defaultdict(list)
    cur_seq = None
    for line in trf_proc.stdout:
        line = line.decode("utf8").strip()
        if line.startswith("@"):
            cur_seq = line[1:]
        else:
            fields = line.split()
            start, end = int(fields[0])-1, int(fields[1])-1 # convert to zero-based coords, inclusive
            repeat_unit = fields[13]
            seqs_to_trfs[cur_seq].append((start, end, repeat_unit))
            # print(seqs[cur_seq][start:end+1], repeat_unit)
    # print(seqs_to_trfs)

    return seqs_to_trfs

def test():
    seqs = {"a":"CGACGACGACGACGA"}
    run_trf(seqs)

if __name__ == '__main__':
    test()