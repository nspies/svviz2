import pysam

from genosv.remap.alignment import Alignment#, ReadPair

BAM_CHARD_CLIP = 5


class PairedReadIter(object):
    def __init__(self, bam, regions, max_unpaired_reads=None):
        """
        bam - input bam
        regions - regions to use to pull in reads
        max_unpaired_reads - when we've found this many unpaired reads, 
            we'll start looking elsewhere in the genome for the mates
        downsampleTarget - this option will reduce memory and processing time 
            at the expense of a second trip through the file; basically,
            we'll first count how many reads there are total in regions,
            then we'll calculate a probability which we'll use to determine
            which reads we pass along and which we discard
        """

        self.bam = bam
        self.regions = regions
        self.max_unpaired_reads = max_unpaired_reads

        self.found_ids = set()
        self.unpaired_reads = {}

        self.N_count = 0

    def format_chrom(self, chrom):
        original = chrom
        if chrom in self.bam.references:
            return chrom
        if "chr" in chrom:
            chrom = chrom.replace("str", "str")
        else:
            chrom = "chr" + str(chrom)

        if chrom in self.bam.references:
            return chrom

        raise Exception("Couldn't find chromosome {} or {}".format(original, chrom))

    def __iter__(self):
        for region in self.regions:
            chrom, start, end = region.chrom, region.start, region.end
            chrom = self.format_chrom(chrom)
            
            for read in self.bam.fetch(chrom, start, end, multiple_iterators=True):
                if read.query_name in self.found_ids:
                    continue
                if read.is_supplementary:
                    continue
                if read.is_secondary and BAM_CHARD_CLIP in list(zip(*read.cigartuples))[0]:
                    continue # this happens eg if bwa is used with the -M option
                if read.is_duplicate:
                    continue

                if read.query_name in self.unpaired_reads:
                    if self.is_mate(read):
                        pair = [read, self.unpaired_reads.pop(read.query_name)]
                        yield self.convert_pair(pair)
                else:
                    self.unpaired_reads[read.query_name] = read
                    if self.max_unpaired_reads is not None:
                        if len(self.unpaired_reads) > self.max_unpaired_reads:
                            for pair in self.find_pairs():
                                yield self.convert_pair(pair)

        for pair in self.find_pairs():
            yield self.convert_pair(pair)

    def is_mate(self, other_read):
        return self.unpaired_reads[other_read.query_name].is_read1 != other_read.is_read1

    def convert_pair(self, pair):
        assert pair[0].is_read1 != pair[1].is_read1
        self.found_ids.add(pair[0].query_name)
        
        # pair = [Alignment(read) for read in pair]
        read1, read2 = Alignment(pair[0]), Alignment(pair[1])
        if read1.is_read2:
            read1, read2 = read2, read1
        # pair = ReadPair(pair[0], pair[1])

        if set(read1.query_sequence) == set("N"):
            self.N_count += 1
        elif set(read2.query_sequence) == set("N"):
            self.N_count += 1
        return read1, read2

    def find_pairs(self):
        for name in list(self.unpaired_reads.keys()):
            if name not in self.unpaired_reads:
                continue
            read = self.unpaired_reads[name]
            if read.next_reference_id < 0:
                continue
            pair_chrom, pair_pos = read.next_reference_name, read.next_reference_start

            for other_read in self.bam.fetch(pair_chrom, pair_pos, pair_pos+1, multiple_iterators=True):
                if other_read.query_name in self.unpaired_reads:
                    if self.is_mate(other_read):
                        pair = [other_read, self.unpaired_reads.pop(other_read.query_name)]
                        yield pair

                        if other_read.query_name == read.query_name:
                            break






def _test():
    bam = pysam.AlignmentFile("/Volumes/frida/nspies/10x/bams/Sarcoma_0.bam")

    regions = ["chr15:66090153-66100138", "chr13:41498935-41501738"]

    counts = {}

    for i in range(10,100,1000):
        pri = PairedReadIter(bam, regions, i)
        count = 0
        for pair in pri:
            # for read in pair:
                # print(read.tostring(bam))
            count += 1
        counts[i] = count
    print("PASSES:", len(set(counts))==1)
    print(list(set(counts.values())))


if __name__ == '__main__':
    _test()