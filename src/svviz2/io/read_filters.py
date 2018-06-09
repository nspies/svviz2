def get_haplotype_filter(hp):
    def filter_haplotype(batch):
        filtered_batch = []
        for read_pair in batch:
            cur_hp = None
            cur_read1 = read_pair.original_read_ends["1"]
            if cur_read1.has_tag("HP"):
                cur_hp = cur_read1.get_tag("HP")
            if cur_hp == hp:
                filtered_batch.append(read_pair)
        return filtered_batch
    return filter_haplotype
 
    
