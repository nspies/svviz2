def get_haplotype(read_or_pair):
    try:
        read = read_or_pair.original_read_ends["1"]
    except AttributeError:
        read = read_or_pair
        
    cur_hp = None
    if read.has_tag("HP"):
        cur_hp = read.get_tag("HP")
        
    return cur_hp
    
def get_haplotype_filter(hp):
    def filter_haplotype(batch):
        filtered_batch = []
        for read_or_pair in batch:
            if get_haplotype(read_or_pair) == hp:
                filtered_batch.append(read_or_pair)
        return filtered_batch
    return filter_haplotype
 
    
