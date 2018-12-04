
def map_realign(batch, realigner, sample):
    # if sample.single_ended:
    #     return map_realign_unpaired(batch, realigner)
    # else:
    return map_realign_pairs(batch, realigner, sample)

def map_realign_pairs(batch, datahub, sample):
    ref_genome_sources = [datahub.local_ref_genome_source]
    alt_genome_sources = [datahub.local_alt_genome_source]

    if datahub.aligner_type == "bwa" and not datahub.args.only_realign_locally:
        ref_genome_sources.append(datahub.genome)

    for genome_source in ref_genome_sources+alt_genome_sources:
        genome_source.set_aligner_params(sample.sequencer, sample.max_base_quality)

    import tqdm
    # for read_or_pair in batch:
    for read_or_pair in tqdm.tqdm(batch):
        # print(read_or_pair)
        #if read_or_pair.name == "D00360:64:HBAP3ADXX:1:2114:9685:53802":
        read_or_pair.realign(ref_genome_sources, alt_genome_sources)

    return batch

# def map_realign_pairs(batch, realigner, sample):
#     results = []
#     for pair in batch:
#         realigned_pair = realigner.realign_pair(pair, sample.read_statistics)
#         results.append(realigned_pair)

#     return results

# def map_realign_unpaired(batch, realigner):
#     import tqdm
#     results = []
#     # for read in tqdm.tqdm(batch):
#     for read in batch:
#         realigned = realigner.realign_single_end(read)
#         results.append(realigned)

#     return results
