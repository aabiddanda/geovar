
import allel
import pandas as pd
import numpy as np 

def vcf_to_freq_table(vcf_file, pop_panel, minor_allele=True):
    """
        Converts a VCF File to a frequency table 
    """
    # NOTE: we only assume that we have two columns in hte 
    pop_df = pd.read_table(pop_panel, usecols=['sample','pop'])
    pop_dict = pop_df.set_index(['sample']).to_dict()['pop']     
    # NOTE: we are assuming that we have only biallelic markers
    vcf_data = allel.read_vcf(vcf_file, alt_number=1)
    samples = vcf_data['samples']
    chrom = vcf_data['variants/CHROM']
    pos = vcf_data['variants/POS']
    ref_alleles = vcf_data['variants/REF']
    alt_alleles = vcf_data['variants/ALT']
    # Generating the population dictionary 
    unique_pops = np.unique([pop_dict[i] for i in pop_dict])
    pop_vector = np.array([pop_dict[i] for i in samples])
    pop_idx_dict = {}
    for p in unique_pops:
        pop_idx_dict[p] = np.where(pop_vector == p)[0]
    # Calculating the Alternative Allele Frequency
    gt = allel.GenotypeArray(vcf_data['calldata/GT'])
    # Calculating the total allele frequency 
    tot_ac_cnt = gt.count_alleles()
    alt_freq = tot_ac_cnt[:,1] / tot_ac_cnt.sum(axis=1)
    global_ac = tot_ac_cnt[:,1]
    global_af = alt_freq
    flip_af = np.repeat(True, alt_freq.size)
    if minor_allele:
        flip_af = (alt_freq > 0.5)
    # Generate the allele freqquency and allele count vectors 
    flip_idx = np.where(flip_af)[0]
    if flip_idx.size > 0:
        for i in flip_idx:
            global_ac[i] = tot_ac_cnt[i,0]
            global_af[i] = 1. - global_af[i]
    print("HERE!")
    # Generate ALT AF per subpopulation 
    allel_cnt_subpops = gt.count_alleles_subpops(pop_idx_dict)
    flip_allele = ~flip_af
    flip_allele = flip_allele.astype(np.int8) 
    f = lambda acnt, flip: acnt[:,flip]/acnt.sum(axis=1)
    # Setting up the final data frame
    af_dict = {i: f(allel_cnt_subpops[i], flip_allele) for i in allel_cnt_subpops}
    af_dict2 = {i: af_dict[i].tolist() for i in af_dict}
    # af_df = pd.DataFrame(af_dict)
    # # Inserting all of the columns that are needed
    # af_df = af_df.insert(0, 'CHR', chrom)
    # af_df = af_df.insert(1, 'SNP', pos)
    # af_df = af_df.insert(2, 'A1', ref_alleles)
    # af_df = af_df.insert(3, 'A2', alt_alleles)
    # af_df = af_df.insert(4, 'MAC', global_ac)
    # af_df = af_df.insert(5, 'MAF', global_af)
    return(af_dict2)

