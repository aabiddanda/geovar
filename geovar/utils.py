"""Utilities for file conversions for GeoVar."""

import allel
import pandas as pd
import numpy as np


def flip_alleles(acnt, flip=False):
    """Flip alleles if based on being the minor allele."""
    assert np.all(acnt >= 0)
    assert acnt.ndim == 2
    flipped = acnt[np.arange(len(acnt)), flip] / acnt.sum(axis=1)
    return flipped


def vcf_to_freq_table(vcf_file, pop_panel, outfile=None, minor_allele=True):
    """Convert a VCF File to a frequency table to be used as input to a GeoVar object.

    Args:
        vcf_file (:obj:`string`): filepath to VCF file (can be bgzipped).
        pop_panel (:obj:`string`): filepath to population panel file.
        outfile (:obj:`string`): file to write output allele frequency table to.
        minor_allele (:obj:`bool`): flag to indicate if we want to polarize to the minor allele.

    """
    # NOTE: we only assume that we have two columns in the file
    pop_df = pd.read_table(
        pop_panel, sep=r"\s", usecols=["sample", "pop"], engine="python"
    )
    pop_dict = pop_df.set_index(["sample"]).to_dict()["pop"]
    # NOTE: we are assuming that we have only biallelic markers
    vcf_data = allel.read_vcf(vcf_file, alt_number=1)
    samples = vcf_data["samples"]
    chrom = vcf_data["variants/CHROM"]
    pos = vcf_data["variants/POS"]
    ref_alleles = vcf_data["variants/REF"]
    alt_alleles = vcf_data["variants/ALT"]
    # Generating the population dictionary
    unique_pops = np.unique([pop_dict[i] for i in pop_dict])
    pop_vector = np.array([pop_dict[i] for i in samples])
    pop_idx_dict = {}
    for p in unique_pops:
        pop_idx_dict[p] = np.where(pop_vector == p)[0]
    # Calculating the Alternative Allele Frequency
    gt = allel.GenotypeArray(vcf_data["calldata/GT"])
    # Calculating the total allele frequency
    tot_ac_cnt = gt.count_alleles()
    alt_freq = tot_ac_cnt[:, 1] / tot_ac_cnt.sum(axis=1)
    global_ac = tot_ac_cnt[:, 1]
    global_af = alt_freq
    flip_af = np.repeat(True, alt_freq.size)
    if minor_allele:
        flip_af = alt_freq > 0.5
    # Generate the allele freqquency and allele count vectors
    flip_idx = np.where(flip_af)[0]
    if flip_idx.size > 0:
        for i in flip_idx:
            global_ac[i] = tot_ac_cnt[i, 0]
            global_af[i] = 1.0 - global_af[i]
            # Swap the alleles here ...
            cur_ref = ref_alleles[i]
            ref_alleles[i] = alt_alleles[i]
            alt_alleles[i] = cur_ref
    # Generate ALT AF per subpopulation
    allel_cnt_subpops = gt.count_alleles_subpops(pop_idx_dict)
    flip_allele = ~flip_af
    flip_allele = flip_allele.astype(np.int8)

    # Setting up the final data frame
    af_dict = {
        i: flip_alleles(allel_cnt_subpops[i], flip_allele) for i in allel_cnt_subpops
    }
    af_df = pd.DataFrame(af_dict)
    # Inserting all of the columns that are needed for the allele frequency
    af_df.insert(0, "CHR", chrom)
    af_df.insert(1, "SNP", pos)
    af_df.insert(2, "A1", ref_alleles)
    af_df.insert(3, "A2", alt_alleles)
    af_df.insert(4, "MAC", global_ac)
    af_df.insert(5, "MAF", global_af)
    if outfile is not None:
        af_df.to_csv(outfile, index=False, sep=" ")
    return af_df


def sep_freq_mat_pops(af_df, known_cols=["CHR", "SNP", "A1", "A2", "MAC", "MAF"]):
    """Convert an allele frequency data frame to a frequency array.

    Args:
        af_df (:obj:`pandas.DataFrame`): allele frequency data frame.
        known_cols (:obj:`list`): list of columns to exclude from being a population.

    """
    # Get columns that are not known
    colnames = af_df.columns
    idx = ~np.isin(colnames, known_cols)
    # Generate frequency matrix and the population names
    freq_mat = af_df[af_df.columns[idx]].values
    pop_names = af_df.columns[idx].tolist()
    return (pop_names, freq_mat)
