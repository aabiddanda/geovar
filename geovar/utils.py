"""Utilities for file conversions for GeoVar."""

import warnings
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
from cyvcf2 import VCF


def flip_alleles(afreq, flip=False):
    """Flip alleles if based on being the minor allele."""
    assert np.all(afreq >= 0.0)
    assert np.all(afreq <= 1.0)
    assert afreq.ndim == 1
    flipped = 1.0 - afreq
    return flipped


def read_pop_panel(pop_panel_file):
    """Read in a population panel file."""
    pop_panel_file_path = Path(pop_panel_file)
    if not pop_panel_file_path.is_file():
        raise ValueError(f"{pop_panel_file} is not a file!")
    else:
        if pop_panel_file_path.suffix == ".csv":
            pop_df = pd.read_csv(pop_panel_file_path, usecols=["sample", "pop"])
        elif pop_panel_file_path.suffix == ".tsv":
            pop_df = pd.read_csv(
                pop_panel_file_path, sep="\t", usecols=["sample", "pop"]
            )
        else:
            pop_df = pd.read_csv(
                pop_panel_file_path, sep=r"\s+", usecols=["sample", "pop"]
            )
        return pop_df


def verify_sample_indices(pop_df, samples):
    """Generate the sample indices."""
    if type(samples) == list:
        samples = np.asarray(samples, dtype=str)
    pop_dict = pop_df.set_index(["sample"]).to_dict()["pop"]
    unique_pops = np.unique([pop_dict[i] for i in pop_dict])
    pop_vector = np.repeat("", samples.size)
    for i, s in enumerate(samples):
        try:
            pop_vector[i] = pop_dict[s]
        except KeyError:
            warnings.warn(f"Sample {s} does not have a population label!", UserWarning)
    pop_vector = np.array(pop_vector)
    pop_idx_dict = {}
    for p in unique_pops:
        pop_idx_dict[p] = np.where(pop_vector == p)[0]
    return unique_pops, pop_idx_dict, pop_dict


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


def vcf_to_freq_table(vcf_file, pop_df, outfile=None, minor_allele=True, **kwargs):
    """Convert a VCF File to a frequency table to be used as input to a GeoVar object.

    Args:
        vcf_file (:obj:`string`): filepath to VCF file (can be bgzipped).
        pop_df (:obj:`pandas.DataFrame`): population data frame in pandas format.
        outfile (:obj:`string`): file to write output allele frequency table to.
        minor_allele (:obj:`bool`): flag to indicate if we want to polarize to the minor allele.

    """
    vcf_filepath = Path(vcf_file)
    if not vcf_filepath.is_file():
        raise ValueError(f"{vcf_file} is not a valid VCF file!")
    vcf = VCF(vcf_filepath, **kwargs)
    unique_pops, pop_idx_dict, pop_dict = verify_sample_indices(pop_df, vcf.samples)
    chrom = []
    pos = []
    ref_alleles = []
    alt_alleles = []
    global_af = []
    global_ac = []
    alt_freq = []
    allele_cnt_subpops = []
    for variant in tqdm(vcf):
        chrom.append(variant.CHROM)
        pos.append(variant.POS)
        ref_alleles.append(variant.REF)
        alt_alleles.append(variant.ALT)
        global_ac.append(variant.AC)
        alt_freq.append(variant.aaf)
        cur_gt = variant.gt_types.copy()
        pop_ac_cnt = [
            cur_gt[pop_idx_dict[i].sum() / pop_idx_dict[i].size] for i in unique_pops
        ]
        allele_cnt_subpops.append(pop_ac_cnt)

    allele_cnt_subpops = np.array(allele_cnt_subpops)
    global_af = alt_freq
    flip_af = np.repeat(False, alt_freq.size)
    if minor_allele:
        flip_af = alt_freq > 0.5
    flip_idx = np.where(flip_af)[0]
    if flip_idx.size > 0:
        for i in flip_idx:
            global_ac[i] = 2 * vcf.samples.size - global_ac[i]
            global_af[i] = 1.0 - global_af[i]
            # Swap the alleles here ...
            cur_ref = ref_alleles[i]
            ref_alleles[i] = alt_alleles[i]
            alt_alleles[i] = cur_ref
            allele_cnt_subpops[i, :] = 1.0 - allele_cnt_subpops
    # Setting up the final data frame
    af_dict = {x: allele_cnt_subpops[:, i] for (i, x) in enumerate(unique_pops)}
    af_df = pd.DataFrame(af_dict)
    # Inserting all of the columns for allele frequencies
    af_df.insert(0, "CHR", chrom)
    af_df.insert(1, "SNP", pos)
    af_df.insert(2, "A1", ref_alleles)
    af_df.insert(3, "A2", alt_alleles)
    af_df.insert(4, "MAC", global_ac)
    af_df.insert(5, "MAF", global_af)
    if outfile is not None:
        af_df.to_csv(outfile, index=False, sep="\t")
    return af_df
