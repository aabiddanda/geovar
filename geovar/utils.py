
import allel
import pandas as import pd

def vcf_to_freq_table(vcf_file, pop_panel):
    """
        Converts a VCF File to a frequency table 
    """
    vcf_data = allel.read_vcf(vcf_file)
    pass


