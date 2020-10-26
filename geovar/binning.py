"""Module to bin allele frequencies to generate GeoVar codes."""

import numpy as np
import pandas as pd
from tqdm import tqdm
from .utils import sep_freq_mat_pops


class GeoVar(object):
    """GeoVar object that performs binning of allele frequencies."""

    def __init__(self, bins=[(0, 0), (0, 0.05), (0.05, 1.0)]):
        """Object to perform binning of allele frequencies.

        Args:
            bins (:obj:`list`): list of tuples
            containing allele frequency categories.

        """
        assert np.all(np.array(bins) <= 1.0)
        assert np.all(np.array(bins) >= 0.0)
        self.bins = bins
        self.freq_mat = None
        self.n_variants = None
        self.n_populations = None
        self.pops = None
        self.geovar_codes = None

    def __str__(self):
        """Representation of fields in a GeoVar object."""
        test_str = "GeoVar\n"
        test_str += "number of variants: %d\n" % self.n_variants
        test_str += "number of pops: %d\n" % self.n_populations
        test_str += "pops: " + ",".join(self.pops) + "\n"
        test_str += "allele freq bins: " + ",".join([str(i) for i in self.bins])
        return test_str

    def add_freq_mat(self, freq_mat_file):
        """Adding an allele frequency table (see example notebook for format).

        Args:
            freq_mat_file (:obj:`string`): filepath to frequency table file
            (see example notebook for formatting).

        """
        af_df = pd.read_table(freq_mat_file, sep=r"\s")
        pops, freq_mat = sep_freq_mat_pops(af_df)
        self.pops = pops
        self.freq_mat = freq_mat
        self.n_variants = freq_mat.shape[0]
        self.n_populations = freq_mat.shape[1]

    def generate_bins(self, bins=[(0, 0), (0, 0.05), (0.05, 1.0)]):
        """Define new bins for each allele frequency categorization.

        Args:
            bins (:obj:`list`): list of tuples specifying bins of allele frequency.

        """
        assert np.all(np.array(bins) < 1.0)
        b = 0.0
        new_bins = []
        for x in bins:
            new_bins.append((b, x))
            b = x
        new_bins.append((b, 1.0))
        self.bins = new_bins

    def geovar_binning(self):
        """Compute the GeoVar codes for each variant across each population."""
        geovar_codes = np.zeros(shape=self.freq_mat.shape, dtype=np.uint16)
        i = 1
        for b in self.bins[1:]:
            idx = np.where((self.freq_mat > b[0]) & (self.freq_mat <= b[1]))
            geovar_codes[idx] = i
            i += 1
        geovar_codes = np.apply_along_axis(
            lambda x: "".join([str(i) for i in x]), 1, geovar_codes
        )
        self.geovar_codes = geovar_codes

    def count_geovar_codes(self):
        """Count the number of unique GeoVar codes in a dataset."""
        assert self.geovar_codes is not None
        uniq_geovar, n_geovar = np.unique(self.geovar_codes, return_counts=True)
        ncat = np.max(np.vstack([list(x) for x in uniq_geovar]).astype(np.uint32))
        return (uniq_geovar, n_geovar, ncat)

    def geovar_codes_streaming(self, freq_mat_file):
        """Version of GeoVar code generation algorithm that streams through file to avoid memory overflow.

        Args:
            freq_mat_file (:obj:`string`): filepath to
            frequency table file (see example notebook for formatting).

        """
        assert self.bins is not None
        geovar_codes = []
        # Setting up the testing bins
        test_bins = np.array([x[1] for x in self.bins])
        with open(freq_mat_file, "r") as f:
            header = f.readline()
            # Take the population labels currently
            pops = np.array(header.split()[6:])
            self.pops = pops
            for line in tqdm(f):
                # Split after the 6th column ...
                maf_vector = np.array(line.split()[6:]).astype(np.float64)
                cur_geovar = np.digitize(maf_vector, test_bins, right=True)
                cur_geovar_code = "".join([str(i) for i in cur_geovar])
                geovar_codes.append(cur_geovar_code)
        # Setting the variables here
        self.geovar_codes = np.array(geovar_codes)
        self.n_variants = self.geovar_codes.size
        self.n_populations = self.pops.size
