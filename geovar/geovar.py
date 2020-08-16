import numpy as np 
import pandas as pd
from utils import sep_freq_mat_pops

class GeoVar(object):

    def __init__(self, freq_mat_file, bins=[(0, 0), (0, 0.05), (0.05, 1.0)]):
        """
            TODO : DOCUMENT this class and the format of the input file here 
        """
        af_df = pd.read_table(freq_mat_file, sep='\s')
        pops, freq_mat = sep_freq_mat_pops(af_df) 
        self.freq_mat = freq_mat
        self.n_variants = freq_mat.shape[0]
        self.n_populations = freq_mat.shape[1]
        self.pops = pops
        self.bins = bins
        self.geovar_codes = None

    def __str__(self):
        """
            Print the relevant parameters of the objects 
        """
        test_str = 'GeoVar\n' 
        test_str += 'number of variants: %d\n' % self.n_variants
        test_str += 'number of pops: %d\n' % self.n_populations
        test_str += 'pops: ' + ','.join(self.pops) + '\n'
        # NOTE : need to print the bins here
        test_str += 'allele freq bins: ' + ','.join([str(i) for i in self.bins]) 
        return(test_str)
    
    def generate_bins(self, endpts):
        """
           TODO : documentation           
        """
        assert(np.all(np.array(endpts) < 1.0))
        b = 0.0
        y = len(endpts)
        bins = []
        for x in endpts:
            bins.append((b, x))
            b = x
        bins.append((b, 1.0))
        self.bins = bins

    def gen_geodist_binning(self):
        """
            TODO : documentation 
        """
        geovar_codes = np.zeros(shape=self.freq_mat.shape, dtype=np.uint16)
        i = 1
        for b in self.bins[1:]:
            idx = np.where((self.freq_mat > b[0]) & (self.freq_mat <= b[1]))
            geovar_codes[idx] = i
            i += 1
        geovar_codes = np.apply_along_axis(lambda x : ''.join([str(i) for i in x]), 1, geovar_codes)
        self.geovar_codes = geovar_codes

    def count_geovar_codes(self):
        """
           TODO : documentation 
        """
        assert(self.geovar_codes is not None)
        uniq_geodist, n_geodist = np.unique(self.geovar_codes, return_counts=True)
        ncat = np.max(np.vstack([list(x) for x in uniq_geodist]).astype(np.uint32))
        return(uniq_geodist, n_geodist, ncat)