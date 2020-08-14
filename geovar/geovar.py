import numpy as np 
import pandas as pd
from tqdm import tqdm

class GeoVar(object):

    def __init__(self, freq_mat_file, bins=[(0, 0), (0, 0.05), (0.05, 1.0)]):
        """
            TODO : DOCUMENT this class and the input file here 
        """
        
        self.freq_mat = freq_mat
        self.n_variants = freq_mat.shape[0]
        self.n_populations = freq_mat.shape[1]
        self.bins = bins

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
        output = np.zeros(shape=self.freq_mat.shape, dtype=np.uint16)
        i = 1
        for b in self.bins[1:]:
            idx = np.where((X > b[0]) & (X <= b[1]))
            output[idx] = i
            i += 1
        return(output)
