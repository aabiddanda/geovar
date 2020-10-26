"""Initialization for GeoVar."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .binning import GeoVar
from .viz import GeoVarPlot
from .utils import vcf_to_freq_table, sep_freq_mat_pops

__all__ = ["GeoVar", "GeoVarPlot", "vcf_to_freq_table", "sep_freq_mat_pops"]

__version__ = "0.1.0"
