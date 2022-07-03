"""Initialization for GeoVar."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .utils import vcf_to_freq_table, sep_freq_mat_pops, read_pop_panel
from .binning import GeoVar
from .viz import GeoVarPlot

__all__ = [
    "GeoVar",
    "GeoVarPlot",
    "vcf_to_freq_table",
    "sep_freq_mat_pops",
    "read_pop_panel",
]

__version__ = "0.1.0"
