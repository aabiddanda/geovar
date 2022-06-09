"""Testing I/O routines for geovar input."""

from geovar.binning import GeoVar
import pytest


def test_bin_boundaries():
    """Testing the bin boundaries."""
    bins_less_than_zero = [(-1, 0), (0, 1)]
    bins_greater_than_one = [(0, 0), (0, 2)]
    with pytest.raises(AssertionError):
        GeoVar(bins=bins_less_than_zero)
    with pytest.raises(AssertionError):
        GeoVar(bins=bins_greater_than_one)


def test_empty_bins():
    """Test that the bins are not empty."""
    with pytest.raises(AssertionError):
        GeoVar(bins=[])
