"""Testing I/O routines for geovar input."""

from geovar.binning import GeoVar
import pytest


@pytest.fixture
def bins_less_than_zero():
    """Bins that are below zero."""
    return [(-1, 0), (0, 1)]


@pytest.fixture
def bins_greater_than_one():
    """Bins that are greater than one."""
    return [(0, 0.5), (0.5, 1.1)]


def test_bin_boundaries(bins_less_than_zero, bins_greater_than_one):
    """Testing the bin boundaries."""
    with pytest.raises(AssertionError):
        GeoVar(bins=bins_less_than_zero)
    with pytest.raises(AssertionError):
        GeoVar(bins=bins_greater_than_one)


def test_empty_bins():
    """Test that the bins are not empty."""
    with pytest.raises(AssertionError):
        GeoVar(bins=[])
