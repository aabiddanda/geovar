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

@pytest.fixture
def bins_alt1():
    """An alternative binning structure."""
    return [(0, 0), (0, 0.01), (0.01, 0.05), (0.05, 1.0)]


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


def test_generate_bins(bins_less_than_zero, bins_greater_than_one):
    """Test the new generation of bins."""
    geov_obj = GeoVar()
    with pytest.raises(AssertionError):
        geov_obj.generate_bins(bins=bins_less_than_zero)
    with pytest.raises(AssertionError):
        geov_obj.generate_bins(bins=bins_greater_than_one)
    geov_obj.generate_bins(bins=[(0, 0), (0, 0.01), (0.01, 1.0)])


@pytest.fixture
def valid_freq_mat():
    """A Valid frequency matrix for input to geovar."""
    return 'geovar/data/new_1kg_nygc.chr22.biallelic_snps.filt.n5000.freq.txt'

@pytest.fixture
def valid_freq_mat_gz():
    """A Valid frequency matrix for input to geovar as a gzipped file."""
    return 'geovar/data/new_1kg_nygc.chr22.biallelic_snps.filt.n5000.freq.txt.gz'

def test_add_freq_mat(valid_freq_mat):
    """Test for addition of a valid frequency matrix file for GeoVar."""
    geov_obj = GeoVar()
    geov_obj.add_freq_mat(valid_freq_mat)


def test_geovar_binning(valid_freq_mat, bins_alt1):
    """Test GeoVar binning under multiple binning schemes."""
    geov_obj = GeoVar()
    geov_obj.add_freq_mat(valid_freq_mat)
    geov_obj.geovar_binning()
    geov_obj.generate_bins(bins_alt1)
    geov_obj.geovar_binning()

def test_count_geovar_codes(valid_freq_mat, bins_alt1):
    geov_obj = GeoVar()
    geov_obj.add_freq_mat(valid_freq_mat)
    geov_obj.geovar_binning()
    u, n_geovar, ncat = geov_obj.count_geovar_codes()
    # Two non-zero categories ...
    assert ncat == 2
    assert all(n_geovar > 0)



def test_geovar_codes_streaming(valid_freq_mat, valid_freq_mat_gz):
    geov_obj = GeoVar()
    geov_obj.geovar_codes_streaming(valid_freq_mat)
    assert geov_obj.n_variants == 5000
    assert geov_obj.pops.size == geov_obj.n_populations
    geov_obj.geovar_codes_streaming(valid_freq_mat_gz)
    assert geov_obj.n_variants == 5000
    assert geov_obj.pops.size == geov_obj.n_populations
