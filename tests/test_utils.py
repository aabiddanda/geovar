"""Testing I/O routines for geovar input."""

import geovar.utils as utils
from hypothesis import given, strategies as st
from hypothesis.extra.numpy import arrays
import numpy as np
import pandas as pd
import pytest


@given(a=arrays(np.int32, 2, elements=st.integers(-10000, 0)))
def test_negative_allele(a):
    """Testing the flipping of allele frequencies."""
    with pytest.raises(AssertionError):
        utils.flip_alleles(a)


@pytest.fixture
def af_test_df():
    """Test great dataframe for allele frequency format."""
    df = pd.DataFrame(
        {
            "CHR": ["chr1", "chr1"],
            "SNP": ["rs1", "rs2"],
            "A1": ["A", "C"],
            "A2": ["G", "T"],
            "MAC": [2, 5],
            "MAF": [0.2, 0.15],
            "POP1": [0.1, 0.2],
            "POP2": [0.05, 0.3],
        }
    )
    return df


@pytest.fixture
def af_test_df_drop(af_test_df):
    """Test dropping some columns."""
    df = af_test_df.drop(["SNP", "MAF"], axis=1)
    return df


def test_sep_freq_mat_pops(af_test_df, expected_pops=["POP1", "POP2"]):
    """Separate frequency matrix of populations."""
    (pops, _) = utils.sep_freq_mat_pops(af_test_df)
    for p in expected_pops:
        assert p in pops


def test_sep_freq_mat_dropped(af_test_df_drop, expected_pops=["POP1", "POP2"]):
    """Separate frequency matrix of populations with a dropped column."""
    (pops, _) = utils.sep_freq_mat_pops(af_test_df_drop)
    for p in expected_pops:
        assert p in pops


@pytest.fixture
def pop_panel_df():
    """Populate population panel dataframe."""
    df = pd.DataFrame({"sample": ["A1", "B1"], "pop": ["A", "A"]})
    return df


@pytest.fixture
def pop_panel_rename(pop_panel_df):
    """Rename a population panel badly."""
    rename_df = pop_panel_df.rename(columns={"pop": "xxx"})
    return rename_df


@pytest.fixture()
def create_pop_panel_csv(tmp_path_factory, pop_panel_df):
    """Create a population panel as a CSV."""
    fn = tmp_path_factory.mktemp("data") / "pop_panel.csv"
    pop_panel_df.to_csv(fn, index=False)
    return fn


@pytest.fixture()
def create_pop_panel_tsv(tmp_path_factory, pop_panel_df):
    """Create a population panel as a TSV."""
    fn = tmp_path_factory.mktemp("data") / "pop_panel.tsv"
    pop_panel_df.to_csv(fn, sep="\t", index=False)
    return fn


@pytest.fixture
def create_pop_panel_txt(tmp_path_factory, pop_panel_df):
    """Create a population panel as a TXT File."""
    fn = tmp_path_factory.mktemp("data") / "pop_panel.txt"
    pop_panel_df.to_csv(fn, sep=" ", index=False)
    return fn


@pytest.fixture()
def create_pop_panel_csv_bad(tmp_path_factory, pop_panel_rename):
    """Create a population panel as a CSV."""
    fn = tmp_path_factory.mktemp("data") / "pop_panel.bad.csv"
    pop_panel_rename.to_csv(fn, index=False)
    return fn


def test_read_pop_panel(
    pop_panel_df, create_pop_panel_csv, create_pop_panel_tsv, create_pop_panel_txt
):
    """Testing the reading of each population panel format."""
    test_files = [create_pop_panel_csv, create_pop_panel_tsv, create_pop_panel_txt]
    for f in test_files:
        df = utils.read_pop_panel(f)
        pd.testing.assert_frame_equal(df, pop_panel_df)


def test_read_bad_pop_panel(create_pop_panel_csv_bad):
    """Testing that population panels must include both samples and pops columns."""
    with pytest.raises(Exception):
        utils.read_pop_panel(create_pop_panel_csv_bad)


@pytest.fixture
def pop_panel_multi_pop_df():
    """Populate population panel dataframe."""
    df = pd.DataFrame({"sample": ["A1", "B1", "C1", "D1"], "pop": ["A", "A", "B", "B"]})
    return df


@pytest.mark.filterwarnings("error")
def test_verify_pops(pop_panel_multi_pop_df):
    """Test that population verification for samples in VCFs works properly."""
    samples_correct = ["A1", "B1", "C1", "D1"]
    samples_subset = ["A1", "C1"]
    samples_excess = ["A1", "E1"]
    utils.verify_sample_indices(pop_panel_multi_pop_df, samples_correct)
    utils.verify_sample_indices(pop_panel_multi_pop_df, samples_subset)
    # This treats warnings as exceptions here...
    with pytest.raises(Exception):
        utils.verify_sample_indices(pop_panel_multi_pop_df, samples_excess)
