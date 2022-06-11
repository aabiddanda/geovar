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
    """Separate frequency matrix of populations."""
    (pops, _) = utils.sep_freq_mat_pops(af_test_df_drop)
    for p in expected_pops:
        assert p in pops
