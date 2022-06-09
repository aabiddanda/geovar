"""Testing I/O routines for geovar input."""

from geovar.utils import flip_alleles
from hypothesis import given, strategies as st
from hypothesis.extra.numpy import arrays
import numpy as np
import pytest


@given(a=arrays(np.int32, 2, elements=st.integers(-10000, 0)))
def test_negative_allele(a):
    """Testing the flipping of allele frequencies."""
    with pytest.raises(AssertionError):
        flip_alleles(a)
