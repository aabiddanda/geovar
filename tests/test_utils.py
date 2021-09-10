"""Testing I/O routines for geovar input."""

import sys
sys.path.append('geovar/')
from utils import *
from hypothesis import given, strategies as st
from hypothesis.extra.numpy import arrays
import pytest

@given(a=arrays(np.int32, 2, elements=st.integers(-10000, 0)))
def test_allele(a):
    """Testing the flipping of allele frequencies."""
    with pytest.raises(AssertionError):
        flip_alleles(a)
