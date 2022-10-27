"""
tests the import of docna
"""

import numpy as np

from docna import Run

def test_safemax():
    assert Run.safemax([0, 1]) == 1

def test_safemaxfail():
    assert Run.safemax([np.nan]) == 0
