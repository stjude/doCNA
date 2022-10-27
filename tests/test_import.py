"""
tests the import of docna
"""

import numpy as np

from doCNA import Run

def test_safemax():
    assert Run.safemax([0, 1]) == 1

def test_safemaxfail():
    assert Run.safemax([]) == 0
