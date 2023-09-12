"""Test clean functions"""

import numpy as np
from physprep.processing import clean


def test_comb_band_stop():
    # Simulate signal
    time = np.array(np.arange(0, 10, 1))
    raw_data = np.sin(50 * 2 * np.pi * time)
    # Set parameters value
    notches = {"f1": 5, "f2": 6}
    nyquist = 10
    q = 1
    # Apply the notch filter
    filtered_timeseries = clean._comb_band_stop(notches, nyquist, raw_data, q)
    # Check that the filtered signal is the same size as the raw signal
    assert filtered_timeseries.size == raw_data.size