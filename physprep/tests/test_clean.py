"""Test clean functions"""

import numpy as np
import pandas as pd

from physprep.processing import clean


def test_preprocessing_workflow():
    pass


def test_remove_padding():
    # Simulate signal
    simulated_data = {
        "time": np.arange(-5, 10, 1),
        "TTL": np.concatenate((np.zeros(5), np.full(5, 4.1), np.zeros(5)), axis=None),
    }
    simulated_data = pd.DataFrame(simulated_data)
    # Remove padding based on the trigger threshold
    cropped_data = clean.remove_padding(simulated_data, trigger_threshold=4)
    assert len(cropped_data) == 5
    # Remove padding based on `start_time`
    cropped_data = clean.remove_padding(simulated_data, start_time=0)
    assert len(cropped_data) == 10
    # Remove padding based on `end_time`
    cropped_data = clean.remove_padding(simulated_data, end_time=5)
    assert len(cropped_data) == 11
    # Remove padding based on `start_time` and `end_time`
    cropped_data = clean.remove_padding(simulated_data, start_time=0, end_time=5)
    assert len(cropped_data) == 6
    # Remove padding when start_time = end_time = trigger_threshold = None
    cropped_data = clean.remove_padding(simulated_data)
    assert len(cropped_data) == len(simulated_data)
    # Remove padding when start_time, end_time and trigger_threshold are not None
    cropped_data = clean.remove_padding(
        simulated_data, start_time=0, end_time=5, trigger_threshold=4
    )
    assert len(cropped_data) == 5


def test_preprocess_signal():
    pass


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
