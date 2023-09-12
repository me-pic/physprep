"""Test time_sqi functions"""

import operator

import numpy as np
import pytest

from physprep.quality import time_sqi


def test_metrics_hr_sqi():
    """Test metrics_hr_sqi function"""
    intervals = [1000, 1000, 1000, 1000, 1000]
    assert time_sqi.metrics_hr_sqi(intervals, metric="mean") == 60.0
    assert time_sqi.metrics_hr_sqi(intervals, metric="median") == 60.0
    assert time_sqi.metrics_hr_sqi(intervals, metric="std") == 0.0
    assert time_sqi.metrics_hr_sqi(intervals, metric="min") == 60.0
    assert time_sqi.metrics_hr_sqi(intervals, metric="max") == 60.0
    # Test invalid metric
    with pytest.raises(ValueError):
        assert time_sqi.metrics_hr_sqi(intervals, metric="invalid") is np.nan


def test_minimal_range_sqi():
    """Test minimal_range_sqi function"""
    signal = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    assert time_sqi.minimal_range_sqi(signal, 5) == 0.4444


def test_threshold_sqi():
    """Test threshold_sqi function"""
    assert time_sqi.threshold_sqi(6, 4, operator.gt) == "Acceptable"
    assert time_sqi.threshold_sqi(6, 4, operator.lt) == "Not acceptable"
    assert time_sqi.threshold_sqi(6, [2, 8]) == "Acceptable"


def rac_sqi(signal, threshold, duration=2):
    """
    Compute the Rate of Amplitude Change (RAC) in the signal for windows
    of length defines by `duration`.

    Parameters
    ----------
    signal : vector
        Signal on which to compute the RAC.
    threshold : float
        Threshold to consider to evaluate the RAC. If the RAC is above
        that threshold, the quality of the signal in the given window is
        considered as bad.
    duration : float
        Duration of the windows on which to compute the RAC.
        Default to 2.

    Returns
    -------
    rac_ratio : float
        Ratio between the number of windows above `threshold` and the overall
        number of windows.

    References
    ----------
    Böttcher, S., Vieluf, S., Bruno, E., Joseph, B., Epitashvili, N., Biondi, A., ...
        & Loddenkemper, T. (2022). Data quality evaluation in wearable monitoring.
        Scientific reports, 12(1), 21412.
    """
    nb_windows = len(signal) // duration
    rac_values = []

    for i in range(nb_windows):
        window_start = i * duration
        window_end = window_start + duration

        window = signal[window_start:window_end]
        highest_value = max(window)
        lowest_value = min(window)

        rac = abs(highest_value - lowest_value) / min(highest_value, lowest_value)
        rac_values.append(rac)

    rac_ratio = np.count_nonzero(np.array(rac_values) > threshold) / len(signal)

    return np.round(rac_ratio, 4)


def test_rac_sqi():
    """Test rac_sqi function"""
    # Test with a signal with a RAC ratio above the threshold
    assert rac_sqi([1, 2, 3, 4, 5, 6, 7, 8, 9], 0.5, 2) == 0.1111
    # Test with a signal with a RAC ratio below the threshold
    assert rac_sqi([1, 2, 3, 4, 5, 6, 7, 8, 9], 0.6, 1) == 0.0
