"""Test time_sqi functions"""

import operator

import numpy as np
import pytest

from physprep.quality import time_sqi

# ==================================================================================
# Signal Quality Indices
# ==================================================================================


def test_sqi_cardiac_overview():
    simulated_data = {
        "Ectopic": 0,
        "Missed": 0,
        "Extra": 0,
        "Long": 0,
        "Short": 0,
        "Peaks": np.arange(0, 10),
    }
    assert time_sqi.sqi_cardiac_overview(simulated_data) == {
        "Ectopic": 0,
        "Missed": 0,
        "Extra": 0,
        "Long": 0,
        "Short": 0,
    }
    simulated_data = {"Ectopic": 0, "Long": 0, "Short": 0, "Peaks": np.arange(0, 10)}
    assert time_sqi.sqi_cardiac_overview(simulated_data) == {
        "Ectopic": 0,
        "Long": 0,
        "Short": 0,
    }
    simulated_data = {"Peaks": np.arange(0, 10)}
    assert time_sqi.sqi_cardiac_overview(simulated_data) == {}


def test_sqi_eda_overview():
    simulated_data = {"a": np.array([]), "b": np.arange(0, 10)}
    assert time_sqi.sqi_eda_overview(len(simulated_data["a"]), 0) == {
        "Quality": "Not acceptable"
    }
    assert time_sqi.sqi_eda_overview(len(simulated_data["b"]), 0) == {
        "Quality": "Acceptable"
    }


# ==================================================================================
# Quality indices : Internal metrics
# ==================================================================================


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


def test_rac_sqi():
    """Test rac_sqi function"""
    # Test with a signal with a RAC ratio above the threshold
    assert time_sqi.rac_sqi([1, 2, 3, 4, 5, 6, 7, 8, 9], 0.5, 2) == 0.1111
    # Test with a signal with a RAC ratio below the threshold
    assert time_sqi.rac_sqi([1, 2, 3, 4, 5, 6, 7, 8, 9], 0.6, 1) == 0.0


def test_threshold_sqi():
    """Test threshold_sqi function"""
    assert time_sqi.threshold_sqi(6, 4, operator.gt) == "Acceptable"
    assert time_sqi.threshold_sqi(6, 4, operator.lt) == "Not acceptable"
    assert time_sqi.threshold_sqi(6, [2, 8]) == "Acceptable"
