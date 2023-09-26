"""Test clean functions"""

import pytest

from physprep.prepare import get_info


def test_order_channels():
    """Test order_channels function"""
    acq_channels = ["PPG_ch", "EDA_ch", "ECG_ch", "TTL"]
    channels = ["cardiac_ppg", "EDA_ch", "cardiac_ecg", "trigger"]
    ch_idx = [1, 3, 4]
    metadata_physio = {
        "cardiac_ecg": {"id": "ECG", "Channel": "ECG_ch"},
        "cardiac_ppg": {"id": "PPG", "Channel": "PPG_ch"},
        "trigger": {"id": "TTL", "Channel": "TTL"},
    }

    ch_names, chsel = get_info.order_channels(acq_channels, metadata_physio)
    assert ch_names == channels
    assert chsel == ch_idx
    with pytest.raises(ValueError):
        ch_names, chsel = get_info.order_channels(acq_channels, {})
    with pytest.raises(ValueError):
        ch_names, chsel = get_info.order_channels(acq_channels[1:], metadata_physio)
