"""Test clean functions"""

import pytest

from physprep.prepare import get_info


def test_order_channels():
    """Test order_channels function"""
    acq_channels = ["PPG_ch", "EDA_ch", "ECG_ch", "TTL"]
    channels = ["PPG", "EDA_ch", "ECG", "TTL"]
    ch_idx = [1, 3, 4]
    metadata_physio = {
        "cardiac_ecg": {"id": "ECG", "channel": "ECG_ch"},
        "cardiac_ppg": {"id": "PPG", "channel": "PPG_ch"},
        "trigger": {"id": "TTL", "channel": "TTL"},
    }

    ch_names, chsel = get_info.order_channels(acq_channels, metadata_physio)
    assert ch_names == channels
    assert chsel == ch_idx
    with pytest.raises(ValueError):
        ch_names, chsel = get_info.order_channels(acq_channels, {})
