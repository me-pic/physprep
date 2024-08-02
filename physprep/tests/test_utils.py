"""Test the utils functions"""

import json
import os

import pandas as pd
import pytest
from bids import BIDSLayout

from physprep import utils


def test_check_filename():
    """
    Test the check_filename function.
    """
    # Instantiate variables
    outdir = "test"
    filename = "test.json"
    extension = ".json"
    # Create directory
    os.mkdir(outdir)
    # Check filename when file does not exist
    filename = utils._check_filename(outdir, filename, extension)
    # Check
    assert filename == "test.json"

    # Create file
    data = {"test": "test"}
    # Save data
    with open(os.path.join(outdir, filename), "w") as tmp:
        json.dump(data, tmp)
    tmp.close()
    # Check filename when file exists, overwrite = False
    with pytest.raises(FileExistsError):
        filename = utils._check_filename(outdir, filename, extension, overwrite=False)
    # Check filename when file exists, overwrite = True
    filename = utils._check_filename(outdir, filename, extension, overwrite=True)
    assert filename == "test.json"

    # Delete file
    os.remove(os.path.join(outdir, filename).strip())
    # Delete directory
    os.rmdir(outdir)


def test_load_json():
    """
    Test the load_json function.
    """
    # Instantiate variables
    filename = "test.json"
    data = {"test": "test"}
    # Save data
    with open(filename, "w") as tmp:
        json.dump(data, tmp)
    tmp.close()
    # Load data
    data_loaded = utils.load_json(filename)
    # Check
    assert data_loaded == data
    # Delete file
    os.remove(filename)


def test_check_input_validity():  # option, valid_options, empty=True):
    # Test empty
    assert utils._check_input_validity("", ["test"], empty=True) == ""  # valid
    assert utils._check_input_validity(" ", ["test"], empty=True) == " "  # valid
    assert utils._check_input_validity("", ["test"], empty=False) is False  # invalid
    # Test string
    assert utils._check_input_validity("test", ["test"]) == "test"  # valid
    # Test int
    assert utils._check_input_validity("1", int) == 1  # valid
    assert utils._check_input_validity("1.0", int) is False  # invalid
    assert utils._check_input_validity("a", int) is False  # invalid
    assert utils._check_input_validity("1", "odd") == 1  # valid
    assert utils._check_input_validity("2", "odd") is False  # invalid
    # Test float
    assert utils._check_input_validity("1.0", float) == 1.0  # valid
    assert utils._check_input_validity("a", float) is False  # invalid


def test_rename_in_bids():
    snake_case = ["resp_signal", "eda_signal"]
    camel_case = ["RespSignal", "EdaSignal"]
    # Test renaming of DataFrame
    df_renamed = utils.rename_in_bids(
        pd.DataFrame({"RESP_Signal": ["test"], "EDA_Signal": ["test"]})
    )
    assert df_renamed.columns.tolist() == snake_case
    df_renamed = utils.rename_in_bids(
        pd.DataFrame({"RespSignal": ["test"], "EdaSignal": ["test"]})
    )
    assert df_renamed.columns.tolist() == snake_case
    df_renamed = utils.rename_in_bids(
        pd.DataFrame({"respSignal": ["test"], "edaSignal": ["test"]})
    )
    assert df_renamed.columns.tolist() == snake_case
    df_renamed = utils.rename_in_bids(
        pd.DataFrame({"resp_signal": ["test"], "eda_signal": ["test"]})
    )
    assert df_renamed.columns.tolist() == snake_case

    # Test renaming of dict
    dict_renamed = utils.rename_in_bids({"RESP_Signal": ["test"], "EDA_Signal": ["test"]})
    assert list(dict_renamed.keys()) == camel_case
    dict_renamed = utils.rename_in_bids({"RespSignal": ["test"], "EdaSignal": ["test"]})
    assert list(dict_renamed.keys()) == camel_case
    dict_renamed = utils.rename_in_bids({"respSignal": ["test"], "edaSignal": ["test"]})
    assert list(dict_renamed.keys()) == camel_case

def test_check_bids_validity():
    # Create temporary fake bids dataset
    tmp_bids_dataset = '../data/bids_test'
    os.mkdir(tmp_bids_dataset)
    with open(os.path.join(tmp_bids_dataset, "dataset_description.json"), "w") as f:
        json.dump({}, f, indent=4)
    layout_bids_dataset = utils._check_bids_validity(tmp_bids_dataset)
    assert type(layout_bids_dataset) == BIDSLayout
    # Create temporary fake non-bids dataset
    tmp_dataset = '../data/test'
    os.mkdir(tmp_dataset)
    layout_dataset = utils._check_bids_validity(tmp_dataset)
    assert os.path.isfile(os.path.join(tmp_dataset, 'dataset_description.json'))
    assert type(layout_dataset) == BIDSLayout
    
    # Delete temporary dataset
    os.remove(os.path.join(tmp_bids_dataset, "dataset_description.json").strip())
    os.remove(os.path.join(tmp_dataset, "dataset_description.json").strip())
    os.rmdir(tmp_bids_dataset)
    os.rmdir(tmp_dataset)