"""Test the utils functions"""

import json
import os

import pytest

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
