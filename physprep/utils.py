"""Utilities for Physprep"""

import json
import os
import pickle


def _check_filename(outdir, filename, extension=None, overwrite=False):
    # Check if file already exist
    if os.path.exists(os.path.join(outdir, filename).strip()):
        if not overwrite:
            raise IOError(
                "Killing the script because the file already exist. "
                "If you want to overwrite the file, set the `overwrite` "
                "flag to `True`"
            )
        else:
            print(
                "WARNING: This file already exist, but will be overwritten. "
                "If you do not want to overwrite the file, kill the script and "
                "set the `overwrite` flag to `False`"
            )

    # Check extension if specified
    if extension is not None:
        root, ext = os.path.splitext(filename)
        if not ext or ext != extension:
            filename = root + extension

    return filename


def _check_input_validity(option, valid_options):
    if type(valid_options) is list:
        if set(["", " "]) <= set(valid_options):
            if option in ["", " "]:
                return option
        if int in valid_options or float in valid_options:
            if "." in option or "," in option:
                return float(option)
            elif option.isdigit() is False:
                print("**Please enter a positive number.")
                return False
            else:
                return int(option)
        else:
            if option not in valid_options:
                print(f"**Please enter a valid option: {', '.join(valid_options)}.")
                return False
            else:
                if option in [
                    "resampling",
                    "resample",
                    "upsampling",
                    "upsample",
                    "downsampling",
                    "downsample",
                ]:
                    option = "signal_resample"
                return option
    if valid_options in [int, "odd"]:
        if option.isdigit() is False:
            print("**Please enter a positive integer.")
            return False
        if valid_options == "odd":
            if option % 2 != 0:
                print("**Please enter an odd number.")
                return False
            else:
                return int(option)
        else:
            return int(option)
    if valid_options == float:
        if "." in option or "," in option:
            return float(option)
        else:
            print("**Please enter a positive float.")
            return False


def _create_ref():
    ref = {}
    ref["authors"] = input("Enter the author(s) name: \n")
    ref["year"] = input("Enter the publication year: \n")
    ref["title"] = input("Enter the publication title: \n")
    publication_type = input("Is the source of information a journal ? [y/n]")
    if publication_type in ["y", "yes"]:
        ref["journal"] = input("Enter the title of the journal: \n")
        ref["volume"] = input("Enter the volume number: \n")
        ref["issue"] = input("Enter the issue number: \n")
    else:
        book = input("Is the source of information a book ? [y/n]")
        if book in ["y", "yes"]:
            ref["publisher"] = input("Enter the name of the publisher: \n")
            ref["location"] = input(
                "Enter the location of the publisher (city and state/country): \n"
            )
    ref["doi"] = input("Enter the DOI: \n")

    return ref


def load_json(filename):
    """
    Parameters
    ----------
    filename : str
        File path of the .json to load.

    Returns
    -------
    data : dict
        Dictionary with the content of the .json passed in argument.
    """
    try:
        with open(filename, "r") as tmp:
            data = json.loads(tmp.read())
        tmp.close()
    except Exception:
        try:
            with open(filename, "rb") as tmp:
                data = pickle.load(tmp)
            tmp.close()
        except Exception:
            pass

    return data


def create_config(outdir, filename, overwrite=False):
    """
    outdir: str, Path
    filename: str
    """
    # Instantiate variables
    step = None
    steps = []
    tmp = {}
    valid_filters = ["butterworth", "fir", "bessel", "savgol", "notch"]
    valid_steps = [
        "filtering",
        "resampling",
        "resample",
        "upsampling",
        "upsample",
        "downsampling",
        "downsample",
        "",
        " ",
    ]

    filename = _check_filename(outdir, filename, extension=".json", overwrite=overwrite)

    while True:
        method = lowcut = highcut = order = desired_sampling_rate = cutoff = False
        tmp = {}

        step = input(
            "Enter a processing step among the following: resampling, "
            "filtering.\nIf you do not want to add a step, just press enter.\n"
        )
        step = _check_input_validity(step.lower(), valid_steps)
        if step not in ["", " "]:
            tmp["step"] = step
            if step in ["filtering", "filter"]:
                while method is False:
                    method = input(
                        "Enter the filter type among the following: "
                        f"{', '.join(valid_filters)}.\n"
                    )
                    tmp["method"] = _check_input_validity(method.lower(), valid_filters)
                if method == "notch":
                    Q = input("Enter the quality factor for the notch filter. \n")
                    tmp["Q"] = _check_input_validity(Q, [int, float])
                if method in ["butterworth", "fir", "bessel"]:
                    while cutoff is False:
                        while lowcut is False:
                            lowcut = input(
                                "Enter the lower cutoff frequency (Hz). "
                                "If you do not want to apply a high pass or band "
                                "pass filter, just press enter. \n"
                            )
                            lowcut = _check_input_validity(lowcut, [int, float, "", " "])
                            if lowcut not in ["", " "]:
                                tmp["lowcut"] = lowcut
                        while highcut is False:
                            highcut = input(
                                "Enter the higher cutoff frequency (Hz). "
                                "If you do not want to apply a low pass filter "
                                "or band pass filter, just press enter. \n"
                            )
                            highcut = _check_input_validity(
                                highcut, [int, float, "", " "]
                            )
                            if highcut not in ["", " "]:
                                tmp["highcut"] = highcut
                        if lowcut in ["", " "] and highcut in ["", " "]:
                            print(
                                "**Please enter either the filter lower cutoff frequency "
                                "and/or the filter higher cutoff frequency"
                            )
                            lowcut = highcut = False
                        else:
                            cutoff = True
                if method in ["savgol", "butterworth", "bessel"]:
                    while order is False:
                        order = input(
                            "Enter the filter order. Must be a positive " "integer.\n"
                        )
                        tmp["order"] = _check_input_validity(order, int)
                if method == "savgol":
                    window_size = input(
                        "Enter the length of the filter window. Must be an odd integer.\n"
                    )
                    tmp["window_size"] = _check_input_validity(window_size, "odd")

            if step in [
                "resampling",
                "resample",
                "upsampling",
                "upsample",
                "downsampling",
                "downsample",
            ]:
                while desired_sampling_rate is False:
                    desired_sampling_rate = input(
                        "Enter the desired sampling frequency "
                        "to resample the signal (in Hz). \n"
                    )
                    desired_sampling_rate = _check_input_validity(
                        desired_sampling_rate, [int, float]
                    )
                    tmp["desired_sampling_rate"] = desired_sampling_rate

            ref = input("Is there a reference related to that step ? [y/n]")
            if ref in ["y", "yes"]:
                tmp["reference"] = _create_ref()
            steps.append(tmp)
        else:
            with open(os.path.join(outdir, filename), "w") as f:
                json.dump(steps, f, indent=4)
            break
