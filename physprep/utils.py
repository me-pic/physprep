"""Utilities for Physprep"""

import json
import os
import pickle


def _check_filename(outdir, filename, extension=None, overwrite=False):
    # Check extension if specified
    if extension is not None:
        root, ext = os.path.splitext(filename)
        if not ext or ext != extension:
            filename = root + extension

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

    return filename


def _check_input_validity(option, valid_options, empty=True):
    if type(valid_options) is list:
        if empty:
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
    publication_type = input("Is the source of information a journal ? [y/n] \n")
    if publication_type in ["y", "yes"]:
        ref["journal"] = input("Enter the title of the journal: \n")
        ref["volume"] = input("Enter the volume number: \n")
        ref["issue"] = input("Enter the issue number: \n")
        ref["page"] = input("Enter the page numbers: \n")
        ref["doi"] = input("Enter the DOI: \n")
    else:
        book = input("Is the source of information a book ? [y/n] \n")
        if book in ["y", "yes"]:
            ref["publisher"] = input("Enter the name of the publisher: \n")
            ref["location"] = input(
                "Enter the location of the publisher (city and state/country): \n"
            )

    ref["url"] = input("Enter the URL of the source: \n")

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


def create_config_preprocessing(outdir, filename, overwrite=False):
    """
    outdir: str, Path
    filename: str
    overwrite: bool
    """
    # Instantiate variables
    steps = []
    valid_filters = ["butterworth", "fir", "bessel", "savgol", "notch"]
    valid_steps = [
        "filtering",
        "resampling",
        "resample",
        "upsampling",
        "upsample",
        "downsampling",
        "downsample",
    ]

    filename = _check_filename(outdir, filename, extension=".json", overwrite=overwrite)

    while True:
        tmp = {}
        step = input(
            "\n Enter a processing step among the following: resampling, "
            "filtering.\nIf you do not want to add a step, just press enter.\n"
        )
        step = _check_input_validity(step.lower(), valid_steps, empty=True)
        if step not in ["", " "]:
            method = lowcut = highcut = order = desired_sampling_rate = cutoff = False
            tmp["step"] = step
            if step in ["filtering", "filter"]:
                while method is False:
                    method = input(
                        "\n Enter the filter type among the following: "
                        f"{', '.join(valid_filters)}.\n"
                    )
                    tmp["method"] = _check_input_validity(method.lower(), valid_filters)
                if method == "notch":
                    Q = input("\n Enter the quality factor for the notch filter. \n")
                    tmp["Q"] = _check_input_validity(Q, [int, float])
                if method in ["butterworth", "fir", "bessel"]:
                    while cutoff is False:
                        while lowcut is False:
                            lowcut = input(
                                "\n Enter the lower cutoff frequency (Hz). "
                                "If you do not want to apply a high pass or band "
                                "pass filter, just press enter. \n"
                            )
                            lowcut = _check_input_validity(
                                lowcut, [int, float], empty=True
                            )
                            if lowcut not in ["", " "]:
                                tmp["lowcut"] = lowcut
                        while highcut is False:
                            highcut = input(
                                "\n Enter the higher cutoff frequency (Hz). "
                                "If you do not want to apply a low pass filter "
                                "or band pass filter, just press enter. \n"
                            )
                            highcut = _check_input_validity(
                                highcut, [int, float], empty=True
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
                            "\n Enter the filter order. Must be a positive " "integer.\n"
                        )
                        tmp["order"] = _check_input_validity(order, int)
                if method == "savgol":
                    window_size = input(
                        "\n Enter the length of the filter window. Must be an odd "
                        "integer.\n"
                    )
                    tmp["window_size"] = _check_input_validity(window_size, "odd")

            if step == "signal_resample":
                while desired_sampling_rate is False:
                    desired_sampling_rate = input(
                        "\n Enter the desired sampling frequency "
                        "to resample the signal (in Hz). \n"
                    )
                    desired_sampling_rate = _check_input_validity(
                        desired_sampling_rate, [int, float]
                    )
                    tmp["desired_sampling_rate"] = desired_sampling_rate

            ref = input("\n Is there a reference related to that step ? [y/n] \n")
            if ref in ["y", "yes"]:
                tmp["reference"] = _create_ref()
            steps.append(tmp)
        else:
            print("\n---Saving configuration file---")
            with open(os.path.join(outdir, filename), "w") as f:
                json.dump(steps, f, indent=4)
            break


def create_config_workflow(outdir, filename, overwrite=False):
    """
    outdir: str, Path
    filename: str
    overwrite: bool
    """
    # Instantiate variables
    signals = {}
    valid_signals = ["PPG", "ECG", "EDA", "RSP"]
    preprocessing_strategy = [
        os.path.splitext(f)[0] for f in os.listdir("./data/preprocessing_strategy/")
    ]
    preprocessing_strategy.append("new")

    filename = _check_filename(outdir, filename, extension=".json", overwrite=overwrite)

    while True:
        signal = preprocessing = False
        while signal is False:
            signal = input(
                "\n Enter the type of signal to process. Currently only the (pre-)"
                f"processing of {', '.join(valid_signals)}. \nIf you do not want to add "
                "another type of signal, just press enter.\n"
            )
            signal = _check_input_validity(signal.upper(), valid_signals, empty=True)

        if signal not in ["", " "]:
            signals[signal] = {}
            # Add preprocessing strategy to the workflow
            while preprocessing is False:
                preprocessing = input(
                    "\n Enter the name of the preprocessing "
                    f"strategy to clean the {signal} signal. Choose among the current "
                    "configuration files by providing the name of the strategy, or "
                    "create a new configuration file. To create a new configuration file "
                    "type `new`. Otherwise, choose among those strategy: "
                    f"{', '.join(preprocessing_strategy[:-1])} \n"
                )
                preprocessing = _check_input_validity(
                    preprocessing, preprocessing_strategy, empty=True
                )

            if preprocessing == "new":
                filename_preprocessing = input(
                    "\n Enter the name of the preprocessing "
                    "strategy. The given name will be used as the name of the json file."
                )
                filename_preprocessing = _check_filename(
                    outdir, filename_preprocessing, extension=".json", overwrite=overwrite
                )
                # Create the preprocessing configuration file
                create_config_preprocessing(
                    outdir, filename_preprocessing, overwrite=overwrite
                )
                # Add preprocessing config file directory to the workflow config file
                signals[signal] = {
                    "preprocessing_strategy": os.path.join(outdir, filename_preprocessing)
                }
            else:
                filename_preprocessing = _check_filename(
                    outdir, preprocessing, extension=".json", overwrite=overwrite
                )
                signals[signal] = {
                    "preprocessing_strategy": os.path.join(outdir, filename_preprocessing)
                }

        else:
            print("\n---Saving configuration file---")
            with open(os.path.join(outdir, filename), "w") as f:
                json.dump(signals, f, indent=4)
            break
