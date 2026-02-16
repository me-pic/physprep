"""Utilities for Physprep"""

import json
import os
import pickle
import pkgutil
import re
import warnings
from pathlib import Path

import pandas as pd
from bids import BIDSLayout
from bids.exceptions import BIDSValidationError
from bids.layout import BIDSFile
from pkg_resources import resource_filename

WORKFLOW_STRATEGIES = ["neuromod"]
PREPROCESSING_STRATEGIES = [
    "neuromod_ecg",
    "neuromod_eda",
    "neuromod_ppg",
    "neuromod_rsp",
]
QA_STRATEGIES = ["neuromod_cardiac", "neuromod_eda", "neuromod_rsp"]


def _check_filename(outdir, filename, extension=None, overwrite=False):
    # Check extension if specified
    if extension is not None:
        root, ext = os.path.splitext(filename)
        if not ext or ext != extension:
            filename = root + extension

    # Check if file already exist
    if os.path.exists(os.path.join(outdir, filename).strip()):
        if not overwrite:
            raise FileExistsError(
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
                if option == "resampling":
                    option = "signal_resample"
                return option
    if valid_options in [int, "odd"]:
        if option.isdigit() is False:
            print("**Please enter a positive integer.")
            return False
        if valid_options == "odd":
            if int(option) % 2 == 0:
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
    # Instantiate variable
    ref = {}
    publication_type = book = False
    # Collect input
    ref["authors"] = input("Enter the author(s) name: \n")
    ref["year"] = input("Enter the publication year: \n")
    ref["title"] = input("Enter the publication title: \n")
    while publication_type is False:
        publication_type = input("Is the source of information a journal ? [y/n] \n")
        publication_type = _check_input_validity(publication_type, ["y", "n"])
    if publication_type == "y":
        ref["journal"] = input("Enter the title of the journal: \n")
        ref["volume"] = input("Enter the volume number: \n")
        ref["issue"] = input("Enter the issue number: \n")
        ref["page"] = input("Enter the page numbers: \n")
        ref["doi"] = input("Enter the DOI: \n")
    else:
        while book is False:
            book = input("Is the source of information a book ? [y/n] \n")
            book = _check_input_validity(book, ["y", "n"])
        if book == "y":
            ref["publisher"] = input("Enter the name of the publisher: \n")
            ref["location"] = input(
                "Enter the location of the publisher (city and state/country): \n"
            )

    ref["url"] = input("Enter the URL of the source: \n")

    return ref


def _check_sub_validity(sub, bids_sub):
    # Make sure sub is a list
    if isinstance(sub, str):
        sub = [sub]
    # Get the common elements between the specified subjects in `sub`,
    # and the ones founds in`bids_sub`
    valid_sub = list(set(sub).intersection(bids_sub))
    invalid_sub = list(set(sub) - set(valid_sub))
    if len(invalid_sub) > 0:
        warnings.warn(
            "The following subject were specified in `sub`, but were not found in "
            f"{bids_sub}: {invalid_sub}. \n Only {valid_sub} will be preprocessed"
        )

    return valid_sub


def _check_ses_validity(ses, bids_ses):
    # Make sure sub is a list
    if isinstance(ses, str):
        ses = [ses]
    # Get the common elements between the specified subjects in `ses`,
    # and the ones founds in`bids_ses`
    valid_ses = list(set(ses).intersection(bids_ses))
    invalid_ses = list(set(ses) - set(valid_ses))
    if len(invalid_ses) > 0:
        warnings.warn(
            "The following sessions were specified in `ses`, but were not found in "
            f"{bids_ses}: {invalid_ses}. \n Only {valid_ses} will be preprocessed"
        )

    return valid_ses


def _check_bids_validity(path, is_derivative=False):
    # Check if path is a BIDS dataset, otherwise create a `dataset_description.json` file
    # Reference: https://github.com/bids-standard/bids-starter-kit/
    if isinstance(path, BIDSLayout):
        return path
    else:
        try:
            layout = BIDSLayout(path, validate=False, is_derivative=is_derivative)
        except BIDSValidationError:
            warnings.warn(
                f"Because {path} is not a BIDS dataset, an empty "
                "`dataset_description.json` file will be created at the root. "
                "MAKE SURE TO FILL THAT FILE AFTERWARD !"
            )
            if is_derivative:
                descrip = pkgutil.get_data(
                    __name__, "data/boilerplates/dataset_description_derivatives.json"
                )
            else:
                descrip = pkgutil.get_data(
                    __name__, "data/boilerplates/dataset_description.json"
                )
            with open(f"{path}/dataset_description.json", "w") as f:
                json.dump(json.loads(descrip.decode()), f, indent=4)
            f.close()
            layout = BIDSLayout(path, validate=False, is_derivative=is_derivative)

        return layout


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


def save_processing(outdir, bids_entities, descriptor, data, metadata, preproc_strategy, save_raw=False):
    """
    outdir: str or pathlib.Path
        Path to the directory where the preprocessed physiological data will be saved.
    bids_entities: str
        BIDS entities to use for saving data
    descriptor; str
        Descriptor that will be used for filename
    data: dict
        Dictionary containing the timeseries for each signal.
    metadata: dict
        Dictionary containing the metadata for each signal
    save_raw: bool
        If True, save the raw timeseries. Otherwise, only save the processed timeseries
    """
    # Save preprocessed signal
    if outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        code = outdir / 'code'
        code.mkdir(parents=True, exist_ok=True)
    else:
        warnings.warn(
            "No output directory specified. Data will be saved in the "
            f"current working directory: {Path.cwd()}\n"
        )
        outdir = Path.cwd()

    # Define BIDS layout for derivatives dataset
    layout_deriv = _check_bids_validity(outdir, is_derivative=True)
    # Add desc entity to dict
    bids_entities["desc"] = descriptor
    # Define pattern to build path for derivatives
    deriv_pattern = "sub-{subject}[/ses-{session}]/{datatype}/sub-{subject}[_ses-{session}][_task-{task}][_run-{run}][_recording-{recording}][_desc-{desc}]_{suffix}.{extension}"

    if not save_raw:
        to_keep = [col for col in data.keys() if "raw" not in col]
    else:
        to_keep = [col for col in data.keys()]

    if isinstance(metadata["SamplingFrequency"], list):
        unique_sf = list(set(metadata["SamplingFrequency"]))
        # Save derivatives with different SamplingFrequency in different files
        for f in unique_sf:
            get_idx = [i for i, s in enumerate(metadata["SamplingFrequency"]) if s == f]
            cols = [
                metadata["Columns"][i]
                for i in get_idx
                if metadata["Columns"][i] in to_keep
            ]

            df = pd.DataFrame({col: data[col] for col in cols})

            bids_entities["recording"] = f"{f}Hz"

            filename = layout_deriv.build_path(
                bids_entities, deriv_pattern, validate=False
            )
            # Make sure directory exists
            Path(BIDSFile(filename).dirname).mkdir(parents=True, exist_ok=True)
            # Save data
            df.to_csv(filename, sep="\t", index=False, compression="gzip")
    else:
        df = pd.DataFrame({col: data[col] for col in to_keep})

        filename = layout_deriv.build_path(bids_entities, deriv_pattern, validate=False)
        # Make sure directory exists
        Path(BIDSFile(filename).dirname).mkdir(parents=True, exist_ok=True)
        # Save data
        df.to_csv(filename, sep="\t", index=False, compression="gzip")

    # Save preproc_strategy
    with open(code / "preprocessing_strategy.json", "w") as json_file:
        json.dump(preproc_strategy, json_file, indent=4)


def save_features(outdir, bids_entities, events):
    # Get BIDS layout
    layout_deriv = _check_bids_validity(outdir, is_derivative=True)
    # Define pattern
    deriv_pattern = "sub-{subject}[/ses-{session}]/{datatype}/sub-{subject}[_ses-{session}][_task-{task}][_run-{run}][_recording-{recording}]_desc-physio_events.tsv"
    # Build path
    filename = layout_deriv.build_path(bids_entities, deriv_pattern, validate=False)
    # Make sure directory exists
    Path(BIDSFile(filename).dirname).mkdir(parents=True, exist_ok=True)
    # Save data
    events.to_csv(filename, sep="\t", index=False)


def save_qa(outdir, bids_entities, qa_output, short=False, report=False):
    # Get BIDS layout
    layout_deriv = _check_bids_validity(outdir, is_derivative=True)
    # Define pattern
    if report:
        bids_entities["desc"] = "qareport"
        bids_entities["suffix"] = "html"
    else:
        if short:
            bids_entities["desc"] = "quality"
        else:
            bids_entities["desc"] = "qualitydesc"
        bids_entities["suffix"] = "json"
    deriv_pattern = "sub-{subject}[/ses-{session}]/{datatype}/sub-{subject}[_ses-{session}][_task-{task}][_run-{run}][_recording-{recording}]_desc-{desc}.{suffix}"
    # Build path
    filename = layout_deriv.build_path(bids_entities, deriv_pattern, validate=False)
    # Make sure directory exists
    Path(BIDSFile(filename).dirname).mkdir(parents=True, exist_ok=True)

    # Save qa
    if report:
        with open(filename, "w") as f:
            f.write(qa_output)
            f.close()
    else:
        with open(filename, "w") as f:
            json.dump(qa_output, f, indent=4)
            f.close()


def create_config_qa(outdir, filename, modality, overwrite=False):
    """
    Generate a configuration file for the qa strategy based on the user inputs.

    Parameters
    ----------
    outdir: str, pathlib.Path
        Saving directory.
    filename: str
        Saving filename.
    modality: str
        Modality to consider. Possible choices: 'PPG', 'ECG', 'EDA', 'RSP'.
    overwrite: bool
        If `True`, overwrite the existing file with the specified `filename` in the
        `outdir` directory. Default is False.
    """
    # Instantiate variables
    tmp = []
    valid_steps = ["sliding", "metric", "feature"]

    # Validate filename
    filename = _check_filename(outdir, filename, extension=".json", overwrite=overwrite)

    # Validate modality:
    if modality not in ["ECG", "PPG", "RSP", "EDA"]:
        raise ValueError(
            f"modality {modality} is not a valid value. Please "
            "among ['ECG', 'PPG', 'RSP', 'EDA']"
        )

    if modality in ["ECG", "PPG"]:
        valid_metrics = [
            "Mean",
            "Median",
            "SD",
            "Min",
            "Max",
            "Skewness",
            "Kurtosis",
            "Quality",
        ]
        valid_features = ["NN_intervals", "HR", "signal"]
    elif modality == "RSP":
        valid_metrics = [
            "Mean",
            "Median",
            "SD",
            "Min",
            "Max",
            "Skewness",
            "Kurtosis",
            "CV",
            "Variability",
            "Quality",
        ]
        valid_features = ["Amplitude", "Rate"]
    elif modality == "EDA":
        valid_metrics = [
            "Mean",
            "Median",
            "SD",
            "Min",
            "Max",
            "Skewness",
            "Kurtosis",
            "Detected_peaks",
            "Quality",
        ]
        valid_features = ["signal", "Tonic", "Phasic", "Peaks"]

    while True:
        step = input(
            "\nEnter a qa step among the following: `sliding` to "
            "specify a window in which to compute some qa metrics. \nIf "
            "you want the metrics to be compute on the whole signal, "
            "do not enter `sliding`. \nIf you want to enter a metric, "
            "enter `metric`. \nIf you do not want to add a step, just "
            "press enter.\n"
        )
        step = _check_input_validity(step.lower(), valid_steps, empty=True)

        if step not in ["", " "]:
            if step == "sliding":
                is_sliding = [True for elem in tmp if "sliding" in elem.keys()]
                if not is_sliding:
                    duration = step_window = False
                    while duration is False:
                        duration = input(
                            "\nEnter the duration of the window (in seconds): \n"
                        )
                        duration = _check_input_validity(
                            duration, [int, float], empty=True
                        )
                    while step_window is False:
                        step_window = input(
                            "\nEnter the step of the window (in seconds), if you want "
                            "to compute the metrics in rolling windows: \n"
                        )
                        if step_window in ["", " "]:
                            step_window = 0
                        step_window = _check_input_validity(
                            step_window, [int, float], empty=True
                        )
                    tmp.append({"sliding": {"duration": duration, "step": step_window}})
                else:
                    print("\nA window was already specified. \n")
            elif step == "metric":
                invalid_metric = True
                while invalid_metric:
                    metric = input(
                        "\nEnter the metric you want to compute among those options: "
                        f"{', '.join(valid_metrics)}.\n"
                    )
                    if metric in valid_metrics:
                        invalid_metric = False
                    else:
                        print("\n WARNING: invalid metric")
                invalid_feature = True
                while invalid_feature:
                    feature = input(
                        f"\nEnter the features on which you want to compute the {metric} "
                        "among the following options: "
                        f"{', '.join(valid_features)}.\n"
                    )
                    if feature in valid_features:
                        invalid_feature = False
                    else:
                        print("\n WARNING: invalid feature")

                tmp.append({"metric": metric, "feature": feature})
        else:
            # Save the configuration file only if there is at least one signal
            if bool(tmp):
                print("\n---Saving configuration file---")
                with open(os.path.join(outdir, filename), "w") as f:
                    json.dump(tmp, f, indent=4)
            break


def create_config_preprocessing(outdir, filename, overwrite=False):
    """
    Generate a configuration file for the preprocessing strategy based on the user inputs.

    Parameters
    ----------
    outdir: str, pathlib.Path
        Saving directory.
    filename: str
        Saving filename.
    overwrite: bool
        If `True`, overwrite the existing file with the specified `filename` in the
        `outdir` directory. Default is False.
    """
    # Instantiate variables
    steps = []
    valid_filters = ["butterworth", "fir", "bessel", "savgol", "notch"]
    valid_noth_method = ["biopac", "bottenhorn"]
    valid_steps = ["filtering", "resampling"]

    filename = _check_filename(outdir, filename, extension=".json", overwrite=overwrite)

    while True:
        tmp = {}
        tmp_params = {}
        step = input(
            "\n Enter a processing step among the following: resampling, "
            "filtering.\nIf you do not want to add a step, just press enter.\n"
        )
        step = _check_input_validity(step.lower(), valid_steps, empty=True)
        if step not in ["", " "]:
            method = lowcut = highcut = order = desired_sampling_rate = cutoff = ref = (
                notch_method
            ) = Q = tr = slices = mb = False
            tmp["step"] = step
            if step in ["filtering", "filter"]:
                while method is False:
                    method = input(
                        "\n Enter the filter type among the following: "
                        f"{', '.join(valid_filters)}.\n"
                    )
                    tmp_params["method"] = _check_input_validity(
                        method.lower(), valid_filters
                    )
                if method == "notch":
                    while Q is False:
                        Q = input("\n Enter the quality factor for the notch filter. \n")
                        Q = _check_input_validity(Q, [int, float])
                        tmp_params["Q"] = Q
                    while notch_method is False:
                        notch_method = input(
                            "\n Enter the notch filter method among the following:  "
                            "biopac, bottenhorn.\n "
                        )
                        notch_method = _check_input_validity(
                            notch_method, valid_noth_method
                        )
                        tmp_params["notch_method"] = notch_method
                    while tr is False:
                        tr = input("\n Enter the tr used to acquired the fMRI data. \n")
                        tr = _check_input_validity(tr, [int, float])
                        tmp_params["tr"] = tr
                    while slices is False:
                        slices = input(
                            "\n Enter the number of slices used to acquired the fMRI "
                            "data.\n "
                        )
                        slices = _check_input_validity(slices, int)
                        tmp_params["slices"] = slices
                    if notch_method == "bottenhorn":
                        while mb is False:
                            mb = input(
                                "\n Enter the multi-band acceleration factor used to "
                                "acquired the fMRI data.\n "
                            )
                            mb = _check_input_validity(mb, int)
                            tmp_params["mb"] = mb
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
                                tmp_params["lowcut"] = lowcut
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
                                tmp_params["highcut"] = highcut
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
                        order = _check_input_validity(order, int)
                        tmp_params["order"] = order
                if method == "savgol":
                    window_size = input(
                        "\n Enter the length of the filter window. Must be an odd "
                        "integer.\n"
                    )
                    window_size = _check_input_validity(window_size, "odd")
                    tmp_params["window_size"] = window_size
            if step == "signal_resample":
                while desired_sampling_rate is False:
                    desired_sampling_rate = input(
                        "\n Enter the desired sampling frequency "
                        "to resample the signal (in Hz). \n"
                    )
                    desired_sampling_rate = _check_input_validity(
                        desired_sampling_rate, [int, float]
                    )
                    tmp_params["desired_sampling_rate"] = desired_sampling_rate

            tmp["parameters"] = tmp_params

            if step in ["filtering", "filter", "signal_resample"]:
                while ref is False:
                    ref = input("\n Is there a reference related to that step ? [y/n] \n")
                    ref = _check_input_validity(ref, ["y", "n"], empty=False)
                if ref == "y":
                    tmp["reference"] = _create_ref()
                steps.append(tmp)
        else:
            print("\n---Saving configuration file---")
            with open(os.path.join(outdir, filename), "w") as f:
                json.dump(steps, f, indent=4)
            break


def create_config_workflow(outdir, filename, overwrite=False):
    """
    Generate a configuration file for the workflow strategy based on the user inputs.

    Parameters
    ----------
    outdir: str, pathlib.Path
        Saving directory.
    filename: str
        Saving filename.
    overwrite: bool
        If `True`, overwrite the existing file with the specified `filename` in the
        `outdir` directory. Default: False.
    """
    # Instantiate variables
    signals = {}
    valid_signals = [
        "cardiac_ppg",
        "PPG",
        "cardiac_ecg",
        "ECG",
        "electrodermal",
        "EDA",
        "respiratory",
        "RSP",
        "RESP",
        "trigger",
        "TTL",
    ]
    preprocessing_strategy = [
        os.path.splitext(f)[0]
        for f in os.listdir("./physprep/data/preprocessing_strategy/")
    ]
    preprocessing_strategy.append("new")

    qa_strategy = [
        os.path.splitext(f)[0] for f in os.listdir("./physprep/data/qa_strategy/")
    ]
    qa_strategy.append("new")

    filename = _check_filename(outdir, filename, extension=".json", overwrite=overwrite)

    while True:
        signal = preprocessing = qa = False
        while signal is False:
            signal = input(
                "\n Enter the type of signal to process. Currently only the (pre-)"
                f"processing of {', '.join(valid_signals)}. \nIf you do not want to add "
                "another type of signal, just press enter.\n"
            )
            signal = _check_input_validity(signal, valid_signals, empty=True)

        if signal not in ["", " "]:
            signals[signal] = {}
            # Associate abbreviation to the signal type
            if signal in ["PPG", "cardiac_ppg"]:
                signals[signal] = {
                    "id": "PPG",
                    "Description": "continuous pulse measurement",
                    "Units": "V",
                }
            elif signal in ["ECG", "cardiac_ecg"]:
                signals[signal] = {
                    "id": "ECG",
                    "Description": "continuous electrocardiogram measurement",
                    "Units": "mV",
                }
            elif signal in ["EDA", "GSR", "electrodermal"]:
                signals[signal] = {
                    "id": "EDA",
                    "Description": "continuous electrodermal measurement",
                    "Units": "microsiemens",
                }
            elif signal in ["RSP", "RESP", "respiratory"]:
                signals[signal] = {
                    "id": "RESP",
                    "Description": "continuous breathing measurement",
                    "Units": "cm H2O",
                }
            elif signal == "trigger":
                signals[signal] = {
                    "id": "TTL",
                    "Description": "continuous measurement of the scanner trigger signal",
                    "Units": "V",
                }

            # Ask for the channel name associated with the signal
            channel = input(
                "\n Enter the name of the channel in your acq file associated with the "
                f"{signal} signal: \n"
            )
            signals[signal].update({"Channel": channel})

            if signal != "trigger":
                # Add preprocessing strategy to the workflow
                while preprocessing is False:
                    preprocessing = input(
                        "\n Enter the name of the preprocessing "
                        f"strategy to clean the {signal} signal. Choose among the "
                        "current configuration files by providing the name of the "
                        "strategy, \nspecify the path to an existing file, "
                        "or create a new configuration file. To create a "
                        "new configuration file type `new`.\n Otherwise, choose among "
                        f"those strategy: {', '.join(preprocessing_strategy[:-1])}.\n"
                    )

                if preprocessing == "new":
                    filename_preprocessing = input(
                        "\n Enter the name of the preprocessing "
                        "strategy. The given name will be used as the name of the json "
                        "file.\n"
                    )
                    filename_preprocessing = _check_filename(
                        outdir,
                        filename_preprocessing,
                        extension=".json",
                        overwrite=overwrite,
                    )
                    # Create the preprocessing configuration file
                    create_config_preprocessing(
                        outdir, filename_preprocessing, overwrite=overwrite
                    )
                    # Add preprocessing config file directory to the workflow config file
                    signals[signal].update(
                        {
                            "preprocessing_strategy": os.path.join(
                                outdir, filename_preprocessing
                            )
                        }
                    )
                else:
                    if preprocessing.split(".")[0] in PREPROCESSING_STRATEGIES:
                        signals[signal].update(
                            {"preprocessing_strategy": preprocessing.split(".")[0]}
                        )
                    else:
                        filename_preprocessing = _check_filename(
                            outdir, preprocessing, extension=".json", overwrite=overwrite
                        )
                        signals[signal].update(
                            {
                                "preprocessing_strategy": os.path.join(
                                    outdir, filename_preprocessing
                                )
                            }
                        )
                # Add qa strategy to the workflow
                while qa is False:
                    qa = input(
                        "\n Enter the name of the qa "
                        f"strategy to clean the {signal} signal. Choose among the "
                        "current configuration files by providing the name of the "
                        "strategy, \nspecify the path to an existing file, "
                        "or create a new configuration file. To create a "
                        "new configuration file type `new`.\n Otherwise, choose among "
                        f"those strategy: {', '.join(qa_strategy[:-1])}.\n"
                    )

                if qa == "new":
                    filename_qa = input(
                        "\n Enter the name of the qa "
                        "strategy. The given name will be used as the name of the json "
                        "file.\n"
                    )
                    filename_qa = _check_filename(
                        outdir,
                        filename_qa,
                        extension=".json",
                        overwrite=overwrite,
                    )
                    # Create the qa configuration file
                    create_config_qa(outdir, filename_qa, overwrite=overwrite)
                    # Add qa config file directory to the workflow config file
                    signals[signal].update(
                        {"qa_strategy": os.path.join(outdir, filename_qa)}
                    )
                else:
                    if qa.split(".")[0] in QA_STRATEGIES:
                        signals[signal].update({"qa_strategy": qa.split(".")[0]})
                    else:
                        filename_qa = _check_filename(
                            outdir, qa, extension=".json", overwrite=overwrite
                        )
                        signals[signal].update(
                            {"qa_strategy": os.path.join(outdir, filename_qa)}
                        )

        else:
            # Save the configuration file only if there is at least one signal
            if bool(signals):
                print("\n---Saving configuration file---")
                with open(os.path.join(outdir, filename), "w") as f:
                    json.dump(signals, f, indent=4)
            break


def get_config(strategy_name, strategy="workflow"):
    """
    Get the strategy configuration file.

    Parameters
    ----------
    strategy_name: str, pathlib.Path
        Name of the workflow_strategy if using a preset or path to the configuration file
        if using a custom workflow strategy.
    strategy: str
        Type of strategy to load. Choose among `workflow`, `preprocessing` or `qa`.
        Default: `workflow`.

    Returns
    -------
    load_strategy: dict
        Dictionary containing the strategy configuration.
    """
    if strategy == "workflow":
        valid_strategies = WORKFLOW_STRATEGIES
        preset_path = "workflow_strategy"
    elif strategy == "preprocessing":
        valid_strategies = PREPROCESSING_STRATEGIES
        preset_path = "preprocessing_strategy"
    elif strategy == "qa":
        valid_strategies = QA_STRATEGIES
        preset_path = "qa_strategy"
    else:
        raise ValueError(
            "The given strategy is not valid. Choose among `workflow` or `preprocessing`."
        )

    # Check if the given strategy is valid
    if strategy_name in valid_strategies:
        strategy_path = resource_filename(
            "physprep", f"data/{preset_path}/{strategy_name}.json"
        )
    elif Path(strategy_name).exists():
        strategy_path = Path(strategy_name)
    else:
        raise ValueError(
            f"The given `strategy_name` {strategy_name} is not valid. "
            f"Please choose among {', '.join(valid_strategies)} or enter a valid path"
        )

    load_strategy = load_json(strategy_path)
    # Specific check for the workflow strategy
    if strategy == "workflow":
        if "trigger" not in load_strategy.keys():
            raise ValueError(
                "The workflow strategy configuration file must contain a key `trigger`."
            )

    return load_strategy


def rename_in_bids(data):
    """
    Rename the columns/keys in BIDS format.

    Parameters
    ----------
    data: dict or DataFrame
        Dictionary or DataFrame containing the data.

    Returns
    -------
    data: dict or DataFrame
        Dictionary or DataFrame containing the data with the new key/column names.
    """
    # If the data is a DataFrame, rename the columns according to the snake_case
    # convention
    bids_names = {}
    names = data.columns if isinstance(data, pd.DataFrame) else data
    if isinstance(data, pd.DataFrame):
        # Rename columns following BIDS convention for tabular files

        for col in names:
            col_snake = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", col)
            # Deal with multiple consecutive uppercase letters
            col_snake = re.sub("([a-z0-9])([A-Z])", r"\1_\2", col_snake).lower()
            col_snake = col_snake.replace("__", "_")
            bids_names.update({col: col_snake})

    # If the data is a dictionary, rename the keys according to the CamelCase convention
    elif isinstance(data, dict):
        # Rename keys following BIDS convention for Key-value files
        for k in names:
            if _is_camel_case(k):
                if k[0].lower() == k[0]:
                    bids_names.update({k: k[0].upper() + k[1:]})
                else:
                    bids_names.update({k: k})
            elif not _is_camel_case(k):
                key_camel = k.split("_")
                key_camel = "".join(map(str.capitalize, key_camel))
                bids_names.update({k: key_camel})

    # Rename columns/keys
    if isinstance(data, pd.DataFrame):
        data.rename(columns=bids_names, inplace=True)
    elif isinstance(data, dict):
        data = dict((bids_names[k], v) for (k, v) in data.items())

    return data


def _is_camel_case(input):
    """
    Helper function to check if the input is in CamelCase format.

    input: str
        Input string to check.
    """
    pattern = r"^[a-zA-Z][a-zA-Z0-9]*$"

    # Use re.match to see if the entire string matches the pattern
    return bool(re.match(pattern, input))


def create_scans_file_eeg(path, sub, ses, overwrite=False):
    """
    Generate the *scans files

    Parameters
    ----------
    path: str or pathlib.Path
        BIDS directory
    sub: str
        sub id (e.g. `01`)
    ses: str
        ses id (e.g. `001`)
    overwrite: bool
        Whether or not overwriting existing *scans files
    """
    filenames, acq_times = [], []

    path = Path(path)

    if (
        os.path.isfile(
            path / f"sub-{sub}" / f"ses-{ses}" / f"sub-{sub}_ses-{ses}_scans.tsv"
        )
        and not overwrite
    ):
        warnings.warn(
            f"sub-{sub}_ses-{ses}_scans.tsv already exists. If you want to"
            "overwrite that file, please specify `overwrite=True`"
        )
    else:
        # Define BIDS Layout
        layout = _check_bids_validity(path)

        # Get the files to list in the *scans file
        files = layout.get(subject=sub, session=ses, suffix="eeg", extension="edf")

        if len(files) == 0:
            raise ValueError(f"No *eeg.edf files found for sub-{sub}_ses-{ses}")
        else:
            for file in files:
                # Adding the filename of `file` to the list
                filenames.append(f"eeg/{file.filename}")
                # Get the acquisition time
                from mne.io import read_raw_edf

                tmp = read_raw_edf(file)
                acq_times.append(tmp.info["meas_date"].strftime("%Y-%m-%dT%H:%M:%S"))

        # Saving the *scans.json sidecar
        scans_json = {
            "filename": {"Description": "Name of the edf file"},
            "acq_time": {
                "LongName": "Acquisition time",
                "Description": "Acquisition time of the particular scan",
            },
        }
        with open(
            path / f"sub-{sub}" / f"ses-{ses}" / f"sub-{sub}_ses-{ses}_scans.json", "w"
        ) as f:
            json.dump(scans_json, f, indent=4)

        # Saving the *scans.tsv file
        scans_tsv = pd.DataFrame({"filename": filenames, "acq_time": acq_times})
        scans_tsv.to_csv(
            path / f"sub-{sub}" / f"ses-{ses}" / f"sub-{sub}_ses-{ses}_scans.tsv",
            sep="\t",
            index=False,
        )
