# -*- coding: utf-8 -*-
# !/usr/bin/env python
"""
Neuromod processing utilities.
"""

import timeit
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from neurokit2 import (
    ecg_peaks,
    eda_process,
    ppg_findpeaks,
    rsp_process,
    signal_fixpeaks,
    signal_rate,
)
from neurokit2.misc import as_vector
from neurokit2.signal.signal_formatpeaks import _signal_from_indices
from systole.correction import correct_peaks, correct_rr
from systole.utils import input_conversion

from physprep.utils import load_json, rename_in_bids, save_processing


def features_extraction_workflow(
    data, metadata, workflow_strategy, outdir=None, filename=None, save=True
):
    """
    Extract features from physiological data.

    Parameters
    ----------
    data : DataFrame or dict
        The physiological timeseries on which to extract features.
    metadata : dict or pathlib.Path
        The metadata associated with the cleaned physiological data, if data has been
        cleaned. Otherwise, the metadata associated with the raw physiological data
        (i.e., the outputted json file from Phys2Bids).
    workflow_strategy : dict
        Dictionary containing the content of the workflow strategy.
    outdir : str or pathlib.Path
        Path to the directory where the preprocessed physiological data will be saved.
    filename : str or pathlib.Path
        Name of the file to save the preprocessed physiological data.
    save : bool
        Specify if the preprocessed signals should be saved.
        Default to True.

    Returns
    -------
    timeseries
    info_dict
    """
    timeseries = {}
    info_dict = {}

    print("Extracting features from physiological data...\n")
    # Load metadata
    if isinstance(metadata, str) or isinstance(metadata, Path):
        metadata = load_json(metadata)

    # if data is a dataframe, convert to dict
    if isinstance(data, pd.DataFrame):
        data = data.to_dict(index="list")

    # Extract features for each signal type in the `workflow_strategy`
    for signal_type in workflow_strategy:
        if signal_type != "trigger":
            # Retrieve SamplingFrequency
            sampling_rate = metadata[signal_type]["SamplingFrequency"]
            # Extract features for each signal type
            if signal_type in data.keys():
                signal = data[signal_type]
            elif workflow_strategy[signal_type]["id"] in data.columns:
                signal = data[workflow_strategy[signal_type]["id"]]
            else:
                raise ValueError(f"Signal type {signal_type} not found in the data.")

            signal = as_vector(signal[f"{signal_type}_clean"])
            print(f"Extracting features for {signal_type}...\n")
            start_time = timeit.default_timer()

            if signal_type.lower() in ["ecg", "cardiac_ecg", "ppg", "cardiac_ppg"]:
                timeserie, info = extract_cardiac(
                    signal,
                    sampling_rate=sampling_rate,
                    data_type=signal_type,
                )
            elif signal_type.lower() in ["respiratory", "rsp", "resp", "breathing"]:
                timeserie, info = extract_respiratory(signal, sampling_rate=sampling_rate)
            elif signal_type.lower() in ["electrodermal", "eda"]:
                timeserie, info = extract_electrodermal(
                    signal, sampling_rate=sampling_rate
                )

            timeseries.update({signal_type: timeserie.to_dict("list")})
            info_dict.update({signal_type: info})

            end_time = np.round(timeit.default_timer() - start_time, 2)
            print(f"{signal_type} features extraction: done in {end_time} sec***\n")

    # Save derivatives
    if save:
        print("Saving extracted features...\n")
        save_processing(outdir, filename, "desc-features", timeseries, info_dict)
        # Save timeseries
        for timeserie in timeseries:
            timeseries[timeserie] = pd.DataFrame(timeseries[timeserie])
        print("Extracted features saved. \n")

    return timeseries, info_dict


# ==================================================================================
# Features extraction functions
# ==================================================================================


def extract_cardiac(signal, sampling_rate=1000, data_type="ppg"):
    """
    Process cardiac signal.

    Extract features of interest for neuromod project.

    Parameters
    ----------
    signal_cleaned : array
        Cleaned PPG/ECG signal.
    sampling_rate : int
        The sampling frequency of `signal` (in Hz, i.e., samples/second).
        Defaults to 1000.
    data_type : str
        Precise the type of signal to be processed (`ppg` or `ecg`).
        Default to 'ppg'.

    Returns
    -------
    timeseries : pd.DataFrame
        A DataFrame containing the cleaned PPG/ECG signals.
        - the raw signal.
        - the cleaned signal.
        - the heart rate as measured based on PPG/ECG peaks.
        - the PPG/ECG peaks marked as "1" in a list of zeros.
    info : dict
        Dictionary containing list of intervals between peaks.
    """
    try:
        # Find peaks
        print("Neurokit processing started")
        if data_type.lower() in ["ecg", "cardiac_ecg"]:
            _, info = ecg_peaks(
                ecg_cleaned=signal,
                sampling_rate=sampling_rate,
                method="promac",
                correct_artifacts=False,
            )
            # Rename Peaks key
            info["Peaks"] = info["ECG_R_Peaks"]
            # Remove duplicated info
            del info["ECG_R_Peaks"]
            del info["sampling_rate"]
            # Correct peaks
            info["CleanPeaksNK"] = signal_fixpeaks(
                peaks=info["Peaks"],
                sampling_rate=sampling_rate,
                interval_min=0.5,
                interval_max=1.5,
                method="neurokit",
            )
        elif data_type.lower() in ["ppg", "cardiac_ppg"]:
            info = ppg_findpeaks(signal, sampling_rate=sampling_rate, method="elgendi")
            info["Peaks"] = info["PPG_Peaks"]
            # Rename Peaks key
            del info["PPG_Peaks"]
            print("Neurokit found peaks")
            info["CleanPeaksNK"] = signal_fixpeaks(
                info["Peaks"],
                sampling_rate=sampling_rate,
                interval_min=0.5,
                interval_max=1.5,
                method="neurokit",
            )
            print("Neurokit fixed peaks")

        else:
            raise ValueError("Please use a valid data type: 'ecg' or 'ppg'")

        print("Formatting peaks into signal")
        peak_list_nk = _signal_from_indices(info["Peaks"], desired_length=len(signal))
        print("Formatting Peaks signal into RR timeseries")
        # peak to intervals
        rr = input_conversion(
            info["Peaks"],
            input_type="peaks_idx",
            output_type="rr_ms",
            sfreq=sampling_rate,
        )

        # correct beat detection
        corrected, (nMissed, nExtra, nEctopic, nShort, nLong) = correct_rr(rr)
        corrected_peaks = correct_peaks(peak_list_nk, n_iterations=4)
        print("systole corrected RR and Peaks series")
        # Compute rate based on peaks
        rate = signal_rate(
            info["CleanPeaksNK"],
            sampling_rate=sampling_rate,
            desired_length=len(signal),
        )

        # sanitize info dict
        info.update(
            {
                "Ectopic": nEctopic,
                "Short": nShort,
                "Long": nLong,
                "Extra": nExtra,
                "Missed": nMissed,
                "CleanRRSystole": corrected.tolist(),
                "SamplingFrequency": sampling_rate,
            }
        )
        # Prepare output
        timeseries = pd.DataFrame(
            {
                f"{data_type.lower()}_clean": signal,
                f"{data_type.lower()}_peaks_nk": peak_list_nk,
                f"{data_type.lower()}_peaks_systole": corrected_peaks["clean_peaks"],
                f"{data_type.lower()}_rate": rate,
            }
        )
        # Renaming cols/keys
        timeseries = rename_in_bids(timeseries)
        info = rename_in_bids(info)

    except Exception:
        print(f"ERROR in {data_type} features extraction procedure")
        traceback.print_exc()
        timeseries = pd.DataFrame({f"{data_type.lower()}": signal})
        info = {"Processed": False}

    return timeseries, info


def extract_electrodermal(signal, sampling_rate=1000, method="neurokit"):
    """
    Process EDA signal.

    Custom processing for neuromod EDA acquisition.

    Parameters
    -----------
    signal : vector
        The EDA timeseries.
    sampling_rate : int
        The sampling frequency of `eda_raw` (in Hz, i.e., samples/second).
        Default to 1000.

    Returns
    -------
    timeseries : pd.DataFrame
        A DataFrame containing the cleaned eda signals.
        - *"EDA_Raw"*: the raw signal.
        - *"EDA_Clean"*: the clean signal.
        - *"EDA_Tonic"*: the tonic component of the EDA signal.
        - *"EDA_Phasic"*: the phasic component of the EDA signal.
    info : dict
        Dictionary containing a list of SCR peaks.
    """
    try:
        timeseries, info = eda_process(signal, sampling_rate=sampling_rate, method=method)
        info.update({"SamplingFrequency": sampling_rate})
        # Remove duplicated info
        del info["sampling_rate"]
        # Renaming cols/keys
        timeseries = rename_in_bids(timeseries)
        info = rename_in_bids(info)
    except Exception:
        print("ERROR in EDA features extraction procedure")
        traceback.print_exc()
        timeseries = pd.DataFrame({"eda": signal})
        info = {"Processed": False}

    for k in info.keys():
        if isinstance(info[k], np.ndarray):
            info[k] = info[k].tolist()

    return timeseries, info


def extract_respiratory(signal, sampling_rate=1000, method="khodadad2018"):
    """
    Parameters
    ----------
    signal : vector
        The RESP timeseries.
    sampling_rate : int
        The sampling frequency of `eda_raw` (in Hz, i.e., samples/second).
        Default to 10000.
    method : str
        Method to use for processing.
        Default to 'khodadad2018'.

    Returns
    -------
    timeseries : pd.DataFrame
        A DataFrame containing the cleaned RSP signals.
    info : dict
        Dictionary containing a list of peaks and troughts.

    References
    ----------
    Khodadad, D., Nordebo, S., Müller, B., Waldmann, A., Yerworth, R., Becher, T., ...
        & Bayford, R. (2018). Optimized breath detection algorithm in electrical impedance
        tomography. Physiological measurement, 39(9), 094001.

    See also
    --------
    https://neuropsychology.github.io/NeuroKit/functions/rsp.html#preprocessing
    """
    try:
        timeseries, info = rsp_process(signal, sampling_rate=sampling_rate, method=method)
        info.update({"SamplingFrequency": sampling_rate})
        # Remove duplicated info
        del info["sampling_rate"]
        # Renaming cols/keys
        timeseries = rename_in_bids(timeseries)
        info = rename_in_bids(info)
    except Exception:
        print("ERROR in RESP features extraction procedure")
        traceback.print_exc()
        timeseries = pd.DataFrame({"resp": signal})
        info = {"Processed": False}

    for k in info.keys():
        if isinstance(info[k], np.ndarray):
            info[k] = info[k].tolist()

    return timeseries, info
