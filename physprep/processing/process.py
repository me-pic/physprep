# -*- coding: utf-8 -*-
# !/usr/bin/env python
"""
Neuromod processing utilities.
"""
# import pickle

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
from utils import load_json


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
        # Retrieve SamplingFrequency
        if isinstance(metadata["SamplingFrequency"], dict):
            sampling_rate = metadata["SamplingFrequency"][signal_type.id]
        else:
            sampling_rate = metadata["SamplingFrequency"]

        # Extract features for each signal type
        signal = as_vector(data[signal_type.id])
        print(f"***{signal_type.id} features extraction: begin***\n")
        start_time = timeit.default_timer()

        if signal_type.id in ["ECG", "PPG"]:
            timeserie, info = extract_cardiac(
                signal,
                sampling_rate=sampling_rate,
                data_type=signal_type.id,
            )
        elif signal_type.id == "RESP":
            timeserie, info = extract_respiratory(signal, sampling_rate=sampling_rate)
        elif signal_type.id == "EDA":
            timeserie, info = extract_electrodermal(signal, sampling_rate=sampling_rate)

        # Remove duplicated info, that info should already be available in the metadata
        # Under the key SamplingFrequency as recommended in the BIDS specification.
        del info["sampling_rate"]

        timeseries.update({signal_type.id: timeserie})
        info_dict.update({signal_type.id: info})

        end_time = timeit.default_timer() - start_time
        print(f"***{signal_type.id} features extraction: done in {end_time} sec***\n")

    # Save derivatives
    if save:
        # Save preprocessed signal
        if outdir is not None:
            outdir = Path(outdir)
        else:
            print(
                "WARNING! No output directory specified. Data will be saved in the "
                f"current working directory: {Path.cwd()}"
            )
            outdir = Path.cwd()
        # If signals have different lengths (i.e. different sampling rates), data will
        # be saved in different files
        vals = timeseries.values()
        ref_length = len(vals[0])
        if all(len(item) == ref_length for item in vals):
            # Save timeseries
            timeseries = pd.DataFrame(timeseries)
            timeseries.to_csv(
                Path(outdir / filename).with_suffix(".tsv.gz"),
                sep="\t",
                index=False,
                compression="gzip",
            )
            # TODO: Save info_dict; check BIDS spec for filename
        else:
            # Save timeseries
            for timeserie in timeseries:
                signal_name = list(timeserie.keys())[0]
                filename_signal = filename.replace("physio", signal_name)
                timeserie = pd.DataFrame(timeserie)
                timeserie.to_csv(
                    Path(outdir / filename_signal).with_suffix(".tsv.gz"),
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
            # TODO: Save info_dict; check BIDS spec for filename

    return timeseries, info_dict


# ==================================================================================
# Features extraction functions
# ==================================================================================


def extract_cardiac(signal, sampling_rate=1000, data_type="PPG"):
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
        Precise the type of signal to be processed. The function currently
        support PPG and ECG signal processing.
        Default to 'PPG'.

    Returns
    -------
    timeseries : dict
        A Dictionary containing the cleaned PPG/ECG signals.
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
        if data_type == "ECG":
            _, info = ecg_peaks(
                ecg_cleaned=signal,
                sampling_rate=sampling_rate,
                method="promac",
                correct_artifacts=False,
            )
            info["ECG_Peaks"] = info["ECG_R_Peaks"]
            info["ECG_Clean_Peaks_NK"] = signal_fixpeaks(
                peaks=info["ECG_Peaks"],
                sampling_rate=sampling_rate,
                interval_min=0.5,
                interval_max=1.5,
                method="neurokit",
            )
        elif data_type == "PPG":
            info = ppg_findpeaks(signal, sampling_rate=sampling_rate, method="elgendi")
            print("Neurokit found peaks")
            info["PPG_Clean_Peaks_NK"] = signal_fixpeaks(
                info["PPG_Peaks"],
                sampling_rate=sampling_rate,
                interval_min=0.5,
                interval_max=1.5,
                method="neurokit",
            )
            print("Neurokit fixed peaks")

        else:
            print("Please use a valid data type: 'ECG' or 'PPG'")

        print("Formatting peaks into signal")
        peak_list_nk = _signal_from_indices(
            info[f"{data_type.upper()}_Peaks"], desired_length=len(signal)
        )
        print("Formatting Peaks signal into RR timeseries")
        # peak to intervals
        rr = input_conversion(
            info[f"{data_type.upper()}_Peaks"],
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
            info[f"{data_type.upper()}_Clean_Peaks_NK"],
            sampling_rate=sampling_rate,
            desired_length=len(signal),
        )

        # sanitize info dict
        info.update(
            {
                f"{data_type}_ectopic": nEctopic,
                f"{data_type}_short": nShort,
                f"{data_type}_long": nLong,
                f"{data_type}_extra": nExtra,
                f"{data_type}_missed": nMissed,
                f"{data_type}_Clean_RR_Systole": corrected.tolist(),
            }
        )
        # Prepare output
        timeseries = {
            f"{data_type}_Clean": signal,
            f"{data_type}_Peaks_NK": peak_list_nk,
            f"{data_type}_Peaks_Systole": corrected_peaks["clean_peaks"],
            f"{data_type}_Rate": rate,
        }

    except Exception:
        print(f"ERROR in {data_type} features extraction procedure")
        traceback.print_exc()
        timeseries = {f"{data_type}": signal}
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
    signals : DataFrame
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
        timeseries.to_dict(orient="list")
    except Exception:
        print("ERROR in EDA features extraction procedure")
        traceback.print_exc()
        timeseries = {"EDA": signal}
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
    timeseries : dict
        A dictionary containing the cleaned RSP signals.
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
        timeseries.to_dict(orient="list")
    except Exception:
        print("ERROR in RESP features extraction procedure")
        traceback.print_exc()
        timeseries = {"RESP": signal}
        info = {"Processed": False}

    for k in info.keys():
        if isinstance(info[k], np.ndarray):
            info[k] = info[k].tolist()

    return timeseries, info
