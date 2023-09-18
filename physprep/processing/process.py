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
            try:
                df, info = extract_cardiac(
                    signal,
                    sampling_rate=sampling_rate,
                    data_type=signal_type.id,
                )
            except Exception:
                print(f"ERROR in {signal_type.id} features extraction procedure")
                traceback.print_exc()
                signals = pd.DataFrame({f"{signal_type.id}": signal})
                info = {"Processed": False}

        elif signal_type.id == "RESP":
            df, info = extract_respiratory(signal, sampling_rate=sampling_rate)
        elif signal_type.id == "EDA":
            df, info = extract_electrodermal(signal, sampling_rate=sampling_rate)

        end_time = timeit.default_timer() - start_time
        print(f"***{signal_type.id} features extraction: done in {end_time} sec***\n")

    # Save derivatives
    if save:
        pass
    # bio_df.to_csv(
    #     os.path.join(outdir, sub, ses, f"{filenames_tsv[idx]}.tsv.gz"),
    #     sep="\t",
    #     index=False,
    # )
    # with open(
    #     os.path.join(outdir, sub, ses, f"{filenames_tsv[idx]}.json"),
    #     "wb",
    # ) as f:
    #     pickle.dump(bio_info, f, protocol=4)
    #     f.close()

    return df, info, signals


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
    signals : DataFrame
        A DataFrame containing the cleaned PPG/ECG signals.
        - the raw signal.
        - the cleaned signal.
        - the heart rate as measured based on PPG/ECG peaks.
        - the PPG/ECG peaks marked as "1" in a list of zeros.
    info : dict
        Dictionary containing list of intervals between peaks.
    """
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
            f"{data_type.upper()}_ectopic": nEctopic,
            f"{data_type.upper()}_short": nShort,
            f"{data_type.upper()}_long": nLong,
            f"{data_type.upper()}_extra": nExtra,
            f"{data_type.upper()}_missed": nMissed,
            f"{data_type.upper()}_Clean_RR_Systole": corrected.tolist(),
        }
    )
    # Prepare output
    signals = pd.DataFrame(
        {
            f"{data_type.upper()}_Clean": signal,
            f"{data_type.upper()}_Peaks_NK": peak_list_nk,
            f"{data_type.upper()}_Peaks_Systole": corrected_peaks["clean_peaks"],
            f"{data_type.upper()}_Rate": rate,
        }
    )

    # Remove duplicated info, that info should already be available in the metadata
    # Under the key SamplingFrequency as recommended in the BIDS specification.
    del info["sampling_rate"]

    return signals, info


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
        signals, info = eda_process(signal, sampling_rate=sampling_rate, method=method)
    except Exception:
        print("ERROR in EDA features extraction procedure")
        traceback.print_exc()
        signals = pd.DataFrame({"EDA": signal})
        info = {"Processed": False}

    for k in info.keys():
        if isinstance(info[k], np.ndarray):
            info[k] = info[k].tolist()

    return signals, info


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
    signals : DataFrame
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
        signals, info = rsp_process(signal, sampling_rate=sampling_rate, method=method)
    except Exception:
        print("ERROR in RESP features extraction procedure")
        traceback.print_exc()
        signals = pd.DataFrame({"RESP": signal})
        info = {"Processed": False}

    for k in info.keys():
        if isinstance(info[k], np.ndarray):
            info[k] = info[k].tolist()

    return signals, info
