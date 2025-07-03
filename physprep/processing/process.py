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
    rsp_peaks,
    signal_fixpeaks,
    signal_rate,
)
from neurokit2.misc import as_vector
from neurokit2.signal.signal_formatpeaks import _signal_from_indices

from physprep.utils import load_json, rename_in_bids, save_processing


def features_extraction_workflow(
    data, metadata, workflow_strategy
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
        (i.e., the outputed json file from Phys2Bids).
    workflow_strategy : dict
        Dictionary containing the content of the workflow strategy.

    Returns
    -------
    features
    """
    features = {}

    print("Extracting features from physiological data...\n")
    # Load metadata
    if isinstance(metadata, str) or isinstance(metadata, Path):
        metadata = load_json(metadata)

    # if data is a dataframe, convert to dict
    if isinstance(data, pd.DataFrame):
        data = data.to_dict("list")
    
    # Extract features for each signal type in the `workflow_strategy`
    for idx, signal_type in enumerate(workflow_strategy):
        if signal_type not in ["trigger", "concurrentWith"]:
            # Retrieve SamplingFrequency
            if isinstance(metadata["SamplingFrequency"], list):
                sampling_rate = metadata["SamplingFrequency"][idx]
            else:
                sampling_rate = metadata["SamplingFrequency"]

            # Extract features for each signal type
            if signal_type not in '\t'.join(data.keys()):
                raise ValueError(f"Signal type {signal_type} not found in the data.")

           
            signal_key = f"{signal_type}_clean"  if f"{signal_type}_clean" in data.keys() else signal_type

            signal = as_vector(data[signal_key])
            print(f"Extracting features for {signal_type}...\n")
            start_time = timeit.default_timer()

            if signal_type.lower() in ["ecg", "cardiac_ecg", "ppg", "cardiac_ppg"]:
                info = extract_cardiac_peaks(
                    signal,
                    sampling_rate=sampling_rate,
                    data_type=signal_type,
                )
            elif signal_type.lower() in ["respiratory", "rsp", "resp", "breathing"]:
                info = extract_respiratory_peaks(signal, sampling_rate=sampling_rate)
            elif signal_type.lower() in ["electrodermal", "eda", "gsr"]:
                info = extract_electrodermal_peaks(
                    signal, sampling_rate=sampling_rate
                )

            features.update({signal_type: info})

            end_time = np.round(timeit.default_timer() - start_time, 2)
            print(f"{signal_type} features extraction: done in {end_time} sec***\n")

    events = convert_2_events(features, metadata)

    return features, events

    
# ==================================================================================
# Features extraction functions
# ==================================================================================


def extract_cardiac_peaks(signal, sampling_rate=1000, data_type="ppg"):
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
            # Correct peaks
            info["CleanPeaksNK"] = signal_fixpeaks(
                peaks=info["ECG_R_Peaks"],
                sampling_rate=sampling_rate,
                interval_min=0.5,
                interval_max=1.5,
                method="neurokit",
            )
            corrected_peaks = local_maxima_peak_correction(
                signal, 
                info["CleanPeaksNK"], 
                window=0.02, 
                sampling_rate=sampling_rate
            )

            # Add peaks to dictionnary
            info = {
                'r_peak': info['ECG_R_Peaks'],
                'r_peak_corrected': corrected_peaks
            }
            print("Formatting peaks into signal")
            peak_list = _signal_from_indices(info['r_peak'], desired_length=len(signal))
            peak_list_corrected = _signal_from_indices(info['r_peak_corrected'], desired_length=len(signal))
        elif data_type.lower() in ["ppg", "cardiac_ppg"]:
            info = ppg_findpeaks(signal, sampling_rate=sampling_rate, method="elgendi")
            # Rename Peaks key
            print("Neurokit found peaks")
            info["CleanPeaksNK"] = signal_fixpeaks(
                info["PPG_Peaks"],
                sampling_rate=sampling_rate,
                interval_min=0.5,
                interval_max=1.5,
                method="neurokit",
            )
            print("Neurokit fixed peaks")
            info = {
                'systolic_peak': info['PPG_Peaks'],
                'systolic_peak_corrected': info['CleanPeaksNK']
            }
            print("Formatting peaks into signal")
            peak_list = _signal_from_indices(info['systolic_peak'], desired_length=len(signal))
            peak_list_corrected = _signal_from_indices(info['systolic_peak_corrected'], desired_length=len(signal))
        else:
            raise ValueError("Please use a valid data type: 'ecg' or 'ppg'")

    except Exception:
        print(f"ERROR in {data_type} features extraction procedure")
        traceback.print_exc()
        info = {"Processed": False}

    return info


def extract_electrodermal_peaks(signal, sampling_rate=1000, method="neurokit"):
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
        _, info = eda_process(signal, sampling_rate=sampling_rate, method=method)
        # Renaming cols
        info = {
            'scr_onset': info['SCR_Onsets'],
            'scr_peak': info['SCR_Peaks'],
            'scr_recovery': info['SCR_Recovery']
        }
    except Exception:
        print("ERROR in EDA features extraction procedure")
        traceback.print_exc()
        info = {"Processed": False}

    for k in info.keys():
        if isinstance(info[k], np.ndarray):
            info[k] = info[k].tolist()

    return info


def extract_respiratory_peaks(signal, sampling_rate=1000, method="khodadad2018"):
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
        _, info = rsp_peaks(signal, sampling_rate=sampling_rate, method=method)
        # Renaming cols/keys
        info = {
            'inhale_max': info['RSP_Peaks'],
            'exhale_max': info['RSP_Troughs']
        }
    except Exception:
        print("ERROR in RESP features extraction procedure")
        traceback.print_exc()
        info = {"Processed": False}

    for k in info.keys():
        if isinstance(info[k], np.ndarray):
            info[k] = info[k].tolist()


    return info


def convert_2_events(features, metadata):
    # Load json dictionnary containing non continuous data
    modalities = list(features.keys())
    duration = 0
    row=0
    df_events = pd.DataFrame(columns=['onset', 'duration', 'trial_type', 'channel'])
    for modality in modalities:
        if isinstance(metadata["SamplingFrequency"], list):
            sampling_rate = metadata["SamplingFrequency"][idx]
        else:
            sampling_rate = metadata["SamplingFrequency"]

        entities = list(features[modality].keys())
        for entity in entities:
            list_value = features[modality][entity]
            if type(list_value) is int:
                list_value = [list_value]
            elif type(list_value) is np.ndarray:
                list_value = list_value.tolist()
            list_value = [value for value in list_value if str(value) != 'nan']
                
            for idx in list_value:
                df_events.loc[row] = pd.Series({
                    'onset':idx/sampling_rate, 
                    'duration':duration, 
                    'trial_type': entity, 
                    'channel': modality
                    })
                row += 1
                        

    return df_events.sort_values(('onset'), ignore_index=True)


def local_maxima_peak_correction(signal, peaks, window=0.02, sampling_rate=1000):
    """
    signal : np.array
        Array containing the continuous timeserie on which the peaks were extracted
    peaks : np.array
        Array containing the index of the extracted peaks
    window : int
        Window to consider before and after the peak (in ms)
    sampling_rate : int
        Sampling rate of the signal
    """
    for idx, peak in enumerate(peaks):
        start = int(peak - window * sampling_rate)
        end = int (peak + window * sampling_rate)

        idx_local_maxima = np.argmax(signal[start:end])
        peaks[idx] = idx_local_maxima + start

    return peaks    