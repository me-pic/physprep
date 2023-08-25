# -*- coding: utf-8 -*-
# !/usr/bin/env python
"""
Neuromod processing utilities.
"""
import glob
import multiprocessing

# dependencies
import os
import pickle

# signal utils
import timeit
import traceback

import click
import neurokit2 as nk
import numpy as np
import pandas as pd

# home brewed cleaning utils
from clean import (
    neuromod_ecg_clean,
    neuromod_eda_clean,
    neuromod_ppg_clean,
    neuromod_rsp_clean,
)

# high-level processing utils
from neurokit2 import ecg_peaks, ppg_findpeaks, signal_fixpeaks, signal_rate
from neurokit2.misc import as_vector
from neurokit2.signal.signal_formatpeaks import _signal_from_indices
from systole.correction import correct_peaks, correct_rr
from systole.utils import input_conversion
from utils import load_json


def process_session(args):
    # Unpack arguments for parallel processing
    source, sub, ses, outdir, multi_echo = args
    multi_echo = bool(multi_echo)
    neuromod_bio_process(source, sub, ses, outdir, False)


# ==================================================================================
# Processing pipeline
# ==================================================================================


# @click.command()
# @click.argument("source", type=str)
# @click.argument("sub", type=str)
# @click.argument("ses", type=str)
# @click.argument("outdir", type=str)
# @click.argument("multi_echo", type=bool)
def neuromod_bio_process(source, sub, ses, outdir, multi_echo):
    """
    Run processing pipeline on specified biosignals.

    Parameters
    ----------
    source : str
        The main directory containing the segmented runs.
    sub : str
        The id of the subject.
    ses : str
        The id of the session.
    outdir : str
        The directory to save the outputs.
    multi_echo : bool
        Indicate if the multi-echo sequence was used or not.
    """
    # Check if `outdir` exists, otherwise create it

    if not os.path.exists(os.path.join(outdir, sub, ses)):
        os.makedirs(os.path.join(outdir, sub, ses))

    # Load tsv files contained in source/sub/ses
    data_tsv, data_json, filenames_tsv = load_segmented_runs(
        source, sub, ses, outdir, remove_padding=True
    )
    sampling_rate = data_json["SamplingFrequency"]

    # initialize returned objects
    for idx, df in enumerate(data_tsv):
        run = filenames_tsv[idx].split("_")[2]
        print(f"---Processing biosignals for {sub} {ses}: {run}---")

        bio_info = {}
        bio_df = pd.DataFrame()

        # initialize each signals
        ppg_raw = df["PPG"]
        rsp_raw = df["RSP"]
        eda_raw = df["EDA"]
        ecg_raw = df["ECG"]

        # ppg
        print("******PPG workflow: begin***")
        start_time = timeit.default_timer()
        ppg, ppg_info = ppg_process(ppg_raw, sampling_rate=sampling_rate)
        bio_info["PPG"] = ppg_info
        bio_df = pd.concat([bio_df, ppg], axis=1)
        end_time = timeit.default_timer() - start_time
        print(f"***PPG workflow: done in {end_time} sec***")

        #  ecg
        print("***ECG workflow: begin***")
        start_time = timeit.default_timer()
        ecg, ecg_info = ecg_process(ecg_raw, sampling_rate=sampling_rate, me=multi_echo)
        bio_info["ECG"] = ecg_info
        bio_df = pd.concat([bio_df, ecg], axis=1)
        end_time = timeit.default_timer() - start_time
        print(f"***ECG workflow: done in {end_time} sec***")

        #  rsp
        print("***Respiration workflow: begin***")
        start_time = timeit.default_timer()
        rsp, rsp_info = rsp_process(
            rsp_raw, sampling_rate=sampling_rate, method="khodadad2018"
        )
        bio_info["RSP"] = rsp_info
        bio_df = pd.concat([bio_df, rsp], axis=1)
        end_time = timeit.default_timer() - start_time
        print(f"***Respiration workflow: done in {end_time} sec***")

        #  eda
        print("***Electrodermal activity workflow: begin***")
        start_time = timeit.default_timer()
        eda, eda_info = eda_process(eda_raw, sampling_rate=sampling_rate)
        bio_info["EDA"] = eda_info
        bio_df = pd.concat([bio_df, eda], axis=1)
        end_time = timeit.default_timer() - start_time
        print(f"***Electrodermal activity workflow: done in {end_time} sec***")

        # save outputs
        print("***Saving processed biosignals: begin***")
        start_time = timeit.default_timer()

        bio_df.to_csv(
            os.path.join(outdir, sub, ses, f"{filenames_tsv[idx]}.tsv.gz"),
            sep="\t",
            index=False,
        )
        with open(
            os.path.join(outdir, sub, ses, f"{filenames_tsv[idx]}.json"),
            "wb",
        ) as fp:
            pickle.dump(bio_info, fp, protocol=4)
            fp.close()
        end_time = timeit.default_timer() - start_time
        print(f"***Saving processed biosignals: done in {end_time} sec***")


# ==================================================================================
# Loading functions
# ==================================================================================


def load_segmented_runs(source, sub, ses, outdir, remove_padding=True):
    """
    Parameters
    ----------
    source : str
        The main directory containing the segmented runs.
    sub : str
        The id of the subject.
    ses : str
        The id of the session.
    outdir : str
        The directory to save the start padding, only if `remove_padding`
        is True.
    remove_padding : bool
        Indicate if the padding should be removed or not from the data.
        If True, the signal included in the start padding is saved in outdir.
        If no padding was used, `remove_padding`should be False.
        Default to True.

    Returns
    -------
    data_tsv : list
        List containing dataframes with the raw signal for each run.
    data_json : dict
        Dictionary containing metadata.
    filenames : list
        List containing the filename for each segmented run without the extension.
    """
    data_tsv, filenames = [], []
    files_tsv = [f for f in os.listdir(os.path.join(source, sub, ses)) if "tsv.gz" in f]
    files_tsv.sort()
    # The sampling rate will be overwritten if there is a json file
    # TODO: Change harded coded sampling rate to a variable from info file
    data_json = {"SamplingFrequency": 10000, "Columns": []}

    for tsv in files_tsv:
        # Signals
        filename = tsv.split(".")[0]
        filenames.append(filename)
        print(f"---Reading data for {sub} {ses}: {filename.split('_')[2]}---")
        # Metadata
        json_file = filename + ".json"
        print("Reading json file")
        data_json = load_json(os.path.join(source, sub, ses, json_file))

        # look at the info file to get the channels if unboundlocalerror
        if data_json["Columns"] == []:
            info = pd.read_json(
                os.path.join(source, sub, f"{sub}_volumes_all-ses-runs.json")
            )
            ch_list = [
                "EDA"
                if "EDA" in ch
                else "ECG"
                if "ECG" in ch
                else "PPG"
                if "PPG" in ch
                else "TTL"
                if "A 5" in ch or "TTL" in ch
                else "RSP"
                for ch in info[ses]["ch_names"]
            ]
            data_json["Columns"] = ch_list
        print("Reading tsv file")
        # Read signals
        tmp_csv = pd.read_csv(
            os.path.join(source, sub, ses, tsv),
            sep="\t",
            compression="gzip",
            names=data_json["Columns"],
        )
        # Remove padding and save it for later
        if remove_padding:
            # Get triggers
            trigger = tmp_csv[tmp_csv["TTL"] > 4]
            # First trigger
            start = list(trigger.index)[0]
            # Last trigger
            end = list(trigger.index)[-1]

            data_tsv.append(tmp_csv[start:end].reset_index(drop=True))
            tmp_csv[:start].to_csv(
                os.path.join(outdir, sub, ses, f"{filename}_noseq.tsv.gz"),
                sep="\t",
                index=False,
            )
        # make a list of the runs
        else:
            data_tsv.append(tmp_csv)

    return data_tsv, data_json, filenames


# ==================================================================================
# Processing functions
# ==================================================================================


def process_cardiac(signal_raw, signal_cleaned, sampling_rate=10000, data_type="PPG"):
    """
    Process cardiac signal.

    Extract features of interest for neuromod project.

    Parameters
    ----------
    signal_raw : array
        Raw PPG/ECG signal.
    signal_cleaned : array
        Cleaned PPG/ECG signal.
    sampling_rate : int
        The sampling frequency of `signal_raw` (in Hz, i.e., samples/second).
        Defaults to 10000.
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
    if data_type in ["ecg", "ECG"]:
        _, info = ecg_peaks(
            ecg_cleaned=signal_cleaned,
            sampling_rate=sampling_rate,
            method="promac",
            correct_artifacts=False,
        )
        info["ECG_Peaks"] = info["ECG_R_Peaks"]
        info["ECG_Clean_Peaks_NK"] = nk.signal_fixpeaks(
            peaks=info["ECG_Peaks"],
            sampling_rate=sampling_rate,
            interval_min=0.5,
            interval_max=1.5,
            method="neurokit",
        )
    elif data_type in ["ppg", "PPG"]:
        info = ppg_findpeaks(
            signal_cleaned, sampling_rate=sampling_rate, method="elgendi"
        )
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
        info[f"{data_type.upper()}_Peaks"], desired_length=len(signal_cleaned)
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
        desired_length=len(signal_cleaned),
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
            f"{data_type.upper()}_Raw": signal_raw,
            f"{data_type.upper()}_Clean": signal_cleaned,
            f"{data_type.upper()}_Peaks_NK": peak_list_nk,
            f"{data_type.upper()}_Peaks_Systole": corrected_peaks["clean_peaks"],
            f"{data_type.upper()}_Rate": rate,
        }
    )

    return signals, info


def ppg_process(ppg_raw, sampling_rate=10000, downsampling_rate=1000):
    """
    Process PPG signal.

    Custom processing function for neuromod PPG acquisition.

    Parameters
    ----------
    ppg_raw : vector
        The raw PPG channel.
    sampling_rate : int
        The sampling frequency of `ppg_raw` (in Hz, i.e., samples/second).
        Default to 10000.
    downsampling_rate : int
        The sampling frequency to use to downsample the signals
        (in Hz, i.e., samples/second). If None, the signals are not downsampled.
        Default to 1000.

    Returns
    -------
    signals : DataFrame
        A DataFrame containing the cleaned ppg signals.
        - *"PPG_Raw"*: the raw signal.
        - *"PPG_Clean"*: the cleaned signal.
        - *"PPG_Rate"*: the heart rate as measured based on PPG peaks.
        - *"PPG_Peaks"*: the PPG peaks marked as "1" in a list of zeros.
    info : dict
        Dictionary containing list of intervals between peaks.
    """
    ppg_signal = as_vector(ppg_raw)

    # Prepare signal for processing
    print("Cleaning PPG")
    ppg_signal, ppg_cleaned = neuromod_ppg_clean(
        ppg_signal, sampling_rate=sampling_rate, downsampling=downsampling_rate
    )
    print("PPG Cleaned")
    # Process clean signal
    try:
        if downsampling_rate is None:
            signals, info = process_cardiac(
                ppg_signal,
                ppg_cleaned,
                sampling_rate=sampling_rate,
                data_type="PPG",
            )
            info["SamplingFrequency"] = sampling_rate
        else:
            signals, info = process_cardiac(
                ppg_signal,
                ppg_cleaned,
                sampling_rate=downsampling_rate,
                data_type="PPG",
            )
            info["SamplingFrequency"] = downsampling_rate
    except Exception:
        print("ERROR in PPG processing procedure")
        traceback.print_exc()
        signals = pd.DataFrame({"PPG_Raw": ppg_signal, "PPG_Clean": ppg_cleaned})
        info = {"Processed": False}

    return signals, info


def ecg_process(
    ecg_raw,
    sampling_rate=10000,
    downsampling_rate=1000,
    method="bottenhorn",
    me=True,
):
    """
    Process ECG signal.

    Custom processing for neuromod ECG acquisition.

    Parameters
    -----------
    ecg_raw : vector
        The raw ECG channel.
    sampling_rate : int
        The sampling frequency of `ecg_raw` (in Hz, i.e., samples/second).
        Default to 10000.
    downsampling_rate : int
        The sampling frequency to use to downsample the signals
        (in Hz, i.e., samples/second). If None, the signals are not downsampled.
        Default to 1000.
    method : str
        The processing pipeline to apply.
        Default to 'bottenhorn'.
    mr : bool
        True if MB-ME sequence was used. Otherwise, considered that the MB-SE
        sequence was used.
        Default to True.

    Returns
    -------
    signals : DataFrame
        A DataFrame containing the cleaned ecg signals.
        - *"ECG_Raw"*: the raw signal.
        - *"ECG_Clean"*: the cleaned signal.
        - *"ECG_Rate"*: the heart rate as measured based on ECG peaks.
    info : dict
        Dictionary containing list of peaks.
    """
    ecg_signal = as_vector(ecg_raw)

    # Prepare signal for processing
    print("Cleaning ECG")
    ecg_signal, ecg_cleaned = neuromod_ecg_clean(
        ecg_signal,
        sampling_rate=sampling_rate,
        method=method,
        me=me,
        downsampling=downsampling_rate,
    )
    print("ECG Cleaned")
    # Process clean signal
    try:
        if downsampling_rate is None:
            signals, info = process_cardiac(
                ecg_signal,
                ecg_cleaned,
                sampling_rate=sampling_rate,
                data_type="ECG",
            )
            info["SamplingFrequency"] = sampling_rate
        else:
            signals, info = process_cardiac(
                ecg_signal,
                ecg_cleaned,
                sampling_rate=downsampling_rate,
                data_type="ECG",
            )
            info["SamplingFrequency"] = downsampling_rate
    except Exception:
        print("ERROR in ECG processing procedure")
        traceback.print_exc()
        signals = pd.DataFrame({"ECG_Raw": ecg_signal, "ECG_Clean": ecg_cleaned})
        info = {"Processed": False}

    return signals, info


def eda_process(eda_raw, sampling_rate=10000, downsampling_rate=1000):
    """
    Process EDA signal.

    Custom processing for neuromod EDA acquisition.

    Parameters
    -----------
    eda_raw : vector
        The raw EDA channel.
    sampling_rate : int
        The sampling frequency of `eda_raw` (in Hz, i.e., samples/second).
        Default to 10000.
    downsampling_rate : int
        The sampling frequency to use to downsample the signals
        (in Hz, i.e., samples/second). If None, the signals are not downsampled.
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
    eda_signal = as_vector(eda_raw)

    # Prepare signal for processing
    print("Cleaning EDA")
    eda_signal, eda_cleaned = neuromod_eda_clean(
        eda_signal, sampling_rate=sampling_rate, downsampling=downsampling_rate
    )
    print("EDA Cleaned")
    # Process clean signal
    try:
        if downsampling_rate is None:
            signals, info = nk.eda_process(
                eda_cleaned, sampling_rate=sampling_rate, method="neurokit"
            )
            info["SamplingFrequency"] = sampling_rate
        else:
            signals, info = nk.eda_process(
                eda_cleaned, sampling_rate=downsampling_rate, method="neurokit"
            )
            info["SamplingFrequency"] = downsampling_rate
        signals["EDA_Raw"] = eda_signal
    except Exception:
        print("ERROR in EDA processing procedure")
        traceback.print_exc()
        signals = pd.DataFrame({"EDA_Raw": eda_signal, "EDA_Clean": eda_cleaned})
        info = {"Processed": False}

    for k in info.keys():
        if isinstance(info[k], np.ndarray):
            info[k] = info[k].tolist()

    return signals, info


def rsp_process(
    rsp_raw, sampling_rate=10000, downsampling_rate=1000, method="khodadad2018"
):
    """
    Parameters
    ----------
    rsp_raw : vector
        The raw RSP channel
    sampling_rate : int
        The sampling frequency of `eda_raw` (in Hz, i.e., samples/second).
        Default to 10000.
    downsampling_rate : int
        The sampling frequency to use to downsample the signals
        (in Hz, i.e., samples/second). If None, the signals are not downsampled.
        Default to 1000.
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
    rsp_signal = as_vector(rsp_raw)

    # Clean and filter respiratory signal
    print("Cleaning RSP")
    rsp_signal, rsp_cleaned = neuromod_rsp_clean(
        rsp_signal, sampling_rate=sampling_rate, downsampling=downsampling_rate
    )
    # Process clean signal
    print("Processing RSP")
    try:
        if downsampling_rate is None:
            signals, info = nk.rsp_process(
                rsp_cleaned, sampling_rate=sampling_rate, method=method
            )
            info["SamplingFrequency"] = sampling_rate
        else:
            signals, info = nk.rsp_process(
                rsp_cleaned, sampling_rate=downsampling_rate, method=method
            )
            info["SamplingFrequency"] = downsampling_rate
        signals["RSP_Raw"] = rsp_signal
        print("RSP Cleaned and processed")
    except Exception:
        print("ERROR in RSP processing procedure")
        traceback.print_exc()
        signals = pd.DataFrame({"RSP_Raw": rsp_signal, "RSP_Clean": rsp_cleaned})
        info = {"Processed": False}

    for k in info.keys():
        if isinstance(info[k], np.ndarray):
            info[k] = info[k].tolist()

    return signals, info


@click.command()
@click.argument("source", type=str)
@click.argument("outdir", type=str)
@click.argument("sub", type=str)
@click.argument("multi_echo", type=bool)
def parallel_neuromod_bio_process(source, sub, outdir, multi_echo):
    num_cpus = multiprocessing.cpu_count()
    num_cpus = round(num_cpus / 3)
    pool = multiprocessing.Pool(processes=num_cpus)
    sessions = [
        os.path.basename(name) for name in glob.glob(os.path.join(source, sub, "ses-*"))
    ]
    print(sessions)
    args = [(source, sub, ses, outdir, multi_echo) for ses in sorted(sessions)]
    pool.map(process_session, args)

    pool.close()
    pool.join()


if __name__ == "__main__":
    # neuromod_bio_process()
    parallel_neuromod_bio_process()
