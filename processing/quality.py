# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""Physiological data quality assessment"""

import os
import glob
import json
import click
import operator
import traceback
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import kurtosis, skew


@click.command()
@click.argument("source", type=str)
@click.argument("sub", type=str)
@click.argument("ses", type=str)
@click.argument("outdir", type=str)
@click.argument("sliding")
def neuromod_bio_sqi(source, sub, ses, outdir, sliding={"duration": 60, "step": 10}):
    """
    Run processing pipeline on specified biosignals.

    Parameters
    ----------
    source : str
        The main directory contaning the segmented runs.
    sub : str
        The id of the subject.
    ses : str
        The id of the session.
    outdir : str
        The directory to save the outputs.
    sliding : dict
    """
    filenames = glob.glob(os.path.join(source, sub, ses, "*_physio*"))
    filenames_signal = [
        f
        for f in filenames
        if "".join(Path(f).suffixes) == ".tsv.gz" and "noseq" not in f
    ]
    filenames_signal.sort()
    # breakpoint()
    for idx, f in enumerate(filenames_signal):
        filename = f.split(".")[0]
        info = load_json(os.path.join(source, sub, ses, filename + ".json"))

        signal = pd.read_csv(os.path.join(source, sub, ses, f), sep="\t")
        summary = {}
        print(f"---QCing on {f.split('/')[-1]}---")
        # Compute metrics on the unsegmented signal
        # if sliding is None:
        print("***Computing quality metrics for PPG signal***")
        try:
            summary["PPG"] = sqi_cardiac(
                signal["PPG_Clean"],
                info["PPG"],
                data_type="PPG",
                sampling_rate=info["PPG"]["sampling_rate"],
            )
        except Exception:
            print("Not able to compute PPG")
            traceback.print_exc()
        print("***Computing quality metrics for EEG signal***")
        try:
            summary["ECG"] = sqi_cardiac(
                signal["ECG_Clean"],
                info["ECG"],
                data_type="ECG",
                sampling_rate=info["ECG"]["sampling_rate"],
            )
        except Exception:
            print("Not able to compute ECG")
            traceback.print_exc()
        print("***Computing quality metrics for EDA signal***")
        try:
            summary["EDA"] = sqi_eda(
                signal, info["EDA"], sampling_rate=info["EDA"]["sampling_rate"]
            )
        except Exception:
            print("Not able to compute EDA")
            traceback.print_exc()
        print("***Computing quality metrics for RSP signal***")
        try:
            summary["RSP"] = sqi_rsp(
                signal, info["RSP"], sampling_rate=info["RSP"]["sampling_rate"]
            )
        except Exception:
            print("Not able to compute RSP")
            traceback.print_exc()
        print("***Generating report***")
        # savefile = Path(source)
        generate_report(
            summary, os.path.join(outdir, sub, ses), f"{filename.split('/')[-1]}"
        )  # .split('/')[-1].split('_')[2]}")
        # Compute metrics on signal segmented in fixed windows
        """
        elif sliding['step'] != 0:
            window_samples = int(sliding['duration'] * sampling_rate)
            num_windows = int(len(signal) / window_samples)

            for i in range(num_windows):
                start = i * window_samples
                end = start + window_samples
                window = signal.iloc([start:end])
                print("***Computing quality metrics for PPG signal***")
                summary["PPG"].append(sqi_cardiac(window["PPG_Clean"], info["PPG"], data_type="PPG"))
                print("***Computing quality metrics for EEG signal***")
                summary["ECG"].append(sqi_cardiac(window["ECG_Clean"], info["ECG"], data_type="ECG"))
                print("***Computing quality metrics for EDA signal***")
                summary["EDA"].append(sqi_eda(signal, info["EDA"]))
                print("***Computing quality metrics for RSP signal***")
                summary["RSP"].append(sqi_rsp(signal, info["RSP"]))
            print("***Generating report***")
            savefile = Path(source)
            generate_report(summary, os.path.join(outdir, sub, ses), f"{sub}_{ses}_task-{filename.parts[-2]}_run-{idx+1}_physio.html")
        # Compute metrics using a sliding window approach    
        else:
            window_samples = int(sliding['duration'] * sampling_rate)
            step_samples = int(sliding['step'] * sampling_rate)
            num_windows = int((len(signal) - window_samples) / step_samples) + 1
            for i in range(num_windows):
                start = i * step_samples
                end = start + window_samples
                window = signal.iloc([start:end])
                print("***Computing quality metrics for PPG signal***")
                summary["PPG"].append(sqi_cardiac(window["PPG_Clean"], info["PPG"], data_type="PPG"))
                print("***Computing quality metrics for EEG signal***")
                summary["ECG"].append(sqi_cardiac(window["ECG_Clean"], info["ECG"], data_type="ECG"))
                print("***Computing quality metrics for EDA signal***")
                summary["EDA"].append(sqi_eda(signal, info["EDA"]))
                print("***Computing quality metrics for RSP signal***")
                summary["RSP"].append(sqi_rsp(signal, info["RSP"]))
            print("***Generating report***")
            savefile = Path(source)
            generate_report(summary, os.path.join(outdir, sub, ses), f"{sub}_{ses}_task-{filename.parts[-2]}_run-{idx+1}_physio.html")      
        """


# ==================================================================================
# Utils
# ==================================================================================


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
    tmp = open(filename)
    data = json.load(tmp)
    tmp.close()

    return data


# ==================================================================================
# Signal Quality Indices
# ==================================================================================


def sqi_cardiac(
    signal_cardiac,
    info,
    data_type="ECG",
    sampling_rate=10000,
    mean_NN=[600, 1200],
    std_NN=300,
):
    """
    Extract SQI for ECG/PPG processed signal

    Parameters
    ----------
    signal_cardiac : DataFrame
        Output from the process.py script.
    info : dict
        Output from the process.py script.
    data_type : str
        Type of the signal. Valid options include 'ecg' and 'ppg'.
        Default to 'ecg'.
    sampling_rate : int
        The sampling frequency of `signal` (in Hz, i.e., samples/second).
        Default to 10000.
    mean_NN : list
        Range to consider to evaluate the quality of the signal based on the
        mean NN interval.
    std_NN : int or float
        Value to consider to evaluate the quality if the signal based on the
        std NN interval.

    Returns
    -------
    summary : DataFrame
        DataFrame containing sqi values.

    Reference
    ---------
    Le, V. K. D., Ho, H. B., Karolcik, S., Hernandez, B., Greeff, H.,
        Nguyen, V. H., ... & Clifton, D. (2022). vital_sqi: A Python
        package for physiological signal quality control.
        Frontiers in Physiology, 2248.
    """
    summary = {}

    # Descriptive indices on NN intervals
    summary["Mean_NN_intervals"] = np.round(
        np.mean(info[f"{data_type}_clean_rr_systole"]), 4
    )
    summary["Median_NN_intervals"] = np.round(
        np.median(info[f"{data_type}_clean_rr_systole"]), 4
    )
    summary["SD_NN_intervals"] = np.round(
        np.std(info[f"{data_type}_clean_rr_systole"], ddof=1), 4
    )
    # Descriptive indices on heart rate
    summary["Mean_HR"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"], metric="mean"
    )
    summary["Median_HR"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"], metric="median"
    )
    summary["SD_HR"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"], metric="std"
    )
    summary["Min_HR"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"], metric="min"
    )
    summary["Max_HR"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"], metric="max"
    )
    # Descriptive indices on overall signal
    summary["Skewness"] = np.round(kurtosis(signal_cardiac), 4)
    summary["Kurtosis"] = np.round(skew(signal_cardiac), 4)
    summary["Ectopic"] = info[f"{data_type}_ectopic"]
    summary["Missed"] = info[f"{data_type}_missed"]
    summary["Extra"] = info[f"{data_type}_extra"]
    summary["Long"] = info[f"{data_type}_long"]
    summary["Short"] = info[f"{data_type}_short"]
    summary["Cumulseconds_rejected"] = info[f"{data_type}_cumulseconds_rejected"]
    summary["%_rejected_segments"] = np.round(
        info[f"{data_type}_%_rejected_segments"], 4
    )

    # Quality assessment based on mean NN intervals and std
    if (
        threshold_sqi(np.mean(info[f"{data_type}_clean_rr_systole"]), mean_NN)
        == "Acceptable"
        and threshold_sqi(
            np.std(info[f"{data_type}_clean_rr_systole"], ddof=1), std_NN, operator.lt
        )
        == "Acceptable"
    ):
        summary["Quality"] = "Acceptable"
    else:
        summary["Quality"] = "Not acceptable"

    return summary


def sqi_eda(signal_eda, info, sampling_rate=10000):
    """
    Extract SQI for EDA processed signal

    Parameters
    ----------
    signal_eda : DataFrame
        Output from the process.py script.
    info : dict
        Output from the process.py script.
    sampling_rate : int
        The sampling frequency of `signal_raw` (in Hz, i.e., samples/second).
        Default to 10000.

    Returns
    -------
    summary : DataFrame
        DataFrame containing sqi values.
    """
    summary = {}
    # Descriptive indices on overall signal
    summary["Minimal_range"] = minimal_range_sqi(
        signal_eda["EDA_Clean"], threshold=0.05
    )
    summary["RAC"] = rac_sqi(signal_eda["EDA_Clean"], threshold=0.2, duration=2)
    summary["Mean_EDA"] = np.round(np.mean(signal_eda["EDA_Clean"]), 4)
    summary["Median_EDA"] = np.round(np.median(signal_eda["EDA_Clean"]), 4)
    summary["SD_EDA"] = np.round(np.std(signal_eda["EDA_Clean"]), 4)
    summary["Min_EDA"] = np.round(np.min(signal_eda["EDA_Clean"]), 4)
    summary["Max_EDA"] = np.round(np.max(signal_eda["EDA_Clean"]), 4)
    # Descriptive indices on SCL
    summary["Mean_SCL"] = np.round(np.mean(signal_eda["EDA_Tonic"]), 4)
    summary["SD_SCL"] = np.round(np.std(signal_eda["EDA_Tonic"]), 4)
    summary["Median_SCL"] = np.round(np.median(signal_eda["EDA_Tonic"]), 4)
    summary["Min_SCL"] = np.round(np.min(signal_eda["EDA_Tonic"]), 4)
    summary["Max_SC"] = np.round(np.max(signal_eda["EDA_Tonic"]), 4)
    # Descriptive indices on SCR
    summary["Mean_SCR"] = np.round(np.mean(signal_eda["EDA_Phasic"]), 4)
    summary["SD_SCR"] = np.round(np.std(signal_eda["EDA_Phasic"]), 4)
    summary["Median_SCR"] = np.round(np.median(signal_eda["EDA_Phasic"]), 4)
    summary["Min_SCR"] = np.round(np.min(signal_eda["EDA_Phasic"]), 4)
    summary["Max_SCR"] = np.round(np.max(signal_eda["EDA_Phasic"]), 4)
    summary["Number_of_detected_onsets"] = len(info["SCR_Onsets"])
    summary["Number_of_detected_peaks"] = len(info["SCR_Peaks"])
    # Descriptive indices on sCR rise time
    summary["Mean_rise_time"] = np.round(
        np.mean(rise_time_sqi(info["SCR_Onsets"], info["SCR_Peaks"], sampling_rate)), 4
    )
    summary["SD_rise_time"] = np.round(
        np.std(rise_time_sqi(info["SCR_Onsets"], info["SCR_Peaks"], sampling_rate)), 4
    )
    summary["Median_rise_time"] = np.round(
        np.median(rise_time_sqi(info["SCR_Onsets"], info["SCR_Peaks"], sampling_rate)),
        4,
    )
    summary["Min_rise_time"] = np.round(
        np.min(rise_time_sqi(info["SCR_Onsets"], info["SCR_Peaks"], sampling_rate)), 4
    )
    summary["Max_rise_time"] = np.round(
        np.max(rise_time_sqi(info["SCR_Onsets"], info["SCR_Peaks"], sampling_rate)), 4
    )

    return summary


def sqi_rsp(signal_rsp, info, sampling_rate=10000, mean_rate=0.5):
    """
    Extract SQI for respiratory processed signal

    Parameters
    ----------
    signal_rsp : DataFrame
        Output from the process.py script.
    info : dict
        Output from the process.py script.
    sampling_rate : int
        The sampling frequency of `signal_raw` (in Hz, i.e., samples/second).
        Default to 10000.
    mean_rate : int or float
        Value to consider to evaluate the quality if the signal based on the
        respiratory rate mean (in Hz).

    Returns
    -------
    summary : DataFrame
        DataFrame containing sqi values.
    """
    summary = {}
    # Descriptive indices on signal amplitude
    summary["Mean_Amp"] = np.round(np.mean(signal_rsp["RSP_Amplitude"]), 4)
    summary["Median_Amp"] = np.round(np.median(signal_rsp["RSP_Amplitude"]), 4)
    summary["SD_Amp"] = np.round(np.std(signal_rsp["RSP_Amplitude"]), 4)
    summary["Min_Amp"] = np.round(np.min(signal_rsp["RSP_Amplitude"]), 4)
    summary["CV_Amp"] = np.round(np.max(signal_rsp["RSP_Amplitude"]), 4)
    summary["variability_Amp"] = np.round(
        np.std(signal_rsp["RSP_Amplitude"]) / np.mean(signal_rsp["RSP_Amplitude"]), 4
    )
    # Descriptive indices on signal rate
    summary["Mean_Rate"] = np.round(np.mean(signal_rsp["RSP_Rate"]), 4)
    summary["Median_Rate"] = np.round(np.median(signal_rsp["RSP_Rate"]), 4)
    summary["SD_Rate"] = np.round(np.std(signal_rsp["RSP_Rate"]), 4)
    summary["Min_Rate"] = np.round(np.min(signal_rsp["RSP_Rate"]), 4)
    summary["Max_Rate"] = np.round(np.max(signal_rsp["RSP_Rate"]), 4)
    summary["CV_Rate"] = np.round(
        np.std(signal_rsp["RSP_Rate"]) / np.mean(signal_rsp["RSP_Rate"]), 4
    )
    # Quality assessment based on the mean respiratory rate
    if (
        threshold_sqi(np.mean(signal_rsp["RSP_Rate"]) / 60, mean_rate, operator.lt)
        == "Acceptable"
    ):
        summary["Quality"] = "Acceptable"
    else:
        summary["Quality"] = "Not acceptable"

    return summary


# ==================================================================================
# Quality indices : Internal metrics
# ==================================================================================


def metrics_hr_sqi(intervals, metric="mean"):
    """
    Compute the mean heart rate from the RR intervals

    Parameters
    ----------
    intervals : vector
        RR intervals.
    metric : str
        Specify the metric to use between 'mean', 'median',
        'std', 'min', 'max'.
        Default to 'mean'.

    Returns
    -------
    metric_rr : float
        Metric related to the heart rate.
    """
    bpm = np.divide(60000, intervals)

    try:
        if metric == "mean":
            metric_rr = np.round(np.mean(bpm), 4)
        elif metric == "median":
            metric_rr = np.round(np.median(bpm), 4)
        elif metric == "std":
            metric_rr = np.round(np.std(bpm), 4)
        elif metric == "min":
            metric_rr = np.round(np.min(bpm), 4)
        elif metric == "max":
            metric_rr = np.round(np.max(bpm), 4)
    except:
        print(f"Invalid metric: {metric}.")

    return metric_rr


def minimal_range_sqi(signal, threshold):
    """
    Compute the ratio between the number of timepoints under a definied `threshold`
    and the signal length.

    Parameters
    ----------
    signal : vector
        Signal on which to compute the minimal range.
    threshold : float
        Threshold to consider to compute the minimal range.

    Returns
    -------
    minimal_range : float
        Ratio between the signal under `threshold` and the overall signal.

    References
    ----------
    Böttcher, S., Vieluf, S., Bruno, E., Joseph, B., Epitashvili, N., Biondi, A., ...
        & Loddenkemper, T. (2022). Data quality evaluation in wearable monitoring.
        Scientific reports, 12(1), 21412.
    """
    nb_minimal = np.count_nonzero(signal < threshold)
    minimal_range = nb_minimal / len(signal)

    return np.round(minimal_range, 4)


def rac_sqi(signal, threshold, duration=2):
    """
    Compute the Rate of Amplitude Change (RAC) in the signal for windows
    of length defines by `duration`.

    Parameters
    ----------
    signal : vector
        Signal on which to compute the RAC.
    threshold : float
        Threshold to consider to evalutate the RAC. If the RAC is above
        that threshold, the quality of the signal in the given window is
        considered as bad.
    duration : float
        Duration of the windows on which to compute the RAC.
        Default to 2.

    Returns
    -------
    rac_ratio : float
        Ratio between the number of windows above `threshold` and the overall
        number of windows.

    References
    ----------
    Böttcher, S., Vieluf, S., Bruno, E., Joseph, B., Epitashvili, N., Biondi, A., ...
        & Loddenkemper, T. (2022). Data quality evaluation in wearable monitoring.
        Scientific reports, 12(1), 21412.
    """
    nb_windows = len(signal) // duration
    rac_values = []

    for i in range(nb_windows):
        window_start = i * duration
        window_end = window_start + duration

        window = signal[window_start:window_end]
        highest_value = max(window)
        lowest_value = min(window)

        rac = abs(highest_value - lowest_value) / min(highest_value, lowest_value)
        rac_values.append(rac)

    rac_ratio = np.count_nonzero(np.array(rac_values) > threshold) / len(signal)

    return np.round(rac_ratio, 4)


def rise_time_sqi(onsets, peaks, sampling_rate):
    """
    Return the delay between onset and peak of each SCR

    Parameters
    ----------
    Onsets :
        List of onsets indices
    Peaks :
        List of peaks indices
    sampling_rate:

    Returns
    -------
    """
    try:
        if len(onsets) != len(peaks):
            raise ValueError("ERROR: onsets and peaks have different length")
        else:
            rise_time = []
            for onset, peak in zip(onsets, peaks):
                rise_time.append((peak - onset) / sampling_rate)

    except Exception as error:
        print(repr(error))

    return rise_time


def threshold_sqi(metric, threshold, op=None):
    """
    Return quality assessment based on a threshold.

    Parameters
    ----------
    metric : int or float
        Metric to consider for the quality assessment.
    threshold : int, float or list
        Value to use to evaluate the quality of the `metric`.
        If a list of two elements is passed, the `metric` needs
        to be within that range to be consider as good.
    op : operator function
        Operator to use for the comparaison between the `metric` and the `threshold`.
        Only considered if `threshold` is a int or a float.
        See https://docs.python.org/2/library/operator.html#module-operator
        for all the possible operators.

    Returns
    -------
    a : str
        Quality Assessment given a metric and a threshold.

    Examples
    --------
    >>> threshold_sqi(6, 4, operator.gt)
    Acceptable
    >>> threshold_sqi(6, 4, operator.lt)
    Not acceptable
    >>> threshold_sqi(6, [2, 8])
    """
    if type(threshold) == list:
        if len(threshold) != 2:
            print(
                "Length of threshold should be 2 to set a range, otherwise pass an int or a float."
            )
        # Check if `metric` is within the range of values defined in `threshold`
        elif min(threshold) <= metric <= max(threshold):
            return "Acceptable"
        else:
            return "Not acceptable"
    else:
        if op(metric, threshold):
            return "Acceptable"
        else:
            return "Not acceptable"


# ==================================================================================
# Signals quality report
# ==================================================================================


def generate_report(summary, save, filename):
    """
    Generate quality assessment report in html format

    Parameters
    ----------
    summary : dict or list of dict
        Dictionnary contaning sqi values for a specified signal.
        List of dictationaries can be passed to include multiple
        signals to the report.
    save : str
        Directory to save the generated report.
    filename : str
        Name of the output file.

    Examples
    --------
    """
    # Generate the report in HTML format
    html_report = """
    <html>
    <head>
    <style>
        table {
        font-family: Arial, sans-serif;
        border-collapse: collapse;
        width: 100%;
        }

        td, th {
        border: 1px solid #dddddd;
        text-align: left;
        padding: 8px;
        }

        th {
        background-color: #dddddd;
        }
    </style>
    </head>
    <body>
    <h2>Signal Quality Report</h2>
    """

    for k in summary.keys():
        html_report += f"""
        <h3>{k} Signal</h3>
        """
        for metric in summary[k].keys():
            if metric == "Quality":
                html_report += f"""
                <br><b>{metric}</b> : {summary[k][metric]}
            """
            else:
                html_report += f"""
                    <br>{metric} : {summary[k][metric]}
            """

    # Complete the HTML report
    html_report += """
    </body>
    </html>
    """

    # Save the HTML report to a file
    print("Saving html report")
    with open(os.path.join(save, f"{filename}_quality.html"), "w") as file:
        file.write(html_report)
        file.close()


if __name__ == "__main__":
    neuromod_bio_sqi()
