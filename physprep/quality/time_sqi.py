# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning

import operator
import traceback

import numpy as np
from scipy.stats import kurtosis, skew

# ==================================================================================
# Signal Quality Indices
# ==================================================================================


def sqi_cardiac_overview(info, data_type="ECG"):
    """
    Report artefacts identified by the library systole during
    the processing of the cardiac signals

    Parameters
    ----------
    info : dict
        Output from the process.py script.
    data_type : str
        Type of the signal. Valid options include 'ECG' and 'PPG'.
        Default to 'ECG'.

    Returns
    -------
    summary : Dict
        Dictionary containing sqi values.

    Reference
    ---------
    Legrand et al., (2022). Systole: A python package for cardiac signal
        synchrony and analysis. Journal of Open Source Software, 7(69), 3832,
        https://doi.org/10.21105/joss.03832
    """
    summary = {}
    # Make sure the data_type string is in capital letters
    data_type = data_type.upper()

    # Descriptive indices on overall signal
    summary["Ectopic"] = info[f"{data_type}_ectopic"]
    summary["Missed"] = info[f"{data_type}_missed"]
    summary["Extra"] = info[f"{data_type}_extra"]
    summary["Long"] = info[f"{data_type}_long"]
    summary["Short"] = info[f"{data_type}_short"]
    summary["Cumulseconds_rejected"] = info[f"{data_type}_cumulseconds_rejected"]
    summary["%_rejected_segments"] = np.round(info[f"{data_type}_%_rejected_segments"], 4)

    return summary


def sqi_cardiac(
    signal_cardiac,
    info,
    data_type="ECG",
    sampling_rate=10000,
    mean_NN=[600, 1200],
    std_NN=300,
    window=None,
):
    """
    Extract SQI for ECG/PPG processed signal.

    Parameters
    ----------
    signal_cardiac : DataFrame
        Output from the process.py script.
    info : dict
        Output from the process.py script.
    data_type : str
        Type of the signal. Valid options include 'ECG' and 'PPG'.
        Default to 'ECG'.
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
    summary : Dict
        Dictionary containing sqi values.

    Reference
    ---------
    Le, V. K. D., Ho, H. B., Karolcik, S., Hernandez, B., Greeff, H.,
        Nguyen, V. H., ... & Clifton, D. (2022). vital_sqi: A Python
        package for physiological signal quality control.
        Frontiers in Physiology, 2248.
    """
    summary = {}
    # Make sure the data_type string is in capital letters
    data_type = data_type.upper()
    if data_type == "PPG":
        peaks = "PPG_Peaks"
    elif data_type == "ECG":
        peaks = "ECG_R_Peaks"
    # Segment info according to the specify window
    if window is not None:
        sub_list = [x for x in info[peaks] if min(window) <= x <= max(window)]
        min_index = info[peaks].index(min(sub_list))
        max_index = info[peaks].index(max(sub_list)) + 1
    else:
        min_index = 0
        max_index = len(info[peaks])

    # Descriptive indices on NN intervals
    summary["Mean_NN_intervals (ms)"] = np.round(
        np.mean(info[f"{data_type}_clean_rr_systole"][min_index:max_index]), 4
    )
    summary["Median_NN_intervals (ms)"] = np.round(
        np.median(info[f"{data_type}_clean_rr_systole"][min_index:max_index]),
        4,
    )
    summary["SD_NN_intervals (ms)"] = np.round(
        np.std(info[f"{data_type}_clean_rr_systole"][min_index:max_index], ddof=1),
        4,
    )
    summary["Min_NN_intervals (ms)"] = np.round(
        np.min(info[f"{data_type}_clean_rr_systole"][min_index:max_index]), 4
    )
    summary["Max_NN_intervals (ms)"] = np.round(
        np.max(info[f"{data_type}_clean_rr_systole"][min_index:max_index]), 4
    )
    # Descriptive indices on heart rate
    summary["Mean_HR (bpm)"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"][min_index:max_index],
        metric="mean",
    )
    summary["Median_HR (bpm)"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"][min_index:max_index],
        metric="median",
    )
    summary["SD_HR (bpm)"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"][min_index:max_index],
        metric="std",
    )
    summary["Min_HR (bpm)"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"][min_index:max_index],
        metric="min",
    )
    summary["Max_HR (bpm)"] = metrics_hr_sqi(
        info[f"{data_type}_clean_rr_systole"][min_index:max_index],
        metric="max",
    )
    # Descriptive indices on overall signal
    summary["Skewness"] = np.round(kurtosis(signal_cardiac), 4)
    summary["Kurtosis"] = np.round(skew(signal_cardiac), 4)

    # Quality assessment based on mean NN intervals and std
    if (
        threshold_sqi(
            np.mean(info[f"{data_type}_clean_rr_systole"][min_index:max_index]),
            mean_NN,
        )
        == "Acceptable"
        and threshold_sqi(
            np.std(
                info[f"{data_type}_clean_rr_systole"][min_index:max_index],
                ddof=1,
            ),
            std_NN,
            operator.lt,
        )
        == "Acceptable"
    ):
        summary["Quality"] = "Acceptable"
    else:
        summary["Quality"] = "Not acceptable"

    return summary


def sqi_eda(signal_eda, info, sampling_rate=10000, window=None):
    """
    Extract SQI for EDA processed signal.

    NOTE: add option to compute the latency between stimulus onset
    and SCR onset

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
    # Segment info according to the specify window
    if window is not None:
        sub_list = [x for x in info["SCR_Peaks"] if min(window) <= x <= max(window)]
        min_index = info["SCR_Peaks"].index(min(sub_list))
        max_index = info["SCR_Peaks"].index(max(sub_list))
    else:
        min_index = 0
        max_index = len(info["SCR_Peaks"])

    # Descriptive indices on overall signal
    # summary["Minimal_range"] = minimal_range_sqi(
    #    signal_eda["EDA_Clean"], threshold=0.05
    # )
    # summary["RAC"] = rac_sqi(signal_eda["EDA_Clean"], threshold=0.2, duration=2)
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
    summary["Max_SCL"] = np.round(np.max(signal_eda["EDA_Tonic"]), 4)
    # Descriptive indices on SCR
    summary["Mean_SCR"] = np.round(np.mean(signal_eda["EDA_Phasic"]), 4)
    summary["SD_SCR"] = np.round(np.std(signal_eda["EDA_Phasic"]), 4)
    summary["Median_SCR"] = np.round(np.median(signal_eda["EDA_Phasic"]), 4)
    summary["Min_SCR"] = np.round(np.min(signal_eda["EDA_Phasic"]), 4)
    summary["Max_SCR"] = np.round(np.max(signal_eda["EDA_Phasic"]), 4)
    summary["Number_of_detected_peaks"] = len(info["SCR_Peaks"][min_index:max_index])
    # Descriptive indices on SCR rise time
    summary["Mean_rise_time"] = np.round(
        np.mean(info["SCR_RiseTime"][min_index:max_index]), 4
    )
    summary["SD_rise_time"] = np.round(
        np.std(info["SCR_RiseTime"][min_index:max_index]), 4
    )
    summary["Median_rise_time"] = np.round(
        np.median(info["SCR_RiseTime"][min_index:max_index]), 4
    )
    try:
        summary["Min_rise_time"] = np.round(
            np.min(info["SCR_RiseTime"][min_index:max_index]), 4
        )
    except Exception:
        traceback.print_exc()
        summary["Min_rise_time"] = np.nan
    try:
        summary["Max_rise_time"] = np.round(
            np.max(info["SCR_RiseTime"][min_index:max_index]), 4
        )
    except Exception:
        traceback.print_exc()
        summary["Max_rise_time"] = np.nan
    # Descriptive indices on SCR Recovery
    summary["Mean_recovery_time"] = np.round(
        np.mean(info["SCR_Recovery"][min_index:max_index]), 4
    )
    summary["SD_recovery_time"] = np.round(
        np.std(info["SCR_Recovery"][min_index:max_index]), 4
    )
    summary["Median_recovery_time"] = np.round(
        np.median(info["SCR_Recovery"][min_index:max_index]), 4
    )
    summary["Min_recovery_time"] = np.round(
        np.min(info["SCR_Recovery"][min_index:max_index]), 4
    )
    summary["Max_recovery_time"] = np.round(
        np.max(info["SCR_Recovery"][min_index:max_index]), 4
    )

    return summary


def sqi_rsp(signal_rsp, sampling_rate=10000, mean_rate=0.5):
    """
    Extract SQI for respiratory processed signal.

    Parameters
    ----------
    signal_rsp : DataFrame
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
        np.std(signal_rsp["RSP_Amplitude"]) / np.mean(signal_rsp["RSP_Amplitude"]),
        4,
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
    Compute the mean heart rate from the RR intervals.

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
    except Exception:
        metric_rr = np.nan
        traceback.print_exc()
        raise ValueError(f"Invalid metric: {metric}.")
    
    return metric_rr


def minimal_range_sqi(signal, threshold):
    """
    Compute the ratio between the number of timepoints under a defined
    `threshold` and the signal length.

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
        Threshold to consider to evaluate the RAC. If the RAC is above
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
        Operator to use for the comparison between the `metric` and the `threshold`.
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
    Acceptable
    """
    if isinstance(threshold, list):
        if len(threshold) != 2:
            print(
                "Length of threshold should be 2 to set a range, "
                "otherwise pass an int or a float."
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
