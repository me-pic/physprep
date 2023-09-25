# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning

import operator

import numpy as np
from scipy.stats import kurtosis, skew

# ==================================================================================
# Signal Quality Indices
# ==================================================================================


def sqi_cardiac_overview(info):
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

    # Descriptive indices on overall signal
    summary["Ectopic"] = info["Ectopic"]
    summary["Missed"] = info["Missed"]
    summary["Extra"] = info["Extra"]
    summary["Long"] = info["Long"]
    summary["Short"] = info["Short"]

    return summary


def sqi_eda_overview(info, feature_quality="ScrOnsets", threshold=0):
    """
    Report SQI on the overall processed EDA signal.

    Parameters
    ----------
    info : dict
        Output from the process.py script.
    feature_quality : str
        Feature to consider to evaluate the quality if the signal.
    threshold : int or float
        Threshold to consider to evaluate the quality if the signal.

    Returns
    -------
    summary : dict
        Dictionary containing sqi values.
    """
    summary = {}
    summary["Quality"] = threshold_sqi(len(info[feature_quality]), threshold, operator.gt)

    return summary


def sqi_cardiac(
    signal_cardiac,
    info,
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
    mean_NN : list
        Range to consider to evaluate the quality of the signal based on the
        mean NN interval.
    std_NN : int or float
        Value to consider to evaluate the quality if the signal based on the
        std NN interval.
    window : list
        List of two elements specifying the window on which to compute the SQI.

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
    peaks = "Peaks"
    # Segment info according to the specify window
    if window is not None:
        start, end = min(window), max(window)
        sub_list = [x for x in info[peaks] if start <= x <= end]
        min_index = np.where(info[peaks] == min(sub_list))[0][0]
        max_index = np.where(info[peaks] == max(sub_list))[0][0] + 1
    else:
        min_index = 0
        max_index = len(info[peaks])

    # Descriptive indices on NN intervals
    summary["Mean_NN_intervals (ms)"] = np.round(
        np.mean(info["CleanRRSystole"][min_index:max_index]), 4
    )
    summary["Median_NN_intervals (ms)"] = np.round(
        np.median(info["CleanRRSystole"][min_index:max_index]),
        4,
    )
    summary["SD_NN_intervals (ms)"] = np.round(
        np.std(info["CleanRRSystole"][min_index:max_index], ddof=1),
        4,
    )
    summary["Min_NN_intervals (ms)"] = np.round(
        np.min(info["CleanRRSystole"][min_index:max_index]), 4
    )
    summary["Max_NN_intervals (ms)"] = np.round(
        np.max(info["CleanRRSystole"][min_index:max_index]), 4
    )
    # Descriptive indices on heart rate
    summary["Mean_HR (bpm)"] = metrics_hr_sqi(
        info["CleanRRSystole"][min_index:max_index],
        metric="mean",
    )
    summary["Median_HR (bpm)"] = metrics_hr_sqi(
        info["CleanRRSystole"][min_index:max_index],
        metric="median",
    )
    summary["SD_HR (bpm)"] = metrics_hr_sqi(
        info["CleanRRSystole"][min_index:max_index],
        metric="std",
    )
    summary["Min_HR (bpm)"] = metrics_hr_sqi(
        info["CleanRRSystole"][min_index:max_index],
        metric="min",
    )
    summary["Max_HR (bpm)"] = metrics_hr_sqi(
        info["CleanRRSystole"][min_index:max_index],
        metric="max",
    )
    # Descriptive indices on overall signal
    summary["Skewness"] = np.round(kurtosis(signal_cardiac), 4)
    summary["Kurtosis"] = np.round(skew(signal_cardiac), 4)

    # Quality assessment based on mean NN intervals and std
    if (
        threshold_sqi(
            np.mean(info["CleanRRSystole"][min_index:max_index]),
            mean_NN,
        )
        == "Acceptable"
        and threshold_sqi(
            np.std(
                info["CleanRRSystole"][min_index:max_index],
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


def sqi_eda(signal_eda, info, window=None):
    """
    Extract SQI for EDA processed signal.

    Parameters
    ----------
    signal_eda : DataFrame
        Output from the process.py script.
    info : dict
        Output from the process.py script.
    window : list
        List of two elements specifying the window on which to compute the SQI.

    Returns
    -------
    summary : DataFrame
        DataFrame containing sqi values.
    """
    summary = {}
    # Segment info according to the specify window
    if window is not None:
        try:
            start, end = min(window), max(window)
            sub_list = [x for x in info["ScrPeaks"] if start <= x <= end]
            min_index = info["ScrPeaks"].index(min(sub_list))
            max_index = info["ScrPeaks"].index(max(sub_list))
        except Exception:
            print("No SCR detected...")
            min_index = max_index = 0
    else:
        min_index = 0
        max_index = len(info["ScrPeaks"])

    # Descriptive indices on overall signal
    # summary["Minimal_range"] = minimal_range_sqi(
    #    signal_eda["EDA_Clean"], threshold=0.05
    # )
    # summary["RAC"] = rac_sqi(signal_eda["EDA_Clean"], threshold=0.2, duration=2)
    summary["Mean_EDA"] = np.round(np.mean(signal_eda["eda_clean"]), 4)
    summary["Median_EDA"] = np.round(np.median(signal_eda["eda_clean"]), 4)
    summary["SD_EDA"] = np.round(np.std(signal_eda["eda_clean"]), 4)
    summary["Min_EDA"] = np.round(np.min(signal_eda["eda_clean"]), 4)
    summary["Max_EDA"] = np.round(np.max(signal_eda["eda_clean"]), 4)
    # Descriptive indices on SCL
    summary["Mean_SCL"] = np.round(np.mean(signal_eda["eda_tonic"]), 4)
    summary["SD_SCL"] = np.round(np.std(signal_eda["eda_tonic"]), 4)
    summary["Median_SCL"] = np.round(np.median(signal_eda["eda_tonic"]), 4)
    summary["Min_SCL"] = np.round(np.min(signal_eda["eda_tonic"]), 4)
    summary["Max_SCL"] = np.round(np.max(signal_eda["eda_tonic"]), 4)
    # Descriptive indices on SCR
    summary["Mean_SCR"] = np.round(np.mean(signal_eda["eda_phasic"]), 4)
    summary["SD_SCR"] = np.round(np.std(signal_eda["eda_phasic"]), 4)
    summary["Median_SCR"] = np.round(np.median(signal_eda["eda_phasic"]), 4)
    summary["Min_SCR"] = np.round(np.min(signal_eda["eda_phasic"]), 4)
    summary["Max_SCR"] = np.round(np.max(signal_eda["eda_phasic"]), 4)
    # Descriptive indices on SCR
    if min_index != max_index and max_index != 0:
        summary["Number_of_detected_peaks"] = len(info["ScrPeaks"][min_index:max_index])
        # Descriptive indices on SCR rise time
        summary["Mean_rise_time"] = np.round(
            np.mean(info["ScrRisetime"][min_index:max_index]), 4
        )
        summary["SD_rise_time"] = np.round(
            np.std(info["ScrRisetime"][min_index:max_index]), 4
        )
        summary["Median_rise_time"] = np.round(
            np.median(info["ScrRisetime"][min_index:max_index]), 4
        )
        summary["Min_rise_time"] = np.round(
            np.min(info["ScrRisetime"][min_index:max_index]), 4
        )
        summary["Max_rise_time"] = np.round(
            np.max(info["ScrRisetime"][min_index:max_index]), 4
        )
        # Descriptive indices on SCR Recovery
        summary["Mean_recovery_time"] = np.round(
            np.mean(info["ScrRecoverytime"][min_index:max_index]), 4
        )
        summary["SD_recovery_time"] = np.round(
            np.std(info["ScrRecoverytime"][min_index:max_index]), 4
        )
        summary["Median_recovery_time"] = np.round(
            np.median(info["ScrRecoverytime"][min_index:max_index]), 4
        )
        summary["Min_recovery_time"] = np.round(
            np.min(info["ScrRecoverytime"][min_index:max_index]), 4
        )
        summary["Max_recovery_time"] = np.round(
            np.max(info["ScrRecoverytime"][min_index:max_index]), 4
        )
    else:
        summary["Number_of_detected_peaks"] = 0
        summary["Mean_rise_time"] = np.nan
        summary["SD_rise_time"] = np.nan
        summary["Median_rise_time"] = np.nan
        summary["Min_rise_time"] = np.nan
        summary["Max_rise_time"] = np.nan
        summary["Mean_recovery_time"] = np.nan
        summary["SD_recovery_time"] = np.nan
        summary["Median_recovery_time"] = np.nan
        summary["Min_recovery_time"] = np.nan
        summary["Max_recovery_time"] = np.nan

    return summary


def sqi_rsp(signal_rsp, mean_rate=0.5):
    """
    Extract SQI for respiratory processed signal.

    Parameters
    ----------
    signal_rsp : DataFrame
        Output from the process.py script.
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
    summary["Mean_Amp"] = np.round(np.mean(signal_rsp["rsp_amplitude"]), 4)
    summary["Median_Amp"] = np.round(np.median(signal_rsp["rsp_amplitude"]), 4)
    summary["SD_Amp"] = np.round(np.std(signal_rsp["rsp_amplitude"]), 4)
    summary["Min_Amp"] = np.round(np.min(signal_rsp["rsp_amplitude"]), 4)
    summary["CV_Amp"] = np.round(np.max(signal_rsp["rsp_amplitude"]), 4)
    summary["variability_Amp"] = np.round(
        np.std(signal_rsp["rsp_amplitude"]) / np.mean(signal_rsp["rsp_amplitude"]),
        4,
    )
    # Descriptive indices on signal rate
    summary["Mean_Rate"] = np.round(np.mean(signal_rsp["rsp_rate"]), 4)
    summary["Median_Rate"] = np.round(np.median(signal_rsp["rsp_rate"]), 4)
    summary["SD_Rate"] = np.round(np.std(signal_rsp["rsp_rate"]), 4)
    summary["Min_Rate"] = np.round(np.min(signal_rsp["rsp_rate"]), 4)
    summary["Max_Rate"] = np.round(np.max(signal_rsp["rsp_rate"]), 4)
    summary["CV_Rate"] = np.round(
        np.std(signal_rsp["rsp_rate"]) / np.mean(signal_rsp["rsp_rate"]), 4
    )
    # Quality assessment based on the mean respiratory rate
    if (
        threshold_sqi(np.mean(signal_rsp["rsp_rate"]) / 60, mean_rate, operator.lt)
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
    else:
        metric_rr = np.nan
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
