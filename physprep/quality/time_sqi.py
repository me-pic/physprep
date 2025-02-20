# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning

import operator

import numpy as np
from scipy.stats import kurtosis, skew
from biosppy.quality import eda_sqi_bottcher

# ==================================================================================
# Signal Quality Indices
# ==================================================================================

MAPPING_METRICS = {
    "Mean": np.mean,
    "SD": np.std,
    "Median": np.median,
    "Min": np.min,
    "Max": np.max,
}


def sqi_cardiac_overview(info, sampling_rate):
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
    # Descriptive indices on overall signal
    peaks = [p for p in info.keys() if 'peak_corrected' in p][0]
    try:
        from systole.utils import input_conversion
        from systole.correction import correct_rr
        rr = input_conversion(
            info[peaks], 
            'peaks_idx', 
            output_type='rr_ms',
            sfreq = sampling_rate)
        
        _, (missed, extra, ectopic, short, long) = correct_rr(rr)

        summary = {
            "Ectopic": ectopic,
            "Missed": missed,
            "Extra": extra,
            "Long": long,
            "Short": short
        }

    except ImportError as e:
        print('Systole not imported... Can not run sqi_cardiac_overview function')
        summary = {}

    return summary


def sqi_eda_overview(feature_quality, threshold=0):
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
    summary["Quality"] = threshold_sqi(feature_quality, threshold, operator.gt)

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

    peaks = [p for p in info.keys() if 'peak_corrected' in p][0]
 
    rr = input_conversion(info[peaks], input_type="peaks_idx", output_type="rr_ms")
    corrected_rr, _ = correct_rr(rr)

    
    # Segment info according to the specify window
    if window is not None:
        start, end = min(window), max(window)
        sub_list = [x for x in info[peaks] if start <= x <= end]
        if len(sub_list) > 0:
            min_index = np.where(info[peaks] == min(sub_list))[0][0]
            max_index = np.where(info[peaks] == max(sub_list))[0][0] + 1
        else:
            min_index = None
    else:
        min_index = 0
        max_index = len(info[peaks])

    # Descriptive indices on overall signal
    summary["Skewness"] = np.round(kurtosis(signal_cardiac), 4)
    summary["Kurtosis"] = np.round(skew(signal_cardiac), 4)

    if min_index is not None:
        # Descriptive indices on NN intervals
        for metric in MAPPING_METRICS:
            summary[f"{metric}_NN_intervals (ms)"] = np.round(
                MAPPING_METRICS[metric](corrected_rr[min_index:max_index]), 4
            )

        # Descriptive indices on heart rate
        for metric in MAPPING_METRICS:
            summary[f"{metric}_HR (bpm)"] = metrics_hr_sqi(
                corrected_rr[min_index:max_index],
                metric=metric.lower(),
            )

        # Quality assessment based on mean NN intervals and std
        if (
            threshold_sqi(
                np.mean(corrected_rr[min_index:max_index]),
                mean_NN,
            )
            == "Acceptable"
            and threshold_sqi(
                np.std(
                    corrected_rr[min_index:max_index],
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
    else:
        # Descriptive indices on NN intervals
        for metric in MAPPING_METRICS:
            summary[f"{metric}_NN_intervals (ms)"] = None

        # Descriptive indices on heart rate
        for metric in MAPPING_METRICS:
            summary[f"{metric}_HR (bpm)"] = None
        summary['Quality'] = "Not acceptable"

    return summary


def sqi_eda(signal_eda, signal_phasic, signal_tonic, info, sampling_rate=10000, window=None):
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
            sub_list = [x for x in info["scr_peak"] if start <= x <= end]
            min_index = info["scr_peak"].index(min(sub_list))
            max_index = info["scr_peak"].index(max(sub_list))
        except Exception:
            print("No SCR detected...")
            min_index = max_index = 0
    else:
        min_index = 0
        max_index = len(info["scr_peak"])

    # Description indices on EDA
    for metric in MAPPING_METRICS:
        summary[f"{metric}_EDA"] = np.round(
            MAPPING_METRICS[metric](signal_eda), 4
        )
    # Descriptive indices on SCL
    if signal_tonic is not None:
        for metric in MAPPING_METRICS:
            summary[f"{metric}_SCL"] = np.round(
                MAPPING_METRICS[metric](signal_tonic), 4
            )
    # Descriptive indices on SCR
    if signal_phasic is not None:
        for metric in MAPPING_METRICS:
            summary[f"{metric}_SCR"] = np.round(
                MAPPING_METRICS[metric](signal_phasic), 4
            )
    # Descriptive indices on SCR
    if min_index != max_index and max_index != 0:
        summary["Number_of_detected_peaks"] = len(info["scr_peak"][min_index:max_index])
    
    rac, qa_bottcher = rac_sqi(signal_eda, sampling_rate)
    if qa_bottcher == 1:
        summary['Quality'] = "Acceptable"
    elif qa_bottcher > 0.5:
        summary['Quality'] = "Mostly acceptable"
    else:
        summary['Quality'] = "Not acceptable"

    return summary


def sqi_rsp(signal_rsp, signal_amplitude, signal_rate, mean_rate=0.5):
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
    for metric in MAPPING_METRICS:
            summary[f"{metric}_Amp"] = np.round(
                MAPPING_METRICS[metric](signal_amplitude), 4
            )
    summary["variability_Amp"] = np.round(
        np.std(signal_amplitude) / np.mean(signal_amplitude),
        4,
    )
    # Descriptive indices on signal rate
    for metric in MAPPING_METRICS:
            summary[f"{metric}_Rate"] = np.round(
                MAPPING_METRICS[metric](signal_rate), 4
            )
    summary["variability_Rate"] = np.round(
        np.std(signal_rate) / np.mean(signal_rate),
        4,
    )

    # Quality assessment based on the mean respiratory rate
    if (
        threshold_sqi(np.mean(signal_rate) / 60, mean_rate, operator.lt)
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
    elif metric == "sd":
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


def rac_sqi(signal, sampling_rate, threshold=0.2, duration=2):
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
    nb_windows = len(signal) // (duration * sampling_rate)
    rac_values, qa_scores = [], []

    for i in range(nb_windows):
        window_start = i * (duration*sampling_rate)
        window_end = window_start + (duration*sampling_rate)

        window = signal[window_start:window_end]
        highest_value, highest_idx = max(window), np.argmax(window)
        lowest_value, lowest_idx = min(window), np.argmin(window)

        if lowest_idx<highest_idx:
            rac = abs(highest_value - lowest_value) / lowest_value
        elif lowest_idx>=highest_idx:
            rac = abs(highest_value - lowest_value) / highest_value
        rac_values.append(rac)
        qa_scores.append((rac < threshold) & (np.mean(window)>0.05))

    rac_ratio = np.count_nonzero(np.array(rac_values) < threshold) / nb_windows
    qa_ratio = np.count_nonzero(np.array(qa_scores)) / nb_windows

    return np.round(rac_ratio, 2), np.round(qa_ratio, 2)


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
