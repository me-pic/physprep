# -*- coding: utf-8 -*-
# !/usr/bin/env python
"""
Neuromod cleaning utilities.
"""
import neurokit2 as nk
import numpy as np
from scipy import signal

from physprep import utils


def preprocessing_workflow(data, metadata, workflow_strategy):
    """
    Apply the preprocessing workflow to a physiological recordings.

    Parameters
    ----------
    data : dataframe
        The raw physiological signal.
    metadata : dict
        The metadata associated with the physiological recording.
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
    clean_signals = {}
    metadata_derivatives = {}
    strategies = {}
    col = []
    sf = []

    # Remove padding in data if any
    if metadata["StartTime"] > 2e-3:
        data = remove_padding(data, trigger_threshold=4)

    # Iterate over content of `workflow_strategy`
    for signal_type in workflow_strategy:
        print(f"Preprocessing {signal_type}...\n")
        if signal_type not in ["trigger", "concurrentWith"]:
            if "preprocessing_strategy" in workflow_strategy[
                signal_type
            ] and workflow_strategy[signal_type]["preprocessing_strategy"] not in [
                "",
                " ",
                None,
            ]:
                # Get the timeseries
                if signal_type in data.columns:
                    raw = data[signal_type]
                elif workflow_strategy[signal_type]["id"] in data.columns:
                    raw = data[workflow_strategy[signal_type]["id"]]
                else:
                    raise ValueError(f"Signal type {signal_type} not found in the data.")
                # Apply the preprocessing strategy
                raw, clean, components, sampling_rate, preproc_strategy = preprocess_signal(
                    raw,
                    workflow_strategy[signal_type]["preprocessing_strategy"],
                    sampling_rate=metadata["SamplingFrequency"],
                )

                sf = [*sf, sampling_rate, sampling_rate]

                clean_signals.update(
                    {f'{signal_type}_raw': raw}
                )
                clean_signals.update(
                    {f'{signal_type}_clean': clean}
                )

                col = [*col, f'{signal_type}_raw', f'{signal_type}_clean']

                if components is not None:
                    clean_signals.update(
                        {f'{signal_type}_tonic': np.array(components['EDA_Tonic'])}
                    )
                    clean_signals.update(
                        {f'{signal_type}_phasic': np.array(components['EDA_Phasic'])}
                    )

                    col = [*col, 'electrodermal_tonic', 'electrodermal_phasic']
                
                strategies.update(
                    {f'{signal_type}': preproc_strategy}
                )
            else:
                print(
                    f"No preprocessing strategy specified for {signal_type}. "
                    "The preprocessing step will be skipped."
                )

    # Checking if there are different sampling frequencies
    if  all(f == sf[0] for f in sf):
        sf = sf[0]

    metadata_derivatives.update(
        {
            "StartTime": data["time"].loc[0],
            "SamplingFrequency": sf,
            "Columns": col,
        }
    )

    return clean_signals, metadata_derivatives, strategies


def remove_padding(data, start_time=None, end_time=None, trigger_threshold=None):
    """
    Remove padding in data if any.

    Parameters
    ----------
    data : dataframe
        The raw physiological timeseries.
    start_time : float
        The time at which the recording started (in seconds). If specified, data will
        be cropped before this time. Default to None.
    end_time : float
        The time at which the recording ended (in seconds). If specified, data will
        be cropped after this time. Default to None.
    trigger_threshold : float
        The threshold to use to detect the trigger. The data before the first
        detected trigger and after the last detected trigger will be removed.
        If `start_time` and/or `end_time` are specified, the trigger detection
        will be skipped and the data will be cropped based on the specified
        times. Default to None.

    Returns
    -------
    data : array or Series
        The raw physiological signal without padding.
    """
    if trigger_threshold is not None:
        # Cropped the data based on the detected triggers
        trigger = data[data["TTL"] > trigger_threshold]
        start_time = list(trigger.index)[0]
        end_time = list(trigger.index)[-1]
    else:
        # Crop the data based on the specified start and end times
        if start_time is not None:
            start_time = list(data[data["time"] == start_time].index)[0]
        if end_time is not None:
            end_time = list(data[data["time"] == end_time].index)[0]

    data = data.loc[start_time:end_time].reset_index(drop=True)

    return data


def preprocess_signal(signal, preprocessing_strategy, sampling_rate=1000):
    """
    Apply a preprocessing workflow to a physiological signal.

    Parameters
    ----------
    signal : array or Series
        The raw physiological signal.
    preprocessing_strategy : dict
        Dictionary containing the preprocessing steps to apply to the signal.
    sampling_rate : float
        The sampling frequency of `signal` (in Hz, i.e., samples/second).

    Returns
    -------
    raw : array or Series
        The raw physiological signal.
    signal : array or Series
        The cleaned physiological signal.
    """
    raw = signal
    components = None
    # Retrieve preprocessing steps
    preprocessing = utils.get_config(preprocessing_strategy, strategy="preprocessing")
    # Iterate over preprocessing steps as defined in the configuration file
    for step in preprocessing:
        print(f"...Applying {step['step']}\n")
        if step["step"] == "filtering":
            if step["parameters"]["method"] == "notch":
                signal = comb_band_stop(signal, sampling_rate, step["parameters"])
            elif step["parameters"]["method"] == "median":
                signal = median_filter(signal, step["parameters"], sampling_rate)
            else:
                signal = nk.signal_filter(
                    signal, sampling_rate=sampling_rate, **step["parameters"]
                )
        elif step["step"] == "signal_resample":
            # Resample the signal
            signal = nk.signal_resample(
                signal, sampling_rate=sampling_rate, **step["parameters"]
            )
            raw = nk.signal_resample(
                raw, sampling_rate=sampling_rate, **step["parameters"]
            )
            sampling_rate = step["parameters"]["desired_sampling_rate"]
        elif step["step"] == "phasic":
            # This step is only to extract phasic and tonic components from EDA signal
            components = nk.eda_phasic(signal, sampling_rate=sampling_rate, **step["parameters"])
        else:
            raise ValueError(
                f"Unknown preprocessing step: {step['step']}. Make sure the "
                "preprocessing strategy is properly defined. For more "
                "details, please refer to the Physprep documentation."
            )
        print(f"   step: {step['step']}, parameters: {step['parameters']} done !\n")

    return raw, signal, components, sampling_rate, preprocessing


def comb_band_stop(data, sampling_rate, params):
    """
    Series of notch filters aligned with the scanner gradient's harmonics.

    Parameters
    ----------
    data : array
        The signal to be filtered.
    sampling_rate : float
        The sampling frequency of `signal` (in Hz, i.e., samples/second).
    params : dict
        The parameters of the scanner sequence (i.e. `tr`, `slices` and `mb` if
        `notch_method` is 'bottenhorn').

    Returns
    -------
    filtered : array
        The filtered signal.

    References
    ----------
    Biopac Systems, Inc. Application Notes: application note 242
        ECG Signal Processing During fMRI
        https://www.biopac.com/wp-content/uploads/app242x.pdf
    Bottenhorn, K. L., Salo, T., Riedel, M. C., Sutherland, M. T.,
        Robinson, J. L., Musser, E. D., & Laird, A. R. (2021).
        Denoising physiological data collected during multi-band,
        multi-echo EPI sequences. bioRxiv, 2021-04.
        https://doi.org/10.1101/2021.04.01.437293

    See also
    --------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.iirnotch.html
    """
    # Setting scanner sequence parameters
    nyquist = np.float64(sampling_rate / 2)
    tr = params["tr"]
    slices = params["slices"]
    Q = params["Q"]
    if params["notch_method"] == "biopac":
        notches = {"slices": slices / tr, "tr": 1 / tr}
    elif params["notch_method"] == "bottenhorn":
        mb = params["mb"]
        notches = {"slices": slices / mb / tr, "tr": 1 / tr}
    # band stopping each frequency specified with notches dict
    for notch in notches:
        for i in np.arange(1, int(nyquist / notches[notch])):
            f0 = notches[notch] * i
            w0 = f0 / nyquist
            b, a = signal.iirnotch(w0, Q)
            data = signal.filtfilt(b, a, data)
    return data


def median_filter(data, params, sampling_rate):
    """
    Series of notch filters aligned with the scanner gradient's harmonics.

    Parameters
    ----------
    data : array
        The signal to be filtered.
    params : dict
        The parameters for the median filtering (i.e. `window_size`).
    sampling_rate : float
        The sampling frequency of `signal` (in Hz, i.e., samples/second).

    Returns
    -------
    filtered : array
        The filtered signal.

    References
    ----------
    Privratsky, A. A., Bush, K. A., Bach, D. R., Hahn, E. M., & Cisler, J. M. (2020). 
        Filtering and model-based analysis independently improve skin-conductance response 
        measures in the fMRI environment: Validation in a sample of women with PTSD. 
        International Journal of Psychophysiology, 158, 86-95.
        https://doi.org/10.1016/j.ijpsycho.2020.09.015

    See also
    --------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.medfilt.html
    """
    window_size = int(params['window_size']*sampling_rate)
    # window_size must be odd
    if window_size % 2 == 0:
        window_size += 1
        
    return signal.medfilt(data, kernel_size=window_size)