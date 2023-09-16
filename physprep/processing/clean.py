# -*- coding: utf-8 -*-
# !/usr/bin/env python
"""
Neuromod cleaning utilities.
"""

# from pathlib import Path

import neurokit2 as nk
import numpy as np
from scipy import signal
import pandas as pd
from physprep import utils


def preprocessing_workflow(data, metadata, workflow_strategy, save=True):
    """
    Parameters
    ----------
    data : dataframe
        The raw physiological signal.
    metadata : dict
        The metadata associated with the physiological recording.
    workflow_strategy : dict
        Dictionary containing the content of the workflow strategy.
    save : bool
        Specify if the preprocessed signals should be saved.
        Default to True.
    """
    clean_signals = {}

    # Remove padding in data if any
    if metadata["StartTime"] > 2e-3:
        data = remove_padding(data)

    # Iterate over content of `workflow_strategy`
    for signal_type in workflow_strategy:
        if signal_type != "trigger":
            if "preprocessing_strategy" in signal and signal[
                "preprocessing_strategy"
            ] not in ["", " ", None]:
                clean = preprocess_signal(
                    data[signal_type.id],
                    signal_type["preprocessing_strategy"],
                    signal_type.id,
                    sampling_rate=metadata["SamplingFrequency"],
                )
                clean_signals.update({signal_type.id: clean})
            else:
                print(f"No preprocessing strategy specified for {signal_type}. The preprocessing step will be skipped.")
    

    if save:
        # Save preprocessed signal
        pass


def remove_padding(data, start_time=None, end_time=None, trigger_threshold=4):
    """
    Remove padding in data if any.

    Parameters
    ----------
    data : array or Series
        The raw physiological signal.
    start_time : float
        The time at which the recording started (in seconds). 
        Default to 0.
    end_time : float
        The time at which the recording ended (in seconds).
        Default to None.
    trigger_threshold : float
        The threshold to use to detect the trigger. The data before the first
        detected trigger and after the last detected trigger will be removed.
        If `start_time` and/or `end_time` are specified, the trigger detection
        will be skipped and the data will be cropped based on the specified
        times. Default to 4.

    Returns
    -------
    data : array or Series
        The raw physiological signal without padding.
    """
    if start_time is None and end_time is None:
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
    signal_type : str
        The type of physiological signal. Can be 'ECG', 'PPG', 'RSP', 'EDA'.
    sampling_rate : float
        The sampling frequency of `signal` (in Hz, i.e., samples/second).

    Returns
    -------
    signal : array or Series
        The cleaned physiological signal.
    """
    # Retrieve preprocessing steps
    preprocessing = utils.get_config(preprocessing_strategy, strategy="preprocessing")
    # Iterate over preprocessing steps as defined in the configuration file
    for step in preprocessing:
        if step["step"] == "filtering":
            if step["parameters"]["method"] == "notch":
                pass
            else:
                signal = nk.signal_filter(
                    signal, sampling_rate=sampling_rate, **step["parameters"]
                )
        elif step["step"] == "signal_resample":
            # Resample the signal
            signal = nk.signal_resample(
                signal, sampling_rate=sampling_rate, **step["parameters"]
            )
        else:
            raise ValueError(
                f"Unknown preprocessing step: {step['step']}. Make sure the "
                "preprocessing strategy is properly defined. For more "
                "details, please refer to the Physprep documentation."
            )

    # Create a dataframe with the raw and cleaned signal
    # Dataframe containing the cleaned physiological signal.
    # timeseries = pd.DataFrame(
    #    {f"{signal_type}_Raw": raw, f"{signal_type}_Clean": signal}
    #)

    return signal


















def neuromod_ppg_clean(ppg_signal, sampling_rate=10000.0, downsampling=None):
    """
    Clean a PPG signal.

    Prepare raw PPG signal for peak detection with specified method.

    Parameters
    ----------
    ppg_signal : list, array or Series
        The raw PPG channel.
    sampling_rate : float
        The sampling frequency of `ppg_signal` (in Hz, i.e., samples/second).
        Default to 10000.
    downsampling : int
        The desired sampling frequency (Hz). If None, the signal is not
        resample.
        Default to None.

    Returns
    -------
    ppg_clean : array
        Vector containing the cleaned PPG signal.
    """
    # Apply band pass filter
    ppg_clean = nk.signal_filter(
        ppg_signal, sampling_rate=sampling_rate, lowcut=0.5, highcut=8, order=3
    )
    # Downsample the signal if specified
    if downsampling is not None:
        ppg_clean = nk.signal_resample(
            ppg_clean,
            sampling_rate=sampling_rate,
            desired_sampling_rate=downsampling,
        )
        ppg_signal = nk.signal_resample(
            ppg_signal,
            sampling_rate=sampling_rate,
            desired_sampling_rate=downsampling,
        )

    return ppg_signal, ppg_clean


# ======================================================================
# Electrocardiogram (ECG)
# =======================================================================


def neuromod_ecg_clean(
    ecg_signal,
    sampling_rate=10000.0,
    method="biopac",
    me=False,
    downsampling=None,
):
    """
    Clean an ECG signal.

    Prepare a raw ECG signal for R-peak detection with the specified method.

    Parameters
    ----------
    ecg_signal : list, array or Series
        The raw ECG channel.
    sampling_rate : float
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Default to 10000.
    method : str
        The processing pipeline to apply between 'biopac' and 'bottenhorn'.
        Default to 'biopac'.
    me : bool
        Specify if the MRI sequence used was the multi-echo (True)
        or the single-echo (False).
        Default to False.
    downsampling : int
        The desired sampling frequency (Hz). If None, the signal is not
        resample.
        Default to None.

    Returns
    -------
    ecg_clean : array
        Vector containing the cleaned ECG signal.
    """
    if me:
        tr = 2.65
        mb = 2
        slices = 70
    else:
        tr = 1.49
        mb = 4
        slices = 60

    if method in ["biopac"]:
        ecg_clean = _ecg_clean_biopac(ecg_signal, sampling_rate, tr=tr, slices=slices)
    if method in ["bottenhorn", "bottenhorn2022"]:
        # Apply comb band pass filter with Bottenhorn correction
        print("... Applying the corrected comb band pass filter.")
        ecg_clean = _ecg_clean_bottenhorn(
            ecg_signal,
            sampling_rate=sampling_rate,
            tr=tr,
            mb=mb,
            slices=slices,
        )
    # Downsample the signal if specified
    if downsampling is not None:
        ecg_clean = nk.signal_resample(
            ecg_clean,
            sampling_rate=sampling_rate,
            desired_sampling_rate=downsampling,
        )
        ecg_signal = nk.signal_resample(
            ecg_signal,
            sampling_rate=sampling_rate,
            desired_sampling_rate=downsampling,
        )

    return ecg_signal, ecg_clean


# =============================================================================
# ECG internal : biopac recommendations
# =============================================================================
def _ecg_clean_biopac(ecg_signal, sampling_rate=10000.0, tr=1.49, slices=60, Q=100):
    """
    Single-band sequence gradient noise reduction.

    This function is a reverse-engineered appropriation of BIOPAC's
    application note 242. It only applies to signals polluted by single-band
    (f)MRI sequence.

    Parameters
    ----------
    ecg_signal : array
        The ECG channel.
    sampling_rate: float
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Default to 10000.
    tr : int
        The time Repetition of the MRI scanner.
        Default to 1.49.
    slices :
        The number of volumes acquired in the tr period.
        Default to 60.
    Q : int
        The filter quality factor.
        Default to 100.

    Returns
    -------
    ecg_clean : array
        The cleaned ECG signal.

    References
    ----------
    Biopac Systems, Inc. Application Notes: application note 242
        ECG Signal Processing During fMRI
        https://www.biopac.com/wp-content/uploads/app242x.pdf
    """
    # Setting scanner sequence parameters
    nyquist = np.float64(sampling_rate / 2)
    notches = {"slices": slices / tr, "tr": 1 / tr}
    # remove baseline wandering
    ecg_clean = nk.signal_filter(
        ecg_signal,
        sampling_rate=int(sampling_rate),
        lowcut=2,
    )
    # Filtering at specific harmonics
    ecg_clean = _comb_band_stop(notches, nyquist, ecg_clean, Q)
    # bandpass filtering
    ecg_clean = nk.signal_filter(
        ecg_clean,
        sampling_rate=sampling_rate,
        lowcut=2,
        highcut=20,
        method="butter",
        order=5,
    )

    return ecg_clean


def _ecg_clean_bottenhorn(
    ecg_signal,
    sampling_rate=10000.0,
    tr=1.49,
    mb=4,
    slices=60,
    Q=100,
    comb=False,
):
    """
    Multiband sequence gradient noise reduction.

    Parameters
    ----------
    ecg_signal : array
        The ECG channel.
    sampling_rate : float
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Default to 10000.
    tr : float
        The time Repetition of the MRI scanner.
        Default to 1.49.
    mb : 4
        The multiband acceleration factor.
        Default to 4.
    slices : int
        The number of volumes acquired in the tr period.
        Default to 60.
    Q : int
        The filter quality factor.
        Default to 100.

    Returns
    -------
    ecg_clean : array
        The cleaned ECG signal.

    References
    ----------
    Bottenhorn, K. L., Salo, T., Riedel, M. C., Sutherland, M. T.,
        Robinson, J. L., Musser, E. D., & Laird, A. R. (2021).
        Denoising physiological data collected during multi-band,
        multi-echo EPI sequences. bioRxiv, 2021-04.
        https://doi.org/10.1101/2021.04.01.437293

    See also
    --------
    https://neuropsychology.github.io/NeuroKit/_modules/neurokit2/signal/signal_filter.html#signal_filter
    """
    # Setting scanner sequence parameters
    nyquist = np.float64(sampling_rate / 2)
    notches = {"slices": slices / mb / tr, "tr": 1 / tr}

    # Remove low frequency artefacts: respiration & baseline wander using high
    # pass butterworth filter (order=2)
    print("... Applying high pass filter.")
    ecg_clean = nk.signal_filter(
        ecg_signal, sampling_rate=sampling_rate, lowcut=2, method="butter"
    )
    # Filtering at fundamental and specific harmonics per Biopac application
    # note #265
    if comb:
        print("... Applying notch filter.")
        ecg_clean = _comb_band_stop(notches, nyquist, ecg_clean, Q)
    # Low pass filtering at 40Hz per Biopac application note #242
    print("... Applying low pass filtering.")
    ecg_clean = nk.signal_filter(ecg_clean, sampling_rate=sampling_rate, highcut=40)

    return ecg_clean


# =============================================================================
# EDA
# =============================================================================
def neuromod_eda_clean(eda_signal, sampling_rate=10000.0, downsampling=None):
    """
    Multiband sequence gradient noise reduction.

    Parameters
    ----------
    eda_signal : array
        The EDA channel.
    sampling_rate : float
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Default to 10000.
    downsampling : int
        The desired sampling frequency (Hz). If None, the signal is not
        resample.
        Default to None.

    Returns
    -------
    eda_clean : array
        The cleaned EDA signal.

    References
    ----------
    Bottenhorn, K. L., Salo, T., Riedel, M. C., Sutherland, M. T.,
        Robinson, J. L., Musser, E. D., & Laird, A. R. (2021).
        Denoising physiological data collected during multi-band,
        multi-echo EPI sequences. bioRxiv, 2021-04.
        https://doi.org/10.1101/2021.04.01.437293

    See also
    --------
    https://neuropsychology.github.io/NeuroKit/functions/eda.html#preprocessing
    """
    # Low pass filtering at 3Hz, order=4
    eda_clean = nk.eda_clean(eda_signal, sampling_rate=sampling_rate)
    # Downsample the signal if specified
    if downsampling is not None:
        eda_clean = nk.signal_resample(
            eda_clean,
            sampling_rate=sampling_rate,
            desired_sampling_rate=downsampling,
        )
        eda_signal = nk.signal_resample(
            eda_signal,
            sampling_rate=sampling_rate,
            desired_sampling_rate=downsampling,
        )

    return eda_signal, eda_clean


# =============================================================================
# RSP
# =============================================================================
def neuromod_rsp_clean(rsp_signal, sampling_rate=10000.0, downsampling=None):
    """
    Clean a PPG signal.

    Parameters
    ----------
    rsp_signal : list, array or Series
        The raw RSP channel.
    sampling_rate : float
        The sampling frequency of `rsp_signal` (in Hz, i.e., samples/second).
        Default to 10000.
    downsampling : int
        The desired sampling frequency (Hz). If None, the signal is not
        resample.
        Default to None.

    Returns
    -------
    rsp_clean : array
        Vector containing the cleaned PPG signal.

    References
    ----------
    Khodadad, D., Nordebo, S., Müller, B., Waldmann, A., Yerworth, R.,
        Becher, T., ... & Bayford, R. (2018). Optimized breath detection
        algorithm in electrical impedance tomography. Physiological
        measurement, 39(9), 094001.
    """
    # Apply bandpass filter
    rsp_clean = nk.rsp_clean(
        rsp_signal, sampling_rate=sampling_rate, method="khodadad2018"
    )
    # Downsample the signal if specified
    if downsampling is not None:
        rsp_clean = nk.signal_resample(
            rsp_clean,
            sampling_rate=sampling_rate,
            desired_sampling_rate=downsampling,
        )
        rsp_signal = nk.signal_resample(
            rsp_signal,
            sampling_rate=sampling_rate,
            desired_sampling_rate=downsampling,
        )

    return rsp_signal, rsp_clean


# =============================================================================
# General functions
# =============================================================================


def _comb_band_stop(notches, nyquist, filtered, Q):
    """
    Series of notch filters aligned with the scanner gradient's harmonics.

    Parameters
    ----------
    notches : dict
        Frequencies to use in the IIR notch filter.
    nyquist : float
        The Nyquist frequency.
    filtered : array
        Data to be filtered.
    Q : int
        The filter quality factor.

    Returns
    -------
    filtered : array
        The filtered signal.

    References
    ----------
    Biopac Systems, Inc. Application Notes: application note 242
        ECG Signal Processing During fMRI
        https://www.biopac.com/wp-content/uploads/app242x.pdf

    See also
    --------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.iirnotch.html
    """
    # band stopping each frequency specified with notches dict
    for notch in notches:
        for i in np.arange(1, int(nyquist / notches[notch])):
            f0 = notches[notch] * i
            w0 = f0 / nyquist
            b, a = signal.iirnotch(w0, Q)
            filtered = signal.filtfilt(b, a, filtered)
    return filtered
