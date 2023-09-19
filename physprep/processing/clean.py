# -*- coding: utf-8 -*-
# !/usr/bin/env python
"""
Neuromod cleaning utilities.
"""

import pickle
from pathlib import Path

import neurokit2 as nk
import numpy as np
import pandas as pd
from scipy import signal

from physprep import utils


def preprocessing_workflow(
    data, metadata, workflow_strategy, outdir=None, filename=None, save=True
):
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

    # Remove padding in data if any
    if metadata["StartTime"] > 2e-3:
        data = remove_padding(data)

    # Iterate over content of `workflow_strategy`
    for signal_type in workflow_strategy:
        if signal_type != "trigger":
            if "preprocessing_strategy" in signal and signal[
                "preprocessing_strategy"
            ] not in ["", " ", None]:
                raw, clean, sampling_rate = preprocess_signal(
                    data[signal_type],
                    signal_type["preprocessing_strategy"],
                    signal_type.id,
                    sampling_rate=metadata["SamplingFrequency"],
                )
                clean_signals.update(
                    {
                        signal_type: {
                            f"{signal_type}_raw": raw,
                            f"{signal_type}_clean": clean,
                        }
                    }
                )
                metadata_derivatives.update(
                    {
                        signal_type: {
                            "StartTime": data["time"].loc[0],
                            "SamplingFrequency": sampling_rate,
                            "Columns": list(f"{signal_type}_raw", f"{signal_type}_clean"),
                        }
                    }
                )
            else:
                print(
                    f"No preprocessing strategy specified for {signal_type}. "
                    "The preprocessing step will be skipped."
                )

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
        for clean_signal in clean_signals:
            filename_signal = filename.replace("physio", f"desc-preproc_{clean_signal}")
            df_signal = pd.DataFrame(clean_signals[clean_signal])
            df_signal.to_csv(
                Path(outdir / filename_signal).with_suffix(".tsv.gz"),
                sep="\t",
                index=False,
                compression="gzip",
            )

            # Save metadata on filtered signals
            with open(Path(outdir / filename_signal).with_suffix(".json"), "w") as f:
                pickle.dump(metadata_derivatives[clean_signal], f, indent=4)
                f.close()

        return clean_signals, metadata_derivatives


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
    # Retrieve preprocessing steps
    preprocessing = utils.get_config(preprocessing_strategy, strategy="preprocessing")
    # Iterate over preprocessing steps as defined in the configuration file
    for step in preprocessing:
        if step["step"] == "filtering":
            if step["parameters"]["method"] == "notch":
                pass  # TODO: implement notch filtering
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
            sampling_rate = step["signal_resample"]["desired_sampling_rate"]
        else:
            raise ValueError(
                f"Unknown preprocessing step: {step['step']}. Make sure the "
                "preprocessing strategy is properly defined. For more "
                "details, please refer to the Physprep documentation."
            )

    return raw, signal, sampling_rate


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
