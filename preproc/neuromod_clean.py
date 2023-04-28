# -*- coding: utf-8 -*-
# !/usr/bin/env python
"""
Neuromod cleaning utilities

NOTE : add downsampling option
NOTE : add neuromod_rsp_clean and neuromod_eda_clean functions
"""
import numpy as np
import pandas as pd
import neurokit2 as nk
from scipy import signal

import sys
sys.path.append("../visu")
#import make_pipeline

def neuromod_bio_clean(tsv=None, data=None, h5=None, sampling_rate=1000, pipeline=False):
    """
    Filter biosignals.

    Prepare biosignals for extraction of characteristics of physiological
    activity with a set of filters and smoothing functions

    Parameters
    ----------
    tsv :
        directory of BIDSified biosignal recording
    h5 :
        directory of h5 file
    data (optional) :
        pandas DataFrame object
    sampling_rate : int
        The sampling frequency of `signal_raw` (in Hz, i.e., samples/second).
        Default to 10000.
    pipeline : bool
        If True the pipeline figure will be generated for each signal, including
        the filter applied to each step and the related parameters, and a 10s window
        of the processed signal following each step.
        Default to False.
    
    Return
    ------
    bio_df : DataFrame
    """

    # check input and sanitize
    #if data and tsv is None:
        #raise ValueError("You have to give at least one of the two \n"
                         #"parameters: tsv or df")

    if tsv is not None:
        data = pd.read_csv(tsv, sep='t', compression='gz')

    if h5 is not None:
        data = pd.read_hdf(h5, key='bio_df')
        sampling_rate = pd.read_hdf(h5, key='sampling_rate')

    # sanitize by columns
    if "RSP" in data.keys():
        rsp = data["RSP"]
    else:
        rsp = None
    if "EDA" in data.keys():
        eda = data["EDA"]
    else:
        eda = None
    if "ECG" in data.keys():
        ecg = data["ECG"]
    elif "EKG" in data.keys():
        ecg = data["EKG"]
    else:
        ecg = None
    if "PPG" in data.keys():
        ppg = data["PPG"]
    else:
        ppg = None

    # Keep unkown columns in data
    cols = ["ECG", "EKG", "PPG", "RSP", "EDA"]
    keep_keys = [key for key in data.keys() if key not in cols]
    if len(keep_keys) != 0:
        keep = data[keep_keys]
    else:
        keep = None

    # initialize output
    bio_df = pd.DataFrame()

    # sanitize input signals
    # PPG_
    if ppg is not None:
        ppg = nk.as_vector(ppg)
        ppg_clean, ppg_pipeline = neuromod_ppg_clean(ppg, sampling_rate=sampling_rate, pipeline=pipeline)

        bio_df = pd.concat([bio_df, ppg_clean], axis=1)
    # ECG
    if ecg is not None:
        ecg = nk.as_vector(ecg)
        ecg_clean, ecg_pipeline = neuromod_ecg_clean(ecg, sampling_rate=sampling_rate, pipeline=pipeline)
        bio_df = pd.concat([bio_df, ecg_clean], axis=1)

    # RSP
    if rsp is not None:
        rsp = nk.as_vector(rsp)
        rsp_clean, rsp_pipeline = nk.rsp_clean(rsp, sampling_rate=sampling_rate) #pipeline=pipeline 
        bio_df = pd.concat([bio_df, rsp_clean], axis=1)

    # EDA
    if eda is not None:
        eda = nk.as_vector(eda)
        eda_clean, eda_pipeline = nk.eda_clean(eda, sampling_rate=sampling_rate) #pipeline=pipeline
        bio_df = pd.concat([bio_df, eda_clean], axis=1)

    return bio_df
# ======================================================================
# Photoplethysmograph (PPG)
# =======================================================================

def neuromod_ppg_clean(ppg_signal, sampling_rate=10000, method='nabian2018'):
    ppg_cleaned = nk.ppg_clean(ppg_signal, sampling_rate=sampling_rate,
                               method=method)

    return ppg_cleaned
# ======================================================================
# Electrocardiogram (ECG)
# =======================================================================

def neuromod_ecg_clean(ecg_signal, trigger_pulse, sampling_rate=10000., method='biopac', me=False, pipeline=False):
    """
    Clean an ECG signal.

    Prepare a raw ECG signal for R-peak detection with the specified method.

    Parameters
    ----------
    ecg_signal : list, array or Series
        The raw ECG channel.
    sampling_rate : int
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Defaults to 10000.
    method : str
        The processing pipeline to apply. Defaults to 'biopac'.
    me : bool
        Specify if the MRI sequence used was the multi-echo (True) or the singlue-echo (False). Defaults to False.
    Returns
    -------
    array
        Vector containing the cleaned ECG signal.
    """
    make_pipeline = dict()
    step_idx = 1
    window_duration = 10
    start = window_duration*sampling_rate
    end = start + window_duration*sampling_rate

    if me:
        tr=2.65
        mb=2
        slices=70
    else:
        tr=1.49
        mb=4
        slices=60

    method = method.lower()  # remove capitalised letters
    if method in ["schmidt", "schmidt2018"]:
        clean = _ecg_clean_schmidt(ecg_signal, sampling_rate)
    else:
        timeseries = pd.DataFrame(columns=['ECG','Trigger', 'seconds'])
        time = np.arange(len(ecg_signal))/sampling_rate
        timeseries['ECG'] = ecg_signal
        timeseries['Trigger'] = trigger_pulse
        timeseries['seconds'] = time
        if method in ["biopac"]:
            clean = _ecg_clean_biopac(timeseries, sampling_rate)
        if method in ['bottenhorn', 'bottenhorn2022']:
            #Remove respiration-related noise using a 2Hz highpass filter
            print('---Cleaning respiration-related noise---')
            ecg_signal_hp, info_filter = butter_highpass_filter(ecg_signal, 2., sampling_rate)
            if pipeline:
                make_pipeline.update({
                    step_idx: info_filter.update({
                        'signal' : ecg_signal_hp[start:end]
                        })
                    })
                step_idx+=1
            #Apply comb band pass filter with Bottenhorn correction
            print('---Applying the corrected comb band pass filter---')
            ecg_bottenhorn, info_filter = _ecg_clean_bottenhorn(ecg_signal_hp, sampling_rate=sampling_rate, tr=tr, mb=mb, slices=slices)
            if pipeline:
                make_pipeline.update({
                    step_idx: info_filter.update({
                        'signal' : ecg_bottenhorn[start:end]
                        })
                    })
                step_idx+=1        
        clean = bandpass_filter(ecg_bottenhorn, f0=24.0, Q=6, low=3, high=34, order=5, sampling_rate=sampling_rate)
        if pipeline:
            make_pipeline.update({
                step_idx: info_filter.update({
                    'signal' : ecg_bottenhorn[start:end]
                    })
                })
       
    return clean, make_pipeline


# =============================================================================
# ECG internal : Schmidt et al. 2016
# =============================================================================
def _ecg_clean_schmidt(ecg_signal, sampling_rate=10000):
    """
    from Schmidt, M., Krug, J. W., & Rose, G. (2016).
    Reducing of gradient induced artifacts on the ECG signal during MRI
    examinations using Wilcoxon filter.
    Current Directions in Biomedical Engineering.
    https://doi.org/10.1515/cdbme-2016-0040
    """
    # enveloppe at least 100 ms in samples
    env = int(0.01*sampling_rate)
    # initialize empty array
    ecg_clean = np.empty(len(ecg_signal), dtype='float64')

    # iterate through each sample
    for current_sample in range(env, (len(ecg_signal)-1)):
        # window to convolve median smoothing operation
        past_samples = current_sample-env
        window = np.array(ecg_signal[past_samples:current_sample])

        # compute the walsh averages with Wilcoxon method
        # initialize empty array half the size of window
        walsh_arr = np.empty(int(len(window)/2))

        # iterate from index 1 to half of window, and -1 to -(half of window)
        # in order to perform reduced computation walsh averages
        for i in range(0, (len(walsh_arr))):
            # -1 is last element and 0 is first
            r_i = -(i+1)
            # the last element is index of middle of window
            if i != len(walsh_arr)-1:
                # mirror elements in window are conjugated until middle index
                walsh_arr[i] = (window[i] + window[r_i])/2

            # middle value in window is last value of walsh array and stays as
            else:
                walsh_arr[i] = window[i]
        # compute the median of walsh averages array
        ecg_clean[current_sample] = np.median(walsh_arr)
    # Bandpass filtering
    ecg_clean = nk.signal_filter(ecg_clean, lowcut=0.05, highcut=45,
                                 method='bessel', order=5)

    return ecg_clean

# =============================================================================
# ECG internal : biopac recommendations
# =============================================================================
def _ecg_clean_biopac(timeseries, sampling_rate=10000., tr=1.49, slices=60, Q=100):
    """
    Clean the ECG signal using a comb band stop filter implemented according to Biopac
    application note on ECG signal processing during fMRI.

    Parameters
    ----------
    timeseries : DataFrame
        DataFrame containing the ECG signal and the trigger timeseries
    sampling_rate : int
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Defaults to 10000.
    tr : float
        TR 
    slice : int
        Number of slices
    Q : float
        Quality factor
    
    Returns
    -------
    filtered : array
        Filtered signal
    info : dict
        Dictionary containing filter info

    Reference
    ---------
    Biopac Systems, Inc. Application Notes: application note 242
    ECG Signal Processing During fMRI
    https://www.biopac.com/wp-content/uploads/app242x.pdf

    See also
    --------
    https://neuropsychology.github.io/NeuroKit/functions/signal.html#signal-filter
    """
    # Setting scanner sequence parameters
    nyquist = np.float64(sampling_rate/2)
    notches = {'slices': slices / tr,
               'tr': 1 / tr}
    # find trigger timing
    triggers = timeseries[timeseries['Trigger'] > 4].index.values
    print(triggers)
    # remove baseline wandering
    filtered = nk.signal_filter(timeseries['ECG'][triggers[1]:triggers[-1]],
                                sampling_rate=int(sampling_rate), lowcut = 2)
    # Filtering at specific harmonics, with trigger timing info
    filtered, info = comb_band_stop(notches, nyquist, filtered, Q, sampling_rate)
    # bandpass filtering
    filtered = nk.signal_filter(filtered, sampling_rate=sampling_rate, lowcut=2,
                                 highcut=20, method='butter', order=5)

    return filtered, info

def _ecg_clean_bottenhorn(ecg_signal, sampling_rate=10000., tr=1.49, mb=4, slices=60, Q=100):
    """
    Clean the ECG signal using a comb band stop filter implementing the correction for multiband 
    factor proposed in Bottenhorn et al. 2021

    Parameters
    ----------
    ecg_signal : list, array or Series
        The ECG channel to filtered.
    sampling_rate : int
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Defaults to 10000.
    tr : float
        TR 
    mb : int
        Multiband acceleration factor
    slice : int
        Number of slices
    Q : float
        Quality factor

    Returns
    -------
    filtered : array
        Filtered signal
    info : dict
        Dictionary containing filter info

    Reference
    ---------
    Bottenhorn, K. L., Salo, T., Riedel, M. C., Sutherland, M. T., Robinson, J. L., Musser, E. D., 
        & Laird, A. R. (2021). Denoising physiological data collected during multi-band, multi-echo 
        EPI sequences. bioRxiv, 2021-04. https://doi.org/10.1101/2021.04.01.437293
    
    https://github.com/62442katieb/mbme-physio-denoising/blob/main/notebooks/denoising_eda.ipynb
    
    See also
    --------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.iirnotch.html 
    https://neuropsychology.github.io/NeuroKit/functions/signal.html#signal-filter 
    """
    # Setting scanner sequence parameters
    nyquist = np.float64(sampling_rate/2)
    notches = {'slices': slices / mb / tr,
               'tr': 1 / tr}
    # remove baseline wandering
    filtered = nk.signal_filter(ecg_signal,
                                sampling_rate=int(sampling_rate), lowcut = 2)
    # Filtering at specific harmonics, with trigger timing info
    filtered, info = comb_band_stop(notches, nyquist, filtered, Q, sampling_rate)
    # bandpass filtering
    filtered = nk.signal_filter(filtered, sampling_rate=sampling_rate, lowcut=2,
                                 highcut=20, method='butter', order=5)

    return filtered, info

# =============================================================================
# EDA
# =============================================================================
def _eda_clean_bottenhorn(eda_signal, sampling_rate=10000., Q=100, mb=4, tr=1.49, slices=60):

    notches = {'slices': slices / mb / tr,'tr': 1 / tr}

    #hp_eda = butter_highpass_filter(scan1['EDA'], 1, fs, order=5)
    bottenhorn_filtered = eda_signal
    for notch in notches:
        bottenhorn_filtered = comb_band_stop(notches[notch], np.float64(sampling_rate/2), bottenhorn_filtered, Q, sampling_rate)


# =============================================================================
# General functions
# =============================================================================

def comb_band_stop(notches, nyquist, filtered, Q, sampling_rate):
    """
    A serie of notch filters aligned with the scanner gradient's harmonics
    
    Parameters
    ----------
    notches : 
    nyquist : 
    filtered : 
    Q : float
        Quality factor
    sampling_rate : int
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Defaults to 10000.

    Returns
    -------
    filtered : 
    dict
        Dictionary containing filter info
    
    Reference
    ---------
    Biopac Systems, Inc. Application Notes: application note 242
    ECG Signal Processing During fMRI
    https://www.biopac.com/wp-content/uploads/app242x.pdf
    """
    # Save argument values 
    info = locals()
    # Compute filter
    ## band stoping each frequency specified with notches dict
    for notch in notches:
        for i in np.arange(1, (nyquist / notches[notch])):
            f0 = notches[notch] * i
            w0 = f0/nyquist
            b,a = signal.iirnotch(w0, Q)
            filtered = signal.filtfilt(b, a, filtered)
    return filtered, {'parameters': info, 'filter': 'Comb band stop'}

def consecutive(data, stepsize=0.000501):
    """
    reference: https://github.com/62442katieb/mbme-physio-denoising/blob/main/notebooks/denoising_eda.ipynb
    """
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def butter_freq_pass(cutoff, sampling_rate, order=5, btype='high'):
    """
    Parameters
    ----------
    cutoff : float
    sampling_rate : int
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Defaults to 10000.
    order
    btype : str

    Returns
    -------
    b
    a

    Reference
    ---------
    https://github.com/62442katieb/mbme-physio-denoising/blob/main/notebooks/denoising_eda.ipynb
    
    See also
    --------

    """
    nyq = 0.5 * sampling_rate
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype=btype, analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, sampling_rate, order=5):
    """
    Parameters
    ----------
    data 
    cutoff : float
    sampling_rate : int
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Defaults to 10000.
    order

    Returns
    -------
    y : 
    dict
        Dictionary containing filter info

    Reference
    ---------
    https://github.com/62442katieb/mbme-physio-denoising/blob/main/notebooks/denoising_eda.ipynb

    See also
    --------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
    """
    # Save argument values 
    info = locals()
    del info['data'] # Remove 'data' from info dictionnary
    # Compute filter
    b, a = butter_freq_pass(cutoff, sampling_rate, order=order, btype='high')
    y = signal.filtfilt(b, a, data)
    return y,  {'parameters': info, 'filter': 'Butterworth high pass'} #, 'signal': y[]}

def butter_lowpass_filter(data, cutoff, sampling_rate, order=5):
    """
    Parameters
    ----------
    data 
    cutoff : float
    sampling_rate : int
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Defaults to 10000.
    order : 

    Returns
    -------
    y : 
    dict
        Dictionary containing filter info

    Reference
    ---------
    https://github.com/62442katieb/mbme-physio-denoising/blob/main/notebooks/denoising_eda.ipynb

    See also
    --------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
    """
    # Save argument values 
    info = locals()
    del info['data'] # Remove 'data' from info dictionnary
    # Compute filter
    b, a = butter_freq_pass(cutoff, sampling_rate, order=order, btype='low')
    y = signal.filtfilt(b, a, data)
    return y, {'parameters': info, 'filter': 'Butterworth low pass'}

def fourier_freq(timeseries, d, fmax):
    """
    Parameters
    ----------
    timeseries
    d
    fmax

    Returns
    -------
    fft
    fft_db
    freq
    limit

    Reference
    ---------
    https://github.com/62442katieb/mbme-physio-denoising/blob/main/notebooks/denoising_eda.ipynb
    
    See also
    --------
    https://numpy.org/doc/stable/reference/generated/numpy.fft.fft.html
    https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html
    """
    fft = np.fft.fft(timeseries)
    freq = np.fft.fftfreq(timeseries.shape[-1], d=d)
    fft_db = 10 * np.log10(abs(fft))
    limit = np.where(freq >= fmax)[0][0]
    return fft, fft_db, freq, limit

def bandpass_filter(data, f0=24.0, Q=6, low=3, high=34, order=5, sampling_rate=10000):
    """
    Parameters
    ----------
    data
    f0
    Q : float
        Quality factor
    low
    high
    order
    sampling_rate : int
        The sampling frequency of `data` (in Hz, i.e., samples/second).
        Defaults to 10000.

    Return
    ------
    ecg_clean
    
    See also
    --------
    """
    nyquist = sampling_rate / 2
    w0 = f0 / nyquist
    sos = signal.butter(order, [low / nyquist, high / nyquist], btype='band', output='sos')
    ecg_clean = signal.sosfiltfilt(sos, data)
    return ecg_clean
