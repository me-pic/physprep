# -*- coding: utf-8 -*-
# !/usr/bin/env python
"""
Neuromod processing utilities
"""
# dependencies
import math

import pandas as pd
import numpy as np
import click
import os
import json
# high-level processing utils
from neurokit2 import eda_process, rsp_process, ecg_peaks,  ppg_findpeaks, ecg_process
from systole.correction import correct_rr, correct_peaks
from heartpy import process, enhance_peaks, exceptions
# signal utils
from systole.utils import input_conversion
from neurokit2.misc import as_vector
from neurokit2 import signal_rate, signal_fixpeaks, signal_filter
from neurokit2.signal.signal_formatpeaks import _signal_from_indices
# home brewed cleaning utils
from neuromod_clean import neuromod_ecg_clean


def neuromod_bio_process(tsv=None, h5=None, df=None, sampling_rate=10000):
    """
    Process biosignals.

    tsv :
        directory of BIDSified biosignal recording
    df (optional) :
        pandas DataFrame object
    """
    #if df and tsv and h5 is None:
        #raise ValueError("You have to give at least one of the two \n"
            #             "parameters: tsv or df")

    if tsv is not None:
        df = pd.read_csv(tsv, sep='\t', compression='gzip')
        print("Reading TSV file")

    if h5 is not None:
        df = pd.read_hdf(h5, key='bio_df')
        sampling_rate = pd.read_hdf(h5, key='sampling_rate')
    # initialize returned objects
    if df is not None:
        print("Reading pandas DataFrame")

        bio_info = {}
        bio_df = pd.DataFrame()

        # initialize each signals
        ppg_raw = df["PPG"]
        rsp_raw = df["RSP"]
        eda_raw = df["EDA"]
        ecg_raw = df["ECG"]

        # ppg
        ppg, ppg_info = neuromod_ppg_process(ppg_raw, sampling_rate=sampling_rate)
        bio_info.update(ppg_info)
        bio_df = pd.concat([bio_df, ppg], axis=1)

        # ecg
       # ecg, ecg_info = neuromod_ecg_process(ecg_raw, sampling_rate=sampling_rate)
        #bio_info.update(ecg_info)
        bio_df = pd.concat([bio_df, ecg_raw], axis=1)

        #  rsp
        rsp, rsp_info = rsp_process(rsp_raw, sampling_rate=sampling_rate,
                                    method='khodadad2018')
        bio_info.update(rsp_info)
        bio_df = pd.concat([bio_df, rsp], axis=1)
        print("Respiration workflow: done")

        #  eda
        eda, eda_info = eda_process(eda_raw, sampling_rate, method='neurokit')
        bio_info.update(eda_info)
        bio_df = pd.concat([bio_df, eda], axis=1)
        print("Electrodermal activity workflow: done")

        # return a dataframe
        bio_df['TTL'] = df['TTL']
        bio_df['time'] = df['time']

    return(bio_df, bio_info)


def neuromod_ppg_process(ppg_raw, sampling_rate=10000):
    """
    Process PPG signal.

    Custom processing function for neuromod PPG acquisition

    Parameters
    -----------
    ppg_raw : vector
        The raw PPG channel.
    sampling_rate : int
        The sampling frequency of `ppg_signal` (in Hz, i.e., samples/second).
        Defaults to 10000.
    Returns
    -------
    signals : DataFrame
        A DataFrame containing the cleaned ppg signals.
        - *"PPG_Raw"*: the raw signal.
        - *"PPG_Clean"*: the cleaned signal.
        - *"PPG_Rate"*: the heart rate as measured based on PPG peaks.
        - *"PPG_Peaks"*: the PPG peaks marked as "1" in a list of zeros.
    info : dict
        containing list of intervals between peaks
    """
    ppg_signal = as_vector(ppg_raw)

    # Clean signal
    ppg_cleaned = signal_filter(ppg_signal, sampling_rate=sampling_rate,
                                lowcut=0.5, highcut=8, order=3)
    print('PPG Cleaned')
    
    # heartpy
    print("HeartPy processing started")
    wd, m = process(ppg_cleaned, sampling_rate, reject_segmentwise=True,
                    interp_clipping=True, report_time=True)
        
    cumsum=0
    rejected_segments=[]
    for i in wd['rejected_segments']:
        cumsum += int(np.diff(i)/sampling_rate)
        rejected_segments.append((int(i[0]),int(i[1])))
    print('Heartpy found peaks')
    
    # Find peaks
    info = ppg_findpeaks(ppg_cleaned, sampling_rate=sampling_rate)
    info['PPG_Peaks'] = info['PPG_Peaks'].tolist()
    peak_list_nk = _signal_from_indices(info['PPG_Peaks'], desired_length=len(ppg_cleaned))
    print('Neurokit found peaks')
    
    # peak to intervals
    rr = input_conversion(info['PPG_Peaks'], input_type='peaks_idx', output_type='rr_ms', sfreq=sampling_rate)
    
    # correct beat detection
    corrected, (nMissed, nExtra, nEctopic, nShort, nLong) = correct_rr(rr) 
    corrected_peaks = correct_peaks(peak_list_nk, n_iterations=4)
    print('systole corrected RR series')
    
    rate = signal_rate(info['PPG_Peaks'], sampling_rate=sampling_rate,
                       desired_length=len(ppg_signal))
    
    
    # sanitize info dict    
    info.update({'PPG_ectopic': nEctopic, 'PPG_short': nShort, 'PPG_long': nLong, 'PPG_extra': nExtra, 'PPG_missed': nMissed,
                 'PPG_clean_rr_systole': corrected.tolist(),'PPG_clean_rr_hp': [float(v) for v in wd['RR_list_cor']],
                 'PPG_rejected_segments': rejected_segments, 
                 'PPG_cumulseconds_rejected': int(cumsum), 
                 'PPG_%_rejected_segments': float(cumsum/(len(ppg_signal)/sampling_rate))})

    # Prepare output  
    signals = pd.DataFrame(
                {"PPG_Raw": ppg_signal, "PPG_Clean": ppg_cleaned, "PPG_Peaks_NK": peak_list_nk,
                 "PPG_Peaks_Systole": corrected_peaks['clean_peaks'],
                 "PPG_Rate": rate}
    )

    return signals, info

def neuromod_ecg_process(ecg_raw, trigger_pulse, sampling_rate=10000, method='bottenhorn'):
    """
    Process neuromod ECG.

    Custom processing for neuromod ECG acquisition.

    Parameters
    -----------
    ecg_raw : vector
        The raw ECG channel.
    sampling_rate : int
        The sampling frequency of `ecg_signal` (in Hz, i.e., samples/second).
        Defaults to 10000.
    method : str
        The processing pipeline to apply. Defaults to 'fmri'
    Returns
    -------
    signals : DataFrame
        A DataFrame containing the cleaned ppg signals.
        - *"ECG_Raw"*: the raw signal.
        - *"ECG_Clean"*: the cleaned signal.
        - *"ECG_Rate"*: the heart rate as measured based on PPG peaks.
    info : dict
        containing list of peaks
    """
    ecg_signal = as_vector(ecg_raw)

    # prepare signal for processing
    ecg_cleaned = neuromod_ecg_clean(ecg_signal, trigger_pulse, sampling_rate=10000, method=method, me=True)

    #heartpy
    print("HeartPy processing started")
    wd, m = process(ecg_cleaned, sampling_rate, reject_segmentwise=True,
                    interp_clipping=True, report_time=True)
    cumsum=0
    rejected_segments=[]
    for i in wd['rejected_segments']:
        cumsum += int(np.diff(i)/sampling_rate)
        rejected_segments.append((int(i[0]),int(i[1])))
    print('Heartpy found peaks')
    
    # Find peaks
    print("Neurokit processing started")
    _, info = ecg_peaks(ecg_cleaned=ecg_cleaned, 
                        sampling_rate=sampling_rate, 
                        method='nabian2018',
                        correct_artifacts=True)
    info['ECG_R_Peaks'] = info['ECG_R_Peaks'].tolist()
    peak_list_nk = _signal_from_indices(info['ECG_R_Peaks'], desired_length=len(ecg_cleaned))
    print('Neurokit found peaks')
    
    # peak to intervals
    rr = input_conversion(info['ECG_R_Peaks'], input_type='peaks_idx', output_type='rr_ms', sfreq=sampling_rate)
    
    # correct beat detection
    corrected, (nMissed, nExtra, nEctopic, nShort, nLong) = correct_rr(rr) 
    corrected_peaks = correct_peaks(peak_list_nk, n_iterations=4)
    print('systole corrected RR series')
    # Compute rate based on peaks
    rate = signal_rate(info['ECG_R_Peaks'], sampling_rate=sampling_rate,
                       desired_length=len(ecg_signal))

    # sanitize info dict    
    info.update({'ECG_ectopic': nEctopic, 'ECG_short': nShort, 'ECG_long': nLong, 'ECG_extra': nExtra, 'ECG_missed': nMissed,
                'ECG_clean_rr_systole': corrected.tolist(),'ECG_clean_rr_hp': [float(v) for v in wd['RR_list_cor']],
                'ECG_rejected_segments': rejected_segments, 
                'ECG_cumulseconds_rejected': int(cumsum), 
                'ECG_%_rejected_segments': float(cumsum/(len(ecg_signal)/sampling_rate))})
    # Prepare output  
    signals = pd.DataFrame(
                {"ECG_Raw": ecg_signal, "ECG_Clean": ecg_cleaned, "ECG_Peaks_NK": peak_list_nk,
                "ECG_Peaks_Systole": corrected_peaks['clean_peaks'],
                "ECG_Rate": rate})

    return signals, info

def neuromod_eda_process():
    return

def load_json(filename):
    tmp = open(filename)
    data = json.load(tmp)
    tmp.close()
    
    return data

def load_segmented_runs(source, sub, ses):
    """
    """
    data_tsv, filenames = [], []
    files_tsv = [f for f in os.listdir(os.path.join(source, sub, ses)) if 'tsv.gz' in f]
    #Remove files ending with _01
    files_tsv = [f for f in files_tsv if '_01.' not in f]
    files_tsv.sort()

    for tsv in files_tsv:
        filename = tsv.split(".")[0]
        filenames.append(filename)
        print(f"---Reading data for {sub} {ses}: run {filename[-2:]}---")
        json = filename+".json"
        print("Reading json file")
        data_json = load_json(os.path.join(source, sub, ses, json))
        print("Reading tsv file")
        data_tsv.append(pd.read_csv(os.path.join(source, sub, ses, tsv),
                                    sep="\t", compression="gzip",
                                    names=data_json["Columns"]))
    
    return data_tsv, filenames
        
@click.command()
@click.argument('source', type=str)
@click.argument('sub', type=str)
@click.argument('ses', type=str)
@click.argument('outdir', type=str)
@click.argument('save', type=bool)
def process_ppg_data(source, sub, ses, outdir, save=True):
    """
    """
    data_tsv, filenames_tsv = load_segmented_runs(source, sub, ses)
    for idx, d in enumerate(data_tsv):
        print(f"---Processing PPG signal for {sub} {ses}: run {filenames_tsv[idx][-2:]}---")
        signals, info = neuromod_ppg_process(d['PPG'], sampling_rate=10000)
        if save:
            print("Saving processed data")
            signals.to_csv(os.path.join(outdir, sub, ses, f"{filenames_tsv[idx]}"+"_ppg_signals"+".tsv"), sep="\t")
            with open(os.path.join(outdir, sub, ses, f"{filenames_tsv[idx]}"+"_ppg_info"+".json"), 'w') as fp:
                json.dump(info, fp)

    return signals, info



@click.command()
@click.argument('source', type=str)
@click.argument('sub', type=str)
@click.argument('ses', type=str)
@click.argument('outdir', type=str)
@click.argument('save', type=bool)
def process_ecg_data(source, sub, ses, outdir, save=True):
    """
    """
    data_tsv, filenames_tsv = load_segmented_runs(source, sub, ses)
    for idx, d in enumerate(data_tsv):
        print(f"---Processing ECG signal for {sub} {ses}: run {filenames_tsv[idx][-2:]}---")
        print('--Cleaning the signal---')
        signals, info = neuromod_ecg_process(d['ECG'], d['TTL'], sampling_rate=10000, method='bottenhorn')
        if save:
            print("Saving processed data")
            signals.to_csv(os.path.join(outdir, sub, ses, f"{filenames_tsv[idx]}"+"_ecg_signals"+".tsv"), sep="\t")
            with open(os.path.join(outdir, sub, ses, f"{filenames_tsv[idx]}"+"_ecg_info"+".json"), 'w') as fp:
                json.dump(info, fp)

    return signals, info


@click.command()
@click.argument('source', type=str)
@click.argument('sub', type=str)
@click.argument('ses', type=str)
@click.argument('outdir', type=str)
@click.argument('save', type=bool)
def process_rsp_data(source, sub, ses, outdir, save =True):
    """
    """
    data_tsv, filenames_tsv = load_segmented_runs(source, sub, ses)
    for idx, d in enumerate(data_tsv):
        print(f"---Processing RSP signal for {sub} {ses}: run {filenames_tsv[idx][-2:]}---")
        signals, _ = rsp_process(d['RSP'], sampling_rate=10000, method='khodadad2018')
        if save:
            print("Saving processed data")
            signals.to_csv(os.path.join(outdir, sub, ses, f"{filenames_tsv[idx]}"+"_rsp_signals"+".tsv"), sep="\t")

    return signals

if __name__ == "__main__":
    #PPG processing pipeline
    #process_ppg_data()
    #ECG processing pipeline
    process_ecg_data()
    #RSP processing pipeline
    #process_rsp_data()
