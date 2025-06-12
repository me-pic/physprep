
# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""
Creating simulated BIDS dataset
"""

import json
import shutil
import numpy as np
import pandas as pd
import neurokit2 as nk
from pathlib import Path

import importlib.resources as pkg_resources
from . import boilerplates



def create_bids_dataset(path, nb_sub=1, nb_ses=1, nb_run=1, duration=60, sampling_rate=1000, noise=0.1, modality=['ECG', 'PPG', 'RSP', 'EDA'], concurrent_with='func', overwrite=False):
    """
    Create simulated BIDS dataset for raw physio data

    Parameters
    ----------
    path : str or pathlib.Path
        Directory to save the fake dataset.
    duration : int
        Length of the timeseries in seconds.
    sampling_rate : int or float
        Sampling rate to use for the simulated timeseries.
    noise : float
        Noise level to add to the simulated timeseries.
    modality : list
        List of physiological modalities to include in the dataset. Possible values
        incluce: `ECG`, `PPG`, `RSP`, `EDA`.
    concurrent_with : str
        Specify in which context the peripheral physiological data were recorded. If they
        were recorded with fMRI data, the value of `concurrent_with` should be `func`. If
        they were recorded with eeg data, the value should be `eeg`.
    """
    valid = ['ECG', 'PPG', 'RSP', 'EDA']
    # Check if path exist, otherwise create it
    path = Path(path)
    if overwrite is False and path.exists():
        # Raise error if directory already exists
        raise FileExistsError(f'{path} already exists.')
    else:
        # Create directory
        path.mkdir(parents=True, exist_ok=True)

    # Verify that the specified modality are valid
    if not all(m in valid for m in modality):
        raise ValueError(f'Values in {modality} must be in {valid}.')

    # Create simulated timeseries corresponding to each modality
    for sub in range(nb_sub):
        # Create sub id
        sub_id = '0'*(2-len(str(nb_sub))) + str(sub+1)
        for ses in range(nb_ses):
            # Create ses id
            ses_id = '0'*(3-len(str(nb_ses))) + str(ses+1)
            for run in range(nb_run):
                run_id = '0'*(2-len(str(nb_run))) + str(run+1)
                # Instantiate empty DataFrame and empty dictionary
                data = pd.DataFrame()
                metadata = {}
                cols = []

                # Add time
                data['time'] = np.linspace(
                    0, 
                    duration*sampling_rate,
                    num=duration*sampling_rate
                )
                # Add SamplingFrequency and StartTime
                metadata.update(
                    {
                        'SamplingFrequency': sampling_rate,
                        'StartTime': 0,
                    }
                )

                for mod in modality:
                    if mod == 'ECG':
                        data['ECG'] = nk.ecg_simulate(
                            duration=duration,
                            sampling_rate=sampling_rate,
                            noise=noise,
                            random_state=0
                        )
                        cols.append('ECG')
                    elif mod == 'PPG':
                        data['PPG'] = nk.ppg_simulate(
                            duration=duration,
                            sampling_rate=sampling_rate,
                            ibi_randomness=noise,
                            random_state=0
                        )
                        cols.append('PPG')
                    elif mod == 'RSP':
                        data['RSP'] = nk.rsp_simulate(
                            duration=duration,
                            sampling_rate=sampling_rate,
                            noise=noise,
                            random_state=0
                        )
                        cols.append('RSP')
                    elif mod == 'EDA':
                        data['EDA'] = nk.eda_simulate(
                            duration=duration,
                            sampling_rate=sampling_rate,
                            noise=noise,
                            random_state=0
                        )
                        cols.append('EDA')
                
                # Create subdirectories
                path_tmp = path / f'sub-{sub_id}' / f'ses-{ses_id}' / concurrent_with 
                path_tmp.mkdir(parents=True, exist_ok=True)
                # Save timeseries
                data.to_csv(path / f'sub-{sub_id}' / f'ses-{ses_id}' / concurrent_with / f'sub-{sub_id}_ses-{ses_id}_task-test_run-{run_id}_physio.tsv.gz', sep='\t')
                # Save related metadata
                metadata.update(
                    {
                        'Columns':cols
                    }
                )
                with open(path / f'sub-{sub_id}' / f'ses-{ses_id}' / concurrent_with / f'sub-{sub_id}_ses-{ses_id}_task-test_run-{run_id}_physio.json', "w") as f:
                    json.dump(metadata, f, indent=4)
                    f.close()
    
    # Saving the dataset_description.json file at the root of `path`
    with pkg_resources.files(boilerplates).joinpath('dataset_description.json').open('rb') as template_file:
        with open(path / 'dataset_description.json', 'wb') as dest_file:
            shutil.copyfileobj(template_file, dest_file)