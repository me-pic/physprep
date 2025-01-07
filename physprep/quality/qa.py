# Migrate computing_sqi function currently in `report.py` here

import traceback
import numpy as np
import neurokit2 as nk
from pathlib import Path

from physprep import utils
from physprep.quality.time_sqi import (
    sqi_cardiac,
    sqi_cardiac_overview,
    sqi_eda,
    sqi_eda_overview,
    sqi_rsp,
)


def computing_sqi(
    workflow,
    timeseries,
    extracted_features,
    metadata
):
    """
    Run processing QC-ing pipeline on specified biosignals, and generate an html report
    containing the quality metrics.

    Parameters
    ----------
    workflow : dict
        Dictionary containing the qc metrics to compute for each modality.
    timeseries: dict
        Dictionary containing the processed timeseries for each modality.
    extracted_features: dict
        Dictionary containing the extracted features for each modality.
    metadata : dict or pathlib.Path
        The metadata associated with the cleaned physiological data, if data has been
        cleaned. Otherwise, the metadata associated with the raw physiological data
        (i.e., the outputed json file from Phys2Bids).
    """
    summary, summary_tmp = ({}, {})

    if isinstance(metadata, str) or isinstance(metadata, Path):
        metadata = load_json(metadata)

    for modality in extracted_features.keys():
        if isinstance(metadata["SamplingFrequency"], list):
            sampling_rate = metadata["SamplingFrequency"][idx]
        else:
            sampling_rate = metadata["SamplingFrequency"]

        # Get sliding
        workflow_qa = utils.get_config(workflow[modality]['qa_strategy'], strategy="qa")
        sliding = [step["sliding"] for step in workflow_qa if "sliding" in step.keys()]

        if not sliding:
            # If sliding not specified, consider the timeserie as one big window
            sliding = {'duration': len(timeseries[modality+"_clean"])/sampling_rate, 'step': 0}
        else:
            sliding = sliding[0]

        window_samples = int(
            sliding["duration"]
            * sampling_rate
        )
        step_samples = int(
            sliding["step"] * sampling_rate
        )
        if sliding['duration'] == len(timeseries[modality+"_clean"])/sampling_rate:
            num_windows = 1
        elif sliding['step'] == 0:
            num_windows = int(len(timeseries[modality+"_clean"]) / window_samples)
        else:
            num_windows = (
                int((len(timeseries[modality+"_clean"]) - window_samples) / step_samples) + 1
            )

        print(f"***Computing quality metrics for {modality} signal***")
        for i in range(num_windows):
            if sliding['duration'] != len(timeseries[modality+"_clean"])/sampling_rate and sliding['step'] == 0:
                start = int(i * window_samples)
            else:
                start = int(i * step_samples)
            start_idx = int(
                start / sampling_rate
            )
            end = int(start + window_samples)
            end_idx = int(end / sampling_rate)
            window = timeseries[modality+"_clean"][start:end]

            if modality in ["ppg", "cardiac_ppg"]:
                summary_tmp[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                    window,
                    extracted_features[modality],
                    window=[start, end],
                )
            elif modality in ["ecg", "cardiac_ecg"]:
                summary_tmp[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                    window,
                    extracted_features[modality],
                    window=[start, end],
                )
            elif modality in ["eda", "gsr", "electrodermal"]:
                if modality+"_phasic" in timeseries.keys():
                    phasic = timeseries[modality+"_phasic"][start:end]
                else:
                    phasic = None
                if modality+"_tonic" in timeseries.keys():
                    tonic = timeseries[modality+"_tonic"][start:end]
                else:
                    tonic = None

                summary_tmp[f"{start_idx}-{end_idx}"] = sqi_eda(
                    window,
                    phasic,
                    tonic,
                    extracted_features[modality],
                    window=[start, end],
                )
            elif modality in ["rsp", "resp", "respiratory"]:
                # Retrieve or compute the respiratory amplitude
                if modality+"_amplitude" in timeseries.keys():
                    amplitude = timeseries[modality+"_amplitude"][start:end]
                else:
                    amplitude = nk.rsp_amplitude(
                        timeseries[modality+"_clean"], 
                        np.array(extracted_features[modality]['inhale_max']),
                        troughs=np.array(extracted_features[modality]['exhale_max'])
                    )
                    amplitude = amplitude[start:end]
                # Retrieve or compute the respiratory rate
                if modality+"_rate" in timeseries.keys():
                    rate = timeseries[modality+"_rate"][start:end]
                else:
                    rate = nk.signal_rate(
                        np.array(extracted_features[modality]['exhale_max']),
                        sampling_rate=sampling_rate,
                        desired_length=len(timeseries[modality+"_clean"]),
                        interpolation_method='monotone_cubic'
                    )
                    rate = rate[start:end]

                summary_tmp[f"{start_idx}-{end_idx}"] = sqi_rsp(
                    window,
                    amplitude,
                    rate
                )

        # Descriptive metrics on the overall run
        if modality in ["ppg", "ecg", "cardiac_ppg", "cardiac_ecg"]:
            summary[modality] = {
                "Overview": sqi_cardiac_overview(
                    extracted_features[modality],
                    sampling_rate
                )
            }
            summary[modality].update(summary_tmp)
        elif modality in ["eda", "electrodermal", "gsr"]:
            summary[modality] = {
                "Overview": sqi_eda_overview(
                    len(extracted_features[modality]["scr_peak"])
                )
            }
            summary[modality].update(summary_tmp)
        elif modality in ["rsp", "resp", "respiratory"]:
            summary[modality] = summary_tmp
            #except Exception:
            #    print(f"Not able to compute {modality}")
            #    traceback.print_exc()

    return summary