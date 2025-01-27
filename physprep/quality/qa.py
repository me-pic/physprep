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
    summary, summary_tmp, summary_short = ({}, {}, {})

    if isinstance(metadata, str) or isinstance(metadata, Path):
        metadata = load_json(metadata)

    for modality in extracted_features.keys():
        if isinstance(metadata["SamplingFrequency"], list):
            sampling_rate = metadata["SamplingFrequency"][idx]
        else:
            sampling_rate = metadata["SamplingFrequency"]

        # Get data
        timeseries_data = timeseries[modality+"_clean"]if modality+"_clean" in timeseries.keys() else timeseries[modality]

        # Get sliding
        workflow_qa = utils.get_config(workflow[modality]['qa_strategy'], strategy="qa")
        sliding = [step["sliding"] for step in workflow_qa if "sliding" in step.keys()]

        if not sliding:
            # If sliding not specified, consider the timeserie as one big window
            sliding = {'duration': len(timeseries_data)/sampling_rate, 'step': 0}
        else:
            sliding = sliding[0]

        window_samples = int(
            sliding["duration"]
            * sampling_rate
        )
        step_samples = int(
            sliding["step"] * sampling_rate
        )
        if sliding['duration'] == len(timeseries_data)/sampling_rate:
            num_windows = 1
        elif sliding['step'] == 0:
            num_windows = int(len(timeseries_data) / window_samples)
        else:
            num_windows = (
                int((len(timeseries_data) - window_samples) / step_samples) + 1
            )

        print(f"***Computing quality metrics for {modality} signal***")
        for i in range(num_windows):
            if sliding['duration'] != len(timeseries_data)/sampling_rate and sliding['step'] == 0:
                start = int(i * window_samples)
            else:
                start = int(i * step_samples)
            start_idx = int(
                start / sampling_rate
            )
            end = int(start + window_samples)
            end_idx = int(end / sampling_rate)
            
            window = timeseries_data[start:end]

            if modality.lower() in ["ppg", "cardiac_ppg"]:
                summary_tmp[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                    window,
                    extracted_features[modality],
                    window=[start, end],
                )
            elif modality.lower() in ["ecg", "cardiac_ecg"]:
                summary_tmp[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                    window,
                    extracted_features[modality],
                    window=[start, end],
                )
            elif modality.lower() in ["eda", "gsr", "electrodermal"]:
                phasic_key = [key for key in timeseries.keys() if "phasic" in key.lower()]
                tonic_key = [key for key in timeseries.keys() if "tonic" in key.lower()]
                if len(phasic_key)==1:
                    phasic = timeseries[phasic_key[0]][start:end]
                else:
                    phasic = None
                if len(tonic_key)==1:
                    tonic = timeseries[tonic_key[0]][start:end]
                else:
                    tonic = None

                summary_tmp[f"{start_idx}-{end_idx}"] = sqi_eda(
                    window,
                    phasic,
                    tonic,
                    extracted_features[modality],
                    window=[start, end],
                )
            elif modality.lower() in ["rsp", "resp", "respiratory"]:
                # Retrieve or compute the respiratory amplitude
                amplitude_key = [key for key in timeseries.keys() if "amplitude" in key.lower()]
                if len(amplitude_key)==1:
                    amplitude = timeseries[amplitude_key][start:end]
                else:
                    amplitude = nk.rsp_amplitude(
                        np.array(timeseries_data), 
                        np.array(extracted_features[modality]['inhale_max']),
                        troughs=np.array(extracted_features[modality]['exhale_max'])
                    )
                    amplitude = amplitude[start:end]
                # Retrieve or compute the respiratory rate
                rate_key = [key for key in timeseries.keys() if "rate" in key.lower()]
                if len(rate_key)==1:
                    rate = timeseries[rate_key][start:end]
                else:
                    rate = nk.signal_rate(
                        np.array(extracted_features[modality]['exhale_max']),
                        sampling_rate=sampling_rate,
                        desired_length=len(timeseries_data),
                        interpolation_method='monotone_cubic'
                    )
                    rate = rate[start:end]

                summary_tmp[f"{start_idx}-{end_idx}"] = sqi_rsp(
                    window,
                    amplitude,
                    rate
                )

        # Generate short summary
        quality_window = 0
        for window in summary_tmp.keys():
            if summary_tmp[window]['Quality'] == "Acceptable":
                quality_window += 1

        summary_short[modality] = {
            "Description": workflow[modality]['Description'],
            "PercentageValid": np.round(quality_window/len(summary_tmp.keys()), 2),
            "QualityAssessment": "Pass" if quality_window/len(summary_tmp.keys()) > 0.8 else "Fail"
        }

        # Descriptive metrics on the overall run
        if modality.lower() in ["ppg", "ecg", "cardiac_ppg", "cardiac_ecg"]:
            summary[modality] = {
                "Overview": sqi_cardiac_overview(
                    extracted_features[modality],
                    sampling_rate
                )
            }
            summary[modality].update(summary_tmp)
        elif modality.lower() in ["eda", "electrodermal", "gsr"]:
            summary[modality] = {
                "Overview": sqi_eda_overview(
                    len(extracted_features[modality]["scr_peak"])
                )
            }
            summary[modality].update(summary_tmp)
        elif modality.lower() in ["rsp", "resp", "respiratory"]:
            summary[modality] = summary_tmp

    summary_short["Description"] = {
        "PercentageValid": f"Percentage of acceptable {sliding['duration']/60}min windows within a run based on modality specific criterion",
        "QualityAssessment": "Quality assessment of the run. Pass if more than 80% of the run is classified as Acceptable, otherwise Fail"
    }

    return summary, summary_short