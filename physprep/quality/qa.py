# Migrate computing_sqi function currently in `report.py` here

import traceback
from pathlib import Path

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
    sliding={"duration": 60, "step": 10}
):
    """
    Run processing QC-ing pipeline on specified biosignals, and generate an html report
    containing the quality metrics.

    Parameters
    ----------
    timeseries: dict
        Dictionary containing the processed timeseries for each modality.
    extracted_features: dict
        Dictionary containing the extracted features for each modality.
    sliding : dict
        Dictionary containing the `duration` (in seconds) of the windows in which to
        calculate the SQI. If `step` is not specified, is equal to 0 or to None,
        the fixed windows approach is used for QC-ing. If `step` is specified, a
        sliding window approach is used for QC-ing, and the value corresponds to the
        step between each window (in seconds). If `duration` is None, the SQI are
        computed on the whole run.
        Default to {'duration': 60, 'step': 10}.
    """
    summary, summary_ppg, summary_ecg, summary_eda, summary_rsp = (
        {},
        {},
        {},
        {},
        {},
    )
    filename = filename.replace("physio", "desc-quality-report")

    # Compute metrics on the unsegmented signal
    if not sliding:
        for modality in extracted_features.keys():
            print(f"***Computing quality metrics for {modality} signal***")
            try:
                if modality in ["cardiac_ppg", "cardiac_ecg"]:
                    summary[modality] = {
                        "Overview": sqi_cardiac_overview(extracted_features[modality])
                    }
                    summary[modality].update(
                        sqi_cardiac(
                            timeseries[f"{modality}_clean"],
                            extracted_features[modality],
                        )
                    )
                if modality == "electrodermal":
                    summary[modality] = {
                        "Overview": sqi_eda_overview(
                            len(extracted_features[modality]["ScrPeaks"])
                        )
                    }
                    summary[modality].update(
                        sqi_eda(
                            timeseries[modality],
                            extracted_features[modality],
                        )
                    )
                if modality == "respiratory":
                    summary[modality] = sqi_rsp(timeseries[modality])
            except Exception:
                print(f"Not able to compute {modality}")
                traceback.print_exc()

        # Generate report
        print("***Generating report***")
        generate_report(
            workflow,
            summary,
            timeseries,
            extracted_features,
            outdir,
            filename,
        )

    # Compute metrics on signal segmented in fixed windows
    elif (
        sliding.get("step") == 0 or not sliding.get("step") or sliding.get("step") is None
    ):
        print("Fixed window")
        for modality in extracted_features.keys():
            try:
                window_samples = int(
                    sliding["duration"]
                    * extracted_features[modality]["SamplingFrequency"]
                )
                num_windows = int(len(timeseries[modality]) / window_samples)

                # Descriptive metrics on the windows
                print(f"***Computing quality metrics for {modality} signal***")
                for i in range(num_windows):
                    start = int(i * window_samples)
                    start_idx = int(
                        start / extracted_features[modality]["SamplingFrequency"]
                    )
                    end = int(start + window_samples)
                    end_idx = int(end / extracted_features[modality]["SamplingFrequency"])
                    window = timeseries[modality].iloc[start:end]
                    if modality == "caridac_ppg":
                        summary_ppg[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                            window[f"{modality}_clean"],
                            extracted_features[modality],
                            window=[start, end],
                        )
                    elif modality == "cardiac_ecg":
                        summary_ecg[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                            window[f"{modality}_clean"],
                            extracted_features[modality],
                            window=[start, end],
                        )
                    elif modality == "electrodermal":
                        summary_eda[f"{start_idx}-{end_idx}"] = sqi_eda(
                            window,
                            extracted_features[modality],
                            window=[start, end],
                        )
                    elif modality == "respiratory":
                        summary_rsp[f"{start_idx}-{end_idx}"] = sqi_rsp(
                            window,
                        )
                # Descriptive metrics on the overall run
                if modality == "cardiac_ppg":
                    summary[modality] = {
                        "Overview": sqi_cardiac_overview(extracted_features[modality])
                    }
                    summary[modality].update(summary_ppg)
                elif modality == "cardiac_ecg":
                    summary[modality] = {
                        "Overview": sqi_cardiac_overview(extracted_features[modality])
                    }
                    summary[modality].update(summary_ecg)
                elif modality == "electrodermal":
                    summary[modality] = {
                        "Overview": sqi_eda_overview(
                            len(extracted_features[modality]["ScrPeaks"])
                        )
                    }
                    summary[modality].update(summary_eda)
                elif modality == "respiratory":
                    summary[modality] = summary_rsp

            except Exception:
                print(f"Not able to compute {modality}")
                traceback.print_exc()

        # Generate report
        print("***Generating report***")
        generate_report(
            workflow,
            summary,
            timeseries,
            extracted_features,
            outdir,
            filename,
            window=True,
        )

    # Compute metrics using a sliding window approach
    else:
        print("Sliding window")
        for modality in extracted_features.keys():
            try:
                window_samples = int(
                    sliding["duration"]
                    * extracted_features[modality]["SamplingFrequency"]
                )
                step_samples = int(
                    sliding["step"] * extracted_features[modality]["SamplingFrequency"]
                )
                num_windows = (
                    int((len(timeseries[modality]) - window_samples) / step_samples) + 1
                )

                print(f"***Computing quality metrics for {modality} signal***")
                for i in range(num_windows):
                    start = int(i * step_samples)
                    start_idx = int(
                        start / extracted_features[modality]["SamplingFrequency"]
                    )
                    end = int(start + window_samples)
                    end_idx = int(end / extracted_features[modality]["SamplingFrequency"])
                    window = timeseries[modality].iloc[start:end]
                    if modality == "cardiac_ppg":
                        summary_ppg[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                            window[f"{modality}_clean"],
                            extracted_features[modality],
                            window=[start, end],
                        )
                    elif modality == "cardiac_ecg":
                        summary_ecg[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                            window[f"{modality}_clean"],
                            extracted_features[modality],
                            window=[start, end],
                        )
                    elif modality == "electrodermal":
                        summary_eda[f"{start_idx}-{end_idx}"] = sqi_eda(
                            window,
                            extracted_features[modality],
                            window=[start, end],
                        )
                    elif modality == "respiratory":
                        summary_rsp[f"{start_idx}-{end_idx}"] = sqi_rsp(
                            window,
                        )

                # Descriptive metrics on the overall run
                if modality == "cardiac_ppg":
                    summary[modality] = {
                        "Overview": sqi_cardiac_overview(extracted_features[modality])
                    }
                    summary[modality].update(summary_ppg)
                elif modality == "cardiac_ecg":
                    summary[modality] = {
                        "Overview": sqi_cardiac_overview(extracted_features[modality])
                    }
                    summary[modality].update(summary_ecg)
                elif modality == "electrodermal":
                    summary[modality] = {
                        "Overview": sqi_eda_overview(
                            len(extracted_features[modality]["ScrPeaks"])
                        )
                    }
                    summary[modality].update(summary_eda)
                elif modality == "respiratory":
                    summary[modality] = summary_rsp
            except Exception:
                print(f"Not able to compute {modality}")
                traceback.print_exc()

        # Generate report
        print("***Generating report***")
        generate_report(
            workflow,
            summary,
            timeseries,
            extracted_features,
            outdir,
            filename,
            window=True,
        )


def qa():
    sub_id = sub.name
    """
    path_files = sorted(sub.glob('*'))
    files = [file for file in path_files if 'quality-metrics' in file.name]
    for file in files:
        print(file)
        qc_data = load_json(file)
        qc_dict = {}
        if dataset in ['movie10', 'friends']:
            dataset_idx = 3
        elif dataset in ['shinobi', 'harrypotter']:
            dataset_idx = 4
        filename = '_'.join(file.name.split('_')[:dataset_idx]) + '_desc-quality.json'
        for modality in qc_data.keys():
            if 'Overview' in qc_data[modality].keys():
                del qc_data[modality]['Overview']
            quality = [qc_data[modality][window]['Quality'] for window in qc_data[modality].keys()]
            percentage_pass = quality.count('Acceptable') / len(quality)
            if percentage_pass > 0.80:
                pass_or_fail = 'Pass'
            else:
                pass_or_fail = 'Fail'
            if modality == 'PPG':
                descript = 'Cardiac pulse signal'
            elif modality == 'ECG':
                descript = 'Electrocardiogram signal'
            elif modality == 'RSP':
                descript = 'Respiratory belt signal'
            elif modality == 'EDA':
                descript = 'Electrodermal signal'
            qc_dict[modality] = {
                'Description': descript,
                'PercentageValid': np.round(percentage_pass, 4),
                'QualityAssessment': pass_or_fail
            }
        qc_dict['Description'] = {
            'PercentageValid': 'Percentage of acceptable 1min windows within a run based on modality specific criterion',
            'QualityAssessment': 'Quality assessment of the run. Pass if more than 80% of the run is classified as Acceptable, otherwise Fail'
        }
        # Saving the output
        p = path_deriv / sub_id / 'func'
        p.mkdir(parents=True, exist_ok=True)
        with open(p / filename, 'w') as f:
            json.dump(qc_dict, f, indent=4)
    """
    path_ses = sorted(sub.glob('ses*'))
    for ses in path_ses:
        print(ses)
        ses_id = ses.name
        path_files = sorted(ses.glob('*'))
        files = [file for file in path_files if 'quality-metrics' in file.name]
        for file in files:
            print(file)
            qc_data = load_json(file)
            qc_dict = {}
            if dataset in ['movie10', 'friends']:
                dataset_idx = 3
            elif dataset in ['shinobi', 'mario']:
                dataset_idx = 4
            filename = '_'.join(file.name.split('_')[:dataset_idx]) + '_desc-quality.json'
            for modality in qc_data.keys():
                if 'Overview' in qc_data[modality].keys():
                    del qc_data[modality]['Overview']
                quality = [qc_data[modality][window]['Quality'] for window in qc_data[modality].keys()]
                percentage_pass = quality.count('Acceptable') / len(quality)
                if percentage_pass > 0.80:
                    pass_or_fail = 'Pass'
                else:
                    pass_or_fail = 'Fail'
                if modality == 'PPG':
                    descript = 'Cardiac pulse signal'
                elif modality == 'ECG':
                    descript = 'Electrocardiogram signal'
                elif modality == 'RSP':
                    descript = 'Respiratory belt signal'
                elif modality == 'EDA':
                    descript = 'Electrodermal signal'
                qc_dict[modality] = {
                    'Description': descript,
                    'PercentageValid': np.round(percentage_pass, 4),
                    'QualityAssessment': pass_or_fail
                }
            qc_dict['Description'] = {
                'PercentageValid': 'Percentage of acceptable 1min windows within a run based on modality specific criterion',
                'QualityAssessment': 'Quality assessment of the run. Pass if more than 80% of the run is classified as Acceptable, otherwise Fail'
            }
            # Saving the output
            p = path_deriv / sub_id / ses_id / 'func'
            p.mkdir(parents=True, exist_ok=True)
            with open(p / filename, 'w') as f:
                json.dump(qc_dict, f, indent=4)
