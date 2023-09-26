# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""Physiological data quality assessment"""

import traceback
from pathlib import Path

from physprep.quality.time_sqi import (
    sqi_cardiac,
    sqi_cardiac_overview,
    sqi_eda,
    sqi_eda_overview,
    sqi_rsp,
)
from physprep.visu.plot_signals import generate_plot


def computing_sqi(
    workflow,
    timeseries,
    extracted_features,
    outdir,
    filename,
    sliding={"duration": 60, "step": 10},
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
    outdir : str or pathlib.Path
        Path to the output derivatives directory.
    filename : str or pathlib.Path
        Name of the file to save the preprocessed physiological data.
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


# ==================================================================================
# Signals quality report
# ==================================================================================


def generate_summary(workflow, filename):
    # Get task info
    task_name = [task for task in filename.split("_") if "task" in task]
    run_number = [run for run in filename.split("_") if "run" in run]
    sub = [sub for sub in filename.split("_") if "sub" in sub][0]
    ses = [ses for ses in filename.split("_") if "ses" in ses][0]

    if len(run_number) == 0:
        run_number = task_name

    # Add meta data
    html_report = f"""
    <h1>Summary</h1>
    <ul>
        <li>Subject ID : {sub.split("-")[1]}</li>
        <li>Session ID : {ses.split("-")[1]}</li>
        <li>Task : {task_name[0].split("-")[1]} (run {run_number[0].split("-")[1]})</li>
    """
    # Add info about recorded modalities
    for modality in workflow:
        if modality != "trigger":
            processing_strategy = workflow[modality]["preprocessing_strategy"]
            html_report += f"""
                <li>{modality} (channel - {workflow[modality]["Channel"]})</li>
                    <ul>
                        <li style="color:rgb(80,80,80);">
                            Description: {workflow[modality]["Description"]}
                        </li>
                        <li style="color:rgb(80,80,80);">
                            Units: {workflow[modality]["Units"]}
                        </li>
                        <li style="color:rgb(80,80,80);">
                            Preprocessing strategy: {processing_strategy}
                        </li>
                    </ul>
                </li>
                """

    html_report += "</ul>"

    return html_report


def generate_report(workflow, summary, data, info, derivatives, filename, window=False):
    # Generate the report in HTML format
    html_report = """
    <!DOCTYPE html>
    <html>
    <head>
    <style>
        table {
        font-family: Arial, sans-serif;
        border-collapse: collapse;
        width: 100%;
        }
        td, th {
        border: 1px solid #dddddd;
        text-align: left;
        padding: 8px;
        }
        th {
        background-color: #dddddd;
        }
    </style>
    <script src="https://cdn.bokeh.org/bokeh/release/bokeh-2.3.3.min.js"
        crossorigin="anonymous"></script>
    </head>
    <body>
    """
    html_report += generate_summary(workflow, filename)
    for k in summary.keys():
        html_report += f"""
        <h1>{k} Signal</h1>
        """
        if window:
            if "Overview" in list(summary[k].keys()):
                dict_overview = summary[k]["Overview"]
                dict_window = summary[k]
                del dict_window["Overview"]
                # Table for overview metrics
                html_report += "<h2>Overview</h2>"
                # Create table headers
                headers = list(dict_overview.keys())
                header_row = "<tr>{}</tr>".format(
                    "".join("<th>{}</th>".format(header) for header in headers)
                )
                # Create table row
                row = "".join(
                    "<td>{}</td>".format(dict_overview.get(col)) for col in headers
                )

                # Merge headers and rows table
                table = f"<table>{header_row + row}</table>"
                html_report += f"<table>{table}</table>"
            else:
                dict_window = summary[k]
            # Table for window-by-window metrics
            html_report += "<h2>Windows</h2>"
            headers = ["Window"] + list(dict_window[list(dict_window.keys())[0]].keys())
            # Create table headers
            header_row = "<tr>{}</tr>".format(
                "".join("<th>{}</th>".format(header) for header in headers)
            )
            # Create table rows
            rows = []
            for w, values in dict_window.items():
                row = "<tr><td>{}</td>{}</tr>".format(
                    w,
                    "".join("<td>{}</td>".format(values.get(col)) for col in headers[1:]),
                )
                rows.append(row)
            # Merge headers and rows tables
            table = "<table>{}</table>".format(header_row + "".join(rows))
            html_report += f"<table>{table}</table>"

        else:
            # Table for overview metrics
            html_report += "<h2>Overview</h2>"
            # Create table headers
            headers = list(summary[k].keys())
            header_row = "<tr>{}</tr>".format(
                "".join("<th>{}</th>".format(header) for header in headers)
            )
            # Create table rows
            row = "".join("<td>{}</td>".format(summary[k].get(col)) for col in headers)

            # Merge headers and rows table
            table = f"<table>{header_row + row}</table>"
            html_report += f"<table>{table}</table>"
        # Add interactive plot
        html_report += "<h2>Plot</h2>"
        # Generate interactive figure
        script, div = generate_plot(data, info, k)
        html_report += f"{script}"
        html_report += f"<div>{div}</div>"

    # Complete the HTML report
    html_report += """
    </body>
    </html>
    """

    # Save the HTML report to a file
    print("Saving html report")
    with open(Path(derivatives / filename).with_suffix(".html"), "w") as file:
        file.write(html_report)
        file.close()
