# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""Physiological data quality assessment"""

import glob
import json
import os
import traceback
from pathlib import Path

import click
import pandas as pd
from time_sqi import (
    sqi_cardiac,
    sqi_cardiac_overview,
    sqi_eda,
    sqi_eda_overview,
    sqi_rsp,
)
from utils import load_json
from visu.plot_signals import generate_plot


@click.command()
@click.argument("source", type=str)
@click.argument("derivatives")
@click.argument("sub", type=str)
@click.argument("ses", type=str)
@click.argument("sliding", type=str)
def neuromod_bio_sqi(source, derivatives, sub, ses, sliding={"duration": 60, "step": 10}):
    """
    Run processing QC-ing pipeline on specified biosignals, and generate an html report
    containing the quality metrics.

    Parameters
    ----------
    source : str
        The directory containing the runs unfiltered.
    derivatives: str
        The directory containing the processed runs.
    sub : str
        The id of the subject.
    ses : str
        The id of the session.
    sliding : dict
        Dictionary containing the `duration` (in seconds) of the windows in which to
        calculate the SQI. If `step` is not specified, is equal to 0 or to None,
        the fixed windows approach is used for QC-ing. If `step` is specified, a
        sliding window approach is used for QC-ing, and the value corresponds to the
        step between each window (in seconds). If `duration` is None, the SQI are
        computed on the whole run.
        Default to {'duration': 60, 'step': 10}.
    """
    sliding = json.loads(sliding)

    filenames = glob.glob(os.path.join(derivatives, sub, ses, "*_physio*"))
    filenames_signal = [
        f
        for f in filenames
        if "".join(Path(f).suffixes) == ".tsv.gz" and "noseq" not in f
    ]
    filenames_signal.sort()

    for idx, f in enumerate(filenames_signal):
        filename = os.path.basename(f).split(".")[0]
        info = load_json(os.path.join(derivatives, sub, ses, filename + ".json"))

        signal = pd.read_csv(os.path.join(derivatives, sub, ses, f), sep="\t")
        summary, summary_ppg, summary_ecg, summary_eda, summary_rsp = (
            {},
            {},
            {},
            {},
            {},
        )
        print(f"---QCing on {f.split('/')[-1]}---")
        # Compute metrics on the unsegmented signal
        if not sliding:
            for modality in list(info.keys()):
                print(f"***Computing quality metrics for {modality} signal***")
                try:
                    if modality in ["PPG", "ECG"]:
                        summary[modality] = sqi_cardiac_overview(
                            info[modality], data_type=modality
                        )
                        summary[modality].update(
                            sqi_cardiac(
                                signal[f"{modality}_Clean"],
                                info[modality],
                                data_type=modality,
                                sampling_rate=info[modality]["sampling_rate"],
                            )
                        )
                    if modality == "EDA":
                        summary["EDA"] = sqi_eda(
                            signal,
                            info["EDA"],
                            sampling_rate=info["EDA"]["sampling_rate"],
                        )
                        summary["EDA"] = sqi_eda_overview(info["EDA"])
                    if modality == "RSP":
                        summary["RSP"] = sqi_rsp(
                            signal, sampling_rate=info["RSP"]["sampling_rate"]
                        )
                except Exception:
                    print(f"Not able to compute {modality}")
                    traceback.print_exc()

            # Generate report
            print("***Generating report***")
            generate_report(
                summary,
                source,
                derivatives,
                sub,
                ses,
                f"{filename.split('/')[-1]}",
            )

        # Compute metrics on signal segmented in fixed windows
        elif (
            sliding.get("step") == 0
            or not sliding.get("step")
            or sliding.get("step") is None
        ):
            print("Fixed window")
            for modality in list(info.keys()):
                try:
                    window_samples = int(
                        sliding["duration"] * info[modality]["sampling_rate"]
                    )
                    num_windows = int(len(signal) / window_samples)

                    # Descriptive metrics on the windows
                    print(f"***Computing quality metrics for {modality} signal***")
                    for i in range(num_windows):
                        start = int(i * window_samples)
                        start_idx = int(start / info[modality]["sampling_rate"])
                        end = int(start + window_samples)
                        end_idx = int(end / info[modality]["sampling_rate"])
                        window = signal.iloc[start:end]
                        if modality == "PPG":
                            summary_ppg[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                                window["PPG_Clean"],
                                info["PPG"],
                                data_type="PPG",
                                sampling_rate=info["PPG"]["sampling_rate"],
                                window=[start, end],
                            )
                        elif modality == "ECG":
                            summary_ecg[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                                window["ECG_Clean"],
                                info["ECG"],
                                data_type="ECG",
                                sampling_rate=info["ECG"]["sampling_rate"],
                                window=[start, end],
                            )
                        elif modality == "EDA":
                            summary_eda[f"{start_idx}-{end_idx}"] = sqi_eda(
                                window,
                                info["EDA"],
                                sampling_rate=info["EDA"]["sampling_rate"],
                                window=[start, end],
                            )
                        elif modality == "RSP":
                            summary_rsp[f"{start_idx}-{end_idx}"] = sqi_rsp(
                                window,
                                sampling_rate=info["RSP"]["sampling_rate"],
                            )
                    # Descriptive metrics on the overall run
                    if modality == "PPG":
                        summary["PPG"] = {
                            "Overview": sqi_cardiac_overview(info["PPG"], data_type="PPG")
                        }
                        summary["PPG"].update(summary_ppg)
                    elif modality == "ECG":
                        summary["ECG"] = {
                            "Overview": sqi_cardiac_overview(info["ECG"], data_type="ECG")
                        }
                        summary["ECG"].update(summary_ecg)
                    elif modality == "EDA":
                        summary["EDA"] = summary_eda
                        summary["EDA"].update({"Overview": sqi_eda_overview(info["EDA"])})
                    elif modality == "RSP":
                        summary["RSP"] = summary_rsp

                except Exception:
                    print(f"Not able to compute {modality}")
                    traceback.print_exc()

            # Generate report
            print("***Generating report***")
            generate_report(
                summary,
                source,
                derivatives,
                sub,
                ses,
                f"{filename.split('/')[-1]}",
                window=True,
            )

        # Compute metrics using a sliding window approach
        else:
            print("Sliding window")
            for modality in list(info.keys()):
                try:
                    window_samples = int(
                        sliding["duration"] * info[modality]["sampling_rate"]
                    )
                    step_samples = int(sliding["step"] * info[modality]["sampling_rate"])
                    num_windows = int((len(signal) - window_samples) / step_samples) + 1

                    print(f"***Computing quality metrics for {modality} signal***")
                    for i in range(num_windows):
                        start = int(i * step_samples)
                        start_idx = int(start / info[modality]["sampling_rate"])
                        end = int(start + window_samples)
                        end_idx = int(end / info[modality]["sampling_rate"])
                        window = signal.iloc[start:end]
                        if modality == "PPG":
                            summary_ppg[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                                window["PPG_Clean"],
                                info["PPG"],
                                data_type="PPG",
                                sampling_rate=info["PPG"]["sampling_rate"],
                                window=[start, end],
                            )
                        elif modality == "ECG":
                            summary_ecg[f"{start_idx}-{end_idx}"] = sqi_cardiac(
                                window["ECG_Clean"],
                                info["ECG"],
                                data_type="ECG",
                                sampling_rate=info["ECG"]["sampling_rate"],
                                window=[start, end],
                            )
                        elif modality == "EDA":
                            summary_eda[f"{start_idx}-{end_idx}"] = sqi_eda(
                                window,
                                info["EDA"],
                                sampling_rate=info["EDA"]["sampling_rate"],
                                window=[start, end],
                            )
                        elif modality == "RSP":
                            summary_rsp[f"{start_idx}-{end_idx}"] = sqi_rsp(
                                window,
                                sampling_rate=info["RSP"]["sampling_rate"],
                            )

                    # Descriptive metrics on the overall run
                    if modality == "PPG":
                        summary["PPG"] = {
                            "Overview": sqi_cardiac_overview(info["PPG"], data_type="PPG")
                        }
                        summary["PPG"].update(summary_ppg)
                    elif modality == "ECG":
                        summary["ECG"] = {
                            "Overview": sqi_cardiac_overview(info["ECG"], data_type="ECG")
                        }
                        summary["ECG"].update(summary_ecg)
                    elif modality == "EDA":
                        summary["EDA"] = summary_eda
                        summary["EDA"].update({"Overview": sqi_eda_overview(info["EDA"])})
                    elif modality == "RSP":
                        summary["RSP"] = summary_rsp
                except Exception:
                    print(f"Not able to compute {modality}")
                    traceback.print_exc()

            # Generate report
            print("***Generating report***")
            generate_report(
                summary,
                source,
                derivatives,
                sub,
                ses,
                f"{filename.split('/')[-1]}",
                window=True,
            )


# ==================================================================================
# Signals quality report
# ==================================================================================


def generate_summary(source, sub, ses, filename):
    # Get task info
    task_name = [task for task in filename.split("_") if "task" in task]
    run_number = [run for run in filename.split("_") if "run" in run]
    if len(run_number) == 0:
        run_number = task_name
    # Get session meta data
    meta_info = [f for f in os.listdir(os.path.join(source, sub)) if ".json" in f]
    if len(meta_info) == 1:
        meta_info = load_json(os.path.join(source, sub, meta_info[0]))
    else:
        meta_info = {}
    source_info = load_json(os.path.join(source, sub, ses, filename + ".json"))
    # Add meta data
    html_report = f"""
    <h1>Summary</h1>
    <ul>
        <li>Subject ID : {sub.split("-")[1]}</li>
        <li>Session ID : {ses.split("-")[1]}</li>
        <li>Task : {task_name[0].split("-")[1]} (run {run_number[0].split("-")[1]})</li>
    """
    # Add info about recorded modalities
    if bool(meta_info):
        ch = source_info["Columns"]
        ch.remove("time")
        ch_names = meta_info[ses]["ch_names"]
        del ch_names[ch.index("TTL")]
        del ch[ch.index("TTL")]
        for c, ch_name in zip(ch, ch_names):
            html_report += f"""
                <li>{c} (channel - {ch_name})
                    <ul>
                        <li style="color:rgb(80,80,80);">
                            Sampling rate: {source_info["SamplingFrequency"]} Hz
                        </li>
                    </ul>
                </li>
                """
        html_report += "</ul>"

    return html_report


def generate_report(summary, source, derivatives, sub, ses, filename, window=False):
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
    html_report += generate_summary(source, sub, ses, filename)
    filename_json = os.path.join(derivatives, sub, ses, filename + ".json")
    for k in summary.keys():
        html_report += f"""
        <h1>{k} Signal</h1>
        <input type="hidden" id=filenameJson value={filename_json}>
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
        script, div = generate_plot(derivatives, sub, ses, filename, k)
        html_report += f"{script}"
        html_report += f"<div>{div}</div>"

    # Complete the HTML report
    html_report += """
    </body>
    </html>
    """

    # Save the HTML report to a file
    print("Saving html report")
    with open(
        os.path.join(derivatives, sub, ses, f"{filename}_quality.html"), "w"
    ) as file:
        file.write(html_report)
        file.close()


if __name__ == "__main__":
    neuromod_bio_sqi()
