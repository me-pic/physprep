# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""Physiological data quality assessment"""

import traceback
from pathlib import Path

from physprep.visu.plot_signals import generate_plot


def generate_summary(workflow, bids_entities):
    # Get task info
    task = bids_entities['task']
    run = bids_entities['run'] if 'run' in bids_entities.keys() else '-'
    sub = bids_entities['subject']
    ses = bids_entities['session']

    # Add meta data
    html_report = f"""
    <h1>Summary</h1>
    <ul>
        <li>Subject ID : {sub}</li>
        <li>Session ID : {ses}</li>
        <li>Task : {task}</li>
        <li>Run : {run}</li>
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


def generate_report(workflow, summary, data, info, metadata_derivatives,derivatives, bids_entities, window=False):
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
    html_report += generate_summary(workflow, bids_entities)
    for k in summary.keys():
        html_report += f"""
        <h1>{k.capitalize()} signal</h1>
        """
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

        # Add interactive plot
        html_report += "<h2>Plot</h2>"
        # Generate interactive figure
        script, div = generate_plot(data, info, metadata_derivatives, k)
        html_report += f"{script}"
        html_report += f"<div>{div}</div>"

    # Complete the HTML report
    html_report += """
    </body>
    </html>
    """
    return html_report
