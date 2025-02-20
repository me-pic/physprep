# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""
Physiological signals plotting functions.
"""
import json
import os
import traceback
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from bokeh.embed import components
from bokeh.io import output_notebook
from bokeh.layouts import column, gridplot
from bokeh.models import BoxAnnotation, ColumnDataSource, RangeTool
from bokeh.plotting import figure, output_file, save
from pandas.core.indexes.datetimes import DatetimeIndex
from systole.plots import plot_rr

from physprep.utils import load_json

output_notebook()


def load_data(outdir, sub, ses):
    path = os.path.join(outdir, sub, ses)
    files = [
        f.split(".")[0] for f in os.listdir(path) if "tsv.gz" in f and "noseq" not in f
    ]
    files.sort()

    data, data_noseq, info = [], [], []
    for f in files:
        print(f)
        data.append(pd.read_csv(os.path.join(outdir, sub, ses, f + ".tsv.gz"), sep="\t"))
        data_noseq.append(
            pd.read_csv(os.path.join(outdir, sub, ses, f + "_noseq.tsv.gz"), sep="\t")
        )
        info.append(load_json(os.path.join(outdir, sub, ses, f + ".json")))

    return data, data_noseq, files, info


def plot_scr(
    signal: np.ndarray = None,
    peaks: np.ndarray = None,
    onsets: np.ndarray = None,
    sfreq: int = None,
    figsize: int = 300,
) -> figure:
    """
    Visualization of phasic component of the EDA signal.
    """
    time = pd.to_datetime(np.arange(0, len(signal)) / sfreq, unit="s")

    source = ColumnDataSource(data={"time": time, "signal": signal})

    p1 = figure(
        title="Skin conductance response",
        sizing_mode="stretch_width",
        height=figsize,
        x_axis_label="Time",
        x_axis_type="datetime",
        y_axis_label="Amplitude (uSiemens)",
        output_backend="webgl",
        tools="pan,wheel_zoom,box_zoom,box_select,reset,save",
        x_range=(time[0], time[-1]),
    )

    p1.line(
        "time",
        "signal",
        source=source,
        legend_label="SCR",
        line_color="#3783a9",
    )

    p1.circle(
        x=time[peaks],
        y=signal[peaks],
        size=5,
        legend_label="SCR - Peaks",
        fill_color="#c56c5e",
        line_color="black",
    )

    p1.circle(
        x=time[onsets],
        y=signal[onsets],
        size=5,
        legend_label="SCR - Onsets",
        fill_color="#6c0073",
        line_color="black",
    )

    cols = (p1,)

    if len(cols) > 1:
        return column(*cols, sizing_mode="stretch_width")
    else:
        return cols[0]


def plot_raw(
    signal: np.ndarray,
    eda_scr: np.ndarray = None,
    eda_scl: np.ndarray = None,
    time: DatetimeIndex = None,
    peaks: np.ndarray = None,
    onsets: np.ndarray = None,
    detector: str = None,
    sfreq: int = None,
    modality: str = "ppg",
    title: str = "Figure",
    show_heart_rate: bool = True,
    show_artefacts: bool = False,
    bad_segments: Optional[List[Tuple[int, int]]] = None,
    decim: int = 10,
    slider: bool = True,
    figsize: int = 300,
    events_params: Optional[Dict] = None,
    **kwargs,
) -> figure:
    """Visualization of PPG or ECG signal with systolic peaks/R wave detection.

    The instantaneous heart rate can be derived in a second row.

    Parameters
    ----------
    time :
        The time index.
    signal :
        The physiological signal (1d numpy array).
    eda_scr :
        The skin conductance response component associated with EDA signal
        (1d numpy array).
    eda_scl :
        The skin conductance level component associate with the EDA signal
        (1d numpy array).
    peaks :
        The peaks or R wave detection (1d boolean array).
    modality :
        The recording modality. Can be `"cardiac_ppg"`, `"cardiac_ecg"`, "electrodermal"
        or `"respiratory"`.
    show_heart_rate :
        If `True`, create a second row and plot the instantanesou heart rate
        derived from the physiological signal
        (calls :py:func:`systole.plots.plot_rr` internally). Defaults to `False`.
    show_artefacts :
        If `True`, the function will call
       :py:func:`systole.detection.rr_artefacts` to detect outliers intervalin the time
        series and outline them using different colors.
    bad_segments :
        Mark some portion of the recording as bad. Grey areas are displayed on the top
        of the signal to help visualization (this is not correcting or transforming the
        post-processed signals). Should be a list of tuples shuch as (start_idx,
        end_idx) for each segment.
    decim :
        Factor by which to subsample the raw signal. Selects every Nth sample (where N
        is the value passed to decim). Default set to `10` (considering that the input
        signal has a sampling frequency of 1000 Hz) to save memory.
    slider :
        If `True`, add a slider to zoom in/out in the signal (only working with
        bokeh backend).
    figsize :
        Figure heights. Default is `300`.
    events_params :
        (Optional) Additional parameters that will be passed to
       :py:func:`systole.plots.plot_events` and plot the events timing in the background.
    kwargs:
        Other keyword arguments passed to the function but unused by the Bokeh backend.

    Returns
    -------
    raw :
        The bokeh figure containing the plot.

    """

    eda_strings = ["eda", "gsr", "electrodermal", "electrodermal activity"]

    if figsize is None:
        figsize = 300

    # if time is None:
    time = pd.to_datetime(np.arange(0, len(signal)) / sfreq, unit="s")

    if peaks is not None:
        source = ColumnDataSource(
            data={
                "time": time[::decim],
                "signal": signal[::decim],
                "peaks": peaks[::decim],
            }
        )
    else:
        source = ColumnDataSource(data={"time": time[::decim], "signal": signal[::decim]})

    if modality.lower() in ["cardiac_ppg", "ppg"]:
        # title = "PPG recording"
        ylabel = "PPG level (a.u.)"
        peaks_label = "Systolic peaks"
        signal_label = "PPG signal"
    elif modality.lower() in ["cardiac_ecg", "ecg"]:
        # title = "ECG recording"
        ylabel = "ECG (mV)"
        peaks_label = "R wave"
        signal_label = "ECG signal"
    elif modality.lower() in ["respiratory", "rsp", "resp"]:
        # title = "Respiration"
        ylabel = "Respiratory signal"
        peaks_label = "End of inspiration"
        signal_label = "Respiratory signal"
    elif modality in ["electrodermal", "eda", "gsr"]:
        # title = "EDA recording"
        ylabel = "EDA (uSiemens)"
        peaks_label = "SCR Peaks"
        signal_label = "EDA signal"
    else:
        ylabel=""
        peaks_label=""
        signal_label=""

    ############
    # Raw plot #
    ############

    raw = figure(
        title=title,
        x_axis_type="datetime",
        sizing_mode="stretch_width",
        height=figsize,
        x_axis_label="Time",
        y_axis_label=ylabel,
        output_backend="webgl",
        x_range=(time[0], time[-1]),
    )

    raw.line(
        "time",
        "signal",
        source=source,
        legend_label=signal_label,
        line_color="#a9373b",
    )

    if eda_scl is not None:
        raw.line(
            "time",
            "signal",
            source=ColumnDataSource(
                data={"time": time[::decim], "signal": eda_scl[::decim]}
            ),
            legend_label="Skin Conductance Level",
            line_color="#37a98b",
        )

    if peaks is not None and modality not in eda_strings:
        raw.circle(
            x=time[peaks],
            y=signal[peaks],
            size=5,
            legend_label=peaks_label,
            fill_color="lightgrey",
            line_color="grey",
        )

    raw.legend.title = "Signal"

    cols = (raw,)

    # Highlight bad segments if provided
    if bad_segments is not None:
        for bads in bad_segments:
            # Plot time range
            event_range = BoxAnnotation(
                left=time[bads[0]],
                right=time[bads[1]],
                fill_alpha=0.2,
                fill_color="grey",
            )
            event_range.level = "underlay"
            raw.add_layout(event_range)

    # Instantaneous heart rate
    ##########################
    if show_heart_rate is True:
        instantaneous_hr = plot_rr(
            peaks.astype(bool),
            input_type="peaks",
            backend="bokeh",
            figsize=figsize,
            slider=False,
            line=True,
            show_artefacts=show_artefacts,
            events_params=events_params,
        )
        instantaneous_hr.x_range = raw.x_range

        cols += (instantaneous_hr,)  # type: ignore

    if peaks is not None and modality in eda_strings:
        scr = plot_scr(eda_scr, peaks.astype(bool), onsets.astype(bool), sfreq)
        scr.x_range = raw.x_range

        cols += (scr,)

    if slider is True:
        select = figure(
            title="Select the time window",
            y_range=raw.y_range,
            y_axis_type=None,
            height=int(figsize * 0.5),
            x_axis_type="datetime",
            tools="",
            toolbar_location=None,
            background_fill_color="#efefef",
        )

        range_tool = RangeTool(x_range=raw.x_range)
        range_tool.overlay.fill_color = "navy"
        range_tool.overlay.fill_alpha = 0.2

        select.line("time", "signal", source=source)
        select.ygrid.grid_line_color = None
        select.add_tools(range_tool)

        cols += (select,)  # type: ignore

    if len(cols) > 1:
        return column(*cols, sizing_mode="stretch_width")
    else:
        return cols[0]


def generate_plot(data, features, metadata, modality):
    """
    Parameters
    ----------
    data : dict
        Dictionary containing the timeserie to plot.
    info : dict
        Dictionary containing the features of the timeserie.
    modality : list
        A list containing the biosignal modalities to plot.
        The options include "ECG", "PPG", "EDA", and "RSP".
    """
    print(f"Plotting {modality} signal: begin")
    if isinstance(metadata['SamplingFrequency'], list):
        sampling_frequency = metadata['SamplingFrequency'][metadata['Columns'].index(f'{modality}_clean')]
    else:
        sampling_frequency = metadata['SamplingFrequency']

    try:
        # Plot cleaned signal during MRI sequence
        if modality.lower() in ["respiratory", "rsp", "resp"]:
            peaks = [1 if idx in features[modality]["inhale_max"] else 0 for idx, i in enumerate([0]*len(data[f"{modality}_clean"]))]

            figure_signal = plot_raw(
                signal=data[f"{modality}_clean"],
                peaks=np.array(peaks).astype(bool),
                sfreq=sampling_frequency,
                modality=modality,
                title=f"{modality} : Scanner on - Clean",
                show_heart_rate=False,
                show_artefacts=True,
            )
        elif modality.lower() in ["electrodermal", "eda", "gsr"]:
            peaks = [1 if idx in features[modality]["scr_peak"] else 0 for idx, i in enumerate([0]*len(data[f"{modality}_clean"]))]
            onsets = [1 if idx in features[modality]["scr_onset"] else 0 for idx, i in enumerate([0]*len(data[f"{modality}_clean"]))]
            if f"{modality}_phasic" in data.keys():
                figure_signal = plot_raw(
                    signal=data[f"{modality}_clean"],
                    eda_scr=data[f"{modality}_phasic"],
                    eda_scl=data[f"{modality}_tonic"],
                    peaks=np.array(peaks).astype(bool),
                    onsets=np.array(onsets).astype(bool),
                    sfreq=sampling_frequency,
                    modality=modality,
                    title=f"{modality} : Scanner on - Clean",
                    show_heart_rate=False,
                    show_artefacts=True,
                )
            else:
                figure_signal = plot_raw(
                    signal=data[f"{modality}_clean"],
                    peaks=np.array(peaks).astype(bool),
                    onsets=np.array(onsets).astype(bool),
                    sfreq=sampling_frequency,
                    modality=modality,
                    title=f"{modality} : Scanner on - Clean",
                    show_heart_rate=False,
                    show_artefacts=True,
                )
        elif modality.lower() in ["cardiac", "ecg", "ppg"]:
            peak_key = [key for key in features[modality].keys() if 'peak_corrected' in key][0]
            peaks = [1 if idx in features[modality][peak_key] else 0 for idx, i in enumerate([0]*len(data[f"{modality}_clean"]))]

            figure_signal = plot_raw(
                signal=data[f"{modality}_clean"],
                peaks=np.array(peaks).astype(bool),
                sfreq=sampling_frequency,
                modality=modality.lower(),
                title=f"{modality} : Scanner on - Clean",
                show_heart_rate=True,
                show_artefacts=True,
            )
        print(f"Plotting {modality} signal: done")
    except Exception:
        print(f"Could not plot {modality} signal")
        traceback.print_exc()

    return components(figure_signal)


def generate_raw_filtered_plots(outdir, sub, ses, modality):
    """
    Generate interactive plots for each modality.

    Parameters
    ----------
    outdir : str
        The directory containing the processed signals.
    task : str
        The id of the task.
    sub : str
        The id of the subject.
    ses : str
        The id of the session.
    modality : list
        A list containing the biosignal modalities to plot.
        The options include "ECG", "PPG", "EDA", and "RSP".
    """
    modality = json.loads(modality)

    data, data_noseq, filenames, info = load_data(outdir, sub, ses)
    outdir = os.path.join(outdir, sub, ses)

    for i, filename in enumerate(filenames):
        figures = []
        idx = 0
        print(f"Generate plots for {filename}")
        for mod in modality:
            tmp = 0
            try:
                # Plot raw signal before MRI sequence
                print(f"Plotting Raw signal with scanner off: {mod} begin")
                figures.append(
                    plot_raw(
                        signal=data_noseq[i][mod],
                        sfreq=10000,
                        modality=mod.lower(),
                        title=f"{mod} : Scanner off - Raw",
                        show_heart_rate=False,
                        show_artefacts=True,
                    )
                )
                idx += 1
                tmp += 1
                print(f"Plotting Raw signal with scanner off: {mod} done")

                # Plot raw signal during MRI sequence
                print(f"Plotting Raw signal with scanner on: {mod} begin")
                figures.append(
                    plot_raw(
                        signal=data[i][f"{mod}_raw"],
                        sfreq=info[i][mod]["Sampling"],
                        modality=mod.lower(),
                        title=f"{mod} : Scanner on - Raw",
                        show_heart_rate=False,
                        show_artefacts=True,
                    )
                )
                idx += 1
                tmp += 1
                print(f"Plotting Raw signal with scanner on: {mod} done")

                # Plot cleaned signal during MRI sequence
                if mod == "respiratory":
                    print(f"Plotting Clean signal with scanner on: {mod} begin")
                    figures.append(
                        plot_raw(
                            signal=data[i][f"{mod}_clean"],
                            peaks=data[i][f"{mod}_peaks"].astype(bool),
                            sfreq=info[i][mod]["SamplingFrequency"],
                            modality="resp",
                            title=f"{mod} : Scanner on - Clean",
                            show_heart_rate=False,
                            show_artefacts=True,
                        )
                    )
                    idx += 1
                    tmp += 1
                    print(f"Plotting Clean signal with scanner on: {mod} done")
                elif mod == "electrodermal":
                    print(f"Plotting Clean signal with scanner on: {mod} begin")
                    figures.append(
                        plot_raw(
                            signal=data[i]["eda_clean"],
                            eda_scr=data[i]["eda_phasic"],
                            eda_scl=data[i]["eda_tonic"],
                            peaks=data[i]["scr_peaks"],
                            onsets=data[i]["scr_onsets"],
                            sfreq=info[i][mod]["SamplingFrequency"],
                            modality="eda",
                            title=f"{mod} : Scanner on - Clean",
                            show_heart_rate=False,
                            show_artefacts=True,
                        )
                    )
                    idx += 1
                    tmp += 1
                    print(f"Plotting Clean signal with scanner on: {mod} done")
                elif mod in ["cardiac_ecg", "cardiac_ppg"]:
                    print(f"Plotting Clean signal with scanner on: {mod} begin")
                    figures.append(
                        plot_raw(
                            signal=data[i][f"{mod}_clean"],
                            peaks=data[i][f"{mod}_peaks_nk"].astype(bool),
                            sfreq=info[i][mod]["SamplingFrequency"],
                            modality=mod.lower(),
                            title=f"{mod} : Scanner on - Clean",
                            show_heart_rate=True,
                            show_artefacts=True,
                        )
                    )
                    idx += 1
                    tmp += 1
                    print(f"Plotting Clean signal with scanner on: {mod} done")
            except Exception:
                idx -= tmp
                del figures[-tmp:]
                print(f"Could not plot {mod} signal")
                traceback.print_exc()

        # Organize figures in grid containing 3 columns
        # (Raw: scanner OFF; Raw: scanner ON; Clean: scanner ON)
        layout = gridplot(
            children=[[figures[i] for i in range(j, j + 3)] for j in range(0, idx - 1, 3)]
        )

        # Save generate figure in html
        output_file(outdir + f"/{filename}_plot.html")
        save(layout)


if __name__ == "__main__":
    generate_plot()
