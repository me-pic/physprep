# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning

"""another util for neuromod phys data conversion."""

import json
import logging
import math
import os

import numpy as np
import pandas as pd
import pprintpp
from bioread import read_file
from neurokit2 import read_acqknowledge

from physprep.prepare.list_sub import list_sub
from physprep.utils import _check_bids_validity, load_json

LGR = logging.getLogger(__name__)


def order_channels(acq_channels, metadata_physio):
    """
    Order channels in the acq file according to the metadata_physio file.

    Parameters
    ----------
    acq_channels : list
        List of channels in the acq file.
    metadata_physio : dict
        Dictionary containing the metadata_physio file.
    """
    ch_names = []
    chsel = []
    for idx, channel in enumerate(acq_channels):
        found = False
        for key in metadata_physio:
            if key != "concurrentWith":
                if channel == metadata_physio[key]["Channel"]:
                    ch_names.append(key)
                    chsel.append(idx + 1)
                    found = True
        if not found:
            ch_names.append(channel)

    if len(chsel) == 0:
        raise ValueError(
            "No correspondence between channels in the acq file and "
            "channels defined in workflow configuration file."
        )

    return ch_names, chsel


def volume_counter(root, sub, metadata_physio, ses=None, tr=1.49, trigger_ch="TTL"):
    """
    Volume counting for each run in a session.

    Parameters
    ----------
    root : str
        Directory containing the biopac data.
        Example: "/home/user/dataset/sourcedata/physio".
    sub : str
        Name of path for a specific subject. Example: "sub-01".
    metadata_physio : dict
        Dictionary containing the metadata_physio file.
    ses : str
        Name of path for a specific session (optional workflow for specific experiment).
        Default to none.
    tr : float
            Value of the TR used in the MRI sequence.
        Default to 1.49.
    trigger_ch : str
        Name of the trigger channel used on Acknowledge.
        Default to 'TTL'.

    Returns
    -------
    ses_runs: dict
        Each key lists the number of volumes/triggers in each run,
        including invalid volumes.
    """
    LGR = logging.getLogger(__name__)
    # Check directory
    if os.path.exists(root) is False:
        raise ValueError("Couldn't find the following directory: ", root)

    # Load metadata_physio
    if isinstance(metadata_physio, os.PathLike) or isinstance(metadata_physio, str):
        metadata_physio = load_json(metadata_physio)

    # List the files that have to be counted
    if ses == "files":
        ses = None
    dirs = list_sub(root, sub, ses)
    ses_runs = {}
    # loop iterating through files in each dict key representing session
    # returned by list_sub for this loop, exp refers to session's name,
    # avoiding confusion with ses argument
    for exp in dirs:
        LGR.info(f"counting volumes in physio file for: {exp}")
        for file in sorted(dirs[exp]):
            # reading acq
            if exp == "files":
                path_to_file = os.path.join(root, sub, file)
            else:
                path_to_file = os.path.join(root, sub, exp, file)
            bio_df, fs = read_acqknowledge(path_to_file)
            print(bio_df)
            # find the correct index of Trigger channel
            if trigger_ch in bio_df.columns:
                trigger_index = list(bio_df.columns).index(trigger_ch)
                # initialize a df with TTL values over 4 (switch either ~0 or ~5)
            else:
                trigger_index = list(bio_df.columns).index("TTL")
                # initialize a df with TTL values over 4 (switch either ~0 or ~5)
            query_df = bio_df[bio_df[bio_df.columns[trigger_index]] > 4]

            # Define session length - this list will be less
            # memory expensive to play with than dataframe
            session = list(query_df.index)

            # maximal TR - the time distance between two adjacent TTL, now
            # given by the ceiling value of the tr (but might be tweaked if
            # needed)
            tr_period = fs * math.ceil(tr)

            # Define session length and adjust with padding
            try:
                start = int(session[0])
            except IndexError:
                LGR.info(f"No trigger channel input apparently; skipping {file}")
                return "No trigger input", bio_df.columns
                continue
            end = int(session[-1])

            # initialize list of sample index to compute nb of volumes per run
            parse_list = []

            # ascertain that session is longer than 3 min
            for idx in range(1, len(session)):
                # define time diff between current successive trigger
                time_delta = session[idx] - session[idx - 1]

                # if the time diff between two trigger values over 4
                # is larger than TR, keep both indexes
                if time_delta > tr_period:
                    parse_start = int(session[idx - 1])
                    parse_end = int(session[idx])
                    # adjust the segmentation with padding
                    # parse start is end of run
                    parse_list += [(parse_start, parse_end)]
            if len(parse_list) == 0:
                runs = round((end - start) / fs / tr + 1)
                if exp not in ses_runs:
                    ses_runs[exp] = [runs]
                else:
                    ses_runs[exp].append([runs])
                continue
            # Create tuples with the given indexes
            # First block is always from first trigger to first parse
            block1 = (start, parse_list[0][0])

            # runs is a list of tuples specifying runs in the session
            runs = []
            # push the resulting tuples (run_start, run_end)
            runs.append(block1)
            for i in range(0, len(parse_list)):
                try:
                    runs.append((parse_list[i][1], parse_list[1 + i][0]))

                except IndexError:
                    runs.append((parse_list[i][1], end))

            # compute the number of trigger/volumes in the run
            for i in range(0, len(runs)):
                runs[i] = round(((runs[i][1] - runs[i][0]) / fs) / tr) + 1
            if exp not in ses_runs:
                ses_runs[exp] = [runs]
            else:
                ses_runs[exp].append(runs)

    ch_names, chsel = order_channels(bio_df.columns, metadata_physio)

    LGR.info(f"Volumes for session :\n{ses_runs}")
    return ses_runs, ch_names, chsel


def duration_counter(trigger_ch, sampling_rate, thr=5):
    """
    Duration counting for each run in a session.

    Parameters
    ----------
    trigger_ch : nd.array
        Trigger timeserie
    sampling_rate : float
        Sampling_rate of `trigger_ch`
    thr : int
        Threshold to using for trigger counting

    Returns
    -------
    runs : List
        Duration of each run based
    """
    idx = np.where(trigger_ch >= thr)[0].tolist()
    runs = []
    start, end = idx[0], None
    for i, tp in enumerate(idx):
        # If index not the first or last in the list
        if tp != idx[-1] and tp != idx[0]:
            # If index not continuous
            if tp != idx[i - 1] + 1:
                end = np.where(trigger_ch[start:tp] < thr)[0][0] + start
                if end is not None and start is not None:
                    runs.append(
                        float(np.floor(100 * (int(end) - start) / sampling_rate) / 100)
                    )
                    start = tp
                    end = None
        elif tp != idx[0]:
            # Make sure we include the last take
            runs.append(float(np.floor(100 * (tp - start) / sampling_rate) / 100))

    return runs


def get_info(root, sub, workflow, ses=None, count_vol=False):
    """
    Get all volumes taken for a sub.
    `get_info` pushes the info necessary to execute the phys2bids multi-run
    workflow to a dictionary.

    Arguments
    ---------
    root : str or pathlib.Path or BIDSLayout
        Root directory of the BIDS dataset containing the data.
    sub : str
        Name of path for a specific subject. Example: "sub-01".
    workflow : path
        Path to workflow strategy config file.
    ses : str
        Name of path for a specific session. Example: "ses-001".
    count_vol : bool
        Specify if you want to count triggers in physio file.
        Default to False.

    Returns
    -------
    ses_runs_vols : dict
        Number of processed runs, number of expected runs, number of
        triggers/volumes per run, sourcedata file location.
    """
    LGR = logging.getLogger(__name__)
    # Get BIDSLayout from root
    layout = _check_bids_validity(root)
    workflow = load_json(workflow)

    if "concurrentWith" in workflow.keys():
        modality = workflow["concurrentWith"]["dataType"]
        trigChannel = workflow["concurrentWith"]["trigChannel"]
        trigThresh = workflow["concurrentWith"]["trigThresh"]
    else:
        print("No `concurrentWith` key found in `workflow`")
        modality = "physio"

    # list matches for a whole subject's dir
    ses_runs_matches = list_sub(
        os.path.join(root, "code"),
        sub,
        ses=ses,
        ext=".tsv",
        show=True,
    )

    # go to fmri matches and get entries for each run of a session
    nb_expected_runs = {}

    # If there is a tsv file matching the acq file and the nii.gz files in root
    ses_info = list_sub(
        os.path.join(layout.root, f"sourcedata/{modality}/"), sub, ses, ext=".acq"
    )

    # iterate through sessions and get _matches.tsv with list_sub dict
    for exp in sorted(ses_runs_matches):
        LGR.info(exp)

        if ses_info[exp] == []:
            LGR.info("No acq file found for this session")
            continue
        elif exp == "files":
            LGR.info("No SES IDs")
            path_to_source = os.path.join(layout.root, f"sourcedata/{modality}/", sub)
        else:
            path_to_source = os.path.join(
                layout.root, f"sourcedata/{modality}/", sub, exp
            )

        if len(ses_runs_matches[exp]) == 0:
            LGR.info("No tsv file matching concurrent recordings found for this session")
            continue
        elif len(ses_runs_matches[exp]) > 1:
            LGR.info(
                "Multiple tsv files matching concurrent recordings found for this session"
            )
            LGR.info(
                f"Will not proceed with {exp}. Please check the tsv files in "
                f"{os.path.join(root, 'sourcedata/physio', sub)}"
            )
            continue

        # Load tsv file outputted from match_acq_bids
        matches = pd.read_csv(
            os.path.join(root, "code", sub, ses, ses_runs_matches[exp][0]), sep="\t"
        )
        # Get all the files related to concurrent recordings
        list_matches = matches.iloc[:, 0].tolist()
        list_matches = [layout.get_file(match) for match in list_matches]

        # initialize a counter and a dictionary
        nb_expected_volumes_run = {}
        tr, concurrent_file = [], []

        # iterate through matches (concurrent recordings)
        for idx, match in enumerate(list_matches):
            entities = match.get_entities()
            metadata = match.get_metadata()
            concurrent_file.append(match.filename)

            if entities["datatype"] == "eeg":
                # For physio runs that are delimited by a continuous trigger,
                # to be able to use phys2bids, we will consider the triggers
                # to be really long tr. So the number of expected timepoints will be
                # one per run, and the tr duration will be equal to the duration
                # of each run
                from mne.io import read_raw_edf

                data_eeg = read_raw_edf(match)
                tr.append(
                    duration_counter(
                        data_eeg[trigChannel][0][0],
                        data_eeg.info["sfreq"],
                        thr=trigThresh,
                    )
                )
                nb_expected_volumes_run[f"{idx+1:02d}"] = [1] * len(tr[idx])
            else:
                # we want to have the TR in a _bold.json to later use it in the
                # volume_counter function
                # Check if metadata are there
                if not metadata:
                    LGR.info(f"No metadata to match : {exp}")
                    continue

                tr.append(metadata["RepetitionTime"])
                # we want to GET THE NB OF VOLUMES in the _bold.json of a given run
                try:
                    nb_expected_volumes_run[f"{idx+1:02d}"] = metadata["dcmmeta_shape"][
                        -1
                    ]
                except KeyError:
                    LGR.info("Cannot access Nifti BIDS metadata")

        # print the thing to show progress
        LGR.info(f"BIDS metadata; number of volumes per run:\n{nb_expected_volumes_run}")
        # push all info in run in dict
        nb_expected_runs[exp] = {}
        # the nb of expected volumes in each run of the session (embedded dict)
        nb_expected_runs[exp] = nb_expected_volumes_run
        nb_expected_runs[exp]["expected_runs"] = len(matches)
        # nb_expected_runs[exp]['processed_runs'] = idx  # counter is used here
        nb_expected_runs[exp]["task"] = entities["task"]
        nb_expected_runs[exp]["tr"] = tr
        nb_expected_runs[exp]["concurrent_with"] = entities["datatype"]
        nb_expected_runs[exp]["concurrent_file"] = concurrent_file

        # save the name
        name = ses_info[exp]
        if name:
            name.reverse()
            nb_expected_runs[exp]["in_file"] = name

        if count_vol:
            ls_run_dict, ls_ch_names, ls_chsel = [], [], []
            # check if biopac file exist, notify the user that we won't
            # count volumes
            try:
                # do not count the triggers in phys file if no physfile
                for idx_n, n in enumerate(name):
                    run_dict = {}
                    if os.path.isfile(os.path.join(path_to_source, n)) is False:
                        tmp_filename = os.path.join(
                            layout.root, "sourcedata", modality, sub, exp, n
                        )
                        LGR.info(
                            "cannot find session directory for "
                            f"sourcedata :\n{tmp_filename}"
                        )
                    else:
                        # count the triggers in physfile otherwise
                        try:
                            if entities["datatype"] == "eeg":
                                bio_data = read_file(
                                    os.path.join(
                                        layout.root,
                                        f"sourcedata/{modality}/",
                                        sub,
                                        exp,
                                        n,
                                    )
                                )
                                trigger_ch = [
                                    (channels.data, channels.samples_per_second)
                                    for channels in bio_data.channels
                                    if channels.name == workflow["trigger"]["Channel"]
                                ]

                                vol_in_biopac = duration_counter(
                                    trigger_ch[0][0],
                                    trigger_ch[0][1],
                                    thr=workflow["trigger"]["trigThresh"],
                                )
                                vol_in_biopac = {exp: vol_in_biopac}
                                # Get the channel names and channel numbers
                                ch_names, chsel = order_channels(
                                    [channel.name for channel in bio_data.channels],
                                    workflow,
                                )
                                LGR.info(
                                    "finished counting duration of segments in physio "
                                    f"file for: {exp}"
                                )
                            else:
                                vol_in_biopac, ch_names, chsel = volume_counter(
                                    os.path.join(layout.root, f"sourcedata/{modality}/"),
                                    sub,
                                    workflow,
                                    ses=exp,
                                    tr=tr,
                                    trigger_ch=workflow["trigger"]["Channel"],
                                )
                                LGR.info(
                                    f"finished counting volumes in physio file for: {exp}"
                                )
                            try:
                                run_dict.update(
                                    {f"run-{idx_n+1:02d}": vol_in_biopac[exp]}
                                )
                            except TypeError:
                                # there were no triggers so stocking a place holder
                                run_dict = vol_in_biopac

                            ls_run_dict.append(run_dict)
                            ls_ch_names.append(list(ch_names))
                            ls_chsel.append(list(chsel))

                        # skip the session if we did not find the file
                        except KeyError:
                            continue
                nb_expected_runs[exp]["recorded_triggers"] = ls_run_dict
                nb_expected_runs[exp]["ch_names"] = ls_ch_names
                nb_expected_runs[exp]["chsel"] = ls_chsel
            except KeyError:
                nb_expected_runs[exp]["recorded_triggers"] = "No triggers found"
                LGR.info(
                    "Directory is empty or file is clobbered/No triggers:\n"
                    f"{os.path.join(path_to_source, sub, exp)}",
                )

                LGR.info(f"skipping :{exp} for task {name}")
        print("~" * 80)

    pprintpp.pprint(nb_expected_runs)

    if os.path.exists(os.path.join(root, 'code', sub)) is False:
        os.mkdir(os.path.join(root, 'code', sub))
    filename = f"{sub}_sessions.json"
    with open(os.path.join(root, 'code', sub, filename), "w") as f:
        json.dump(nb_expected_runs, f, indent=4)

    return nb_expected_runs


if __name__ == "__main__":
    get_info()
