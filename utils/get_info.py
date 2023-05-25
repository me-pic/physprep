# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning

"""another util for neuromod phys data conversion."""

import os
import glob
import json
import math
import click
import logging
import pprintpp
import fnmatch
import pandas as pd
from list_sub import list_sub
from neurokit2 import read_acqknowledge
LGR = logging.getLogger(__name__)


def volume_counter(root, sub, ses=None, tr=1.49, trigger_ch="TTL"):
    """
    Volume counting for each run in a session.

    Parameters
    ----------
    root : str
        Directory containing the biopac data. Example: "/home/user/dataset/sourcedata/physio".
    subject : str
        Name of path for a specific subject. Example: "sub-01".
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
        Each key lists the number of volumes/triggers in each run, including invalid volumes.
    """
    LGR = logging.getLogger(__name__)
    # Check directory
    if os.path.exists(root) is False:
        raise ValueError("Couldn't find the following directory: ", root)

    # List the files that have to be counted
    dirs = list_sub(root, sub, ses)
    ses_runs = {}
    # loop iterating through files in each dict key representing session returned by list_sub
    # for this loop, exp refers to session's name, avoiding confusion with ses argument
    for exp in dirs:
        LGR.info(f"counting volumes in physio file for: {exp}")
        for file in sorted(dirs[exp]):
            # reading acq
            bio_df, fs = read_acqknowledge(os.path.join(root, sub, exp, file))
            # find the correct index of Trigger channel
            if trigger_ch in bio_df.columns:
                trigger_index = list(bio_df.columns).index(trigger_ch)

            # initialize a df with TTL values over 4 (switch either ~0 or ~5)
            query_df = bio_df[bio_df[bio_df.columns[trigger_index]] > 4]

            # Define session length - this list will be less
            # memory expensive to play with than dataframe
            session = list(query_df.index)

            # maximal TR - the time distance between two adjacent TTL, now given by the ceiling value of the tr (but might be tweaked if needed)
            tr_period = fs * math.ceil(tr)

            # Define session length and adjust with padding
            start = int(session[0])
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
    LGR.info(f"Volumes for session :\n{ses_runs}")
    return ses_runs, bio_df.columns


@click.command()
@click.argument("root", type=str)
@click.argument("sub", type=str)
@click.option("--ses", type=str, default=None, required=False)
@click.option("--count_vol", type=bool, default=False, required=False)
@click.option("--show", type=bool, default=True, required=False)
@click.option("--save", type=str, default=None, required=False)
@click.option("--tr", type=float, default=None, required=False)
@click.option("--tr_channel", type=str, default=None, required=False)
@click.option("--scanning_sheet", type=str)
def call_get_info(
    root, sub, ses=None, count_vol=False, show=True, save=None, tr=None, tr_channel=None, scanning_sheet=None,
):
    """
    Call `get_info` function only if `get_info.py` is called as CLI

    For parameters description, please refer to the documentation of the `get_info` function
    """
    LGR = logging.getLogger(__name__)
    get_info(root, sub, ses, count_vol, show, save, tr, tr_channel, scanning_sheet)


def get_info(
    root, sub, ses=None, count_vol=False, show=True, save=None, tr=None, tr_channel=None, scanning_sheet=None,
):
    """
    Get all volumes taken for a sub.
    `get_info` pushes the info necessary to execute the phys2bids multi-run
    workflow to a dictionary. It can save it to `_volumes_all-ses-runs.json`
    in a specified path, or be printed in your terminal.
    The examples given in the Arguments section assume that the data followed this structure :
    home/
    └── users/
        └── dataset/
            └── sourcedata/
                └── physio/
                    ├── sub-01/
                    |   ├── ses-001/
                    |   |   └── file.acq
                    |   ├── ses-002/
                    |   |   └── file.acq
                    |   └── ses-0XX/
                    |       └── file.acq
                    └── sub-XX/
                        ├── ses-001/
                        |   └── file.acq
                        └── ses-0XX/
                            └── file.acq
    Arguments
    ---------
    root : str
        Root directory of dataset containing the data. Example: "/home/user/dataset/".
    sub : str
        Name of path for a specific subject. Example: "sub-01".
    ses : str
        Name of path for a specific session. Example: "ses-001".
    count_vol : bool
        Specify if you want to count triggers in physio file.
        Default to False.
    show : bool
        Specify if you want to print the dictionary.
        Default to True.
    save : str
        Specify where you want to save the dictionary in json format.
        If not specified, the output will be saved where you run the script.
        Default to None.
    tr : float
        Value of the TR used in the MRI sequence.
        Default to None.
    trigger_ch : str
        Name of the trigger channel used on Acknowledge.
        Defaults to None.

    Returns
    -------
    ses_runs_vols : dict
        Number of processed runs, number of expected runs, number of triggers/volumes per run,
        sourcedata file location.

    Examples
    --------
    In script
    >>> ses_runs_vols = get_info(root="/home/user/dataset/", sub="sub-01", ses="ses-001", count_vol=True, save="/home/user/dataset/info/", tr=2.0)
    In terminal
    >>> python get_info.py /home/user/dataset/ sub-01 --ses ses-001 --count_vol True --save /home/user/dataset/info/ --tr 2.0 --tr_channel 'Custom, HLT100C - A 5'
    """
    LGR = logging.getLogger(__name__)
    # list matches for a whole subject's dir
    ses_runs_matches = list_sub(
        os.path.join(root, "sourcedata/physio/"), sub, ses=ses, ext=".tsv", show=show
    )

    # go to fmri matches and get entries for each run of a session
    nb_expected_runs = {}

    # If there is a tsv file matching the acq file and the nii.gz files in root
    ses_info = list_sub(os.path.join(root, "sourcedata/physio/"), sub, ses, ext=".acq")
    # iterate through sessions and get _matches.tsv with list_sub dict
    for exp in sorted(ses_runs_matches):
        LGR.info(exp)
        if ses_info[exp] == []:
            LGR.info("No acq file found for this session")
            continue

        # initialize a counter and a dictionary
        nb_expected_volumes_run = {}
        tasks = []
        matches = glob.glob(os.path.join(root, sub, exp, "func", "*bold.json"))
        # we want only want to keep 1 of the two files per run
        if "mario" in matches[0]:
            matches = fnmatch.filter(matches, "*-mag*")
        matches.sort()
        # iterate through _bold.json
        for idx, filename in enumerate(matches):                
            task = filename.rfind(f"{exp}_") + 8
            task_end = filename.rfind("_")
            tasks += [filename[task:task_end]]

            # read metadata
            with open(filename) as f:
                bold = json.load(f)
            # we want to have the TR in a _bold.json to later use it in the volume_counter function
            tr = bold["RepetitionTime"]
            # we want to GET THE NB OF VOLUMES in the _bold.json of a given run
            try:
                nb_expected_volumes_run[f"{idx+1:02d}"] = bold["dcmmeta_shape"][-1]
            except KeyError:
                pprintpp.pprint(f"skipping {exp} because .json info non-existant")
                if scanning_sheet is not None:
                    pprintpp.pprint(f"checking scanning sheet for {sub}/{exp}")
                    df_sheet=pd.read_csv(scanning_sheet)
                    vol_idx=df_sheet[df_sheet[sub]==f"p{sub[4:]}_friends{exp[4:]}"].index+idx
                    vols=int(df_sheet["#volumes"].iloc[vol_idx])
                    nb_expected_volumes_run[f'{idx+1:02d}']=vols
                # log that we are unable to run the thing
                else:
                    LGR.info(f"Cannot access Nifti BIDS metadata")
                    continue


        # print the thing to show progress
        LGR.info(f"Nifti BIDS metadata; number of volumes per run:\n{nb_expected_volumes_run}")
        # push all info in run in dict
        nb_expected_runs[exp] = {}
        # the nb of expected volumes in each run of the session (embedded dict)
        nb_expected_runs[exp] = nb_expected_volumes_run
        nb_expected_runs[exp]["expected_runs"] = len(matches)
        # nb_expected_runs[exp]['processed_runs'] = idx  # counter is used here
        nb_expected_runs[exp]["task"] = tasks
        nb_expected_runs[exp]["tr"] = tr

        # save the name
        name = ses_info[exp]
        if name:
            name.reverse()
            nb_expected_runs[exp]["in_file"] = name

        if count_vol:
            run_dict = {}
            # check if biopac file exist, notify the user that we won't
            # count volumes
            try:
                # do not count the triggers in phys file if no physfile
                if (
                    os.path.isfile(
                        os.path.join(root, "sourcedata/physio", sub, exp, name[0])
                    )
                    is False
                ):
                    LGR.info(f"cannot find session directory for sourcedata :\n"
                        f"{os.path.join(root, 'sourcedata/physio', sub, exp, name[0])}")
                else:
                    # count the triggers in physfile otherwise
                    try:
                        vol_in_biopac, ch_names = volume_counter(
                            os.path.join(root, "sourcedata/physio/"),
                            sub,
                            ses=exp,
                            tr=tr,
                            trigger_ch=tr_channel,
                        )
                        LGR.info(f"finished counting volumes in physio file for: {exp}")

                        for i, run in enumerate(vol_in_biopac[exp]):
                            run_dict.update({f"run-{i+1:02d}": run})

                        nb_expected_runs[exp]["recorded_triggers"] = run_dict
                        nb_expected_runs[ses]["ch_names"] = ch_names

                    # skip the session if we did not find the _bold.json
                    except KeyError:
                        continue
            except KeyError:
                nb_expected_runs[exp]["recorded_triggers"] = "No triggers found"
                LGR.info(
                    "Directory is empty or file is clobbered/No triggers:\n"
                    f"{os.path.join(root, 'sourcedata/physio', sub, exp)}",
                )

                LGR.info(f"skipping :{exp} for task {filename}")
        print("~" * 80)

    if show:
        pprintpp.pprint(nb_expected_runs)

    if save is not None:
        nb_expected_runs = pprintpp.pformat(nb_expected_runs)
        if os.path.exists(os.path.join(save, sub)) is False:
            os.mkdir(os.path.join(save, sub))
        filename = f"{sub}_volumes_all-ses-runs.json"
        with open(os.path.join(save, sub, filename), "w") as fp:
            fp.write(nb_expected_runs)
    return nb_expected_runs


if __name__ == "__main__":
    log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    call_get_info()
