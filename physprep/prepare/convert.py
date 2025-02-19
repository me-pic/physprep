# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""
Physiological data conversion to BIDS.
"""

import gc
import logging
import os
import shutil

from bids import BIDSLayout
from phys2bids.phys2bids import phys2bids

from physprep import utils


def convert(root, sub, ses=None, pad=0, overwrite=False):
    """
    Phys2Bids conversion for one subject data.

    Parameters
    ----------
    root : str
        Root of the BIDS directory.
    subject : str
        Name of path for a specific subject. Example: "sub-01".
    ses : list
        Name of path for specific session(s). Example: "ses-01" or
        ["ses-001", "ses-002"].
        Default to None.

    See also
    --------
    https://phys2bids.readthedocs.io/en/latest/index.html
    """
    logger = logging.getLogger(__name__)
    # fetch info
    logger.info(f"Reading fetcher in:\n{os.path.join(root, 'code', sub)}")
    fetcher = f"{sub}_sessions.json"
    layout = BIDSLayout(root)
    info = utils.load_json(os.path.join(root, "code", sub, fetcher))
    # define sessions
    if ses is None:
        ses = sorted(list(info.keys()))
        # Remove the sessions that are already processed
        if overwrite is False:
            existing = [f"ses-{s}" for s in layout.get_sessions()]
            setA = set(ses)
            # Get new set with elements that are only in sessions but not in
            # existing
            ses = sorted(list(setA.difference(existing)))

    elif isinstance(ses, list) is False:
        ses = [ses]

    # iterate through info
    for col in sorted(ses):
        # If the directory does not contain ses-* sub-folders
        if not layout.get_sessions():
            ses_id = "001"
            indir = os.path.join(root, "sourcedata", info[col]["concurrent_with"], sub)
        else:
            ses_id = col[-3:]
            indir = os.path.join(
                root, "sourcedata", info[col]["concurrent_with"], sub, col
            )

        if info[col] is None:
            logger.info(f"Empty session : {col}")
            continue

        for idx, file in enumerate(info[col]["in_file"]):
            # Define ch_name and trigger idx
            # skip empty sessions
            if "ch_names" not in info[col]:
                logger.info(
                    "No channel names provided nor found; "
                    "can't find the trigger to segment"
                )
                continue
            else:
                chtrig = info[col]["ch_names"][idx].index("trigger") + 1

            # Iterate through files in each session and run phys2bids
            logger.info(f"Converting : {col}")
            if info[col]["concurrent_with"] == "eeg":
                tr = []
                for tr_idx in info[col]["recorded_triggers"][idx][f"run-0{idx+1}"]:
                    for t in info[col]["tr"][idx]:
                        if t - tr_idx < 0.02:
                            tr.append(t)
                nte = [1] * len(tr)
            else:
                tr = info[col]["tr"][idx]
                nte = info[col]["recorded_triggers"][idx][f"run-0{idx+1}"]

            phys2bids(
                file,
                info=False,
                indir=indir,
                outdir=os.path.join(root, sub, col, info[col]["concurrent_with"]),
                heur_file=None,
                sub=sub[-2:],
                ses=ses_id,
                chtrig=chtrig,
                chsel=info[col]["chsel"][idx],
                num_timepoints_expected=nte,
                tr=tr,
                thr=4,
                pad=pad,
                ch_name=info[col]["ch_names"][idx],
                yml="",
                debug=True,
                quiet=False,
            )
            # Rename files
            files_concurrent = [
                file
                for file in layout.get(
                    session=col[-3:], suffix=info[col]["concurrent_with"]
                )
                if "json" not in file.filename
            ]
            file_concurrent = [
                file
                for file in files_concurrent
                if info[col]["concurrent_file"][idx] in file.filename
            ][0]
            bids_entities = {
                "subject": file_concurrent.get_entities()["subject"],
                "session": file_concurrent.get_entities()["session"],
                "datatype": file_concurrent.get_entities()["datatype"],
                "task": file_concurrent.get_entities()["task"],
                "run": file_concurrent.get_entities()["run"],
                "suffix": "physio",
                "extension": "json",
            }
            pattern = "sub-{subject}[/ses-{session}]/{datatype}/sub-{subject}[_ses-{session}][_task-{task}][_run-{run}][_{suffix}].{extension}"
            # Renaming the *_physio.json file
            renamed = layout.build_path(bids_entities, pattern, validate=False)
            os.rename(
                os.path.join(
                    root,
                    sub,
                    col,
                    info[col]["concurrent_with"],
                    file.split(".")[0] + ".json",
                ),
                renamed,
            )
            # Renaming the *_physio.tsv.gz file
            bids_entities["extension"] = "tsv.gz"
            renamed = layout.build_path(bids_entities, pattern, validate=False)
            os.rename(
                os.path.join(
                    root,
                    sub,
                    col,
                    info[col]["concurrent_with"],
                    file.split(".")[0] + ".tsv.gz",
                ),
                renamed,
            )
        gc.collect()

        # Make sure folder exists
        outdir_phys2bids = os.path.join(root, "code", sub, col, "phys2bids")
        os.makedirs(outdir_phys2bids, exist_ok=True)

        # Move code/ folder to appropriate directory
        shutil.move(
            os.path.join(root, sub, col, info[col]["concurrent_with"], "code"),
            outdir_phys2bids,
        )

        print("~" * 30)
