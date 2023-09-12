# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""Neuromod phys data rename converted files."""

import glob
import logging
import os

import click
import numpy as np
import pandas as pd


@click.command()
@click.argument("indir", type=click.Path(exists=True), required=True)
@click.argument("sub", type=str, required=True)
@click.option("--ses", type=str, required=False, default=None)
@click.option("--neuromod", type=bool, required=False, default=False)
@click.option("--min_volumes", type=int, required=False, default=350)
def co_register_physio(indir, sub, ses=None, neuromod=False, min_volumes=350):
    """
    Comply to BIDS and co-register functional acquisitions.

    Rename valid files and remove invalid files for 1 subject's directory

    Parameters:
    ------------
    indir : path
        directory to save data and to retrieve acquisition info (`.json file`)
    subject : string
        name of path for a specific subject (e.g.'sub-03')
    ses : list
        specific session numbers can be listed (e.g. ['ses-001', 'ses-002']
    neuromod : bool
        Set to True if working with Neuromod data. Default is False.
    min_volumes : int
        Minimum number of  volumes expected for a run. If a run has less volumes than this
        number, it will be removed. Default is 350.
    Returns:
    --------
    BIDS-compliant /func directory for physio files
    """
    logger = logging.getLogger(__name__)
    # fetch info
    info = pd.read_json(os.path.join(indir, sub, f"{sub}_volumes_all-ses-runs.json"))
    # define sessions
    if ses is None:
        ses = info.columns
    elif isinstance(ses, list) is False:
        ses = [ses]

    # iterate through sessions
    for s in ses:
        logger.info(f"renaming files in session : {s}")

        # list files in the session
        tsv = glob.glob(os.path.join(indir, sub, s, "*.tsv.gz"))
        tsv.sort()

        if tsv is None or len(tsv) == 0:
            print(f"no physio file for {s}")
            continue

        json = glob.glob(os.path.join(indir, sub, s, "*.json"))
        json.sort()

        log = glob.glob(os.path.join(indir, sub, s, "code", "conversion", "*.log"))
        log.sort()

        png = glob.glob(os.path.join(indir, sub, s, "code", "conversion", "*.png"))
        png.sort()

        # sanitize list of triggers
        triggers = list(info[s]["recorded_triggers"].values())
        triggers = list(np.concatenate(triggers).flat)

        if neuromod:
            # remove padding trigger for videogames tasks
            if any(game in indir for game in ["mario", "shinobi"]):
                triggers = [trigger - 1 if trigger > 200 else trigger for trigger in triggers]

        if len(info[s]["task"]) is not info[s]["expected_runs"]:
            logger.info("Number of tasks does not match expected number of runs")
            continue

        if info[s]["recorded_triggers"].values is None:
            logger.info(
                "No recorded triggers information - check physio files " f"for {s}"
            )
            continue

        if len(info[s]["task"]) == 0:
            logger.info(f"No task name listed ; skipping {s}")
            continue

        if len(info[s]["task"]) == 1:
            this_one = triggers.index(info[s]["01"])
            to_be_del = list(range(0, len(triggers)))
            to_be_del.remove(this_one)

        # if input is normal, then check co-registration
        else:
            to_be_del = []
            # remove files that don't contain enough volumes
            for idx, volumes in enumerate(triggers):
                if volumes < min_volumes:
                    to_be_del.append(idx)

        # these can be safely removed
        for idx in to_be_del:
            os.remove(tsv[idx])
            os.remove(json[idx])
            os.remove(log[idx])

        # these are to be kept
        triggers = np.delete(triggers, to_be_del)
        tsv = np.delete(tsv, to_be_del)
        json = np.delete(json, to_be_del)

        # check if number of volumes matches neuroimaging JSON sidecar
        for idx, volumes in enumerate(triggers):
            i = f"{idx+1:02d}"
            logger.info(info[s][i])
            filename_tsv = f"{sub}_{s}_{info[s]['task'][idx]}_physio.tsv.gz"
            filename_json = f"{sub}_{s}_{info[s]['task'][idx]}_physio.json"
            if volumes != info[s][i]:
                logger.info(
                    f"Recorded triggers info for {s} does not match with "
                    f"BOLD sidecar ({volumes} != {info[s][i]})\n"
                    f"Skipping {s}"
                )
                break

            else:
                os.rename(
                    tsv[idx],
                    os.path.join(indir, sub, s, filename_tsv),
                )
                os.rename(
                    json[idx],
                    os.path.join(indir, sub, s, filename_json),
                )


if __name__ == "__main__":
    log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    co_register_physio()
