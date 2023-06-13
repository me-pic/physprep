# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""Neuromod phys data rename converted files."""

import glob
import pandas as pd
import numpy as np
import click
import os
import logging


@click.command()
@click.argument("indir", type=click.Path(exists=True), required=True)
@click.argument("sub", type=str, required=True)
@click.option("sessions", type=str, nargs="*")
def co_register_physio(indir, sub, sessions=None):
    """
    Comply to BIDS and co-register functional acquisitions.

    Rename valid files and remove invalid files for 1 subject's directory

    Parameters:
    ------------
    indir : path
        directory to save data and to retrieve acquisition info (`.json file`)
    subject : string
        name of path for a specific subject (e.g.'sub-03')
    sessions : list
        specific session numbers can be listed (e.g. ['ses-001', 'ses-002']
    Returns:
    --------
    BIDS-compliant /func directory for physio files
    """
    logger = logging.getLogger(__name__)
    # fetch info
    info = pd.read_json(f"{indir}{sub}/{sub}_volumes_all-ses-runs.json")
    # define sessions
    if sessions is None:
        sessions = info.columns
    elif isinstance(sessions, list) is False:
        sessions = [sessions]

    # iterate through sesssions
    for ses in sessions:
        logger.info(f"renaming files in session : {ses}")

        # list files in the session
        tsv = glob.glob(f"{indir}{sub}/{ses}/*.tsv.gz")
        tsv.sort()

        if tsv is None or len(tsv) == 0:
            print(f"no physio file for {ses}")
            continue

        json = glob.glob(f"{indir}{sub}/{ses}/*.json")
        json.sort()

        log = glob.glob(f"{indir}{sub}/{ses}/code/conversion/*.log")
        log.sort()

        png = glob.glob(f"{indir}{sub}/{ses}/code/conversion/*.png")
        png.sort()

        # sanitize list of triggers
        triggers = list(info[ses]["recorded_triggers"].values())
        triggers = list(np.concatenate(triggers).flat)

        # check sanity of info - expected runs is number of runs in BOLD sidecar
        # if info[ses]['expected_runs'] is not info[ses]['processed_runs']:
        #    print(f"Expected number of runs {info[ses]['expected_runs']} "
        #          "does not match info from neuroimaging metadata")

        if len(info[ses]["task"]) is not info[ses]["expected_runs"]:
            logger.info("Number of tasks does not match expected number of runs")
            continue

        if info[ses]["recorded_triggers"].values is None:
            logger.info(
                f"No recorded triggers information - check physio files for {ses}"
            )
            continue

        if len(info[ses]["task"]) == 0:
            logger.info(f"No task name listed ; skipping {ses}")
            continue

        if len(info[ses]["task"]) == 1:
            this_one = triggers.index(info[ses]["01"])
            to_be_del = list(range(0, len(triggers)))
            to_be_del.remove(this_one)

        # if input is normal, then check co-registration
        else:
            to_be_del = []

            # remove files that don't contain enough volumes
            for idx, volumes in enumerate(triggers):
                # NOTE: this should not be hardcoded
                if volumes < 400:
                    to_be_del.append(idx)

        # these can be safely removed
        for idx in to_be_del:
            os.remove(tsv[idx])
            os.remove(json[idx])
            os.remove(log[idx])
            os.remove(png[idx])

        # theses are to be kept
        triggers = np.delete(triggers, to_be_del)
        tsv = np.delete(tsv, to_be_del)
        json = np.delete(json, to_be_del)
        log = np.delete(log, to_be_del)
        png = np.delete(png, to_be_del)

        # check if number of volumes matches neuroimaging JSON sidecar
        for idx, volumes in enumerate(triggers):
            i = f"{idx+1:02d}"
            logger.info(info[ses][i])
            if volumes != info[ses][i]:
                logger.info(
                    f"Recorded triggers info for {ses} does not match with "
                    f"BOLD sidecar ({volumes} != {info[ses][i]})\n"
                    f"Skipping {ses}"
                )
                break

            else:
                os.rename(
                    tsv[idx],
                    f"{indir}{sub}/{ses}/{sub}_{ses}_{info[ses]['task'][idx]}_physio.tsv.gz",
                )
                os.rename(
                    json[idx],
                    f"{indir}{sub}/{ses}/{sub}_{ses}_{info[ses]['task'][idx]}_physio.json",
                )


if __name__ == "__main__":
    log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    co_register_physio()
