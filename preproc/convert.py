# -*- coding: utf-8 -*-
# !/usr/bin/env python -W ignore::DeprecationWarning
"""Neuromod phys data conversion."""

from phys2bids.phys2bids import phys2bids
import argparse
import sys
import pandas as pd
import gc
import os

import sys
sys.path.append("../utils")
from CLI import _get_parser3

def neuromod_phys2bids(sourcedata, scratch, sub, ses=None, tr=1.49):
    """
    Phys2Bids conversion for one subject data

    Parameters:
    ------------
    sourcedata : path
        main directory containing the biopac data (e.g. /to/dataset/info)
    scratch : path
        directory to save data and to retrieve acquisition info (`.json file`)
    subject : string
        name of path for a specific subject (e.g.'sub-03')
    sessions : list
        specific session numbers can be listed (e.g. ['ses-001', 'ses-002']
    tr : float
        tr value used for mri acquisition
    ch_name : list
        specify the name of the channels in the acqknowledge file

    Returns:
    --------
    phys2bids output
    """
    # fetch info
    info = pd.read_json(os.path.join(scratch, sub, f"{sub}_volumes_{ses}-runs.json"))
    # define sessions
    if ses is None:
        ses = info.columns
    elif isinstance(ses, list) is False:
        ses = [ses]
    # Define ch_name
        if ch_name is None:
            print("Warning: you did not specify a value for ch_name, the values that will be use are the following: ")
            ch_name = ["EDA", "PPG", "ECG", "TTL", "RSP"]
            print(ch_name)
            print("Please make sure, those values are the right ones !")
            chtrig=4
        else:
            chtrig = ch_name.index('TTL')

    # iterate through info
    for col in ses:
        # skip empty sessions
        if info[col] is None:
            continue
        print(col)

        # Iterate through files in each session and run phys2bids
        filename = info[col]["in_file"]
        if filename is list:
            for i in range(len(filename) - 1):
                phys2bids(
                    filename[i],
                    info=False,
                    indir=os.path.join(sourcedata, sub, col),
                    outdir=os.path.join(scratch, sub, col),
                    heur_file=None,
                    sub=sub[-2:],
                    ses=col[-3:],
                    chtrig=chtrig,
                    chsel=None,
                    num_timepoints_expected=info[col]["recorded_triggers"][
                        f"run-0{i+1}"
                    ],
                    tr=info[col]["tr"],
                    thr=4,
                    pad=9,
                    ch_name=info[col]["ch_name"],
                    yml="",
                    debug=False,
                    quiet=False,
                )
        else:
            try:
                phys2bids(
                    filename,
                    info=False,
                    indir=os.path.join(sourcedata, "physio", sub, col),
                    outdir=os.path.join(scratch, sub, col),
                    heur_file=None,
                    sub=sub[-2:],
                    ses=col[-3:],
                    chtrig=chtrig,
                    chsel=None,
                    num_timepoints_expected=info[col]["recorded_triggers"]["run-01"],
                    tr=info[col]["tr"],
                    thr=4,
                    pad=9,
                    ch_name=info[col]["ch_name"],
                    yml="",
                    debug=False,
                    quiet=False,
                )
            except AttributeError:
                filename.sort()
                for i in range(len(filename)):
                    print(i)
                    phys2bids(
                        filename[i],
                        info=False,
                        indir=os.path.join(sourcedata, sub, col),
                        outdir=os.path.join(scratch, sub, col),
                        heur_file=None,
                        sub=sub[-2:],
                        ses=col[-3:],
                        chtrig=chtrig,
                        chsel=None,
                        num_timepoints_expected=info[col]["recorded_triggers"][
                            f"run-0{i+1}"
                        ],
                        tr=info[col]["tr"],
                        thr=4,
                        pad=9,
                        ch_name=info[col]["ch_name"],
                        yml="",
                        debug=False,
                        quiet=False,
                    )

            except TypeError:
                print(f"No input file for {col}")
                continue
        gc.collect()
        print("~" * 30)


def _main(argv=None):
    options = _get_parser2().parse_args(argv)
    neuromod_phys2bids(**vars(options))


if __name__ == "__main__":
    _main(sys.argv[1:])
