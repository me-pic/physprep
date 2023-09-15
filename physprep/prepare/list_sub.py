# -*- coding: utf-8 -*-
# !/usr/bin/env python
"""CLI for physio utils."""

import json
import logging
import os

import pprintpp

LGR = logging.getLogger(__name__)


def list_sub(root, sub, ses=None, ext=".acq", save=None, show=False):
    """
    List a subject's files.

    Returns a dictionary entry for each session in a subject's directory.
    Each entry is a list of files for a given subject/ses directory.
    If ses is given, only one dictionary entry is returned.

    Arguments
    ---------
    root : str
        Root directory of dataset containing the data. Example: "home/user/dataset/".
    sub : str
        Name of path for a specific subject. Example: "sub-01".
    ses : str
        Name of path for a specific session. Example: "ses-001".
        Default to None.
    ext : str
        Specify the extension of the files are we looking for.
        Default to '.acq'; biosignals from biopac.
    save : str
        Specify where you want to save the lists.
        If not specified, the output will not be saved.
        Default to None.
    show : bool
        If True, prints the output dict.
        Default to False.

    Returns
    -------
    ses_runs : dict
        Dictionary containing the sessions id in the subject's folder, and the name
        of the acqknowledge file.
        Returned if `ses` is not specified.
    files : dict
        Dictionary of acqknowledge filenames for each session.
        Returned if `ses` is specified.

    Examples
    --------
    In script
    >>> ses_runs = list_sub(root="/home/user/dataset/", sub="sub-01")
    In terminal
    >>> python list_sub.py /home/user/dataset/ sub-01 --ses ses-001 --ext acq --show True
    """
    # Check the subject's
    path_sub = os.path.join(root, sub)
    if os.path.isdir(path_sub) is False:
        raise ValueError("Couldn't find the subject's path \n", os.path.join(root, sub))
    file_list = []
    ses_runs = {}
    ses_list = os.listdir(path_sub)
    # list files in only one session
    if ses is not None:
        dir = os.path.join(path_sub, ses)
        # if the path exists, list .acq files
        if os.path.exists(dir):
            for filename in os.listdir(dir):
                if filename.endswith(ext):
                    file_list += [filename]
            if show:
                print("list of sessions in subjet's directory: ", ses_list)
                print("list of files in the session:", file_list)

            # return a dictionary entry for the specified session
            files = {str(ses): file_list}

            if save is not None:
                with open(
                    os.path.join(save, sub, ses, "list_sub_acq_files.json"),
                    "w",
                ) as fp:
                    json.dump(files, fp)

            return files
        else:
            print("list of sessions in subjet's directory: ", ses_list)
            raise Exception("Session path you gave does not exist")

    # list files in all sessions (or here, exp for experiments)
    elif os.path.isdir(os.path.join(path_sub, ses_list[0])) is True:
        for exp in ses_list:
            if exp.endswith(".json"):
                continue
            # re-initialize the list
            file_list = []
            # iterate through directory's content
            for filename in os.listdir(os.path.join(path_sub, exp)):
                if filename.endswith(ext):
                    file_list += [filename]

            # save the file_list as dict item
            ses_runs[exp] = file_list

        # display the lists (optional)
        if show:
            pprintpp.pprint(ses_runs)

        if save is not None:
            with open(os.path.join(save, sub, "list_sub_acq_files.json"), "w") as fp:
                json.dump(ses_runs, fp)

        return ses_runs
    # list files in a sub directory without sessions
    else:
        # push filenames in a list
        for filename in os.listdir(path_sub):
            if filename.endswith(ext):
                file_list += [filename]
        # store list
        ses_runs["files"] = file_list

        if save is not None:
            with open(os.path.join(save, sub, "list_sub_acq_files.json"), "w") as fp:
                json.dump(ses_runs, fp)

        # return a dictionary of sessions each containing a list of files
        return ses_runs
