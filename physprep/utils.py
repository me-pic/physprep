"""Utilities for Physprep"""

import json
import pickle


def load_json(filename):
    """
    Parameters
    ----------
    filename : str
        File path of the .json to load.

    Returns
    -------
    data : dict
        Dictionary with the content of the .json passed in argument.
    """
    try:
        with open(filename, "r") as tmp:
            data = json.loads(tmp.read())
        tmp.close()
    except Exception:
        try:
            with open(filename, "rb") as tmp:
                data = pickle.load(tmp)
            tmp.close()
        except Exception:
            pass

    return data


def create_config():
    filename = input("Enter the name of your config file \n")
    tmp = None
    order = "0"
    valid_filter_type = ["butterworth", "fir", "bessel", "savgol", "", " "]
    while tmp not in ["", " "]:
        while True:
            tmp = input(
                "Enter the filter type among the following: butterworth, fir, bessel, "
                "savgol"
            )
            if tmp not in valid_filter_type:
                print(
                    "Please enter a valid filter type (i.e., butterworth, fir, bessel, "
                    "savgol)"
                )
                continue
            else:
                method = tmp
                break

        if tmp in ["", " "]:
            break
        else:
            lowcut = input(
                "Enter the lower cutoff frequency (Hz). If you do not want to apply "
                "a high pass or band pass filter, just press enter"
            )
            highcut = input(
                "Enter the higher cutoff frequency (Hz). If you do not want to apply "
                "a low pass filter or band pass filter, just press enter"
            )
            while True:
                order = input("Enter the filter order. Must be a positive integer.")
                if order.isdigit() is False:
                    print("Please enter a positive integer")
                    continue
                else:
                    break

    return filename, method, lowcut, highcut, order
