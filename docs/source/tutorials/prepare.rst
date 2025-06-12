Put those data in BIDS
----------------------

`Physprep` does not include functions to automatically BIDSified raw physiological 
timeseries in its workflow. The reasoning behind this choice is that the BIDSification of 
raw (source) physiological data can sometimes fail, and it is the user responsability to 
ensure that this step is done properly before running the `Physprep` workflow. However, 
`Physprep` includes some helper functions that could help you BIDSifying your data ! This 
tutorial shows how to use those functions.

Preparing your files
^^^^^^^^^^^^^^^^^^^^

You will first need to general a file matching your source physiological data (e.g. your 
`.acq` files to your `bold` or `eeg` data). For that, we will use the `match_acq_bids` 
function, which take into account that there is a sub-<sub>[_ses-<ses>]_scans.tsv file in 
your sub-<sub>/[ses-<ses>/] subdirectories (see `BIDS documentation <https://bids-specification.readthedocs.io/en/stable/modality-agnostic-files.html#scans-file>`_ 
for more details on the scans files).

.. code-block:: python

    from physprep.prepare.match_acq_bids import match_acq_bids

    # Change the paths to reflect your own
    bids_path = '/path/to/your/bids/dataset'
    source_path = '/path/where/to/find/your/acq/files'

    match_acq_bids(bids_path, biopac_path)

Let's see what you have created...

.. code-block:: python

    import pandas as pd

    # Change the path to reflect your own
    matches = pd.read_csv('/path/sourcedata/eeg_or_func/sub-*/ses-*/*matches.tsv', sep='\t')
    print(matches)


.. note::
    `match_acq_bids`` uses the timestamp in the `.acq` and `.edf` (if you have a EEG dataset) or `bold` (if you have a fMRI dataset) files to match the `.acq` file with the concurrent `*eeg.edf` or `*bold.nii.gz` file. However, it is possible that this function fails if, for example, those files were recorded on computers with unsynchronized clock. Because of that, users NEED to validate the `.acq` files in `sourcedata/eeg/` (for EEG dataset) or `sourcedata/func` (for fMRI dataset), and make sure that the content of the `*matches.tsv`is adequate.

Retrieving info for BIDSification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The bidsification for the `.acq` files relies on the python library `phys2bids <https://physiopy.github.io/libraries/phys2bids/>`_. In order to properly segment and bidsify our `.acq` recordings with `physprep`, we need to retrieve some info from our data using the `get_info` function.

.. code-block:: python

    import os
    from physprep.prepare.get_info import get_info

    # Change the path to reflect your own
    workflow_strategy = '/path/to/your/workflow_strategy.json'
    help(get_info)
    get_info(path, sub, workflow_strategy, ses=ses, count_vol=True)

Let's see what we have just created...

.. code-block:: python

    # Change the path to reflect your own
    info_sub = load_json('/path/code/sub-*/sub-*_sessions.json')
    print(info_sub)


BIDSifying the physio data
^^^^^^^^^^^^^^^^^^^^^^^^^^

We are now ready to BIDSify the data ! For that, we'll use the `convert` function.

.. code-block:: python

    from physprep.prepare.convert import convert

    convert(path, sub, ses=ses, pad=0, overwrite=False, validate=True)

If the output of `convert` is correct, we can know use `Physprep CLI <./cli.html>`_ !