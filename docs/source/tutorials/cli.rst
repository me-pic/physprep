One line to (almost) run them all
---------------------------------

`Physprep` is a Command Line Interface (CLI) that automatically (pre)processed and assessed
the quality of peripheral physiological data. It takes in input a configuration file and 
a BIDS directory containing raw physiological data, and return cleaned timeseries,
extracted features, quality assessment and summary quality report. 

For this example, we will create a BIDS dataset with simulated data using a function
in the physprep `data` submodule.

.. code-block:: python

    from physprep.data import data

    # Let's see what parameters this function can take
    help(data.create_bids_dataset)
    # Let's create the dataset
    data.create_bids_dataset('/path/to/save/those/data')


Once the dataset has been created, we can run the command in our terminal ! 

.. warning::
    Make sure you have `installed Physprep <./../installation.md>`_ in your virtual environment before.

For this example, we will use the default configuration file provided in `Physprep`, i.e.
the `neuromod` json file. You can also replace that first argument by your configuration
file that you have created previously (see `using configuration files <./../usage.html>`_). We 
will also specified the path where our raw physio data are, the subject id, the 
session id, the path where we want to save the derivatives, and the `--save_report` flag 
to save a summary descriptive report of our preprocessed signals.

.. note::
    To run the pipeline on all subjects and sessions in your raw physio directory, do not 
    use the `--sub`` and `--ses arguments`

.. code-block:: bash

    physprep 'neuromod' '/path/where/the/simulated/data/have/been/saved' --sub '01' --ses '001' --derivatives_dir '/path/to/save/derivatives/' --save_report
