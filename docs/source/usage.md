# Usage Notes

```bash
usage: physprep [--help] [--workflow_strategy WORKFLOW_STRATEGY] [--indir_bids INDIR_BIDS] [--sub SUB]
                [--ses SES] [--indir_raw_physio INDIR_RAW_PHYSIO] [--skip_match_acq_bids]
                [--skip_convert] [--padding PADDING]

Preprocess raw physiological data acquired in MRI, extract features, and
generate quality report.

Required:
  workflow strategy        Name of the workflow_strategy if using a preset.
                           It is also possible to use a custom file by
                           providing the path to a JSON file containing
                           workflow strategy. In that case, please check
                           Physprep documentation to make sure your file is
                           properly formatted.
  indir_bids               Path to the directory containing the BIDS-like
                           dataset.
  sub                      Subject label.

Options:
  --ses TEXT               Session label.
  --indir_raw_physio PATH  Path to the directory containing the raw
                           physiological data. Specify if raw
                           physiological data is not in the BIDS
                           directory. For more details, about the BIDS
                           data structure, please refer to the
                           documentation.
  --skip_match_acq_bids    If specified, the workflow will not match the
                           acq files with the bold files. If acq files are
                           already organized properly, this flag can be
                           specified. For more details, see the
                           documentation of the mathc_acq_bids.py script.
  --skip_convert           If specified, the workflow will not convert the
                           physiological data recordings in BIDS format.
                           This implies that the data is already segmented
                           in runs, organized in a BIDS-like structure
                           (i.e., one tsv.gz and one json file per run),
                           and named following the BIDS recommandations.
  --padding INTEGER        Time (in seconds) of padding to add at the
                           beginning and end of each run. Default to 9.
  --help                   Show this message and exit.
```

## Using configuration files

Physprep pipeline uses configuration files to define the preprocessing steps for each
modality (presets can be found in `physprep/data`). The configuration files are organized
in two levels:
- `physprep/data/workflow_strategy` defining the type of physiological signal to process
and the preprocessing strategy for each type of signal.
- `physprep/data/preprocessing_strategy` defining the preprocessing steps and parameters
related to them.

### Workflow strategy

The workflow configuration file is defined as the following:

```
{
    "<name_of_the_physiological_signal>": {
        "channel": "<name_of_the_channel>",
        "preprocessing_strategy": "<name_of_the_preprocessing_strategy>",
    },
}
```

See examples in `physprep/data/workflow_strategy`.

Every item within the file corresponds to a supported physiological signal, such as "ECG",
"PPG", "EDA", or "RSP". For each modality to (pre-)process, users have the possibility
to specify the channel name in the acquisition file, and the filename/directory of
the preprocessing strategy.

An helper function is provided in the `utils` module to create a customized workflow
configuration file:

```python
from physprep.utils import create_config_workflow

# Check the function documentation
help(create_config_workflow)

# Call the function specifying the parameters value. You will be asked to enter different
# parameters from a set of currently supported values.
create_config_workflow(
    "/path/to/your/config/file",
    "your_amazing_workflow_strategy.json",
    overwrite=True
)
```

### Preprocessing strategy

```
[
    {
        "step": "<name_of_the_step>",
        "parameters": {
            "<step_parameters>": "<value>",
            ...
        },
        "reference": {
            "<reference_parameters>": "value",
            ...
        }
    },
]
```

See presets in `physprep/data/preprocessing_strategy`.

Users also have the possibility to create their own preprocessing configuration file using
the `create_config_preprocessing` function from the utils module:

```python
from physprep.utils import create_config_preprocessing

# Check the function documentation
help(create_config_preprocessing)

# Call the function specifying the parameters value. You will be asked to enter different
# parameters from a set of currently supported values.
create_config_preprocessing(
    "/path/to/your/config/file",
    "your_amazing_preprocessing_strategy.json",
    overwrite=True
)
```
