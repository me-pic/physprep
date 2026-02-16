# Usage Notes

```bash
usage: physprep [workflow_strategy WORKFLOW_STRATEGY] [indir_bids INDIR_BIDS] [--sub SUB]
                [--ses SES] [--derivatives_dir DERIVATIVES_DIR] [--save_report]
                [--track_carbon] [--country_code] [--help]

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

Options:
  --sub TEXT               Subject label (e.g., `01`).
  --ses TEXT               Session label (e.g., `001`).
  --derivatives_dir PATH   Path to the derivatives directory. If not specified,
                           derivatives will be save in the `indir_bids` directory.
  --save_report            If specified, an quality report will be generated and
                           saved for each run.
  --track_carbon           If specified, carbon tracker will be used to track power
                           use of the pipeline.
  --country_code           Country ISO code used by carbon trackers. Will only be
                           used if --track_carbon flag is specified.
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
        "id": "<signal id>",
        "Description": "description of the signal",
        "Units": "recording units",
        "Channel": "<name_of_the_channel>",
        "preprocessing_strategy": "<name_of_the_preprocessing_strategy>",
        "qa_strategy": "<name_of_the_qa_strategy>",
    },
}
```

See example in `physprep/data/workflow_strategy`.

Every item within the file corresponds to a supported physiological signal, such as "ECG",
"PPG", "EDA", or "RSP". For each modality to (pre-)process, users have the possibility
to specify the channel name in the acquisition file, the filename/directory of
the preprocessing strategy, and the filename/directory of the qa strategy.

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

### QA strategy

```
[
    {
        "sliding": {
            "duration": <duration_in_seconds>,
            "step": <duration_in_seconds>
        }
    },
    {
        "metric": "<metric>",
        "feature": "<feature_name>"
    },
    {
        "metric": "Quality",
        "feature": "<feature_name>",
        "threshold": {"feature": "criterion"}
    }
]
```

See presets in `physprep/data/qa_strategy`.

Users also have the possibility to create their own qa configuration file using
the `create_config_qa` function from the utils module:

```python
from physprep.utils import create_config_qa

# Check the function documentation
help(create_config_qa)

# Call the function specifying the parameters value. You will be asked to enter different
# parameters from a set of currently supported values.
create_config_qa(
    "/path/to/your/config/file",
    "your_amazing_qa.json",
    "ECG",
    overwrite=True
)
```

(carbon-tracker)=
## Carbon tracker

`physprep` workflow integrates [`codecarbon`](https://mlco2.github.io/codecarbon/#) to track the energy usage and carbon emissions related to the preprocessing of your physiological data. The carbon emissions are estimated based on the energy consumption and the carbon intensity of electricity for the country specified in the `--country_code` parameter (default to "CAN" - Canada). However, since the carbon intensity of the electricity is based on the national average, a more precise estimate could be calculated by taking the carbon intensity value of a specific province/state/region (see [electricitymaps](https://app.electricitymaps.com)). To track the energy consumption with physprep, specify the flag `--track_carbon` in the CLI. The `codecarbon` output will be saved under `code/emissions.csv`.
