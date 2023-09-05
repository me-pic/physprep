# Usage Notes

```bash
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

```bash
from utils import create_config_workflow

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

```bash
from utils import create_config_preprocessing

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