# Workflow

0. Setup the configuration files: either use the presets (see
 `physprep/data/workflow_strategy`, `physprep/data/preprocessing_strategy`, 
 `physprep/data/qa_strategy`) or create
 your {doc}`own configuration files</usage/>`.
1. Prepare the data: segment continuous physiological recordings into BIDSified runs based
 on the MR metadata. This part is not automatically included in the workflow, but Physprep
 provides functions to BIDSified your dataset.
2. Preprocess the data: filter and resample the data based on the parameters specified in
 the `preprocessing_strategy` configuration file.
3. Extract the features: detect or compute different features for each type of supported
 physiological signals (e.g. peaks, troughs, amplitude).
4. Generate a report: provide descriptive metrics as specified in the `qa_strategy`
 configuration file, and interactive plots associated with the processed signals.
