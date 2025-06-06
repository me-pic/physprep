# Outputs

Two main kinds of output can be produced by the `physprep` workflow: raw physiolocal data
organized in BIDS format (according to the BIDS specification v1.10.0), and physiological 
derivatives.

## Raw physiological data 

If the `skip_convert` flag is not specified when running the workflow, all of the output 
will be save at the participant level under `sub-<sub_id>/ses-<ses_id>/<func/eeg>`. If the 
physiological data were acquired concurrently with fmri data, the output will be save 
under the `func` subfolder. If they were acquired concurrently with eeg data, the output 
will be save under the `eeg` subfolder.

### Data and metadata files

For each session (or run) the following files wille be generated:

- `[matches]_physio.tsv.gz` : raw segmented physiological timeseries.
- `[matches]_physio.json` : metadata file containing the tsv columns names, start time, and 
signal sampling frequency information.

## Physiological derivatives

If `derivatives_dir` is specified when running the workflow, all the derivatives will be
save at the participant level in that directory. If not specified, the derivatives will
be save at the participant level under `<indir_bids>/derivatives/physprep`.

### Data files

For each session (or run) the following files wille be generated:

- `[matches]_desc-preproc_physio.tsv.gz` : preprocessed timeseries.
- `[matches]_desc-quality.json` : quality assessment on the preprocessed timeseries.
- `[matches]_desc-physio_events.tsv` : extracted features.