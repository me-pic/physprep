"""
Physprep workflow.

Preprocess raw physiological data acquired in MRI, extract features, and generate quality
report.
"""

from pathlib import Path

import click

from physprep import utils
from physprep.processing import clean, process
from physprep.quality import qa, report


@click.command()
@click.argument(
    "workflow_strategy",
    type=click.Path(),
)
@click.argument(
    "indir_bids",
    type=click.Path(),
)
@click.option(
    "--sub",
    type=str,
    default=None,
    required=False,
    help="Subject id. Use only to process the data of that specific subject. "
    "For example: if you specify --sub 01, only sub-01 data will be processed.",
)
@click.option(
    "--ses",
    type=str,
    default=None,
    required=False,
    help="Session label. Use only to process the data of that specific session. "
    "For example: if you specify --ses 001, only data from ses-001 will be "
    "processed. If specify, but --sub not specified, data from the specified "
    "session (e.g. ses-001) across all subjects will be processed.",
)
@click.option(
    "--derivatives_dir",
    type=click.Path(),
    default=None,
    help="Path to the derivatives directory.",
)
@click.option(
    "--save_report",
    is_flag=True,
    help="If specified, an quality report will be generated and saved for each run.",
)
def main(
    workflow_strategy,
    indir_bids,
    sub=None,
    ses=None,
    derivatives_dir=None,
    save_report=False,
):
    """Physprep workflow.

    Preprocess raw physiological data acquired in MRI, extract features, and generate
    quality report.
    \b

    Parameters

    ----------

    workflow_strategy : str or pathlib.Path

        Name of the workflow_strategy if using a preset. It is also possible to use a
        custom file by providing the path to a JSON file containing workflow strategy.
        In that case, please check Physprep documentation to make sure your file is
        properly formatted.

    indir_bids : str or pathlib.Path

        Path to the directory containing the BIDS-like dataset.

    sub : str, optional

        Subject id. E.g. '01'

    ses : str, optional

        Session id, by default `None`. E.g. '001'

    indir_raw_physio : str or pathlib.Path

        Path to the directory containing the raw physiological data. Specify if raw
        physiological data is not in the BIDS directory.

    derivatives_dir : str or pathlib.Path

        Path to the output directory. If `None`, a `derivatives/` directory will be
        created in the root of the directory specified by `indir_bids` (i.e., source BIDS
        dataset), as per [BIDS specification](https://bids-specification.readthedocs.io/en/stable/common-principles.html#storage-of-derived-datasets).

    save_report : bool, optional

        If specified, an quality report will be generated and saved for each run.

    """
    # Set up directories
    # Check if directories exist
    indir_bids = Path(indir_bids)
    if not indir_bids.exists():
        raise FileNotFoundError(f"{indir_bids} does not exist.")
    # Check if indir_bids is already a bids dataset
    layout = utils._check_bids_validity(indir_bids)

    # Create output directories
    if derivatives_dir is None:
        derivatives_dir = indir_bids / "derivatives" / "physprep"
    else:
        derivatives_dir = Path(derivatives_dir)
    derivatives_dir.mkdir(parents=True, exist_ok=True)

    # Get workflow info as defined in the configuration file `workflow_strategy`
    workflow = utils.get_config(workflow_strategy, strategy="workflow")

    # Clean & process physiological data
    # Change to iterate through files and not sessions + files by using BIDSLayout:
    # Defines parameters to get the physio files
    info_layout = {"extension": "tsv.gz", "suffix": "physio"}
    if sub is not None:
        info_layout["subject"] = utils._check_sub_validity(sub, layout.get_subjects())
    if ses is not None:
        info_layout["session"] = utils._check_ses_validity(ses, layout.get_sessions())
    # Get directory for files containing physio timeseries
    files = layout.get(**info_layout)
    # Make sure `files` is not an empty list
    if len(files) == 0:
        raise FileNotFoundError(
            "No files found. Please make sure you specified the correct directory, "
            "and if applicable the correct values for `sub` and/or `ses`."
        )

    for file in files:
        # Load metada
        metadata = file.get_metadata()
        if not bool(metadata):
            raise FileNotFoundError(f"No metadata file associated with {file}.")
        # Load data
        print(f"\nLoading {file}...\n")
        data = file.get_df()
        print("Data loaded.\n")

        # Preprocess data
        print("Preprocessing data...\n")
        preprocessed_signals, metadata_derivatives = clean.preprocessing_workflow(
            data, metadata, workflow
        )
        print("Saving preprocessed signals...\n")
        utils.save_processing(
            derivatives_dir,
            file.get_entities(),
            "preproc",
            preprocessed_signals,
            metadata_derivatives,
        )
        print("Preprocessing done.\n")

        # Extract features
        print("Extracting features...\n")
        features, events = process.features_extraction_workflow(
            preprocessed_signals, metadata_derivatives, workflow
        )
        print("Saving extracted features...\n")
        utils.save_features(derivatives_dir, file.get_entities(), events)
        print("Features extraction done.\n")

        # Generate quality report
        print("Assessing quality of the data...\n")
        qa_metrics, qa_short = qa.computing_sqi(
            workflow, preprocessed_signals, features, metadata_derivatives
        )
        print("Saving quality assessment...\n")
        utils.save_qa(derivatives_dir, file.get_entities(), qa_metrics)
        utils.save_qa(derivatives_dir, file.get_entities(), qa_short, short=True)
        print("Data quality assessed.\n")

        if save_report:
            print("Generating QC report... \n")
            qa_report = report.generate_report(
                workflow,
                qa_metrics,
                preprocessed_signals,
                features,
                metadata_derivatives,
                derivatives_dir,
                file.get_entities(),
            )
            print("QC report generated. \n")
            utils.save_qa(derivatives_dir, file.get_entities(), qa_report, report=True)


if __name__ == "__main__":
    main()
