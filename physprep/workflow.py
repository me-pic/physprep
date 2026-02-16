"""
Physprep workflow.

Preprocess raw physiological data acquired in MRI, extract features, and generate quality
report.
"""

from pathlib import Path
from codecarbon import OfflineEmissionsTracker

import sys
import click

from physprep import utils
from physprep.processing import clean, process
from physprep.quality import qa, report

import logging
logger = logging.getLogger(__name__)

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
@click.option(
    "--track_carbon",
    is_flag=True,
    help="If specified, carbon tracker will be used to track power use of the pipeline."
)
@click.option(
    "--country_code",
    type=str,
    default="CAN",
    help="Country ISO code used by carbon trackers. Will only be used if --track_carbon flag is specified."
)
def main(
    workflow_strategy,
    indir_bids,
    sub=None,
    ses=None,
    derivatives_dir=None,
    save_report=False,
    track_carbon=False,
    country_code="CAN"
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
    code_dir = Path(derivatives_dir / "code")
    code_dir.mkdir(parents=True, exist_ok=True)

    log_file = code_dir / "log.txt"

    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
        force=True
    )

    logger.info("Command: %s", " ".join([sys.executable] + sys.argv))

    if track_carbon:
        logger.info("CodeCarbon started")
        logger.info("Using country_iso_code=%s", country_code)
        logger.info("Saving report at: %s", code_dir)

        tracker = OfflineEmissionsTracker(
            output_dir=code_dir, country_iso_code=country_code
        )
        tracker.start()

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
        error_msg = "No files found. Please make sure you specified the correct directory, and if applicable the correct values for `sub` and/or `ses`."
        logger.error(error_msg)
        raise FileNotFoundError(
            error_msg
        )

    for file in files:
        # Load metada
        metadata = file.get_metadata()
        if not bool(metadata):
            error_msg = f"No metadata file associated with {file}."
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        # Load data
        logger.info("\nLoading %s...\n", file)
        data = file.get_df()
        logger.info("Data loaded.\n")


        # Preprocess data
        logger.info("Preprocessing data...\n")
        preprocessed_signals, metadata_derivatives, preproc_strategy = clean.preprocessing_workflow(
            data, metadata, workflow
        )
        logger.info("Saving preprocessed signals...\n")
        utils.save_processing(
            derivatives_dir,
            file.get_entities(),
            "preproc",
            preprocessed_signals,
            metadata_derivatives,
            preproc_strategy
        )
        logger.info("Preprocessing done.\n")

        # Extract features
        logger.info("Extracting features...\n")
        features, events = process.features_extraction_workflow(
            preprocessed_signals, metadata_derivatives, workflow
        )
        logger.info("Saving extracted features...\n")
        utils.save_features(derivatives_dir, file.get_entities(), events)
        logger.info("Features extraction done.\n")

        # Generate quality report
        logger.info("Assessing quality of the data...\n")
        qa_metrics, qa_short = qa.computing_sqi(
            workflow, preprocessed_signals, features, metadata_derivatives
        )
        logger.info("Saving quality assessment...\n")
        utils.save_qa(derivatives_dir, file.get_entities(), qa_metrics)
        utils.save_qa(derivatives_dir, file.get_entities(), qa_short, short=True)
        logger.info("Data quality assessed.\n")

        if save_report:
            logger.info("Generating QC report... \n")
            qa_report = report.generate_report(
                workflow,
                qa_metrics,
                preprocessed_signals,
                features,
                metadata_derivatives,
                derivatives_dir,
                file.get_entities(),
            )
            logger.info("QC report generated. \n")
            utils.save_qa(derivatives_dir, file.get_entities(), qa_report, report=True)
    
    if track_carbon:
        emissions = tracker.stop()
        logger.info("CodeCarbon stopped")
        logger.info("Saving report at %s", code_dir)
        logger.info("Carbon emissions: %s kg", emissions)


if __name__ == "__main__":
    main()
