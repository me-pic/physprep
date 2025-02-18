"""
Physprep workflow.

Preprocess raw physiological data acquired in MRI, extract features, and generate quality
report.
"""

from pathlib import Path

import click

from physprep import utils
from physprep.prepare import convert, get_info, match_acq_bids, rename
from physprep.processing import clean, process
from physprep.quality import report, qa


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
    "--indir_raw_physio",
    type=click.Path(),
    default=None,
    help="Path to the directory containing the raw physiological data. Specify if raw "
    "physiological data is not in the BIDS directory. For more details, about the BIDS "
    "data structure, please refer to the documentation.",
)
@click.option(
    "--derivatives_dir",
    type=click.Path(),
    default=None,
    help="Path to the derivatives directory.",
)
@click.option(
    "--skip_match_acq_bids",
    is_flag=True,
    help="If specified, the workflow will not match the acq files with the bold files. "
    "If acq files are already organized properly, this flag can be specified. For more "
    "details, see the documentation of the mathc_acq_bids.py script.",
)
@click.option(
    "--skip_convert",
    is_flag=True,
    help="If specified, the workflow will not convert the physiological data recordings "
    "in BIDS format. This implies that the data is already segmented in runs, organized "
    "in a BIDS-like structure (i.e., one tsv.gz and one json file per run), and named "
    "following the BIDS recommandations.",
)
@click.option(
    "--heur",
    type=str,
    help="File needed to convert raw data into BIDS format. For more details, check the "
    "phys2bids documentation.",
)
@click.option(
    "--padding",
    type=int,
    default=9,
    help="Time (in seconds) of padding to add at the beginning and end of each run. "
    "Default to 9.",
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
    indir_raw_physio=None,
    derivatives_dir=None,
    skip_match_acq_bids=False,
    skip_convert=False,
    heur = None,
    padding=9,
    save_report=False
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

    skip_match_acq_bids : bool, optional

        If specified, the workflow will not match the acq files with the bold files. Use
        if acq files are already organized properly.

    skip_convert : bool, optional

        If specified, the workflow will not convert the physiological data recordings
        in BIDS format.

    heur : str, optional

        File needed to convert raw data into BIDS format. For more details, check the
        phys2bids documentation.

    padding : int, optional

        Time (in seconds) of padding to add at the beginning and end of each run. This
        parameter is used if `skip_convert` is set to False. See Phys2BIDS documentation
        for more details. By default, 9.

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

    if indir_raw_physio is not None:
        indir_raw_physio = Path(indir_raw_physio)
        if not indir_raw_physio.exists():
            raise FileNotFoundError(f"{indir_raw_physio} does not exist.")


    # Create output directories
    #segmented_dir = indir_bids / sub
    #segmented_dir.mkdir(parents=True, exist_ok=True)
    if derivatives_dir is None:
        derivatives_dir = indir_bids / 'derivatives' / 'physprep' 
    else : 
        derivatives_dir = Path(derivatives_dir)
    derivatives_dir.mkdir(parents=True, exist_ok=True)

    # Get workflow info as defined in the configuration file `workflow_strategy`
    workflow = utils.get_config(workflow_strategy, strategy="workflow")
    
    # Match acq files with bold files if specified
    if not skip_match_acq_bids:
        match_acq_bids(indir_bids, indir_raw_physio)
    """
    if not skip_convert:
        # Get information about the physiological recordings
        info_sessions = get_info.get_info(
            indir_bids,
            sub,
            ses,
            count_vol=True,
            save=indir_bids / 'sourcedata',
            tr_channel=workflow['trigger']['channel'],
        )
        # Convert physiological data to BIDS format with phys2bids
        # Session-level
        ls_ses = layout.get_sessions()
        if ses is not None :
            ls_ses = list(ses)
            
        if len(ls_ses) == 1:
            convert.convert(
                raw_dir, segmented_dir, sub, ses=ses, info=info_sessions, pad=padding
            )
        # Subject-level: if ses-*, conversion will be done for all sessions
        else:
            convert.convert(raw_dir, segmented_dir, sub, info=info_sessions, pad=padding)
        # Rename physiological data to BIDS format
        rename.co_register_physio(segmented_dir, sub, ses=ses)
    """
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
            workflow,
            preprocessed_signals,
            features,
            metadata_derivatives
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
                file.get_entities()
            )
            print("QC report generated. \n")
            utils.save_qa(derivatives_dir, file.get_entities(), qa_report, report=True)


if __name__ == "__main__":
    main()
