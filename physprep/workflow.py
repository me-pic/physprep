"""
Physprep workflow.

Preprocess raw physiological data acquired in MRI, extract features, and generate quality
report.
"""

from pathlib import Path

import click
import pandas as pd

from physprep import utils
from physprep.prepare import convert, get_info, match_acq_bids, rename
from physprep.processing import clean, process
from physprep.quality import report


@click.command()
@click.argument(
    "workflow_strategy",
    type=click.Path(),
)
@click.argument(
    "indir_bids",
    type=click.Path(),
)
@click.argument("sub", type=str)
@click.option("--ses", type=str, default=None, required=False, help="Session label.")
@click.option(
    "--indir_raw_physio",
    type=click.Path(),
    default=None,
    help="Path to the directory containing the raw physiological data. Specify if raw "
    "physiological data is not in the BIDS directory. For more details, about the BIDS "
    "data structure, please refer to the documentation.",
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
    "--padding",
    type=int,
    default=9,
    help="Time (in seconds) of padding to add at the beginning and end of each run. "
    "Default to 9.",
)
def main(
    workflow_strategy,
    indir_bids,
    sub,
    ses=None,
    indir_raw_physio=None,
    skip_match_acq_bids=False,
    skip_convert=False,
    padding=9,
    n_jobs=None,
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

    outdir : str or pathlib.Path
        Path to the directory where the processed physiological data will be saved.

    sub : str

        Subject label.

    ses : str, optional

        Session label, by default `None`.

    indir_raw_physio : str or pathlib.Path

        Path to the directory containing the raw physiological data. Specify if raw
        physiological data is not in the BIDS directory.

    skip_match_acq_bids : bool, optional

        If specified, the workflow will not match the acq files with the bold files. Use
        if acq files are already organized properly.

    skip_convert : bool, optional

        If specified, the workflow will not convert the physiological data recordings
        in BIDS format.

    padding : int, optional

        Time (in seconds) of padding to add at the beginning and end of each run. This
        parameter is used if `skip_convert` is set to False. See Phys2BIDS documentation
        for more details. By default, 9.
    """
    # TODO: add dataset_description.json file in derivatives/ directory
    # Set up directories
    # Check if directories exist
    indir_bids = Path(indir_bids)
    if not indir_bids.exists():
        raise FileNotFoundError(f"{indir_bids} does not exist.")
    if indir_raw_physio is not None:
        indir_raw_physio = Path(indir_raw_physio)
        if not indir_raw_physio.exists():
            raise FileNotFoundError(f"{indir_raw_physio} does not exist.")
    # Create output directories
    if ses is not None and not isinstance(ses, list):
        ls_ses = [ses]
        raw_dir = indir_bids / "sourcedata" / sub
        segmented_dir = indir_bids / sub
        derivatives_dir = indir_bids / "derivatives" / "physprep" / sub
    elif ses is None:
        ls_ses = sorted(Path(indir_bids / sub).glob("ses-*"))
        # If no ses-* subdirectory in sub
        if len(ls_ses) == 0:
            raw_dir = indir_bids / "sourcedata" / sub / "func"
            segmented_dir = indir_bids / sub / "func"
            derivatives_dir = indir_bids / "derivatives" / "physprep" / sub
        # If ses-* subdirectories in sub
        else:
            ses = ls_ses
            raw_dir = indir_bids / "sourcedata" / sub
            segmented_dir = indir_bids / sub
            derivatives_dir = indir_bids / "derivatives" / "physprep" / sub

    raw_dir.mkdir(parents=True, exist_ok=True)
    segmented_dir.mkdir(parents=True, exist_ok=True)
    derivatives_dir.mkdir(parents=True, exist_ok=True)

    # Get workflow info as defined in the configuration file `workflow_strategy`
    workflow = utils.get_config(workflow_strategy, strategy="workflow")

    # Match acq files with bold files if specified
    if not skip_match_acq_bids:
        match_acq_bids(indir_bids, indir_raw_physio)
    if not skip_convert:
        # Get information about the physiological recordings
        info_sessions = get_info.get_info(
            indir_bids,
            sub,
            ses,
            count_vol=True,
            save=indir_bids / "sourcedata",
            tr_channel=workflow["trigger"]["channel"],
        )
        # Convert physiological data to BIDS format with phys2bids
        # Session-level
        if len(ls_ses) == 1:
            convert.convert(
                raw_dir, segmented_dir, sub, ses=ses, info=info_sessions, pad=padding
            )
        # Subject-level: if ses-*, conversion will be done for all sessions
        else:
            convert.convert(raw_dir, segmented_dir, sub, info=info_sessions, pad=padding)
        # Rename physiological data to BIDS format
        rename.co_register_physio(segmented_dir, sub, ses=ses)

    # Clean & process physiological data
    if len(ls_ses) >= 1:
        for s in ls_ses:
            runs = sorted(s.glob("func/*_physio.*"))
            # Remove duplicated elements in runs with same filename but different
            # extension
            runs = list(set([run.parent / run.stem for run in runs]))
            # Need to run it twice because of the tsv.gz extension
            runs = list(set([run.parent / run.stem for run in runs]))
            runs.sort()

            for run in runs:
                filename = run.stem
                print(f"\nLoading {filename}...\n")
                # Load data
                metadata = utils.load_json(run.with_suffix(".json"))
                data = pd.read_csv(
                    run.with_suffix(".tsv.gz"), sep="\t", names=metadata["Columns"]
                )
                print("Data loaded.\n")
                print("Preprocessing data...\n")
                # Preprocess data
                preprocessed_signals, metadata_derivatives = clean.preprocessing_workflow(
                    data, metadata, workflow, Path(derivatives_dir / s.stem), filename
                )
                print("Preprocessing done.\n")
                print("Extracting features...\n")
                # Extract features
                timeseries, features = process.features_extraction_workflow(
                    preprocessed_signals,
                    metadata_derivatives,
                    workflow,
                    Path(derivatives_dir / s.stem),
                    filename,
                )
                print("Features extracted.\n")
                print("Generating quality report...\n")
                # Generate quality report
                report.computing_sqi(segmented_dir.parent, derivatives_dir.parent, sub, s)
                print("Quality report generated.\n")
    else:
        pass


if __name__ == "__main__":
    main()
