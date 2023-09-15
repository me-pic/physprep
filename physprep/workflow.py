"""
Physprep workflow.

Preprocess raw physiological data acquired in MRI, extract features, and generate quality
report.
"""

from pathlib import Path

import click

from physprep import utils
from physprep.prepare import convert, get_info, match_acq_bids, rename
from physprep.processing import clean  # , process

# from physprep.quality import report


@click.command()
@click.option(
    "workflow_strategy",
    type=click.Path(),
    help="Name of the workflow_strategy if using preset (Ex: 'neuromod'). For all "
    "available presets, please check Physprep documentation. Otherwise, path to custom "
    "workflow strategy file (Ex: '/path/to/my_config_file.json').",
)
@click.argument(
    "indir_bids",
    type=click.Path(),
    help="Path to the bids-like dataset directory.",
)
@click.argument(
    "outdir",
    type=click.Path(),
    help="Path to the directory where the processed physiological data will be saved.",
)
@click.argument("sub", type=str, help="Subject label.")
@click.option("--ses", type=str, default=None, required=False, help="Session label.")
@click.option(
    "indir_raw_physio",
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
    outdir,
    sub,
    ses=None,
    indir_raw_physio=None,
    skip_match_acq_bids=False,
    skip_convert=False,
    padding=9,
):
    """
    Physprep workflow.

    Preprocess raw physiological data acquired in MRI, extract features, and generate
    quality report.

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
    if ses is not None:
        ls_ses = [ses]
        raw_dir = indir_bids / "sourcedata" / sub / ses / "func"
        segmented_dir = indir_bids / sub / ses / "func"
        derivatives_dir = indir_bids / "derivatives" / "physprep" / sub / ses
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
            raw_dir = indir_bids / "sourcedata"
            segmented_dir = indir_bids
            derivatives_dir = indir_bids / "derivatives" / "physprep" / sub

    raw_dir.mkdir(parents=True, exists_ok=True)
    segmented_dir.mkdir(parents=True, exists_ok=True)
    derivatives_dir.mkdir(parents=True, exists_ok=True)

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
    for signal in workflow:
        if signal != "trigger":
            if "preprocessing_strategy" in signal and signal[
                "preprocessing_strategy"
            ] not in ["", " ", None]:
                preprocessing = utils.get_config(
                    signal["preprocessing_strategy"], strategy="preprocessing"
                )
                raw_signal, clean_signal = clean.preprocessing_workflow(
                    preprocessing, segmented_dir, sub, ses=ses
                )
                print(raw_signal, clean_signal)
                pass  # TODO
    # Generate quality report
