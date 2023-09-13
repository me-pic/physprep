"""
Physprep workflow.

Preprocess raw physiological data acquired in MRI, extract features, and generate quality
report.
"""

import click

from physprep import utils
from physprep.prepare import convert, get_info, match_acq_bids, rename

# from physprep.processing import clean, process
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
    "indir_mri",
    type=click.Path(),
    help="Path to the directory containing the raw MRI data.",
)
@click.argument(
    "indir_physio",
    type=click.Path(),
    help="Path to the directory containing the raw physiological data.",
)
@click.argument(
    "outdir",
    type=click.Path(),
    help="Path to the directory where the processed physiological data will be saved.",
)
@click.argument("sub", type=str, help="Subject label.")
@click.option("--ses", type=str, default=None, required=False, help="Session label.")
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
def main(
    workflow_strategy,
    indir_mri,
    indir_physio,
    outdir,
    sub,
    ses=None,
    skip_match_acq_bids=False,
    skip_convert=False,
):
    """
    Physprep workflow.

    Preprocess raw physiological data acquired in MRI, extract features, and generate
    quality report.

    Parameters
    ----------
    indir_mri : str or pathlib.Path
        Path to the directory containing the raw MRI data.
    indir_physio : str or pathlib.Path
        Path to the directory containing the raw physiological data.
    outdir : str or pathlib.Path
        Path to the directory where the processed physiological data will be saved.
    sub : str
        Subject label.
    ses : str, optional
        Session label, by default `None`.
    workflow_strategy : str or pathlib.Path
        Name of the workflow_strategy if using a preset. It is also possible to use a
        custom file by providing the path to a JSON file containing workflow strategy.
        In that case, please check Physprep documentation to make sure your file is
        properly formatted.
    """

    # Get workflow info as defined in the configuration file `workflow_strategy`
    workflow_strategy = utils.get_workflow_strategy(workflow_strategy)

    # Match acq files with bold files if specified
    if not skip_match_acq_bids:
        match_acq_bids(indir_mri, indir_physio)
    if not skip_convert:
        # Get information about the physiological data
        get_info()
        # Convert physiological data to BIDS format with phys2bids
        convert()
        # Rename physiological data to BIDS format
        rename()
    # Clean physiological data

    # Process physiological data

    # Generate quality report

    pass
