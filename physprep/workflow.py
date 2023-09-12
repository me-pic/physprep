"""
Physprep workflow.

Preprocess raw physiological data acquired in MRI, extract features, and generate quality
report.
"""

import click

# from physprep import utils
# from physprep.prepare import convert, get_info, rename
# from physprep.processing import clean, process
# from physprep.quality import report


@click.command()
@click.argument(
    "workflow_strategy",
    type=click.Path(),
    help="Path to a JSON file containing the workflow strategy.",
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
def main(workflow_strategy, indir_mri, indir_physio, outdir, sub, ses=None):
    """
    Physprep workflow.

    Preprocess raw physiological data acquired in MRI, extract features, and generate
    quality report.

    Parameters
    ----------
    workflow_strategy : str or pathlib.Path
        Path to a JSON file containing the workflow strategy.
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
    """
    # Get information about the physiological data

    # Convert physiological data to BIDS format with phys2bids

    # Rename physiological data to BIDS format

    # Clean physiological data

    # Process physiological data

    # Generate quality report

    pass
