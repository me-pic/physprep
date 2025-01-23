import datetime
import logging
import os

import bioread
import click
import pandas
from pathlib import Path
from pytz import timezone


@click.command()
@click.argument("bids_path", type=click.Path())
@click.argument("biopac_path", type=click.Path())
@click.option("--overwrite", is_flag=True)
def match_all_bolds(bids_path, biopac_path, overwrite=True):
    """
    Match the Acqknowldge files (.acq) with the bold files (.nii.gz).
    The correspondence between the acq and the nii files is saved in
    a tsv file.

    Parameters
    ----------
    bids_path : str
        BIDS data directory.
    biopac_path : str
        Biopac data directory.
    """
    bids_path = Path(bids_path)
    biopac_path = Path(biopac_path)
    tz = timezone("Canada/Eastern")
    # get acq file starts and end datetimes
    acqk_files = sorted(list(biopac_path.glob('*.acq')))
    acqk_files_startends = []
    for acqk in acqk_files:
        try:
            acq_h = bioread.read_headers(str(acqk))
            acq_start = acq_h.earliest_marker_created_at
            if acq_start is None:
                logging.error(f"no start marker in: {acqk}")
                continue
            acq_end = acq_start
            if len(acq_h.time_index):
                acq_end = acq_start + datetime.timedelta(seconds=acq_h.time_index[-1])
            acqk_files_startends.append((acqk, acq_start, acq_end))
        except Exception:
            logging.error(f"read error for file: {acqk}")

    sourcedata = bids_path / "sourcedata" / "physio"
    sourcedata.mkdir(parents=True, exist_ok=True)
    sessions_list = sorted(bids_path.glob("sub-*/ses-*"))
    for session in sessions_list:
        session_bids_path = session.relative_to(bids_path)
        session_sourcedata = sourcedata / session_bids_path
        session_sourcedata.mkdir(parents=True, exist_ok=True)
        sub_ses_prefix = str(session_bids_path).replace(os.path.sep, "_")

        try:
            scans = pandas.read_csv(
                session / (sub_ses_prefix + "_scans.tsv"),
                delimiter="\t",
                parse_dates=["acq_time"],
            )
        except Exception as e:
            logging.error(f"read error for file: {session / (sub_ses_prefix + '_scans.tsv')}")
            raise ValueError(f"Please verify that {session / (sub_ses_prefix + '_scans.tsv')} exits")

            
        list_matches_out = session_sourcedata / (
            sub_ses_prefix + "_physio_fmri_matches.tsv"
        )
        breakpoint()
        if list_matches_out.exists() and not overwrite:
            continue
        matches = []
        for idx, scan in scans.iterrows():
            acq_files = []
            if any(suffix in scan.filename for suffix in ["bold.nii.gz", "eeg.edf"]):
                acq_time = tz.localize(scan.acq_time.to_pydatetime())

                for acqk_wtiming in acqk_files_startends:
                    # Takes into account that the start time of the bold or eeg 
                    # acquisition is somewhere between the start and end of the acq 
                    # recording
                    if acqk_wtiming[1] < acq_time and acqk_wtiming[2] > acq_time:
                        acq_files.append(acqk_wtiming)
                    # Try to see if the acq recording started after the bold or eeg 
                    # acquisition
                    else:
                        # Get the end time of `scan` file for eeg data
                        if 'eeg.edf' in scan.filename:
                            from mne.io import read_raw_edf

                            tmp = read_raw_edf(session / scan.filename)
                            tmp_date = pandas.Timestamp(tmp.info['meas_date'])
                            acq_time_end = tmp_date.tz_convert(tz) + datetime.timedelta(seconds=tmp.times[-1])
                            # Takes into account that the start time of the acq recording
                            # is somewhere between the start and end of the eeg acquisition
                            if acq_time < acqk_wtiming[1] and acq_time_end > acqk_wtiming[1]:
                                acq_files.append(acqk_wtiming)

                if len(acq_files) == 0:
                    logging.error(f"No acq file found for: {scan.filename}")
                    matches.append((session / scan.filename, None))
                else:
                    if len(acq_files) > 1:
                        if not all(
                            [acq[1] == acq_files[0][1] for acq in acq_files[1:]]
                        ):  # duplicated files
                            logging.warning(
                                "More that one acq file found for: "
                                f"{scan.filename} \n {acq_files}"
                            )
                    bname = os.path.basename(acq_files[0][0])
                    dest_path = session_sourcedata / bname
                    matches.append(
                        (
                            session / scan.filename,
                            dest_path.relative_to(bids_path),
                        )
                    )
                    if not dest_path.exists() and acq_files[0][0].exists():
                        logging.info(f"moving {acq_files[0][0]} to {dest_path}")
                        os.rename(acq_files[0][0], dest_path)
        list_matches_out.write_text("\n".join([f"{m[0]}\t{m[1]}" for m in matches]))


if __name__ == "__main__":
    logging.basicConfig(level=logging.getLevelName('INFO'))
    match_all_bolds()
