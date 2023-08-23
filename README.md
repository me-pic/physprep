# Physprep workflow

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

This standalone repository acts as use-case documentation for physiological data processing steps. The proposed workflow integrates community-based Ptyhon librairies such as [phys2bids](https://github.com/physiopy/phys2bids) and [neurokit2](https://github.com/neuropsychology/NeuroKit). 

The repo is separated in three main modules, and provides a setp-by-step tutorial for each of them:

`utils\` 
1. `list_sub.py`: list all the physiological files for a given subject (and a given session).
2. `get_info.py`: retrieve physiological files information.
3. `match_acq_bids.py`: match Acqknowledge files (.acq) with the fMRI Nifti files (.nii.gz).

`processing\`
1. `convert.py`: use [phys2bids](https://github.com/physiopy/phys2bids) to segment the acqknowledge files in runs following the BIDS format.
2. `clean.py`: implement functions to filter the physiological signals, and to remove the artifacts induced by the MRI.
3. `process.py`: build a processing pipeline based on `clean.py` functions.

`quality\`
1. `report.py`: provide a summary of the quality of the processed signal.
2. `time.py`: provide quantitative metrics computed on the time domain
3. `frequency.py`: provide quantitative metrics computed on the frequency domain

`visu\` :construction_worker:

## Acqknowlegments
Thanks to the generous data donation from a subject of the [Courtois-Neuromod project](https://www.cneuromod.ca/), research communities like [PhysioPy](https://physiopy.github.io/) will benefit from common data access to test and optimize their physio data preparation workflows, using [BIDS](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/06-physiological-and-other-continuous-recordings.html) format. 
