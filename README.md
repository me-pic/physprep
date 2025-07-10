# PhysPrep workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15664327.svg)](https://doi.org/10.5281/zenodo.15664327)
[![Static Badge](https://img.shields.io/badge/Apache--2.0-license?style=flat&label=license&color=green)](https://github.com/courtois-neuromod/physprep/blob/main/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/physprep/badge/?version=latest)](https://physprep.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Description

Peripheral biosignals are valuable for capturing fluctuations in cognitive and affective states, making it possible to characterize brain-body interactions in neuroimaging setup. The relevance of these measures when acquiring neuroimaging data and the subsequent increase in the amount of data processing highlight the importance of establishing reliable processing and reporting practices for biosignals acquired in magnetic resonance (MR) context, but a standardized head-to-tail workflow has yet to be proposed. Physprep is an open-source tool to prepare biosignals for analysis. This workflow presents itself as user-friendly and flexible by providing customization options through configuration files. The (pre-)processing of cardiac signals (ECG and PPG), electrodermal signals (EDA), and respiratory signals is currently supported.

## Quick start

Clone the project and create a virtual environment with Python >=3.9 and <=3.11 !
```
git clone git@github.com:courtois-neuromod/physprep.git

python3 -m venv <venv_name>
source <venv_name>/bin/activate

cd physprep
pip install .
```

## Acknowlegments

We want to give credits to the developers of the libraries on which PhysPrep mainly relies, namely:

&ensp; &nbsp; &nbsp; &nbsp; :raised_hands:  [Phys2Bids](https://github.com/physiopy/phys2bids) and the [Physiopy community](https://physiopy.github.io/) <br>
&ensp; &nbsp; &nbsp; &nbsp; :sparkles:  [NeuroKit2](https://github.com/neuropsychology/NeuroKit) <br>
