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

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="25%"><a href="https://github.com/me-pic"><img src="https://avatars.githubusercontent.com/u/77584086?v=4?s=200" width="200px;" alt="Marie-Eve Picard"/><br /><sub><b>Marie-Eve Picard</b></sub></a><br /><a href="#ideas-me-pic" title="Ideas, Planning, & Feedback">🤔</a> <a href="#code-me-pic" title="Code">💻</a> <a href="#test-me-pic" title="Tests">⚠️</a><a href="#maintenance-me-pic" title="Maintenance">🚧</a></td>
      <td align="center" valign="top" width="25%"><a href="https://github.com/sangfrois"><img src="https://avatars.githubusercontent.com/u/38385719?v=4?s=200" width="200px;" alt="François Lespinasse"/><br /><sub><b>François Lespinasse</b></sub></a><br /><a href="#ideas-sangfrois" title="Ideas, Planning, & Feedback">🤔</a> <a href="#code-sangfrois" title="Code">💻</a></td>
      <td align="center" valign="top" width="25%"><a href="https://github.com/emdupre"><img src="https://avatars.githubusercontent.com/u/15017191?v=4?s=200" width="200px;" alt="Elizabeth DuPre"/><br /><sub><b>Elizabeth DuPre</b></sub></a><br /><a href="#bug-emdupre" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="25%"><a href="https://github.com/bpinsard"><img src="https://avatars.githubusercontent.com/u/1155388?v=4?s=200" width="200px;" alt="Basile Pinsard"/><br /><sub><b>Basile Pinsard</b></sub></a><br /><a href="#code-bpinsard" title="Code">💻</a></td>
    </tr>
  </tbody>
</table>

## Acknowlegments

Please cite `physprep` using the Zenodo DOI:

```bibtex
@software{PhysPrep,
    author = {Marie-Eve Picard and François Lespinasse and Elizabeth DuPre and Basile Pinsard},
    license = {Apache-2.0},
    title = {{physprep}},
    url = {https://github.com/courtois-neuromod/physprep/},
    doi = {https://doi.org/10.5281/zenodo.15858672}
}
```

We also want to give credits to the developers of the libraries on which PhysPrep mainly relies, namely:

:raised_hands:  [Phys2Bids](https://github.com/physiopy/phys2bids) and the [Physiopy community](https://physiopy.github.io/) <br>

Please consider citing `phys2bids` using the Zenodo DOI:

```bibtex
@software{Phys2bids,
    author = {The phys2bids developers},
    license = {Apache-2.0},
    title = {{phy2bids}},
    url = {https://github.com/physiopy/phys2bids},
    doi = {https://doi.org/10.5281/zenodo.3586045}
}
```

:sparkles:  [NeuroKit2](https://github.com/neuropsychology/NeuroKit) <br>

Please consider citing `neurokit2` using the Zenodo DOI:

```bibtex
@software{Neurokit2,
    author = {Dominique Makowski and Tam Pham and Zen J. Lau and Jan C. Brammer},
    license = {MIT},
    title = {{neurokit2}},
    url = {https://github.com/neuropsychology/NeuroKit/},
    doi = {https://doi.org/10.5281/zenodo.13331945}
}
```
