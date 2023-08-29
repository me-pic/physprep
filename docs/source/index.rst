.. Physprep documentation master file, created by
   sphinx-quickstart on Tue Aug 29 13:20:35 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Physprep: Physiological Data Preparation Workflow Adapted for MR Settings
=========================================================================

About
-----

Peripheral biosignals are valuable for capturing fluctuations in cognitive and affective
states, making it possible to characterize brain-body interactions in neuroimaging setup.
The relevance of these measures when acquiring neuroimaging data and the subsequent
increase in the amount of data processing highlight the importance of establishing
reliable processing and reporting practices for biosignals acquired in magnetic resonance
(MR) context, but a standardized head-to-tail workflow has yet to be proposed. Physprep
is an open-source tool to prepare biosignals for analysis. This workflow presents itself
as user-friendly and flexible by providing customization options through configuration
files. The (pre-)processing of cardiac signals (ECG and PPG), electrodermal signals (EDA),
and respiratory signals is currently supported.

.. note::

   This project is under active development.

Physprep was originally developed to preprocess physiological data acquired concurrently
with large-scale individual functional MRI as part of the `CNeuromMod project
<https://docs.cneuromod.ca/en/latest/>`. This workflow combines well-maintained,
community-based Python libraries such as
`Phys2bids <https://phys2bids.readthedocs.io/en/latest/>` and `NeuroKit2
<https://neuropsychology.github.io/NeuroKit/>`, and relies on the BIDS standard.

Contents
--------

.. toctree::
   :maxdepth: 1

   installation.md
   usage.md
   workflow.md
   outputs.md
   authors.md
   api.rst

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
