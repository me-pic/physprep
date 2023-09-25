# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

from physprep import __copyright__, __packagename__

sys.path.insert(0, os.path.abspath(".."))
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = __packagename__
copyright = __copyright__
author = "Marie-Eve Picard"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_rtd_theme",
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_context = {
    "display_github": True,
    "github_user": "me-pic",
    "github_repo": "physprep",
    "github_version": "docs/setup",
    "conf_py_path": "https://github.com/courtois-neuromod/physprep/tree/main/docs/source",
}
html_static_path = ["_static"]

# -- Options for myst_parser -------------------------------------------------
myst_enable_extensions = ["colon_fence"]
