# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path
import grobnercrystals
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'grobnercrystals'
copyright = '2025, Abigail Price'
author = 'Abigail Price'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.githubpages',
    ]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ["custom.css"]