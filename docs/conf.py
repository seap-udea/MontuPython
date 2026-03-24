# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "MontuPython"
copyright = "2023-present, Jorge I. Zuluaga"
author = "Jorge I. Zuluaga"
import re as _re
_ver_file = os.path.join(os.path.dirname(__file__), "..", "montu", "version.py")
release = _re.search(r"version\s*=\s*['\"]([^'\"]+)['\"]", open(_ver_file).read()).group(1)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "nbsphinx",
    "sphinx_rtd_theme",
    "sphinx.ext.mathjax",
]

# MyST Parser configuration
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    "html_image",
    "colon_fence",
]

napoleon_numpy_docstring = True
napoleon_google_docstring = False

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
pygment_style = "sphinx"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "collapse_navigation": False,
    "titles_only": False,
    "external_links": [
        {"name": "Source", "url": "https://github.com/seap-udea/MontuPython"},
    ],
}

html_js_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js",
    "https://unpkg.com/@jupyter-widgets/html-manager@^0.20.0/dist/embed-amd.js",
    "https://cdn.plot.ly/plotly-latest.min.js",
]
