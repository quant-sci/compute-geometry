"""Sphinx configuration for the compute-geometry documentation."""

import os
import sys
from datetime import date

# Make the cgeom package importable for autodoc.
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------

project = "compute-geometry"
author = "Kleyton da Costa"
copyright = f"{date.today().year}, {author}"

# Pull the version straight from the installed package metadata.
try:
    from importlib.metadata import version as _version

    release = _version("compute-geometry")
except Exception:  # pragma: no cover - fallback when not installed
    release = "0.1.2"
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "myst_parser",
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- Autodoc / autosummary ---------------------------------------------------

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# Napoleon understands the Google-style docstrings used across the codebase.
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_use_rtype = False

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# -- MyST --------------------------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "smartquotes",
]
myst_heading_anchors = 3

# -- HTML output (Furo theme) ------------------------------------------------

html_theme = "furo"
html_title = "compute-geometry"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "../public/logo.svg"
html_favicon = "../public/logo.svg"

html_theme_options = {
    "sidebar_hide_name": True,
    "navigation_with_keys": True,
    "source_repository": "https://github.com/kleyt0n/compute-geometry",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#0466c8",
        "color-brand-content": "#0466c8",
    },
    "dark_css_variables": {
        "color-brand-primary": "#5b8def",
        "color-brand-content": "#5b8def",
    },
}
