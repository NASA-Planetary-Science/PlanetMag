# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

project = 'PlanetMag'
copyright = '2023, California Institute of Technology'
author = 'Corey J. Cochrane and Marshall J. Styczinski'
release = 'v1.0.0'
github_url = 'https://github.com/coreyjcochrane/planetmag'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinxcontrib.matlab',
              'sphinx.ext.autodoc',
              'myst_parser',
              'sphinx.ext.napoleon']
source_suffix = ['.rst', '.md']
matlab_src_dir = os.path.abspath('..')
sys.path.insert(0, os.path.abspath('../'))
autodoc_member_order = 'alphabetical'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'  # Install with pip install sphinx-rtd-theme
html_static_path = ['_static']
html_css_files = ['css/custom.css']
html_logo = '../misc/PlanetMag_logoDocs.png'
html_favicon = '../misc/PlanetMag_logo.ico'

html_theme_options = {
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': '#2980B9',  # Default is #2980B9
    'logo_only': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': -1
}
