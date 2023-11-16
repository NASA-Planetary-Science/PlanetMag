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

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinxcontrib.apidoc',
              'sphinxcontrib.matlab',
              'myst_parser']
source_suffix = ['.rst', '.md']

_HERE = os.path.dirname(__file__)
_ROOT_DIR = os.path.abspath(os.path.join(_HERE, '..'))
#_PACKAGE_DIR = os.path.abspath(os.path.join(_HERE, '../PlanetMag'))  # Not needed until Python conversion

sys.path.insert(0, _ROOT_DIR)
#sys.path.insert(0, _PACKAGE_DIR)

matlab_src_dir = _ROOT_DIR
matlab_auto_link = 'all'

apidoc_module_dir = _ROOT_DIR
apidoc_output_dir = 'stubs'
apidoc_template_dir = 'templates'
apidoc_excluded_paths = ['config*', 'setup.py']
apidoc_separate_modules = True
apidoc_module_first = True

templates_path = ['templates']
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
