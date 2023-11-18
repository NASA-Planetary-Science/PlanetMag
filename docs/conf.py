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
              'sphinx.ext.mathjax',
              'myst_parser']
source_suffix = ['.rst', '.md']
mathjaxVer = 2

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


# -- Options for LaTeX math formatting-----------
# https://www.sphinx-doc.org/en/master/latex.html

imgmath_image_format = 'svg'
imgmath_use_preview = True
imgmath_latex_preamble = r'\usepackage[notextcomp]{stix}' + \
                         r'\usepackage[version=4]{mhchem}' + \
                         r'\usepackage{siunitx}' + \
                         r'\usepackage{upgreek}' + \
                         r'\sisetup{group-separator={\,}, group-minimum-digits={5}, group-digits={integer}}'
latex_elements = {
    'extrapackages': r'\usepackage[notextcomp]{stix}' +
                     r'\usepackage[version=4]{mhchem}' +
                     r'\usepackage{siunitx}' +
                     r'\usepackage{upgreek}',
    'preamble': r'\sisetup{group-separator={\,}, group-minimum-digits={5}, group-digits={integer}}'
}

latex_packages = ['upgreek', 'mhchem']
if mathjaxVer == 3:
    mathjax_path='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js'
    mathjax3_config = {
      'loader': {'load': [f'[tex]/{pkg}' for pkg in latex_packages]},
      'tex': {'packages': {'[+]': [latex_packages]}}
    }
else:
    mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML'
    latex_packages.remove('upgreek')
    latex_packages += ['unicode']
    upgreekDefs = r"""
\def\upalpha{{\unicode{945}}}
\def\upbeta{{\unicode{946}}}
\def\upgamma{{\unicode{947}}}
\def\updelta{{\unicode{948}}}
\def\upepsilon{{\unicode{949}}}
\def\upzeta{{\unicode{950}}}
\def\upeta{{\unicode{951}}}
\def\uptheta{{\unicode{952}}}
\def\upiota{{\unicode{953}}}
\def\upkappa{{\unicode{954}}}
\def\uplambda{{\unicode{955}}}
\def\upmu{{\unicode{956}}}
\def\upnu{{\unicode{957}}}
\def\upxi{{\unicode{958}}}
\def\upomicron{{\unicode{959}}}
\def\uppi{{\unicode{960}}}
\def\uprho{{\unicode{961}}}
\def\upsigma{{\unicode{963}}}
\def\uptau{{\unicode{964}}}
\def\upupsilon{{\unicode{965}}}
\def\upphi{{\unicode{966}}}
\def\upchi{{\unicode{967}}}
\def\uppsi{{\unicode{968}}}
\def\upomega{{\unicode{969}}}
\def\upAlpha{{\unicode{913}}}
\def\upBeta{{\unicode{914}}}
\def\upGamma{{\unicode{915}}}
\def\upDelta{{\unicode{916}}}
\def\upEpsilon{{\unicode{917}}}
\def\upZeta{{\unicode{918}}}
\def\upEta{{\unicode{919}}}
\def\upTheta{{\unicode{920}}}
\def\upIota{{\unicode{921}}}
\def\upKappa{{\unicode{922}}}
\def\upLambda{{\unicode{923}}}
\def\upMu{{\unicode{924}}}
\def\upNu{{\unicode{925}}}
\def\upXi{{\unicode{926}}}
\def\upOmicron{{\unicode{927}}}
\def\upPi{{\unicode{928}}}
\def\upRho{{\unicode{929}}}
\def\upSigma{{\unicode{931}}}
\def\upTau{{\unicode{932}}}
\def\upUpsilon{{\unicode{933}}}
\def\upPhi{{\unicode{934}}}
\def\upChi{{\unicode{935}}}
\def\upPsi{{\unicode{936}}}
\def\upOmega{{\unicode{937}}}"""
    mathjax2_config = {
        'extensions': ['tex2jax.js'] , #+ [f'TeX/{pkg}.js' for pkg in latex_packages],
        'TeX': {'extensions': [f'{pkg}.js' for pkg in latex_packages],
                'Macros': {'upphi': r'\unicode{x03d5}'}},
        'HTML-CSS': {'fonts': ['STIX-Web']}
    }
    # mathjax2_config = {
    #     'loader': {'load': [f'[tex]/{pkg}' for pkg in latex_packages]},
    #     'tex': {'packages': {'[+]': [latex_packages]}}
    # }
