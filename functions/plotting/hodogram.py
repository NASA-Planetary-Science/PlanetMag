"""
    Functions for loading data and printing individual hodograms for each moon.
"""

# Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
# Created by Corey J. Cochrane and Marshall J. Styczinski
# Maintained by Marshall J. Styczinski
# Contact: corey.j.cochrane@jpl.nasa.gov
########################################

import numpy as np
import os
import logging
from hdf5storage import loadmat  # Install hdf5storage with pip

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Set up log messages
printFmt = '[%(levelname)s] %(message)s'
logHP = logging.getLogger('healpy')
stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter(printFmt))
logHP.setLevel(logging.WARNING)
log = logging.getLogger('PlanetMag')
log.setLevel(logging.DEBUG)
log.addHandler(stream)

mpl.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'text.latex.preamble': r'\usepackage{stix}'
})
meta = {
    'Creator': 'Created with PlanetMag'
}

mpEnd = 'noMP'
models = {
    'Jupiter': {
        0: f'JRM33C2020{mpEnd}',
        1: f'VIP4C1981{mpEnd}',
        2: f'O6K1997{mpEnd}',
        3: f'KS2005{mpEnd}',
        4: f'JRM09C2020{mpEnd}',
        5: f'JRM09C1981{mpEnd}',
        6: f'VIP4K1997{mpEnd}',
        7: f'JRM33C2020{mpEnd}'
    },
    'Saturn': {
        0: f'Cassini11{mpEnd}',
        1: f'B2010{mpEnd}',
        2: f'Cassini11{mpEnd}'
    },
    'Uranus': {
        0: f'AH5{mpEnd}',
        1: f'Q3{mpEnd}',
        2: f'AH5{mpEnd}'
    },
    'Neptune': {
        0: f'O8{mpEnd}',
        1: f'O8{mpEnd}'
    }
}


def loaddat(planet, moonName, model=None, datDir='out', figDir='figures'):
    """
    Load time series data from disk and format into horizontal and vertical for hodogram plotting.

    Parameters
    ----------
    planet : str
        Name of parent planet for hodogram data to load.
    moonName : str
        Name of moon for hodogram data to load.
    model : int, default=None
        Model index number as defined in GetModelOpts. Passing None selects the default.
    datDir : str, default='out'
        Data directory from which to load FFT data.
    figDir : str, default='figures'
        Output directory for figure files.

    Returns
    -------
    xx : float, shape 1xN
        Selected :math:`\\vec{B}` component for plotting on the :math:`x` axis in hodogram plot.
        The selected component depends on the moon, so that the most interesting components are
        plotted.
    yy : float, shape 1xN
        Selected :math:`\\vec{B}` component for plotting on the :math:`x` axis in hodogram plot.
        The selected component depends on the moon, so that the most interesting components are
        plotted.
    xlabel : str
        Text to display for label on x axis.
    ylabel : str
        Text to display for label on y axis.
    title : str
        Text to display at top for plot title.
    figFile : str
        Output file name for figure to print, including extension.
    """

    if model is None:
        model = 0

    fBase = f'evalB{moonName}'
    datFile = os.path.join(datDir, f'{fBase}{models[planet][model]}.mat')
    figFile = os.path.join(figDir, f'{moonName}Hodogram{models[planet][model]}.png')

    data = loadmat(datFile)
    Bvec_nT = data['BvecMoon']

    xx = Bvec_nT[0,:]
    xlabel = '$B_x$ (nT), IAU frame'
    if planet == 'Saturn':
        yy = Bvec_nT[2,:]
        ylabel = '$B_z$ (nT), IAU frame'
    else:
        yy = Bvec_nT[1,:]
        ylabel = '$B_y$ (nT), IAU frame'

    title = f'{moonName} hodogram, {models[planet][model][:-4]} model'

    return xx, yy, xlabel, ylabel, title, figFile


def plotHodogramSingle(xx, yy, xlabel, ylabel, title, figFile, deftFigsize=None):
    """
    Plots a single hodogram figure for given input data and labels.

    Parameters
    ----------
    xx : float, shape 1xN
        :math:`x` axis data for hodogram plot.
    yy : float, shape 1xN
        :math:`y` axis data for hodogram plot.
    xlabel : str
        Text to display for label on x axis.
    ylabel : str
        Text to display for label on y axis.
    title : str
        Text to display at top for plot title.
    figFile : str
        Output file name for figure to print, including extension.
    deftFigsize : tuple of float, shape 1x2, default=None
        Default figure size in screen inches (the default of matplotlib).
    """

    if deftFigsize is None:
        deftFigsize = (4, 2)

    fig = plt.figure(figsize=deftFigsize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])

    ax.plot(xx, yy, color='xkcd:forest green')
    ax.set_aspect(1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    plt.tight_layout()
    fig.savefig(figFile, bbox_inches='tight', dpi=300, metadata=meta)
    log.debug(f'Hodogram saved to file: {figFile}')
    plt.close()

    return


if __name__ == '__main__':
    xx, yy, xlabel, ylabel, title, figFile = loaddat('Jupiter', 'Europa')
    plotHodogramSingle(xx, yy, xlabel, ylabel, title, figFile)
