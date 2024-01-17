"""
    Functions for plotting Fast Fourier Transform (FFT) data.
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


def loaddat(moonName, datDir='out', figDir='figures'):
    """
    Reloads FFT data from disk into numpy arrays.

    Parameters
    ----------
    moonName : str
        Name of moon to use for reloading data.
    datDir : str, default='out'
        Data directory from which to load FFT data.
    figDir : str, default='figures'
        Output directory for figure files.

    Returns
    -------
    T_h : float, shape 1xN
        Periods for excitation spectrum data in hours, i.e.\ the x axis values.
    BxFFT_nT : float, shape 1xN
        Magnitude of the :math:`\hat{x}` component of the excitation spectrum power for each period.
    ByFFT_nT : float, shape 1xN
        Magnitude of the :math:`\hat{y}` component of the excitation spectrum power for each period.
    BzFFT_nT : float, shape 1xN
        Magnitude of the :math:`\hat{z}` component of the excitation spectrum power for each period.
    xlabel : str
        Text to display for label on x axis.
    ylabel : str
        Text to display for label on y axis.
    title : str
        Text to display at top for plot title.
    figFile : str
        Output file name for figure to print, including extension.
    """

    fBase = f'{moonName}FTdata'
    datFile = os.path.join(datDir, f'{fBase}.mat')

    data = loadmat(datFile)
    magModelDescrip = data['magModelDescrip'][0]
    coordType = data['coordType'][0]
    BxFFT_nT = np.squeeze(data['B1vec1'])
    ByFFT_nT = np.squeeze(data['B1vec2'])
    BzFFT_nT = np.squeeze(data['B1vec3'])
    T_h = np.squeeze(data['T_h'])

    xlabel = 'Excitation period (h)'
    ylabel = 'Component amplitude (nT)'

    title = f'{moonName} excitation spectrum, {magModelDescrip}'
    figFile = os.path.join(figDir, f'{moonName}FFT{magModelDescrip}{coordType}.png')

    return T_h, BxFFT_nT, ByFFT_nT, BzFFT_nT, xlabel, ylabel, title, figFile


def plotFFTSingle(T_h, BxFFT_nT, ByFFT_nT, BzFFT_nT, xlabel, ylabel, title, figFile,
    deftFigsize=None, LW=None):
    """
    Plot Fast Fourier Transform (FFT) data for a single excitation spectrum.

    Parameters
    ----------
    T_h : float, shape 1xN
        Periods for excitation spectrum data in hours, i.e.\ the x axis values.
    BxFFT_nT : float, shape 1xN
        Magnitude of the :math:`\hat{x}` component of the excitation spectrum power for each period.
    ByFFT_nT : float, shape 1xN
        Magnitude of the :math:`\hat{y}` component of the excitation spectrum power for each period.
    BzFFT_nT : float, shape 1xN
        Magnitude of the :math:`\hat{z}` component of the excitation spectrum power for each period.
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
    LW : float, default=None
        Linewidth in pt.
    """

    if deftFigsize is None:
        deftFigsize = (5, 3.5)
    if LW is None:
        LW = 0.5

    fig = plt.figure(figsize=deftFigsize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])

    ax.plot(T_h, np.abs(BxFFT_nT), color='xkcd:green', label='$B_x$', linewidth=LW)
    ax.plot(T_h, np.abs(ByFFT_nT), color='xkcd:blue', label='$B_y$', linewidth=LW)
    ax.plot(T_h, np.abs(BzFFT_nT), color='xkcd:brick red', label='$B_z$', linewidth=LW)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([2, np.max(T_h)])
    ax.legend()

    plt.tight_layout()
    fig.savefig(figFile, bbox_inches='tight', dpi=300, metadata=meta)
    log.debug(f'FFT saved to file: {figFile}')
    plt.close()

    return


if __name__ == '__main__':
    T_h, BxFFT_nT, ByFFT_nT, BzFFT_nT, xlabel, ylabel, title, figFile = loaddat('Europa')
    plotFFTSingle(T_h, BxFFT_nT, ByFFT_nT, BzFFT_nT, xlabel, ylabel, title, figFile)
