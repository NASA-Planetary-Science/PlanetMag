"""
Run this script after running RunPublication() with Matlab.
"""

# Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
# Created by Corey J. Cochrane and Marshall J. Styczinski
# Maintained by Marshall J. Styczinski
# Contact: corey.j.cochrane@jpl.nasa.gov
########################################

import numpy as np
import os
import logging
from healpy.projector import CartesianProj  # Install healpy with pip; must use wsl if installing
from healpy.pixelfunc import vec2pix, npix2nside  # healpy on Windows
from hdf5storage import loadmat  # Install hdf5storage with pip

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
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

vecCompChoice = 'mag'
ppgc = 540  # Number of points per great circle -- sets angular resolution
nLonMap = ppgc + 1
nLatMap = int(ppgc/2) + 1
latlonFontSize = 14
cLabelPad = 5
cLabelFontSize = 10
titleFontSize = 16
nLatTicks = 5
nLonTicks = 5
latMin = -90
latMax = 90
DO_360 = True  # Whether to plot from 0 to 360 (True) or from -180 to 180 (False).
MATCH_PREV = False  # Whether to use contour levels meant to match previous presentation slides
deftFigsize = (6,3)

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
prevCont = {
    'Jupiter': np.arange(3, 13.1, 1).astype(np.int_),
    'Saturn': np.arange(0.2, 0.56, 0.05),
    'Uranus': np.arange(0.1, 1.11, 0.1),
    'Neptune': np.arange(0.1, 0.81, 0.1)
}

latTicks = np.linspace(-90, 90, nLatTicks, endpoint=True, dtype=np.int_)
latMap_deg = np.linspace(-90, 90, nLatMap)
if DO_360:
    lonMap_deg = np.linspace(0, 360, nLonMap)
else:
    lonMap_deg = np.linspace(-180, 180, nLonMap)


def FixLons(lons):
    """
    Configure longitudes to be from 0 to 360.

    Parameters
    ----------
    lons : float, array-like
        Longitude(s) to adjust to be between 0 and 360.

    Returns
    -------
    lons : float, array-like
        Adjusted longitudes.
    """

    if np.all(lons == 360):
        return 360 * np.ones_like(lons)
    else:
        return lons[lons != 360] % 360


def LonHemisphere(longitude, EAST=True):
    """
    Retrieve the E/W hemisphere letter for degree unit labels.

    Parameters
    ----------
    longitude : float, array-like
        Longitude for which to retrieve hemisphere letter.
    EAST : bool, default=True
        Whether to treat longitude values as all east (used in maps running from 0--360).

    Returns
    -------
    hemisphere : str
        A string to accompany degree units in map labels to label the hemisphere. Either nothing
        (for 0), E, or W.
    """

    if EAST:
        longitude = FixLons(longitude)
    if longitude > 0:
        hemisphere = 'E'
    elif longitude < 0:
        hemisphere = 'W'
    elif longitude == 0:
        hemisphere = ''
    else:
        hemisphere = 'E'
    return hemisphere


def EastFormatted(longitude, EAST=True):
    """
    Formatting function to pass to axis labelers for longitudes.

    Parameters
    ----------
    longitude : float, array-like
        Longitude for which to generate tick labels.
    EAST : bool, default=True
        Whether to treat longitude values as all east (used in maps running from 0--360).

    Returns
    -------
    fmtString : str
        A string to represent the passed longitude for a given axis tick.
    """

    fmtString = u'{longitude:{num_format}}{degree}{hemisphere}'
    return fmtString.format(longitude=abs(longitude), num_format='g',
                            hemisphere=LonHemisphere(longitude, EAST=EAST),
                            degree=u'\u00B0')


def GetSign(val):
    """
    Get + or - for text labels of numbers that may be negative.

    Parameters
    ----------
    val : float
        Value for which to get + or - sign.

    Returns
    -------
    sign : str
        A string to prepend to text labels that may be negative.

    """
    if val < 0:
        sign = u'\u2013'
    else:
        sign = ''
    return sign


def Cformat(field):
    """
    Formatting function for magnetic field contour labels.

    Parameters
    ----------
    field : float, array-like
        Magnetic field component or magnitude value for contour to label.

    Returns
    -------
    fmtString : str
        A string to represent the contour label.
    """

    fmtString = u'{sign}{field:{num_format}}'
    return fmtString.format(field=abs(field), sign=GetSign(field), num_format='g')


def LatHemisphere(latitude):
    """
    Retrieve the N/S hemisphere letter for degree unit labels.

    Parameters
    ----------
    latitude : float, array-like
        Latitude for which to retrieve hemisphere letter.

    Returns
    -------
    hemisphere : str
        A string to accompany degree units in map labels to label the hemisphere. Either nothing
        (for 0), N, or S.
    """
    if latitude == 0:
        hemisphere = ''
    elif latitude > 0:
        hemisphere = 'N'
    else:
        hemisphere = 'S'
    return hemisphere


def LatFormatted(latitude):
    """
    Formatting function to pass to axis labelers for latitudes.

    Parameters
    ----------
    latitude : float, array-like
        Latitude for which to generate tick labels.

    Returns
    -------
    fmtString : str
        A string to represent the passed latitude for a given axis tick.
    """

    fmtString = u'{latitude:{num_format}}{degree}{hemisphere}'
    return fmtString.format(latitude=abs(latitude), num_format='g',
                            hemisphere=LatHemisphere(latitude),
                            degree=u'\u00B0')


# Generate lat/lon formatters
if DO_360:
    LonFormatter = tick.FuncFormatter(lambda v, pos: EastFormatted(v, EAST=True))
    lonMin = 0
    lonMax = 360
    lonTicks = np.linspace(lonMin, lonMax, nLonTicks, dtype=np.int_)
else:
    LonFormatter = tick.FuncFormatter(lambda v, pos: EastFormatted(v, EAST=False))
    lonMin = -180
    lonMax = 180
    lonTicks = np.linspace(lonMin, lonMax, nLonTicks, dtype=np.int_)

LatFormatter = tick.FuncFormatter(lambda v, pos: LatFormatted(v))
ConFormatter = tick.FuncFormatter(lambda v, pos: Cformat(v))


def SetMap(ax):
    """
    Format axis labels as a cylindrical projection map.

    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes
        Axis object to format as a map.
    """

    ax.set_xticks(lonTicks)
    ax.set_yticks(latTicks)
    ax.tick_params(axis='both', which='major', labelsize=latlonFontSize)
    ax.xaxis.set_major_formatter(LonFormatter)
    ax.yaxis.set_major_formatter(LatFormatter)

    return


def loadSurfaceMap(planet, model=None, ioDir='publication', fBase='surfMap_', vecComp='mag'):
    """
    Load data from disk for a magnetic field surface map.

    Parameters
    ----------
    planet : str
        Planet for which to load saved data.
    model : int or str, default=None
        Index number for model to load, as needed to complete the file name, or 'all'.
    ioDir : str, default='publication'
        Directory name relative to run dir where data files are found.
    fBase : str, default='surfMap_'
        Filename pattern base to use for loading data files.
    vecComp : str, default='mag'
        Vector component to use for map plotting. Options are 'x', 'y', 'z', 'r', 'theta', 'phi',
        and 'mag'.

    Returns
    -------
    Bout_G : list of float, array-like
        Magnetic field map of the selected vector component at evaluated HEALpix pixels.
    titles : list of str
        Title strings for placing at the top of each plot axis object.
    inpFile : list of str
        File paths sans extension used for each surface map.
    """

    if model is None:
        model = 0

    if type(model) == int:
        inpFile = [os.path.join(ioDir, f'{fBase}{planet}_{models[planet][model]}')]
    elif model == 'all':
        inpFile = [os.path.join(ioDir, f'{fBase}{planet}_{models[planet][thisModel]}')
                   for thisModel in range(1, np.size(list(models[planet].values())))]
    else:
        raise ValueError(f'model option "{model}" not supported. Pass None, "all", or an integer.')

    Bout_G = []
    titles = []
    for file in inpFile:
        data = loadmat(f'{file}.mat')
        magModelDescrip = data['magModelDescrip'][0, 0]
        Brtp_G = data['Brtp_G']
        Bxyz_G = data['Bxyz_G']
        Bmag_G = data['Bmag_G']

        vecCompStr = vecComp
        if vecComp == 'mag':
            Bcomp_G = Bmag_G
            vecCompStr = r'\mathrm{mag}'
        elif vecComp == 'x':
            Bcomp_G = Bxyz_G[0, :]
        elif vecComp == 'y':
            Bcomp_G = Bxyz_G[1, :]
        elif vecComp == 'z':
            Bcomp_G = Bxyz_G[2, :]
        elif vecComp == 'r':
            Bcomp_G = Brtp_G[0, :]
        elif vecComp in ['t', 'th', 'theta']:
            Bcomp_G = Brtp_G[1, :]
            vecCompStr = r'\theta'
        elif vecComp in ['p', 'phi']:
            Bcomp_G = Brtp_G[2, :]
            vecCompStr = r'\phi'
        else:
            raise ValueError(f'vecComp {vecComp} not recognized. Options are "mag", "r", '
                             '"theta", "phi", "x", "y", "z".')

        Bout_G.append(np.squeeze(Bcomp_G))
        titles.append(f'{planet} {magModelDescrip} $B_{vecCompStr}$ (G)')

    return Bout_G, titles, inpFile


def healpixMap(ax, Bplot_G, vecComp, title, levels=None, cmap=None):
    """

    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes
        Axis object for the given map.
    Bplot_G : float, array-like
        Magnetic field map of the selected vector component at evaluated HEALpix pixels.
    vecComp : str
        Vector component to use for map plotting. Options are 'x', 'y', 'z', 'r', 'theta', 'phi',
        and 'mag'.
    title : str
        Title string for placing at the top the plot.
    levels : float, array-like, default=None
        Contour levels to pass to ax.contour method.
    cmap : str, default=None
        Colormap identifier string to pass to ax.pcolormesh method.
    """

    nside = npix2nside(np.size(Bplot_G))
    cp = CartesianProj(flipconv='geo')
    cp.set_proj_plane_info(xsize=nLonMap, ysize=nLatMap,
                           lonra=np.array([lonMin, lonMax]), latra=None)
    vecFunc = lambda vec: vec2pix(nside, vec[0], vec[1], vec[2])
    pix = vecFunc(cp.xy2vec(cp.ij2xy()))
    plotData = Bplot_G[pix]

    dmax = dmin = None
    if cmap is None:
        if vecComp == 'mag':
            cmap = 'viridis'
        else:
            cmap = 'seismic'
            dmax = np.max(np.abs(Bplot_G))
            dmin = -dmax

    SetMap(ax)
    asymMap = ax.pcolormesh(lonMap_deg, latMap_deg, plotData, vmin=dmin, vmax=dmax,
                            shading='auto', cmap=cmap, rasterized=True)
    asymContours = ax.contour(lonMap_deg, latMap_deg, plotData,
                              levels=levels, colors='black')
    ax.clabel(asymContours, fmt=Cformat, fontsize=cLabelFontSize, inline_spacing=cLabelPad)
    ax.set_title(title, size=titleFontSize)
    ax.set_aspect(1)
    return


def plotAllIndividual():
    """
    Plot magnetic field surface maps from data saved to disk, separately.
    """

    for planet in ['Jupiter', 'Saturn', 'Uranus', 'Neptune']:
        Bout_G, titles, figFiles = loadSurfaceMap(planet, model='all', vecComp=vecCompChoice)
        if MATCH_PREV:
            doLevels = prevCont[planet]
        else:
            doLevels = None

        for Bplot_G, title, file in zip(Bout_G, titles, figFiles):
            outFile = f'{file}.pdf'
            fig = plt.figure(figsize=deftFigsize)
            grid = GridSpec(1, 1)
            ax = fig.add_subplot(grid[0, 0])

            healpixMap(ax, Bplot_G, vecCompChoice, title, levels=doLevels)
            plt.tight_layout()
            fig.savefig(outFile, format='pdf', dpi=300, metadata=meta)
            log.debug(f'Surface field map figure saved to file: {outFile}')
            plt.close()

    return


def plotDefaultTogether(fBase='surfMap_'):
    """
    Plot magnetic field surface maps for each planet, all together in one figure.

    Parameters
    ----------
    fBase : str, default='surfMap_'
        Filename pattern base to use for loading data files.
    """

    outFile = os.path.join('publication', f'{fBase}AllDefault_B{vecCompChoice}.pdf')
    fig = plt.figure(figsize=(2*deftFigsize[0], 2*deftFigsize[1]))
    grid = GridSpec(2, 2)
    planets = ['Jupiter', 'Saturn', 'Uranus', 'Neptune']
    axes = np.array([fig.add_subplot(grid[i, j]) for i in range(2) for j in range(2)])

    for ax, planet in zip(axes, planets):
        if MATCH_PREV:
            doLevels = prevCont[planet]
        else:
            doLevels = None

        Bplot_G, title, _ = loadSurfaceMap(planet, model=0, vecComp=vecCompChoice)
        healpixMap(ax, Bplot_G[0], vecCompChoice, title[0], cmap=None, levels=doLevels)

    plt.tight_layout()
    fig.savefig(outFile, bbox_inches='tight', format='pdf', dpi=300, metadata=meta)
    log.debug(f'Surface field map figure saved to file: {outFile}')
    plt.close()

    return


if __name__ == '__main__':
    plotDefaultTogether()
    plotAllIndividual()
