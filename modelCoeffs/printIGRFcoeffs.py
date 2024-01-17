# Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
# Created by Corey J. Cochrane and Marshall J. Styczinski
# Maintained by Marshall J. Styczinski
# Contact: corey.j.cochrane@jpl.nasa.gov
########################################

import numpy as np
import os

def printIGRFcoeffs(fNameIn='IGRF13.shc', fNameOut='coeffsEarthIGRF13', datDir='modelCoeffs'):
    """
    Convert IGRF-13 spherical harmonic coefficients described in IGRF13.shc to a PlanetMag-compatible
    format.

    Parameters
    ----------
    fNameIn : str, default='IGRF13.shc'
        File name for table of spherical harmonics
    fNameOut : str, default='coeffsEarthIGRF13'
        File name base for output .csv files (one for g, one for h)
    datDir : str, default='modelCoeffs'
        Relative directory path where fNameIn is found and to which fNameOut will be printed
    """

    fPathIn = os.path.join(datDir, fNameIn)
    fPathOut = os.path.join(datDir, fNameOut)

    gHeader = 'Contains \'g\' Schmidt semi-normalized spherical harmonic coefficients for IGRF13 internal field model of Earth in gauss (G), organized as g_mn, so that each row represents a single n and each column a single m. See https://doi.org/10.1186/s40623-020-01288-x'
    gCols = '       g0n,       g1n,       g2n,       g3n,       g4n,       g5n,       g6n,       g7n,       g8n,       g9n,      g10n,      g11n,      g12n,      g13n'
    hHeader = 'Contains \'g\' Schmidt semi-normalized spherical harmonic coefficients for IGRF13 internal field model of Earth in gauss (G), organized as g_mn, so that each row represents a single n and each column a single m. See https://doi.org/10.1186/s40623-020-01288-x'
    hCols = '       h0n,       h1n,       h2n,       h3n,       h4n,       h5n,       h6n,       h7n,       h8n,       h9n,      h10n,      h11n,      h12n,      h13n'

    data = np.loadtxt(fPathIn, skiprows=5)
    g = np.zeros((13,14))
    h = np.zeros((13,14))
    nlist = np.asarray(data[:,0], dtype=np.int_)
    mlist = np.asarray(data[:,1], dtype=np.int_)
    gm = mlist[mlist>=0]
    hm = abs(mlist[mlist<0])
    gn = nlist[mlist>=0]
    hn = nlist[mlist<0]
    cIGRF2020 = 27
    glin = data[mlist>=0, cIGRF2020-1] / 1e5
    hlin = data[mlist<0, cIGRF2020-1] / 1e5
    for n,m,thisg in zip(gn,gm,glin):
        g[n-1,m] = thisg
    for n,m,thish in zip(hn,hm,hlin):
        h[n-1,m] = thish

    with open(f'{fPathOut}g.csv', 'w') as f:
        f.write(gHeader + '\n')
        f.write(gCols + '\n')
        for i, row in enumerate(g):
            rowStr = ','.join([f'{val:>10.6f}' if val != 0 else '         0' for val in row])
            f.write(rowStr + '\n')

    with open(f'{fPathOut}h.csv', 'w') as f:
        f.write(hHeader + '\n')
        f.write(hCols + '\n')
        for i, row in enumerate(h):
            rowStr = ','.join([f'{val:>10.6f}' if val != 0 else '         0' for val in row])
            f.write(rowStr + '\n')

if __name__ == '__main__':
    printIGRFcoeffs()
