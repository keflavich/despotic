"""
This module contains some utility functions for momentum-driven winds
"""

import numpy as np
from scipy.special import erf
import scipy.constants as physcons
me = physcons.m_e * 1e3
mH = physcons.m_p * 1e3 + me
h = physcons.h * 1e7
c = physcons.c * 1e2
G = physcons.G * 1e3
kB = physcons.k * 1e7
e = physcons.e * c/1e1

def sxMach(mach):
    """
    Returns the dispersion in log column density versus Mach number

    Parameters:
       mach: float or arraylike
          Mach number

    Returns:
       sx: float or array
          dispersion of log column density
    """
    alpha = 2.5
    rfac = 0.5*(3.0-alpha)/(2.0-alpha) * \
           (1.0-np.asarray(mach)**(2*(2-alpha))) / \
           (1.0-np.asarray(mach)**(2*(3-alpha)))
    return np.sqrt(np.log(1.0+rfac*np.asarray(mach)**2/4.))

def zetaM(xcrit, sx):
    """
    Returns the mass fraction at column density x < xcrit

    Parameters:
       xcrit: float or arraylike
          maximum column density
       sx: float or arraylike
          dispersion of column densities

    Returns:
       zetaM: float or array
          mass fraction with x < xcrit
    """
    return 0.5*(1.0-erf((-2*np.asarray(xcrit)+np.asarray(sx)**2) /
                        (2.**1.5*np.asarray(sx))))

def zetaA(xcrit, sx):
    """
    Returns the area fraction at column density x < xcrit

    Parameters:
       xcrit: float or arraylike
          maximum column density
       sx: float or arraylike
          dispersion of column densities

    Returns:
       zetaA: float or array
          area fraction with x < xcrit
    """
    return 0.5*(1.0+erf((2*np.asarray(xcrit)+np.asarray(sx)**2) /
                        (2.**1.5*np.asarray(sx))))

def pM(x, sx):
    """
    Returns the mass PDF evaluated at column density x

    Parameters:
       x: float or arraylike
          dimensionless log column density
       sx: float or arraylike
          dispersion of column densities

    Returns:
       pM: float or array
          mass PDF of column densities evaluated at x
    """
    return 1.0/np.sqrt(2.0*np.pi*np.asarray(sx)**2) * \
        np.exp(-(np.asarray(x)-np.asarray(sx)**2/2.)**2 /
               (2.*np.asarray(sx)**2))

def pA(x, sx):
    """
    Returns the area PDF evaluated at column density x

    Parameters:
       x: float or arraylike
          dimensionless log column density
       sx: float or arraylike
          dispersion of column densities

    Returns:
       pA: float or array
          area PDF of column densities evaluated at x
    """
    return 1.0/np.sqrt(2.0*np.pi*np.asarray(sx)**2) * \
        np.exp(-(x+np.asarray(sx)**2/2.)**2/(2.*np.asarray(sx)**2))

def tX(abd, Omega, wl, muH=1.4):
    """
    Returns the line timescale tX

    Parameters:
       abd: float or array
          abundance relative to hydrogen
       muH: float or array
          gas mass per H nucleus, in units of H masses
       Omega: float or array
          transition oscillator strength
       wl: float or array
          transition wavelength in cm

    Returns:
       tX: float or array
          transition time parameter tX
    """
    mX = muH*mH/abd
    tX = Omega*wl*e**2 / (4*G*mX*me*c)
    return tX

