"""
This module implements various shielding functions for
photodissociation of species by the interstellar radiation field. Each
function defined here returns a shielding factor f_shield by which the
free-space photodissociation rate is reduced by line shielding. The
naming convention for the function is

fShield_X_Y

where X is the species name and Y is a string describing the source
paper / authors for the shielding function. For example, the CO
shielding function of van Dishoeck & Black (1988) is fShield_CO_vDB.
"""

########################################################################
# Copyright (C) 2013 Mark Krumholz
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################

import numpy as np
from scipy.interpolate import RectBivariateSpline

######################################################################
# Self-shielding of H2, following Draine & Bertoldi (1996, ApJ, 468,
# 269)
######################################################################

def fShield_H2_DB(NH2, sigma):
    """
    This function returns the shielding factor for H2
    photodissociation as a function of the H2 column density and
    velocity dispersion, based on the approximation formula of Draine
    & Bertoldi (1996)

    Parameters
    ----------
    NH2 : float or array
        H2 column density in cm^-2
    sigma : float or array
        velocity dispersion in cm s^-1

    Returns
    -------
    fShield : float or array
         the shielding factor for the input NH2 and sigma
    """

    x = NH2 / 5e14
    b5 = np.sqrt(2.0) * sigma / 1.0e5
    return 0.965/(1.0+x/b5)**2 + 0.035/(1.0+x)**0.5 * \
        np.exp(-8.5e-4*(1.0+x)**0.5)


######################################################################
# Shielding of CO by CO and H2, following van Dishoeck & Black (1988,
# ApJ, 334, 771)
######################################################################

# Tabulated data from van Dishoeck & Black
_logNCOvDB = np.array([0, 13, 14, 15, 16, 17, 18, 19]) # CO column densities
_logNH2vDB = np.array([0, 19, 20, 21, 22, 23])         # H2 column densities
_ThetavDB = np.array([ \
        [1.0, 9.681e-1, 7.764e-1, 3.631e-1, 7.013e-2, 1.295e-2, 1.738e-3, 9.985e-5], \
            [8.215e-1, 7.916e-1, 6.160e-1, 2.749e-1, 5.351e-2, 1.065e-2, 1.519e-3, 8.818e-5], \
            [7.160e-1, 6.900e-1, 5.360e-1, 2.359e-1, 4.416e-2, 8.769e-3, 1.254e-3, 7.558e-5], \
            [3.500e-1, 3.415e-1, 2.863e-1, 1.360e-1, 2.500e-2, 4.983e-3, 7.151e-4, 3.796e-5], \
            [4.973e-2, 4.877e-2, 4.296e-2, 2.110e-2, 4.958e-3, 9.245e-4, 1.745e-4, 8.377e-6], \
            [1.310e-4, 1.293e-4, 1.160e-4, 6.346e-5, 1.822e-5, 6.842e-6, 3.622e-6, 3.572e-7], \
            ])  # Tabulated shielding factors

# Extend van Dishoeck & Black's table 3 dex in each direction by
# linear extrapolation
_logNCOvDB1 = np.zeros(9)
_logNCOvDB1[:-1] = _logNCOvDB
_logNCOvDB1[-1] = 22
_exCO = (_logNCOvDB1[-1] - _logNCOvDB1[-2]) / (_logNCOvDB1[-2] - _logNCOvDB1[-3])
_logNH2vDB1 = np.zeros(7)
_logNH2vDB1[:-1] = _logNH2vDB
_logNH2vDB1[-1] = 25
_exH2 = (_logNH2vDB1[-1] - _logNH2vDB1[-2]) / (_logNH2vDB1[-2] - _logNH2vDB1[-3])
_ThetavDB1 = np.zeros((7,9))
_ThetavDB1[:6,:8] = _ThetavDB
for _i, _logNH2 in enumerate(_logNH2vDB):
    _ThetavDB1[_i,-1] = _ThetavDB1[_i,-2] * (_ThetavDB1[_i,-2]/_ThetavDB1[_i,-3])**_exCO
for _j, _logNH2 in enumerate(_logNCOvDB):
    _ThetavDB1[-1,_j] = _ThetavDB1[-2,_j] * (_ThetavDB1[-2,_j]/_ThetavDB1[-3,_j])**_exH2
_ThetavDB1[-1,-1] = np.sqrt(
    _ThetavDB1[-1,-2] * (_ThetavDB1[-1,-2]/_ThetavDB1[-1,-3])**_exCO *
    _ThetavDB1[-2,-1] * (_ThetavDB1[-2,-1]/_ThetavDB1[-3,-1])**_exH2)

# Create interpolation functions at various orders
_logThetavDBinterp1 = RectBivariateSpline(_logNH2vDB1, _logNCOvDB1, \
                                              np.log(_ThetavDB1), \
                                              kx=1, ky=1)
_logThetavDBinterp2 = RectBivariateSpline(_logNH2vDB1, _logNCOvDB1, \
                                              np.log(_ThetavDB1), \
                                              kx=2, ky=2)
_logThetavDBinterp3 = RectBivariateSpline(_logNH2vDB1, _logNCOvDB1, \
                                              np.log(_ThetavDB1), \
                                              kx=3, ky=3)

# van Dishoeck & Black's shielding function
def fShield_CO_vDB(NCO, NH2, order=1):
    """
    This function returns the shielding factor for CO
    photodissociation as a function of CO and H2 column densities,
    based on the model of van Dishoeck & Black (1987)

    Parameters
    ----------
    NCO : float or array
        CO column density in cm^-2
    NH2 : float or array
        H2 column density in cm^-2
    order : integer
        order of spline interpolation on van Dishoeck & Black's table;
        1 is the safest choice, but 2 and 3 are also provided, and may
        be more accurate in some ranges

    Returns
    -------
    fShield : array
        the shielding factor for the input NCO and NH2 values

    Raises
    ------
    ValueError if order is not 1, 2, or 3
    """
    logNCO = np.log10(np.clip(NCO, 10.**_logNCOvDB[0], 10.**_logNCOvDB[-1]))
    logNH2 = np.log10(np.clip(NH2, 10.**_logNH2vDB[0], 10.**_logNH2vDB[-1]))
    if order==1:
        return np.exp(_logThetavDBinterp1(logNH2, logNCO))
    elif order==2:
        return np.exp(_logThetavDBinterp2(logNH2, logNCO))
    elif order==3:
        return np.exp(_logThetavDBinterp3(logNH2, logNCO))
    else:
        raise ValueError, "order must be 1, 2, or 3"
