"""
This module defines the composition class, which carries information
about the bulk chemical composition of clouds and provides methods for
performing computations related to this data.
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
from despoticError import despoticError

# Physical constants describing H2
thetaRot = 85.3    # Rotation constant in K
thetaVib = 5984.   # Vibrational constant in K

class composition:
    """
    A class describing the chemical composition of an interstellar
    cloud, and providing methods to perform calculations using those
    properties.

    Attributes
    ----------
    xHI : float
        abundance of HI per H nucleus
    xoH2 : float
        abundance of ortho-H2 per H nucleus (note that the maximum
        possible value of xoH2 is 0.5, since it is per H nucleus)
    xpH2 : float
        abundance of paa-H2 per H nucleus  (note that the maximum
        possible value of xoH2 is 0.5, since it is per H nucleus)
    xHe : float
        abundance of He per H nucleus
    xe : float
        abundance of free electrons per H nucleus
    xHplus : float
        abundance of H+ per H nucleus
    mu : float
        mean mass per free particle, in units of H mass
    muH : float
        mean mass per H nucleus, in units of H mass
    qIon : float
        energy added to the gas per primary CR / x-ray ionization
    cv : float
        dimensionless specific heat per H nucleus at constant volume;
        the usual specific heat per unit volume may be obtained by
        multiplying this by nH * kB, and the specific heat per unit
        mass may be obtained by multiplying by nH * muH * kB

    Class methods
    -------------
    computeDerived -- compute derived quantities: mu, muH, qIon
    computeCv -- compute cv


    """

########################################################################
# Method to initialize
########################################################################
    def __init__(self):
        """
        This method initializes the class.

        Parameters
        ----------
        None

        Returns
        -------
        Nothing
        """

        # Initial values when class is created
        self.xHI = 0.0
        self.xoH2 = 0.0
        self.xpH2 = 0.0
        self.xHe = 0.0
        self.xe = 0.0
        self.xHplus = 0.0
        self.mu = 0.0           # Mean mass per particle
        self.muH = 0.0          # Mean mass per H nucleus
        self.qIon = 0.0
        self.cv = 0.0           # Specific heat per H nucleus


########################################################################
# Method to compute mu, muH, qIon from stored composition
########################################################################
    def computeDerived(self, nH):
        """
        Compute the derived quantities mu, muH, qIon

        Parameters
        ----------
        nH : float
            volume density in H cm^-3

        Returns
        -------
        Nothing

        Remarks
        -------
        For the purposes of this procedure, we treat electrons as
        massless.
        """

        # Mean particle masses
        self.mu = (self.xHI + self.xHplus + \
                       2.0*(self.xpH2 + self.xoH2) + \
                       4.0*self.xHe) / \
            (self.xHI + self.xHplus + self.xpH2 + \
                 self.xoH2 + self.xHe + self.xe)
        self.muH = self.xHI + self.xHplus + \
            2.0*(self.xpH2 + self.xoH2) + \
            4.0*self.xHe

        # Heating rate in eV for HI; fit from Draine (2011)
        qHI = self.xHI * (6.5 + 26.5*(self.xe/(self.xe+0.07))**0.5)
        # Heating rate in H2; fit to Glassgold, Galli, & Padovani
        # (2012), table 6
        lognH = np.log10(nH)
        if lognH < 2.0:
            qH2 = 10.0
        elif lognH < 4.0:
            qH2 = 10.0 + 3.0*(lognH-2.0)/2.0
        elif lognH < 7.0:
            qH2 = 13.0 + 4.0*(lognH-4.0)/3.0
        elif lognH < 10.0:
            qH2 = 17.0 + 1.0*(lognH-7.0)/3.0
        else:
            qH2 = 18.0
        qH2 *= 2.0 * (self.xpH2 + self.xoH2)

        # Total heating rate
        self.qIon = (qHI + qH2) * 1.6e-12


########################################################################
# Method to compute cv(T)
########################################################################
    def computeCv(self, T, noSet=False, Jmax=40):
        """
        Compute the specific heat per H nucleus, in erg K^-1 H^-1

        Parameters
        ----------
        T : float or array
            temperature in K
        noSet : Boolean
            if True, the value of cv stored in the class is not
            altered, but the calculated cv is still returned
        Jmax : int
            maximum J to be used in evaluating the rotational
            partition function; should be set to a value such that T
            << J(J+1) * thetaRot, there thetaRot = 85.3 K. Defaults to
            40.

        Returns
        -------
        cv : float or array
            value of cv
        """

        # Translational part
        cvtrans = 1.5 * (self.xHI + self.xHplus + self.xpH2 + \
                             self.xoH2 + self.xHe + self.xe)

        # Vibrational part
        x = -thetaVib/T
        cvvib = (self.xoH2 + self.xpH2) * x**2 * \
            np.exp(x) / (1.0-np.exp(x))**2

        # para-H2 rotational part. The value we want is given by
        # d/dT [ T^2/Z dZ/dT ] = 
        #     (1/Z) d/dT (T^2 dZ/dT) - [(T/Z) dZ/dT]^2

        # compute rotational partition function of pH2, zpH2. Also
        # compute d/dT (zpH2) and d/dT (T^2 dzpH2/dT).
        if self.xpH2 > 0:
            x = -thetaRot/T
            j = np.arange(0,Jmax+0.1,2)
            zpH2 = np.tensordot( 2*j+1, np.exp(np.outer(j*(j+1), x)), axes=1 )
            d_zpH2_dT = (thetaRot/T**2) * \
                np.tensordot( (2*j+1)*(j+1)*j, \
                               np.exp(np.outer(j*(j+1), x)), axes=1 )
            d_T2_dzpH2_dT_dT = x**2 * \
                np.tensordot( j**2*(j+1)**2*(2*j+1), \
                               np.exp(np.outer(j*(j+1), x)), axes=1 )
            cvpH2rot = self.xpH2 * \
                (d_T2_dzpH2_dT_dT / zpH2 - (T * d_zpH2_dT / zpH2)**2)
        else:
            cvpH2rot = 0.0

        # do ortho-H2; this is the same as para-H2, except that the
        # partition function is multiplied by an extra factor of 3 for
        # degeneracy, and that J is odd instead of even
        if self.xoH2 > 0:
            x = -thetaRot/T
            j = np.arange(1,Jmax+0.1,2)
            zoH2 = 3.0*np.tensordot( 2*j+1, np.exp(np.outer(j*(j+1), x)), axes=1 )
            d_zoH2_dT = 3.0*(thetaRot/T**2) * \
                np.tensordot( (2*j+1)*(j+1)*j, \
                               np.exp(np.outer(j*(j+1), x)), axes=1 )
            d_T2_dzoH2_dT_dT = 3.0*x**2 * \
                np.tensordot( j**2*(j+1)**2*(2*j+1), \
                               np.exp(np.outer(j*(j+1), x)), axes=1 )
            cvoH2rot = self.xoH2 * \
                (d_T2_dzoH2_dT_dT / zoH2 - (T * d_zoH2_dT / zoH2)**2)
        else:
            cvoH2rot = 0.0

        # Total cv; if T is an array, return an array; if not, return
        # a scalar
        if noSet:
            if hasattr(T, '__iter__'):
                return cvtrans + cvvib + cvpH2rot + cvoH2rot
            else:
                return (cvtrans + cvvib + cvpH2rot + cvoH2rot)[0]
        else:
            if hasattr(T, '__iter__'):
                self.cv = cvtrans + cvvib + cvpH2rot + cvoH2rot
            else:
                self.cv = (cvtrans + cvvib + cvpH2rot + cvoH2rot)[0]
            return self.cv
