"""
This module contains the function lineProfLTE, which is capable of
computing the line profile for a cloud in LTE from a specified cloud
propeties and for a specified emitter.
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
from scipy.integrate import odeint
from emitterData import emitterData
from despoticError import despoticError

# Define some global physical constants in cgs units
import scipy.constants as physcons
kB = physcons.k/physcons.erg
c = physcons.c/physcons.centi
mH = physcons.m_p/physcons.gram
h = physcons.h*1e7
sigma=physcons.sigma/physcons.erg*physcons.centi**2
a = 4*sigma/c
G = physcons.G*1e3

# Small number to avoid divide by zeros
small = 1e-50

########################################################################
# This routine computes the line profile for a specified transition,
# assuming that the gas is in LTE
########################################################################
def lineProfLTE(emdat, u, l, R, denProf, TProf, \
                    vProf=0.0, sigmaProf=0.0, \
                    offset=0.0, TCMB=2.73, vOut=None, vLim=None, \
                    nOut=100, dv=None, mxstep=10000):
        """
        Return the brightness temperature versus velocity for a
        specified line, assuming the level populations are in LTE.

        Parameters
        ----------
        em : class emitterData
            emitterData object describing the emitting species for
            which the computation is to be made
        u : int
            upper state of line to be computed
        l : int
            lower state of line to be computed
        R : float
            cloud radius in cm
        denProf : float or callable
            If denProf is a float, this give the density in particles
            cm^-3 of the emitting species, which is taken to be
            uniform. denProf can also be a function giving the density
            as a function of radius; see remarks below for details.
        TProf : float or callable
            same as denProf, but giving the temperature in K
        vProf : float or callable, optional
            same as vProf, but giving the bulk radial
            velocity in cm/s; if omitted, bulk velocity is set to 0
        sigmaProf : float or callable, optional
            same as denProf, but giving the non-thermal
            velocity dispersion in cm/s; if omitted, non-thermal
            velocity dispersion is set to 0


        Returns
        -------
        TB : array
             brightness temperature as a function of velocity (in K)
        vOut : array
               velocities at which TB is computed (in cm s^-1)

        Additional parameters
        ---------------------
        offset : float, optional
            fractional distance from cloud center at which
            measurement is made; 0 = at cloud center, 1 = at
            cloud edge; valid values are 0 - 1
        vOut : sequence, optional
            sequence of velocities (relative to line center at 0) at
            which the output is to be returned
        vLim : sequence (2), optional
            maximum and minimum velocities relative to line center at
            which to compute TB
        nOut : int, optional
            number of velocities at which to output
        dv : float, optional
            velocity spacing at which to produce output
        TCMB : float, optional
            CMB temperature used as a background to the cloud, in
            K. Defaults to 2.73.
        mxstep : int, optional
            maximum number of steps in the ODE solver; default is
            10,000

	Raises
	------
	despoticError is the specified upper and lower state have no
	radiative transition between them, or if offset is not in the
	range 0 - 1

        Remarks
        -------
        The functions denProf, TProf, vProf, and sigmaProf, if
        specified, should accept one floating argument, and return one
        floating value. The argument r is the radial position within
        the cloud in normalized units, so that the center is at r = 0
        and the edge at r = 1. The return value should be the density,
        temperature, velocity, or non-thermal velocity dispesion at
        that position, in cgs units. 
        """

        # Step 1: safety check
        if emdat.EinsteinA[u,l] == 0.0:
            raise depoticError, 'no radiative transition from state ' \
                +str(u)+' to state '+str(l)+' found'
        if offset < 0.0 or offset > 1.0:
            raise despoticError, 'offset must be in the range 0 - 1'

        # Step 2: set up the helper class to compute normalization
        # constants
        te = _transferEqn(emdat, u, l, R, denProf, TProf, vProf, \
				  sigmaProf, offset)

        # Step 3: construct list of velocities at which to output
        if vOut == None:
            if dv == None:
                if vLim == None:
                    # No input given, so take velocity limits to be
                    # offset from line center by max of 5*sigmaTot + abs(v0)
                    vLim = [-5*te.sigmaTot - abs(te.v0), \
                                 5*te.sigmaTot + abs(te.v0)]
                # Compute vOut from vLim and nOut
                vOut = np.arange(vLim[0], vLim[1]*(1.0+1e-6), \
                                  (vLim[1]-vLim[0])/(nOut-1))
            else:
                # dv is non-zero, so set velocities from dv and nOut
                vOut = np.arange(-dv*(nOut/2.0), \
                                   dv*(nOut/2.0+1.0e-6), dv)


        # Step 4: integrate the ODE at the given velocities
        iOut = np.zeros(len(vOut))
        for i, v in enumerate(vOut):

            # Frequency normalized to line-center value
            f = 1 + v/c

            # Normalized CMB intensity at this wavelength
            ICMB = (2*h*emdat.freq[u,l]**3/c**2) / \
                (np.exp(h*f*emdat.freq[u,l]/(kB*TCMB))-1.0) / te.I0

            # Limits of integration
            sLim = [-np.sqrt(1.0-offset**2), np.sqrt(1.0-offset**2)]

            # Evaluate the integral
            iOut[i] = odeint(te.rhs, ICMB, sLim, \
                                 mxstep=mxstep, args = (f,))[1]

            # Subtract off the CMB
            iOut[i] -= ICMB

        # Step 5: convert intensity to brightness temperature; be
        # careful to handle 0 or negative intensities correctly
        TB = (h*emdat.freq[u,l]/kB) / \
            np.log(1.0 + 2.0*h*emdat.freq[u,l]**3 / \
			   (c**2*abs(iOut)*te.I0+small))
        TB[iOut < 0] *= -1.0
        TB[iOut == 0.0] = 0.0

        # Step 6: return
        return TB, vOut


########################################################################
# These are helper classes into which we load all the information we
# need about the various parameters in the transfer equation, so that
# we don't have to keep passing them around.
########################################################################
class _normFunc:
    def __init__(self, func, norm):
        self.func = func
        self.norm = norm
    def f(self, arg):
        return self.func(arg) / self.norm

class _normFunc2:
    def __init__(self, func, norm):
        self.func = func
        self.norm = norm
    def f(self, arg):
        return self.func(arg*self.norm)

# Dummy function that just returns 1.0
def _unity(r):
    return 1.0

class _transferEqn:

    # Function to return the RHS of the transfer equation
    def rhs(self, I, x, f):

        # Compute normalized radius
        r = np.sqrt(x**2 + self.offset**2)

        # Compute line shape function
        sigmaf = np.sqrt(self.betas**2*self.T.f(r) + \
                          self.betaNT**2*self.sigma.f(r)**2)
        f0 = 1 - self.beta*self.u.f(r)*np.sin(x/(r+small))
        phif = 1.0/np.sqrt(2*np.pi*sigmaf**2) * \
            np.exp(-(f-f0)**2/(2*sigmaf**2))

        # Return RHS
        return self.d.f(r) * (self.prefac/self.Z.f(self.T.f(r))) * \
            (np.exp(-self.Theta/self.T.f(r)) - self.tau0* \
                 (1.0-np.exp(-self.Theta/self.T.f(r)))*I) * phif


    def rhs1(self, x, I, f):

        # Compute normalized radius
        r = np.sqrt(x**2 + self.offset**2)

        # Compute line shape function
        sigmaf = np.sqrt(self.betas**2*self.T.f(r) + \
                          self.betaNT**2*self.sigma.f(r)**2)
        f0 = 1 - self.beta*self.u.f(r)*np.sin(x/(r+small))
        phif = 1.0/np.sqrt(2*np.pi*sigmaf**2) * \
            np.exp(-(f-f0)**2/(2*sigmaf**2))

        # Return RHS
        return self.d.f(r) * (self.prefac/self.Z.f(self.T.f(r))) * \
            (np.exp(-self.Theta/self.T.f(r)) - self.tau0* \
                 (1.0-np.exp(-self.Theta/self.T.f(r)))*I) * phif


    # Initialization function
    def __init__(self, emdat, u, l, R, denProf, TProf, vProf, \
                     sigmaProf, offset):

        # Get normalizations, and set up function pointers
        if isinstance(denProf, float):
            d0 = denProf
            self.d = _normFunc(_unity, 1.0)
        else:
            d0 = denProf(1.0)
            self.d = _normFunc(denProf, d0)
        if isinstance(TProf, float):
            T0 = TProf
            self.T = _normFunc(_unity, 1.0)
        else:
            T0 = TProf(1.0)
            self.T = _normFunc(TProf, T0)
        if isinstance(vProf, float):
            self.v0 = vProf
            self.u = _normFunc(_unity, 1.0)
        else:
            self.v0 = vProf(1.0)
            self.u = _normFunc(vProf, self.v0)
        if isinstance(sigmaProf, float):
            sigma0 = sigmaProf
            self.sigma = _normFunc(_unity, 1.0)
        else:
            sigma0 = sigmaProf(1.0)
            self.sigma = _normFunc(sigmaProf, sigma0)

        # Compute derived quantities that we need to store
        cs0 = np.sqrt(kB*T0/(emdat.molWgt*mH))
        self.beta = self.v0/c
        self.betas = cs0/c
        self.betaNT = sigma0 / c
        wavelength = c / emdat.freq[u, l]
        self.Theta = h*c / (wavelength*kB*T0)
        self.tau0 = emdat.EinsteinA[u,l] * wavelength**3 * d0 * R / (2*c)
        self.I0 = emdat.EinsteinA[u,l] * d0 * h * R
        self.prefac = d0 * emdat.levWgt[u] * np.exp(-emdat.levTemp[l]/T0) \
            / (4*np.pi)
        self.Z = _normFunc2(emdat.partFunc, T0)
        self.offset = offset
        self.sigmaTot = np.sqrt(cs0**2+self.v0**2)
