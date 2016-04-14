"""
This module provides a generic driver for chemistry network
calculations. Each chemistry network must define three methods:
__init__ uses an input cloud to initialize the network, dxdt returns
the instantaneous rate of change of abundance for all elements in the
network, and applyAbundances writes the abundances in the network back
to the cloud. The chemEvol procedure uses these to update the
chemistry network for a cloud.
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
from despotic.despoticError import despoticError
from .abundanceDict import abundanceDict
from copy import deepcopy
import scipy.constants as physcons
kB = physcons.k*1e7

def chemEvol(cloud, tFin, tInit=0.0, nOut=100, dt=None,
             tOut=None, network=None, info=None,
             addEmitters=False, evolveTemp='fixed',
             isobaric=False, tempEqParam=None,
             dEdtParam=None):
    """
    Evolve the abundances of a cloud using the specified chemical
    network.

    Parameters
       cloud : class cloud
          cloud on which computation is to be performed
       tFin : float
          end time of integration, in sec
       tInit : float
          start time of integration, in sec
       nOut : int
          number of times at which to report the temperature; this
          is ignored if dt or tOut are set
       dt : float
          time interval between outputs, in sec; this is ignored if
          tOut is set
       tOut : array
          list of times at which to output the temperature, in s;
          must be sorted in increasing order
       network : chemical network class
          a valid chemical network class; this class must define the
          methods __init__, dxdt, and applyAbundances; if None, the
          existing chemical network for the cloud is used
       info : dict
          a dict of additional initialization information to be passed
          to the chemical network class when it is instantiated
       addEmitters : Boolean
          if True, emitters that are included in the chemical
          network but not in the cloud's existing emitter list will
          be added; if False, abundances of emitters already in the
          emitter list will be updated, but new emiters will not be
          added to the cloud
       evolveTemp : 'fixed' | 'gasEq' | 'fullEq' | 'evol'
          how to treat the temperature evolution during the chemical
          evolution; 'fixed' = treat tempeature as fixed; 'gasEq' = hold
          dust temperature fixed, set gas temperature to instantaneous
          equilibrium value; 'fullEq' = set gas and dust temperatures to
          instantaneous equilibrium values; 'evol' = evolve gas
          temperature in time along with the chemistry, assuming the
          dust is always in instantaneous equilibrium
       isobaric : Boolean
          if set to True, the gas is assumed to be isobaric during the
          evolution (constant pressure); otherwise it is assumed to be
          isochoric; note that (since chemistry networks at present are
          not allowed to change the mean molecular weight), this option
          has no effect if evolveTemp is 'fixed'
       tempEqParam : None | dict
          if this is not None, then it must be a dict of values that
          will be passed as keyword arguments to the cloud.setTempEq,
          cloud.setGasTempEq, or cloud.setDustTempEq routines; only used
          if evolveTemp is not 'fixed'
       dEdtParam : None | dict
          if this is not None, then it must be a dict of values that
          will be passed as keyword arguments to the cloud.dEdt
          routine; only used if evolveTemp is 'evol'

    Returns
       time : array
          array of output times, in sec
       abundances : class abundanceDict
          an abundanceDict giving the abundances as a function of time
       Tg : array
          gas temperature as a function of time; returned only if
          evolveTemp is not 'fixed'
       Td : array
          dust temperature as a function of time; returned only if
          evolveTemp is not 'fixed' or 'gasEq'

    Raises
       despoticError, if network is None and the cloud does not already
       have a defined chemical network associated with it
    """

    # Check if we have been passed a new chemical network. If so,
    # initialize it and associate it with the cloud, unless it is the
    # same type as the current network; if not, make sure the cloud
    # has a network associated with it before proceeding.
    if network is not None:
        if not hasattr(cloud, 'chemnetwork'):
            cloud.chemnetwork = network(cloud=cloud, info=info)
        elif not isinstance(cloud.chemnetwork, network):
            cloud.chemnetwork = network(cloud=cloud, info=info)
    elif not hasattr(cloud, 'chemnetwork'):
        raise despoticError(
            'if network is None, cloud must have' +
            ' an existing chemnetwork')

    # Set up output times
    if tOut==None:
        if dt==None:
            tOut = tInit + np.arange(nOut+1)*float(tFin-tInit)/nOut
        else:
            tOut = np.arange(tInit, (tFin-tInit)*(1+1e-10), dt)

    # Sanity check on output times: eliminate any output times
    # that are not between tInit and tFin, and make sure final time is
    # tFin
    tOut1 = tOut[tOut >= tInit]
    tOut1 = tOut1[tOut1 <= tFin]
    if tOut1[-1] < tFin:
        tOut1 = np.append(tOut1, tFin)

    # Set the isobar
    if isobaric:
        isobar = cloud.Tg * cloud.nH
    else:
        isobar = -1

    # See how we're handling the temperature
    if evolveTemp == 'fixed':

        # Simplest case: fixed temperature, so just evolve the
        # chemical network alone
        xOut = odeint(cloud.chemnetwork.dxdt, cloud.chemnetwork.x,
                      tOut1)

    elif evolveTemp == 'gasEq':

        # We're evolving the temperature as well as the chemical
        # abundances, but we're doing so assuming the temperature is
        # always in equilibrium; we therefore use our wrapper class,
        # defined below
        dxdtwrap = _dxdt_wrapper(cloud, isobar, gasOnly=True, 
                                 tempEqParam=tempEqParam)
        xOut = odeint(dxdtwrap.dxdt_Teq, cloud.chemnetwork.x, tOut1)

        # Go back and compute the equilibrium gas temperature at each
        # of the requested output times
        Tg = np.zeros(len(tOut1))
        for i in range(Tg.size):
            cloud.chemnetwork.x = xOut[i,:]
            cloud.chemnetwork.applyAbundances()
            if tempEqParam is not None:
                cloud.setGasTempEq(**tempEqParam)
            else:
                cloud.setGasTempEq()
            Tg[i] = cloud.Tg

    elif evolveTemp == 'fullEq':

        # Same as the gasEq case, but now we set but the gas and dust
        # temperature to equilibrium
        dxdtwrap = _dxdt_wrapper(cloud, isobar, gasOnly=False,
                                 tempEqParam=tempEqParam)
        xOut = odeint(dxdtwrap.dxdt_Teq, cloud.chemnetwork.x, tOut1)

        # Go back and compute the equilibrium gas and dust
        # temperatures at each of the requested output times
        Tg = np.zeros(len(tOut1))
        Td = np.zeros(len(tOut1))
        for i in range(Tg.size):
            cloud.chemnetwork.x = xOut[i,:]
            cloud.chemnetwork.applyAbundances()
            if tempEqParam is not None:
                cloud.setTempEq(**tempEqParam)
            else:
                cloud.setTempEq()
            Tg[i] = cloud.Tg
            Td[i] = cloud.Td

    elif evolveTemp == 'evol':

        # Evolve the chemistry, gas, and dust temperatures
        # simultaneously
        dxdtwrap = _dxdt_wrapper(cloud, isobar,
                                 tempEqParam=tempEqParam,
                                 dEdtParam=dEdtParam)
        xTInit = np.append(cloud.chemnetwork.x, cloud.Tg)
        xTOut = odeint(dxdtwrap.dxTdt, xTInit, tOut1)
        xOut = xTOut[:,:-1]
        Tg = xTOut[:,-1]

        # Go back and compute equilibrium dust temperatures at each
        # output time
        Td = np.zeros(Tg.size)
        for i in range(Tg.size):
            cloud.chemnetwork.x = xOut[i,:]
            cloud.chemnetwork.applyAbundances()
            if tempEqParam is not None:
                cloud.setDustTempEq(**tempEqParam)
            else:
                cloud.setDustTempEq()
            Td[i] = cloud.Td

    else:
        raise despoticError(
            'chemEvol: invalid option ' + str(evolveTemp) +
            'for evolveTemp')

    # Write final results to chemnetwork
    cloud.chemnetwork.x = xOut[-1,:]

    # Write final abundances back to cloud
    cloud.chemnetwork.applyAbundances(addEmitters=addEmitters)

    # If the final time was not one of the requested output times,
    # chop it off the data to be returned
    if np.sum(tFin == tOut):
        xOut = xOut[:-1,:]
        if evolveTemp != 'fixed':
            Tg = Tg[:-1]
        if evolveTemp == 'evol' or evolveTemp == 'fullEq':
            Td = Td[:-1]

    # Return output
    if evolveTemp == 'fixed':
        return tOut, abundanceDict(cloud.chemnetwork.specList,
                                   np.transpose(xOut)), 
    elif evolveTemp == 'gasEq':
        return tOut, abundanceDict(cloud.chemnetwork.specList,
                                   np.transpose(xOut)), Tg
    else:
        return tOut, abundanceDict(cloud.chemnetwork.specList,
                                   np.transpose(xOut)), Tg, Td


# This is a helper class that wraps the dxdt function of a chemical
# network so that we can set the temperature (and possibly the
# density) to equilibrium while using it
class _dxdt_wrapper(object):

    # Initialization: just store a pointer to the cloud and assorted
    # other information that tells us how to call the temperature
    # setting routines
    def __init__(self, cloud, isobar, gasOnly=True, tempEqParam=None,
                 dEdtParam=None):
        self.cloud = cloud
        self.isobar = isobar
        self.gasOnly = gasOnly
        self.tempEqParam = tempEqParam
        self.dEdtParam = deepcopy(dEdtParam)
        if self.dEdtParam is not None:
            self.dEdtParam['fixedLevPop'] = True
            self.dEdtParam['sumOnly'] = True
            self.dEdtParam['gasOnly'] = True

    # Function to evolve the chemical network while also setting the
    # temperature and, if necessary, the density, to instantaneous
    # equilibrium
    def dxdt_Teq(self, xin, time):

        # Write abundances to the cloud
        self.cloud.chemnetwork.applyAbundances()

        # Set the temperature to the equilibrium
        if self.gasOnly:
            if self.tempEqParam is None:
                self.cloud.setGasTempEq()
            else:
                self.cloud.setGasTempEq(**self.tempEqParam)
        else:
            if self.tempEqParam is None:
                self.cloud.setTempEq()
            else:
                self.cloud.setTempEq(**self.tempEqParam)

        # If we are isobaric, update the density
        if self.isobar > 0:
            self.cloud.nH = isobar / self.cloud.Tg

        # Call the dxdt routine, and return its value
        dxdt = self.cloud.chemnetwork.dxdt(xin, time)

        idx = np.argmin(xin/np.abs(dxdt))
        print ('time = {:e}, Tg = {:f}, xCO = {:e}, xlim = {:e}' +
               ', x/xdotlim = {:e}, speclim = {:s}').format(
                   time, self.cloud.Tg, xin[4], xin[idx],
                   xin[idx]/dxdt[idx],
                   self.cloud.chemnetwork.specList[idx])

        return dxdt

    # Function to evolve the chemical network along with the
    # temperature; the input vector xTin contains the chemical
    # abundances first, then gas temperature in the last slot
    def dxTdt(self, xTin, time):

        # Output holder
        dxTdt_out = np.zeros(xTin.shape)

        # Update the cloud temperature
        self.cloud.Tg = xTin[-1]

        # If we're isobaric, update density
        if self.isobar > 0:
            self.cloud.nH = self.isobar / self.cloud.Tg

        # Write abundances to the cloud
        self.cloud.chemnetwork.x = xTin[:-1]
        self.cloud.chemnetwork.applyAbundances()

        # Compute the new dust temperature
        if self.tempEqParam is None:
            self.cloud.setDustTempEq()
        else:
            self.cloud.setDustTempEq(**self.tempEqParam)

        # Compute cooling rates; note that the level populations do
        # not need to be recomputed here, because they were already
        # computed in solving for the dust temperature value, which
        # depends on the line heating rate
        if self.dEdtParam is None:
            rates = self.cloud.dEdt(fixedLevPop=True, gasOnly=True, 
                                    sumOnly=True)
        else:
            rates = self.cloud.dEdt(**self.dEdtParam)
        dEdtGas = rates['dEdtGas']

        # Compute specific heat c_v at current temperature
        self.cloud.comp.computeCv(self.cloud.Tg)

        # Get temperature derivative from dE/dt by dividing by c_v
        if self.isobar > 0:
            dxTdt_out[-1] = dEdtGas/((self.cloud.comp.cv+1.0)*kB)
        else:
            dxTdt_out[-1] = dEdtGas/(self.cloud.comp.cv*kB)

        # Get time derivatives of chemical abundances
        dxTdt_out[:-1] = self.cloud.chemnetwork.dxdt(xTin[:-1], time)

        idx = np.argmin(xTin/np.abs(dxTdt_out))
        if idx != len(xTin)-1:
            print ('time = {:e}, Tg = {:f}, xCO = {:e}, xlim = {:e}' +
                   ', x/xdotlim = {:e}, speclim = {:s}').format(
                       time, self.cloud.Tg, xTin[5], xTin[idx],
                       xTin[idx]/dxTdt_out[idx],
                       self.cloud.chemnetwork.specList[idx])
        else:
            print ('time = {:e}, Tg = {:f}, xCO = {:e}' +
                   ', T/Tdot = {:e}').format(
                       time, self.cloud.Tg, xTin[5], 
                       xTin[-1]/dxTdt_out[-1])

        # Return
        return dxTdt_out
