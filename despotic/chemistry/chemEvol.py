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
from abundanceDict import abundanceDict

def chemEvol(cloud, tFin, tInit=0.0, nOut=100, dt=None,
                 tOut=None, network=None, info=None,
                 addEmitters=False):
    """
    Evolve the abundances of a cloud using the specified chemical
    network.

    Parameters
    ----------
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

    Returns
    -------
    abundances : class abundanceDict
        an abundanceDict giving the abundances as a function of time
    time : array of floats
        array of output times, in sec

    Raises
    ------
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
        raise despoticError, 'if network is None, cloud must have' + \
            ' an existing chemnetwork'

    # Set up output times
    if tOut==None:
        if dt==None:
            tOut = tInit + np.arange(nOut+1)*float(tFin-tInit)/nOut
        else:
            tOut = np.arange(tInit, (tFin-tInit)*(1+1e-10), dt)

    # Sanity check on output times: eliminate any output times
    # that are not between tInit and tFin
    tOut1 = tOut[tOut >= tInit]
    tOut1 = tOut1[tOut1 <= tFin]

    # Now evolve cloud for requested amount of time
    xOut = odeint(cloud.chemnetwork.dxdt, cloud.chemnetwork.x,
                  tOut1)

    # Write results to chemnetwork
    cloud.chemnetwork.x = xOut[-1,:]

    # If necessary, continue integrating up to tEvol
    if tOut1[-1] < tFin:
        xOut1 = odeint(cloud.chemnetwork.dxdt, cloud.chemnetwork.x,
                       np.array([tFin-tOut1[-1]]))
        cloud.chemnetwork.x = xOut1[-1,:]

    # Write abundances back to cloud
    cloud.chemnetwork.applyAbundances(addEmitters=addEmitters)

    # Return output
    return abundanceDict(cloud.chemnetwork.specList,
                         np.transpose(xOut)), tOut



def setChemEq(cloud, tEqGuess=None, network=None, info=None,
              addEmitters=False, tol=1e-6, maxTime=1e16,
              verbose=False, convList=None):
    """
    Set the chemical abundances for a cloud to their equilibrium
    values, computed using a specified chemical netowrk.

    Parameters
    ----------
    cloud : class cloud
        cloud on which computation is to be performed
    tEqGuess : float
        a guess at the timescale over which equilibrium will be
        achieved; if left unspecified, the code will attempt to
        estimate this time scale on its own
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
    tol : float
        tolerance requirement on the equilibrium solution
    convList : list
        list of species to include when calculating tolerances to
        decide if network is converged; species not listed are not
        considered. If this is None, then all species are considered
        in deciding if the calculation is converged.
    verbose : Boolean
        if True, diagnostic information is printed as the calculation
        proceeds

    Returns
    -------
    converged : Boolean
        True if the calculation converged, False if not

    Raises
    ------
    despoticError, if network is None and the cloud does not already
    have a defined chemical network associated with it

    Remarks
    -------
    The final abundances are written to the cloud whether or not the
    calculation converges.
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
        raise despoticError, 'if network is None, cloud must have' + \
            ' an existing chemnetwork'

    # Get initial timesale estimate if we were not given one. Picking
    # this is very tricky, because chemical networks often involve
    # reactions with a wide range of timescales, and the rates of
    # change starting from arbitrary initial conditions may be highly
    # non-representative of those found elsewhere in parameter
    # space. To make a guess, we compute dx/dt at the initial
    # conditions and at a point slightly perturbed from them, and then
    # use the most restrictive timestep we find.
    if tEqGuess == None:

        # Print status
        if verbose:
            print "setChemEquil: estimating characteristic " + \
                "equilibration timescale..."

        # Compute current time derivatives
        xdot1 = cloud.chemnetwork.dxdt(cloud.chemnetwork.x, 0.0)
        dt = cloud.chemnetwork.x / (xdot1+__small)
        x2 = cloud.chemnetwork.x + dt*xdot1
        xdot2 = cloud.chemnetwork.dxdt(x2, 0.0)

        # Use larger of the two xdot's to define a timescale
        tEqGuess = max(np.amax(cloud.chemnetwork.x/(xdot1+__small)),
                       np.amax(x2/(xdot2+__small)))

        # Make sure tEqGuess doesn't exceed maxTime
        if tEqGuess > maxTime:
            tEqGuess = maxTime

    # Print status
    if verbose:
        print "setChemEquil: estimated equilibration timescale = " + \
            str(tEqGuess) + " sec"

    # Decide which species we will consider in determining if things
    # are converged
    if convList == None:
        convList = cloud.chemnetwork.specList
    convArray = np.array([cloud.chemnetwork.specList.index(spec)
                          for spec in convList])

    # Now evolve in time for estimated equilibrium timescale and check
    # convergence; if not converged, increase time and keep running
    # until we converge or maximum time is reached.
    err = np.zeros(convArray.size)+10.0*tol
    t = 0.0
    tEvol = tEqGuess
    lastCycle = False
    while True:

        # Evolve for specified time
        xOut = odeint(cloud.chemnetwork.dxdt, cloud.chemnetwork.x,
                      np.array([t, t+tEvol/2.0, t+tEvol]))

        # Compute residual
        err = abs(xOut[-2,convArray]/(xOut[-1,convArray]+__small)-1)

        # Write results to chemnetwork
        cloud.chemnetwork.x = xOut[-1,:]

        # Print status
        if verbose:
            print "setChemEquil: evolved from t = " + str(t) + \
                " to "+str(t+tEvol)+" sec, residual = " + \
                str(np.amax(err)) + " for species " + \
                cloud.chemnetwork.specList[convArray[np.argmax(err)]]

        # Check for convergence
        if np.amax(err) < tol:
            break

        # Update time and timestep, or break if we've exceed maximum
        # allowed time
        if t + tEvol < maxTime:
            t += tEvol
            tEvol *= 2.0
        else:
            if lastCycle:
                break
            else:
                tEvol = maxTime-t
                lastCycle = True

    # Print status
    if verbose:
        if np.amax(err) < tol:
            print "setChemEquil: abundances converged: " + \
                str(abundanceDict(cloud.chemnetwork.specList,
                                  cloud.chemnetwork.x))
        else:
            print "setChemEquil: reached maximum time of " + \
                str(maxTime) + " sec without converging"

    # Write results to the cloud
    cloud.chemnetwork.applyAbundances(addEmitters=addEmitters)

    # Report on whether we converged
    return np.amax(err) < tol
