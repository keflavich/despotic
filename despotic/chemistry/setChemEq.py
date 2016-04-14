"""
This module contains methods that can be used to compute the
equilibrium chemical state of a chemical network, by integrating the
chemical state forward in time until it converges. There is no
guarantee, however, that equilibria are unique, or even that they
exist.
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
from .chemEvol import chemEvol
import scipy.constants as physcons
kB = physcons.k*1e7

# Small numerical value
__small = 1e-100

def setChemEq(cloud, tEqGuess=None, network=None, info=None,
              addEmitters=False, tol=1e-6, maxTime=1e16,
              verbose=False, smallabd=1e-15, convList=None,
              evolveTemp='fixed', isobaric=False, tempEqParam=None,
              dEdtParam=None, maxTempIter=50):
    """
    Set the chemical abundances for a cloud to their equilibrium
    values, computed using a specified chemical netowrk.

    Parameters
       cloud : class cloud
          cloud on which computation is to be performed
       tEqGuess : float
          a guess at the timescale over which equilibrium will be
          achieved; if left unspecified, the code will attempt to
          estimate this time scale on its own
       network : chemNetwork class
          a valid chemNetwork class; this class must define the
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
       evolveTemp : 'fixed' | 'iterate' | 'iterateDust' | 'gasEq' | 'fullEq' | 'evol'
          how to treat the temperature evolution during the chemical
          evolution:

          * 'fixed' = treat tempeature as fixed
          * 'iterate' = iterate between setting the gas temperature and
            chemistry to equilibrium
          * 'iterateDust' = iterate between setting the gas and dust
            temperatures and the chemistry to equilibrium
          * 'gasEq' = hold dust temperature fixed, set gas temperature to
            instantaneous equilibrium value as the chemistry evolves
          * 'fullEq' = set gas and dust temperatures to instantaneous
            equilibrium values while evolving the chemistry network
          * 'evol' = evolve gas temperature in time along with the
            chemistry, assuming the dust is always in instantaneous
            equilibrium

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
       tol : float
          tolerance requirement on the equilibrium solution
       convList : list
          list of species to include when calculating tolerances to
          decide if network is converged; species not listed are not
          considered. If this is None, then all species are considered
          in deciding if the calculation is converged.
       smallabd : float
          abundances below smallabd are not considered when checking for
          convergence; set to 0 or a negative value to consider all
          abundances, but beware that this may result in false
          non-convergence due to roundoff error in very small abundances
       maxTempIter : int
          maximum number of iterations when iterating between chemistry
          and temperature; only used if evolveTemp is 'iterate' or
          'iterateDust'
       verbose : Boolean
          if True, diagnostic information is printed as the calculation
          proceeds

    Returns
       converged : Boolean
          True if the calculation converged, False if not

    Raises
       despoticError, if network is None and the cloud does not already
       have a defined chemical network associated with it

    Remarks
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
        raise despoticError(
            'if network is None, cloud must have' +
            ' an existing chemnetwork')

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
            print("setChemEquil: estimating characteristic " + 
                  "equilibration timescale...")

        # Compute current time derivatives
        xdot1 = cloud.chemnetwork.dxdt(cloud.chemnetwork.x, 0.0)
        dt1 = np.amin(np.abs((cloud.chemnetwork.x+max(smallabd,0)) /
                            (xdot1+__small)))
        x2 = cloud.chemnetwork.x + dt1*xdot1
        xdot2 = cloud.chemnetwork.dxdt(x2, 0.0)
        dt2 = np.amin(np.abs((x2+max(smallabd,0)) / (xdot2+__small)))

        # Use larger of the two xdot's to define a timescale, but
        # divide by 10 for safety, with a minimum of 10^7 sec
        tEqGuess = max(dt1, dt2)
        tEqGuess /= 10.0
        tEqGuess = max(tEqGuess, 1e7)

        # If we're evolving the gas temperature too, estimate a
        # timescale for its evolution
        if evolveTemp == 'evol':
            if dEdtParam is None:
                rates = cloud.dEdt(gasOnly=True, sumOnly=True)
            else:
                dEdtParam1 = deepcopy(dEdtParam)
                dEdtParam1['gasOnly'] = True
                dEdtParam1['sumOnly'] = True
                rates = cloud.dEdt(**dEdtParam1)
            if isobaric:
                dTdt = rates['dEdtGas'] / \
                       ((cloud.comp.computeCv(cloud.Tg)+1)*kB)
            else:
                dTdt = rates['dEdtGas'] / \
                       (cloud.comp.computeCv(cloud.Tg)*kB)
            tEqGuess = max(tEqGuess, cloud.Tg/(np.abs(dTdt)+__small))

        # Make sure tEqGuess doesn't exceed maxTime
        if tEqGuess > maxTime:
            tEqGuess = maxTime

    # Print status
    if verbose:
        print("setChemEquil: estimated equilibration timescale = " + 
              str(tEqGuess) + " sec")

    # Decide which species we will consider in determining if
    # abundances are converged
    if convList == None:
        convList = cloud.chemnetwork.specList
    convArray = np.array([cloud.chemnetwork.specList.index(spec)
                          for spec in convList])

    # If we're isobaric, save the isobar
    if isobaric:
        isobar = cloud.Tg * cloud.nH

    # Outer loop, if we're iterating between temperature and chemistry
    if evolveTemp == 'iterate' or evolveTemp == 'iterateDust':
        tempConverge = False
    else:
        tempConverge = True
    itCount = 0
    while True:

        # Evolve the chemistry in time for estimated equilibrium
        # timescale and check convergence; if not converged, increase
        # time and keep running until we converge or maximum time is
        # reached.
        err = np.zeros(convArray.size)+10.0*tol
        t = 0.0
        tEvol = tEqGuess
        lastCycle = False
        while True:

            # Evolve for specified time
            if evolveTemp != 'iterate' and evolveTemp != 'iterateDust':
                out = chemEvol(cloud, t+tEvol, tInit=t, nOut=3,
                               evolveTemp=evolveTemp, isobaric=isobaric,
                               tempEqParam=tempEqParam,
                               dEdtParam=dEdtParam)
            else:
                out = chemEvol(cloud, t+tEvol, tInit=t, nOut=3,
                               evolveTemp='fixed', 
                               addEmitters=addEmitters)
            xOut = np.array(out[1].values())

            # Compute residual
            err = abs(xOut[convArray,-2]/(xOut[convArray,-1]+__small)-1)

            # If smallabd is set, exclude species will small abundances
            # from the calculation
            if smallabd > 0.0:
                err[xOut[convArray,-1] < smallabd] = 0.1*tol

            # Add temperature to residual if we're evolving it
            if evolveTemp == 'evol':
                TOut = xOut[2]
                err = np.append(err, abs(TOut[-2]/(TOut[-1]+__small)-1))

            # Print status
            if verbose:
                if evolveTemp != 'evol' or np.argmax(err) < len(err-1):
                    print(
                        "setChemEquil: evolved from t = " + str(t) + 
                        " to "+str(t+tEvol)+" sec, residual = " + 
                        str(np.amax(err)) + " for species " + 
                        cloud.chemnetwork.specList[convArray[np.argmax(err)]])
                else:
                    print(
                        "setChemEquil: evolved from t = " + str(t) + 
                        " to "+str(t+tEvol)+" sec, residual = " + 
                        str(np.amax(err)) + " for temperature")

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
                ad = abundanceDict(cloud.chemnetwork.specList,
                                   cloud.chemnetwork.x)
                print(
                    "setChemEquil: abundances converged: " + str(ad))
            else:
                print(
                    "setChemEquil: reached maximum time of " + 
                    str(maxTime) + " sec without converging")

        # Floor small negative abundances to avoid numerical problems
        idx = np.where(np.logical_and(cloud.chemnetwork.x <= 0.0,
                                      np.abs(cloud.chemnetwork.x) < smallabd))
        cloud.chemnetwork.x[idx] = smallabd

        # If we failed to converge on the chemistry, bail out now
        if np.amax(err) >= tol:
            return False

        # Are we iterating on temperature?
        if not tempConverge:

            # Yes, so update temperature
            Tglast = cloud.Tg
            Tdlast = cloud.Td
            if evolveTemp == 'iterate':
                if tempEqParam is None:
                    cloud.setGasTempEq()
                else:
                    cloud.setGasTempEq(**tempEqParam)
            else:
                if tempEqParam is None:
                    cloud.setTempEq()
                else:
                    cloud.setTempEq(**tempEqParam)

            # If we're isobaric, also update the density
            if isobaric:
                cloud.nH = isobar / cloud.Tg

            # Check for temperature convergence
            resid = max(abs((cloud.Tg-Tglast)/cloud.Tg),
                        abs((cloud.Td-Tdlast)/cloud.Td))
            if resid < tol:
                tempConverge = True
            else:
                tempConverge = False

            # Print status
            if verbose:
                print("setChemEquil: updated temperatures to " +
                      "Tg = {:f}, Td = {:f}, residual = {:e}"). \
                    format(cloud.Tg, cloud.Td, resid)
                if tempConverge:
                    print("Temperature converged!")

        # Break if we've also converged on the temperature
        if tempConverge:
            break

        # Update iteration counter and see if we have gone too many
        # times
        itCount += 1
        if itCount > maxTempIter:
            break

    # Write results to the cloud if we converged
    if tempConverge:
        cloud.chemnetwork.applyAbundances(addEmitters=addEmitters)

    # Report on whether we converged
    return tempConverge
