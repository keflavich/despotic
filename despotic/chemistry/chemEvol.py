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
