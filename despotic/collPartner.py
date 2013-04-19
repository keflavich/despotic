"""
This module defines the class collPartner, which stores information
about a single collision partner. This is a helper class, and is
always instantiated by giving it a file object that points to a
particular place in a Leiden database-formatted file, which is the
start of the listing for that collision partner. The class defines a
method that returns collision rates at a specified temperature.
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

from numpy import *
from scipy.interpolate import interp1d
from despoticError import *

class collPartner:
    """
    Store information about a particular collision partner for a given
    species.

    Attributes
    ----------
    nlev : int
        number of energy levels for the emitting species
    ntrans : int
        number of collisional transitions in the data table
    ntemp : int
        number of temperatures in the data table
    tempTable : array(ntemp)
        list of temperatues at which collision rate coefficients are
        given
    colUpper : int array(ntrans)
        list of upper states for collisions
    colLower : int array(ntrans)
        list of lower states for collisions
    colRate : array(ntrans, ntemp)
        table of downward collision rate coefficients, in cm^3 s^-1
    colRateInterp : list(ntrans) of functions
        each function in the list takes one variable, the temperature,
        as an argument, and returns the collision rate coefficient for
        the corresponding transition at the given temperature; only
        downard transitions are included

    Methods
    -------
    __init__ -- initialization method
    colRates -- method to return collision rates for all downward
        transitions at a given temperature or list of temperatures
    collRateMatrix -- method to return the collision rate matrix for
        all transitions for this collision patner as a function of
        temperature
    """

########################################################################
# Class initialization method
########################################################################
    def __init__(self, fp, nlev, extrap=False):
        """
        Initialize a collPartner object

        Parameters
        ----------
        fp : file
            a file object that points to the start of the collision
            rate data for one species in a LAMDA file
        extrap : Boolean
            if True, when the interpolation functions for this
            collision partner are created, data points will be added
            at low and high temperatures so that it is possible to
            extrapolate off the table. Extrapolation is done by
            assuming the the collision rate follows a powerlaw in
            temperature

        Returns
        -------
        Nothing
        """

        # Store number of levels
        self.nlev = nlev

        # Read number of transitions listed
        line = ''
        while (line.strip() == ''):
            line = fp.readline()
            if line.strip()[0] == '!':
                line = ''
        self.ntrans = int(line.strip())
        line = ''
        while (line.strip() == ''):
            line = fp.readline()
            if line.strip()[0] == '!':
                line = ''

        # Read temperature list
        self.ntemp = int(line.strip())
        self.tempTable = zeros(self.ntemp)
        line = ''
        while (line.strip() == ''):
            line = fp.readline()
            if line.strip()[0] == '!':
                line = ''
        for i, t in enumerate(line.split()):
            self.tempTable[i] = float(t)

        # Read table of collision rate coefficients
        self.colUpper = zeros(self.ntrans, dtype='int')
        self.colLower = zeros(self.ntrans, dtype='int')
        self.colRate = zeros((self.ntrans, self.ntemp))
        for i in range(self.ntrans):
            line = ''
            while (line.strip() == ''):
                line = fp.readline()
                if line.strip()[0] == '!':
                    line = ''
            linesplit = line.split()
            self.colUpper[i] = int(linesplit[1])-1  # Correct to 0 offset
            self.colLower[i] = int(linesplit[2])-1  # Correct to 0 offset
            for j in range(self.ntemp):
                self.colRate[i,j] = float(linesplit[j+3])

        # Mask zeros in collision rates to avoid numerical problems
        self.colRate[self.colRate==0.0] = 1.0e-30

        # Generate interpolating functions from table
        self._buildInterp(extrap)


########################################################################
# Method to return collision rates for every downward transition in the
# table at a specified temperature; if keyword transition is specied,
# the value is returned for a specified set of transitions
########################################################################
    def colRates(self, temp, transition=None):
        """
        Return interpolated collision rates for all transitions at a
        given temperature or list of temperatures

        Parameters
        ----------
        temp : float or array of floats
            temperature(s) at which collision rates are computed, in K
        transition : int array(2, N)
            list of upper and lower states for which collision rates
            are to be computed; default behavior is to computer for
            all known transitions

        Returns
        -------
        rates : array(ntrans) or array(ntrans, ntemp)
            collision rates at the specified temperature(s)
        """

        if transition==None:
            if hasattr(temp, '__iter__'):
                # Temperature is a list or array
                rates=zeros((self.ntrans, len(temp)))
                for i in range(self.ntrans):
                    rates[i,:]=self.colRateInterp[i](log(temp))
            else:
                # Temperature is a scalar
                rates=zeros(self.ntrans)
                for i in range(self.ntrans):
                    rates[i]=self.colRateInterp[i](log(temp))
            return exp(rates)
        else:
            u=transition[0]
            l=transition[1]
            if hasattr(temp, '__iter__'):
                # Temperature is a list or array
                rates=zeros((len(u), len(temp)))
                for i, ulev in enumerate(u):
                    idx=where(self.colUpper == u and self.colLower == l)[0]
                    if len(idx) == 1:
                        rates[i,:] = \
                            exp(self.colRateInterp[idx](log(temp)))
            else:
                # Temperature is a scalar
                rates=zeros(len(u))
                for i, ulev in enumerate(u):
                    idx=where(self.colUpper == u and self.colLower == l)[0]
                    if len(idx) == 1:
                        rates[i] = \
                            exp(self.colRateInterp[idx](log(temp)))
            return rates


########################################################################
# Method to return the collision rate matrix
########################################################################
    def colRateMatrix(self, temp, levWgt, levTemp):
        """
        Return interpolated collision rates for all transitions at a
        given temperature, stored as an nlev x nlev matrix.

        Parameters
        ----------
        temp : float
            temperature at which collision rates are computed, in K
        levWgt : array of float
            array of statistical weights for each level
        levTemp : array of float
            array of level energies, measured in K

        Returns
        -------
        k : array(nlev, nlev)
            collision rates at the specified temperature; element i,j
            of the matrix gives the collision rate from state i to
            state j
        """
        k = zeros((self.nlev, self.nlev))
        # Downward transitions
        k[self.colUpper, self.colLower] += \
            self.colRates(temp)
        # Upward transitions
        k[self.colLower, self.colUpper] += \
            k[self.colUpper, self.colLower] * \
            (levWgt[self.colUpper] / levWgt[self.colLower]) * \
            exp( -(levTemp[self.colUpper]-levTemp[self.colLower]) / \
                      temp )
        return k

########################################################################
# Method to turn extrapolation on/off
########################################################################
    def _buildInterp(self, extrap):
        """
        Build collision rate interpolation functions

        Parameters
        ----------
        extrap : Boolean
            true turns on extrapolation, false turns it off

        Returns
        -------
        Nothing
        """
        # Generate interpolating functions from table
        self.colRateInterp = []
        if self.ntemp > 1:
            if extrap==False:   # Standard case, no extrapolation allowed
                for i in range(self.ntrans):
                    self.colRateInterp.append( \
                        interp1d(log(self.tempTable), \
                                     log(self.colRate[i,:]), \
                                     kind='linear') )
            else:  # Extrapolation allowed
                # In this case we create a new collision rate table, which
                # contains linearly extrapolated data points at very high
                # an low temperatures, and use that to construct our
                # inteprolating functions
                logTempExtrap = zeros(self.ntemp+2)
                logTempExtrap[0] = -30.0
                logTempExtrap[1:-1] = log(self.tempTable)
                logTempExtrap[-1] = 30.0
                logColRateExtrap = zeros(self.ntemp+2)
                for i in range(self.ntrans):
                    logColRateExtrap[1:-1] = log(self.colRate[i,:])
                    logColRateExtrap[0] = logColRateExtrap[1] + \
                        (logTempExtrap[0] - logTempExtrap[1]) * \
                        (logColRateExtrap[1] - logColRateExtrap[2]) / \
                        (logTempExtrap[1] - logTempExtrap[2])
                    logColRateExtrap[-1] = logColRateExtrap[-2] + \
                        (logTempExtrap[-1] - logTempExtrap[-2]) * \
                        (logColRateExtrap[-2] - logColRateExtrap[-3]) / \
                        (logTempExtrap[-2] - logTempExtrap[-3])
                    self.colRateInterp.append( \
                        interp1d(logTempExtrap, logColRateExtrap, \
                                     kind='linear') )
        else:
            # Special case: the table only gives a single
            # temperature. In this case, just construct an
            # interpolating function that assumes the rate coefficient
            # is constant
            logTempExtrap = zeros(3)
            logTempExtrap[0] = -30.0
            logTempExtrap[1] = log(self.tempTable[0])
            logTempExtrap[2] = 30.0
            logColRateExtrap = zeros(3)
            for i in range(self.ntrans):
                logColRateExtrap[:] = log(self.colRate[i,0])
                self.colRateInterp.append( \
                    interp1d(logTempExtrap, logColRateExtrap, \
                                 kind='linear') )
        

########################################################################
# End of collPartner class
########################################################################

