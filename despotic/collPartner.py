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

import numpy as np
from scipy.interpolate import interp1d
from .despoticError import despoticError
import pyximport
pyximport.install(setup_args={'include_dirs':np.get_include()})
from .collPartner_helper import colRates_all_scalar, \
    colRates_all_vector, colRates_some_scalar, \
    colRates_some_vector

class collPartner:
    """
    A utility class to store information about a particular collision
    partner for a given species, and to interpolate collision rates
    from those data

    Parameters
       fp : file
          a file object that points to the start of the collision
          rate data for one species in a LAMDA file
       nlev : int
          number of levels for the emitting species
       extrap : Boolean
          if True, then computing the collision rate with a
          temperature that is outside the table will result in the
          maximum or minimum value in the table being returned; if
          False, it will raise an error

    Class attributes
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
    """

    #####################################################################
    # Class initialization method
    #####################################################################
    def __init__(self, fp, nlev, extrap=True):
        """
        Initialize a collPartner object

        Parameters
           fp : file
              a file object that points to the start of the collision
              rate data for one species in a LAMDA file
           nlev : int
              number of levels for the emitting species
           extrap : Boolean
              if True, then computing the collision rate with a
              temperature that is outside the table will result in the
              maximum or minimum value in the table being returned; if
              False, it will raise an error
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
        self.tempTable = np.zeros(self.ntemp)
        line = ''
        while (line.strip() == ''):
            line = fp.readline()
            if line.strip()[0] == '!':
                line = ''
        for i, t in enumerate(line.split()):
            self.tempTable[i] = float(t)

        # Read table of collision rate coefficients
        self.colUpper = np.zeros(self.ntrans, dtype='int')
        self.colLower = np.zeros(self.ntrans, dtype='int')
        self.colRate = np.zeros((self.ntrans, self.ntemp))
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

        # Store whether extrapolation is allowed
        self.extrap = extrap

        # Compute log of temperatures and collision rates; these will
        # be used for interpolation
        self.logTempTable = np.log(self.tempTable)
        self.logColRate = np.log(self.colRate)


    #####################################################################
    # Method to return collision rates for every downward transition in
    # the table at a specified temperature; if keyword transition is
    # specified, the value is returned for a specified set of transitions
    #####################################################################
    def colRates(self, temp, transition=None):
        """
        Return interpolated collision rates for all transitions at a
        given temperature or list of temperatures

        Parameters
           temp : float | array
              temperature(s) at which collision rates are computed, in K
           transition : array of int, shape (2, N)
              list of upper and lower states for which collision rates
              are to be computed; default behavior is to computer for
              all known transitions

        Returns
           rates : array, shape (ntrans) | array, shape (ntrans, ntemp)
              collision rates at the specified temperature(s)
        """

        # If extrapolation is not allowed, make sure all temperatures
        # we've gotten are in the allowed range
        if not self.extrap:
            if np.any(np.array(temp) < self.tempTable[0]) or \
               np.any(np.array(temp) > self.tempTable[-1]):
                raise despoticError(
                    "temperature out of bounds in collPartner")

        # Handle the four possible cases; all transitions or some,
        # single scalar temperature or array of temperatures
        if transition is None:
            if hasattr(temp, '__iter__'):
                rates = colRates_all_vector(
                    self.logColRate, self.logTempTable, array(temp))
            else:
                rates = colRates_all_scalar(
                    self.logColRate, self.logTempTable, temp)
            return rates
        else:
            if hasattr(temp, '__iter__'):
                rates = colRates_some_vector(
                    self.colUpper, self.colLower, self.logColRate,
                    self.logTempTable, array(temp), transition)
            else:
                rates = colRates_some_scalar(
                    self.colUpper, self.colLower, self.logColRate,
                    self.logTempTable, temp, transition)
            return rates


    #####################################################################
    # Method to return the collision rate matrix
    #####################################################################
    def colRateMatrix(self, temp, levWgt, levTemp):
        """
        Return interpolated collision rates for all transitions at a
        given temperature, stored as an nlev x nlev matrix.

        Parameters
           temp : float
              temperature at which collision rates are computed, in K
           levWgt : array
              array of statistical weights for each level
           levTemp : array
              array of level energies, measured in K

        Returns
           k : array, shape (nlev, nlev)
              collision rates at the specified temperature; element i,j
              of the matrix gives the collision rate from state i to
              state j
        """
        k = np.zeros((self.nlev, self.nlev))
        # Downward transitions
        k[self.colUpper, self.colLower] += self.colRates(temp)
        # Upward transitions
        k[self.colLower, self.colUpper] \
            += k[self.colUpper, self.colLower] * \
            (levWgt[self.colUpper] / levWgt[self.colLower]) * \
            np.exp( -(levTemp[self.colUpper]-levTemp[self.colLower]) / \
                 temp )
        return k

########################################################################
# End of collPartner class
########################################################################

