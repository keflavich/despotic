"""
This module defines the chemNetwork class. This is a purely abstract
class that defines the required elements of all chemistry
networks. Chemistry network classes should be derived from it.
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

from despotic.chemistry import abundanceDict
from despotic.despoticError import despoticError

class chemNetwork(object):
    """
    This is a purely abstract class that defines the methods that all
    chemistry networks are required to implement. Chemistry networks
    should be derived from it, and should override its
    methods. Attempting to instantiate this directly will lead to an
    error.

    Parameters
       cloud : class cloud
          a cloud object to which this network should be attached
       info : dict
          a dict of additional information to be passed to the network
          on instantiation

    Class attributes
       specList : list
          list of strings giving the names of the species being
          treated in the chemical network
       x : array
          array of abundances of the species in specList
       cloud : class cloud
          a cloud object to which this chemical network is attached;
          can be None
    """

    # specList is a list of strings, with the ith element giving the
    # name of the ith species in the network
    specList = None

    # x is a numpy array, with the ith element giving the abundace of
    # the ith species relative to hydrogen nuclei. x must be a rank 1
    # array with the same number of elements as specList.
    x = None

    # cloud is an object of class cloud to which this chemical network
    # is attached. Chemical networks need not be attached to clouds.
    cloud = None

    # The init method for a chemical network must accept two keywords:
    # cloud, which receives an object of class cloud, and info, which
    # received a dict containing optional additional arguments. This
    # routine should initialize specList and x.
    def __init__(self, cloud=None, info=None):
        raise despoticError(
            "chemNetwork is an abstract class, " + 
            "and should never be instantiated directly. Only " + 
            "instantiate classes derived from it.")

    # The dxdt method takes two arguments: xin is a numpy array in the
    # same format as x, giving the input abundances; time is a float
    # giving the time. It should return a numpy array with the same
    # number of elements as xin, giving the time derivative of the
    # abundances of the species in CGS units.
    def dxdt(self, xin, time):
        """
        This routine returns the time rate of change of the abundances
        for all species in the network.

        Parameters
           xin : array
              array of starting abundances
           time : float
              current time in sec

        Returns
           dxdt : array
              the time derivative of all species abundances
        """
        raise despoticError(
            "chemNetwork is an abstract class, " + 
            "and should never be instantiated directly. Only " + 
            "instantiate classes derived from it.")

    # The applyAbundances method is responsible for taking the
    # abundances stored in the array x and writing them to the
    # emitters in the cloud. The optional argument addEmitters, which
    # should default to False, specifies that emitters contained in
    # the network but not part of the emitter list for the cloud
    # should be added.
    def applyAbundances(self, addEmitters=False):
        """
        This method writes abundances from the chemical network back
        to the cloud to which this network is attached.

        Parameters
           addEmitters : bool
              if True, and the network contains emitters that are not
              part of the parent cloud, then the network will attempt
              to add them using cloud.addEmitter. Otherwise this
              routine will change the abundances of whatever emitters
              are already attached to the cloud, but will not add new
              ones.

        Returns:
           Nothing
        """

        raise despoticError(
            "chemNetwork is an abstract class, " + 
            "and should never be instantiated directly. Only " + 
            "instantiate classes derived from it.")

    # This defines abundances as a property of the chemNetwork, and
    # provides a convenient way of getting the abundances in a nice,
    # user-readable form using abundanceDict. Derived classes
    # generally will not need to override this definition, but they
    # are free to do so, as none of the chemical network driver
    # classes rely on these definitions.
    @property
    def abundances(self):
        """
        The current abundances of every species in the chemical
        network, stored as an abundanceDict.
        """
        self._abundances = abundanceDict(self.specList, self.x)
        return self._abundances

    @abundances.setter
    def abundances(self, other):
        self._abundances = abundanceDict(self.specList, self.x)
        self._abundances.update(other)
        self.applyAbundances()
