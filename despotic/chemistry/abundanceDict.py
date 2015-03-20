"""
This module defines the abundanceDict class.
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
import collections
from despotic.despoticError import despoticError

class abundanceDict(collections.MutableMapping,dict):
    """
    An abundanceDict object is a wrapper around an array of
    abundances, and maps between human-readable chemical names
    (e.g. CO) and numeric indices in the array. This mapping is
    created when the dict is first initialized, and is immutable
    thereafter. However, the underlying array to which the mapping
    applies can be altered.
    """

########################################################################
# Initialization method
########################################################################
    def __init__(self, specList, x):
        """
        This method initializes the species list for this abundance
        dict, and also points to the numeric abundance array it wraps
        around.

        Parameters
        ----------
        specList : list of strings
             list of species names for this abundanceDict
        x : array of rank 1 or 2
             array of abundances; the length of the first dimension of
             x must be equal to the length of specList

        Returns
        -------
        Nothing
        """

        # Make sure input specList and x are properly formatted
        if np.rank(x) < 1:
            raise despoticError, "x must be " + \
                "a numpy array of rank >= 1"
        elif x.shape[0] != len(specList):
            raise despoticError, "first dimension of " + \
                "x must be same length as specList"

        self.x = x
        self.__specDict = {}
        for i, s in enumerate(specList):
            self.__specDict[s] = i


########################################################################
# getitem, setitem methods operate on the associated numpy array
########################################################################
    def __getitem__(self, key):
        """
        __getitem__ works just as for an ordinary dict
        """
        if key not in self.__specDict:
            raise KeyError
        return self.x[self.__specDict[key]]

    def __setitem__(self, key, value):
        """
        __setitem__ sets the value in the array x corresponding to the
        input species name.
        """
        if key not in self.__specDict:
            raise despoticError, "cannot add new species to " + \
                "abundanceDict"
        self.x[self.__specDict[key]] = value

########################################################################
# disallow deletions from the __specDict key
########################################################################
    def __delitem__(self, key):
        """
        raises an error, since abundanceDicts are
        immutable
        """
        raise despoticError, "cannot delete species from " + \
            "abundanceDict"

    def clear(self):
        """
        raises an error, since abundanceDicts are
        immutable
        """
        raise despoticError, "cannot delete species from " + \
            "abundanceDict"

    def pop(self, key):
        """
        raises an error, since abundanceDicts are
        immutable
        """
        raise despoticError, "cannot delete species from " + \
            "abundanceDict"

    def popitem(self):
        """
        raises an error, since abundanceDicts are
        immutable
        """
        raise despoticError, "cannot delete species from " + \
            "abundanceDict"

    def update(self):
        """
        raises an error, since abundanceDicts are
        immutable
        """
        raise despoticError, "cannot delete species from " + \
            "abundanceDict"

########################################################################
# define how to print abundanceDict objects
########################################################################
    def __repr__(self):
        """
        define how to print abundanceDict objects
        """
        stringRep = "{"
        for i, k in enumerate(self.__specDict):
            stringRep += "'" + k + "': " + \
                str(self.x[self.__specDict[k]])
            if i != len(self.__specDict)-1:
                stringRep += ", "
        stringRep += "}"
        return stringRep

########################################################################
# the methods below act just like they do on an ordinary dict
# whose keys are __specDict and whose values are the elements of x
########################################################################
    def __iter__(self):
        """
        __iter__ works just as for an ordinary dict
        """
        return dict.__iter__(self.__specDict)

    def __len__(self):
        """
        __len__ works just as for an ordinary dict
        """
        return self.x.size

    def __contains__(self, key):
        """
        __contains__ works just as for an ordinary dict
        """
        return dict.__contains__(self.__specDict, key)

    def keys(self):
        """
        keys works just as for an ordinary dict
        """
        return self.__specDict.keys()

    def values(self):
        """
        values returns a list of numpy arrays corresponding to the
        rows of x
        """
        return [self[k] for k in self.keys()]

    def has_key(self, key):
        """
        has_key works just as for an ordinary dict
        """
        return self.__specDict.has_key(key)

    def copy(self):
        """
        copy works just as for an ordinary dict
        """
        specList = ['']*len(self.__specDict)
        for k in self.__specDict.keys():
            specList[self.__specDict[k]] = k
        newAD = abundanceDict(specList, self.x)
        return newAD
