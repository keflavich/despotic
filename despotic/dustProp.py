"""
This module defines the dustProp class, which carries information
about the properties of the dust in a cloud.
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

class dustProp:
    """
    A class describing the properties of dust grains.

    Attributes
    ----------
    sigma10 : float
        dust opacity to thermal radiation at 10 K, in cm^2 H^-1
    sigmaPE : float
        dust opacity averaged over 8 - 13.6 eV, the range that
        dominates grain photoelectric heatin
    sigmaISRF : float
        dust opacity averaged the range that dominates grain starlight
        heating
    Zd : float
        dust abundance normalized to solar neighborhood value
    beta : float
        dust spectral index in the mm, sigma ~ nu^beta
    alphaGD : grain-gas coupling coefficient

    Methods
    ------
    None
    """

    # Initialize to generic MW values
    sigma10 = 2.0e-25
    sigmaPE = 1.0e-21
    sigmaISRF = 3.0e-22
    Zd = 1.0
    beta = 2.0
    alphaGd = 3.2e-34

