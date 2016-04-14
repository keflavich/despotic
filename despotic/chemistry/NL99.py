"""
This module implements the reduced carbon-oxygen chemistry network of
Nelson & Langer (1999, ApJ, 524, 923).
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
import string
from despotic.despoticError import despoticError
from .shielding import fShield_CO_vDB
from .abundanceDict import abundanceDict
from .chemNetwork import chemNetwork
from .reactions import cr_reactions, photoreactions, reaction_matrix
import scipy.constants as physcons
import warnings

import os

# Check if we're trying to compile on readthedocs, in which case we
# need to disable a ton of this stuff because otherwise we get into
# problems with the lack of numpy and scipy
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

########################################################################
# Physical and numerical constants
########################################################################
if not on_rtd:
    kB = physcons.k*1e7
    mH = (physcons.m_p+physcons.m_e)*1e3
    _small = 1e-100

########################################################################
# List of species used in this chemistry network; note that the
# network does not track H or H+
########################################################################
if not on_rtd:
    specList = ['He+', 'H3+', 'OHx', 'CHx', 'CO', 'C', 'C+', 'HCO+', 'O',
                'M+']
    specListExtended = specList + ['H2', 'He', 'M', 'e-']


########################################################################
# Data on cosmic ray reactions
# Reactions are, in order:
# (1)       cr + H2      -> H3+  + e-
# (2)       cr + He      -> He+  + e-
########################################################################
if not on_rtd:
    _cr_react = [
        { 'spec' : ['H2', 'H3+', 'e-'], 'stoich' : [-1, 1, 1], 'rate': 2.0 },
        { 'spec' : ['He', 'He+', 'e-'], 'stoich' : [-1, 1, 1], 'rate': 1.1 }
    ]
    _cr = cr_reactions(specListExtended, _cr_react)

########################################################################
# Data on photoreactions
# Reactions are, in order:
# (1) h nu + C -> C+ + e
# (2) h nu + CHx -> C
# (3) h nu + CO -> C + O
# (4) h nu + OHx -> O
# (5) h nu + M -> M+ + e
# (6) h nu + HCO+ -> CO
########################################################################
if not on_rtd:
    _ph_react = [
        { 'spec' : ['C', 'C+', 'e-'], 'stoich' : [-1, 1, 1],
          'rate' : 3.0e-10*1.7, 'av_fac' : 3.0 },
        { 'spec' : ['CHx', 'C'], 'stoich' : [-1, 1],
          'rate' : 1.0e-9*1.7, 'av_fac' : 1.5 },
        { 'spec' : ['CO', 'C', 'O'], 'stoich' : [-1, 1, 1],
          'rate' : 1.0e-10*1.7, 'av_fac' : 1.7, 
          'shield_fac' : fShield_CO_vDB },
        { 'spec' : ['OHx', 'O'], 'stoich' : [-1, 1],
          'rate' : 5.0e-10*1.7, 'av_fac' : 1.9 },
        { 'spec' : ['HCO+', 'CO'], 'stoich' : [-1, 1],
          'rate' : 1.5e-10*1.7, 'av_fac' : 2.5 }
    ]
    _ph = photoreactions(specListExtended, _ph_react)

########################################################################
# Data on two-body reactions
# Reactions are, in order:
# (1)  H3+  + C   -> CHx  + H2
# (2)  H3+  + O   -> OHx  + H2
# (3)  H3+  + CO  -> HCO+ + H2
# (4)  He+  + H2  -> He
# (5)  He+  + CO  -> C+   + O  + He
# (6)  C+   + H2  -> CHx
# (7)  C+   + OHx -> HCO+
# (8)  O    + CHx -> CO
# (9)  C    + OHx -> CO
# (10) He+  + e-  -> He
# (11) H3+  + e-  -> H2
# (12) C+   + e-  -> C
# (13) HCO+ + e-  -> CO
# (14) M+   + e-  -> M
# (15) H3+  + M   -> M+  + H2
########################################################################
if not on_rtd:
    _twobody_react = [
        { 'spec'   : [  'H3+',   'C',  'CHx', 'H2'       ],
          'stoich' : [    -1 ,   -1 ,     1 ,   1        ] },
        { 'spec'   : [  'H3+',   'O',  'OHx', 'H2'       ],
          'stoich' : [    -1 ,   -1 ,     1 ,   1        ] },
        { 'spec'   : [  'H3+',  'CO', 'HCO+', 'H2'       ],
          'stoich' : [    -1 ,   -1 ,     1 ,   1        ] },
        { 'spec'   : [  'He+',  'H2',   'He'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [  'He+',  'CO',   'C+',  'O', 'He' ],
          'stoich' : [    -1 ,   -1 ,     1 ,   1 ,   1  ] },
        { 'spec'   : [   'C+',  'H2',  'CHx'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [   'C+', 'OHx', 'HCO+'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [    'O', 'CHx',   'CO'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [    'C', 'OHx',   'CO'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [  'He+',  'e-',   'He'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [  'H3+',  'e-',   'H2'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [   'C+',  'e-',    'C'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [ 'HCO+',  'e-',   'CO'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [   'M+',  'e-',    'M'             ],
          'stoich' : [    -1 ,   -1 ,     1              ] },
        { 'spec'   : [  'H3+',   'M',   'M+', 'H2'       ],
          'stoich' : [    -1 ,   -1 ,     1 ,   1        ] } ]
    _twobody = reaction_matrix(specListExtended, _twobody_react)

    # Two-body reaciton rate coefficients
    _k2 = np.array([ 
        2.0e-9, 8.0e-10, 1.7e-9, 7.0e-15, 1.6e-9, 4.0e-16, 1.0e-9, 
        2.0e-10, 5.8e-12, 9.0e-11, 1.9e-6, 1.4e-10, 3.3e-5, 
        3.8e-10, 2.0e-9])
    _k2Texp = np.array([ 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.64, -0.54, 
        -0.61, -1.0, -0.65, 0.0])

def _twobody_ratecoef(T):
    return _k2 * T**_k2Texp


########################################################################
# Set some default abundances
########################################################################
_xHedefault = 0.1
_xCdefault = 2.0e-4
_xOdefault = 4.0e-4
_xMdefault = 2.0e-7
_xH2 = 0.5

########################################################################
# Define the NL99 class
########################################################################
class NL99(chemNetwork):
    """
    This class the implements the chemistry network of Nelson & Langer
    (1999, ApJ, 524, 923).

    Parameters
       cloud : class cloud
          a DESPOTIC cloud object from which initial data are to be
          taken
       info : dict
          a dict containing additional parameters

    Remarks
       The dict info may contain the following key - value pairs:

       'xC' : float
          the total C abundance per H nucleus; defaults to 2.0e-4
       'xO' : float
          the total H abundance per H nucleus; defaults to 4.0e-4
       'xM' : float
          the total refractory metal abundance per H
          nucleus; defaults to 2.0e-7
       'sigmaDustV' : float
          V band dust extinction cross
          section per H nucleus; if not set, the default behavior
          is to assume that sigmaDustV = 0.4 * cloud.dust.sigmaPE
       'AV' : float
          total visual extinction; ignored if sigmaDustV is set
       'noClump' : bool
          if True, the clumping factor is set to 1.0; defaults to False
    """

    ####################################################################
    # Method to initialize
    ####################################################################
    def __init__(self, cloud=None, info=None):
        """
        Parameters
           cloud : class cloud
              a DESPOTIC cloud object from which initial data are to be
              taken
           info : dict
              a dict containing additional parameters

        Remarks
           The dict info may contain the following key - value pairs:

           'xC' : float
              the total C abundance per H nucleus; defaults to 2.0e-4
           'xO' : float
              the total H abundance per H nucleus;
              defaults to 4.0e-4
           'xM' : float
              the total refractory metal abundance per H
              nucleus; defaults to 2.0e-7
           'sigmaDustV' : float
              V band dust extinction cross
              section per H nucleus; if not set, the default behavior
              is to assume that sigmaDustV = 0.4 * cloud.dust.sigmaPE
           'AV' : float
              total visual extinction; ignored if
              sigmaDustV is set
           'noClump' : bool
              if True, the clump factor is set to
              1.0; defaults to False
        """

        # List of species for this network; provide a pointer here so
        # that it can be accessed through the class
        self.specList = specList
        self.specListExtended = specListExtended

        # Store the input info dict
        self.info = info

        # Array to hold abundances
        self.x = np.zeros(10)

        # Total metal abundance
        if info is None:
            self.xM = _xMdefault
        else:
            if 'xM' in info:
                self.xM = info['xM']
            else:
                self.xM = _xMdefault

        # Extract information from the cloud if one is given
        if cloud is None:

            # No cloud given, so set some defaults
            self.cloud = None

            # Physical properties
            self._xHe = _xHedefault
            self._ionRate = 2.0e-17
            self._NH = _small
            self._temp = _small
            self._chi = 1.0
            self._nH = _small
            self._AV = 0.0
            if info is not None:
                if 'AV' in info:
                    self._AV = info['AV']

            # Set initial abundances
            if info is None:
                self.x[6] = _xCdefault
                self.x[8] = _xOdefault
            else:
                if 'xC' in info:
                    self.x[6] = info['xC']
                else:
                    self.x[6] = _xCdefault
                if 'xO' in info:
                    self.x[8] = info['xO']
                else:
                    self.x[8] = _xOdefault
            self.x[9] = self.xM

        else:

            # Cloud is given, so get information out of it
            self.cloud = cloud

            # Sanity check: make sure cloud is pure H2
            if cloud.comp.xH2 != 0.5:
                raise despoticError(
                    "NL99 network only valid " + 
                    "for pure H2 composition")

            # Sanity check: make sure cloud contains some He, since
            # network will not function properly at He abundance of 0
            if cloud.comp.xHe == 0.0:
                raise despoticError(
                    "NL99 network requires " + 
                    "non-zero He abundance")

            # Set abundances

            # Make a case-insensitive version of the emitter list for
            # convenience
            try:
                emList = dict(zip(map(string.lower, 
                                      cloud.emitters.keys()), 
                                  cloud.emitters.values()))
            except:
                # This somewhat more bulky construction is required in
                # python 3
                lowkeys = [k.lower() for k in cloud.emitters.keys()]
                lowvalues = list(cloud.emitters.values())
                emList = dict(zip(lowkeys, lowvalues))

            # OH and H2O
            if 'oh' in emList:
                self.x[2] += emList['oh'].abundance
            if 'ph2o' in emList:
                self.x[2] += emList['ph2o'].abundance
            if 'oh2o' in emList:
                self.x[2] += emList['oh2o'].abundance
            if 'p-h2o' in emList:
                self.x[2] += emList['p-h2o'].abundance
            if 'o-h2o' in emList:
                self.x[2] += emList['o-h2o'].abundance

            # CO
            if 'co' in emList:
                self.x[4] = emList['co'].abundance

            # Neutral carbon
            if 'c' in emList:
                self.x[5] = emList['c'].abundance

            # Ionized carbon
            if 'c+' in emList:
                self.x[6] = emList['c+'].abundance

            # HCO+
            if 'hco+' in emList:
                self.x[7] = emList['hco+'].abundance

            # Sum input abundances of C, C+, CO, HCO+ to ensure that
            # all carbon is accounted for. If there is too little,
            # assume the excess is C+. If there is too much, throw an
            # error.
            if info is None:
                xC = _xCdefault
            elif 'xC' in info:
                xC = info['xC']
            else:
                xC = _xCdefault
            xCtot = self.x[4] + self.x[5] + self.x[6] + self.x[7]
            if xCtot < xC:
                # Print warning if we're altering existing C+
                # abundance.
                if 'c' in emList:
                    print("Warning: input C abundance is " + 
                        str(xC) + ", but total input C, C+, CHx, CO, " + 
                        "HCO+ abundance is " + str(xCtot) + 
                        "; increasing xC+ to " + str(self.x[6]+xC-xCtot))
                self.x[6] += xC - xCtot
            elif xCtot > xC:
                # Throw an error if input C abundance is smaller than
                # what is accounted for in initial conditions
                raise despoticError(
                    "input C abundance is " + 
                    str(xC) + ", but total input C, C+, CHx, CO, " + 
                    "HCO+ abundance is " + str(xCtot))

            # O
            if 'o' in emList:
                self.x[8] = emList['o'].abundance
            elif info is None:
                self.x[8] = _xOdefault - self.x[2] - self.x[4] - \
                    self.x[7]
            elif 'xO' in info:
                self.x[8] = info['xO'] - self.x[2] - self.x[4] - \
                    self.x[7]
            else:
                self.x[8] = _xOdefault - self.x[2] - self.x[4] - \
                    self.x[7]

            # As with C, make sure all O is accounted for, and if not
            # park the extra in OI
            if info is None:
                xO = _xOdefault
            elif 'xC' in info:
                xO = info['xO']
            else:
                xO = _xOdefault
            xOtot = self.x[2] + self.x[4] + self.x[7] + self.x[8]
            if xOtot < xO:
                # Print warning if we're altering existing O
                # abundance.
                if 'o' in emList:
                    print("Warning: input O abundance is " + 
                          str(xO) + ", but total input O, OHx, CO, " + 
                          "HCO+ abundance is " + str(xOtot) + 
                          "; increasing xO to " + str(self.x[8]+xO-xOtot))
                self.x[8] += xO - xOtot
            elif xOtot > xO:
                # Throw an error if input O abundance is smaller than
                # what is accounted for in initial conditions
                raise despoticError(
                    "input C abundance is " + 
                    str(xO) + ", but total input O, OHx, CO, " + 
                    "HCO+ abundance is " + str(xOtot))


        # Initial electrons = metals + C+ + HCO+
        xeinit = self.xM + self.x[6] + self.x[7]

        # Initial He+
        self.x[0] = self.xHe*self.ionRate / \
            (self.nH*(_k2[9]*self.temp**_k2Texp[9]*xeinit+_k2[3]*_xH2))

        # Initial H3+
        self.x[1] = _xH2*self.ionRate / \
            (self.nH*(_k2[10]*self.temp**_k2Texp[10]*xeinit+_k2[2]*self.x[8]))

        # Initial M+
        self.x[9] = self.xM

    ####################################################################
    # Define some properties so that, if we have a cloud, quantities
    # that are stored in the cloud point back to it
    ####################################################################
    @property
    def nH(self):
        """
        volume density of H nuclei
        """
        if self.cloud is None:
            return self._nH
        else:
            return self.cloud.nH

    @nH.setter
    def nH(self, value):
        if self.cloud is None:
            self._nH = value
        else:
            self.cloud.nH = value

    @property    
    def temp(self):
        """
        gas kinetic temperature
        """
        if self.cloud is None:
            return self._temp
        else:
            return self.cloud.Tg

    @temp.setter
    def temp(self, value):
        if self.cloud is None:
            self._temp = value
        else:
            self.cloud.Tg = value

    @property
    def sigmaNT(self):
        """
        non-thermal velocity dispersion
        """
        if self.cloud is None:
            return 0.0
        else:
            return self.cloud.sigmaNT

    @property
    def cfac(self):
        """
        clumping factor; cannot be set directly, calculated from temp
        and sigmaNT
        """
        if self.cloud is None:
            return 1.0
        else:
            if self.info is None:
                cs2 = kB * self.cloud.Tg / (self.cloud.comp.mu * mH)
                return np.sqrt(1.0 + 0.75*self.sigmaNT**2/cs2)
            elif 'noClump' in self.info:
                if self.info['noClump'] == True:
                    return 1.0
                else:
                    cs2 = kB * self.cloud.Tg / (self.cloud.comp.mu * mH)
                    return np.sqrt(1.0 + 0.75*self.cloud.sigmaNT**2/cs2)
            else:
                cs2 = kB * self.cloud.Tg / (self.cloud.comp.mu * mH)
                return np.sqrt(1.0 + 0.75*self.cloud.sigmaNT**2/cs2)

    @cfac.setter
    def cfac(self, value):
        raise despoticError(
            "cannot set cfac directly; set sigmaNT or temp instead")

    @property
    def xHe(self):
        """
        He abundance
        """
        if self.cloud is None:
            return self._xHe
        else:
            return self.cloud.comp.xHe

    @xHe.setter
    def xHe(self, value):
        if self.cloud is None:
            self._xHe = value
        else:
            self.cloud.comp.xHe = value

    @property
    def ionRate(self):
        """
        primary ionization rate from cosmic rays and x-rays
        """
        if self.cloud is None:
            return self._ionRate
        else:
            return self.cloud.rad.ionRate

    @ionRate.setter
    def ionRate(self, value):
        if self.cloud is None:
            self._ionRate = value
        else:
            self.cloud.rad.ionRate = value

    @property
    def chi(self):
        """
        ISRF strength, normalized to solar neighborhood value
        """
        if self.cloud is None:
            return self._chi
        else:
            return self.cloud.rad.chi

    @chi.setter
    def chi(self, value):
        if self.cloud is None:
            self._chi = value
        else:
            self.cloud.rad.chi = value

    @property
    def NH(self):
        """
        column density of H nuclei
        """
        if self.cloud is None:
            return self._NH
        else:
            return self.cloud.colDen / 2.0

    @NH.setter
    def NH(self, value):
        if self.cloud is None:
            self._NH = value
        else:
            self.cloud.colDen = 2.0*value

    @property
    def AV(self):
        """
        visual extinction in mag
        """
        if self.cloud is None:
            if self.info is None:
                return self._AV
            elif 'AV' in self.info:
                return self.info['AV']
            else:
                return self._AV
        else:
            if self.info is None:
                return 0.4 * self.cloud.dust.sigmaPE * self.NH
            elif 'sigmaDustV' in self.info:
                # Note factor to convert from mag to true
                # dimensionless units
                return self.NH * self.info['sigmaDustV'] / \
                    np.log(100**0.2)
            elif 'AV' in self.info:
                return self.info['AV']
            else:
                return 0.4 * self.cloud.dust.sigmaPE * self.NH

    @AV.setter
    def AV(self, value):
        if self.cloud is None:
            if self.info is None:
                self._AV = value
            elif 'AV' in self.info:
                self.info['AV'] = value
            else:
                self._AV = value
        else:
            if self.info is None:
                raise despoticError(
                    "cannot set AV directly unless it is part of info")
            elif 'AV' not in self.info:
                raise despoticError(
                    "cannot set AV directly unless it is part of info")
            else:
                self.info['AV'] = value
                
    ####################################################################
    # Override the abundances property of the base chemNetwork class
    # so that we return the derived abundances as well as the
    # variables ones. For the setter, let users set abundances, but if
    # they try to set ones that are derived, issue a warning.
    ####################################################################

    @property
    def abundances(self):
        """
        abundances of all species in the chemical network
        """
        self._abundances = abundanceDict(self.specListExtended,
                                         self.extendAbundances())
        return self._abundances

    @abundances.setter
    def abundances(self, value):
        # Loop over key-value pairs in the and update them, but skip
        # any that are in the extended species list but not the main
        # list; issue warning if this happens
        warn = False
        abd = abundanceDict(self.specList, self.x)
        for k, v in other.items():
            if (k in specListExtended) and \
               (k not in specList):
                warn = True
                continue
            abd[k] = v
        if warn:
            warnings.warn("For NL99 network, cannot set abundances"
            " of derived species H2, He, M, e-; abundances set only "
            " for other species")


    ####################################################################
    # Method to get derived abundances from ones being stored; this
    # adds slots for H2, HeI, MI, and e
    ####################################################################
    def extendAbundances(self, xin=None):
        """
        Compute abundances of derived species not directly followed in
        the network.

        Parameters
           xin : array
              abundances of species directly tracked in the network;
              if left as None, the abundances stored internally to the
              network are used

        Returns
           x : array
              abundances, including those of derived species
        """

        # Object we'll be returning
        xgrow = np.zeros(14)

        # Copy abundances if passed in; otherwise user stored ones
        if xin is None:
            xgrow[:10] = self.x
        else:
            xgrow[:10] = xin

        # H2 abundances is hardwired for NL99 network
        xgrow[10] = _xH2

        # He abundance = total He abundance - He+ abundance
        xgrow[11] = self.xHe - xgrow[0]

        # Neutral metal abundance = total metal abundance - ionized
        # metal abundance
        xgrow[12] = self.xM - xgrow[9]

        # e abundance = He+ + H3+ + C+ + HCO+ + M+
        xgrow[13] = xgrow[0] + xgrow[1] + xgrow[6] + xgrow[7] \
                    + xgrow[9]

        # Return
        return xgrow


    ####################################################################
    # Method to return the time derivative of all chemical rates
    ####################################################################
    def dxdt(self, xin, time):
        """
        This method returns the time derivative of all abundances for
        this chemical network.

        Parameters
           xin : array(10)
              current abundances of all species
           time : float
              current time; not actually used, but included as an
              argument for compatibility with odeint

        Returns
           dxdt : array(10)
              time derivative of x
        """

        # Get abundances of derived quantities
        xgrow = self.extendAbundances(xin)

        # Get clumping factor and effective density
        cfac = self.cfac
        n = self.nH * cfac

        # Use the cosmic ray and photoreaction rate calculators to get
        # their contributions to dxdt
        xdot = _cr.dxdt(xgrow, n, self.ionRate)
        import pdb; pdb.set_trace()
        xdot += _ph.dxdt(xgrow, n, self.chi, self.AV, 
                         [[self.NH*xgrow[4], self.NH/2.0]])


        # Add two-body reactions
        xdot += _twobody.dxdt(xgrow, n, _twobody_ratecoef(self.temp))

        # Return results
        return np.ravel(xdot)[:len(specList)]


    ####################################################################
    # Method to write the currently stored abundances to a cloud
    ####################################################################
    def applyAbundances(self, addEmitters=False):
        """
        This method writes the abundances produced by the chemical
        network to the cloud's emitter list.

        Parameters
           addEmitters : Boolean
              if True, emitters that are included in the chemical
              network but not in the cloud's existing emitter list will
              be added; if False, abundances of emitters already in the
              emitter list will be updated, but new emiters will not be
              added to the cloud

        Returns
           Nothing

        Remarks
           If there is no cloud associated with this chemical network,
           this routine does nothing and silently returns.
        """

        # SAFETY check: make sure we have an associated cloud to which
        # we can write
        if self.cloud == None:
            return

        # Make a case-insensitive version of the emitter list for
        # convenience
        try:
            emList = dict(zip(map(string.lower, 
                                  self.cloud.emitters.keys()), 
                              self.cloud.emitters.values()))
        except:
            # This somewhat more bulky construction is required in
            # python 3
            lowkeys = [k.lower() for k in self.cloud.emitters.keys()]
            lowvalues = list(self.cloud.emitters.values())
            emList = dict(zip(lowkeys, lowvalues))

        # Save rtios of ^12C to ^13C, and ^16O to ^18O
        if '13co' in emList and 'co' in emList:
            c13_12 = emList['13co'].abundance / \
                emList['co'].abundance
        if 'c18o' in emList and 'co' in emList:
            o18_16 = emList['c18o'].abundance / \
                emList['co'].abundance

        # OH, assuming OHx is half OH
        if 'oh' in emList:
            emList['oh'].abundance = self.x[2]/2.0
        elif addEmitters:
            try:
                self.cloud.addEmitter('oh', self.x[2]/2.0)
            except despoticError:
                print('Warning: unable to add OH; cannot find LAMDA file')

        # H2O, assuming OHx is half H2O, and that oH2O and pH2O are
        # in the same ratio as H2
        fp = self.cloud.comp.xpH2 / self.cloud.comp.xH2
        if 'ph2o' in emList:
            emList['ph2o'].abundance = self.x[2]/2.0*fp
        elif 'p-h2o' in emList:
            emList['p-h2o'].abundance = self.x[2]/2.0*fp
        elif addEmitters:
            try:
                self.cloud.addEmitter('ph2o', self.x[2]/2.0*fp)
            except despoticError:
                print('Warning: unable to add p-H2O; cannot find LAMDA file')
        if 'oh2o' in emList:
            emList['oh2o'].abundance = self.x[2]/2.0*(1-fp)
        elif 'o-h2o' in emList:
            emList['o-h2o'].abundance = self.x[2]/2.0*(1-fp)
        elif addEmitters:
            try:
                self.cloud.addEmitter('oh2o', self.x[2]/2.0*(1-fp))
            except despoticError:
                print('Warning: unable to add o-H2O; cannot find LAMDA file')

        # CO
        if 'co' in emList:
            emList['co'].abundance = self.x[4]
        elif addEmitters:
            try:
                self.cloud.addEmitter('co', self.x[4])
            except despoticError:
                print('Warning: unable to add CO; cannot find LAMDA file')

        # if we have 13CO or C18O, make their abundances match that of CO
        # multiplied by the appropriate isotopic abundances
        if '13co' in emList:
            emList['13co'].abundance = self.x[4]*c13_12
        if 'c18o' in emList:
            emList['c18o'].abundance = self.x[4]*o18_16

        # C
        if 'c' in emList:
            emList['c'].abundance = self.x[5]
        elif addEmitters:
            try:
                self.cloud.addEmitter('c', self.x[5])
            except despoticError:
                print('Warning: unable to add C; cannot find LAMDA file')

        # C+
        if 'c+' in emList:
            emList['c+'].abundance = self.x[6]
        elif addEmitters:
            try:
                self.cloud.addEmitter('c+', self.x[6])
            except despoticError:
                print('Warning: unable to add C+; cannot find LAMDA file')

        # HCO+
        if 'hco+' in emList:
            emList['hco+'].abundance = self.x[7]
        elif addEmitters:
            try:
                self.cloud.addEmitter('hco+', self.x[7])
            except despoticError:
                print('Warning: unable to add HCO+; cannot find LAMDA file')

        # O
        if 'o' in emList:
            emList['o'].abundance = self.x[8]
        elif addEmitters:
            try:
                self.cloud.addEmitter('o', self.x[8])
            except despoticError:
                print('Warning: unable to add O; cannot find LAMDA file')

