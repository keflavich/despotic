"""
This module implements the reduced carbon-oxygen chemistry network of
Nelson & Langer (1999, ApJ, 524, 923), with the added H2 chemistry of
Glover & MacLow (2007, ApJS, 169, 239), as combined in Glover & Clark
(2012, MNRAS, 421, 9).
"""

########################################################################
# Copyright (C) 2015 Mark Krumholz
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
from .shielding import fShield_CO_vDB, fShield_H2_DB
from .abundanceDict import abundanceDict
from .chemNetwork import chemNetwork
from .reactions import cr_reactions, photoreactions, gr_reactions
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
# List of species used in this chemistry network; specList is the set
# of species that are actually evolved, while specListExtended
# includes both species that are evolved and those whose abundances
# are derived from conservation laws.
########################################################################
if not on_rtd:
    specList = ['H+', 'H2', 'H3+', 'He+', 'OHx', 'CHx', 'CO', 'C', 
                'C+', 'HCO+', 'O', 'M+']
    specListExtended = specList + ['H', 'He', 'M', 'e-']


########################################################################
# Data on cosmic ray reactions
# Reactions are, in order:
# (1)       cr + H       -> H+   + e-
# (2)       cr + He      -> He+  + e-
# (3)       cr + H2      -> H3+  + H
########################################################################
if not on_rtd:
    _cr_react = [
        { 'spec' : ['H',  'H+',  'e-'], 'stoich' : [-1, 1, 1], 'rate': 1.0 },
        { 'spec' : ['He', 'He+', 'e-'], 'stoich' : [-1, 1, 1], 'rate': 1.1 },
        { 'spec' : ['H2', 'H3+', 'H' ], 'stoich' : [-1, 1, 1], 'rate': 2.0 }
    ]
    _cr = cr_reactions(specListExtended, _cr_react)

########################################################################
# Data on photoreactions
# Reactions are, in order:
# (1) h nu + H2   -> H  + H
# (2) h nu + CO   -> C  + O
# (3) h nu + CI   -> C+ + e
# (4) h nu + CHx  -> C  + H
# (5) h nu + OHx  -> O  + H
# (6) h nu + M    -> M+ + e
# (7) h nu + HCO+ -> CO + H+
########################################################################
if not on_rtd:
    _ph_react = [
        { 'spec' : ['H2', 'H'], 'stoich' : [-1, 2], 
          'rate' : 3.3e-11*1.7, 'av_fac' : 3.74, 
          'shield_fac' : fShield_H2_DB },
        { 'spec' : ['CO', 'C', 'O'], 'stoich' : [-1, 1, 1], 
          'rate' : 1.2e-10*1.7, 'av_fac' : 3.53, 
          'shield_fac' : fShield_CO_vDB },
        { 'spec' : ['C', 'C+', 'e-'], 'stoich' : [-1, 1, 1], 
          'rate' : 1.8e-10*1.7, 'av_fac' : 3.0 },
        { 'spec' : ['CHx', 'C', 'H'], 'stoich' : [-1, 1, 1], 
          'rate' : 1.0e-9, 'av_fac' : 1.5 },
        { 'spec' : ['OHx', 'O', 'H'], 'stoich' : [-1, 1, 1],
          'rate' : 5.0e-10, 'av_fac' : 1.7 },
        { 'spec' : ['M', 'M+', 'e-'], 'stoich' : [-1, 1, 1],
          'rate' : 2.0e-10*1.7, 'av_fac' : 1.9 },
        { 'spec' : ['HCO+', 'CO', 'H+'], 'stoich' : [-1, 1, 1],
          'rate' : 1.5e-10, 'av_fac' : 2.5 } ]
    _ph = photoreactions(specListExtended, _ph_react)

########################################################################
# Data on two-body reactions
# Reactions are, in order:
# (1)  H3+  + C   -> CHx  + H2
# (2)  H3+  + O   -> OHx  + H2
# (3)  H3+  + CO  -> HCO+ + H2
# (4)  He+  + H2  -> He   + H  + H+
# (5)  He+  + CO  -> C+   + O  + He
# (6)  C+   + H2  -> CHx  + H
# (7)  C+   + OHx -> HCO+
# (8)  O    + CHx -> CO   + H
# (9)  C    + OHx -> CO   + H
# (10) He+  + e-  -> He
# (11) H3+  + e-  -> H2   + H
# (12) H3+  + e-  -> 3H
# (13) C+   + e-  -> C
# (14) HCO+ + e-  -> CO   + H
# (15) H+   + e-  -> H
# (16) H2   + H   -> 3H
# (17) H2   + H2  -> H2   + 2H
# (18) H    + e-  -> H+   + 2e-
# (19) He+  + H2  -> H+   + He  + H
# (20) H    + H   -> H2              (grain catalyzed)
# (21) H+   + e-  -> H               (grain catalyzed)
# (22) M+   + e-  -> M
# (23) H3+  + M   -> M+   + H   + H2
########################################################################
if not on_rtd:
    _twobody_react = [
        { 'spec'   : [ 'H3+' , 'C'  , 'CHx' , 'H2'      ], 
          'stoich' : [  -1   ,  -1  ,  1    ,  1        ] },
        { 'spec'   : [ 'H3+' , 'O'  , 'OHx' , 'H2'      ], 
          'stoich' : [  -1   ,  -1  ,  1    ,  1        ] },
        { 'spec'   : [ 'H3+' , 'CO' , 'HCO+', 'H2'      ], 
          'stoich' : [  -1   ,  -1  ,  1    ,  1        ] },
        { 'spec'   : [ 'He+' , 'H2' , 'He'  , 'H' , 'H+'], 
          'stoich' : [  -1   ,  -1  ,  1    ,  1  ,  1  ] },
        { 'spec'   : [ 'He+' , 'CO' , 'C+'  , 'He', 'O' ],
          'stoich' : [  -1   ,  -1  ,  1    ,  1  ,  1  ] },
        { 'spec'   : [ 'C+'  , 'H2' , 'CHx' , 'H'       ],
          'stoich' : [  -1   ,  -1  ,  1    ,  1        ] },
        { 'spec'   : [ 'C+'  , 'OHx', 'HCO+'            ],
          'stoich' : [  -1   ,  -1  ,  1                ] },
        { 'spec'   : [ 'O'   , 'CHx', 'CO'  , 'H'       ],
          'stoich' : [  -1   ,  -1  ,  1    ,  1        ] },
        { 'spec'   : [ 'C'   , 'OHx', 'CO'  , 'H'       ],
          'stoich' : [  -1   ,  -1  ,  1    ,  1        ] },
        { 'spec'   : [ 'He+' , 'e-' , 'He'              ],
          'stoich' : [  -1,     -1,    1                ] },
        { 'spec'   : [ 'H3+' , 'e-' , 'H2'  , 'H'       ],
          'stoich' : [  -1   ,  -1  ,  1    ,  1        ] },
        { 'spec'   : [ 'H3+' , 'e-' , 'H'               ],
          'stoich' : [  -1   ,  -1  ,  3                ] },
        { 'spec'   : [ 'C+'  , 'e-' , 'C'               ],
          'stoich' : [  -1   ,  -1  ,  1                ] },
        { 'spec'   : [ 'HCO+', 'e-' , 'CO'  , 'H'       ],
          'stoich' : [  -1   ,  -1  ,  1    ,  1        ] },
        { 'spec'   : [ 'H+'   , 'e-' , 'H'              ],
          'stoich' : [  -1   ,  -1  ,  1                ] },
        { 'spec'   : [ 'H2'  , 'H'  , 'H'               ],
          'stoich' : [  -1   ,  -1  ,  3                ] },
        { 'spec'   : [ 'H2'  , 'H2' , 'H'               ],
          'stoich' : [  -2   ,  1   ,  2                ] },
        { 'spec'   : [ 'H'   , 'e-' , 'H+'  , 'e-'      ],
          'stoich' : [  -1   ,  -1  ,  1    ,  2        ] },
        { 'spec'   : [ 'He+' , 'H2' , 'H+'  , 'He', 'H' ],
          'stoich' : [  -1   ,  -1  ,  1    ,  1  ,  1  ] },
        { 'spec'   : [ 'H'   , 'H2'                     ],
          'stoich' : [  -2   ,  1                       ],
          'grain'  : True                                 },
        { 'spec'   : [ 'H+'  , 'e-' , 'H'               ],
          'stoich' : [  -1   ,  -1  ,  1                ],
          'grain'  : True                                 },
        { 'spec'   : [ 'M+'  , 'e-' , 'M'               ],
          'stoich' : [  -1   ,  -1  ,  1                ] },
        { 'spec'   : [ 'H3+' , 'M'  , 'M+',   'H2', 'H' ],
          'stoich' : [  -1   ,  -1  ,  1  ,    1  ,  1  ] }
    ]
    _twobody = gr_reactions(specListExtended, _twobody_react)

def _twobody_ratecoef(T, Td, nH, nH2, ne, chi):
    """
    This returns the rate coefficients for all two-body reactions

    Parameters
    ----------
    T : float
       Gas temperature in K
    Td : float
       Dust temperature in K
    nH : float
       Number density of H in cm^-3
    nH2 : float
       Number density of H2 in cm^-3
    ne : float
       Number density of e- in cm^-3
    chi : float
       Radiation field normalized to Solar neighborhood value; this
       should already take into account the effects of dust
       attenuation

    Returns
    -------
    ratecoef : array(23)
       Rate coefficients for the 23 reactions in the network
    """

    # Sanitize inputs
    T = np.maximum(T, _small)
    Td = np.maximum(Td, _small)
    nH = np.maximum(nH, _small)
    nH2 = np.maximum(nH2, _small)
    ne = np.maximum(ne, _small)
    chi = np.maximum(chi, _small)

    # Some derived quantities used below
    logT = np.log10(T)
    logT4 = logT - 4.0
    T4 = T/1e4
    lnT = np.log(T)
    lnTe = np.log(8.6173e-5*T)  # ln of T measured in eV
    n = nH + 2*nH2
    xH = nH / n
    xH2 = 1.0 - xH

    ratecoef = np.zeros(len(_twobody_react))
    ratecoef[0] = 2.0e-9
    ratecoef[1] = 8.0e-10
    ratecoef[2] = 1.7e-9
    ratecoef[3] = 7.0e-15
    ratecoef[4] = 1.4e-9*(T/300.)**-0.5
    ratecoef[5] = 4.0e-16
    ratecoef[6] = 1.0e-9
    ratecoef[7] = 2.0e-10
    ratecoef[8] = 5.0e-12*T**0.5
    ratecoef[9] = 1.0e-11/np.sqrt(T) * \
                  (11.19 - 1.676*logT - 0.2852*logT**2 +
                   0.04433*logT**3)
    ratecoef[10] = 2.34e-8 * (T/300.)**-0.52
    ratecoef[11] = 4.36e-8 * (T/300.)**-0.52
    ratecoef[12] = 4.67e-12 * (T/300.)**-0.6
    ratecoef[13] = 2.76e-7 * (T/300.)**-0.64
    ratecoef[14] = 2.753e-14*(315614./T)**1.5 * \
                   (1.0+(115188./T)**0.407)**-2.242

    kHl = 6.67e-12*np.sqrt(T)*np.exp(-(1.0+63590./T))
    kHh = 3.52e-9*np.exp(-43900./T)
    ncrH = 10.**(3.0 - 0.416*logT4 - 0.327*logT4**2)
    ncrH2 = 10.**(4.845 - 1.3*logT4 + 1.62*logT4**2)
    ncr = 1.0/(xH/ncrH + xH2/ncrH2)
    ratecoef[15] = np.exp(
        (n/ncr) / (1.0+n/ncr) * np.log(kHh) +
        1.0 / (1.0+n/ncr) * np.log(kHl))
    kH2l = 5.996e-30*T**4.1881 / (1.0+6.761e-6*T)**5.6881 * \
           np.exp(-54657.4/T)
    kH2h = 1.9e-9*np.exp(-53300./T)
    ratecoef[16] = np.exp(
        (n/ncr) / (1.0+n/ncr) * np.log(kH2h) +
        1.0 / (1.0+n/ncr) * np.log(kH2l))
    ratecoef[17] = np.exp(
        - 32.71396786 + 13.5365560*lnTe - 5.73932875*lnTe**2
        + 1.56315498*lnTe**3 - 0.2877056*lnTe**4
        + 0.0348255977*lnTe**5 - 2.63197617e-3*lnTe**6
        + 1.11954395e-4*lnTe**7 - 2.03914985e-6*lnTe**8)
    ratecoef[18] = 3.7e-14*np.exp(-35./T)
    fA = 1.0/(1.0+1.0e4*np.exp(-600/Td))
    ratecoef[19] = 3.0e-18*T**0.5*fA / \
                   (1.0+0.04*(T+Td)**0.5 + 0.002*T+8.0e-6*T**2)
    psi = chi*np.sqrt(T)/ne
    ratecoef[20] = 12.25e-14 / \
                   (1.0 + 8.074e-6*psi**1.378 *
                    (1.0 + 508.7*T**0.01586*psi**(-0.4723-1.102e-5*lnT)))
    ratecoef[21] = 3.8e-10*T**-0.65
    ratecoef[22] = 2.0e-9

    # Return
    return ratecoef


########################################################################
# Set some default abundances
########################################################################
_xHedefault = 0.1
_xCdefault = 2.0e-4
_xOdefault = 4.0e-4
_xMdefault = 2.0e-7
_xHdefault = 1.0
_Zddefault = 1.0

########################################################################
# Define the NL99 class
########################################################################
class NL99_GC(chemNetwork):
    """
    This class the implements the CO chemistry network of Nelson & Langer
    (1999, ApJ, 524, 923) coupled to the H2 chemistry network of
    Glover & MacLow (2007, ApJS, 169, 239), as combined by Glover &
    Clark (2012, MNRAS, 421, 9).

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
       'Zd' : float
          dust abundance in solar units; defaults to 1.0
       'sigmaDustV' : float
          V band dust extinction cross
          section per H nucleus; if not set, the default behavior
          is to assume that sigmaDustV = 0.4 * cloud.dust.sigmaPE
       'AV' : float
          total visual extinction; ignored if sigmaDustV is set
       'noClump' : bool
          if True, the clumping factor is set to 1.0; defaults to False
       'sigmaNT' : float
          non-thermal velocity dispersion
       'temp' : float
          gas kinetic temperature
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

        Returns
           Nothing

        Remarks
           The dict info may contain the following key - value pairs:

           'xC' : float
              the total C abundance per H nucleus; defaults to 2.0e-4
           'xO' : float
              the total H abundance per H nucleus; defaults to 4.0e-4
           'xM' : float
              the total refractory metal abundance per H
              nucleus; defaults to 2.0e-7
           'Zd' : float
          dust abundance in solar units; defaults to 1.0
           'sigmaDustV' : float
              V band dust extinction cross
              section per H nucleus; if not set, the default behavior
              is to assume that sigmaDustV = 0.4 * cloud.dust.sigmaPE
           'AV' : float
              total visual extinction; ignored if sigmaDustV is set
           'noClump' : bool
              if True, the clumping factor is set to 1.0; defaults to False
           'sigmaNT' : float
              non-thermal velocity dispersion
           'temp' : float
              gas kinetic temperature
        """

        # List of species for this network; provide a pointer here so
        # that it can be accessed through the class
        self.specList = specList
        self.specListExtended = specListExtended

        # Store the input info dict
        self.info = info

        # Array to hold abundances; wrap it in an abundanceDict for
        # convenience
        self.x = np.zeros(len(specList))
        abd = abundanceDict(specList, self.x)

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
            self._sigmaNT = _small
            self._Zd = _Zddefault
            if info is not None:
                if 'AV' in info.keys():
                    self._AV = info['AV']
                if 'sigmaNT' in info.keys():
                    self._sigmaNT = info['sigmaNT']
                if 'Zd' in info.keys():
                    self._Zd = info['Zd']
                if 'Tg' in info.keys():
                    self._temp = info['Tg']
                if 'Td' in info.keys():
                    self._Td = info['Td']

            # Set initial abundances
            if info is None:
                # If not specied, start all hydrogen as H, all C as C+,
                # all O as OI
                abd['C+'] = _xCdefault
                abd['O'] = _xOdefault
            else:
                if 'xC' in info.keys():
                    abd['C+'] = info['xC']
                else:
                    abd['C+'] = _xCdefault
                if 'xO' in info.keys():
                    abd['O'] = info['xO']
                else:
                    abd['O'] = _xOdefault
            abd['M+'] = self.xM

        else:

            # Cloud is given, so get information out of it
            self.cloud = cloud

            # Sanity check: make sure cloud contains some He, since
            # network will not function properly at He abundance of 0
            if cloud.comp.xHe == 0.0:
                raise despoticError(
                    "NL99_GC network requires " + 
                    "non-zero He abundance")

            # Set abundances

            # Make a case-insensitive version of the emitter list for
            # convenience
            try:
                # This construction is elegant, but it relies on the
                # existence of a freestanding string.lower function,
                # which has been removed in python 3
                emList = dict(zip(map(string.lower, 
                                      cloud.emitters.keys()), 
                                  cloud.emitters.values()))
            except:
                # This somewhat more bulky construction is required in
                # python 3
                lowkeys = [k.lower() for k in cloud.emitters.keys()]
                lowvalues = list(cloud.emitters.values())
                emList = dict(zip(lowkeys, lowvalues))

            # Hydrogen
            abd['H+'] = cloud.comp.xHplus
            abd['H2'] = cloud.comp.xpH2 + cloud.comp.xoH2

            # OH and H2O
            if 'oh' in emList:
                abd['OHx'] += emList['oh'].abundance
            if 'ph2o' in emList:
                abd['OHx'] += emList['ph2o'].abundance
            if 'oh2o' in emList:
                abd['OHx'] += emList['oh2o'].abundance
            if 'p-h2o' in emList:
                abd['OHx'] += emList['p-h2o'].abundance
            if 'o-h2o' in emList:
                abd['OHx'] += emList['o-h2o'].abundance

            # CO
            if 'co' in emList:
                abd['CO'] = emList['co'].abundance

            # Neutral carbon
            if 'c' in emList:
                abd['C'] = emList['c'].abundance

            # Ionized carbon
            if 'c+' in emList:
                abd['C+'] = emList['c+'].abundance

            # HCO+
            if 'hco+' in emList:
                abd['HCO+'] = emList['hco+'].abundance

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
            xCtot = abd['CO'] + abd['C'] + abd['C+'] + abd['HCO+']
            if xCtot < xC:
                # Print warning if we're altering existing C+
                # abundance.
                if 'c' in emList:
                    print("Warning: input C abundance is " + 
                          str(xC) + ", but total input C, C+, CHx, CO, " + 
                          "HCO+ abundance is " + str(xCtot) + 
                          "; increasing xC+ to " + str(abd['C+']+xC-xCtot))
                abd['C+'] += xC - xCtot
            elif xCtot > xC:
                # Throw an error if input C abundance is smaller than
                # what is accounted for in initial conditions
                raise despoticError(
                    "input C abundance is " + 
                    str(xC) + ", but total input C, C+, CHx, CO, " + 
                    "HCO+ abundance is " + str(xCtot))

            # O
            if 'o' in emList:
                abd['O'] = emList['o'].abundance
            elif info is None:
                abd['O'] = _xOdefault - abd['OHx'] - abd['CO'] - \
                    abd['HCO+']
            elif 'xO' in info:
                abd['O'] = info['xO'] - abd['OHx'] - abd['CO'] - \
                    abd['HCO+']
            else:
                abd['O'] = _xOdefault - abd['OHx'] - abd['CO'] - \
                    abd['HCO+']

            # As with C, make sure all O is accounted for, and if not
            # park the extra in OI
            if info is None:
                xO = _xOdefault
            elif 'xC' in info:
                xO = info['xO']
            else:
                xO = _xOdefault
            xOtot = abd['OHx'] + abd['CO'] + abd['HCO+'] + abd['O']
            if xOtot < xO:
                # Print warning if we're altering existing O
                # abundance.
                if 'o' in emList:
                    print("Warning: input O abundance is " + 
                          str(xO) + ", but total input O, OHx, CO, " + 
                          "HCO+ abundance is " + str(xOtot) + 
                          "; increasing xO to " + str(abd['O']+xO-xOtot))
                abd['O'] += xO - xOtot
            elif xOtot > xO:
                # Throw an error if input O abundance is smaller than
                # what is accounted for in initial conditions
                raise despoticError(
                    "input O abundance is " + 
                    str(xO) + ", but total input O, OHx, CO, " + 
                    "HCO+ abundance is " + str(xOtot))

            # Finally, make sure all H nuclei are accounted for and
            # that all hydrogenic abundances are >= 0
            abd1 = self.abundances
            xH = abd1['H+'] + abd1['OHx'] + abd1['CHx'] + abd1['HCO+'] + \
                 abd1['H'] + 2*abd1['H2'] + 3*abd1['H3+']
            if (abs(xH-1.0) > 1.0e-8):
                raise despoticError(
                    "input hydrogen abundances " +
                    "add up to xH = "+str(xH)+" != 1!")
            if (abd1['H+'] < 0) or \
               (abd1['OHx'] < 0) or (abd1['CHx'] < 0) or \
               (abd1['HCO+'] < 0) or (abd1['H'] < 0) or \
               (abd1['H2'] < 0) or (abd1['H3+'] < 0):
                raise despoticError(
                    "abundances of some " + 
                    "hydrogenic species are < 0; abundances " + 
                    "are: " + repr(abd1))

        # Get rate coefficients in starting state; we will use these
        # to put some species close to equilibrium initially
        abd1 = self.abundances
        cfac = self.cfac
        n = self.nH * cfac
        rcoef = _twobody_ratecoef(
            self.temp, self.cloud.Td, n, abd1['H2']*n, abd1['e-']*n, 
            np.exp(-self.AV)*self.chi)

        # Set initial He+ to equilibrium value between creation by
        # cosmic rays and destruction by recombination with free
        # electrons, H2, and CO
        abd['He+'] \
            = self.xHe * self.ionRate / \
            (n * (abd1['H2'] * (rcoef[3]+rcoef[18]) +
                  abd1['CO'] * rcoef[4] +
                  abd1['e-'] * rcoef[9]))

        # Set initial H3+ in the same way as for He+; then correct H2
        # abundance to conserve total hydrogen
        abd['H3+'] \
            = abd1['H2'] * self.ionRate / \
            (n * (abd1['C'] * rcoef[0] +
                  abd1['O'] * rcoef[1] +
                  abd1['CO'] * rcoef[2] +
                  abd1['e-'] * (rcoef[10] + rcoef[11])))
        abd['H2'] -= 1.5*abd['H3+']

        # Initial M+
        abd['M+'] = self.xM


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
    def sigmaNT(self):
        """
        non-thermal velocity dispersion
        """
        if self.cloud is None:
            return self._sigmaNT
        else:
            return self.cloud.sigmaNT

    @sigmaNT.setter
    def sigmaNT(self, value):
        if self.cloud is None:
            self._sigmaNT = sigmaNT
        else:
            self.cloud.sigmaNT = sigmaNT

    @property
    def temp(self):
        if self.cloud is None:
            return self._temp
        else:
            return self.cloud.Tg

    @temp.setter
    def temp(self, value):
        """
        gas kinetic temperature
        """
        if self.cloud is None:
            self._temp = value
        else:
            self.cloud.Tg = value

    @property
    def cfac(self):
        """
        clumping factor; cannot be set directly, calculated from temp
        and sigmaNT
        """
        # Return 1.0 if cloud is not set, or if noClump is in info and
        # is set to True
        if self.cloud is None:
            return 1.0
        if self.info is not None:
            if 'noClump' in self.info.keys():
                if self.info['noClump']:
                    return 1.0

        # If we're here, compute mu and use that to get the sound
        # speed
        cs2 = kB * self.cloud.Tg / (self.mu() * mH)
        return np.sqrt(1.0 + 0.75*self.cloud.sigmaNT**2/cs2)

    @cfac.setter
    def cfac(self, value):
        raise despoticError("cannot set cfac directly")

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
    def Zd(self):
        """
        dust abundance normalized to the solar neighborhood value
        """
        if self.cloud is None:
            return self._Zd
        else:
            return self.cloud.dust.Zd

    @Zd.setter
    def Zd(self, value):
        if self.cloud is None:
            self._Zd = Zd
        else:
            self.cloud.dustProp.Zd = value

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
    # Return mean particle mass for a given chemical composition
    ####################################################################
    def mu(self, xin=None):
        """
        Return mean particle mass in units of H mass

        Parameters
           xin : array
              Chemical composition for which computation is to be
              done; if left as None, current chemical composition is
              used

        Returns
           mu : float
              Mean mass per free particle, in units of H mass
        """
        if xin is None:
            abd = self.abundances
        else:
            abd = self.extendAbundances(xin, outdict=True)
        mu = (abd['H+'] + abd['H'] + 2*abd['H2'] + 3*abd['H3+'] +
              4*(abd['He+']+abd['He'])) / \
            (abd['H+'] + abd['H'] + abd['H2'] + abd['H3+'] +
             abd['He+'] + abd['He'])
        return mu

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
        self._abundances = self.extendAbundances(outdict=True)
        return self._abundances

    @abundances.setter
    def abundances(self, other):
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
            warnings.warn("For NL99_GO network, cannot set abundances"
            " of derived species H, He, M, e-; abundances set only "
            " for other species")
        self.applyAbundances()

    ####################################################################
    # Method to get derived abundances from ones being stored; this
    # adds slots for H2, HeI, MI, and e
    ####################################################################
    def extendAbundances(self, xin=None, outdict=False):
        """
        Compute abundances of derived species not directly followed in
        the network.

        Parameters
           xin : array
              abundances of species directly tracked in the network;
              if left as None, the abundances stored internally to the
              network are used
           outdict : bool
              if True, the values are returned as an abundanceDict; if
              False, they are returned as a pure numpy array

        Returns
           x : array | abundanceDict
              abundances, including those of derived species
        """

        # Make abundance dict of the input species list
        if xin is None:
            abd = abundanceDict(specList, self.x)
        else:
            abd = abundanceDict(specList, xin)

        # Make object we'll be returning
        abd_out = abundanceDict(specListExtended, 
                                np.zeros(len(specListExtended)))
        if xin is None:
            abd_out.x[:len(specList)] = self.x
        else:
            abd_out.x[:len(specList)] = xin

        # Derive HI abundance from other types of H
        abd_out['H'] = 1.0 - abd['H+'] - 2*abd['H2'] - 3*abd['H3+'] \
                       - abd['OHx'] - abd['CHx'] - abd['HCO+']

        # He abundance = total He abundance - He+ abundance
        abd_out['He'] = self.xHe - abd['He+']

        # Neutral metal abundance = total metal abundance - ionized
        # metal abundance
        abd_out['M'] = self.xM - abd['M+']

        # e abundance = H+ + H3+ + He+ + C+ + HCO+ + M+
        abd_out['e-'] = abd['H+'] + abd['H3+'] + abd['He+'] + \
                        abd['C+'] + abd['HCO+'] + abd['M+']

        # Return
        if outdict:
            return abd_out
        else:
            return abd_out.x


    ####################################################################
    # Method to return the time derivative of all chemical rates
    ####################################################################
    def dxdt(self, xin, time):
        """
        This method returns the time derivative of all abundances for
        this chemical network.

        Parameters
           xin : array(12)
              current abundances of all species
           time : float
              current time; not actually used, but included as an
              argument for compatibility with odeint

        Returns
           dxdt : array(12)
              time derivative of x
        """

        # Get abundances of derived quantities
        abd = self.extendAbundances(xin, outdict=True)

        # Get mean particle mass and velocity dispersion; for particle
        # mass, for simplicity we ignore all species but H, He, H2,
        # and e-, and their ions
        sigma_tot = np.sqrt(self.sigmaNT**2 + kB*self.temp 
                            / (self.mu(xin)*mH))

        # Get clumping factor and effective density
        cfac = self.cfac
        n = self.nH * cfac

        # Use the cosmic ray and photoreaction rate calculators to get
        # their contributions to dxdt
        xdot = _cr.dxdt(abd.x, n, self.ionRate)
        xdot += _ph.dxdt(abd.x, n, self.chi, self.AV, 
                         [[self.NH*abd['H2'], sigma_tot],
                          [self.NH*abd['CO'], self.NH*abd['H2']]])

        # Get rates for two-body reactions and add their contribution
        # to xdot
        ratecoef = _twobody_ratecoef(
            self.temp, self.cloud.Td, n, abd['H2']*n, abd['e-']*n, 
            np.exp(-self.AV)*self.chi)
        xdot += _twobody.dxdt(abd.x, n, ratecoef, self.Zd)

        # Return results, omitting derived species
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

        # Get the current abundances as an abundanceDict
        abd = self.abundances

        # Hydrogen; leave the ortho-to-para ratio unchanged, or set it
        # to 0.25 if the initial value is undefined
        if self.cloud.comp.xH2 != 0:
            fortho = self.cloud.comp.xoH2 / self.cloud.comp.xH2
        else:
            fortho = 1./3.
        self.cloud.comp.xHI = abd['H']
        self.cloud.comp.xHplus = abd['H+']
        self.cloud.comp.xoH2 = abd['H2'] * fortho
        self.cloud.comp.xpH2 = abd['H2'] * (1.0 - fortho)
        self.cloud.comp.xe = abd['e-']

        # Hydrogen safety check: DESPOTIC doesn't like it if all
        # hydrogen isn't accounted for, and it doesn't know about H2+
        # or H3+; in certain odd parts of parameter space this can
        # cause problems, which we fix by marginally increasing the
        # abundance of the most abundance hydrogenic species in the
        # composition of the cloud
        try:
            self.cloud.comp._check_abundance()
        except ValueError:
            self.cloud.comp.xHplus = 1.0 - self.cloud.comp.xHI - \
                                     2.0*self.cloud.comp.xH2

        # Make a case-insensitive version of the emitter list for
        # convenience
        try:
            # This construction is elegant, but it relies on the
            # existence of a freestanding string.lower function,
            # which has been removed in python 3
            emList = dict(zip(map(string.lower, 
                                  self.cloud.emitters.keys()), 
                              self.cloud.emitters.values()))
        except:
            # This somewhat more bulky construction is required in
            # python 3
            lowkeys = [k.lower() for k in self.cloud.emitters.keys()]
            lowvalues = list(self.cloud.emitters.values())
            emList = dict(zip(lowkeys, lowvalues))

        # Save ratios of ^12C to ^13C, and ^16O to ^18O
        if '13co' in emList and 'co' in emList:
            c13_12 = emList['13co'].abundance / \
                emList['co'].abundance
        if 'c18o' in emList and 'co' in emList:
            o18_16 = emList['c18o'].abundance / \
                emList['co'].abundance

        # OH, assuming OHx is half OH
        if 'oh' in emList:
            emList['oh'].abundance = abd['OHx']/2.0
        elif addEmitters:
            try:
                self.cloud.addEmitter('oh', abd['OHx']/2.0)
            except despoticError:
                print('Warning: unable to add OH; cannot find LAMDA file')

        # H2O, assuming OHx is half H2O, and that oH2O and pH2O are in
        # the same ratio as oH2 and pH2
        if 'ph2o' in emList:
            emList['ph2o'].abundance = abd['OHx']/2.0*(1-fortho)
        elif 'p-h2o' in emList:
            emList['p-h2o'].abundance = abd['OHx']/2.0*(1-fortho)
        elif addEmitters:
            try:
                self.cloud.addEmitter('ph2o', 
                                      abd['OHx']/2.0*(1-fortho))
            except despoticError:
                print('Warning: unable to add p-H2O; cannot find LAMDA file')
        if 'oh2o' in emList:
            emList['oh2o'].abundance = abd['OHx']/2.0*fortho
        elif 'o-h2o' in emList:
            emList['o-h2o'].abundance = abd['OHx']/2.0*fortho
        elif addEmitters:
            try:
                self.cloud.addEmitter('oh2o', abd['OHx']/2.0*fortho)
            except despoticError:
                print('Warning: unable to add o-H2O; cannot find LAMDA file')

        # CO
        if 'co' in emList:
            emList['co'].abundance = abd['CO']
        elif addEmitters:
            try:
                self.cloud.addEmitter('co', abd['CO'])
            except despoticError:
                print('Warning: unable to add CO; cannot find LAMDA file')

        # if we have 13CO or C18O, make their abundances match that of CO
        # multiplied by the appropriate isotopic abundances
        if '13co' in emList:
            emList['13co'].abundance = abd['CO']*c13_12
        if 'c18o' in emList:
            emList['c18o'].abundance = abd['CO']*o18_16

        # C
        if 'c' in emList:
            emList['c'].abundance = abd['C']
        elif addEmitters:
            try:
                self.cloud.addEmitter('c', abd['C'])
            except despoticError:
                print('Warning: unable to add C; cannot find LAMDA file')

        # C+
        if 'c+' in emList:
            emList['c+'].abundance = abd['C+']
        elif addEmitters:
            try:
                self.cloud.addEmitter('c+', abd['C+'])
            except despoticError:
                print('Warning: unable to add C+; cannot find LAMDA file')

        # HCO+
        if 'hco+' in emList:
            emList['hco+'].abundance = abd['HCO+']
        elif addEmitters:
            try:
                self.cloud.addEmitter('hco+', abd['HCO+'])
            except despoticError:
                print('Warning: unable to add HCO+; cannot find LAMDA file')

        # O
        if 'o' in emList:
            emList['o'].abundance = abd['O']
        elif addEmitters:
            try:
                self.cloud.addEmitter('o', abd['O'])
            except despoticError:
                print('Warning: unable to add O; cannot find LAMDA file')

