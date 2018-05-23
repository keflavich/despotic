"""
This module implements the H - C - O chemistry network of Gong,
Ostriker, & Wolfire, 2017, ApJ, 843, 38.
"""

########################################################################
# Copyright (C) 2018 Mark Krumholz
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
from .shielding import fShield_CO_VvDB, fShield_C_TH, fShield_H2_DB
from .abundanceDict import abundanceDict
from .chemNetwork import chemNetwork
from .reactions import cr_reactions, photoreactions, gr_reactions
import scipy.constants as physcons
import warnings

########################################################################
# Physical and numerical constants
########################################################################
kB = physcons.k/physcons.erg
mH = (physcons.m_p+physcons.m_e)/physcons.gram
_small = 1e-100

########################################################################
# List of species used in this chemistry network; specList is the set
# of species that are actually evolved, while specListExtended
# includes both species that are evolved and those whose abundances
# are derived from conservation laws.
########################################################################
specList = ['H2', 'H+', 'H2+', 'H3+', 'He+', 'O+', 'C+', 'CO',
            'HCO+', 'Si+', 'CHx', 'OHx']
specListExtended = specList + ['H', 'He', 'C', 'O', 'Si', 'e-']

########################################################################
# Data on cosmic ray reactions
# Reactions are, in order:
# (1)       cr + H       -> H+   + e-
# (2)       cr + H2      -> H2+  + e-
# (3)       cr + He      -> He+  + e-
# (4)       cr + C       -> C+   + e-
# (5)       cr + CO + H  -> HCO+ + e-
# (6) gamma_cr + C       -> C+   + e-
# (7) gamma_cr + CO      -> C    + O
# (8) gamma_cr + Si      -> Si+  + e-
# Note: gamma_cr = cosmic-ray induced secondary photon
########################################################################
_cr_react = [
    # Note that we handle reactions (1) and (2),
    # cr + H -> H+ + e-
    # cr + H2 -> H2+ + e-
    # elsewhere because the rates include both primary and secondary
    # ionizations, which the CR rate machinery doesn't handle. The
    # reactions included below are only primary ionizations.
    # cr + He -> He+ + e-
    { 'spec'      : ['He', 'He+', 'e-'],
      'stoich'    : [  -1,     1,    1],
      'rate'      : 1.1 },
    # cr + C -> C+ + e-
    { 'spec'      : ['C', 'C+', 'e-'],
      'stoich'    : [ -1,    1,    1],
      'rate'      : 3.85 },
    # cr + CO -> HCO+ + e-
    { 'spec'      : ['CO', 'HCO+', 'e-'],
      'stoich'    : [  -1,      1,    1],
      'rate'      : 6.52 },
    # gamma_cr + C -> C+ + e-
    { 'spec'      : ['C', 'C+', 'e-'],
      'stoich'    : [ -1,    1,    1],
      'rate'      : 560. },
    # gamma_cr + CO -> C + O
    { 'spec'      : ['CO', 'C', 'O'],
      'stoich'    : [  -1,   1,   1],
      'rate'      : 90. },
    # gamma_cr + Si -> Si+ + e-
    { 'spec'      : ['Si', 'Si+', 'e-'],
      'stoich'    : [  -1,     1,    1],
      'rate'      : 8400. }
    ]
_cr = cr_reactions(specListExtended, _cr_react)

########################################################################
# Data on photoreactions
# Reactions are, in order:
# (1) gamma + C    -> C+  + e-
# (2) gamma + CHx  -> C   + H
# (3) gamma + CO   -> C   + O
# (4) gamma + OHx  -> O   + H
# (5) gamma + Si   -> Si+ + e-
# (6) gamma + H2   -> H   + H
########################################################################
_ph_react = [
    # gamma + C -> C+ + e-
    { 'spec'       : ['C', 'C+', 'e-'],
      'stoich'     : [ -1,    1,    1],
      'rate'       : 3.5e-10,
      'av_fac'     : 3.76,
      'shield_fac' : fShield_C_TH },
    # gamma + CHx -> C + H
    { 'spec'       : ['CHx', 'C', 'H'],
      'stoich'     : [   -1,   1,   1],
      'rate'       : 9.1e-10,
      'av_fac'     : 2.12 },
    # gamma + CO -> C + O
    { 'spec'       : ['CO', 'C', 'O'],
      'stoich'     : [  -1,   1,   1],
      'rate'       : 2.4e-10,
      'av_fac'     : 3.88,
      'shield_fac' : fShield_CO_VvDB },
    # gamma + OHx -> O + H
    { 'spec'       : ['OHx', 'O', 'H'],
      'stoich'     : [   -1,   1,   1],
      'rate'       : 3.8e-10,
      'av_fac'     : 2.66 },
    # gamma + Si -> Si+ + e-
    { 'spec'       : ['Si', 'Si+', 'e-'],
      'stoich'     : [   -1,    1,    1],
      'rate'       : 4.5e-9,
      'av_fac'     : 2.61 },
    # gamma + H2 -> H + H
    { 'spec'       : ['H2', 'H'],
      'stoich'     : [  -1,   2],
      'rate'       : 5.7e-11,
      'av_fac'     : 4.18,
      'shield_fac' : fShield_H2_DB }
    ]
_ph = photoreactions(specListExtended, _ph_react)


########################################################################
# Data on two-body reactions
# Reactions are, in order:
# (1) H3+ + C -> CHx + H2
# (2) H3+ + O -> OHx + H2
# (3) H3+ + O + e- -> H2 + O + H
# (4) O+ + H2 -> OHx + H
# (5) O+ + H2 + e- > O + H + H
# (6) H3+ + CO -> HCO+ + H2
# (7) He+ + H2 -> H+ + He + H
# (8) He+ + CO -> C+ + O + He
# (9) C+ + H2 -> CHx + H
# (10) C+ + H2 + e- -> C + H + H
# (11) C+ + OHx -> HCO+
# (12) CHx + O -> CO + H
# (13) OHx + C -> CO + H
# (14) He+ + e- -> He
# (15) H3+ + e- -> H2 + H
# (16) H3+ + e- -> 3H
# (17) C+ + e- -> C
# (18) HCO+ + e- -> CO + H
# (19) H2+ + H2 -> H3+ + H
# (20) H2+ + H -> H+ + H2
# (21) H+ + e- -> H
# (22) H2 + H -> 3H
# (23) H2 + H2 -> H2 + 2H
# (24) H + e- -> H+ + e-
# (25) He+ + H2 -> H2+ + He
# (26) CHx + H -> H2 + C
# (27) OHx + O -> 2O + H
# (28) Si+ + e- -> Si
# (29) He+ + OHx -> O+ + He + H
# (30) H+ + O -> O+ + H
# (31) O+ + H -> H+ + O
# (32) H + H -> H2 (grain catalyzed)
# (33) H+ + e- -> H (grain catalyzed)
# (34) C+ + e- -> C (grain catalyzed)
# (35) He+ + e- -> He (grain catalyzed)
# (36) Si+ + e- -> Si (grain catalyzed)
########################################################################
_twobody_react = [
    # (1) H3+ + C -> CHx + H2
    { 'spec'   : ['H3+', 'C', 'CHx', 'H2'],
      'stoich' : [   -1,  -1,     1,    1] },
    # (2) H3+ + O -> OHx + H2
    { 'spec'   : ['H3+', 'O', 'OHx', 'H2'],
      'stoich' : [   -1,  -1,     1,    1] },
    # (3) H3+ + O + e- -> H2 + O + H
    { 'spec'   : ['H3+', 'O', 'H2', 'O', 'H'],
      'stoich' : [   -1,  -1,    1,   1,   1] },
    # (4) O+ + H2 -> OHx + H
    { 'spec'   : [ 'O+', 'H2', 'OHx', 'H'],
      'stoich' : [   -1,   -1,     1,   1] },
    # (5) O+ + H2 + e- > O + H + H
    { 'spec'   : ['O+', 'H2', 'O', 'H'],
      'stoich' : [  -1,   -1,   1,   2] },
    # (6) H3+ + CO -> HCO+ + H2
    { 'spec'   : [ 'H3+', 'CO', 'HCO+', 'H2'],
      'stoich' : [    -1,   -1,      1,    1] },
    # (7) He+ + H2 -> H+ + He + H
    { 'spec'   : [ 'He+', 'H2', 'H+', 'He', 'H' ],
      'stoich' : [    -1,   -1,    1,    1,   1] },
    # (8) He+ + CO -> C+ + O + He
    { 'spec'   : [ 'He+', 'CO', 'C+', 'O', 'He'],
      'stoich' : [    -1,   -1,    1,   1,    1] },
    # (9) C+ + H2 -> CHx + H
    { 'spec'   : [ 'C+', 'H2', 'CHx', 'H'],
      'stoich' : [   -1,   -1,     1,   1] },
    # (10) C+ + H2 + e- -> C + H + H
    { 'spec'   : [ 'C+', 'H2', 'C', 'H'],
      'stoich' : [   -1,   -1,   1,   2] },
    # (11) C+ + OHx -> HCO+
    { 'spec'   : [ 'C+', 'OHx', 'HCO+' ],
      'stoich' : [   -1,    -1,      1] },
    # (12) CHx + O -> CO + H
    { 'spec'   : [ 'CHx', 'O', 'CO', 'H'],
      'stoich' : [    -1,  -1,    1,   1] },
    # (13) OHx + C -> CO + H
    { 'spec'   : [ 'OHx', 'C', 'CO', 'H'],
      'stoich' : [    -1,  -1,    1,   1] },
    # (14) He+ + e- -> He
    { 'spec'   : [ 'He+', 'e-', 'He'],
      'stoich' : [    -1,   -1,    1] },
    # (15) H3+ + e- -> H2 + H
    { 'spec'   : [ 'H3+', 'e-', 'H2', 'H' ],
      'stoich' : [    -1,   -1,    1,   1 ] },
    # (16) H3+ + e- -> 3H
    { 'spec'   : [ 'H3+', 'e-', 'H'],
      'stoich' : [    -1,   -1,   3] },
    # (17) C+ + e- -> C
    { 'spec'   : [ 'C+', 'e-', 'C' ],
      'stoich' : [   -1,   -1,   1 ] },
    # (18) HCO+ + e- -> CO + H
    { 'spec'   : [ 'HCO+', 'e-', 'CO', 'H'],
      'stoich' : [     -1,   -1,    1,   1] },
    # (19) H2+ + H2 -> H3+ + H
    { 'spec'   : [ 'H2+', 'H2', 'H3+', 'H'],
      'stoich' : [    -1,   -1,     1,   1] },
    # (20) H2+ + H -> H+ + H2
    { 'spec'   : [ 'H2+', 'H', 'H+', 'H2'],
      'stoich' : [    -1,  -1,    1,    1] },
    # (21) H+ + e- -> H
    { 'spec'   : [ 'H+', 'e-', 'H'],
      'stoich' : [   -1,   -1,   1] },
    # (22) H2 + H -> 3H
    { 'spec'   : [ 'H2', 'H', 'H'],
      'stoich' : [   -1,  -1,   2] },
    # (23) H2 + H2 -> H2 + 2H
    { 'spec'   : [ 'H2', 'H2', 'H'],
      'stoich' : [   -2,    1,   2] },
    # (24) H + e- -> H+ + e-
    { 'spec'   : [ 'H', 'e-', 'H+', 'e-'],
      'stoich' : [  -1,   -1,    1,    1] },
    # (25) He+ + H2 -> H2+ + He
    { 'spec'   : [ 'He+', 'H2', 'H2+', 'He'],
      'stoich' : [    -1,   -1,     1,    1] },
    # (26) CHx + H -> H2 + C
    { 'spec'   : [ 'CHx', 'H', 'H2', 'C' ],
      'stoich' : [    -1,  -1,   1,   1] },
    # (27) OHx + O -> 2O + H
    { 'spec'   : [ 'OHx', 'O', 'O', 'H'],
      'stoich' : [    -1,  -1,   2,   1] },
    # (28) Si+ + e- -> Si
    { 'spec'   : [ 'Si+', 'e-', 'Si'],
      'stoich' : [    -1,   -1,    1] },
    # (29) He+ + OHx -> O+ + He + H
    { 'spec'   : [ 'He+', 'OHx', 'O+', 'He', 'H' ],
      'stoich' : [    -1,    -1,    1,    1,   1 ] },
    # (30) H+ + O -> O+ + H
    { 'spec'   : [ 'H+', 'O', 'O+', 'H'],
      'stoich' : [   -1,  -1,    1,   1] },
    # (31) O+ + H -> H+ + O
    { 'spec'   : [ 'O+', 'H', 'H+', 'O' ],
      'stoich' : [   -1,  -1,    1,   1 ] },
    # (32) H + H -> H2 (grain catalyzed)
    { 'spec'   : [ 'H', 'H2'],
      'stoich' : [  -2,    1],
      'grain'  : True },
    # (33) H+ + e- -> H (grain catalyzed)
    { 'spec'   : [ 'H+', 'e-', 'H'],
      'stoich' : [   -1,   -1,   1],
      'grain'  : True },
    # (34) C+ + e- -> C (grain catalyzed)
    { 'spec'   : [ 'C+', 'e-', 'C'],
      'stoich' : [   -1,   -1,   1],
      'grain'  : True },
    # (35) He+ + e- -> He (grain catalyzed)
    { 'spec'   : [ 'He', 'e-', 'He' ],
      'stoich' : [   -1,   -1,    1 ],
      'grain'  : True },
    # (36) Si+ + e- -> Si (grain catalyzed)
    { 'spec'   : [ 'Si+', 'e-', 'Si' ],
      'stoich' : [    -1,   -1,    1 ],
      'grain'  : True }
]
_twobody = gr_reactions(specListExtended, _twobody_react)

def _twobody_ratecoef(T, n, nH, nH2, ne, chi):
    """
    This returns the rate coefficients for all two-body reactions

    Parameters
    ----------
    T : float
       Gas temperature in K
    n : float
       Number density of H nuclei in cm^-3
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
    ratecoef : array
       Rate coefficients for the reactions in the network
    """
    
    # Sanitize inputs
    T = np.maximum(T, _small)
    nH = np.maximum(nH, _small)
    nH2 = np.maximum(nH2, _small)
    ne = np.maximum(ne, _small)
    chi = np.maximum(chi, _small)

    # Some derived quantities
    Te = 8.6173e-5*T     # Temperature in units of eV/k_B
    logT = np.log10(T)
    logT4 = logT - 4.0
    lnT = np.log(T)
    lnTe = np.log(Te)
    xH = nH / n
    xH2 = nH2 / n
    psi = 1.7*np.sqrt(T)/ne

    # Rate coefficients for reactions
    ratecoef = np.zeros(len(_twobody_react))
    
    # (1) H3+ + C -> CHx + H2
    ci = np.array([3.4e-8, 6.97e-9, 1.31e-7, 1.51e-4])
    Ti = np.array([7.62, 1.38, 2.66e1, 8.11e3])
    ratecoef[0] = 1.04e-9*(300.0/T)**0.00231 + \
                  T**1.5 * np.sum(ci*np.exp(-Ti/T))

    # (2) H3+ + O -> OHx + H2
    kH2OpH2 = 6.0e-10
    kH2Ope = 5.3e-6 / np.sqrt(T)
    r = kH2OpH2 * nH2 / (kH2OpH2 * nH2 + kH2Ope * ne)
    ratecoef[1] = 1.99e-9 * T**-0.190 * r
    
    # (3) H3+ + O + e- -> H2 + O + H
    ratecoef[2] = 1.99e-9 * T**-0.190 * (1.0 - r)

    # (4) O+ + H2 -> OHx + H
    ratecoef[3] = 1.6e-9 * r
    
    # (5) O+ + H2 + e- > O + H + H
    ratecoef[4] = 1.6e-9 * (1.0 - r)
    
    # (6) H3+ + CO -> HCO+ + H2
    ratecoef[5] = 1.7e-9
    
    # (7) He+ + H2 -> H+ + He + H
    ratecoef[6] = 1.26e-13 * np.exp(-22.5/T)
    
    # (8) He+ + CO -> C+ + O + He
    ratecoef[7] = 1.6e-9
    
    # (9) C+ + H2 -> CHx + H
    ratecoef[8] = 2.31e-13 * T**-1.3 * np.exp(-23./T)
    
    # (10) C+ + H2 + e- -> C + H + H
    ratecoef[9] = 0.99e-13 * T**-1.3 * np.exp(-23./T)
    
    # (11) C+ + OHx -> HCO+
    ratecoef[10] = 9.15e-10 * (0.62 + 45.41/np.sqrt(T))
    
    # (12) CHx + O -> CO + H
    ratecoef[11] = 7.7e-11
    
    # (13) OHx + C -> CO + H
    ratecoef[12] = 7.95e-10 * T**-0.339 * np.exp(0.108/T)
    
    # (14) He+ + e- -> He
    ratecoef[13] = 1.0e-11/np.sqrt(T) * (11.19 - 1.676*logT
                                         - 0.2852*logT**2
                                         + 0.04433*logT**3)
    
    # (15) H3+ + e- -> H2 + H
    ratecoef[14] = 4.54e-7 * T**-0.52
    
    # (16) H3+ + e- -> 3H
    ratecoef[15] = 8.46e-7 * T**-0.52
    
    # (17) C+ + e- -> C
    alpha = np.sqrt(T/6.67e-3)
    beta = np.sqrt(T/1.943e6)
    gamma = 0.7849+0.1597*np.exp(-49550./T)
    krr = 2.995e-9 / (alpha*(1.0+alpha)**(1.0-gamma) *
                      (1.0+beta)**(1.0+gamma))
    kdr = T**-1.5 * (6.346e-9 * np.exp(-12.17/T) +
                     9.793e-9*np.exp(-73.8/T) +
                     1.634e-6*np.exp(-15230./T))
    ratecoef[16] = krr + kdr
    
    # (18) HCO+ + e- -> CO + H
    ratecoef[17] = 1.06e-5 * T**-0.64
    
    # (19) H2+ + H2 -> H3+ + H
    ratecoef[18] = 2.84e-9 * T**0.042 * np.exp(-T/46600.)
    
    # (20) H2+ + H -> H+ + H2
    ratecoef[19] = 6.4e-10
    
    # (21) H+ + e- -> H
    ratecoef[20] = 2.753e-14 * (315614./T)**1.5 * \
                   (1.0 + (115188./T)**0.407)**-2.242
    
    # (22) H2 + H -> 3H
    kHl = 6.67e-12*np.sqrt(T)*np.exp(-(1.0+63590./T))
    kHh = 3.52e-9*np.exp(-43900./T)
    ncrH = 10.**(3.0 - 0.416*logT4 - 0.327*logT4**2)
    ncrH2 = 10.**(4.845 - 1.3*logT4 + 1.62*logT4**2)
    ncr = 1.0/(xH/ncrH + xH2/ncrH2)
    if kHh > 0.0 and kHl > 0.0:
        ratecoef[21] = np.exp(
            (n/ncr) / (1.0+n/ncr) * np.log(kHh) +
            1.0 / (1.0+n/ncr) * np.log(kHl))
    
    # (23) H2 + H2 -> H2 + 2H
    kH2l = 5.996e-30*T**4.1881 / (1.0+6.761e-6*T)**5.6881 * \
           np.exp(-54657.4/T)
    kH2h = 1.9e-9*np.exp(-53300./T)
    if kH2l > 0 and kH2h > 0:
        ratecoef[22] = np.exp(
            (n/ncr) / (1.0+n/ncr) * np.log(kH2h) +
            1.0 / (1.0+n/ncr) * np.log(kH2l))
    
    # (24) H + e- -> H+ + e-
    ratecoef[23] = np.exp(
        -3.271396786e1
        +1.35365560e1*lnTe 
        -5.73932875*lnTe**2
        +1.56315498*lnTe**3
        -2.877056e-1*lnTe**4
        +3.48255977e-2*lnTe**5
        -2.63197617e-3*lnTe**6
        +1.11954395e-4*lnTe**7
        -2.03914985e-6*lnTe**8)
    
    # (25) He+ + H2 -> H2+ + He
    ratecoef[24] = 7.2e-15
    
    # (26) CHx + H -> H2 + C
    ratecoef[25] = 2.81e-11*T**0.26
    
    # (27) OHx + O -> 2O + H
    ratecoef[26] = 3.5e-11
    
    # (28) Si+ + e- -> Si
    ratecoef[27] = 1.46e-10*T**-0.62
    
    # (29) He+ + OHx -> O+ + He + H
    ratecoef[28] = 1.35e-9*(0.62+45.41/np.sqrt(T))
    
    # (30) H+ + O -> O+ + H
    ratecoef[29] = (1.1e-11*T**0.517+4.0e-10*T**0.00669)*np.exp(-227./T)
    
    # (31) O+ + H -> H+ + O
    ratecoef[30] = 4.99e-11*T**0.405+7.5e-10*T**-0.458
    
    # (32) H + H -> H2 (grain catalyzed)
    ratecoef[31] = 3.0e-17
    
    # (33) H+ + e- -> H (grain catalyzed)
    ratecoef[32] = 12.25e-14 / (
        1.0 + 8.074e-6*psi**1.378 *
        (1. + 508.7*T**0.01586*psi**(-0.4723-1.102e-5*lnT)))
    
    # (34) C+ + e- -> C (grain catalyzed)
    ratecoef[33] = 45.58e-14 / (
        1.0 + 6.089e-3*psi**1.128 *
        (1.0 + 433.1*T**0.4845*psi**(-0.8120-1.333e-4*lnT)))
    
    # (35) He+ + e- -> He (grain catalyzed)
    ratecoef[34] = 5.572e-14 / (
        1.0 + 3.185e-7*psi**1.512 *
        (1.0 + 5115.*T**3.903e-7*psi**(-0.4956-5.494e-7*lnT)))
    
    # (36) Si+ + e- -> Si (grain catalyzed)
    ratecoef[35] = 2.166e-14 / (
        1.0 + 5.678e-8*psi**1.874 *
        (1.0 + 43750.*T**1.635e-6*psi**(-0.8964-7.538e-5*lnT)))

    # Return
    return ratecoef


########################################################################
# Set some default abundances
########################################################################
_xHedefault = 0.1
_xCdefault = 1.6e-4
_xOdefault = 3.2e-4
_xSidefault = 1.7e-6
_xHdefault = 1.0
_Zddefault = 1.0


########################################################################
# Define the GOW chemistry class
########################################################################
class GOW(chemNetwork):
    """
    This module implements the H - C - O chemistry network of Gong,
    Ostriker, & Wolfire, 2017, ApJ, 843, 38.

    Parameters
       cloud : class cloud
         a DESPOTIC cloud object from which initial data are to be
         taken
       info : dict
         a dict containing additional parameters

    Remarks
       The dict info may contain the following key - value pairs:

       'xC' : float
          the total C abundance per H nucleus; defaults to 1.6e-4
       'xO' : float
          the total H abundance per H nucleus; defaults to 3.2e-4
       'xSi' : float
          the total Si abundance per H nucleus; defaults to 1.7e-6
       'Zd' : float
          dust abundance in solar units; defaults to 1.0
       'sigmaDustV' : float
          V band dust extinction cross section per H nucleus; if not
          set, the default behavior is to assume that sigmaDustV = 0.4
          * cloud.dust.sigmaPE
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

           'xHe' : float
               the total H abundance per H nucleus; defaults to 0.1
           'xC' : float
               the total C abundance per H nucleus; defaults to 1.6e-4
           'xO' : float
               the total H abundance per H nucleus; defaults to 3.2e-4
           'xSi' : float
               the total Si abundance per H nucleus; defaults to 1.7e-6
           'Zd' : float
               dust abundance in solar units; defaults to 1.0
           'sigmaDustV' : float
               V band dust extinction cross section per H nucleus; if not
               set, the default behavior is to assume that sigmaDustV = 0.4
               * cloud.dust.sigmaPE
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

        # Set total elemental abundances
        if info is None:
            self.xHe = _xHedefault
            self.xO = _xOdefault
            self.xC = _xCdefault
            self.xSi = _xSidefault
        else:
            if 'xHe' in info.keys():
                self.xHe = info['xHe']
            else:
                self.xHe = _xHedefault
            if 'xO' in info.keys():
                self.xO = info['xO']
            else:
                self.xO = _xOdefault
            if 'xC' in info.keys():
                self.xC = info['xC']
            else:
                self.xC = _xCdefault
            if 'xSi' in info.keys():
                self.xSi = info['xSi']
            else:
                self.xSi = _xSidefault

        # Extract information from the cloud if one is given
        if cloud is None:

            # No cloud given, so set some defaults
            self.cloud = None

            # Physical properties
            self._xHe = _xHedefault
            self._ionRate = 2.0e-16
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

            # Start all H as H, all He as He, all C as C+, all O as O,
            # all Si as Si+
            if info is None:
                # If not specied, start all hydrogen as H, all C as C+,
                # all O as OI, all Si as Si+
                abd['H'] = 1.0
                abd['He'] = self.xHe
                abd['C+'] = self.xC
                abd['Si+'] = self.xSi

        else:

            # Cloud is given, so get information out of it
            self.cloud = cloud

            # Sanity check: make sure cloud contains some He, since
            # network will not function properly at He abundance of 0
            if cloud.comp.xHe == 0.0:
                raise despoticError(
                    "GOW network requires " + 
                    "non-zero He abundance")
        
            # Set abundances from emitter list

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

            # C+
            if 'c+' in emList:
                abd['C+'] = emList['c+'].abundance

            # HCO+
            if 'hco+' in emList:
                abd['HCO+'] = emList['hco+'].abundance

            # He+
            if 'he+' in emList:
                abd['He+'] = emList['he+'].abundance

            # Si+
            if 'si+' in emList:
                abd['Si+'] = emList['si+'].abundance

            # O+
            if 'o+' in emList:
                abd['O+'] = emList['o+'].abundance

            # Sanity check: make sure input abundances do not exceed
            # total elemental abundance for any element; if it is,
            # throw error

            # H
            if 2*abd['H2'] + abd['H+'] + 2*abd['H2+'] + 3*abd['H3+'] + \
               abd['HCO+'] + abd['CHx'] + abd['OHx'] > 1:
                raise ValueError(
                    "total H elemental abundance (from H+, H2, H2+, "
                    "H3+, HCO+, CHx, OHx) > 1 per H nucleus; "
                    "exiting due to inconsistent elemental "
                    "abundances")

            # C
            xCtot = abd['C+'] + abd['CO'] + abd['HCO+'] + abd['CHx']
            if xCtot > self.xC:
                raise ValueError(
                    "total C elemental abundance (from C+, CO, "
                    "HCO+, CHx) x_C = {:e}".format(xCtot) +
                    " exceeds input C elemental abundance {:e}".
                    format(self.xC) + "; exiting due to inconsistent "
                    "elemental abundances")

            # O
            xOtot = abd['O+'] + abd['CO'] + abd['HCO+'] + abd['OHx']
            if xOtot > self.xO:
                raise ValueError(
                    "total O elemental abundance (from O+, CO, "
                    "HCO+, OHx) x_C = {:e}".format(xOtot) +
                    " exceeds input O elemental abundance {:e}".
                    format(self.xO) + "; exiting due to inconsistent "
                    "elemental abundances")
 
            # Si
            xSitot = abd['Si+']
            if xSitot > self.xSi:
                raise ValueError(
                    "total Si elemental abundance (from Si+) "
                    "x_Si = {:e}".format(xSitot) +
                    " exceeds input Si elemental abundance {:e}".
                    format(self.xSi) + "; exiting due to inconsistent "
                    "elemental abundances")



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
    # Method to get derived abundances from ones being stored
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
        abd_out['H'] = 1.0 - abd['H+'] - 2*abd['H2'] - 2*abd['H2+']\
                       - 3*abd['H3+'] - abd['OHx'] - abd['CHx'] \
                       - abd['HCO+']

        # He abundance = total He abundance - He+ abundance
        abd_out['He'] = self.xHe - abd['He+']

        # O abundance
        abd_out['O'] = self.xO - abd['O+'] - abd['CO'] - abd['HCO+'] \
                       - abd['OHx']

        # C abundance
        abd_out['C'] = self.xC - abd['C+'] - abd['CO'] - abd['HCO+'] \
                       - abd['CHx']

        # Neutral Si abundance = total Si abundance - ionized
        # Si abundance
        abd_out['Si'] = self.xSi - abd['Si+']

        # e abundance
        abd_out['e-'] = abd['H+'] + abd['H2+'] + abd['H3+'] + \
                        abd['He+'] + abd['C+'] + abd['O+'] + \
                        abd['HCO+'] + abd['Si+']

        # Return
        if outdict:
            return abd_out
        else:
            return abd_out.x
        

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
        # abundance of the most abundant hydrogenic species in the
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

                
    ####################################################################
    # Method to return the time derivative of all chemical rates
    ####################################################################
    def dxdt(self, xin, time):
        """
        This method returns the time derivative of all abundances for
        this chemical network.

        Parameters
           xin : array
              current abundances of all species
           time : float
              current time; not actually used, but included as an
              argument for compatibility with odeint

        Returns
           dxdt : array
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
                         [[self.NH*abd['C'], self.NH*abd['H2']],
                          [self.NH*abd['CO'], self.NH*abd['H2']],
                          [self.NH*abd['H2'], sigma_tot]])

        # Do H and H2 CR ionization here; we handle these separately
        # because the rates are non-trivial functions of the total H
        # and H2 abundances, so they can't be handled through the
        # standard CR reaction machinery
        Gamma_H = (2.3*abd['H2'] + 1.5*abd['H'])*self.ionRate
        xdot[0,1] += Gamma_H * abd['H']      # H+
        xdot[0,0] -= 2*Gamma_H * abd['H2']   # H2
        xdot[0,2] += 2*Gamma_H * abd['H2']   # H2+

        # Get rates for two-body reactions and add their contribution
        # to xdot
        ratecoef = _twobody_ratecoef(
            self.temp, n, abd['H']*n, abd['H2']*n,
            abd['e-']*n, np.exp(-1.87*self.AV)*self.chi)
        xdot += _twobody.dxdt(abd.x, n, ratecoef, self.Zd)

        # Return results, omitting derived species
        return np.ravel(xdot)[:len(specList)]

    ####################################################################
    # Method to return the time derivative of all chemical rates due
    # to only one type of reaction
    ####################################################################
    def dxdt_type(self, xin, rtype='2'):
        """
        This method returns the time derivative of all abundances for
        this chemical network.

        Parameters
           xin : array
              current abundances of all species
           rtype : '2' | 'CR' | 'ph'
              reaction type to return: two-body reactions, CR
              reactions, or photoionization reactions

        Returns
           dxdt : array
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

        # Get appropriate reaction type
        if rtype == 'CR' :
            xdot = _cr.dxdt(abd.x, n, self.ionRate)
            Gamma_H = (2.3*abd['H2'] + 1.5*abd['H'])*self.ionRate
            xdot[0,1] += Gamma_H * abd['H']      # H+
            xdot[0,0] -= 2*Gamma_H * abd['H2']   # H2
            xdot[0,2] += 2*Gamma_H * abd['H2']   # H2+
        elif rtype == 'ph':
            xdot = _ph.dxdt(abd.x, n, self.chi, self.AV, 
                            [[self.NH*abd['C'], self.NH*abd['H2']],
                             [self.NH*abd['CO'], self.NH*abd['H2']],
                             [self.NH*abd['H2'], sigma_tot]])
        elif rtype == '2':
            ratecoef = _twobody_ratecoef(
                self.temp, n, abd['H']*n, abd['H2']*n,
                abd['e-']*n, np.exp(-1.87*self.AV)*self.chi)
            print(ratecoef)
            xdot = _twobody.dxdt(abd.x, n, ratecoef, self.Zd)
        else:
            raise ValueError('unknown reaction type '+str(rtype))

        # Return results, omitting derived species
        return np.ravel(xdot)[:len(specList)]
