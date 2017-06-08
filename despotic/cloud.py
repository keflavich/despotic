"""
This module defines the class cloud, which is the basic class of
DESPOTIC. The cloud class defines the properties of a cloud, both
physical and chemical, and provides methods to perform a variety of
computations on that cloud.
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
from scipy.optimize import root
from scipy.optimize import brentq
from scipy.integrate import odeint
from .emitter import emitter
from .composition import composition
from .radiation import radiation
from .dustProp import dustProp
from .chemistry import abundanceDict, chemEvol, setChemEq
from .despoticError import despoticError
from copy import deepcopy

# Define some global physical constants in cgs units
import scipy.constants as physcons
kB = physcons.k*1e7
c = physcons.c*1e2
mH = physcons.m_p*1e3
h = physcons.h*1e7
sigma=physcons.sigma*1e3
a = 4*sigma/c
G = physcons.G*1e3
Ryd = h*c*physcons.Rydberg/100.
TLyA = 3*Ryd/(4*kB)
TLyB = 8*Ryd/(9*kB)

# Small numerical value
small = 1e-50

# Hardwired temperature limits for solvers
Tlo = 1.0
Thi = 1.0e4

class cloud(object):
    """
    A class describing the properties of an interstellar cloud, and
    providing methods to perform calculations using those properties.

    Parameters
       fileName : string
          name of file from which to read cloud description
       verbose : Boolean
          print out information about the cloud as we read it

    Class attributes
       nH : float
          number density of H nuclei, in cm^-3
       colDen : float
          center-to-edge column density of H nuclei, in cm^-2
       sigmaNT : float
          non-thermal velocity dispersion, in cm s^-1
       dVdr : float
          radial velocity gradient, in s^-1 (or cm s^-1 cm^-1)
       Tg : float
          gas kinetic temperature, in K
       Td : float
          dust temperature, in K
       comp : class composition
          a class that stores information about the chemical
          composition of the cloud
       dust : class dustProp
          a class that stores information about the properties of the
          dust in a cloud
       rad : class radiation
          the radiation field impinging on the cloud
       emitters : dict
          keys of the dict are names of emitting species, and values
          are objects of class emitter
       chemnetwork : chemical network class (optional)
          a chemical network that is to be used to perform
          time-dependent chemical evolution calcualtions for this cloud
    """


    ####################################################################
    # Method to initialize
    ####################################################################
    def __init__(self, fileName=None, verbose=False):
        """
        Parameters
           fileName : string
              name of file from which to read cloud description
           verbose : Boolean
              print out information about the cloud as we read it

        Returns
           Nothing
        """

        # Initial values when class is created
        self.nH = 0.0
        self.colDen = 0.0
        self.sigmaNT = 0.0
        self.dVdr = 0.0
        self.Tg = 0.0
        self.Td = 0.0
        self.comp = composition()
        self.rad = radiation()
        self.dust = dustProp()
        self.emitters = {}

        # Read if file name is given
        if fileName != None:
            self.read(fileName, verbose=verbose)


    ####################################################################
    # Method to return the list of emitters for a cloud
    ####################################################################
    def emitterList(self):
        return self.emitters.keys()

    ####################################################################
    # Method to get and set abundances in a cloud
    ####################################################################

    @property
    def abundances(self):
        """
        This property contains the abundances of all emitting species,
        stored as an abundanceDict
        """
        specList = self.emitters.keys()
        x = np.array([self.emitters[s].abundance for s in
                      specList])
        return abundanceDict(specList, x)

    @abundances.setter
    def abundances(self, abd):
        """
        Set abundances of all emitting species; if these species exist
        in the chemical network, they will be updated too

        Parameters
           abd : dict or abundanceDict
              dict containing species names as keys, and abundances of
              those species as values; species in abd that are not
              already in the emitter list will be added
        """
        for k in abd.keys():
            if k in self.emitters.keys():
                self.emitters[k].abundance = abd[k]
            else:
                self.addEmitter(k, abd[k])

    @property
    def chemabundances(self):
        """
        The property contains the abundances of all species in the
        chemical network, stored as an abundanceDict
        """
        if self.chemnetwork is None:
            raise despoticError(
                "cloud: cannot get chemabundances when chemnetwork is None")
        return self.chemnetwork.abundances

    @chemabundances.setter
    def chemabundances(self, abd):
        """
        Set abundances of all species in the chemical network, and
        update the emitter abundances to match

        Parameters
           abd : dict or abundanceDict
              dict containing species names as keys, and abundances of
              those species as values; species in abd that are not
              already in the emitter list will be added
        """
        if self.chemnetwork is None:
            raise despoticError(
                "cloud: cannot set chemabundances when chemnetwork is None")
        self.chemnetwork.abundances = abd
        self.chemnetwork.applyAbundances()

    ####################################################################
    # Method to set velocity dispersion to virial value
    ####################################################################
    def setVirial(self, alphaVir=1.0, setColDen=False, setnH=False, \
                      NTonly=False):
        """
        Set sigmaNT, colDen, or nH to the value required to give a
        virial ratio of unity

        Parameters
           alphaVir : float
              virial ratio to be used in computation; defaults to 1
           setColDen : Boolean
              if True, sigmaNT and nH are fixed, and colDen is
              altered to give the desired virial ratio
           setnH : Boolean
              if True, sigmaNT and colDen are fixed, and nH is altered
              to give the desired virial ratio
           NTonly : Boolean
              if True, the virial ratio is computed considering only the
              non-thermal component of the velocity dispersion

        Returns
           Nothing

        Remarks
           By default the routine fixes nH and colDen and computes
           sigmaNT, but this can be overridden by specifying either
           setColDen or setnH. It is an error to set both of these to
           True.
        """

        # Safety check
        if setColDen==True and setnH==True:
            raise despoticError(
                "setVirial: cannot use both setColDen and setnH")

        # Thermal velocity disperison squared
        if NTonly == False:
            sigmaThSqr = kB*self.Tg / (self.comp.mu*mH)
        else:
            sigmaThSqr = 0.0

        # Default case
        if setColDen==False and setnH==False:
            # Total velocity dispersion squared, from definition of
            # alphaVir
            sigmaTotSqr = (3.0*np.pi*alphaVir*G/20.0) * \
                self.colDen**2/self.nH * self.comp.muH*mH
            # Set non-thermal part
            if sigmaTotSqr > sigmaThSqr:
                self.sigmaNT = np.sqrt(sigmaTotSqr - sigmaThSqr)
            else:
                self.sigmaNT = 0.0
                print("setVirial warning: setting sigmaNT = 0, " +
                      "virial ratio still exceeds desired value")
        elif setnH==True:
            # Set nH from colDen and sigmaNT
            sigmaTotSqr = np.sqrt(self.sigmaNT**2 + sigmaThSqr)
            self.nH = (3.0*np.pi*alphaVir*G/20.0) * \
                self.colDen**2/sigmaTotSqr * self.comp.muH*mH
        else:
            # Set colDen from nH and sigmaNT
            sigmaTotSqr = np.sqrt(self.sigmaNT**2 + sigmaThSqr)
            self.colDen = np.sqrt( \
                sigmaTotSqr*self.nH / (3.0*np.pi*alphaVir*G/20.0) / \
                    (self.comp.muH*mH) )
                                    

    ####################################################################
    # Method to read cloud properties from a file
    ####################################################################
    def read(self, fileName, verbose=False):
        """
        Read the composition from a file

        Pamameters
           fileName : string
              string giving the name of the composition file
           verbose : Boolean
              print out information about the cloud as it is read
        
        Returns
           Nothing

        Remarks
           For the format of cloud files, see the documentation
        """

        # Read file
        try:
            fp = open(fileName, 'r')
            if verbose:
                print("Reading from file "+fileName+"...")
        except IOError:
            raise despoticError("cannot open file "+fileName)
        for line in fp:

            # Skip empty and comment lines
            if line=='\n':
                continue
            if line.strip()[0] == "#":
                continue

            # Break line up based on equal sign
            linesplit = line.split("=")
            if len(linesplit) < 2:
                raise despoticError("Error parsing input line: "+line)
            if linesplit[1] == '':
                raise despoticError("Error parsing input line: "+line)

            # Trim trailing comments from portion after equal sign
            linesplit2 = linesplit[1].split('#')

            # Proceed based on the token that precedes the equal sign
            if linesplit[0].upper().strip() == 'NH':

                self.nH = float(linesplit2[0])
                if verbose:
                    print("Setting nH = "+str(self.nH))

            elif linesplit[0].upper().strip() == 'COLDEN':

                self.colDen = float(linesplit2[0])
                if verbose:
                    print("Setting column density = " +
                          str(self.colDen) + " H cm^-2")

            elif linesplit[0].upper().strip() == 'SIGMANT':

                self.sigmaNT = float(linesplit2[0])
                if verbose:
                    print("Setting sigmaNT = " + 
                          str(self.sigmaNT) + " cm s^-1")

            elif linesplit[0].upper().strip() == 'DVDR':

                self.dVdr = float(linesplit2[0])
                if verbose:
                    print("Setting dVdr = " + 
                          str(self.dVdr) + " cm s^-1 cm^-1")

            elif linesplit[0].upper().strip() == 'TG':

                self.Tg = float(linesplit2[0])
                if verbose:
                    print("Setting Tg = "+str(self.Tg) + " K")

            elif linesplit[0].upper().strip() == 'TD':

                self.Td = float(linesplit2[0])
                if verbose:
                    print("Setting Td = "+str(self.Td) + " K")

            elif linesplit[0].upper().strip() == 'ALPHAGD':

                self.dust.alphaGD = float(linesplit2[0])
                if verbose:
                    print("Setting alpha_GD = " + 
                          str(self.dust.alphaGD) + 
                          " erg cm^3 K^-3/2")

            elif linesplit[0].upper().strip() == 'SIGMAD10':

                self.dust.sigma10 = float(linesplit2[0])
                if verbose:
                    print("Setting sigma_d,10 = " + 
                          str(self.dust.sigma10) + 
                          " cm^2 g^-1")

            elif linesplit[0].upper().strip() == 'SIGMADPE':

                self.dust.sigmaPE = float(linesplit2[0])
                if verbose:
                    print("Setting sigma_d,PE = " + 
                          str(self.dust.sigmaPE) + 
                          " cm^2 H^-1")

            elif linesplit[0].upper().strip() == 'SIGMADISRF':

                self.dust.sigmaISRF = float(linesplit2[0])
                if verbose:
                    print("Setting sigma_d,ISRF = " + 
                          str(self.dust.sigmaISRF) + 
                          " cm^2 H^-1")

            elif linesplit[0].upper().strip() == 'ZDUST':

                self.dust.Zd = float(linesplit2[0])
                if verbose:
                    print("Setting Z'_d = " + 
                          str(self.dust.Zd))

            elif linesplit[0].upper().strip() == 'BETADUST':

                self.dust.beta = float(linesplit2[0])
                if verbose:
                    print("Setting beta_dust = "+str(self.dust.beta))

            elif linesplit[0].upper().strip() == 'XHI':

                self.comp.xHI = float(linesplit2[0])
                if verbose:
                    print("Setting xHI = "+str(self.comp.xHI))

            elif linesplit[0].upper().strip() == 'XPH2':

                self.comp.xpH2 = float(linesplit2[0])
                if verbose:
                    print("Setting xpH2 = "+str(self.comp.xpH2))

            elif linesplit[0].upper().strip() == 'XOH2':

                self.comp.xoH2 = float(linesplit2[0])
                if verbose:
                    print("Setting xoH2 = "+str(self.comp.xoH2))

            elif linesplit[0].upper().strip() == 'H2OPR':

                self.comp.H2OPR = float(linesplit2[0])
                if verbose:
                    print("Setting H2 ortho-para ratio = "+
                          str(self.comp.H2OPR))

            elif linesplit[0].upper().strip() == 'XH2':

                self.comp.xH2 = float(linesplit2[0])
                if self.comp.H2OPR is None:
                    self.comp.H2OPR = 0.25
                    print("Warning: H2 OPR unspecified, assuming 0.25")
                if verbose:
                    print("Setting xpH2 = "+str(self.comp.xpH2))
                    print("Setting xoH2 = "+str(self.comp.xoH2))

            elif linesplit[0].upper().strip() == 'XHE':

                self.comp.xHe = float(linesplit2[0])
                if verbose:
                    print("Setting xHe = "+str(self.comp.xHe))

            elif linesplit[0].upper().strip() == 'XE':

                self.comp.xe = float(linesplit2[0])
                if verbose:
                    print("Setting xe = "+str(self.comp.xe))

            elif linesplit[0].upper().strip() == 'XH+':

                self.comp.xe = float(linesplit2[0])
                if verbose:
                    print("Setting xH+ = "+str(self.comp.xe))

            elif linesplit[0].upper().strip() == 'TCMB':

                self.rad.TCMB = float(linesplit2[0])
                if verbose:
                    print("Setting T_CMB = "+str(self.rad.TCMB)+" K")

            elif linesplit[0].upper().strip() == 'TRADDUST':

                self.rad.TradDust = float(linesplit2[0])
                if verbose:
                    print("Setting T_radDust = " +
                          str(self.rad.TradDust)+" K")

            elif linesplit[0].upper().strip() == 'RADDUTDILUTION':

                self.rad.fdDilute = float(linesplit2[0])
                if verbose:
                    print("Setting radDust dilution factor = " + 
                          str(self.rad.fdDilute))

            elif linesplit[0].upper().strip() == 'IONRATE':

                self.rad.ionRate = float(linesplit2[0])
                if verbose:
                    print("Setting primary ionization rate = " + 
                          str(self.rad.ionRate)+" s^-1 H^-1")

            elif linesplit[0].upper().strip() == 'CHI':

                self.rad.chi = float(linesplit2[0])
                if verbose:
                    print("Setting chi = " + 
                          str(self.rad.chi))

            elif linesplit[0].upper().strip() == 'EMITTER':

                # Emitter lines are complicated. There are two
                # required elements, a name and an abundance, that
                # must come first. There are also four optional
                # elements: energySkip, noExtrap, file:FileName,
                # and URL:url

                # Split up the tokens after the equal sign
                linesplit3 = linesplit2[0].split()

                # Make sure the number of tokens is acceptable
                if len(linesplit3) < 2 or len(linesplit3) > 6:
                    raise despoticError("Error parsing input line: "+line)

                # Do we have optional tokens?
                if len(linesplit3) == 2:

                    # Handle case of just two tokens
                    if verbose:
                        print("Adding emitter "+linesplit3[0]+ 
                              " with abundance "+linesplit3[1])
                    self.addEmitter(linesplit3[0], \
                                        float(linesplit3[1]))

                else:

                    # We have optional tokens; initialize various
                    # options to their defaults, then alter them based
                    # on the tokens we've been given
                    energySkip=False
                    extrap=True
                    emitterFile=None
                    emitterURL=None
                    for token in linesplit3[2:]:
                        if token.upper().strip() == 'ENERGYSKIP':
                            energySkip=True
                        elif token.upper().strip() == 'EXTRAPOLATE':
                            # Allowed to maintain backward compatibility
                            pass
                        elif token.upper().strip() == 'NOEXTRAP':
                            extrap=False
                        elif token.upper().strip()[0:5] == 'FILE:':
                            emitterFile=token[5:].strip()
                        elif token.upper().strip()[0:4] == 'URL:':
                            emitterURL=token[4:].strip()
                        else:
                            raise despoticError(
                                'unrecognized token "' +
                                token.strip()+'" in line: '
                                + line)

                    # Now print message and add emitter
                    if verbose:
                        msg = "Adding emitter "+linesplit3[0]+ \
                                  " with abundance "+linesplit3[1]
                        if energySkip:
                            msg += "; setting energySkip"
                        if not extrap:
                            msg += "; disallowing extrapolation"
                        if emitterFile != None:
                            msg += "; using file name "+emitterFile
                        if emitterURL != None:
                            msg += "; using URL "+emitterURL
                        print(msg)
                    self.addEmitter(linesplit3[0],
                                    float(linesplit3[1]),
                                    energySkip=energySkip,
                                    extrap=extrap,
                                    emitterFile=emitterFile,
                                    emitterURL=emitterURL)

            else:
                # Line does not correspond to any known keyword, so
                # throw an error
                raise despoticError("unrecognized token " +
                    linesplit[0].strip() + " in file " + fileName)

        # Close file
        fp.close()

        # Check that the hydrogen adds up. If not, raise error
        if self.comp.xHI + self.comp.xHplus + \
                2.0*(self.comp.xpH2 + self.comp.xoH2) != 1:
            raise despoticError(
                "total hydrogen abundance xHI + xH+ + 2 xH2 != 1")

        # Set derived properties based on composition, temperature
        self.comp.computeDerived(self.nH)
        if self.Tg > 0.0:
            self.comp.computeCv(self.Tg)

        # If verbose, print results for derived quantities
        if verbose:
            print("Derived quantities:")
            print("   ===> mean mass per particle = " + 
                  str(self.comp.mu) + " mH")
            print("   ===> mean mass per H = " + 
                  str(self.comp.muH) + " mH")
            print("   ===> energy added per ionization = " + 
                  str(self.comp.qIon/1.6e-12) + " eV")
            if self.Tg > 0.0:
                print("   ===> c_v/(k_B n_H mu_H) = " + 
                      str(self.comp.cv))


    ####################################################################
    # Method to add an emitter
    ####################################################################
    def addEmitter(self, emitName, emitAbundance, emitterFile=None,
                   emitterURL=None, energySkip=False, extrap=True):
        """
        Method to add an emitting species

        Pamameters
           emitName : string
              name of the emitting species
           emitAbundance : float
              abundance of the emitting species relative to H
           emitterFile : string
              name of LAMDA file containing data on this species; this
              option overrides the default
           emitterURL : string
              URL of LAMDA file containing data on this species; this
              option overrides the default
           energySkip : Boolean
              if set to True, this species is ignored when calculating
              line cooling rates
           extrap : Boolean
              If set to True, collision rate coefficients for this species
              will be extrapolated to temperatures outside the range given
              in the LAMDA table. If False, no extrapolation is perfomed,
              and providing temperatures outside the range in the table
              produces an error

        Returns
           Nothing
        """
        self.emitters[emitName] \
            = emitter(emitName, emitAbundance, 
                      emitterFile=emitterFile, 
                      emitterURL=emitterURL, 
                      extrap=extrap, energySkip=energySkip)


    ####################################################################
    # Method to calculate instantaneous values of all heating, cooling
    # terms
    ####################################################################
    def dEdt(self, c1Grav=0.0, thin=False, LTE=False, 
             fixedLevPop=False, noClump=False, 
             escapeProbGeom='sphere', PsiUser=None, 
             sumOnly=False, dustOnly=False, gasOnly=False, 
             dustCoolOnly=False, dampFactor=0.5, 
             verbose=False, overrideSkip=False):
        """
        Return instantaneous values of heating / cooling terms

        Parameters
           c1Grav : float
              if this is non-zero, the cloud is assumed to be
              collapsing, and energy is added at a rate
              Gamma_grav = c1 mu_H m_H cs^2 sqrt(4 pi G rho)
           thin : Boolean
              if set to True, cloud is assumed to be opticall thin
           LTE : Boolean
             if set to True, gas is assumed to be in LTE
           fixedLevPop : Boolean
             if set to True, level populations and escape
             probabilities are not recomputed, so the cooling rate is
             based on whatever values are stored
           escapeProbGeom : string, 'sphere' or 'LVG' or 'slab'
             specifies the geometry to be assumed in calculating
             escape probabilities
           noClump : Boolean
             if set to True, the clumping factor used in estimating
             rates for n^2 processes is set to unity
           dampFactor : float
              damping factor to use in level population calculations;
              see emitter.setLevPopEscapeProb
           PsiUser : callable
              A user-specified function to add additional heating /
              cooling terms to the calculation. The function takes the
              cloud object as an argument, and must return a two-element
              array Psi, where Psi[0] = gas heating / cooling rate,
              Psi[1] = dust heating / cooling rate. Positive values
              indicate heating, negative values cooling, and units are
              assumed to be erg s^-1 H^-1.
           sumOnly : Boolean
              if true, rates contains only four entries: dEdtGas and
              dEdtDust give the heating / cooling rates for the
              gas and dust summed over all terms, and maxAbsdEdtGas and
              maxAbsdEdtDust give the largest of the absolute values of
              any of the contributing terms for dust and gas
           gasOnly : Boolean
              if true, the terms GammaISRF, GammaDustLine, LambdaDust, \
              and PsiUserDust are omitted from rates. If both gasOnly
              and sumOnly are true, the dict contains only dEdtGas
           dustOnly : Boolean
              if true, the terms GammaPE, GammaCR, LambdaLine,
              GamamDLine, and PsiUserGas are omitted from rates. If both
              dustOnly and sumOnly are true, the dict contains only
              dEdtDust. Important caveat: the value of dEdtDust returned
              in this case will not exactly match that returned if
              dustOnly is false, because it will not contain the
              contribution from gas line cooling radiation that is
              absorbed by the dust
           dustCoolOnly : Boolean
              as dustOnly, but except that now only the terms
              LambdaDust, PsiGD, and PsiUserDust are computed
           overrideSkip : Boolean
              if True, energySkip directives are ignored, and cooling
              rates are calculated for all species

        Returns
           rates : dict
             A dict containing the values of the various heating and
             cooling rate terms; all quantities are in units of erg s^-1
             H^-1, and by convention positive = heating, negative =
             cooling; for dust-gas exchange, positive indicates heating
             of gas, cooling of dust

           Elements of the dict are as follows by default, but can be
           altered by the additional keywords listed below:

           GammaPE : float
              photoelectric heating rate
           GammaCR : float
              cosmic ray heating rate
           GammaGrav : float
              gravitational contraction heating rate
           LambdaLine : dict
              cooling rate from lines; dictionary keys correspond to
              species in the emitter list, values give line cooling
              rate for that species
           LambdaLyA : float
              cooling rate via Lyman alpha emission
           LambdaLyB : float
              cooling rate via Lyman beta emission
           Lambda2p : float
              cooling rate via 2 photon emission
           PsiGD : float
              dust-gas energy exchange rate
           GammaDustISRF : float
              dust heating rate due to the ISRF
           GammaDustCMB : float
              dust heating rate due to CMB
           GammaDustIR : float
              dust heating rate due to IR field
           GammaDustLine : float
              dust heating rate due to absorption of line radiation
           LambdaDust : float
              dust cooling rate due to thermal emission
           PsiUserGas : float
              gas heating / cooling rate from user-specified
              function; only included if PsiUser != None
           PsiUserDust : float
              gas heating / cooling rate from user-specified
              function; only included is PsiUser != None
        """

        # Make sure composition-derived quantities are initialized
        if self.comp.mu == 0.0:
            self.comp.computeDerived(self.nH)

        # Get clumping factor
        if noClump == False:
            cs2 = kB * self.Tg / (self.comp.mu * mH)
            cfac = np.sqrt(1.0 + 0.75*self.sigmaNT**2/cs2)
        else:
            cfac = 1.0

        # Gas terms
        if dustOnly == False and dustCoolOnly == False:

            # Photoelectric heating rate
            GammaPE = 4.0e-26 * self.rad.chi * self.dust.Zd * \
                np.exp(-self.colDen * self.dust.sigmaPE / 2.0)

            # Gravitational heating rate (Masunaga & Inutsuka 1998)
            GammaGrav = c1Grav * kB * self.Tg / (self.comp.mu * mH) * \
                (4 * np.pi * G * self.nH * self.comp.muH * mH) * \
                self.comp.muH * mH

            # Cosmic ray heating rate
            GammaCR = self.rad.ionRate * self.comp.qIon

            # Cooling rate via collisional excitation of neutral
            # hydrogen by electrons; these rate coefficients come from
            # Osterbrock & Ferland, table 3.16
            if self.Tg < 1e3:
                upsilon2s = 0.29
                upsilon2p = 0.51
                upsilon3s = 0.066
                upsilon3p = 0.12
                upsilon3d = 0.063
            elif self.Tg < 1.5e3:
                wgt = np.log(self.Tg/1e3)/np.log(1.5)
                upsilon2s = 0.29*(1.0-wgt) + 0.32*wgt
                upsilon2p = 0.51*(1.0-wgt) + 0.60*wgt
                upsilon3s = 0.066*(1.0-wgt) + 0.071*wgt
                upsilon3p = 0.12*(1.0-wgt) + 0.13*wgt
                upsilon3d = 0.063*(1.0-wgt) + 0.068*wgt
            elif self.Tg < 2e3:
                wgt = np.log(self.Tg/1.5e3)/np.log(2.0/1.5)
                upsilon2s = 0.32*(1.0-wgt) + 0.35*wgt
                upsilon2p = 0.60*(1.0-wgt) + 0.69*wgt
                upsilon3s = 0.071*(1.0-wgt) + 0.077*wgt
                upsilon3p = 0.13*(1.0-wgt) + 0.14*wgt
                upsilon3d = 0.068*(1.0-wgt) + 0.073*wgt
            else:
                upsilon2s = 0.35
                upsilon2p = 0.69
                upsilon3s = 0.077
                upsilon3p = 0.14
                upsilon3d = 0.073
            # Cooling rates; note the factor of 2 in the denominator
            # is the statistical weight of the 1s state of neutral
            # hydrogen
            fac = 8.629e-6/(2*np.sqrt(self.Tg))
            exfacLyA = np.exp(-TLyA/self.Tg)
            exfacLyB = np.exp(-TLyB/self.Tg)
            Lambda2p = fac * exfacLyA * upsilon2s * self.comp.xHI * \
                       self.comp.xe * self.nH * kB * TLyA
            LambdaLyA = fac * exfacLyA * upsilon2p * self.comp.xHI * \
                       self.comp.xe * self.nH * kB * TLyA
            LambdaLyB = fac * exfacLyB * \
                        (upsilon3s + upsilon3p + upsilon3d) * \
                        self.comp.xHI * self.comp.xe * self.nH * \
                        kB * TLyB

            # Line cooling rate, and heating rate of dust by lines
            LambdaLine = {}
            GammaDLine = 0.0

            # Iterate over emitting species
            for em in self.emitters.values():

                # Skip energetically unimportant emitters
                if em.energySkip and not overrideSkip:
                    continue

                # Calculate level populations and escape probabilities, using
                # specified assumptions
                if fixedLevPop == False:
                    if LTE==True:   # LTE
                        em.setLevPopLTE(self.Tg) 
                        if thin==False:  # LTE, not optically thin
                            em.setEscapeProb(self)
                    elif thin==True:
                        # Optically thin but not in LTE
                        em.setLevPop(self, thin=thin, 
                                     noClump=noClump)
                    else:
                        # Neither optically thin nor in LTE; note that
                        # we try multiple times with progressively
                        # smaller damping factors if need be
                        attempts = 0
                        dFac = dampFactor
                        while em.setLevPopEscapeProb( 
                                self, escapeProbGeom = escapeProbGeom, 
                                noClump = noClump, 
                                verbose = verbose, 
                                dampFactor = dFac) == False:
                            # If we're here, we failed to converge the
                            # level populations, so try again with a
                            # smaller damping factor; allow two
                            # retries before giving up
                            dFac = dFac / 2.0
                            attempts = attempts + 1
                            if attempts > 5:
                                raise despoticError("convergence " +
                                    "failed for "+em.name)
                            else:
                                print("Warning: convergence failed " + 
                                    "for "+em.name+", RETRYING " + 
                                    "with damping factor = " + str(dFac))
                else:
                    # Safety check
                    if em.levPopInitialized == False:
                        raise despoticError(
                            "for emitter " + em.name + ": " + 
                            "cannot use fixedLevPop in dEdt" + 
                            " if any level populations are uninitialized")
                    if em.escapeProbInitialized == False and \
                            thin == False:
                        raise despoticError(
                            "for emitter " + em.name + ": " +
                            "cannot use fixedLevPop in dEdt" +
                            " if any escape probabilities are uninitialized")

                # Calculate cooling rates per H for all lines
                lineLum = em.luminosityPerH(self.rad, thin=thin)

                # Calculate total luminosities
                LambdaLine[em.name] = sum(lineLum)

                # Calculate dust heating rate due to lines
                if gasOnly == False:
                    betaDLine = 1.0 / \
                        (1.0 + 0.375 * self.colDen*self.dust.sigma10* \
                             (em.data.radFreq/(10.0*kB/h))**self.dust.beta)
                    GammaDLine += sum((1.0 - betaDLine)*lineLum)

        # End gas terms

        # Dust terms
        if gasOnly == False:

            # Optically thin dust cooling rate
            LambdaDThin = self.dust.sigma10 * \
                (self.Td/10.)**self.dust.beta * \
                c * a * self.Td**4

            # Optically thick dust cooling rate
            LambdaDThick = c * a * self.Td**4 / self.colDen

            # Actual dust cooling rate
            LambdaD = min(LambdaDThin, LambdaDThick)

            if dustCoolOnly == False:

                # ISRF heating rate
                GammaISRF = 3.9e-24 * self.rad.chi * self.dust.Zd * \
                    np.exp(-self.dust.sigmaISRF * self.colDen/2.0)

                # CMB heating rate
                GammaDCMB = self.dust.sigma10 * \
                    (self.rad.TCMB/10.)**self.dust.beta * \
                    c * a * self.rad.TCMB**4

                # IR heating rate
                GammaDIR = self.dust.sigma10 * \
                    (self.rad.TradDust/10.)**self.dust.beta * \
                    c * a * self.rad.TradDust**4

                # Optically thick IR heating limit
                GammaDIRThick = c * a * self.rad.TradDust**4 / self.colDen

                # Actual heating rate
                GammaDIR = min(GammaDIR, GammaDIRThick)

        # End dust terms

        # Grain-gas energy exchange rate
        PsiGD = self.dust.alphaGD * cfac * self.nH * \
            np.sqrt(self.Tg) * (self.Td - self.Tg)

        # User terms
        if PsiUser != None:
            PsiUserVal = PsiUser(self)
        else:
            PsiUserVal = np.zeros(2)

        # Build dict of results
        rates = {}
        if sumOnly == False:
            if dustOnly == False and dustCoolOnly == False:
                rates['GammaPE'] = GammaPE
                rates['GammaGrav'] = GammaGrav
                rates['GammaCR'] = GammaCR
                rates['LambdaLine'] = LambdaLine
                rates['LambdaLyA'] = LambdaLyA
                rates['LambdaLyB'] = LambdaLyB
                rates['Lambda2p'] = Lambda2p
                if PsiUser != None:
                    rates['PsiUserGas'] = PsiUserVal[0]
            if gasOnly == False:
                if dustCoolOnly == False:
                    rates['GammaDustISRF'] = GammaISRF
                    rates['GammaDustCMB'] = GammaDCMB
                    rates['GammaDustIR'] = GammaDIR
                rates['LambdaDust'] = LambdaD
                if PsiUser != None:
                    rates['PsiUserDust'] = PsiUserVal[1]
            if gasOnly == False and dustOnly == False and \
                    dustCoolOnly == False:
                rates['GammaDustLine'] = GammaDLine
            rates['PsiGD'] = PsiGD
        else:
            if dustOnly == False and dustCoolOnly == False:
                rates['dEdtGas'] = GammaPE + GammaGrav + GammaCR - \
                    sum(LambdaLine.values()) - LambdaLyA \
                    - LambdaLyB - Lambda2p + PsiGD + PsiUserVal[0]
                rates['maxAbsdEdtGas'] = \
                    max(abs(GammaPE), abs(GammaGrav), abs(GammaCR),
                        abs(sum(LambdaLine.values())),
                        abs(LambdaLyA), abs(LambdaLyB), abs(Lambda2p),
                        abs(PsiGD), abs(PsiUserVal[0]))
            if gasOnly == False:
                rates['dEdtDust'] = - LambdaD - PsiGD + PsiUserVal[1]
                rates['maxAbsdEdtDust'] = \
                    max(abs(LambdaD), abs(PsiGD), abs(PsiUserVal[1]))
                if dustCoolOnly == False:
                    rates['dEdtDust'] += GammaISRF + GammaDCMB + \
                        GammaDIR
                    rates['maxAbsdEdtDust'] = \
                        max(rates['maxAbsdEdtDust'], abs(GammaISRF),
                            abs(GammaDCMB), abs(GammaDIR))
                    if dustOnly == False:
                        rates['dEdtDust'] += GammaDLine
                        rates['maxAbsdEdtDust'] = \
                            max(rates['maxAbsdEdtDust'], \
                                    abs(GammaDLine))

        # Return
        return rates


    ####################################################################
    # Method to calculate equilibrium dust temperature for fixed gas
    # temperature
    ####################################################################
    def setDustTempEq(self, PsiUser=None, Tdinit=None,
                      noLines=False, noClump=False,
                      verbose=False, dampFactor=0.5):
        """
        Set Td to equilibrium dust temperature at fixed Tg

        Parameters
           Tdinit : float
              initial guess for gas temperature
           PsiUser : callable
              A user-specified function to add additional heating /
              cooling terms to the calculation. The function takes the
              cloud object as an argument, and must return a two-element
              array Psi, where Psi[0] = gas heating / cooling rate,
              Psi[1] = dust heating / cooling rate. Positive values
              indicate heating, negative values cooling, and units are
              assumed to be erg s^-1 H^-1.
           noLines : Boolean
              If True, line heating of the dust is ignored. This can
              make the calculation significantly faster.
           noClump : Boolean
             if set to True, the clumping factor used in estimating
             rates for n^2 processes is set to unity
           dampFactor : float
             damping factor to use in level population calculations;
             see emitter.setLevPopEscapeProb
           verbose : Boolean
             if True, diagnostic information is printed

        Returns
           success : Boolean
              True if dust temperature calculation converged, False if
              not
        """

        # Make sure composition-derived quantities are initialized
        if self.comp.mu == 0.0:
            self.comp.computeDerived(self.nH)

        # Step 1: initialize dust temperature if initial value is given
        if Tdinit != None:
            self.Td = Tdinit

        # Step 2: compute heating and cooling rates initially; this is
        # useful so that we can compute the terms that don't depend on
        # the dust temperature just once and store them for the rest
        # of the calculation
        rates = self.dEdt(PsiUser=PsiUser, dustOnly=noLines,
                          noClump=noClump, verbose=verbose,
                          dampFactor=dampFactor)
        GammaSum = rates['GammaDustISRF'] + rates['GammaDustIR'] + \
            rates['GammaDustCMB']
        GammaSumMax = max(rates['GammaDustISRF'],
                          rates['GammaDustIR'],
                          rates['GammaDustCMB'])
        if noLines == False:
            GammaSum += rates['GammaDustLine']
            GammaSumMax = max(GammaSumMax, rates['GammaDustLine'])

        # Step 3: if we don't have an initial guess, make one by
        # assuming that dust-gas coupling is negligible, and that the
        # cloud is not optically thick, which are often the case
        if self.Td == 0.0:
            self.Td = (GammaSum /
                       (c*a*self.dust.sigma10*0.1**self.dust.beta)) ** \
                (1.0/(4.0+self.dust.beta))

        # Step 4: get initial scaling
        rates = self.dEdt(dustOnly=True, sumOnly=True,
                          PsiUser=PsiUser,
                          fixedLevPop=True,
                          verbose=verbose)
        lumScale = rates['maxAbsdEdtDust']

        # Step 5: solve for T_d using Brent's method
        try:
            self.Td = np.exp(
                brentq(_dustTempResid, np.log(Tlo), 
                       np.log(Thi), maxiter=200,
                       args=(self, PsiUser,
                             GammaSum, GammaSumMax,
                             lumScale, verbose)))
        except (despoticError, RuntimeError):
            pass

        # Check that we're really converged; if not, try again
        # starting from guessed starting position
        if abs(_dustTempResid(
                np.log(self.Td), self, PsiUser, GammaSum,
                GammaSumMax, lumScale, False)) > 1.0e-3:
            self.Td = \
                (GammaSum /
                 (c*a*self.dust.sigma10*0.1**self.dust.beta)) ** \
                (1.0/(4.0+self.dust.beta))
            try:
                self.Td = np.exp(
                    brentq(_dustTempResid, 0, np.log(1e5),
                           maxiter=200,
                           args=(self, PsiUser,
                                 GammaSum, GammaSumMax,
                                 lumScale, verbose)))
            except (despoticError, RuntimeError):
                # If we're here, we failed to converge
                return False

        # Test for success and return appropriate value
        if abs(_dustTempResid(
                np.log(self.Td), self, PsiUser, GammaSum,
                GammaSumMax, lumScale, False)) > 1.0e-3:
            return False
        else:
            return True


    ####################################################################
    # Method to calculate equilibrium gas temperature for dust gas
    # temperature
    ####################################################################
    def setGasTempEq(self, c1Grav=0.0, thin=False, noClump=False,
                     LTE=False, Tginit=None, fixedLevPop=False,
                     escapeProbGeom='sphere', PsiUser=None,
                     verbose=False):
        """
        Set Tg to equilibrium gas temperature at fixed Td

        Parameters
           c1Grav : float
              if this is non-zero, the cloud is assumed to be
              collapsing, and energy is added at a rate
              Gamma_grav = c1 mu_H m_H cs^2 sqrt(4 pi G rho)
           thin : Boolean
              if set to True, cloud is assumed to be opticall thin
           LTE : Boolean
              if set to True, gas is assumed to be in LTE
           Tginit : float
              initial guess for gas temperature
           fixedLevPop : Boolean
              if set to True, level populations are held fixed
              at the starting value, rather than caclculated
              self-consistently from the temperature
           escapeProbGeom : 'sphere' | 'LVG' | 'slab'
              specifies the geometry to be assumed in computing escape
              probabilities
           noClump : Boolean
              if set to True, the clumping factor used in estimating
              rates for n^2 processes is set to unity
           PsiUser : callable
              A user-specified function to add additional heating /
              cooling terms to the calculation. The function takes the
              cloud object as an argument, and must return a two-element
              array Psi, where Psi[0] = gas heating / cooling rate,
              Psi[1] = dust heating / cooling rate. Positive values
              indicate heating, negative values cooling, and units are
              assumed to be erg s^-1 H^-1.
           verbose : Boolean
              if True, print status messages while running

        Returns
           success : Boolean
              True if the calculation converges, False if it does not
        """

        # Make sure composition-derived quantities are initialized
        if self.comp.mu == 0.0:
            self.comp.computeDerived(self.nH)

        # Initialize gas temperatures if necessary
        if Tginit != None:
            self.Tg = Tginit
        else:
            if self.Tg==0.0:
                self.Tg = 10.0

        # Get initial scaling
        rates = self.dEdt(c1Grav=c1Grav,
                          thin=thin, LTE=LTE,
                          escapeProbGeom=escapeProbGeom,
                          gasOnly=True, noClump=noClump,
                          sumOnly=True, PsiUser=PsiUser)
        lumScale = rates['maxAbsdEdtGas']

        # Solve for Tg
        try:
            resid1 = _gasTempResid(
                np.log(Tlo), self, c1Grav, thin, LTE,
                escapeProbGeom, PsiUser, noClump,
                lumScale, verbose)
            Thi_in = Thi
            while True:
                resid2 = _gasTempResid(
                    np.log(Thi_in), self, c1Grav, thin, LTE,
                    escapeProbGeom, PsiUser, noClump,
                    lumScale, verbose)
                if resid1 * resid2 > 0:
                    Thi_in *= 10.0
                else:
                    break
            self.Tg = np.exp(
                brentq(_gasTempResid, np.log(Tlo), 
                       np.log(Thi_in),
                       args=(self, c1Grav, thin, LTE,
                             escapeProbGeom, PsiUser,
                             noClump, lumScale, verbose)))
        except (despoticError, RuntimeError):
            # If we're here, we failed to converge
            return False

        # Check for success and return appropriate value
        if abs(_gasTempResid(
                np.log(self.Tg), self, c1Grav, thin, LTE,
                escapeProbGeom, PsiUser,
                noClump, lumScale, verbose)) > 1.0e-3:
            return False
        else:
            return True


    ####################################################################
    # Method to calculate equilibrium gas and dust temperatures
    # simultaneously
    ####################################################################
    def setTempEq(self, c1Grav=0.0, thin=False, noClump=False,
                  LTE=False, Tinit=None, fixedLevPop=False,
                  escapeProbGeom='sphere', PsiUser=None,
                  verbose=False, tol=1e-4):
        """
        Set Tg and Td to equilibrium gas and dust temperatures

        Parameters
           c1Grav : float
              coefficient for rate of gas gravitational heating
           thin : Boolean
              if set to True, cloud is assumed to be opticall thin
           LTE : Boolean
              if set to True, gas is assumed to be in LTE
           Tinit : array(2)
              initial guess for gas and dust temperature
           noClump : Boolean
              if set to True, the clumping factor used in estimating
              rates for n^2 processes is set to unity
           fixedLevPop : Boolean
              if set to true, level populations are held fixed
              at the starting value, rather than caclculated
              self-consistently from the temperature
           escapeProbGeom : 'sphere' | 'LVG' | 'slab'
              specifies the geometry to be assumed in computing escape
              probabilities
           PsiUser : callable
              A user-specified function to add additional heating /
              cooling terms to the calculation. The function takes the
              cloud object as an argument, and must return a two-element
              array Psi, where Psi[0] = gas heating / cooling rate,
              Psi[1] = dust heating / cooling rate. Positive values
              indicate heating, negative values cooling, and units are
              assumed to be erg s^-1 H^-1.
           verbose : Boolean
              if True, the code prints diagnostic information as it runs

        Returns
           success : Boolean
              True if the iteration converges, False if it does not
        """

        # Make sure composition-derived quantities are initialized
        if self.comp.mu == 0.0:
            self.comp.computeDerived(self.nH)

        # As an initial guess, set the dust temperature to equilibrium
        # at fixed gas temperature, then set the gas temperature to
        # equilibrium at fixed dust temperature
        if Tinit is None:
            Tinit = [None, None]
        self.setDustTempEq(noClump=noClump, Tdinit=Tinit[1], 
                           PsiUser=PsiUser, verbose=verbose)
        self.setGasTempEq(
            c1Grav=c1Grav, thin=thin, noClump=noClump,
            LTE=LTE, Tginit=Tinit[0], fixedLevPop=fixedLevPop,
            escapeProbGeom=escapeProbGeom, PsiUser=PsiUser,
            verbose=verbose)
        Tinit = np.array([self.Tg, self.Td])

        # Get luminosity scaling
        rates = self.dEdt(c1Grav=c1Grav, thin=thin, LTE=LTE, 
                          escapeProbGeom=escapeProbGeom, 
                          sumOnly=True, PsiUser=PsiUser, 
                          noClump=noClump, 
                          verbose=verbose)
        lumScale = np.zeros(2)
        lumScale[0] = rates['maxAbsdEdtGas']
        lumScale[1] = rates['maxAbsdEdtDust']

        # Iterate to get equilibrium temperatures
        res = root(_gdTempResid, np.log(Tinit),
                   args=(self, c1Grav, thin,
                         LTE, escapeProbGeom, 
                         PsiUser, noClump, lumScale, 
                         verbose), 
                   method='hybr', options = { 'xtol' : tol })

        # If we failed to converge, make one more try from a
        # significantly larger starting temperature; larger starting
        # temperatures seem to provide easier convergence sometimes
        if not res.success:
            Tnew = Tinit * 5.0
            Tnew[Tnew < 50.] = 50.
            res = root(_gdTempResid, np.log(Tnew), 
                       args=(self, c1Grav, thin, 
                             LTE, escapeProbGeom, 
                             PsiUser, noClump, lumScale, 
                             verbose), 
                       method='hybr', options = { 'xtol' : tol })

        # Make sure we've converged.
        if res.success == False:
            return False

        # Store final result
        self.Tg = np.exp(res.x[0])
        self.Td = np.exp(res.x[1])
        return True


    ####################################################################
    # Method to calculate time-dependent evolution of gas temperature
    # for given starting conditions; note that we assume that the dust
    # is always in thermal equilibrium, since its equilibration time
    # is small compared to that of the gas.
    ####################################################################
    def tempEvol(self, tFin, tInit=0.0, c1Grav=0.0, noClump=False,
                 thin=False, LTE=False, fixedLevPop=False,
                 escapeProbGeom='sphere', nOut=100, dt=None,
                 tOut=None, PsiUser=None, isobaric=False,
                 fullOutput=False, dampFactor=0.5,
                 verbose=False):
        """
        Calculate the evolution of the gas temperature in time

        Parameters
           tFin : float
              end time of integration, in sec
           tInit : float
              start time of integration, in sec
           c1Grav : float
              coefficient for rate of gas gravitational heating
           thin : Boolean
              if set to True, cloud is assumed to be opticall thin
           LTE : Boolean
              if set to True, gas is assumed to be in LTE
           isobaric : Boolean
              if set to True, cooling is calculated isobarically;
              otherwise (default behavior) it is computed
              isochorically
           fixedLevPop : Boolean
              if set to true, level populations are held fixed
              at the starting value, rather than caclculated
              self-consistently from the temperature
           noClump : Boolean
              if set to True, the clumping factor used in estimating
              rates for n^2 processes is set to unity
           escapeProbGeom : string, 'sphere' or 'LVG' or 'slab'
              specifies the geometry to be assumed in computing escape
              probabilities
           nOut : int
              number of times at which to report the temperature; this
              is ignored if dt or tOut are set
           dt : float
              time interval between outputs, in s; this is ignored if
              tOut is set
           tOut : array
              list of times at which to output the temperature, in s;
              must be sorted in increasing order
           PsiUser : callable
              A user-specified function to add additional heating /
              cooling terms to the calculation. The function takes the
              cloud object as an argument, and must return a two-element
              array Psi, where Psi[0] = gas heating / cooling rate,
              Psi[1] = dust heating / cooling rate. Positive values
              indicate heating, negative values cooling, and units are
              assumed to be erg s^-1 H^-1.
           fullOutput : Boolean
              if True, the full cloud state is returned at every time,
              as opposed to simply the gas temperature
           dampFactor : float
              damping factor to use in calculating level populations
              (see emitter for details)


        Returns
           if fullOutput == False:

           Tg : array
              array of gas temperatures at specified times, in K
           time : array
              array of output times, in sec

           if fullOutput == True:

           cloudState : list
              each element of the list is a deepcopy of the cloud
              state at the corresponding time; there is one list
              element per output time
           time : array of floats
              array of output times, in sec

        Remarks
           If the settings for nOut, dt, or tOut are such that some of
           the output times requested are past tEvol, the cloud will only
           be evolved up to time tEvol. Similarly, if the last output
           time is less than tEvol, the cloud will still be evolved up to
           time tEvol, and the gas temperature Tg will be set to its
           value at that time.   
        """

        # Make sure composition-derived quantities are initialized
        if self.comp.mu == 0.0:
            self.comp.computeDerived(self.nH)

        # Set up array of output times if necessary
        if tOut==None:
            if dt==None:
                tOut = tInit + np.arange(nOut+1)*float(tFin-tInit)/nOut
            else:
                tOut = np.arange(tInit, (tFin-tInit)*(1+1e-10), dt)

        # Sanity check on output times: eliminate any output times
        # that are not between tInit and tFin
        tOut1 = tOut[tOut >= tInit]
        tOut1 = tOut1[tOut1 <= tFin]

        # If we're isobaric, record the isobar we're on; otherwise set
        # the isobar value to -1 to flag that we're isochoric
        if isobaric:
            isobar = self.nH * self.Tg
        else:
            isobar = -1

        # Integrate the ODE to the requested times; if fullOuptut is
        # set, we need to manually stop the integration at the
        # requested times so that we can dump detailed output.
        if fullOutput == False:
            Tgout = \
                odeint(_gasTempDeriv, self.Tg, tOut1, \
                           args=(self, c1Grav, thin, LTE, \
                                     escapeProbGeom, PsiUser, \
                                     noClump, isobar, dampFactor, \
                                     verbose))
        else:
            cloudList = [deepcopy(cloud)]
            for i, t in enumerate(tOut1[1:]):
                tvec = [tOut[i], t]
                Tgtemp = \
                    odeint(_gasTempDeriv, self.Tg, tvec, \
                               args=(self, c1Grav, thin, LTE, \
                                         escapeProbGeom, PsiUser, \
                                         noClump, isobar, \
                                         dampFactor, verbose))
                cloudList.append(deepcopy(cloud))

        # If necessary, continue integrating up to tEvol
        if tOut1[-1] < tFin:
            odeint(_gasTempDeriv, self.Tg, \
                                 np.array([tFin-tOut1[-1]]), \
                                 args=(self, c1Grav, thin, LTE, \
                                           escapeProbGeom, PsiUser, \
                                           noClump, isobar, \
                                           dampFactor, verbose))


        # Return gas temperature history
        if fullOutput == False:
            return Tgout, tOut
        else:
            return cloudList, tOut


    ####################################################################
    # Method to return the continuum-subtrated luminosity / intensity /
    # brightness temperature of lines from the specified emitter
    ####################################################################
    def lineLum(self, emitName, LTE=False, noClump=False,
                transition=None, thin=False, intOnly=False,
                TBOnly=False, lumOnly=False,
                escapeProbGeom='sphere', dampFactor=0.5,
                noRecompute=False, abstol=1.0e-8,
                verbose=False):
        """
        Return the frequency-integrated intensity of various lines

        Parameters
           emitName : string
              name of the emitter for which the calculation is to be
              performed
           LTE : Boolean
              if True, and level populations are unitialized, they will
              be initialized to their LTE values; if they are
              initialized, this option is ignored
           noClump : Boolean
              if set to True, the clumping factor used in estimating
              rates for n^2 processes is set to unity
           transition : list of two arrays
              if left as None, luminosity is computed for all
              transitions; otherwise only selected transitions are
              computed, with transition[0] = array of upper states
              transition[1] = array of lower states
           thin : Boolean
              if True, the calculation is done assuming the cloud is
              optically thin; if level populations are uninitialized,
              and LTE is not set, they will be computed assuming the
              cloud is optically thin
           intOnly : Boolean
              if true, the output is simply an array containing the
              frequency-integrated intensity of the specified lines;
              mutually exclusive with TBOnly and lumOnly
           TBOnly : Boolean
              if True, the output is simply an array containing the
              velocity-integrated brightness temperatures of the
              specified lines; mutually exclusive with intOnly and
              lumOnly
           lumOnly : Boolean
              if True, the output is simply an array containing the
              luminosity per H nucleus in each of the specified lines;
              mutually eclusive with intOnly and TBOonly
           escapeProbGeom : 'sphere' | 'LVG' | 'slab'
              sets problem geometry that will be assumed in calculating
              escape probabilities; ignored if the escape probabilities
              are already initialized
           dampFactor : float
              damping factor to use in level population calculations;
              see emitter.setLevPopEscapeProb
           noRecompute : False
              if True, level populations and escape probabilities are
              not recomputed; instead, stored values are used

        Returns
           res : list or array

           if intOnly, TBOnly, and lumOnly are all False, each element
           of the list is a dict containing the following fields:

           'freq' : float
              frequency of the line in Hz
           'upper' : int
              index of upper state, with ground state = 0 and states
              ordered by energy
           'lower' : int
              index of lower state
           'Tupper' : float
              energy of the upper state in K (i.e. energy over kB)
           'Tex' : float
              excitation temperature relating the upper and lower levels
           'intIntensity' : float
              frequency-integrated intensity of the line, with the CMB
              contribution subtracted off; units are erg cm^-2 s^-1 sr^-1 
           'intTB' : float
              velocity-integrated brightness temperature of the line,
              with the CMB contribution subtracted off; units are K km
              s^-1
           'lumPerH' : float
              luminosity of the line per H nucleus; units are erg s^-1
              H^-1
           'tau' : float
              optical depth in the line, not including dust
           'tauDust' : float
              dust optical depth in the line

        if intOnly, TBOnly, or lumOnly are True: res is an array
        containing the intIntensity, TB, or lumPerH fields of the dict
        described above
        """

        # Make sure composition-derived quantities are initialized
        if self.comp.mu == 0.0:
            self.comp.computeDerived(self.nH)

        # Step 1. Safety check and initial setup
        if not emitName in self.emitters:
            raise despoticError('unknown emitter '+emitName)
        em=self.emitters[emitName]
        if transition==None:
            u = em.data.radUpper
            l = em.data.radLower
        else:
            u = transition[0]
            l = transition[1]

        # Step 2. Unless we've been asked not to, compute level
        # populations and escape probabilities
        if not noRecompute:
            if LTE==True:
                em.setLevPopLTE(self.Tg)
                if thin==False:
                   em.setEscapeProb(self, escapeProbGeom=escapeProbGeom) 
            elif thin==True:
                em.setLevPop(self, thin=True, noClump=noClump)
            else:
                em.setLevPopEscapeProb(
                    self, noClump=noClump, dampFactor=dampFactor, 
                    escapeProbGeom=escapeProbGeom, 
                    abstol=abstol, verbose=verbose)
        # Safety check: make sure we're initialized
        if em.levPopInitialized == False:
            raise despoticError('cannot use noRecompute if level' + 
                ' popuplations are uninitialized')
        if em.escapeProbInitialized == False and thin == False:
            raise despoticError('cannot use noRecompute if escape' + 
                ' probabilities are uninitialized')

        # Step 3. Compute luminosity per H
        lumPerH = self.emitters[emitName]. \
            luminosityPerH(self.rad, transition=transition,
                           thin=thin)
        if lumOnly == True:
            return lumPerH

        # Step 4. Compute frequency-integrated intensity (units of erg
        # cm^-2 s^-1 sr^-1), including dust attenuation
        tauDust = self.colDen*self.dust.sigma10 * \
            (em.data.freq[u,l]/(10.0*kB/h))**self.dust.beta
        intIntensity = lumPerH * self.colDen / (4.0*np.pi) / \
            (1.0 + 0.375*tauDust)
        if intOnly == True:
            return intIntensity

        # Step 5. Convert to velocity-integrated brightness
        # temperature; note the division by 10^5 to convert from cm
        # s^-1 to km s^-1; also not the special handling of negative
        # integrated intensities, corresponding to lines where there
        # is absorption of the background CMB. By convention we assign
        # these negative brightness temperatures, with a magnitude
        # equal to the brightness temperature of the absolute value of
        # the intensity
        intIntensityMask = intIntensity.copy()
        intIntensityMask[intIntensity <= 0] = small
        TB = h*self.emitters[emitName].data.freq[u,l]/kB / \
            np.log(1+2*h*em.data.freq[u,l]**3 / (c**2*intIntensityMask)) * \
                    c / (em.data.freq[u,l]) \
                    / 1e5
        TB[intIntensity == 0] = 0.0
        mask = np.where(intIntensity < 0)
        TB[mask] = \
            -h*self.emitters[emitName].data.freq[u[mask],l[mask]]/kB / \
            np.log(1+2*h*em.data.freq[u[mask],l[mask]]**3 / \
                    (-c**2*intIntensity[mask])) * \
                    c / (em.data.freq[u[mask],l[mask]]) \
                    / 1e5
        if TBOnly == True:
            return TB

        # Step 6. Construct output dict
        outDict = []
        if thin:
            tau = lumPerH * 0.0
        else:
            tau = em.opticalDepth(transition=transition, \
                                  escapeProbGeom=escapeProbGeom)
        for i, T in enumerate(TB):
            line = { \
                'freq' : em.data.freq[u[i], l[i]], \
                    'upper' : u[i], \
                    'lower' : l[i], \
                    'Tupper' : em.data.levTemp[u[i]], \
                    'Tex' : em.data.dT[u[i],l[i]] / \
                    np.log( em.data.levWgt[u[i]]*em.levPop[l[i]] / \
                             (em.data.levWgt[l[i]]*em.levPop[u[i]]) ), \
                    'lumPerH' : lumPerH[i], \
                    'intIntensity' : intIntensity[i], \
                    'intTB' : TB[i], \
                    'tau' : tau[i], \
                    'tauDust' : tauDust[i] }
            outDict.append(line)
        return outDict


    #####################################################################
    # Method to calculate time-dependent evolution of chemical
    # abundances; note that this is just a wrapper routine for
    # despotic.chemistry.setChemEq.setChemEq
    #####################################################################
    def setChemEq(self, tEqGuess=None, network=None, info=None,
                  addEmitters=False, tol=1e-6, maxTime=1e16,
                  verbose=False, smallabd=1e-15, convList=None, 
                  evolveTemp='fixed', isobaric=False,
                  tempEqParam=None, dEdtParam=None, maxTempIter=50):
        """
        Set the chemical abundances for a cloud to their equilibrium
        values, computed using a specified chemical network.

        Parameters
           tEqGuess : float
              a guess at the timescale over which equilibrium will be
              achieved; if left unspecified, the code will attempt to
              estimate this time scale on its own
           network : chemNetwork object
              the chemNetwork object to use; if None, the existing
              chemnetwork member of the class (if it exists) is used
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
              evolution; 'fixed' = treat tempeature as fixed; 'iterate' =
              iterate between setting the gas temperature and chemistry to
              equilibrium; 'iterateDust' = iterate between setting the gas
              and dust temperatures and the chemistry to equilibrium;
              'gasEq' = hold dust temperature fixed, set gas temperature to
              instantaneous equilibrium value as the chemistry evolves;
              'fullEq' = set gas and dust temperatures to instantaneous
              equilibrium values while evolving the chemistry network;
              'evol' = evolve gas temperature in time along with the
              chemistry, assuming the dust is always in instantaneous
              equilibrium
           isobaric : Boolean
              if set to True, the gas is assumed to be isobaric during
              the evolution (constant pressure); otherwise it is assumed
              to be isochoric; note that (since chemistry networks at
              present are not allowed to change the mean molecular
              weight), this option has no effect if evolveTemp is 'fixed'
           tempEqParam : None | dict
              if this is not None, then it must be a dict of values that
              will be passed as keyword arguments to the cloud.setTempEq,
              cloud.setGasTempEq, or cloud.setDustTempEq routines; only
              used if evolveTemp is not 'fixed'
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
           despoticError, if network is None and the cloud does not
           already have a defined chemical network associated with it

        Remarks
           The final abundances are written to the cloud whether or
           not the calculation converges.
        """

        return setChemEq(self, tEqGuess=tEqGuess, network=network,
                         info=info, tol=tol, maxTime=maxTime,
                         verbose=verbose, smallabd=smallabd, 
                         convList=convList, addEmitters=addEmitters,
                         evolveTemp=evolveTemp,
                         isobaric=isobaric, tempEqParam=tempEqParam,
                         dEdtParam=dEdtParam,
                         maxTempIter=maxTempIter)

########################################################################
# Method to calculate time-dependent evolution of chemical abundances;
# note that this is just a wrapper routine for
# despotic.chemistry.chemEvol.chemEvol
########################################################################
    def chemEvol(self, tFin, tInit=0.0, nOut=100, dt=None,
                 tOut=None, network=None, info=None,
                 addEmitters=False, evolveTemp='fixed',
                 isobaric=False, tempEqParam=None,
                 dEdtParam=None):
        """
        Evolve the chemical abundances of this cloud in time.

        Parameters
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
               list of times at which to output the temperature, in sec;
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
           evolveTemp : 'fixed' | 'gasEq' | 'fullEq' | 'evol'
              how to treat the temperature evolution during the
              chemical evolution; 'fixed' = treat tempeature as fixed;
              'gasEq' = hold dust temperature fixed, set gas temperature
              to instantaneous equilibrium value; 'fullEq' = set gas and
              dust temperatures to instantaneous equilibrium values;
              'evol' = evolve gas temperature in time along with the
              chemistry, assuming the dust is always in instantaneous
              equilibrium
           isobaric : Boolean
              if set to True, the gas is assumed to be isobaric during
              the evolution (constant pressure); otherwise it is assumed
              to be isochoric; note that (since chemistry networks at
              present are not allowed to change the mean molecular
              weight), this option has no effect if evolveTemp is 'fixed'
           tempEqParam : None | dict
              if this is not None, then it must be a dict of values that
              will be passed as keyword arguments to the cloud.setTempEq,
              cloud.setGasTempEq, or cloud.setDustTempEq routines; only
              used if evolveTemp is not 'fixed'
           dEdtParam : None | dict
              if this is not None, then it must be a dict of values that
              will be passed as keyword arguments to the cloud.dEdt
              routine; only used if evolveTemp is 'evol'

        Returns
           time : array of floats
               array of output times, in sec
           abundances : class abundanceDict
               an abundanceDict giving the abundances as a function of time
           Tg : array
              gas temperature as a function of time; returned only if
              evolveTemp is not 'fixed'
           Td : array
              dust temperature as a function of time; returned only if
              evolveTemp is not 'fixed' or 'gasEq'

        Raises
           despoticError, if network is None and the cloud does not already
           have a defined chemical network associated with it
        """

        return chemEvol(self, tFin, tInit=tInit, nOut=nOut,
                        dt=dt, tOut=tOut, network=network, info=info,
                        addEmitters=addEmitters,
                        evolveTemp=evolveTemp,
                        isobaric=isobaric, tempEqParam=tempEqParam,
                        dEdtParam=dEdtParam)


########################################################################
# End of class gasProp
########################################################################

def turb_heating_generator(lengthscale=3e18, turbulence=True):
    """
    Generator for a turbulent heating function with a given length scale
    (cloud objects do not generally contain size scale information)

    Parameters
    ----------
    lengthscale : float
        The length scale in centimeters
    turbulence : bool
        A boolean flag to turn the turbulent heating on or off.  If False,
        the heating will be zero

    Examples
    --------
    >>> turb_heating = turb_heating_generator(3e18)
    >>> gmc.setTempEq(escapeProbGeom='LVG', PsiUser=turb_heating)

    """
    def turb_heating(cloud, lengthscale=lengthscale):
        """ Turbulent heating rate depends on cloud linewidth
        (sigma_nonthermal) and driving scale of the turbulence
        DESPOTIC wants units of erg/s/H (per hydrogen), so the turbulent
        heating rate n sigma^3 / L is divided by n to get just sigma^3/L

        MacLow 1999, 2002, 2004 gives exactly:
            3e-27 erg cm^-3 s^-1 (n/1 cm^3) (v/10 km/s)**3 (L/100 pc)**-1

        Parameters
        ----------
        cloud : cloud
            A cloud object.  Must have sigmaNT defined, where sigmaNT is the 1D
            nonthermal velocity dispersion
        lengthscale : float
            The length scale in cm.  Most simply interpreted as the driving
            length scale
        """
        if turbulence:
            # sigmaNT is assumed to be the 1D velocity dispersion, so it is
            # multiplied by sqrt(3)
            # it has units cm/s
            gamturb = (1.4 * mH *
                       (0.5*3**1.5 * (cloud.sigmaNT)**3 / (lengthscale)))
            return [gamturb, 0]
        else:
            return [0,0]
    return turb_heating


try:
    from PyQt4 import QtGui, QtCore
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
    import numpy
    import sys

    class cloud_gui(QtGui.QWidget):

        def __init__(self, parent=None):
            '''A class which implements an interactive GUI to run the 'cloud' class.'''

            super(cloud_gui, self).__init__(parent)

            self.initUI()

        def initUI(self):
            '''Initialize the user interface'''

            import matplotlib
            matplotlib.use('Qt4Agg')
            import pylab

            #instantiating the cloud object
            self.cloud = cloud()

            # set up the plotting figure and axes (axs[0,0] is the one in the top left corner)
            self.figure, self.axs = pylab.subplots(4, 1, sharex = False, sharey = False, figsize=(8,8) )
            pylab.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95, wspace=0.0, hspace=0.0)            

            # this is the Canvas Widget that displays the `figure`
            # it takes the `figure` instance as a parameter to __init__
            self.canvas = FigureCanvas(self.figure)

            # this is the Navigation widget
            # it takes the Canvas widget and a parent
            self.toolbar = NavigationToolbar(self.canvas, self)

            #---------laying out the line QT edit feilds----------------------       
            self.lbl_emitter = QtGui.QLabel('emitter')
            self.qle_emitter = QtGui.QLineEdit('CO', self)

            self.lbl_x_emitter = QtGui.QLabel('x_emitter')
            self.qle_x_emitter = QtGui.QLineEdit('1e-4', self)
            
            self.lbl_nH = QtGui.QLabel('nH') 
            self.qle_nH = QtGui.QLineEdit('1e3', self)
            
            self.lbl_colDen = QtGui.QLabel('colDen')
            self.qle_colDen = QtGui.QLineEdit('1e22', self)        

            self.lbl_sigmaNT = QtGui.QLabel('sigmaNT')
            self.qle_sigmaNT = QtGui.QLineEdit('20e3', self)        

            self.lbl_Tg = QtGui.QLabel('Tg')
            self.qle_Tg = QtGui.QLineEdit('30.0', self)        

            self.lbl_xoH2 = QtGui.QLabel('xoH2')
            self.qle_xoH2 = QtGui.QLineEdit('0.1', self)        

            self.lbl_xpH2 = QtGui.QLabel('xpH2')
            self.qle_xpH2 = QtGui.QLineEdit('0.4', self)        

            self.lbl_xHe = QtGui.QLabel('xHe')
            self.qle_xHe = QtGui.QLineEdit('0.1', self)        
            #-------------------------------------------------------------------
            
            # The button connected to `plot` method
            self.button = QtGui.QPushButton('Plot')
            self.button.clicked.connect(self.plot)
            
            grid = QtGui.QGridLayout()
            grid.setSpacing(10)
            
            #location 1,0 on the grid spanning 1 row and 3 columns
            grid.addWidget(self.toolbar, 1, 0, 1, 3) 

            #adding the parameter widgets to the UI
            grid.addWidget(self.lbl_emitter, 1, 4)
            grid.addWidget(self.qle_emitter, 1, 5, 1, 2)

            grid.addWidget(self.lbl_x_emitter, 2, 4)
            grid.addWidget(self.qle_x_emitter, 2, 5, 1, 2)

            grid.addWidget(self.lbl_nH, 3, 4)
            grid.addWidget(self.qle_nH, 3, 5, 1, 2)

            grid.addWidget(self.lbl_colDen, 4, 4)
            grid.addWidget(self.qle_colDen, 4, 5, 1, 2)
            
            grid.addWidget(self.lbl_sigmaNT, 5, 4)
            grid.addWidget(self.qle_sigmaNT, 5, 5, 1, 2)

            grid.addWidget(self.lbl_Tg, 6, 4)
            grid.addWidget(self.qle_Tg, 6, 5, 1, 2)

            grid.addWidget(self.lbl_xoH2, 7, 4)
            grid.addWidget(self.qle_xoH2, 7, 5, 1, 2)

            grid.addWidget(self.lbl_xpH2, 8, 4)
            grid.addWidget(self.qle_xpH2, 8, 5, 1, 2)

            grid.addWidget(self.lbl_xHe, 9, 4)
            grid.addWidget(self.qle_xHe, 9, 5, 1, 2)

            #adding the plot button to the widget        
            grid.addWidget(self.button , 15, 5)
            
            grid.addWidget(self.canvas , 2, 0, 15, 3)

            self.setLayout(grid) 
            
            #setting the default ranges and plot data and labels and removing all xticks
            for ax in self.axs:
                ax.set_xticklabels([])
                ax.set_xlim(0, 20)
            
            #setting the label of the bottom x-axis and the ticklabels
            self.axs[-1].set_xlabel(r'J$_{\rm upper}$')
            self.axs[-1].set_xticklabels(numpy.int32(numpy.linspace(0, 20, 5)))
            
            #setting the default limits of the y-axes
            pylab.setp(self.axs[0], 'yscale', 'log'   , 'ylim', [1e-8 , 1e8], 'ylabel', r'T$_B$') 
            pylab.setp(self.axs[1], 'yscale', 'log'   , 'ylim', [1e-15, 1e3], 'ylabel', r'I$_{\rm int}$') 
            pylab.setp(self.axs[2], 'yscale', 'linear', 'ylim', [0.0  , 1e3], 'ylabel', r'$\tau$') 
            pylab.setp(self.axs[3], 'yscale', 'linear', 'ylim', [0.0  , 1e3], 'ylabel', r'$\tau_{\rm dust}$')

            #making the y-ticks more organized
            self.axs[0].set_yticks(self.axs[0].get_yticks()[2:-2])
            self.axs[1].set_yticks(self.axs[1].get_yticks()[2:-2])
            self.axs[2].set_yticks(self.axs[2].get_yticks()[1:-1])
            self.axs[3].set_yticks(self.axs[3].get_yticks()[1:-1])
            
            #making dummpy plots to be used later for faster plotting
            self.plt0_em,  = self.axs[0].plot([], [], 'ro')
            self.plt0_ab,  = self.axs[0].plot([], [], 'bo')
            self.plt0_all, = self.axs[0].plot([], [], 'k--')
            
            self.plt1_em,  = self.axs[1].plot([], [], 'ro')
            self.plt1_ab,  = self.axs[1].plot([], [], 'bo')
            self.plt1_all, = self.axs[1].plot([], [], 'k--')
            
            self.plt2_em,  = self.axs[2].plot([], [], 'ro')
            self.plt2_ab,  = self.axs[2].plot([], [], 'bo')
            self.plt2_all, = self.axs[2].plot([], [], 'k--')        
            
            self.plt3_em,  = self.axs[3].plot([], [], 'ro')
            self.plt3_ab,  = self.axs[3].plot([], [], 'bo')
            self.plt3_all, = self.axs[3].plot([], [], 'k--')
            
            self.setGeometry(300, 300, 1000, 700)
            self.setWindowTitle('Despotic-cloud')    
            self.show()
            
        def plot(self):

            #setting the parameters from the gui to the cloud object
            self.cloud.nH        = numpy.float(getattr(self, 'qle_' + 'nH').text())
            self.cloud.colDen    = numpy.float(getattr(self, 'qle_' + 'colDen').text())
            self.cloud.sigmaNT   = numpy.float(getattr(self, 'qle_' + 'sigmaNT').text())
            self.cloud.Tg        = numpy.float(getattr(self, 'qle_' + 'Tg').text())
            self.cloud.comp.xoH2 = numpy.float(getattr(self, 'qle_' + 'xoH2').text())
            self.cloud.comp.xpH2 = numpy.float(getattr(self, 'qle_' + 'xpH2').text())
            self.cloud.comp.xHe  = numpy.float(getattr(self, 'qle_' + 'xHe').text())

            emitter_str = numpy.string_(getattr(self, 'qle_' + 'emitter').text())         
            emitter_x   = numpy.float(getattr(self, 'qle_' + 'x_emitter').text())

            self.cloud.addEmitter(emitter_str, emitter_x)  
            
            #solving for the emission
            lines = self.cloud.lineLum(emitter_str)

            #plotting the info of the emitter
            
            #the upper level of the transition
            u = numpy.array([l['upper'] for l in lines])
            
            #brightness temperature (in K.km.s-1)
            intTB = numpy.array([l['intTB'] for l in lines])

            self.plt0_all.set_xdata(u)
            self.plt0_all.set_ydata(numpy.fabs(intTB))
            self.axs[0].set_ylim([numpy.fabs(intTB).min()/10.0, numpy.fabs(intTB).max()*10.0])
            
            ind_em = numpy.where(intTB > 0.0)[0]
            ind_ab = numpy.where(intTB < 0.0)[0]
            if ind_em.size > 0:                
                self.plt0_em.set_xdata(u[ind_em])
                self.plt0_em.set_ydata(intTB[ind_em])
            else:
                self.plt0_em.set_xdata([])
                self.plt0_em.set_ydata([])
                
            if ind_ab.size > 0:
                self.plt0_ab.set_xdata(u[ind_ab])
                self.plt0_ab.set_ydata(numpy.fabs(intTB[ind_ab]))
            else:
                self.plt0_ab.set_xdata([])
                self.plt0_ab.set_ydata([])

            #intensity after subtracting the CMB contributuin (in erg.cm-2.s-1)
            intIntensity = numpy.array([l['intIntensity'] for l in lines])
            mn, mx = [max(numpy.fabs(intIntensity).min()/10.0, 1e-15), numpy.fabs(intIntensity).max()*10.0] 
            self.axs[1].set_ylim([mn, mx])
     
            self.plt1_all.set_xdata(u)
            self.plt1_all.set_ydata(numpy.fabs(intIntensity))
            ind_em = numpy.where(intIntensity > 0.0)[0]
            ind_ab = numpy.where(intIntensity < 0.0)[0]
            if ind_em.size > 0:                
                self.plt1_em.set_xdata(u[ind_em])
                self.plt1_em.set_ydata(intIntensity[ind_em])
            else:
                self.plt1_em.set_xdata([])
                self.plt1_em.set_ydata([])
                
            if ind_ab.size > 0:
                self.plt1_ab.set_xdata(u[ind_ab])
                self.plt1_ab.set_ydata(numpy.fabs(intIntensity[ind_ab]))
            else:
                self.plt1_ab.set_xdata([])
                self.plt1_ab.set_ydata([])

            #optical depth
            tau = numpy.array([l['tau'] for l in lines])
            mn, mx = [numpy.fabs(tau).min()/2, numpy.fabs(tau).max()*2]
            self.axs[2].set_ylim(mn, mx)
            self.axs[2].set_yticks(numpy.linspace(mn, mx, 4))
            self.axs[2].set_yticks(self.axs[2].get_yticks()[0:-1])
            
            self.plt2_all.set_xdata(u)
            self.plt2_all.set_ydata(numpy.fabs(tau))
            ind_em = numpy.where(tau > 0.0)[0]
            ind_ab = numpy.where(tau < 0.0)[0]
            if ind_em.size > 0:                
                self.plt2_em.set_xdata(u[ind_em])
                self.plt2_em.set_ydata(tau[ind_em])
            else:
                self.plt2_em.set_xdata([])
                self.plt2_em.set_ydata([])
                
            if ind_ab.size > 0:
                self.plt2_ab.set_xdata(u[ind_ab])
                self.plt2_ab.set_ydata(numpy.fabs(tau[ind_ab]))
            else:
                self.plt2_ab.set_xdata([])
                self.plt2_ab.set_ydata([])

            #optical depth dust
            tauDust = numpy.array([l['tauDust'] for l in lines])
            self.axs[3].set_ylim([numpy.fabs(tauDust).min()/2, numpy.fabs(tauDust).max()*2])
            
            mn, mx = [numpy.fabs(tauDust).min()/2, numpy.fabs(tauDust).max()*2]
            self.axs[3].set_ylim(mn, mx)
            self.axs[3].set_yticks(numpy.linspace(mn, mx, 4))
            self.axs[3].set_yticks(self.axs[3].get_yticks()[0:-1])

            self.plt3_all.set_xdata(u)
            self.plt3_all.set_ydata(numpy.fabs(tauDust))
            ind_em = numpy.where(tauDust > 0.0)[0]
            ind_ab = numpy.where(tauDust < 0.0)[0]
            if ind_em.size > 0:                
                self.plt3_em.set_xdata(u[ind_em])
                self.plt3_em.set_ydata(tauDust[ind_em])
            else:
                self.plt3_em.set_xdata([])
                self.plt3_em.set_ydata([])
                
            if ind_ab.size > 0:
                self.plt3_ab.set_xdata(u[ind_ab])
                self.plt3_ab.set_ydata(numpy.fabs(tauDust[ind_ab]))
            else:
                self.plt3_ab.set_xdata([])
                self.plt3_ab.set_ydata([])

            self.canvas.draw()   
            
    def run_cloud_gui():
        
        app = QtGui.QApplication(sys.argv)
        gui = cloud_gui()
        sys.exit(app.exec_())

    #fail
except Exception as E:

    def run_cloud_gui():
        print('Cannot use DESPOTIC in gui mode. Failed to import PyQt4, or matplotlib.backends.backend_qt4agg')
        print("The error was: "+str(E))
    

########################################################################
# Helper function to return the residuals for calculation of gas and
# dust temperatures. The function takes as its primary argument a
# two-element vector giving (Tgas, Tdust). Additional
# arguments are:
#    cloud -- the calling cloud
#    TCMB -- CMB temperature
#    TradDust -- radiation field temperature seen by the dust
#    zetaCR -- cosmic ray ionization rate
#    GammaPE0 -- photoelectric heating rate in *unattenuated* gas;
#               actual heating rate will be reduced due to dust
#               opacity
#    c1Grav -- coefficient to describe heating rate due to
#              gravitational compression; rate = c1Grav * c_s**2 *
#              sqrt(4 pi G rho) * mu_H * m_H
#    lumScale -- luminosity scale in the problem, used to normalize
#                the dEdt values calculated by this routine. This is a
#                two-element vector, giving scales for gas and dust
#                independently, since they can be quite different.
#    thin -- Boolean; if true, gas is assumed to be optically thin
#    LTE -- Boolean; if true, gas is assumed to be in LTE
#    escapeProbGoem -- string; specified geoemetry to be assumed in
#         computing escape probabilities
#    PsiUser -- an optional user-specified function that gives extra
#        heating and cooling
#    noClump -- Boolean; if true, clumping factors are turned off
########################################################################
def _gdTempResid(logTgd, cloud, c1Grav, thin, LTE,
                 escapeProbGeom, PsiUser, noClump, lumScale,
                 verbose):

    # Insert current dust and gas temperatures into cloud structure;
    # floor to CMB temperature to prevent numerical badness in case
    # the iterative solver has wantered off to negative values
    cloud.Tg = max(np.exp(logTgd[0]), cloud.rad.TCMB)
    cloud.Td = max(np.exp(logTgd[1]), cloud.rad.TCMB)
    if verbose:
        print("")
        print("***")
        print("setTempEq: calculating residual at Tg = " + 
              str(cloud.Tg) + " K, Td = " + str(cloud.Td) + " K...")

    # Get net heating / cooling rates
    rates = cloud.dEdt(c1Grav=c1Grav, thin=thin, LTE=LTE,
                       escapeProbGeom=escapeProbGeom,
                       sumOnly=True, PsiUser=PsiUser,
                       noClump=noClump,
                       verbose=verbose)

    # Print status
    if verbose:
        print("setTempEq: dE_gas/dt = " + str(rates['dEdtGas']) + 
              " erg s^-1 H^-1, dE_dust/dt = " + 
              str(rates['dEdtDust']) + " erg s^-1 H^-1, " + 
              "residuals = " + str(rates['dEdtGas']/lumScale[0]) + 
              " " + str(rates['dEdtDust']/lumScale[1]))

    # Return result
    return np.array([rates['dEdtGas'], rates['dEdtDust']]) \
        / lumScale


########################################################################
# Helper function to return the residuals for calculation of dust
# temperatures at fixed Tgas. The function takes the arguments:
#    Td -- dust temperature (float)
#    cloud -- the calling cloud
#    lumScale -- characteristic luminosity values
#    PsiUser -- an optional user-specified function that gives extra
#        heating and cooling
#    GammaSum -- sum of heating terms that don't depend on T_d
#    GammaSumMax -- maximum of absolute values of the heating terms
#        that go into GammaSum, needed for normalization
#    verbose -- verbose or not
########################################################################
def _dustTempResid(logTd, cloud, PsiUser, GammaSum, GammaSumMax, 
                   lumScale, verbose):

    # Insert current dust temperature into gas structure, with a floor
    # equal to the CMB floor to prevent numerical problems if the
    # rootfinder wanders into negative values
    cloud.Td = max(np.exp(logTd), small)

    # Get net heating / cooling rates
    rates = cloud.dEdt(dustCoolOnly=True, sumOnly=True,
                       PsiUser=PsiUser,
                       fixedLevPop=True,
                       verbose=verbose)

    # Print status if verbose
    if verbose:
        print("_dustTempResid called with Td = "+str(cloud.Td) + 
              ", dEdt = " + str(GammaSum+rates['dEdtDust']) + 
              ", residual = " + 
              str((GammaSum+rates['dEdtDust'])/lumScale))

    # Return dE/dt with correct normalization
    return (GammaSum+rates['dEdtDust'])/lumScale


########################################################################
# Helper function to return the residuals for calculation of gas
# temperatures at fixed Tdust. The function takes the arguments:
#    Tg -- gas temperature (float)
#    cloud -- the calling cloud
#    c1Grav -- coefficient to describe heating rate due to
#              gravitational compression; rate = c1Grav * c_s**2 *
#              sqrt(4 pi G rho) * mu_H * m_H
#    thin -- Boolean; if true, gas is assumed to be optically thin
#    LTE -- Boolean; if true, gas is assumed to be in LTE
#    escapeProbGoem -- string; specified geoemetry to be assumed in
#        computing escape probabilities
#    PsiUser -- user-specified function to add extra heating / cooling
#        terms
#    noClump -- Boolean; if true, clumping factors are turned off
#    verbose -- verbose or not
########################################################################
def _gasTempResid(logTg, cloud, c1Grav, thin, LTE,
                  escapeProbGeom, PsiUser, noClump, lumScale,
                  verbose):

    # Insert current dust temperature into gas structure, using CMB
    # temperature as a floor
    cloud.Tg = max(np.exp(logTg), cloud.rad.TCMB)

    # Get net heating / cooling rates
    rates = cloud.dEdt(c1Grav=c1Grav,
                       thin=thin, LTE=LTE,
                       escapeProbGeom=escapeProbGeom,
                       gasOnly=True, noClump=noClump,
                       sumOnly=True, PsiUser=PsiUser)

    # If verbose, print status
    if verbose:
        print ("_gasTempResid called with Tg = "+str(cloud.Tg) + 
               ", dEdt = " + str(rates['dEdtGas']) + 
               ", residual = " + 
               str(rates['dEdtGas']/lumScale))

    # Return dE/dt with correct normalization
    return rates['dEdtGas']/lumScale


########################################################################
# Helper function to return dEdt for gas and dust in a form that the
# ODE integrator is happy with. Arguments are:
#    Tg -- current gas temperature
#    time -- current time
#    cloud -- the calling cloud
#    c1Grav -- coefficient to describe heating rate due to
#              gravitational compression; rate = c1Grav * c_s**2 *
#              sqrt(4 pi G rho) * mu_H * m_H
#    thin -- Boolean; if true, gas is assumed to be optically thin
#    LTE -- Boolean; if true, gas is assumed to be in LTE
#    escapeProbGoem -- string; specified geoemetry to be assumed in
#         computing escape probabilities
#    PsiUser -- user-specified function to add extra heating / cooling
#        terms
#    noClump -- Boolean; if true, clumping factors are turned off
#    isobar -- float; if > 0, this gives the fixed value of n*T; if
#        <= 0, this indicates the calculation is done isochorically,
#        i.e. n = constant
########################################################################
def _gasTempDeriv(Tg, time, cloud, c1Grav, thin, LTE, \
                      escapeProbGeom, PsiUser, noClump, isobar, \
                      dampFactor, verbose):

    # Update temperature
    cloud.Tg = Tg[0]

    # If we're isobaric, update density
    if isobar > 0:
        cloud.nH = isobar / cloud.Tg

    if verbose:
        print("t = "+str(time)+": Tg = "+str(cloud.Tg) + 
              ", nH = "+str(cloud.nH))

    # Compute new dust temperature
    cloud.setDustTempEq(PsiUser=PsiUser, verbose=verbose,
                        dampFactor=dampFactor)

    # Call dEdt to get time rate of change of energy for gas; note
    # that the level populations do not need to be recomputed here,
    # because they were already computed in solving for the dust
    # temperature value, which depends on the line heating rate
    rates = cloud.dEdt(c1Grav=c1Grav,
                       thin=thin, LTE=LTE,
                       escapeProbGeom=escapeProbGeom,
                       gasOnly=True, noClump=noClump,
                       sumOnly=True, PsiUser=PsiUser,
                       fixedLevPop=True)
    dEdtGas = rates['dEdtGas']

    # Compute specific heat c_v at current temperature
    cloud.comp.computeCv(cloud.Tg)

    # Convert from dE/dt to dTemp/dt by dividing by the c_V
    if isobar > 0:
        if verbose:
            print("dE/dt = "+str(dEdtGas)+", dT/dt = " + 
                  str(dEdtGas/((cloud.comp.cv+1.0)*kB)))
        return dEdtGas/((cloud.comp.cv+1.0)*kB)
    else:
        if verbose:
            print("dE/dt = "+str(dEdtGas)+", dT/dt = " + 
                  str(dEdtGas/((cloud.comp.cv)*kB)))
        return dEdtGas/(cloud.comp.cv*kB)

