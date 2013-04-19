"""
This module provides the emitter class, which describes a single
emitter species, and defines methods for performing calculations on
this species. It also stores a master list of emitting species for
which the corresponding emitterData has been created; this is the dict
knownEmitterData.
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
import os
from despoticError import despoticError
from emitterData import emitterData

########################################################################
# Global variables
########################################################################

# Define some global physical constants in cgs units
import scipy.constants as physcons
kB = physcons.k/physcons.erg
c = physcons.c/physcons.centi
mH = (physcons.m_p+physcons.m_e)/physcons.gram
h = physcons.h*1e7
hc = h*c

# Small numerical value, to avoid divide by 0 issues
small = 1e-50

# Machine epsilon
machineeps = np.finfo(float).eps

# This dict holds emitter data that we've already read in, so that we
# don't need to read it a second time if we create another emitter
# object with the same data
knownEmitterData = {}

########################################################################
# class emitter definition
########################################################################
class emitter:
    """
    Class to store the properties of a single emitting species, and
    preform computations on those properties. Note that all quantities
    are stored in cgs units.

    Attributes
    ----------
    name : string
        name of emitting species
    abundance : float
        abundance of emitting species relative to H nuclei
    data : class emitterData
        physical data for this emitter
    levPop : array(nlev)
        array giving populations of each level
    levPopInitialized : Boolean
        flag for whether levPop is initialized or uninitialized
    escapeProb : array(nlev, nlev)
        2d array giving escape probability for photons emitted in a
        line connecting the given level pair 
    escapeProbInitialized : Boolean
        flag for whether escapeProb is initialized or uninitialized
    energySkip : Boolean
        flag that this species should be skipped when computing line
        cooling rates

    Methods
    -------------
    __init__ -- initialize
    __deepcopy__ -- override the default deepcopy operation
    __getstate__ -- routine to provide pickling support
    __setstate__ -- routine to provide pickling support
    setLevPopLTE -- set level populations to their LTE values at the
        given temperature
    setThin -- reset all escape probabilities to unity, the value for
        an optically thin gas
    setLevPop -- calculate the level populations from the currently
        stored escape probabilities, together with the input gas
        physical conditions and composition
    setEscapeProb -- calculate the escape probabilties from the
        currently-stored level populations and the specified gas
        properties.
    setLevPopEscapeProb -- compute the level populations and escape
        probabilities simultaneously for the specified gas physical
        properties and composition
    opticalDepth -- return the optical depth associated with the
        currently-stored escape probabilities
    luminosityPerH -- compute the luminosity per H nucleus due to line
        emission
    setExtrap -- turns extrapolation on or off for this emitter
    """

########################################################################
# Class initialization method
########################################################################
    def __init__(self, emitName, emitAbundance, extrap=False, \
                     energySkip=False, emitterFile=None, \
                     emitterURL=None):

        """
        Initialization routine

        Parameters
        ----------
        emitName : string
            name of this species
        emitAbundance : float
            abundance of species relative to H
        emitterFile : string
            name of LAMDA file containing data on this species; this
            option overrides the default
        emitterURL : string
            URL of LAMDA file containing data on this species; this
            option overrides the default
        energySkip : Boolean
            if True, this species is skipped when computing line
            cooling rates
        extrap : Boolean
            if True, collision rate coefficients for this species will
            be extrapolated to temperatures outside the range given in
            the LAMDA tables

        Returns
        -------
        Nothing
        """

        # Basic data
        self.name = emitName
        self.abundance = emitAbundance
        self.energySkip = energySkip

        # Get emitter data
        if emitName in knownEmitterData:
            # Emitter data in known list, so just set a pointer
            self.data = knownEmitterData[emitName]
        else:
            # Emitter data not in known list, so read and then set a
            # pointer
            self.data = emitterData(emitName, emitterFile=emitterFile, \
                                        emitterURL=emitterURL, \
                                        extrap=extrap)
            knownEmitterData[emitName] = self.data
 
        # Initialize level populations and escape probabilities, and
        # set flag
        self.levPopInitialized = False
        self.levPop = np.zeros(self.data.nlev)
        self.escapeProbInitialized = False
        self.escapeProb = np.zeros((self.data.nlev, self.data.nlev)) \
            + 1.0


########################################################################
# Override deepcopy method
# We override the default deepcopy because we don't want to deepcopy
# emitterData or collPartner objects that store physical constants,
# and should be shared by all instances of an emitter class that
# represent the same atom or molecule. This maintains the generic
# behavior of this module, i.e. there is a centrally-maintained list
# of emitter data, and individual emitter instances only contain
# pointers to it, they don't each replicate it.
########################################################################
    def __deepcopy__(self, memo={}):
        """
        This defines the deepcopy method for emitting species; a
        custom deepcopy is required to handle the emitterData.

        Parameters
        ----------
        memo : dict
            The dict of objects that have been copied

        Returns
        -------
        Nothing
        """

        # Create a new emitter with the same name, abundance, and
        # energySkip setting as this one
        newcopy = emitter(self.name, self.abundance, \
                              energySkip = self.energySkip)

        # Copy all the attributes to the new emitter except the
        # emitterData
        newcopy.levPopInitialized = self.levPopInitialized
        newcopy.escapeProbInitialized = self.escapeProbInitialized
        newcopy.levPop = self.levPop.copy()
        newcopy.escapeProb = self.escapeProb.copy()

        # Add the new object we've created to the memo
        memo[id(self)] = newcopy

        # Return the new object
        return newcopy


########################################################################
# getstate and setstate methods
# We need custom getstate and setstate methods to handle pickling of
# emitters, because the emitterData requires special handling. These
# data should not be pickled, as they are centrally stored and shared
# by all instances of the same emitter. When pickle dump is called, the
# data are removed from the list of attributes to be written, and when
# load is called they are re-read if necessary. However, there are
# certain attributes of the emitterData that we do need to store: the
# extrapolation state and the file name, since these cannot be
# determined just from re-reading the LAMDA file. These are added to the
# dictionary before pickling, and removed after pickling.
########################################################################
    def __getstate__(self):
        """
        This method provides custom pickling behavior for this class;
        custom behavior is required to handle the corresponding
        emitterData object. This method removes the emitterData from
        the list of items to be pickled, and adds the information
        needed to re-read this data when the class is unpickled.

        Parameters
        ----------
        None

        Returns
        -------
        odict : dict
           a dict for this object that is suitable for pickling, with
           certain information added and some removed
        """
        odict = self.__dict__.copy() # copy the dict since we change it
        odict['extrap'] = self.data.extrap
        odict['lamdaFile'] = self.data.lamdaFile
        del odict['data']            # remove entry for emitter data
        return odict

    def __setstate__(self, idict):
        """
        This method provides custom pickling behavior for this class;
        custom behavior is required to handle the corresponding
        emitterData object. This method re-reads the requisite
        emitterData, then removes the extraneous information that was
        used to describe it from the dict for the class.

        Parameters
        ----------
        idict : dict
            dict describing the class as read from the pickle

        Returns
        -------
        Nothing
        """

        # First extract the extraneous parts from the dict
        extrap = idict['extrap']
        lamdaFile = idict['lamdaFile']
        del idict['extrap']
        del idict['lamdaFile']

        # Copy remaining attributes to self object
        self.__dict__.update(idict)

        # Now we need to reconstruct the emitterData for this object

        # Check if the emitter data for this emitter is known
        if self.name in knownEmitterData:

            # Emitter data in known list, so just set a pointer
            self.data = knownEmitterData[self.name]

        else:

            # Emitter data not in known list. This is tricky. We will
            # first try looking for the file at the stored
            # location. However, if this fails, we will go looking for
            # a generic LAMDA file containing the data for this
            # emitter in the standard location. This additional test
            # is done to maximize portability, so that it should be
            # possible to take a pickled object that was created on
            # one platform and unpickle it on another platform, as
            # long as the other platform has access the the LAMDA
            # data.

            try:
                # Look for data file in stored location
                self.data = emitterData(self.name, \
                                            emitterFile=lamdaFile, \
                                            extrap=extrap)
            except despoticError:
                # We didn't find the file in the stored location, so
                # try looking for it as we would any other data file,
                # without specifying a location
                self.data = emitterData(self.name, \
                                            extrap=extrap)

            # Store emitter in list
            knownEmitterData[self.name] = self.data


########################################################################
# Method to turn extrapolation on or off
########################################################################
    def setExtrap(self, extrap):
        """
        Turn extrapolation on or off for this emitter

        Parameters
        ----------
        extrap : Boolean
           true = extrapolation on, false = extrapolation off

        Returns
        -------
        Nothing

        Remarks
        -------
        Since emitter data is shared by all emitters of the same
        species, turning extrapolation on or off affects all instances
        of emitters of this species, not just this instance.
        """

        # Check if new setting is different from current one; if not,
        # do nothing
        if extrap == self.data.extrap:
            return

        # Rebuild interpolation functions and store new settings
        for p in self.data.partners.values():
            p.buildInterp(extrap)
        self.data.extrap = extrap


########################################################################
# Method to set level populations to LTE values
########################################################################
    def setLevPopLTE(self, temp):
        """
        Set the level populations of this species to their LTE values

        Parameters
        ----------
        temp : float
            temperature in K

        Returns
        -------
        Nothing
        """

        # Partition function
        z = self.data.partFunc(temp)
        # Level populations
        self.levPop = self.data.levWgt * \
            np.exp(-self.data.levTemp/temp) / z
        # Floor
        self.levPop[self.levPop < small] = small
        # Flag that level populations are now set to valid value
        self.levPopInitialized = True


########################################################################
# Method to set escape probabilities to unity
########################################################################
    def setThin(self):
        """
        Set the escape probabilities for this species to unity

        Parameters
        ----------
        None

        Returns
        -------
        Nothing
        """
        self.escapeProb[:,:] = 1.0
        self.escapeProbInitialized = True


########################################################################
# Method to set level populations in statistical equilibrium from the
# stored escape probabilities. An optional flag tells the code to
# assume that the gas is optically thin, treat all escape
# probabilities as unity.
########################################################################
    def setLevPop(self, thisCloud, thin=False, noClump=False, \
                      diagOnly=False, verbose=False):
        """
        Compute the level populations for this species using the
        stored escape probabilities

        Parameters
        ----------
        thisCloud : class cloud
            a cloud object containing the physical and chemical
            properties of this cloud
        thin : Boolean
            if True, the stored escape probabilities are ignored, and
            the cloud is assumed to be optically thin (equivalent to
            assuming all escape probabilities are 1)
        noClump : Boolean
           if set to True, the clumping factor used in estimating
           rates for n^2 processes is set to unity

        Returns
        -------
        infoDict : dict
            A dictionary containing a variety of diagnostic
            information. See user manual for details.

        Other parameters
        ----------------
        diagOnly : Boolean
            if true, diagnostic information is returned, but no
            attempt is made to solve the equations or calculate the
            level popuplations (useful for debugging)
        """

        # Initialize dict of diagnostic information
        infoDict = {}

        # Get collision rate matrix
        q = self.data.collRateMatrix(thisCloud.nH, \
                                         thisCloud.comp, \
                                         thisCloud.Tg)
        infoDict['qNoClump'] = q

        # Apply clumping factor
        if noClump == False:
            cs2 = kB * thisCloud.Tg / (thisCloud.comp.mu * mH)
            cfac = np.sqrt(1.0 + 0.75*thisCloud.sigmaNT**2/cs2)
            q *= cfac
        infoDict['q'] = q

        # Get photon occupation number
        if thisCloud.rad.TCMB > 0.0:
            ngamma = np.zeros((self.data.nlev, self.data.nlev))
            ngamma[self.data.radUpper, self.data.radLower] = \
                thisCloud.rad.ngamma(self.data.radTemp)
        else:
            ngamma = 0.0
        infoDict['ngamma'] = ngamma

        # Construct transition rate matrix:
        #    M_ij = (q_ji + beta_ji (1 + ngamma,ji) A_ji + \
        #       (g_i/g_j) beta_ij n_gamma,ij A_ij) / 
        #       Sum_l [q_il + beta_il (1 + n,gamma,il) A_il + \
        #          (g_l/g_i) beta_li n_gamma,li A_li]
        if thin == False:
            m = np.transpose(q) + \
                np.transpose((1.0+ngamma) * self.escapeProb * \
                              self.data.EinsteinA) + \
                              self.data.wgtRatio * ngamma * \
                              self.escapeProb * \
                              self.data.EinsteinA
        else:
            m = np.transpose(q) + \
                np.transpose((1.0+ngamma) * self.data.EinsteinA) + \
                self.data.wgtRatio * ngamma * self.data.EinsteinA
        msave=m.copy() # save a deep copy so we can do row-by-row manipulation
        for i in range(self.data.nlev):
            m[i,:]=m[i,:]/(np.sum(msave[:,i])+small)

        # Store diagnostic information
        infoDict['inRateCoef'] = msave

        # We want to solve the matrix equation M x = x for the level
        # popluations x, with the constraint that Sum(x) = 1.0. The
        # most efficient way to do this is to rewrite the equation as
        # (M - I) x = 0, and then add a row to the matrix M - I
        # giving the constraint equation, then use a least-squares
        # solve.

        # Step 1. Build M - I
        m1=m-np.identity(self.data.nlev)

        # Step 2. Build a matrix m2 with one extra row in which every
        # element of x is just multiplied by 1 (i.e. the row is all
        # 1's)
        m2=np.zeros((self.data.nlev+1, self.data.nlev))
        m2[:-1,:]=m1
        m2[-1,:]=1.0

        # Step 3. Construct the right hand side of the equation to be
        # solved, (M - I) x = 0, plus the extra RHS for the equation
        # Sum(x) = 1.
        b=np.zeros(self.data.nlev+1)
        b[-1]=1.0
        infoDict['m'] = m2

        # Step 4. Check condition number of matrix. If it is
        # ill-conditioned, eliminate levels from the
        # calculation. First remove levels for which the population in
        # negligible in LTE with either the gas or the radiation
        # field. Then go through and recursively eliminate levels with
        # negligible transition rates into them.
        conditionFlag = False
        if np.linalg.cond(m2) > 1.0e8:

            # Flag that we've eliminated some levels
            conditionFlag = True

            # Track levels we're keeping; initially this is all
            levKeep = np.arange(self.data.nlev, dtype='int')

            # Calculate populations in LTE with gas and radiation
            levPopLTEGas = self.data.levWgt * \
                np.exp(-self.data.levTemp/thisCloud.Tg) / \
                self.data.partFunc(thisCloud.Tg)
            levPopLTERad = self.data.levWgt * \
                np.exp(-self.data.levTemp/thisCloud.rad.TCMB) / \
                self.data.partFunc(thisCloud.rad.TCMB)

            # Kill levels that are below the minimum for both
            levPopMax = np.maximum(levPopLTEGas, levPopLTERad)
            levKeep = np.delete(levKeep, np.where(levPopMax < machineeps))
            levDel = np.setdiff1d(np.arange(self.data.nlev, dtype='int'), levKeep)

            # Build new matrix with some levels removed
            m = msave
            m=np.delete(m, levDel, axis=0)
            m=np.delete(m, levDel, axis=1)
            msave1=m.copy()
            m1=m.copy()
            if len(levKeep) > 1:
                for i in range(len(levKeep)):
                    m[i,:]=m[i,:]/(np.sum(m1[:,i])+small)
            m1=m-np.identity(len(levKeep))
            m2=np.zeros((len(levKeep)+1, len(levKeep)))
            m2[:-1,:]=m1
            m2[-1,:]=1.0
            b = np.zeros(len(levKeep)+1)
            b[-1] = 1.0

            # Store diagnostic information
            infoDict['levDel'] = levDel
            infoDict['inRateCoefRed'] = msave1
            infoDict['mRed'] = m2

            # Check condition number again. If it's still too high,
            # find levels with large ratios of coefficients describing
            # transitions in to transitions out, and eliminate
            # them. Repeat this iteratively until the condition number
            # is acceptable.
            if np.linalg.cond(m2) > 1.0/(1.0e5*machineeps):

                # Get upper limits on level populations
                sumRateOut = np.sum(msave1, axis=0)
                sumRateIn = np.sum(msave1, axis=1)
                levPopLim = np.sumRateIn / sumRateOut

                while np.linalg.cond(m2) > 1.0/(1.0e5*machineeps) and \
                        np.argmin(levPopLim) < 1.0e-6:

                    # Find smallest level population limit left; flag
                    # it to make sure we don't hit it again
                    idx = np.argmin(levPopLim)
                    levPopLim[idx] = 1e100

                    # Delete the offending level
                    idxdel = idx - len(levDel[levDel < idx])
                    if verbose:
                        print idx
                        print idxdel
                        print levKeep
                        print levDel

                    levKeep = np.delete(levKeep, idxdel)
                    levDel = np.setdiff1d(np.arange(self.data.nlev, \
                                                  dtype='int'), \
                                           levKeep)

                    # Rebuild M and b with levels removed
                    m = msave
                    m=np.delete(m, levDel, axis=0)
                    m=np.delete(m, levDel, axis=1)
                    msave1=m.copy()
                    m1=m.copy()
                    if len(levKeep) > 1:
                        for i in range(len(levKeep)):
                            m[i,:]=m[i,:]/(np.sum(m1[:,i])+small)
                    m1=m-np.identity(len(levKeep))
                    m2=np.zeros((len(levKeep)+1, len(levKeep)))
                    m2[:-1,:]=m1
                    m2[-1,:]=1.0
                    b = np.zeros(len(levKeep)+1)
                    b[-1] = 1.0

                    # Store diagnostic information
                    infoDict['levDel2'] = levDel
                    infoDict['inRateCoefRed2'] = msave1
                    infoDict['mRed2'] = m2

        # Step 5. Do the non-linear least squares solve to get the
        # level populations. Handle special case of only on level
        # present by just setting the population of that level to
        # unity.
        if diagOnly == True:
            return infoDict
        if len(b) > 2:
            self.levPop, res, rank, s = np.linalg.lstsq(m2, b)
        else:
            self.levPop = array([1.0])

        # Step 5a. If we eliminated levels due to ill-conditioning,
        # put them back in now, with populations of 0
        if conditionFlag:
            levPop = self.levPop.copy()
            self.levPop = np.zeros(self.data.nlev)
            self.levPop[levKeep] = levPop

        # Step 6. Floor to avoid numerical issues
        self.levPop[self.levPop < small] = small

        # Flag that level populations are now set
        self.levPopInitialized = True

        # Return diagnostic information
        return infoDict


########################################################################
# Method to compute the escape probabilities from the stored level
# populations. By default this computation happens for all
# transitions, but optional keywords can specify a specific [upper,
# lower] combination.
########################################################################
    def setEscapeProb(self, thisCloud, transition=None, \
                          escapeProbGeom='sphere'):
        """
        Compute escape probabilities from stored level populations

        Parameters
        ----------
        thisCloud : class cloud
            a cloud object containing the physical and chemical
            properties of this cloud
        transition : sequence of two sequences of int
                     transition[0] = array of upper states
                     transition[1] = array of lower states
        escapeProbGeom : string, 'sphere' or 'LVG' or 'slab'
             sets problem geometry that will be assumed in calculating
             escape probabilities

        Returns
        -------
        Nothing

        Raises
        ------
        despoticError, if the level populations are not initialized or
        the escape probability geometry is not 'sphere', 'slab', or
        'LVG'
        """

        if self.levPopInitialized == False:
            raise despoticError, "cannot compute escape probabilities before level populations are initialized; call one of the setLevPop routines first"
        if transition==None:
            u = self.data.radUpper
            l = self.data.radLower
        else:
            u = transition[0]
            l = transition[1]
        if escapeProbGeom == "sphere":
            sigmaTot = np.sqrt(thisCloud.sigmaNT**2 + \
                                kB*thisCloud.Tg/(self.data.molWgt*mH))
            tau = (self.data.levWgt[u]/self.data.levWgt[l]) * \
                self.data.EinsteinA[u,l] * (c/self.data.freq[u,l])**3 * \
                3.0 / (16.0*(2.0*np.pi)**1.5*sigmaTot) * \
                self.abundance * thisCloud.colDen * self.levPop[l] * \
                (1.0 - self.levPop[u]*self.data.levWgt[l] / \
                     (self.levPop[l]*self.data.levWgt[u]))
            self.escapeProb[u,l] = 1.0 / (1.0+0.5*tau)
            # The following is RADEX's approximation
            #self.escapeProb[u,l] = 1.5/tau * (1-2/tau**2+(2/tau+2/tau**2)*exp(-tau))
        elif escapeProbGeom == "LVG":
            tau = (self.data.levWgt[u]/self.data.levWgt[l]) * \
                self.data.EinsteinA[u,l] * (c/self.data.freq[u,l])**3 / \
                (8.0*np.pi*abs(thisCloud.dVdr)) * \
                self.abundance * thisCloud.nH * self.levPop[l] * \
                (1.0 - self.levPop[u]*self.data.levWgt[l] / \
                     (self.levPop[l]*self.data.levWgt[u]))
            # Note that, in computing escape probabilities, we need to
            # handle the case of very small tau with care to avoid
            # roundoff problems and make sure we correctly limit to
            # beta -> 1 as tau -> 0
            idx = tau > 1e-6
            self.escapeProb[u[idx], l[idx]] = \
                (1.0 - np.exp(-tau[idx])) / tau[idx]
            idx = tau <= 1e-6
            self.escapeProb[u[idx], l[idx]] = \
                1.0 - tau[idx]/2.0
        elif escapeProbGeom == "slab":
            sigmaTot = np.sqrt(thisCloud.sigmaNT**2 + \
                                kB*thisCloud.Tg/(self.data.molWgt*mH))
            tau = (self.data.levWgt[u]/self.data.levWgt[l]) * \
                self.data.EinsteinA[u,l] * (c/self.data.freq[u,l])**3 / \
                (4.0*(2.0*np.pi)**1.5*sigmaTot) * \
                self.abundance * thisCloud.colDen * self.levPop[l] * \
                (1.0 - self.levPop[u]*self.data.levWgt[l] / \
                     (self.levPop[l]*self.data.levWgt[u]))
            idx = tau > 1e-6
            self.escapeProb[u[idx], l[idx]] = \
                (1.0 - np.exp(-3.0*tau[idx])) / (3.0*tau[idx])
            idx = tau <= 1e-6
            self.escapeProb[u[idx], l[idx]] = \
                1.0 - 3.0*tau[idx]/2.0
        else:
            raise despoticError, "unknown escapeProbGeom " + \
                str(escapeProbGeom)

        # Make negative escape probabilities positive to avoid
        #numerical problems
        self.escapeProb[self.escapeProb < 0.0] = \
            -self.escapeProb[self.escapeProb < 0.0]

        self.escapeProb[l,u] = self.escapeProb[u,l]
        self.escapeProbInitialized = True


########################################################################
# Method to compute the escape probabilities and level populations
# simultaneously.
########################################################################
    def setLevPopEscapeProb(self, thisCloud, escapeProbGeom='sphere', \
                                noClump=False, verbose=False, \
                                reltol=1e-6, abstol=1e-8, \
                                maxiter=200, veryverbose=False, \
                                dampFactor=0.5):
        """
        Compute escape probabilities and level populations
        simultaneously.

        Parameters
        ----------
        thisCloud : class cloud
            a cloud object containing the physical and chemical
            properties of this cloud
        escapeProbGeom : string, 'sphere' or 'LVG' or 'slab'
           sets problem geometry that will be assumed in calculating
           escape probabilities
        noClump : Boolean
           if set to True, the clumping factor used in estimating
           rates for n^2 processes is set to unity

        Returns
        -------
        success: Boolean
           true if iteration converges, false if it does not

        Additional Parameters
        ---------------------
        verbose : Boolean
           if true, diagnostic information is printed
        veryverbose : Boolean
           if true, a very large amount of diagnostic information is
           printed; probably useful only for debugging
        reltol : float
           relative tolerance; convergence is considered to have
           occured when |f_i(n+1) - f_i(n)| / 
           max(f_i(n+1), f_i(n)) < reltol
        abstol : float
           absolute tolerance; convergence is considered to have
           occured when |f_i(n)+1 - f_i(n)| < abstol
        maxiter : int
           maximum number of iterations to allow
        dampFactor : float
           a number in the range (0, 1] that damps out changes in level
           populations at each iteration. A value of 1 means no
           damping, while a value of 0 means the level populations
           never change.

        Raises
        ------
        despoticError, if the escape probability geometry is not
        'sphere', 'slab', or 'LVG'

        Remarks
        -------
        Convergence occurs when either the relative or the absolute
        tolerance condition is satisfied. To disable either relative
        or absolute tolerance checking, just set the appropriate
        tolerance <= 0. However, be warned that in many circumstances
        disabling absolute tolerances will gaurantee non-convergence,
        because truncation errors tend to produce large relative
        residuals for high energy states whose populations are very
        low, and no amount of iterating will reduce these errors
        substantially.
        """

        # Get collision rate matrix with clumping factor correction
        qcoltrans = self.data.collRateMatrix(thisCloud.nH, \
                                                 thisCloud.comp, \
                                                 thisCloud.Tg)
        if noClump == False:
            cs2 = kB * thisCloud.Tg / (thisCloud.comp.mu * mH)
            cfac = np.sqrt(1.0 + 0.75*thisCloud.sigmaNT**2/cs2)
            qcoltrans *= cfac
        qcol = np.transpose(qcoltrans)

        # Get photon occupation number
        if thisCloud.rad.TCMB > 0.0:
            ngamma = np.zeros((self.data.nlev, self.data.nlev))
            ngamma[self.data.radUpper, self.data.radLower] = \
                thisCloud.rad.ngamma(self.data.radTemp)
        else:
            ngamma = 0.0

        # Get radiative transition rate matrix without escape
        # probability corrections
        qrad = np.transpose((1.0+ngamma) * self.data.EinsteinA) + \
                self.data.wgtRatio * ngamma * self.data.EinsteinA

        # If current level populations are unitialized, set them to
        # their LTE values as an initial guess
        if self.levPopInitialized == False:
            self.setLevPopLTE(thisCloud.Tg)
            #self.setLevPop(thisCloud, thin=True, noClump=noClump)

        # Create workspace matrices, and initialize the last rows
        # M and b, which hold the summation constraint
        m = np.zeros((self.data.nlev+1, self.data.nlev))
        m1 = np.zeros((self.data.nlev, self.data.nlev))
        b = np.zeros(self.data.nlev+1)
        m[-1,:] = 1.0
        b[-1] = 1.0

        # prepare to iterate -- copy starting level population vector,
        # initialize norms and counter, create matrix of the correct
        # size
        levPopOld = self.levPop.copy()
        ctr = 1
        relnorm = 1.0
        absnorm = 1.0
        if veryverbose:
            lastresid = 1.0

        # iterate until converged
        while relnorm > reltol and absnorm > abstol and ctr < maxiter:

            # damp changes from last cycle
            self.levPop = dampFactor*self.levPop + \
                (1.0-dampFactor)*levPopOld

            # store last iterates
            levPopOld = self.levPop.copy()

            # compute new escape probabilities
            self.setEscapeProb(thisCloud, \
                                   escapeProbGeom=escapeProbGeom)

            # calculate new transition rate matrix
            m1 = qcol + self.escapeProb*qrad
            m[:-1,:] = np.transpose(np.transpose(m1)/(np.sum(m1,axis=0)+small)) - \
                np.identity(self.data.nlev)

            # check condition number
            if np.linalg.cond(m) <= 1.0/(1.0e5*machineeps):

                # Condition number ok, so solve
                self.levPop, res, rank, s = np.linalg.lstsq(m, b)

            else:

                # Condition numberis not ok, so we need to eliminate
                # some levels

                # Track levels we're keeping; initially this is all
                levKeep = np.arange(self.data.nlev, dtype='int')
                
                # Calculate populations in LTE with gas and radiation
                levPopLTEGas = self.data.levWgt * \
                    np.exp(-self.data.levTemp/thisCloud.Tg) / \
                    self.data.partFunc(thisCloud.Tg)
                levPopLTERad = self.data.levWgt * \
                    np.exp(-self.data.levTemp/thisCloud.rad.TCMB) / \
                    self.data.partFunc(thisCloud.rad.TCMB)

                # Kill levels that are below the minimum for both
                levPopMax = np.maximum(levPopLTEGas, levPopLTERad)
                levKeep = np.delete(levKeep, np.where(levPopMax < machineeps))
                levDel = np.setdiff1d(np.arange(self.data.nlev, dtype='int'), levKeep)

                # Build reduced m
                mred1=np.delete(m1, levDel, axis=0)
                mred1=np.delete(mred1, levDel, axis=1)
                mred = np.zeros((len(levKeep)+1, len(levKeep)))
                mred[:-1,:] = np.transpose(np.transpose(mred1) / \
                                            (np.sum(mred1,axis=0)+small)) - \
                                            np.identity(len(levKeep))
                mred[-1,:] = 1.0

                # Check condition number again
                if np.linalg.cond(mred) > 1.0/(1.0e5*machineeps):

                    # Condition number still bad, so progressively
                    # eliminate levels based on population upper
                    # limits.

                    # Get level limits
                    sumRateOut = np.sum(mred1, axis=0)
                    sumRateIn = np.sum(mred1, axis=1)
                    levPopLim = sumRateIn / sumRateOut

                    # Loop until condition number is acceptable
                    while np.linalg.cond(mred) > 1.0/(1.0e5*machineeps) \
                            and np.argmin(levPopLim) < 1.0e-6:

                        # Delete the level with the most stringent
                        # limit; reset the limit for this level to
                        # something large
                        idx = np.argmin(levPopLim)
                        idxDel = idx - len(levDel[levDel < idx])
                        levPopLim[idx] = 1e100
                        levKeep = np.delete(levKeep, idxDel)
                        levDel = np.setdiff1d( \
                            np.arange(self.data.nlev, dtype='int'), \
                                levKeep)

                        # Build a new reduced matrix
                        mred1=np.delete(m1, levDel, axis=0)
                        mred1=np.delete(mred1, levDel, axis=1)
                        mred = np.zeros((len(levKeep)+1, len(levKeep)))
                        mred[:-1,:] = \
                            np.transpose(np.transpose(mred1) / \
                                          (np.sum(mred1,axis=0)+small)) - \
                                          np.identity(len(levKeep))
                        mred[-1,:] = 1.0

                # We finally have a well-conditioned matrix, so
                # now construct the RHS and solve
                b = np.zeros(len(levKeep)+1)
                b[-1] = 1.0
                levPop, res, rank, s = np.linalg.lstsq(mred, b)
                self.levPop[levKeep] = levPop
                self.levPop[levDel] = 0.0

            # floor level populations to avoid divide by 0 errors
            self.levPop[self.levPop < small] = small

            # get absolute and relative residuals
            resid = abs(levPopOld-self.levPop)
            relResid = resid / np.maximum(levPopOld, self.levPop)
            absnorm = np.amax(resid)
            relnorm = np.amax(relResid)
            # print status message
            if veryverbose:
                print 'setLevPopEscapeProb for ' + self.name + \
                    ': iter ' + str(ctr) + ': abs resid = ' + \
                    str(absnorm) + ', state; ' + str(argmax(resid)) + \
                    ', rel resid = ' + str(relnorm) + \
                    ', state = ' + str(argmax(relResid))
            # update counter
            ctr = ctr+1

        # make sure we've converged
        if relnorm > reltol and absnorm > abstol:
            if verbose:
                print "setLevPopEscapeProb: level populations " + \
                    "DID NOT CONVERGE for " + self.name + \
                    " in "+str(ctr)+" iterations; " + \
                    "relative residual = " + str(relnorm) + \
                    ", absolute residual = " + str(absnorm)
            return False

        # print status message
        if verbose:
            print "setLevPopEscapeProb: level populations " + \
                "converged for " + self.name + \
                " in "+str(ctr)+" iterations; " + \
                "relative residual = " + str(relnorm) + \
                ", absolute residual = " + str(absnorm)

        # flag that level populations and escape probabilities are now
        # valid
        self.levPopInitialized = True
        self.escapeProbInitialized = True

        # return success
        return True


########################################################################
# Method to return the optical depth associated with all transitions
# (or with the specified transitions if the transition keyword is
# specified) based on the stored escape probabilities.
########################################################################
    def opticalDepth(self, transition=None, escapeProbGeom='sphere'):
        """
        Return the optical depths of various lines, computed from the
        stored escape probabilities.

        Parameters
        ----------
        transition : list of two arrays of shape (M)
                     transition[0] = array of upper states
                     transition[1] = array of lower states
        escapeProbGeom : string, 'sphere' or 'LVG' or 'slab'
             sets problem geometry that will be assumed in calculating
             escape probabilities

        Returns
        -------
        tau : array, shape (M)
              optical depths in specified lines; by default
              M = len(radUpper)

        Raises
        ------
        despoticError, if escape probabilities are not initialized or
        the escape probability geometry is not 'sphere', 'slab', or
        'LVG'
        """

        if self.escapeProbInitialized == False:
            raise despoticError, "Error: cannot compute optical depths before escape probabilities are initialized"
        if transition==None:
            u = self.data.radUpper
            l = self.data.radLower
        else:
            u = transition[0]
            l = transition[1]
        if escapeProbGeom == 'sphere':
            return 2.0*(1.0/self.escapeProb[u,l] - 1.0)
        elif escapeProbGeom == 'LVG':
            tau=np.zeros(len(u))
            for i, beta in enumerate(self.escapeProb[u,l]):
                tau[i] = brentq(_betaTauLVG, 1e-10, 1.0/beta, \
                                    args=(beta,))
            return tau
        elif escapeProbGeom == 'slab':
            tau=np.zeros(len(u))
            for i, beta in enumerate(self.escapeProb[u,l]):
                tau[i] = brentq(_betaTauSlab, 1e-10, 3.0/beta, \
                                    args=(beta,))
            return tau
        else:
            raise despoticError, "unknown escapeProbGeom " + \
                str(escapeProbGeom)


########################################################################
# Method to return the luminosity per H with all transitions
# (or with the specified transition if the transition keyword is
# specified) based on the stored escape probabilities. Note that the
# value returned is net luminosity, i.e. emission minus CMB
# absorption. Also note that this can be negative, indicating that the
# gas absorbs more energy than the CMB than it emits.
########################################################################
    def luminosityPerH(self, rad, transition=None, total=False, \
                           thin=False):
        """
        Return the luminosities of various lines, computed from the
        stored level populations and escape probabilities.

        Parameters
        ----------
        rad : class radiation
            radiation field impinging on the cloud
        transition : list of two arrays of shape (M)
            transition[0] = array of upper states
            transition[1] = array of lower states
        total : Boolean
            if True, the routine returns a single float rather than an
            array; this float is the sum of the luminosities per H
            nucleus of all lines

        Returns
        -------
        lum : array, shape (M), or float
              luminosities per H in specified lines; by default
              M = len(radUpper)

        Raises
        ------
        despoticError, if the level populations are not initialized,
        or if the escape probabilities are not initialized and thin is
        not True
        """

        # Safety check
        if self.levPopInitialized == False:
            raise despoticError, "Error: cannot compute luminosities before level populations are initialized; call one of the setLevPop routines first"

        # Get list of lines
        if transition==None:
            u = self.data.radUpper
            l = self.data.radLower
        else:
            u = transition[0]
            l = transition[1]

        # Get photon occupation number
        if rad.TCMB > 0.0:
            Tnu = self.data.freq[u,l]*h/kB
            Tnu[Tnu == 0.0] = rad.TCMB*1e20
            ngamma = rad.ngamma(Tnu)
        else:
            ngamma = 0.0

        # Get net luminosity
        lum = ((1.0+ngamma)*self.levPop[u] - \
                   self.data.wgtRatio[u,l]*ngamma*self.levPop[l]) * \
                   self.data.EinsteinA[u,l] * h*self.data.freq[u,l] * \
                   self.abundance
        if thin==False:
            lum *= self.escapeProb[u,l]
        if total==True:
            return np.sum(lum)
        else:
            return lum


########################################################################
# End of class emitterSpecies
########################################################################


########################################################################
# Helper function to return (1-exp(-tau))/tau, the functional form
# giving the escape probability vs. optical depth in the LVG
# approximation.
########################################################################
def _betaTauLVG(tau, beta):
    return beta-(1.0-np.exp(-tau))/tau

########################################################################
# Helper function to return (1-exp(-3 tau))/(3 tau), the functional form
# giving the escape probability vs. optical depth in the slab
# approximation.
########################################################################
def _betaTauSlab(tau, beta):
    return beta-(1.0-np.exp(-3.0*tau))/(3.0*tau)

