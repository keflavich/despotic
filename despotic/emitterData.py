"""
This module defines the class emitterData that stores information
about the physical properties of emitting species: energy levels,
Einstein coefficients, collision partners and their associated rate
coefficients, etc.
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

########################################################################
# Global variables
########################################################################

import numpy as np
import os
import datetime as dt
from copy import deepcopy
from .despoticError import despoticError
from .collPartner import collPartner
from .fetchLamda import fetchLamda
import warnings

# Define some global physical constants in cgs units
import scipy.constants as physcons
kB = physcons.k*1e7
c = physcons.c*1e2
mH = (physcons.m_p+physcons.m_e)*1e3
h = physcons.h*1e7
hc = h*c
G = physcons.G*1e3

# List of known collision partners, and their corresponding notations
# in the Leiden database
knownPart=['HI', 'pH2', 'oH2', 'He', 'e', 'H+']

########################################################################
# class emitterData
########################################################################
class emitterData(object):
    """
    Class to store the physical properties of a single emitting
    species, and preform computations on those properties. Note that
    all quantities are stored in cgs units.

    Parameters
       emitName : string
          name of this species
       emitterFile : string
          name of LAMDA file containing data on this species; this
          option overrides the default
       emitterURL : string
          URL of LAMDA file containing data on this species; this
          option overrides the default
       extrap : Boolean
          if True, collision rate coefficients for this species will
          be extrapolated to temperatures outside the range given in
          the LAMDA tables
       noRefresh : Boolean
          if True, the routine will not attempt to automatically
          fetch updated versions of files from the web

    Class attributes
       name : string
          name of emitting species
       lamdaFile : string
          name of file from which species was read
       molWgt : float
          molecular weight of species, in units of H masses
       nlev : int
          number of energy levels of the species
       levEnergy : array, shape(nlev)
          energies of levels
       levTemp : array, shape(nlev)
          energies of levels in K (i.e. levEnergy / kB)
       levWgt : array shape(nlev)
          degeneracies (statistical weights) of levels
       nrad : int
          number of radiative transition this species has
       radUpper : integer array, shape (nrad)
          array containing upper states for radiative transitions
       radLower : integer array, shape (nrad)
          array containing lower states for radiative transitions
       radFreq : array, shape (nrad)
          array of frequencies of radiative transitions
       radTemp : array, shape (nrad)
          same as radFreq, but multiplied by h/kB to give units of K
       radTUpper : array, shape (nrad)
          array of temperatures (E/kB) of upper radiative states
       radA : array, shape (nrad)
          array of Einstin A coefficients of radiative transitions
       partners : dict
          listing collision partners; keys are partner names, values are
          objects of class collPartner
       EinsteinA : array, shape (nlev, nlev)
          2d array of nlev x nlev giving Einstein A's for radiative
          transitions connecting each level pair
       freq : array, shape (nlev, nlev)
          2d array of nlev x nlev giving frequency of radiative
          transitions connecting each level pair
       dT : array, shape (nlev, nlev)
          2d array of nlev x nlev giving energy difference between each
          level pair in K; it is positive for i > j, and negative for i
          < j
       wgtRatio : array, shape (nlev, nlev)
          2d array giving ratio of statistical weights of each level
          pair
       extrap : Boolean
          if True, collision rate coefficients for this emitter are
          allowed to be extrapolated off the data table
    """

    ####################################################################
    # Initialization routine
    ####################################################################
    def __init__(self, emitName, emitterFile=None, emitterURL=None,
                 extrap=True, noRefresh=False, noWarn=False):
        """
        Initialization routine

        Parameters
           emitName : string
              name of this species
           emitterFile : string
              name of LAMDA file containing data on this species; this
              option overrides the default
           emitterURL : string
              URL of LAMDA file containing data on this species; this
              option overrides the default
           extrap : Boolean
              if True, collision rate coefficients for this species will
              be extrapolated to temperatures outside the range given in
              the LAMDA tables
           noRefresh : Boolean
              if True, the routine will not attempt to automatically
              fetch updated versions of files from the web
           noWarn : Boolean
              if True, suppress warning messages about missing
              collision rate data
        """

        # Save extrapolation and warning state
        self.__extrap = extrap
        self.noWarn = noWarn

        # Check if we have been given a file name for this emitter. If
        # we have, try to read it, and raise an error if we fail.
        if emitterFile != None:
            try:
                fp = open(emitterFile, 'r')
                self.readLamda(fp, extrap=extrap)
                fp.close()
                self.lamdaFile = emitterFile     # Store filename
                return
            except IOError:
                # If we're here, try looking in the sub-directory
                # LAMDA of $DESPOTIC_HOME
                if 'DESPOTIC_HOME' in os.environ:
                    # Use environment variable if available
                    emitterPath = \
                        os.path.join(os.environ['DESPOTIC_HOME'], \
                                         'LAMDA')
                else:
                    # No path given, so just use relative location
                    emitterPath = 'LAMDA'
                # Try with new path
                try:
                    fileName = os.path.join(emitterPath, emitterFile)
                    fp = open(fileName, 'r')
                    self.readLamda(fp, extrap=extrap)
                    fp.close()
                    self.lamdaFile = fileName     # Store filename
                    return
                except IOError:
                    raise despoticError(
                        "unable to open file "+emitterFile)

        # Set name
        self.name = emitName

        # If we're here, we weren't given a file name. Check if we
        # have been given a URL for this emitter. If we have, try to
        # read it, and raise an error if we fail. If we downloaded a
        # file successfully, read it and return
        if emitterURL != None:
            emitterFile = fetchLamda(emitterURL)
            if emitterFile == None:
                raise despoticError(
                    "unable to download file "+emitterURL)
            else:
                try:
                    fp = open(emitterFile, 'r')
                    self.readLamda(fp, extrap=extrap)
                    fp.close()
                    self.lamdaFile = emitterFile     # Store filename
                    return
                except IOError:
                    raise despoticError(
                        "unable to open file "+emitterFile)

        # If we're here, were weren't given an emitter file name or
        # URL, so we have to look for the data ourselves.

        # Step 1: build possible file names. For a given emitter
        # species, there are several possible names in LAMDA, since
        # in some cases only certain forms of the data are available
        # (for example only extrapolated data exist for some species)
        emitterFile = []
        emitterFile.append(emitName+'.dat')      # Usual extension
        emitterFile.append(emitName+'@xpol.dat') # Extrapolation flag
        emitterFile.append(emitName+'atom.dat')  # Atom flag
        if emitName.lower() != emitName:         # Lowercase versions
            namesCopy = deepcopy(emitterFile)
            for name in namesCopy:
                emitterFile.append(name.lower())

        # Step 2: set path to data files
        if 'DESPOTIC_HOME' in os.environ:
            # Use environment variable if available
            emitterPath = \
                os.path.join(os.environ['DESPOTIC_HOME'], \
                                 'LAMDA')
        else:
            # No path given, so just use relative location
            emitterPath = 'LAMDA'

        # Step 3: try to open and read a local copy of LAMDA file,
        # using the various possible file names; return if we are
        # successful. Also, refresh old files in the process.
        for fname in emitterFile:
            fileName = os.path.join(emitterPath, fname)
            try:
                fp = open(fileName, 'r')
                # If we're here, we've found and opened the
                # file. Before reading, though, make sure it's not out
                # of date. If it is, try to grab a new version off the
                # web unless instructed otherwise.
                if noRefresh == False:
                    mtime = dt.datetime.\
                        fromtimestamp(os.stat(fileName).st_mtime)
                    if mtime < dt.datetime.today()-dt.timedelta(182.5):
                        fp.close()
                        fetchLamda(os.path.basename(fileName), \
                                       path='', \
                                       fileName=fileName)
                        fp = open(fileName, 'r')
                # Now read
                self.readLamda(fp, extrap=extrap)
                fp.close()
                self.lamdaFile = fileName        # Store filename
                return
            except IOError:
                pass

        # Step 4: if we're here, we failed to find a local copy of the
        # LAMDA file, so look for it on the web; if we succeed, read
        # the resulting file and return
        for fname in emitterFile:
            downloadName = fetchLamda(fname)
            if downloadName != None:
                try:
                    fp = open(downloadName, 'r')
                    self.readLamda(fp, extrap=extrap)
                    fp.close()
                    self.lamdaFile = downloadName    # Store filename
                    return
                except IOError:
                    raise despoticError(
                        "downloaded file " + 
                        downloadName + " from LAMDA, but " + 
                        "cannot open it")

        # Step 5: if we're here, we've failed to find the file on the
        # web, so give up and raise an error
        raise despoticError(
            'cannot locate LAMDA file for emitter ' + emitName + 
            '; try setting URL or file name by hand')



    ####################################################################
    # Routine to read LAMDA data files
    ####################################################################
    def readLamda(self, fp, extrap=False):
        """
        Read a LAMDA-formatted file

        Parameters
           fp : file
              pointer to the start of a LAMDA-formatted file
           extrap : Boolean
              if True, collision rate coefficients for this species will
              be extrapolated to temperatures outside the range given in
              the LAMDA tables

        Returns
           Nothing
        """

        # Read header information: molecule name and weight
        line = ''
        while (line.strip() == ''):
            line = fp.readline()
            if line.strip()[0] == '!':
                line = ''
        line = ''
        while (line.strip() == ''):
            line = fp.readline()
            if line.strip()[0] == '!':
                line = ''
        self.molWgt = float(line)

        # Now read level data
        line = ''
        while (line.strip() == ''):
            line = fp.readline()
            if line.strip()[0] == '!':
                line = ''
        self.nlev = int(line)
        self.levEnergy = np.zeros(self.nlev)
        self.levWgt = np.zeros(self.nlev)
        for i in range(self.nlev):
            line = ''
            while (line.strip() == ''):
                line = fp.readline()
                if line.strip()[0] == '!':
                    line = ''
            linesplit = line.split()
            self.levEnergy[i] = float(linesplit[1])*hc
            self.levWgt[i] = float(linesplit[2])
        self.levTemp = self.levEnergy / kB

        # Read radiative transition data
        line = ''
        while (line.strip() == ''):
            line = fp.readline()
            if line.strip()[0] == '!':
                line = ''
        self.nrad = int(line)
        self.radUpper = np.zeros(self.nrad, dtype='int')
        self.radLower = np.zeros(self.nrad, dtype='int')
        self.radFreq = np.zeros(self.nrad)
        self.radTUpper = np.zeros(self.nrad)
        self.radA = np.zeros(self.nrad)
        for i in range(self.nrad):
            line = ''
            while (line.strip() == ''):
                line = fp.readline()
                if line.strip()[0] == '!':
                    line = ''
            linesplit = line.split()
            self.radUpper[i] = int(linesplit[1]) - 1   # Correct to 0 offset
            self.radLower[i] = int(linesplit[2]) - 1   # Correct to 0 offset
            self.radA[i] = float(linesplit[3])
            self.radFreq[i] = float(linesplit[4])*1e9
            self.radTUpper[i] = float(linesplit[5])

        # Read collision partners
        line = ''
        while (line.strip() == ''):
            line = fp.readline()
            if line.strip()[0] == '!':
                line = ''
        npart=int(line)
        self.partners={}
        for i in range(npart):
            # First read name
            line = ''
            while (line.strip() == ''):
                line = fp.readline()
                if line.strip()[0] == '!':
                    line = ''
            idstring = line.strip().split()[0]
            # Assign particle ID based on LAMDA codes
            if idstring=='1':
                # Generic H2; read into para, store in ortho below
                partName = 'pH2'
            elif idstring=='2':
                partName = 'pH2'
            elif idstring=='3':
                partName = 'oH2'
            elif idstring=='4':
                partName = 'e'
            elif idstring=='5':
                partName = 'HI'
            elif idstring=='6':
                partName = 'He'
            elif idstring=='7':
                partName = 'H+'
            else:
                raise despoticError(
                    "unknown collision partner ID " + idstring)

            self.partners[partName] = \
                collPartner(fp, self.nlev, extrap=extrap)
            # For code 1, set ortho-H2 rates equal to para-H2 rates
            if idstring=='1':
                self.partners['oH2'] = self.partners['pH2']

        # Close file
        fp.close()

        # Do we have values of only ortho-H2 or para-H2? If so, assume
        # one is equal to the other.
        if ('oH2' in self.partners) and not ('pH2' in self.partners):
            self.partners['pH2'] = self.partners['oH2']
        if ('pH2' in self.partners) and not ('oH2' in self.partners):
            self.partners['oH2'] = self.partners['pH2']

        # Generate 2d arrays containing Einstein coefficients, frequencies,
        # and degeneracy ratios. These are useful for computation.
        self.EinsteinA = np.zeros((self.nlev, self.nlev))
        self.EinsteinA[self.radUpper, self.radLower] = self.radA
        self.freq = np.zeros((self.nlev, self.nlev))
        self.freq[self.radUpper, self.radLower] = self.radFreq
        self.freq[self.radLower, self.radUpper] = self.radFreq
        self.radTemp = (h/kB) * self.radFreq
        self.dT = np.zeros((self.nlev, self.nlev))
        self.dT[self.radUpper, self.radLower] = (h/kB) * self.radFreq
        self.dT[self.radLower, self.radUpper] = -(h/kB) * self.radFreq
        self.wgtRatio = np.zeros((self.nlev, self.nlev))
        for i in range(self.nlev):
            self.wgtRatio[i,:] = self.levWgt[i]/self.levWgt


    ####################################################################
    # Extrapolation decorator
    ####################################################################
    @property
    def extrap(self):
        return self.__extrap

    @extrap.setter
    def extrap(self, ext):
        # Change value for all collision partners
        for p in self.partners.values():
            p.extrap = ext
        self.__extrap = ext

    ####################################################################
    # Method to compute partition function
    ####################################################################
    def partFunc(self, temp):
        """
        Compute the partition function for this species at the given
        temperature.

        Parameters
           temp : float | array
              gas kinetic temperature

        Returns
           Z : float | array
              the partition function Z(T) for this species
        """

        # If temp is an list, return an array. If it's a number, just
        # return a number.
        if hasattr(temp, '__iter__'):
            return np.inner(self.levWgt, 
                            np.exp(-np.outer(1.0/temp, self.levTemp)))
        else:
            return sum(self.levWgt*np.exp(-self.levTemp/temp))


    ####################################################################
    # Method to return the collision rate matrix for this emitter as a
    # function of density, composition, and temperature
    ####################################################################
    def collRateMatrix(self, nH, comp, temp):
        """
        This routine computes the matrix of collision rates (not rate
        coefficients) between every pair of levels.

        Parameters
           nH : float
              number density of H nuclei
           comp : class composition
              bulk composition of the gas
           temp : float
              gas kinetic temperature

        Returns
           q : array, shape (nlev, nlev)
              array in which element ij is the rate of collisional
              transitions from state i to state j, in s^-1
        """

        # ensure that there are sane abundances defined for colliders, else the
        # solver will fail with a rather opaque message
        comp._check_abundance()

        # Initialize collision rate matrix
        q = np.zeros((self.nlev, self.nlev))

        # Loop over collision partners
        for p in knownPart:

            # Set appropriate density of collision parnter, and make
            # sure that we have collision partner data that we need
            if p=='HI':
                if comp.xHI==0:
                    continue
                else:
                    n = nH * comp.xHI
            elif p=='pH2':
                if comp.xpH2==0:
                    continue
                else:
                    n = nH * comp.xpH2
            elif p=='oH2':
                if comp.xoH2==0:
                    continue
                else:
                    n = nH * comp.xoH2
            elif p=='He':
                if comp.xHe==0:
                    continue
                else:
                    n = nH * comp.xHe
                    if not 'He' in self.partners:
                        # If He is not available, assume it is the
                        # same as pH2 divided by the square root of
                        # the ratio of the reduced mass for this
                        # molecule plus He to the reduced mass for
                        # this molecule plus H_2
                        mredHe = self.molWgt*4.0 / (self.molWgt + 4.0)
                        mredH2 = self.molWgt*2.0 / (self.molWgt + 2.0)
                        n = n / np.sqrt(mredHe/mredH2)
            elif p=='e':
                if comp.xe==0:
                    continue
                else:
                    n = nH * comp.xe
            elif p=='H+':
                if comp.xHplus==0:
                    continue
                else:
                    n = nH * comp.xHplus

            # Compute downward transition rates, handling He as a
            # proxy for H2 as a special case
            if (p != 'He') or ('He' in self.partners):
                try:
                    q += n * self.partners[p].\
                        colRateMatrix(temp, self.levWgt, self.levTemp)
                except (ValueError, despoticError):
                    raise despoticError(
                        "Temperature T = "+str(temp) + 
                        " K is outside tabulated range " + 
                        str(self.partners[p].tempTable[0])+" - " + 
                        str(self.partners[p].tempTable[-1]) + 
                        " K for species "+self.name+", partner " + 
                        p)
                except KeyError:
                    if not self.noWarn:
                        warnings.warn("collision rates not "
                                      "available for "+self.name+
                                      " with species "
                                      +p+" (fractional abundance "
                                      "{:e}".format(n/nH)+") will be "
                                      "omitted from calculation")
            else:
                p1='pH2'
                try:
                    q += n * self.partners[p1].\
                        colRateMatrix(temp, self.levWgt, self.levTemp)
                except (ValueError, despoticError):
                    raise despoticError("Temperature T = "+str(temp) + 
                        " K is outside tabulated range " + 
                        str(self.partners[p1].tempTable[0])+" - " + 
                        str(self.partners[p1].tempTable[-1]) + 
                        " K for species "+self.name+", partner " + 
                        p)
                except KeyError:
                    if not self.noWarn:
                        warnings.warn("collision rates not "
                                      "available for "+self.name+
                                      " with species "
                                      +p+" (fractional abundance "
                                      "{:e}".format(n/nH)+") will be "
                                      "omitted from calculation")

        # Return matrix
        return q

    ####################################################################
    # Method to derive the tX parameter defined by Krumholz+ (2016)
    ####################################################################
    def tX(self, mX, trans=None):
        """
        Returns the tX line strength parameter of Krumholz+ (2016)

        Parameters
           mX : float
              total mass per particle of this species, in g
           trans : int, array, or None
              if set, tX is returned only for the specified
              transitions in the transition list; default is that it
              is returned for all transitions

        Returns
           tX : float or array
              transition strength parameter for the specified
              transitions, in seconds
        """
        if trans is None:
            gi = self.levWgt[self.radUpper]
            gj = self.levWgt[self.radLower]
            lam = c/self.radFreq
            A = self.radA
        else:
            gi = self.levWgt[self.radUpper][trans]
            gj = self.levWgt[self.radLower][trans]
            lam = c/self.radFreq[trans]
            A = self.radA[trans]
        tX = (gi/gj) * A*lam**3 / (32*np.pi**2*G*mX)
        return tX

    ####################################################################
    # Methods to derive the Xthin and alphathin parameters defined by
    # Krumholz+ (2016)
    ####################################################################
    def Xthin(self, abd, trans=None):
        """
        Returns the Xthin parameter of Krumholz+ (2016)

        Parameters
           abd : float
              abundance of the species, relative to H
           trans : int, array, or None
              if set, Xthin is returned only for the specified
              transitions in the transition list; default is that it
              is returned for all transitions

        Returns
           Xthin : float or array
              Xthin parameter for the specified transitions, in cm^-2
              / (K km s^-1)
        """
        if trans is None:
            lam = c/self.radFreq
            A = self.radA
            T = h*self.radFreq/kB
        else:
            lam = c/self.radFreq[trans]
            A = self.radA[trans]
            T = h*self.radFreq[trans]/kB
        Xthin = 8.0*np.pi/(A*T*lam**3*abd)*1e5   # 1e5 is to convert to km/s
        return Xthin

    def alphathin(self, mX, trans=None):
        """
        Returns the Xthin parameter of Krumholz+ (2016)

        Parameters
           mX : float
              total mass per particle of this species, in g
           trans : int, array, or None
              if set, Xthin is returned only for the specified
              transitions in the transition list; default is that it
              is returned for all transitions

        Returns
           Xthin : float or array
              Xthin parameter for the specified transitions, in cm^-2
              / (K km s^-1)
        """
        if trans is None:
            lam = c/self.radFreq
            A = self.radA
            T = h*self.radFreq/kB
        else:
            lam = c/self.radFreq[trans]
            A = self.radA[trans]
            T = h*self.radFreq[trans]/kB
        Xthin = 8.0*np.pi*mX/(A*T*lam**3)*1e5   # 1e5 is to convert to km/s
        return Xthin
    

########################################################################
# End of class emitterData
########################################################################




