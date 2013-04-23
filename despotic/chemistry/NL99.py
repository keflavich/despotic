"""
This module implements the reduced carbon-oxygen chemistry network of
Nelson & Langer (1999, ApJ, 524, 923)
"""

import numpy as np
import string
from despotic.despoticError import despoticError
from shielding import fShield_CO_vDB
import scipy.constants as physcons

########################################################################
# Physical constants
########################################################################
kB = physcons.k/physcons.erg
mH = (physcons.m_p+physcons.m_e)/physcons.gram

########################################################################
# Define a mapping from species name to numerical index.
########################################################################
specDict = { 'He+' : 0, 
              'H3+' : 1, 
              'OHx' : 2, 
              'CHx' : 3, 
              'CO' : 4, 
              'C' : 5, 
              'C+' : 6, 
              'HCO+' : 7, 
              'O' : 8, 
              'M+' : 9 }
specList = ['He+', 'H3+', 'OHx', 'CHx', 'CO', 'C', 'C+', 'HCO+', 'O',
            'M+']


########################################################################
# Data on photoreactions
# Reactions are, in order:
# h nu + CI -> C+ + e
# h nu + CHx -> CI + H
# h nu + CO -> CI + O
# h nu + OHx -> OI + H
# h nu + M -> M+ + e
# h nu + HCO+ -> CO + H
########################################################################
_kph = np.array([ 
        3.0e-10, 1.0e-9, 1.0e-10, 5.0e-10, 
        2.0e-10, 1.5e-10])
_avfac = np.array([3.0, 1.5, 3.0, 1.7, 1.9, 2.5])
_inph = np.array([5, 3, 4, 2, 12, 7], dtype='int')
_outph1 = np.array([6, 5, 5, 8, 9, 4], dtype='int')
_outph2 = np.array([10, 10, 8, 10, 10, 10], dtype='int')


########################################################################
# Data on two-body reactions
# Reactions are, in order:
# (0) H3+ + CI -> CHx + H2
# (1) H3+ + OI -> OHx + H2
# (2) H3+ + CO -> HCO+ + H2
# (3) He+ + H2 -> He + H + H+
# (4) He+ + CO -> C+ + O + He
# (5) C+ + H2 -> CHx + H
# (6) C+ + OHx -> HCO+
# (7) OI + CHx -> CO + H
# (8) CI + OHx -> CO + H
# (9) He+ + e -> He + h nu
# (10) H3+ + e -> H2 + H
# (11) C+ + e -> CI + h nu
# (12) HCO+ + e -> CO + H
# (13) M+ + e -> M + h nu
# (14) H3+ + M -> M+ + H + H2
########################################################################
_k2 = np.array([ 
        2.0e-9, 8.0e-10, 1.7e-9, 7.0e-15, 1.6e-9, 4.0e-16, 1.0e-9, 
        2.0e-10, 5.8e-12, 9.0e-11, 1.9e-6, 1.4e-10, 3.3e-5, 
        3.8e-10, 2.0e-9])
_k2Texp = np.array([ 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.64, -0.54, 
        -0.61, -1.0, -0.65, 0.0])
_in2bdy1 = np.array([1, 1, 1, 0, 0, 6, 6, 8, 5, 0, 1, 6, 7, 9, 1], 
                    dtype='int') 
_in2bdy2 = np.array([5, 8, 4, 10, 4, 10, 2, 3, 2, 13, 13, 13, 13, 
                     13, 12], dtype='int')
_out2bdy1 = np.array([3, 2, 7, 10, 6, 3, 7, 4, 4, 10, 10, 5, 4, 10, 
                      9], dtype='int')
_out2bdy2 = np.array([10, 10, 10, 10, 8, 10, 10, 10, 10, 10, 10, 
                      10, 10, 10, 10], dtype='int')



########################################################################
# Set some default abundances
########################################################################
_xCdefault = 2.0e-4
_xOdefault = 4.0e-4
_xMdefault = 2.0e-7
_xH2 = 0.5

########################################################################
# Define the NL99 class
########################################################################
class NL99:
    """
    This class the implements the chemistry network of Nelson & Langer
    (1999, ApJ, 524, 923).
    """

########################################################################
# Method to initialize
########################################################################
    def __init__(self, cloud=None, info=None):
        """
        Parameters
        ----------
        cloud : class cloud
             a DESPOTIC cloud object from which initial data are to be
             taken
        info : dict
             a dict containing additional parameters

        Returns
        -------
        Nothing

        Raises
        ------
        despoticError, if the dict info contains non-allowed entries

        Remarks
        -------
        The dict info may contain the following keys:

        'xC' : total C abundance per H nucleus
        'xO' : total H abundance per H nucleus
        'xM' : total refractory metal abundance per H nucleus
        'sigmaDustV' : V band dust extinction cross section per H
        nucleus
        'AV' : total visual extinction; mutually exclusive with
        sigmaDustV
        'noClump' : a Boolean; if True, the clump factor is set to 1.0

        Default values are xC = 2e-4, xO = 4e-4, xM = 2e-7, sigmaDustV
        = 0.4 * cloud.dust.sigmaPE
        """

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
            self.xHe = 0.1
            self.ionRate = 2.0e-17
            self.NH = 0.0
            self.temp = 0.0
            self.chi = 1.0
            self.nH = 0.0
            self.AV = 0.0

            # Set initial abundances
            self.emList = None
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
            if gmc.XH2 != 0.5:
                raise despoticError, "NL99 network only valid " + \
                    "for pure H2 composition"

            # Physical properties
            self.xHe = cloud.comp.xHe
            self.ionRate = cloud.rad.ionRate
            self.NH = cloud.colDen / 2.0
            self.temp = cloud.Tg
            self.chi = cloud.rad.chi
            if info is None:
                cs2 = kB * cloud.Tg / (cloud.comp.mu * mH)
                cfac = np.sqrt(1.0 + 0.75*cloud.sigmaNT**2/cs2)
            elif 'noClump' in info:
                if info['noClump'] == True:
                    cfac = 1.0
                else:
                    cs2 = kB * cloud.Tg / (cloud.comp.mu * mH)
                    cfac = np.sqrt(1.0 + 0.75*cloud.sigmaNT**2/cs2)
            else:
                cs2 = kB * cloud.Tg / (cloud.comp.mu * mH)
                cfac = np.sqrt(1.0 + 0.75*cloud.sigmaNT**2/cs2)
            self.nH = cloud.nH * cfac
            if info is None:
                self.AV = 0.4 * cloud.dust.sigmaPE * self.NH
            elif 'sigmaDustV' in info:
                # Note factor to convert from mag to true
                # dimensionless units
                self.AV = self.NH * info['sigmaDustV'] / \
                    np.log(100**0.2)
            else:
                self.AV = 0.4 * cloud.dust.sigmaPE * self.NH

            # Set abundances

            # Make a case-insensitive version of the emitter list for
            # convenience
            self.emList = dict(zip(map(string.lower, 
                                       cloud.emitters.keys()), 
                                   cloud.emitters.values()))

            # OH and H2O
            if 'oh' in self.emList:
                self.x[2] += self.emList['oh'].abundance
            if 'ph2o' in self.emList:
                self.x[2] += self.emList['ph2o'].abundance
            if 'oh2o' in self.emList:
                self.x[2] += self.emList['oh2o'].abundance
            if 'p-h2o' in self.emList:
                self.x[2] += self.emList['p-h2o'].abundance
            if 'o-h2o' in self.emList:
                self.x[2] += self.emList['o-h2o'].abundance

            # CO
            if 'co' in self.emList:
                self.x[4] = self.emList['co'].abundance

            # Neutral carbon
            if 'c' in self.emList:
                self.x[5] = self.emList['c'].abundance

            # Ionized carbon
            if 'c+' in self.emList:
                self.x[6] = self.emList['c+'].abundance

            # HCO+
            if 'hco+' in self.emList:
                self.x[7] = self.emList['hco+'].abundance

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
                if self.x[6] != 0.0:
                    print "Warning: input C abundance is " + \
                        str(xC) + ", but total input C, C+, CO, " + \
                        "HCO+ abundance is " + str(xCtot) + \
                        "; increasing xC+ to " + str(xC-xCtot)
                self.x[6] = xC - xCtot
            elif xCtot > xC:
                # Throw an error if input C abundance is smaller than
                # what is accounted for in initial conditions
                raise despoticError, "input C abundance is " + \
                    str(xC) + ", but total input C, C+, CO, " + \
                    "HCO+ abundance is " + str(xCtot)

            # O
            if 'o' in self.emList:
                self.x[8] = self.emList['o'].abundance
            elif info is None:
                self.x[8] = _xOdefault
            elif 'xO' in info:
                self.x[8] = info['xO']
            else:
                self.x[8] = _xOdefault

            # Ratios of ^12C to ^13C, and ^16O to ^18O
            if '13co' in self.emList and 'co' in self.emList:
                self.c13_12 = self.emList['13co'].abundance / \
                    self.emList['co'].abundance
            if 'c18o' in self.emList and 'co' in self.emList:
                self.o18_16 = self.emList['c18o'].abundance / \
                    self.emList['co'].abundance

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


########################################################################
# Method to return the time derivative of all chemical rates
########################################################################
    def dxdt(self, xin, time):
        """
        This method returns the time derivative of all abundances for
        this chemical network.

        Parameters
        ----------
        xin : array(10)
             current abundances of all species
        time : float
             current time; not actually used, but included as an
             argument for compatibility with odeint

        Returns
        -------
        dxdt : array(10)
             time derivative of x
        """

        # Vector to store results; it is convenient for this to have
        # some phantom slots; slots 10, 11, 12, and 13 store
        # abundances of H2, HeI, MI, and e, respectively
        xdot = np.zeros(14)
        xgrow = np.zeros(14)
        xgrow[:10] = xin
        xgrow[10] = _xH2
        xgrow[11] = self.xHe - xin[0]
        xgrow[12] = self.xM - xin[9]
        xgrow[13] = xin[0] + xin[1] + xin[6] + xin[7] + xin[9]

        # Cosmic ray / x-ray ionization reactions
        xdot[0] = xgrow[11]*self.ionRate
        xdot[1] = self.ionRate

        # Photon reactions
        ratecoef = self.chi/1.6*np.exp(-_avfac*self.AV)*_kph
        rate = ratecoef*xgrow[_inph]
        # Apply CO line shielding factor
        rate[2] = rate[2] * fShield_CO_vDB(xgrow[4]*self.NH, self.NH/2.0) 
        for i, n in enumerate(_inph):
            xdot[_inph[i]] -= rate[i]
            xdot[_outph1[i]] += rate[i]
            xdot[_outph2[i]] += rate[i]

        # Two-body reactions
        rate = _k2*self.temp**_k2Texp*self.nH*xgrow[_in2bdy1]*xgrow[_in2bdy2]
        for i, n in enumerate(_in2bdy1):
            xdot[_in2bdy1[i]] -= rate[i]
            xdot[_in2bdy2[i]] -= rate[i]
            xdot[_out2bdy1[i]] += rate[i]
            xdot[_out2bdy2[i]] += rate[i]

        # Return results
        return xdot[:10]


########################################################################
# Method to take an array of abundances and turn them into a dict of
# abundances for easier human reading
########################################################################
    def xToDict(self, abundances=None):
        """
        This method turns a numerical array of abundances into a dict
        of abundances, for easier human reading
        
        Parameters
        ----------
        abundances : array or None
             abundances to be transformed

        Returns
        -------
        xDict : dict
             dict of abundances

        Raises
        ------
        despoticError if abundances is not None or a 1D or 2D array
        with final dimension of length 10

        Remarks
        -------
        If abundances is None, then xDict is a dict where each entry
        describes the abundance of the corresponding element in
        self.x. If abundances is a 1D array, xDict is a dict where
        each entry gives the abundance of the corresponding element in
        the input abundances array. If abundances is a 2D array, xDict
        is a dict where each entry is an array whose length is equal
        to the length of the first dimension of the array.
        """

        # Prepare output dict
        xDict = {}

        # Proceed based on input abundances
        if abundances==None:

            # None, so return data from stored abundances
            for i, s in enumerate(specList):
                xDict[s] = x[i]

        else:

            # Abundances have been input
            if rank(abundances) < 1 or rank(abundances) > 2:

                # Throw error: input array must be 1d or 2d
                raise despoticError, "abundances must be " + \
                    "array of rank 1 or 2"

            elif rank(abundances) == 1:

                # 1D input; check length
                if abundances.shape[-1] != 10:
                    raise despoticError, "last dimension of " + \
                        "abundances must be of size 10"

                # Construct output
                for i, s in enumerate(specList):
                    xDict[s] = abundances[i]

            else:

                # 2D input; check length of last dimension
                if abundances.shape[-1] != 10:
                    raise despoticError, "last dimension of " + \
                        "abundances must be of size 10"

                # Construct output
                for i, s in enumerate(specList):
                    xDict[s] = abundances[:,i]

        # Return output
        return xDict


########################################################################
# Define abundances as a property, which just gets the current values
# in the x array and returns them in a nicely-formatted dict
########################################################################
    @property
    def abundances(self):
        if hasattr(self, '_abundances'):
            return self._abundances

    @abundances.setter
    def abundances(self, 

########################################################################
# Method to take an dict of abundances and turn them into an array in
# the ordering used by this chemistry network
########################################################################
    def xFromDict(self, abundances=None):
        """
        This method takes 
        """


########################################################################
# Method to write the currently stored abundances to a cloud
########################################################################
    def applyAbundances(abundances=None, cloud=None):
        """
        This method writes abundances to a cloud's emitter list

        Parameters
        ----------
        abundances : None or array(10)
             array of abundances
