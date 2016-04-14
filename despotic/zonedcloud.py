"""
This module defines the zonedcloud class. A zonedcloud is a collection
of cloud objects at different extinctions and densities, which are
taken to be exposed to the same external radiation field.
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
from .cloud import cloud
from .dustProp import dustProp
from .despoticError import despoticError
from copy import deepcopy
import scipy.constants as physcons
mH = physcons.m_p*1e3
kB = physcons.k*1e7
G = physcons.G*1e3
h = physcons.h*1e7
c = physcons.c*1e2

class zonedcloud(object):
    """
    A class consisting of an interstellar cloud divided into different
    column density / extinction zones.


    Parameters
       fileName : string
          name of file from which to read cloud description
       colDen : array
          Array of column densities marking zone centers
       AV : array
          Array of visual extinction values (in mag) marking zone
          centers; AV is converted to column density using a V-band
          cross section equal to 0.4 * sigmaPE; this argument
          ignored if colDen is not None
       nZone : int
          Number of zones into which to divide the cloud, from 0
          to the maximum column density found in the cloud
          description file fileName; ignored if colDen or AV is
          not None
       geometry : 'sphere' | 'slab'
          geometry to assume for the cloud, either 'sphere'
          (onion-like) or 'slab' (layer cake-like)
       verbose : Boolean
          print out information about the cloud as we read it
    """

    ####################################################################
    # Method to initialize
    ####################################################################
    def __init__(self, fileName=None, colDen=None, AV=None, nZone=16, 
                 geometry='sphere', verbose=False):
        """
        This creates a zoned cloud, with zones at different column
        densities. The user must set one of colDen, AV, or fileName.

        Parameters
           fileName : string
              name of file from which to read cloud description
           colDen : array
              Array of column densities marking zone centers
           AV : array
              Array of visual extinction values (in mag) marking zone
              centers; AV is converted to column density using a V-band
              cross section equal to 0.4 * sigmaPE; this argument
              ignored if colDen is not None
           nZone : int
              Number of zones into which to divide the cloud, from 0
              to the maximum column density found in the cloud
              description file fileName; ignored if colDen or AV is
              not None
           geometry : 'sphere' | 'slab'
              geometry to assume for the cloud, either 'sphere'
              (onion-like) or 'slab' (layer cake-like)
           verbose : Boolean
              print out information about the cloud as we read it

        Returns
           Nothing
        """

        # Make sure the user has given us column density, AV, or file
        # name
        if (colDen is None) and (AV is None) and (fileName is None):
            raise despoticError(
                "zonedcloud: must specify one of: " + 
                "colDen, AV, or fileName")

        # Create a base cloud
        basecloud = cloud(fileName=fileName, verbose=verbose)

        # Get column densities to all zones
        if colDen is not None:

            # Column density explicitly specified
            self._colDen = np.zeros(len(colDen)+1)
            self._colDen[1:-1] = 0.5*(colDen[1:]+colDen[:-1])
            self._colDen[-1] = 2*colDen[-1]-colDen[-2]

        elif AV is not None:

            # Get column density from AV
            colDen = AV / (0.4*basecloud.dust.sigmaPE)
            self._colDen = np.zeros(len(colDen)+1)
            self._colDen[1:-1] = 0.5*(colDen[1:]+colDen[:-1])
            self._colDen[-1] = 2*colDen[-1]-colDen[-2]

        else:

            # Get column density from column density of base cloud,
            # divided by number of zones
            self._colDen = basecloud.colDen / nZone * (np.arange(nZone+1))

        # Make sure column densities are sorted, positive, and start
        # with 0
        self._colDen.sort()
        if self._colDen[0] < 0:
            raise despoticError(
                "zonedcloud: column densities must be > 0!")
        if self._colDen[0] > 0:
            np.insert(self._colDen, 0, 0.0)

        # Create list of clouds
        self.zones = [deepcopy(basecloud) for n in
                      range(len(self._colDen)-1)]

        # Set column densities of zones
        self.colDen = 0.5*(self._colDen[:-1] + self._colDen[1:])

        # Force dust and radiation field properties to be the same in
        # each cloud; that way if one is changed all will be changed
        for z in self.zones:
            z.dust = self.zones[0].dust
            z.rad = self.zones[0].rad

        # Store the geometry
        self.geometry = geometry

    ####################################################################
    # Properties that allow access to attributes of zones. Each of
    # these properties matches an attribute of the cloud
    # class. Depending on the attribute, the properties can be used to
    # get/set the attributes as a single value, a list, or an array.
    ####################################################################

    @property
    def nH(self):
        return np.array([z.nH for z in self.zones])

    @nH.setter
    def nH(self, nH):
        # nH can be a float or listlike
        if hasattr(nH, '__iter__'):
            if len(nH) != len(self.zones):
                raise despoticError(
                    "zonedcloud: could not " + 
                    "broadcast "+str(len(nH))+" items to " + 
                    str(len(self.zones)) + " zones")
            for z, x in zip(self.zones, nH):
                z.nH = x
        else:
            for z in self.zones:
                z.nH = nH

    @property
    def colDen(self):
        return np.array([z.colDen for z in self.zones])

    @colDen.setter
    def colDen(self, colDen):

        # colDen must be listlike, in increasing order
        if len(colDen) != len(self.zones):
            raise despoticError(
                "zonedcloud: could not " + 
                "broadcast "+str(len(colDen))+" items to " + 
                str(len(self.zones)) + " zones")
        for i in range(len(colDen)-1):
            if (colDen[i] >= colDen[i+1]) or (colDen[i] <= 0):
                raise despoticError(
                    "zonedcloud: column " + 
                    "densities must be positive and " + 
                    "monotonically increasing")

        # Set new column densities
        for z, x in zip(self.zones, colDen):
            z.colDen = x

        # Get new edge column densities
        self._colDen[1:-1] = 0.5*(colDen[1:]+colDen[:-1])
        self._colDen[-1] = 2*colDen[-1]-colDen[-2]

    @property
    def sigmaNT(self):
        return np.array([z.sigmaNT for z in self.zones])

    @sigmaNT.setter
    def sigmaNT(self, sigmaNT):
        # sigmaNT can be a float or listlike
        if hasattr(sigmaNT, '__iter__'):
            if len(sigmaNT) != len(self.zones):
                raise despoticError(
                    "zonedcloud: could not " + 
                    "broadcast "+str(len(sigmaNT))+" items to " + 
                    str(len(self.zones)) + " zones")
            for z, x in zip(self.zones, sigmaNT):
                z.sigmaNT = x
        else:
            for z in self.zones:
                z.sigmaNT = sigmaNT

    @property
    def dVdr(self):
        return np.array([z.dVdr for z in self.zones])

    @dVdr.setter
    def dVdr(self, dVdr):
        # dVdr can be a float or listlike
        if hasattr(dVdr, '__iter__'):
            if len(dVdr) != len(self.zones):
                raise despoticError(
                    "zonedcloud: could not " + 
                    "broadcast "+str(len(dVdr))+" items to " + 
                    str(len(self.zones)) + " zones")
            for z, x in zip(self.zones, dVdr):
                z.dVdr = x
        else:
            for z in self.zones:
                z.dVdr = dVdr

    @property
    def Tg(self):
        return np.array([z.Tg for z in self.zones])

    @Tg.setter
    def Tg(self, Tg):
        # Tg can be a float or listlike
        if hasattr(Tg, '__iter__'):
            if len(Tg) != len(self.zones):
                raise despoticError(
                    "zonedcloud: could not " + 
                    "broadcast "+str(len(Tg))+" items to " + 
                    str(len(self.zones)) + " zones")
            for z, x in zip(self.zones, Tg):
                z.Tg = x
        else:
            for z in self.zones:
                z.Tg = Tg

    @property
    def Td(self):
        return np.array([z.Td for z in self.zones])

    @Td.setter
    def Td(self, Td):
        # Td can be a float or listlike
        if hasattr(Td, '__iter__'):
            if len(Td) != len(self.zones):
                raise despoticError(
                    "zonedcloud: could not " + 
                    "broadcast "+str(len(Td))+" items to " + 
                    str(len(self.zones)) + " zones")
            for z, x in zip(self.zones, Td):
                z.Td = x
        else:
            for z in self.zones:
                z.Td = Td

    @property
    def comp(self):
        return [z.comp for z in self.zones]

    @comp.setter
    def comp(self, comp):
        # comp can be a single item or or listlike; if it is a single
        # item, we need to make copies to ensure that different zones'
        # compositions can be different
        if hasattr(comp, '__iter__'):
            if len(comp) != len(self.zones):
                raise despoticError(
                    "zonedcloud: could not " + 
                    "broadcast "+str(len(comp))+" items to " + 
                    str(len(self.zones)) + " zones")
            for z, x in zip(self.zones, comp):
                z.comp = deepcopy(x)
        else:
            for z in self.zones:
                z.comp = deepcopy(comp)

    @property
    def mu(self):
        # Return mean mass per free particle in each zone, in units of
        # H mass
        if self.zones[0].comp.mu == 0.0:
            for z in self.zones:
                z.comp.computeDerived(z.nH)
        return np.array([c.mu for c in self.comp])

    @mu.setter
    def mu(self, mu):
        # Disallow direct setting of mu
        raise despoticError("zonedcloud: cannot directly set mu")

    @property
    def muH(self):
        # Return mean mass per H nucleus in each zone, in units of
        # H mass
        if self.zones[0].comp.mu == 0.0:
            for z in self.zones:
                z.comp.computeDerived(z.nH)
        return np.array([c.muH for c in self.comp])

    @mu.setter
    def muH(self, muH):
        # Disallow direct setting of muH
        raise despoticError("zonedcloud: cannot directly set muH")

    # Note: dust and rad are the same for all zones
    @property
    def dust(self):
        return self.zones[0].dust

    @dust.setter
    def dust(self, dust):
        self.zones[0].dust = dust

    @property
    def rad(self):
        return self.zones[0].rad

    @rad.setter
    def rad(self, rad):
        self.zones[0].rad = rad

    # Emitters can be different for every cloud because they have
    # different abundances, but the same species must be included in
    # the emitters for every zone
    @property
    def emitters(self):
        return [z.emitters for z in self.zones]

    @emitters.setter
    def emitters(self, emitters):
        # emitters can be a single element or a list; if it is
        # listlike, each element of the list must have the same keys;
        # if it is a single item, we need to make copies for the
        # different zones so that the emitter abundances and level
        # populations can be different in different zones
        if hasattr(emitters, '__iter__'):
            for x in emitters:
                if set(x.keys()) != set(emitters[0].keys()):
                    raise despoticError(
                        "zonedcloud: all zones must have the " +
                        "same emitters")
            if len(emitters) != len(self.zones):
                raise despoticError(
                    "zonedcloud: could not " + 
                    "broadcast "+str(len(emitters))+" items to " + 
                    str(len(self.zones)) + " zones")
            for z, x in zip(self.zones, emitters):
                z.emitters = deepcopy(x)
        else:
            for z in self.zones:
                z.emitters = deepcopy(emitters)

    @property
    def chemnetwork(self):
        return [z.chemnetwork for z in self.zones]

    @chemnetwork.setter
    def chemnetwork(self, chemnetwork):
        # chemnetwork can be a single item or listlike; if it is a
        # list, we need to make sure that all items are of the same
        # type; if it is a single item, we need to make sure that each
        # zone has a separate copy in order to ensure that the
        # networks can be evolved independently
        if hasattr(chemnetwork, '__iter__'):
            for x in chemnetwork:
                if type(x) != type(chemnetwork[0]):
                    raise despoticError(
                        "zonedcloud: all zones must have the " + 
                        "same chemnetwork type")
            if len(chemnetwork) != len(self.zones):
                raise despoticError(
                    "zonedcloud: could not " + 
                    "broadcast "+str(len(chemnetwork))+" items to " + 
                    str(len(self.zones)) + " zones")
            for z, x in zip(self.zones, chemnetwork):
                z.chemnetwork = deepcopy(x)
        else:
            for z in self.zones:
                z.chemnetwork = deepcopy(chemnetwork)

    ####################################################################
    # These methods return the mass and radial location of each zone
    ####################################################################

    def radius(self, edge=False):
        """
        Return the radius of each zone.

        Parameters
           edge : Boolean
              if True, the value returned is the radii of the zone
              edges; otherwise it is the radii of the zone centers

        Returns
           rad : array
              radii of zone centers (default) or edges (if edge is
              True)

        Remarks
           if the geometry is 'slab', the values returned are the
           depths into the slab rather than the radii
        """

        dr = (self._colDen[1:] - self._colDen[:-1]) / self.nH
        redge = np.cumsum(dr[::-1])[::-1]
        redge = np.append(redge, 0.0)
        if edge:
            return redge
        else:
            return 0.5*(redge[1:]+redge[:-1])
            

    def mass(self, edge=False):
        """
        Returns the mass in each zone.

        Parameters
           edge : Boolean
              if True, the value returned gives the cumulative mass at
              each zone edge, starting from the outer edge; otherwise
              the value returned is the mass of each zone

        Returns
           mass : array
              mass of each zone (if edge is False), or cumulative mass
              to each zone edge (if edge is True)

        Remarks
           if the geometry is 'slab', the masses are undefined, and
           this return returns None
        """

        # In slab geometry, mass is undefined, so return None
        if self.geometry == 'slab':
            return None

        # In spherical geometry, get radial exent and density of each
        # zone, and use it to get mass of each zone
        nH = self.nH
        rho = nH*self.muH*mH
        dr = (self._colDen[1:] - self._colDen[:-1]) / nH
        redge = np.cumsum(dr[::-1])[::-1]
        redge = np.append(redge, 0.0)
        mass = 4./3.*np.pi*(redge[:-1]**3-redge[1:]**3)*rho

        # Return either mass or cumulative sum of masses
        if edge:
            return np.cumsum(mass)
        else:
            return mass

    ####################################################################
    # These methods get and set abundances; each comes in two
    # versions, one of which handles quantities that are mass-weighted
    # over the entire cloud, and the other of which deals with
    # zone-by-zone values
    ####################################################################

    @property
    def abundances_zone(self):
        """
        Return abundances of all emitting species in all zones
        """
        return [z.abundances for z in self.zones]

    @abundances_zone.setter
    def abundances_zone(self, other):
        """
        Set abundances in all zones
        """
        if hasattr(other, '__iter__'):
            for z, o in zip(self.zones, other):
                z.abundances = o
        else:
            for z in self.zones:
                z.abundances = other

    @property
    def abundances(self):
        """
        Returns abundances of all emitting species, mass-weighted over
        cloud
        """
        mass = self.mass()
        abd = self.abundances_zone
        abd_sum = mass[0]*abd[0]
        for m, a in zip(mass[1:], abd[1:]):
            abd_sum += m*a
        return abd_sum / np.sum(mass)

    @abundances.setter
    def abundances(self, other):
        """
        Set all abundances
        """
        if not isinstance(other, dict):
            for z, o in zip(self.zones, other):
                z.abundances = o
        else:
            for z in self.zones:
                z.abundances = other

    @property
    def chemabundances_zone(self):
        """
        Return abundances of all emitting species in all zones
        """
        return [z.chemabundances for z in self.zones]

    @chemabundances_zone.setter
    def chemabundances_zone(self, other):
        """
        Set abundances of all emitting species in all zones
        """
        if not isinstance(other, dict):
            for z, o in zip(self.zones, other):
                z.chemabundances = o
        else:
            for z in self.zones:
                z.chemabundances = other

    @property
    def chemabundances(self):
        """
        Returns abundances of all species in the chemical network,
        mass-weighted over the zonedcloud
        """
        mass = self.mass()
        abd = self.chemabundances_zone
        abd_sum = mass[0]*abd[0]
        for m, a in zip(mass[1:], abd[1:]):
            abd_sum += m*a
        return abd_sum / np.sum(mass)

    @chemabundances.setter
    def chemabundances(self, other):
        """
        Set all abundances in the chemical network
        """
        if hasattr(other, '__iter__'):
            for z, o in zip(self.zones, other):
                z.chemabundances = o
        else:
            for z in self.zones:
                z.chemabundances = other

    ####################################################################
    # These methods wrap the corresponding method in the cloud class,
    # and apply them to all zones
    ####################################################################

    def setVirial(self, alphaVir=1.0, NTonly=False):
        """
        This routine sets the velocity dispersion in all zones to the
        virial value

        Parameters
           alphaVir : float
              virial ratio to be used in computation; defaults to 1
           NTonly : Boolean
              if True, the virial ratio is computed considering only the
              non-thermal component of the velocity dispersion

        Returns
           Nothing
        """
        
        # Thermal velocity disperison squared
        if NTonly == False:
            sigmaThSqr = np.sum(self.mass()*kB*self.Tg / (self.mu*mH)) / \
                         np.sum(self.mass())
        else:
            sigmaThSqr = 0.0

        # Get total velocity dispersion
        m = np.sum(self.mass())
        r = self.radius(edge=True)[0]
        sigmaTotSqr = G*m*alphaVir / (5.0*r)

        # Set non-thermal part
        if sigmaTotSqr > sigmaThSqr:
            self.sigmaNT = np.sqrt(sigmaTotSqr - sigmaThSqr)
        else:
            self.sigmaNT = 0.0
            print("setVirial warning: setting sigmaNT = 0, " + 
                  "virial ratio still exceeds desired value")


    ####################################################################
    # Method to add an emitter; this is applied to every zone
    ####################################################################
    def addEmitter(self, emitName, emitAbundance, emitterFile=None,
                   emitterURL=None, energySkip=False, extrap=True):
        """
        Add an emitting species

        Pamameters
           emitName : string
              name of the emitting species
           emitAbundance : float or listlike
              abundance of the emitting species relative to H; if this
              is listlike, it must have the same number of elements as
              the number of zones
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
              if set to True, collision rate coefficients for this species
              will be extrapolated to temperatures outside the range given
              in the LAMDA table. If False, no extrapolation is perfomed,
              and providing temperatures outside the range in the table
              produces an error

        Returns
           Nothing
        """
        if hasattr(emitAbundance, '__iter__'):
            if len(emitAbundance) != len(self.zones):
                raise despoticError(
                    "zonedcloud: could not " + 
                    "broadcast "+str(len(emitAbundance))+" items to " + 
                    str(len(self.zones)) + " zones")
            for z, x in zip(self.zones, emitAbundance):
                z.addEmitter(emitName, x, emitterFile=emitterFile,
                             emitterURL=emitterURL, extrap=extrap,
                             energySkip=energySkip)
        else:
            for z in self.zones:
                z.addEmitter(emitName, emitAbundance, 
                             emitterFile=emitterFile,
                             emitterURL=emitterURL, extrap=extrap,
                             energySkip=energySkip)

    ####################################################################
    # Method to return the continuum-subtrated luminosity / intensity /
    # brightness temperature of lines from the specified emitter
    ####################################################################
    def lineLum(self, emitName, LTE=False, noClump=False,
                transition=None, thin=False, intOnly=False,
                TBOnly=False, lumOnly=False,
                escapeProbGeom=None, dampFactor=0.5,
                noRecompute=False, abstol=1.0e-8,
                verbose=False):
        """
        Return the frequency-integrated intensity of various lines,
        summed over all zones

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
           intOnly: Boolean
              if true, the output is simply an array containing the
              frequency-integrated intensity of the specified lines;
              mutually exclusive with TBOnly and lumOnly
           TBOnly: Boolean
              if True, the output is simply an array containing the
              velocity-integrated brightness temperatures of the
              specified lines; mutually exclusive with intOnly and
              lumOnly
           lumOnly: Boolean
              if True, the output is simply an array containing the
              luminosity per H nucleus in each of the specified lines;
              mutually eclusive with intOnly and TBOonly
           escapeProbGeom : 'sphere' | 'slab' | 'LVG'
              sets problem geometry that will be assumed in calculating
              escape probabilities; ignored if the escape probabilities
              are already initialized; if left as None, escapeProbGeom
              = self.geometry
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

        # Set geometry
        if escapeProbGeom is None:
            escapeProbGeom = self.geometry

        # Compute zone by zone luminosity
        zoneLum = [z.lineLum(emitName, LTE=LTE, noClump=noClump,
                             transition=transition, thin=thin,
                             intOnly=intOnly, TBOnly=TBOnly,
                             lumOnly=lumOnly,
                             escapeProbGeom=escapeProbGeom,
                             dampFactor=dampFactor,
                             noRecompute=noRecompute) 
                   for z in self.zones]

        # Sum over zones
        if intOnly:

            # Integrated intensities just add
            zoneLum = np.array(zoneLum)
            lineLum = np.sum(zoneLum, axis=0)
            return lineLum

        elif TBOnly:

            # For velocity-integrated brightness temperatures, we need
            # to back out the intensity, add the intensities, then
            # turn back into a brightness temperature
            intTB = np.array(zoneLum)
            if transition is None:
                u = self.zones[0].emitters[emitName].data.radUpper
                l = self.zones[0].emitters[emitName].data.radLower
            else:
                u = transition[0]
                l = transition[1]
            freq = self.zones[0].emitters[emitName].data.freq[u,l]
            intIntensity \
                = 2*h*freq**3/(c**2 * (np.exp(h*c/(kB*intTB*1e5)) - 1.0))
            intIntensity = np.sum(intIntensity, axis=0)
            intTB = h*c/kB / \
                    np.log(1.0 + 2.0*h*freq**3/(c**2*intInensity)) / 1e5
            return intTB

        elif lumOnly:

            # For luminosity per H, we need to mass-weight
            lumPerH = np.array(zoneLum)
            mass = self.mass()
            return np.sum(mass*np.transpose(lumPerH), axis=1) \
                / np.sum(mass)

        else:

            # Returning full output, so use each of the summing
            # procedures above
            lineLum = zoneLum[0]

            # Intensities just add
            intIntensity = np.array(
                [[x['intIntensity'] for x in z] for z in zoneLum])
            intIntensity = np.sum(intIntensity, axis=0)
            for l, i in zip(lineLum, intIntensity):
                l['intIntensity'] = i

            # Derive velocity-integrated brightness temperature from
            # intensity
            freq = np.array([x['freq'] for x in zoneLum[0]])
            intTB = h*c/kB / \
                    np.log(1.0 + 2.0*h*freq**3/ (c**2*intIntensity)) \
                    / 1e5
            for l, t in zip(lineLum, intTB):
                l['intTB'] = t

            # Derive luminosity per unit mass by mass-weighing zones
            lumPerH = np.array(
                [[x['lumPerH'] for x in z] for z in zoneLum])
            mass = self.mass()
            lumPerH = np.sum(mass*np.transpose(lumPerH), axis=1) \
                      / np.sum(mass)
            for l, t in zip(lineLum, lumPerH):
                l['lumPerH'] = t

            # Excitation temperatues and optical depths really only
            # make sense zone by zone, so we return these as arrays
            for i, l in enumerate(lineLum):
                l['Tex'] = np.array([z[i]['Tex'] for z in zoneLum])
                l['tau'] = np.array([z[i]['tau'] for z in zoneLum])

            # Return
            return lineLum


    ####################################################################
    # Method to calculate instantaneous values of all heating, cooling
    # terms
    ####################################################################
    def dEdt(self, c1Grav=0.0, thin=False, LTE=False, 
             fixedLevPop=False, noClump=False, 
             escapeProbGeom=None, PsiUser=None, 
             sumOnly=False, dustOnly=False, gasOnly=False, 
             dustCoolOnly=False, dampFactor=0.5, 
             verbose=False, overrideSkip=False,
             zones=False):
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
           escapeProbGeom : 'sphere' | 'slab' | 'LVG'
              sets problem geometry that will be assumed in calculating
              escape probabilities; ignored if the escape probabilities
              are already initialized; if left as None, escapeProbGeom
              = self.geometry
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
           zones : Boolean
              if True, heating and cooling rates are returned for each
              zone; if False, the values returned are mass-weighted over
              the entire cloud

        Returns
           rates : dict
             A dict containing the values of the various heating and
             cooling rate terms; all quantities are in units of erg s^-1
             H^-1, and by convention positive = heating, negative =
             cooling; for dust-gas exchange, positive indicates heating
             of gas, cooling of dust. By default these quantities are
             mass-weighted over the entire cloud, but if zones is True
             then they are returned as arrays for each zone

           Elements of the dict are as follows by default, but can be
           altered by the additional keywords listed below:

           GammaPE : float or array
              photoelectric heating rate
           GammaCR : float or array
              cosmic ray heating rate
           GammaGrav : float or array
              gravitational contraction heating rate
           LambdaLine : dict or array
              cooling rate from lines; dictionary keys correspond to
              species in the emitter list, values give line cooling
              rate for that species
           PsiGD : float or array
              dust-gas energy exchange rate
           GammaDustISRF : float or array
              dust heating rate due to the ISRF
           GammaDustCMB : float or array
              dust heating rate due to CMB
           GammaDustIR : float or array
              dust heating rate due to IR field
           GammaDustLine : float or array
              dust heating rate due to absorption of line radiation
           LambdaDust : float or array
              dust cooling rate due to thermal emission
           PsiUserGas : float or array
              gas heating / cooling rate from user-specified
              function; only included if PsiUser != None
           PsiUserDust : float or array
              gas heating / cooling rate from user-specified
              function; only included is PsiUser != None
        """

        # Set geometry
        if escapeProbGeom is None:
            escapeProbGeom = self.geometry

        # Do computation for each zone
        dEdt_zones \
            = [z.dEdt(c1Grav=c1Grav, thin=thin, LTE=LTE, 
                      fixedLevPop=fixedLevPop, noClump=noClump, 
                      escapeProbGeom=escapeProbGeom, PsiUser=PsiUser, 
                      sumOnly=sumOnly, dustOnly=dustOnly, gasOnly=gasOnly, 
                      dustCoolOnly=dustCoolOnly, dampFactor=0.5, 
                      verbose=verbose, overrideSkip=overrideSkip)
               for z in self.zones]

        # Form the sum or make the arrays
        dEdt = {}
        if not zones:
            mass = self.mass()
        for k in dEdt_zones[0].keys():
            if zones:
                dEdt[k] = np.array([z[k] for z in dEdt_zones])
            else:
                if k != 'LambdaLine':
                    dEdt[k] = np.sum(mass*
                                     np.array([z[k] for z in
                                               dEdt_zones])) / \
                        np.sum(mass)
                else:
                    lambdaLine = {}
                    for s in dEdt_zones[0][k].keys():
                        lambdaLine[s] \
                            = np.sum(mass*
                                     np.array([z[k][s] for z in
                                               dEdt_zones])) / \
                            np.sum(mass)
                    dEdt[k] = lambdaLine

        # Return
        return dEdt

    ####################################################################
    # Methods that just wrap the corresponding cloud methods, and
    # apply them zone by zone
    ####################################################################
    def setDustTempEq(self, **kwargs):
        """
        Set Td to equilibrium dust temperature at fixed Tg

        Parameters
           kwargs : dict
              these arguments are passed through to the corresponding
              function for each zone

        Returns
           success : Boolean
              True if dust temperature calculation converged, False if
              not
        """
        conv = []
        for i, z in enumerate(self.zones):
            if 'verbose' in kwargs.keys():
                if kwargs['verbose']:
                    print("Finding equilibrium for zone " +
                          str(i+1) + " / " + str(len(self.zones)))
            conv.append(z.setDustTempEq(**kwargs))
        return np.all(conv)


    def setGasTempEq(self, **kwargs):
        """
        Set Tg to equilibrium gas temperature at fixed Td

        Parameters
           kwargs : dict
              these arguments are passed through to the corresponding
              function for each zone

        Returns
           success : Boolean
              True if the calculation converges, False if it does not

        Remarks
           if the key escapeProbGeom is not in kwargs, then
           escapeProbGeom will be set to self.geometry
        """
        if 'escapeProbGeom' not in kwargs.keys():
            kwargs['escapeProbGeom'] = self.geometry
        conv = []
        for i, z in enumerate(self.zones):
            if 'verbose' in kwargs.keys():
                if kwargs['verbose']:
                    print("Finding equilibrium for zone " +
                          str(i+1) + " / " + str(len(self.zones)))
            conv.append(z.setGasTempEq(**kwargs))
        return np.all(conv)


    def setTempEq(self, **kwargs):
        """
        Set Tg and Td to equilibrium gas and dust temperatures

        Parameters
           kwargs : dict
              these arguments are passed through to the corresponding
              function for each zone

        Returns
           success : Boolean
              True if the calculation converges, False if it does not

        Remarks
           if the key escapeProbGeom is not in kwargs, then
           escapeProbGeom will be set to self.geometry
        """
        if 'escapeProbGeom' not in kwargs.keys():
            kwargs['escapeProbGeom'] = self.geometry
        conv = []
        for i, z in enumerate(self.zones):
            if 'verbose' in kwargs.keys():
                if kwargs['verbose']:
                    print("Finding equilibrium for zone " +
                          str(i+1) + " / " + str(len(self.zones)))
            conv.append(z.setTempEq(**kwargs))
        return np.all(conv)


    def setChemEq(self, **kwargs):
        """
        Set the chemical abundances for a cloud to their equilibrium
        values, computed using a specified chemical network.

        Parameters
           kwargs : dict
              these arguments are passed through to the corresponding
              function for each zone

        Returns
           success : Boolean
              True if the calculation converges, False if it does not

        Remarks
           if the key escapeProbGeom is not in kwargs, then
           escapeProbGeom will be set to self.geometry
        """
        if 'tempEqParam' not in kwargs.keys():
            kwargs['tempEqParam'] \
                = { 'escapeProbGeom' : self.geometry }
        elif 'escapeProbGeom' not in kwargs['tempEqParam'].keys():
            kwargs['tempEqParam']['escapeProbGeom'] = self.geometry
        conv = []
        for i, z in enumerate(self.zones):
            if 'verbose' in kwargs.keys():
                if kwargs['verbose']:
                    print("Finding equilibrium for zone " +
                          str(i+1) + " / " + str(len(self.zones)))
            conv.append(z.setChemEq(**kwargs))
        return np.all(conv)
