"""
This module defines the pwind class; this class is mostly just a thin
wrapper over the underlying pwind c++ code, with some extra
functionality to integrate with despotic for the purposes of
calculating atomic properties, and to provide vector broadcasting.
"""

########################################################################
# Copyright (C) 2016 Mark Krumholz
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

from .pwind_interface import libpwind
from .pwind_util import pM, pA
import numpy as np
from ctypes import c_double

########################################################################
# Physical constants, in cgs units
########################################################################
import scipy.constants as physcons
me = physcons.m_e * 1e3
mH = physcons.m_p * 1e3 + me
h = physcons.h * 1e7
c = physcons.c * 1e2
G = physcons.G * 1e3
kB = physcons.k * 1e7
e = physcons.e * c/1e1

########################################################################
# Define the main pwind class
########################################################################

class pwind(object):

    ################################################################
    # Constructor
    ################################################################
    def __init__(self, Gamma, mach, 
                 driver='ideal', potential='point',
                 expansion='intermediate',
                 geometry='sphere', fcrit=1.0,
                 epsabs=1.0e-4, epsrel=1.0e-2,
                 theta=None, phi=None, theta_in=None,
                 tau0=None, uh=None, interpabs=1.0e-2,
                 interprel=1.0e-2):
        """
        Creates a generic wind object to compute observable property
        of winds

        Parameters
           Gamma : float
              Eddington factor
           mach : float
              Mach number
           driver : 'ideal' | 'radiation' | 'hot'
              wind driving mechanism; allowed values are 'ideal'
              (ideal momentum-driven wind), 'radiation'
              (radiation-driven wind), and 'hot' (wind driven by hot
              gas entrainment)
           potential : 'point' | 'isothermal'
              gravitational potential confining the wind; allowed
              values are 'point' and 'isothermal', for point and
              isothermal potentials, respectively
           expansion : 'area' | 'intermediate' | 'solid angle'
              cloud expansion law: clouds can maintain constant area
              ('area'), maintain constant solid angle ('solid angle'),
              or have intermediate behavior ('intermediate') where the
              area increases with distance as r and the solid angle
              decreases as 1/r
           geometry : 'sphere' | 'cone' | 'cone sheath'
              geometry of the wind; allowed values are sphere (covers
              all space), conical, and conical sheath (meaning that
              the wind is bounded by an inner and outer cone)
           fcrit : float
              material is only considered to be launched into the wind
              if x < f_crit x_crit; must be <= 1.0
           epsabs : float
              absolute error tolerance for numerical integrations
           epsrel : float
              relative error tolerance for numerical integrations
           theta : float, in the range (0, pi/2)
              opening angle of the outer edge of the wind for cone or
              cone sheath geometry; ignored for all other geometries
           phi : float, in the range [-pi/2, pi/2]
              inclination of the wind cone central axis relative to
              the plane of the sky for either cone or cone sheath
              geometry; ignored for all other geometries; phi = 0
              corresponds to a wind cone in the plane of the sky, phi > 0
              corresponds to the varpi > 0 side of the wind pointed
              away from the observer
           theta_in : float, in range (0, theta)
              opening angle of the inner edge of the wind for cone
              sheath geometry, ignored for all other geometries
           tau0 : float
              for radiation-driven winds, the optical depth at the
              mean surface density; must be specified for
              radiation-driven winds, ignored for all other drivers
           uh : float
              for hot gas-driven winds, the hot gas speed relative to
              the escape speed; must be specified for hot gas-driven
              winds, ignored for all other drivers
           interpabs : float
              absolute error tolerance of the iterpolation tables used
              for hot gas-driven winds; ignored for all other drivers;
              note that using smaller tolerances quickly becomes very
              expensive in memory and computation time, so values much
              below 10^-2 are not recommended
           interprel : float
              same as interpabs, but giving relative rather tahn
              absolute error tolerance
        """

        # Record information
        self.Gamma_ = Gamma
        self.mach_ = mach
        self.potential_ = potential
        self.expansion_ = expansion
        self.geometry_ = geometry
        self.fcrit_ = fcrit
        self.theta_ = theta
        self.phi_ = phi
        self.theta_in_ = theta_in
        self.driver_ = driver
        self.epsabs_ = epsabs
        self.epsrel_ = epsrel
        self.tau0_ = tau0
        self.uh_ = uh
        self.interpabs_ = interpabs
        self.interprel_ = interprel

        # Build the c++ class we are wrapping
        self.__pw = None
        self.__geom = None
        self.__init_geom()
        self.__init_lib()

    ################################################################
    # Destructor
    ################################################################
    def __del__(self):
        if self.__pw is not None:
            libpwind.pwind_free(self.__pw)
            libpwind.pwind_geom_free(self.__geom)

    ################################################################
    # Method to instantiate the c++ geometry class
    ################################################################
    def __init_geom(self):
        if self.__geom is not None:
            libpwind.pwind_geom_free(self.__geom)
            self.__geom = None
        if self.geometry == 'sphere' :
            self.__geom = libpwind.pwind_geom_sphere_new()
        elif self.geometry == 'cone' :
            if self.theta is None or self.phi is None:
                raise ValueError("pwind: for cone geometry, must "
                                 "set theta and phi")
            if self.theta <=0 or self.theta >= np.pi/2.0 or \
               self.phi < -np.pi/2 or self.phi > np.pi/2:
                raise ValueError("pwind: for cone geometry, must "
                                 "have theta in (0, pi/2) and phi "
                                 "in [-pi/2, pi/2]")
            self.__geom = libpwind.pwind_geom_cone_new(self.theta, self.phi)
        elif self.geometry == 'cone sheath' :
            if self.theta is None or self.phi is None:
                raise ValueError("pwind: for cone sheaht geometry, must "
                                 "set theta, phi, and theta_in")
            if self.theta <=0 or self.theta >= np.pi/2.0 or \
               self.phi < -np.pi/2 or self.phi > np.pi/2 or \
               self.theta_in <= 0 or self.theta_in >= self.theta:
                raise ValueError("pwind: for cone sheath geometry, must "
                                 "have theta in (0, pi/2), phi "
                                 "in [-pi/2, pi/2], and theta_in "
                                 "in (0, theta)")
            self.__geom = libpwind.pwind_geom_cone_sheath_new(
                self.theta, self.theta_in, self.phi)
        else:
            raise ValueError("pwind: unknown geometry '" +
                             str(self.geometry) + "'")
        
    ################################################################
    # Method to instantiate the c++ class wind
    ################################################################
    def __init_lib(self):
        if self.Gamma_ <= 0.0:
            raise ValueError("pwind: Gamma must be > 0")
        if self.mach_ <= 0.0:
            raise ValueError("pwind: mach must be > 0")
        if self.fcrit_ > 1.0:
            raise ValueError("pwind: fcrit must be <= 1")
        if self.__pw is not None:
            libpwind.pwind_free(self.__pw)
            self.__pw = None
        if self.driver_ == 'ideal':
            if self.potential_ == 'point':
                if self.expansion_ == 'area':
                    self.__pw = libpwind.pwind_ideal_pa_new(
                        self.Gamma, self.mach, self.__geom,
                        self.epsabs_, self.epsrel_, self.fcrit_)
                elif self.expansion_ == 'intermediate':
                    self.__pw = libpwind.pwind_ideal_pi_new(
                        self.Gamma, self.mach, self.__geom,
                        self.epsabs_, self.epsrel_, self.fcrit_)
                elif self.expansion_ == 'solid angle':
                    self.__pw = libpwind.pwind_ideal_ps_new(
                        self.Gamma, self.mach, self.__geom,
                        self.epsabs_, self.epsrel_, self.fcrit_)
                else: raise(ValueError("pwind: unknown expansion "
                                      "value "+str(self.expansion_)))
            elif self.potential_ == 'isothermal':
                if self.expansion_ == 'area':
                    self.__pw = libpwind.pwind_ideal_ia_new(
                        self.Gamma, self.mach, self.__geom,
                        self.epsabs_, self.epsrel_, self.fcrit_)
                elif self.expansion_ == 'intermediate':
                    self.__pw = libpwind.pwind_ideal_ii_new(
                        self.Gamma, self.mach, self.__geom,
                        self.epsabs_, self.epsrel_, self.fcrit_)
                elif self.expansion_ == 'solid angle':
                    self.__pw = libpwind.pwind_ideal_is_new(
                        self.Gamma, self.mach, self.__geom,
                        self.epsabs_, self.epsrel_, self.fcrit_)
                else: raise(ValueError("pwind: unknown expansion "
                                      "value "+str(self.expansion_)))
            else: raise(ValueError("pwind: unknown potential " +
                                  str(self.potential_)))
        elif self.driver_ == 'radiation':
            if self.tau0 is None:
                raise(ValueError('pwind: for radiation-driven wind, '
                                 'must explicitly set tau0'))
            if self.Gamma*self.tau0 <= 1.0:
                raise ValueError('pwind: for radiation-driven wind, ',
                                 'Gamma*tau0 must be > 1')
            if self.potential_ == 'point':
                if self.expansion_ == 'area':
                    self.__pw = libpwind.pwind_rad_pa_new(
                        self.Gamma, self.mach, self.tau0,
                        self.__geom, self.epsabs_, self.epsrel_, self.fcrit_)
                elif self.expansion_ == 'intermediate':
                    self.__pw = libpwind.pwind_rad_pi_new(
                        self.Gamma, self.mach, self.tau0,
                        self.__geom, self.epsabs_, self.epsrel_, self.fcrit_)
                elif self.expansion_ == 'solid angle':
                    self.__pw = libpwind.pwind_rad_ps_new(
                        self.Gamma, self.mach, self.tau0,
                        self.__geom, self.epsabs_, self.epsrel_, self.fcrit_)
                else: raise(ValueError("pwind: unknown expansion "
                                      "value "+str(self.expansion_)))
            elif self.potential_ == 'isothermal':
                if self.expansion_ == 'area':
                    self.__pw = libpwind.pwind_rad_ia_new(
                        self.Gamma, self.mach, self.tau0,
                        self.__geom, self.epsabs_, self.epsrel_, self.fcrit_)
                elif self.expansion_ == 'intermediate':
                    self.__pw = libpwind.pwind_rad_ii_new(
                        self.Gamma, self.mach, self.tau0,
                        self.__geom, self.epsabs_, self.epsrel_, self.fcrit_)
                elif self.expansion_ == 'solid angle':
                    self.__pw = libpwind.pwind_rad_is_new(
                        self.Gamma, self.mach, self.tau0,
                        self.__geom, self.epsabs_, self.epsrel_, self.fcrit_)
                else: raise(ValueError("pwind: unknown expansion "
                                      "value "+str(self.expansion_)))
            else: raise(ValueError("pwind: unknown potential " +
                                  str(self.potential_)))
        elif self.driver_ == 'hot':
            if self.uh is None:
                raise(ValueError('pwind: for hot gas-driven wind, '
                                 'must explicitly set uh'))
            if self.uh <= 0.0:
                raise ValuError('pwind: for hot gas-driven wind, '
                                'uh must be > 0')
            if self.potential_ == 'point':
                if self.expansion_ == 'area':
                    self.__pw = libpwind.pwind_hot_pa_new(
                        self.Gamma, self.mach, self.uh,
                        self.__geom, self.epsabs_, self.epsrel_,
                        self.interpabs_, self.interprel_, self.fcrit_)
                elif self.expansion_ == 'intermediate':
                    self.__pw = libpwind.pwind_hot_pi_new(
                        self.Gamma, self.mach, self.uh,
                        self.__geom, self.epsabs_, self.epsrel_,
                        self.interpabs_, self.interprel_, self.fcrit_)
                elif self.expansion_ == 'solid angle':
                    self.__pw = libpwind.pwind_hot_ps_new(
                        self.Gamma, self.mach, self.uh,
                        self.__geom, self.epsabs_, self.epsrel_,
                        self.interpabs_, self.interprel_, self.fcrit_)
                else: raise(ValueError("pwind: unknown expansion "
                                      "value "+str(self.expansion_)))
            elif self.potential_ == 'isothermal':
                if self.expansion_ == 'area':
                    self.__pw = libpwind.pwind_hot_ia_new(
                        self.Gamma, self.mach, self.uh,
                        self.__geom, self.epsabs_, self.epsrel_,
                        self.interpabs_, self.interprel_, self.fcrit_)
                elif self.expansion_ == 'intermediate':
                    self.__pw = libpwind.pwind_hot_ii_new(
                        self.Gamma, self.mach, self.uh,
                        self.__geom, self.epsabs_, self.epsrel_,
                        self.interpabs_, self.interprel_, self.fcrit_)
                elif self.expansion_ == 'solid angle':
                    self.__pw = libpwind.pwind_hot_is_new(
                        self.Gamma, self.mach, self.uh,
                        self.__geom, self.epsabs_, self.epsrel_,
                        self.interpabs_, self.interprel_, self.fcrit_)
                else: raise(ValueError("pwind: unknown expansion "
                                       "value "+str(self.expansion_)))
            else: raise(ValueError("pwind: unknown potential " +
                                   str(self.potential_)))
        else: raise(ValueError("pwind: unknown driver "+
                               str(self.driver_)))

    ################################################################
    # Define cloud properties; these are things where, if we change
    # them, we need to recompute some dependent quantity
    ################################################################
    @property
    def driver(self):
        return self.driver_
    @driver.setter
    def driver(self, val):
        self.driver_ = val
        self.__init_lib()
    @property
    def expansion(self):
        return self.expansion_
    @expansion.setter
    def expansion(self, val):
        self.expansion_ = val
        self.__init_lib()
    @property
    def potential(self):
        return self.potential_
    @potential.setter
    def potential(self, val):
        self.potential_ = val
        self.__init_lib()
    @property
    def Gamma(self):
        return self.Gamma_
    @Gamma.setter
    def Gamma(self, val):
        self.Gamma_ = Gamma
        self.__init_lib()
    @property
    def mach(self):
        return self.mach_
    @mach.setter
    def mach(self, val):
        self.mach_ = mach
        libpwind.set_mach(val, self.__pw)
    @property
    def geometry(self):
        return self.geometry_
    @geometry.setter
    def geometry(self, val):
        if val == self.geometry_:
            return
        self.geometry_ = val
        self.__init_geom()
        libpwind.set_geometry(self.__pw, self.__geom)        
    @property
    def theta(self):
        return self.theta_
    @theta.setter
    def theta(self, val):
        self.theta_ = val
        if self.geometry == 'cone' or self.geometry == 'cone sheath':
            self.__init_geom()
            libpwind.set_geometry(self.__pw, self.__geom)
    @property
    def phi(self):
        return self.phi_
    @phi.setter
    def phi(self, val):
        self.phi_ = val
        if self.geometry == 'cone' or self.geometry == 'cone sheath':
            self.__init_geom()
            libpwind.set_geometry(self.__pw, self.__geom)
    @property
    def theta_in(self):
        return self.theta_in_
    @theta_in.setter
    def theta_in(self, val):
        self.theta_in_ = val
        if self.geometry == 'cone sheath':
            self.__init_geom()
            libpwind.set_geometry(self.__pw, self.__geom)
    @property
    def epsabs(self):
        return self.epsabs_
    @epsabs.setter
    def epsabs(self, val):
        self.epsabs_ = val
        libpwind.set_epsabs(val, self.__pw)
    @property
    def epsrel(self):
        return self.epsrel_
    @epsrel.setter
    def epsrel(self, val):
        self.epsrel_ = val
        libpwind.set_epsrel(val, self.__pw)
    @property
    def fcrit(self):
        return self.fcrit_
    @fcrit.setter
    def fcrit(self, val):
        self.fcrit_ = val
        libpwind.set_fcrit(val, self.__pw)
    @property
    def tau0(self):
        return self.tau0_
    @tau0.setter
    def tau0(self, val):
        self.tau0_ = val
        if self.driver == 'radiation':
            self.__init_lib()
    @property
    def uh(self):
        return self.uh_
    @uh.setter
    def uh(self, val):
        self.uh_ = val
        if self.driver == 'hot':
            self.__init_lib()
    @property
    def interpabs(self):
        return self.interpabs_
    @interpabs.setter
    def interpabs(self, val):
        self.interpabs_ = val
        if self.driver == 'hot':
            self.__init_lib()
    @property
    def interprel(self):
        return self.interprel_
    @interprel.setter
    def interprel(self, val):
        self.interprel_ = val
        if self.driver == 'hot':
            self.__init_lib()
    @property
    def sx(self):
        return libpwind.sx(self.__pw)
    @sx.setter
    def sx(self, val):
        raise NotImplementedError(
            "pwind: can't set sx directly")
    @property
    def xcrit(self):
        return libpwind.xcrit(self.__pw)
    @xcrit.setter
    def xcrit(self, val):
        raise NotImplementedError(
            "pwind: can't set xcrit directly")
    @property
    def umax(self):
        return libpwind.umax(self.__pw)
    @umax.setter
    def umax(self, val):
        raise NotImplementedError(
            "pwind: can't set umax directly")
    @property
    def zetaM(self):
        return libpwind.zetaM(self.__pw)
    @zetaM.setter
    def zetaM(self, val):
        raise NotImplementedError(
            "pwind: can't set zetaM directly")
    @property
    def zetaA(self):
        return libpwind.zetaA(self.__pw)
    @zetaA.setter
    def zetaA(self, val):
        raise NotImplementedError(
            "pwind: can't set zetaA directly")

    ################################################################
    # Wind kinematics methods
    ################################################################
    def y(self, a):
        """
        Return the cloud area function y(a) for this wind

        Parameters:
           a : float or arraylike
              dimensionless radius

        Returns:
           y : float or arraylike
              cloud area function
        """
        if hasattr(a, '__iter__'):
            result = np.zeros(np.asarray(a).shape)
            for i, a1 in enumerate(np.asarray(a).flat):
                result.flat[i] = libpwind.y(a1, self.__pw)
            return result
        else:
            return libpwind.y(a, self.__pw)

    def dyda(self, a):
        """
        Return the derivative of the cloud area function dy/da for
        this wind

        Parameters:
           a : float or arraylike
              dimensionless radius

        Returns:
           dyda : float or arraylike
              derivative of cloud area function
        """
        if hasattr(a, '__iter__'):
            result = np.zeros(np.asarray(a).shape)
            for i, a1 in enumerate(np.asarray(a).flat):
                result.flat[i] = libpwind.dyda(a1, self.__pw)
            return result
        else:
            return libpwind.dyda(a, self.__pw)

    def m(self, a):
        """
        Return the potential shape function m(a) for this wind

        Parameters:
           a : float or arraylike
              dimensionless radius

        Returns:
           m : float or arraylike
              potential shape function
        """
        if hasattr(a, '__iter__'):
            result = np.zeros(np.asarray(a).shape)
            for i, a1 in enumerate(np.asarray(a).flat):
                result.flat[i] = libpwind.m(a1, self.__pw)
            return result
        else:
            return libpwind.m(a, self.__pw)
        
    def U2(self, x, a):
        """
        Return the square velocity U_a^2(x) for this wind

        Parameters:
           x : float or arraylike
              dimensionless column density
           a : float or arraylike
              dimensionless radius

        Returns:
           U2 : float or arraylike
              dimensionless square velocity; a value of nan is returned
              for any combinations of x and a that are forbidden for
              this wind
        """
        bcast = np.broadcast(x, a)
        result = np.zeros(bcast.shape)
        i = 0
        for (x_, a_) in bcast:
            # The c++ routines don't do bounds checking for reasons of
            # speed, but we will for safety
            if a_ < 1.0:
                result.flat[i] = np.nan
            elif a_ >= libpwind.amax_abs(self.__pw):
                result.flat[i] = np.nan
            else:
                xlim = np.zeros(2)
                libpwind.xlimits(a_, self.__pw, xlim)
                if x_ < xlim[0] or x_ > xlim[1]:
                    result.flat[i] = np.nan
                else:
                    result.flat[i] = libpwind.U2(x_, a_, self.__pw)
            i = i+1
        if (not hasattr(x, '__iter__')) and \
           (not hasattr(a, '__iter__')):
            result = float(result)
        return result

    def dU2dx(self, x, a):
        """
        Return d/dx(U_a^2(x)) for this wind

        Parameters:
           x : float or arraylike
              dimensionless column density
           a : float or arraylike
              dimensionless radius

        Returns:
           dU2dx : float or arraylike
              dimensionless d/dx(U^2); a value of nan is returned
              for any combinations of x and a that are forbidden for
              this wind
        """
        bcast = np.broadcast(x, a)
        result = np.zeros(bcast.shape)
        i = 0
        for (x_, a_) in bcast:
            # The c++ routines don't do bounds checking for reasons of
            # speed, but we will for safety
            if a_ < 1.0:
                result.flat[i] = np.nan
            elif a_ >= libpwind.amax_abs(self.__pw):
                result.flat[i] = np.nan
            else:
                xlim = np.zeros(2)
                libpwind.xlimits(a_, self.__pw, xlim)
                if x_ < xlim[0] or x_ > xlim[1]:
                    result.flat[i] = np.nan
                else:
                    result.flat[i] = libpwind.dU2dx(x_, a_, self.__pw)
            i = i+1
        if (not hasattr(x, '__iter__')) and \
           (not hasattr(a, '__iter__')):
            result = float(result)
        return result

    def dU2da(self, x, a):
        """
        Return d/da(U_a^2(x)) for this wind

        Parameters:
           x : float or arraylike
              dimensionless column density
           a : float or arraylike
              dimensionless radius

        Returns:
           dU2da : float or arraylike
              dimensionless d/da(U^2); a value of nan is returned
              for any combinations of x and a that are forbidden for
              this wind
        """
        bcast = np.broadcast(x, a)
        result = np.zeros(bcast.shape)
        i = 0
        for (x_, a_) in bcast:
            # The c++ routines don't do bounds checking for reasons of
            # speed, but we will for safety
            if a_ < 1.0:
                result.flat[i] = np.nan
            elif a_ >= libpwind.amax_abs(self.__pw):
                result.flat[i] = np.nan
            else:
                xlim = np.zeros(2)
                libpwind.xlimits(a_, self.__pw, xlim)
                if x_ < xlim[0] or x_ > xlim[1]:
                    result.flat[i] = np.nan
                else:
                    result.flat[i] = libpwind.dU2dx(x_, a_, self.__pw)
            i = i+1
        if (not hasattr(x, '__iter__')) and \
           (not hasattr(a, '__iter__')):
            result = float(result)
        return result

    def X(self, ur, a):
        """
        Return the column density X_a(ur) for this wind

        Parameters:
           ur : float
              dimensionless radial velocity
           a : float
              dimensionless radius

        Returns:
           X : array
              dimensionless column densities; a value of nan is
              returned for combinations of ur and a that are forbitten
              for this wind
        """
        bcast = np.broadcast(ur, a)
        result = np.zeros(bcast.shape)
        i = 0
        for (ur_, a_) in bcast:
            # The c++ routines don't do bounds checking for reasons of
            # speed, but we will for safety
            if a_ <= 1.0:
                result.flat[i] = np.nan
            elif a_ >= libpwind.amax_abs(self.__pw):
                result.flat[i] = np.nan
            else:
                alim = np.zeros(4)
                nlim = libpwind.alimits(ur_, 0.0, 0.0, self.__pw, alim)
                alim.resize(nlim)
                result.flat[i] = np.nan
                for j in range(alim.size/2):
                    if a_ >= alim[2*j] and a_ <= alim[2*j+1]:
                        result.flat[i] \
                            = libpwind.X(ur_, a_, self.__pw)
                        break
            i+=1
        if (not hasattr(ur, '__iter__')) and \
           (not hasattr(a, '__iter__')):
            result = float(result)
        return result

    ################################################################
    # Some utility functions
    ################################################################
    def U(self, x, a):
        """
        Return the radial velocity U_a(x) for this wind

        Parameters:
           x : float or arraylike
              dimensionless column density
           a : float or arraylike
              dimensionless radius

        Returns:
           U : float or arraylike
              dimensionless radial velocity
        """
        return np.sqrt(self.U2(x,a))
    
    def rho(self, x, a):
        """
        Return the density divided by Mdot/(4 pi r0^2 v0) as a
        function of x and a

        Parameters:
           x : float or arraylike
              dimensionless column density
           a : float or arraylike
              dimensionless radius

        Returns:
           rho : float or arraylike
              dimensionless density
        """
        return pM(x, self.sx) / \
            (pA(x, self.sx) * self.U(x,a) * self.y(a))

    def drhodx(self, x, a):
        """
        Return the differential contribution to the mean density,
        drho_mean/dx, from material with initial surface density x at
        radius a

        Parameters:
           x : float or arraylike
              dimensionless column density
           a : float or arraylike
              dimensionless radius

        Returns:
           drhodx : float or arraylike
              dimensionless differential density
        """
        return pM(x, self.sx) / (a**2 * self.U(x,a))
        

    def dfcdx(self, x, a):
        """
        Return the differential covering factor dfc/dx as a
        function of x and a

        Parameters:
           x : float or arraylike
              dimensionless column density
           a : float or arraylike
              dimensionless radius

        Returns:
           dfcdx : float or arraylike
              differential covering factor
        """
        return pA(x, self.sx)*self.y(a)/a**2

    def fc(self, a):
        """
        Return the covering factor f_c as a function of radius

        Parameters:
           a : float or arraylike
              dimensionless radius

        Returns:
           fc : float or array
              covering fraction
        """
        return self.zetaA*self.y(a)/a**2

    def f_area(self):
        """
        Return fraction of the unit sphere covered by the wind for the
        current geometry; distinct from the covering fraction f_c,
        which is the fraction of the area within the wind geometric
        region that is covered

        Parameters:
           None

        Returns:
           f_area : float
              area fraction
        """
        if self.geometry == 'sphere':
            return 1.0
        elif self.geometry == 'cone':
            return 1.0-np.cos(self.theta)
        elif self.geometry == 'cone sheath':
            return np.cos(self.theta_in) - np.cos(theta_out)

    def dfcda(self, a):
        """
        Return derivative of the covering factor df_c / da

        Parameters:
           a : float or arraylike
              dimensionless radius

        Returns:
           dfcda : float or arraylike
              value of df_c/da
        """
        return self.zetaA * (a*self.dyda(a)-2*self.y(a))/a**3
        
    def s_crit(self, varpi, varpi_t, u=0.0):
        """
        Return the distances along of the line of sight where a given
        line of sight enters and exits the wind, or passes
        through the plane of the sky

        Parameters:
           varpi : float
              dimensionless impact parameter along the wind axis
           varpi_t : float
              dimensionless impact parameter transverse to the wind axis
           u : float
              velocity of interest; if u > 0, only distances on the
              far side wind are returned; if u < 0, only
              distances on the near side are returned; if u == 0.0,
              both far and near side distances are returned

        Returns:
           s_crit : array
              dimensionless locations of wind entry or exit, or
              passage through the midplane, ordered from smallest to
              largest value of s
        """
        s_crit = np.zeros(8)
        ns = libpwind.s_crit(varpi, varpi_t, u, self.__pw, s_crit)
        s_crit.resize(ns)
        return s_crit
 
    def a_crit(self, varpi, varpi_t, u=0.0):
        """
        Return the radii where a given line of sight enters and exits
        the wind, or passes through the plane of the sky

        Parameters:
           varpi : float
              dimensionless impact parameter along the wind axis
           varpi_t : float
              dimensionless impact parameter transverse to the wind axis
           u : float
              velocity of interest; if u > 0, only distances on the
              far side wind are returned; if u < 0, only
              distances on the near side are returned; if u == 0.0,
              both far and near side distances are returned

        Returns:
           a_crit : array(0), array(2) or array(4)
              dimensionless radii of wind entry or exit, or
              passage through the midplane, ordered from smallest to
              largest value of s
        """
        a_crit = np.zeros(8)
        na = libpwind.a_crit(varpi, varpi_t, u, self.__pw, a_crit)
        a_crit.resize(na)
        return a_crit

    def pdot(self, a, fg=None, tctw=None):
        """
        This returns the total momentum flux through the shell at
        radius a, normalised to the driving momentum flux.

        Parameters:
           a : float or array
              dimensionless radius
           fg : float or array
              gas fraction; must be in the range (0, 1]
           tctw : float or array
              ratio of the crossing time of the launch region to the
              wind evacuation time

        Returns:
           pdot : float or array
              momentum flux normalised to driving momentum flux

        Notes:
           If both fg and tctw are left as None, then the quantity
           returned will be the approximate wind momentum flux,
           derived taking the mass outflow rate to be Mdot = fg M_0 /
           zeta_M tc. If both are set, the result is the exact
           momentum flux. It is an error if one is None and the other
           is not.
        """
        bcast = np.broadcast(a, fg, tctw)
        pdot = np.zeros(bcast.shape)
        i = 0
        for (a_, fg_, tctw_) in bcast:
            if fg_ is None and tctw_ is None:
                pdot.flat[i] = libpwind.pdot_approx(a_, self.__pw)
            elif fg_ is not None and tctw_ is not None:
                pdot.flat[i] = libpwind.pdot_exact(a, fg, tctw,
                                                   self.__pw)
            else:
                raise ValueError(
                    "pwind.pdot: fg and tctw must both be None, or"
                    " neither be None")
            i = i+1
        if not hasattr(a, '__iter__') and \
           not hasattr(fg, '__iter__') and \
           not hasattr(tctw, '__iter__'):
            pdot = float(pdot)
        return pdot
    
    ################################################################
    # Functions that return observable quantities
    ################################################################

    def tau(self, u, tXtw=None, abd=None, Omega=None, wl=None,
            muH=1.4, tw=None, fj=1.0, boltzfac=0.0, u_trans=None,
            correlated=True, fw=None, varpi=0.0, varpi_t=0.0, a0=1.0,
            a1=np.finfo(c_double).max):
        """
        This returns the optical depth through the wind.

        Parameters:
           u : float or arraylike
              dimensionless line of sight velocity
           tXtw : float or arraylike
              ratio tX/tw for the wind
           abd : float or arraylike
              abundance of absorbers relative to H
           Omega : float or arraylike
              oscillator strength for transition
           wl : float or arraylike
              wavelength for transition, in cm
           muH : float or arraylike
              gas mass per H nucleus, in units of H masses
           tw : float or arraylike
              mass removal timescale, in sec
           fj : float or array
              fraction of the emitters in the lower state of the
              transition
           boltzfac : float or array
              Boltzmann factor exp(-E_ij/kB T_ex) for the two states of
              the transition, where T_ex = excitation temperature
           u_trans : arraylike
              velocity offset between transitions for multiple
              transition computations
           correlated : bool
              if True, assume correlated winds; if False, uncorrelated
           fw : float or arraylike
              covering factor of wind at launch point; if left as
              None, defaults to zeta_A; only used if correlated is True
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis
           a0 : float or arraylike
              minimum radius from which to integrate optical depth
           a1 : float or arraylike
              maximum radius from which to integrate optical depth

        Returns:
           tau : float or array
              optical depth at specified velocity

        Notes:
           The user can choose to specify either tXtw, or the
           combination of parameters (abd, Omega, wl, tw). In the
           latter case, the user can also optionally change the gas
           mass per H nucleus by setting muH. If tXtw is set, the
           other parameters will be ignored. If it is not set, then an
           error is raised if any of abd, Omega, wl, or tw is not
           specified.

           If u_trans is set, this changes the interpretation of the
           other variables. If u_trans is not None, then the trailing
           dimension of either tXtw or (Omega, wl) must have the same
           number of elements as u_trans, and will be interpreted as
           giving the value of tXtw, or equivalently the oscillator
           strength and wavelength, for each transition.
        """
        if correlated:
            return self.tau_c(u, tXtw=tXtw, abd=abd, Omega=Omega, wl=wl,
                              muH=muH, tw=tw, fj=fj, boltzfac=boltzfac,
                              u_trans=u_trans, fw=fw, varpi=varpi,
                              varpi_t=varpi_t, a0=a0, a1=a1)
        else:
            return self.tau_uc(u, tXtw=tXtw, abd=abd, Omega=Omega, wl=wl,
                               muH=muH, tw=tw, fj=fj, boltzfac=boltzfac,
                               u_trans=u_trans, varpi=varpi,
                               varpi_t=varpi_t, a0=a0, a1=a1)


    def tau_uc(self, u, tXtw=None, abd=None, Omega=None, wl=None,
               muH=1.4, tw=None, fj=1.0, boltzfac=0.0, u_trans=None,
               varpi=0.0, varpi_t=0.0, a0=1.0,
               a1=np.finfo(c_double).max):
        """
        This returns the optical depth assuming that the wind is
        uncorrelated.

        Parameters:
           u : float or arraylike
              dimensionless line of sight velocity
           tXtw : float or arraylike
              ratio tX/tw for the wind
           abd : float or arraylike
              abundance of absorbers relative to H
           Omega : float or arraylike
              oscillator strength for transition
           wl : float or arraylike
              wavelength for transition, in cm
           muH : float or arraylike
              gas mass per H nucleus, in units of H masses
           tw : float or arraylike
              mass removal timescale, in sec
           fj : float or array
              fraction of the emitters in the lower state of the
              transition
           boltzfac : float or array
              Boltzmann factor exp(-E_ij/kB T_ex) for the two states of
              the transition, where T_ex = excitation temperature
           u_trans : arraylike
              velocity offset between transitions for multiple
              transition computations
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis
           a0 : float or arraylike
              minimum radius from which to integrate optical depth
           a1 : float or arraylike
              maximum radius from which to integrate optical depth

        Returns:
           tau : float or array
              optical depth at specified velocity

        Notes:
           The user can choose to specify either tXtw, or the
           combination of parameters (abd, Omega, wl, tw). In the
           latter case, the user can also optionally change the gas
           mass per H nucleus by setting muH. If tXtw is set, the
           other parameters will be ignored. If it is not set, then an
           error is raised if any of abd, Omega, wl, or tw is not
           specified.

           If u_trans is set, this changes the interpretation of the
           other variables. If u_trans is not None, then the trailing
           dimension of either tXtw or (Omega, wl) must have the same
           number of elements as u_trans, and will be interpreted as
           giving the value of tXtw, or equivalently the oscillator
           strength and wavelength, for each transition.
        """
        # Check inputs, compute tXtw if necessary
        if tXtw is None:
            if (abd is None) or (Omega is None) or \
               (wl is None) or (tw is None):
                raise ValueError(
                    "pwind.tau_uc: must set either tXtw or "
                    "(abd, Omega, wl, tw)")
            tXtw = tX(np.asarray(abd), np.asarray(Omega),
                      np.asarray(wl), np.asarray(muH)) / \
                      np.asarray(tw)
        # Call c++ code to compute result; note the different
        # broadcasting depending on whether the trailing dimension of
        # tXtw is to be summed over transitions or not
        if u_trans is None:
            bcast = np.broadcast(np.asarray(u),
                                 np.asarray(tXtw),
                                 np.asarray(fj),
                                 np.asarray(boltzfac),
                                 np.asarray(varpi),
                                 np.asarray(varpi_t),
                                 np.asarray(a0),
                                 np.asarray(a1))
        else:
            # Need to do a bit of footwork here so that we can
            # broadcast over every dimension of tXtw except the final
            # one
            if np.asarray(tXtw).ndim == 1:
                tXtw_tmp = np.empty((1,), dtype=object)
                tXtw_tmp[0] = tXtw
            else:
                tXtw_tmp = np.empty(np.asarray(tXtw)[...,0].shape,
                                    dtype=object)
                for i in range(tXtw_tmp.size):
                    idx = np.unravel_index(i, tXtw_tmp.shape)
                    tXtw_tmp[idx] = np.asarray(tXtw)[idx]
            bcast = np.broadcast(np.asarray(u),
                                 tXtw_tmp,
                                 np.asarray(fj),
                                 np.asarray(boltzfac),
                                 np.asarray(varpi),
                                 np.asarray(varpi_t),
                                 np.asarray(a0),
                                 np.asarray(a1))
        tau_uc = np.zeros(bcast.shape)
        i = 0
        for (u_, tXtw_, fj_, boltzfac_, varpi_, varpi_t_,
             a0_, a1_) in bcast:
            if u_trans is None:
                tau_uc.flat[i] \
                    = libpwind.tau_uc(u_, tXtw_, fj_, boltzfac_,
                                      varpi_, varpi_t_, a0_, a1_,
                                      self.__pw)
            else:
                tau_uc.flat[i] \
                    = libpwind.tau_uc_vec(u_,
                                          np.array(u_trans,
                                                   dtype=c_double),
                                          np.array(tXtw_,
                                                   dtype=c_double),
                                          fj_, boltzfac_, len(u_trans),
                                          varpi_, varpi_t_, a0_, a1_,
                                          self.__pw)
            i = i+1
        if u_trans is None:
            if (not hasattr(u, '__iter__')) and \
               (not hasattr(tXtw, '__iter__')) and \
               (not hasattr(fj, '__iter__')) and \
               (not hasattr(boltzfac, '__iter__')) and \
               (not hasattr(varpi, '__iter__')) and \
               (not hasattr(varpi_t, '__iter__')) and \
               (not hasattr(a0, '__iter__')) and \
               (not hasattr(a1, '__iter__')):
                tau_uc = float(tau_uc)
        else:
            if (not hasattr(u, '__iter__')) and \
               (not hasattr(np.asarray(tXtw)[...,0], '__iter__')) and \
               (not hasattr(fj, '__iter__')) and \
               (not hasattr(boltzfac, '__iter__')) and \
               (not hasattr(varpi, '__iter__')) and \
               (not hasattr(varpi_t, '__iter__')) and \
               (not hasattr(a0, '__iter__')) and \
               (not hasattr(a1, '__iter__')):
                tau_c = float(tau_uc)
        return tau_uc
                
    def tau_c(self, u, tXtw=None, abd=None, Omega=None, wl=None,
              muH=1.4, tw=None, fw=None, fj=1.0, boltzfac=0.0,
              u_trans=None, varpi=0.0, varpi_t=0.0, a0=1.0,
              a1=np.finfo(c_double).max):
        """
        This returns the optical depth assuming that the wind is
        correlated.

        Parameters:
           u : float or arraylike
              dimensionless line of sight velocity
           tXtw : float or arraylike
              ratio tX/tw for the wind
           abd : float or arraylike
              abundance of absorbers relative to H
           Omega : float or arraylike
              oscillator strength for transition
           wl : float or arraylike
              wavelength for transition, in cm
           muH : float or arraylike
              gas mass per H nucleus, in units of H masses
           tw : float or arraylike
              mass removal timescale, in sec
           fw : float or arraylike
              covering factor of wind at launch point; if left as
              None, defaults to zeta_A
           fj : float or array
              fraction of the emitters in the lower state of the
              transition
           boltzfac : float or array
              Boltzmann factor exp(-E_ij/kB T_ex) for the two states of
              the transition, where T_ex = excitation temperature
           u_trans : arraylike
              velocity offset between transitions for multiple
              transition computations
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis
           a0 : float or arraylike
              minimum radius from which to integrate optical depth
           a1 : float or arraylike
              maximum radius from which to integrate optical depth

        Returns:
           tau : float or array
              optical depth at specified velocity

        Notes:
           The user can choose to specify either tXtw, or the
           combination of parameters (abd, Omega, wl, tw). In the
           latter case, the user can also optionally change the gas
           mass per H nucleus by setting muH. If tXtw is set, the
           other parameters will be ignored. If it is not set, then an
           error is raised if any of abd, Omega, wl, or tw is not
           specified.

           If u_trans is set, this changes the interpretation of the
           other variables. If u_trans is not None, then the trailing
           dimension of either tXtw or (Omega, wl) must have the same
           number of elements as u_trans, and will be interpreted as
           giving the value of tXtw, or equivalently the oscillator
           strength and wavelength, for each transition.
        """
        # Check inputs, compute tXtw if necessary
        if tXtw is None:
            if (abd is None) or (Omega is None) or \
               (wl is None) or (tw is None):
                raise ValueError(
                    "pwind.tau_uc: must set either tXtw or "
                    "(abd, Omega, wl, tw)")
            tXtw = tX(np.asarray(abd), np.asarray(Omega),
                      np.asarray(wl), np.asarray(muH)) / \
                      np.asarray(tw)
        # Call c++ code to compute result; note the different
        # broadcasting depending on whether the trailing dimension of
        # tXtw is to be summed over transitions or not
        if u_trans is None:
            bcast = np.broadcast(np.asarray(u), np.asarray(tXtw),
                                 np.asarray(fj), np.asarray(boltzfac),
                                 np.asarray(fw), np.asarray(varpi),
                                 np.asarray(varpi_t), np.asarray(a0),
                                 np.asarray(a1))
        else:
            # Need to do a bit of footwork here so that we can
            # broadcast over every dimension of tXtw except the final
            # one
            if np.asarray(tXtw).ndim == 1:
                tXtw_tmp = np.empty((1,), dtype=object)
                tXtw_tmp[0] = tXtw
            else:
                tXtw_tmp = np.empty(np.asarray(tXtw)[...,0].shape,
                                    dtype=object)
                for i in range(tXtw_tmp.size):
                    idx = np.unravel_index(i, tXtw_tmp.shape)
                    tXtw_tmp[idx] = np.asarray(tXtw)[idx]
            bcast = np.broadcast(np.asarray(u), tXtw_tmp,
                                 np.asarray(fj), np.asarray(boltzfac),
                                 np.asarray(fw), np.asarray(varpi),
                                 np.asarray(varpi_t), np.asarray(a0),
                                 np.asarray(a1))
        tau_c = np.zeros(bcast.shape)
        i = 0
        for (u_, tXtw_, fj_, boltzfac_, fw_, varpi_, varpi_t_,
             a0_, a1_) in bcast:
            if fw_ is None:
                fw_ = self.zetaA
            if u_trans is None:
                tau_c.flat[i] \
                    = libpwind.tau_c(u_, tXtw_, fj_, boltzfac_,
                                     fw_, varpi_, varpi_t_, a0_, a1_,
                                     self.__pw)
            else:
                tau_c.flat[i] \
                    = libpwind.tau_c_vec(u_,
                                         np.array(u_trans,
                                                  dtype=c_double),
                                         np.array(tXtw_,
                                                  dtype=c_double),
                                         fj_, boltzfac_, len(u_trans),
                                         fw_, varpi_, varpi_t_, a0_, a1_,
                                         self.__pw)
            i = i+1
        if u_trans is None:
            if (not hasattr(u, '__iter__')) and \
               (not hasattr(tXtw, '__iter__')) and \
               (not hasattr(fj, '__iter__')) and \
               (not hasattr(boltzfac, '__iter__')) and \
               (not hasattr(fw, '__iter__')) and \
               (not hasattr(varpi, '__iter__')) and \
               (not hasattr(varpi_t, '__iter__')) and \
               (not hasattr(a0, '__iter__')) and \
               (not hasattr(a1, '__iter__')):
                tau_c = float(tau_c)
        else:
            if (not hasattr(u, '__iter__')) and \
               (not hasattr(np.asarray(tXtw)[...,0], '__iter__')) and \
               (not hasattr(fj, '__iter__')) and \
               (not hasattr(boltzfac, '__iter__')) and \
               (not hasattr(fw, '__iter__')) and \
               (not hasattr(varpi, '__iter__')) and \
               (not hasattr(varpi_t, '__iter__')) and \
               (not hasattr(a0, '__iter__')) and \
               (not hasattr(a1, '__iter__')):
                tau_c = float(tau_c)
        return tau_c


    def Phi_c(self, u, fw=None, varpi=0.0, varpi_t=0.0, a0=1.0,
              a1=np.finfo(c_double).max):
        """
        Return the correlated absorption function Phi_c

        Parameters:
           u : float or arraylike
              dimensionless line of sight velocity
           fw : float or arraylike
              covering factor of wind at launch point; if left as
              None, defaults to zeta_A
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis
           a0 : float or arraylike
              minimum radius from which to integrate
           a1 : float or arraylike
              maximum radius from which to integrate

        Returns:
           Phi_c : float or array
              correlated absorption function
        """
        bcast = np.broadcast(np.asarray(u), np.asarray(fw),
                             np.asarray(varpi),
                             np.asarray(varpi_t),
                             np.asarray(a0),
                             np.asarray(a1))
        Phi_c = np.zeros(bcast.shape)
        i = 0
        for (u_, fw_, varpi_, varpi_t_, a0_, a1_) in bcast:
            if fw_ is None:
                fw_ = self.zetaA
            Phi_c.flat[i] = libpwind.Phi_c(u_, fw_, varpi_, varpi_t_,
                                           a0_, a1_, self.__pw)
            i = i+1
        if (not hasattr(u, '__iter__')) and \
           (not hasattr(fw, '__iter__')) and \
           (not hasattr(varpi, '__iter__')) and \
           (not hasattr(varpi_t, '__iter__')) and \
           (not hasattr(a0, '__iter__')) and \
           (not hasattr(a1, '__iter__')):
            Phi_c = float(Phi_c)
        return Phi_c
    
    def Phi_uc(self, u, varpi=0.0, varpi_t=0.0, a0=1.0,
               a1=np.finfo(c_double).max):
        """
        Return the uncorrelated absorption function Phi_c

        Parameters:
           u : float or arraylike
              dimensionless line of sight velocity
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind
              axis
           a0 : float or arraylike
              minimum radius from which to integrate Phi_uc
           a1 : float or arraylike
              maximum radius from which to integrate Phi_uc

        Returns:
           Phi_uc : float or array
              correlated absorption function
        """
        bcast = np.broadcast(np.asarray(u),
                             np.asarray(varpi),
                             np.asarray(varpi_t),
                             np.asarray(a0),
                             np.asarray(a1))
        Phi_uc = np.zeros(bcast.shape)
        i = 0
        for (u_, varpi_, varpi_t_, a0_, a1_) in bcast:
            Phi_uc.flat[i] = libpwind.Phi_uc(u_, varpi_, varpi_t_,
                                             a0_, a1_, self.__pw)
            i = i+1
        if (not hasattr(u, '__iter__')) and \
           (not hasattr(varpi, '__iter__')) and \
           (not hasattr(varpi_t, '__iter__')) and \
           (not hasattr(a0, '__iter__')) and \
           (not hasattr(a1, '__iter__')):
            Phi_uc = float(Phi_uc)
        return Phi_uc

    def Xi(self, u, varpi=0.0, varpi_t=0.0):
        """
        Return the optically thin emission line shape function
        Xi.

        Parameters:
           u : float or arraylike
              dimensionless line of sight velocity
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis

        Returns:
           Xi : float or array
              subcritical line emission shape function
        """
        bcast = np.broadcast(np.asarray(u), np.asarray(varpi),
                             np.asarray(varpi_t))
        Xi = np.zeros(bcast.shape)
        i = 0
        for (u_, varpi_, varpi_t_) in bcast:
            #print u_, varpi_, varpi_t_
            Xi.flat[i] = libpwind.Xi(u_, varpi_, varpi_t_, self.__pw)
            i = i+1
        if (not hasattr(u, '__iter__')) and \
           (not hasattr(varpi, '__iter__')) and \
           (not hasattr(varpi_t, '__iter__')):
            Xi = float(Xi)
        return Xi

    def xi(self, varpi=0.0, varpi_t=0.0):
        """
        Return the optically thin integrated intensity function xi.

        Parameters:
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis

        Returns:
           xi : float or array
              subcritical line emission intensity function
        """
        bcast = np.broadcast(np.asarray(varpi),
                             np.asarray(varpi_t))
        xi = np.zeros(bcast.shape)
        i = 0
        for (varpi_, varpi_t_) in bcast:
            xi.flat[i] = libpwind.xi(varpi_, varpi_t_, self.__pw)            
            i = i+1
        if (not hasattr(varpi, '__iter__')) and \
           (not hasattr(varpi_t, '__iter__')):
            xi = float(xi)
        return xi

    
    def eta(self, u, tXtw=None, abd=None, Omega=None, wl=None,
            muH=1.4, tw=None, fj=1.0, boltzfac=0.0, correlated=True,
            fw=None, varpi=0.0, varpi_t=0.0, thin=False):
        """
        Return the LTE emission profile function eta

        Parameters:
           u : float or array
              dimensionless line of sight velocity
           tXtw : float or array
              ratio of timescales tX and tw
           abd : float or arraylike
              abundance of absorbers relative to H
           Omega : float or arraylike
              oscillator strength for transition
           wl : float or arraylike
              wavelength for transition, in cm
           muH : float or arraylike
              gas mass per H nucleus, in units of H masses
           tw : float or arraylike
              mass removal timescale, in sec
           fj : float or array
              fraction of the emitters in the lower state of the
              transition
           boltzfac : float or array
              Boltzmann factor exp(-E_ij/kB T) for the two states of
              the transition
           correlated : bool
              if True, assume correlated winds; if False, uncorrelated
           fw : float or arraylike
              covering factor of wind at launch point; if left as
              None, defaults to zeta_A; only used if correlated is True
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis
           thin : bool or array of bool
              if True, the escape probability is set to 1

        Returns:
           eta : float or array
              the LTE emission function eta

        Notes:
           The user can choose to specify either tXtw, or the
           combination of parameters (abd, Omega, wl, tw). In the
           latter case, the user can also optionally change the gas
           mass per H nucleus by setting muH. If tXtw is set, the
           other parameters will be ignored. If it is not set, then an
           error is raised if any of abd, Omega, wl, or tw is not
           specified.
        """
        
        # Check inputs, compute tXtw if necessary
        if tXtw is None:
            if (abd is None) or (Omega is None) or \
               (wl is None) or (tw is None):
                raise ValueError(
                    "pwind.eta: must set either tXtw or "
                    "(abd, Omega, wl, tw)")
            tXtw = tX(np.asarray(abd), np.asarray(Omega),
                      np.asarray(wl), np.asarray(muH)) / np.asarray(tw)
        # Call c++ code to compute result
        bcast = np.broadcast(np.asarray(u), np.asarray(tXtw),
                             np.asarray(fj), np.asarray(boltzfac),
                             np.asarray(correlated), np.asarray(fw),
                             np.asarray(varpi), np.asarray(varpi_t),
                             np.asarray(thin))
        eta = np.zeros(bcast.shape)
        i = 0
        for (u_, tXtw_, fj_, boltzfac_, correlated_, fw_, varpi_,
             varpi_t_, thin_) in bcast:
            if fw_ is None:
                fw_ = self.zetaA
            eta.flat[i] = libpwind.eta(u_, tXtw_, fj_, boltzfac_,
                                       correlated_, fw_, varpi_,
                                       varpi_t_, thin_, self.__pw)
            i = i+1
        if (not hasattr(u, '__iter__')) and \
           (not hasattr(tXtw, '__iter__')) and \
           (not hasattr(fj, '__iter__')) and \
           (not hasattr(boltzfac, '__iter__')) and \
           (not hasattr(correlated, '__iter__')) and \
           (not hasattr(fw, '__iter__')) and \
           (not hasattr(varpi, '__iter__')) and \
           (not hasattr(varpi_t, '__iter__')) and \
           (not hasattr(thin, '__iter__')):
            eta = float(eta)
        return eta

    def temp_LTE(self, u, T,
                 # Emitter data
                 emit=None,
                 # Parameters to determine tXtw
                 tw=None, tXtw=None, abd=None, Omega=None, wl=None,
                 muH=1.4, 
                 # Parameters to determine level populations
                 fj=1.0, boltzfac=0.0, trans=0,
                 # Geometric and radiation transfer parameters
                 varpi=0.0, varpi_t=0.0, correlated=True, fw=None,
                 thin=False,
                 # Output parameters
                 TA=False):
        """
        Return the brightness or antenna temperature of line emission
        for a species in LTE observed at velocity u

        Parameters:
           u : float or array
              dimensionless line of sight velocity
           T : float or array
              wind kinetic temperature, in K
           emit : emitter or sequence of emitter
              an emitter object for the transition of interest
           tw : float or arraylike
              mass removal timescale, in sec
           tXtw : float or array
              ratio of timescales tX and tw
           abd : float or arraylike
              abundance of absorbers relative to H
           Omega : float or arraylike
              oscillator strength for transition
           wl : float or arraylike
              wavelength for transition, in cm
           muH : float or arraylike
              gas mass per H nucleus, in units of H masses
           fj : float or array
              fraction of the emitters in the lower state of the
              transition
           boltzfac : float or array
              Boltzmann factor exp(-E_ij/kB T) for the two states of
              the transition
           trans : int or array of int
              transition for which to compute the emission
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis
           correlated : bool
              if True, assume correlated winds; if False, uncorrelated
           fw : float or arraylike
              covering factor of wind at launch point; if left as
           thin : bool or array of bool
              if True, the escape probability is set to 1
           TA : bool
              if True, quantity returned is antenna temperature;
              default is to return brightness temperature

        Returns:
           Temp : float or array
              brightness temperature (or antenna temperature if TA is
              set to True) in K as a function of the velocity u

        Note:
           The user can provide the required input in the following
           ways, which are checked in order.

           1. Set emit and tw, and optionally trans; in this case
              tXtw, abd, Omega, wl, fj, and boltzfac are all ignored;
              the abundance and all properties of the line are taken
              from emit, and the Boltzmann factors are computed
              directly from the emitter object and the input
              temperature.
           2. Set tXtw and wl, and optionally fj and boltzfac; in this
              case emit, tw, abd, and Omega are ignored.
           3. Set abd, Omega, wl, and tw, and optionally fj and
              boltzfac; in this case tXtw is computed from (abd,
              Omega, wl, tw).

           It is an error if not enough parameters are specified to
           satisfy conditions 1, 2, or 3.
        """
        # Check inputs and compute derived quantities if needed
        if emit is not None and tw is not None:
            # Get tXtw, fj, boltzfac from emit and tw
            bcast = np.broadcast(np.asarray(T), np.asarray(emit),
                                 np.asarray(tw), np.asarray(muH),
                                 np.asarray(trans))
            tXtw = np.zeros(bcast.shape)
            fj = np.zeros(bcast.shape)
            boltzfac = np.zeros(bcast.shape)
            wl = np.zeros(bcast.shape)
            i = 0
            for (T_, emit_, tw_, muH_, trans_) in bcast:
                # Get tXtw
                tXtw.flat[i] = emit_.tX(muH_*mH / emit_.abundance,
                                        trans=trans_) / tw_
                # Get wavelength and Boltzmann factor
                up = emit_.data.radUpper[trans_]
                lo = emit_.data.radLower[trans_]
                Eij = (emit_.data.levEnergy[up] -
                       emit_.data.levEnergy[lo])
                wl.flat[i] = h*c / Eij
                boltzfac.flat[i] = np.exp(-Eij / (kB*T_))
                # Get fj
                emit_.setLevPopLTE(T_)
                fj.flat[i] = emit_.levPop[lo]
                # Increment
                i = i+1
        elif tXtw is None:
            # Get tXtw from abd, Omega, wl, tw
            if (abd is None) or (Omega is None) or \
               (wl is None) or (tw is None):
                raise ValueError(
                    "pwind.temp_LTE: insufficient inputs specified")
            tXtw = tX(np.asarray(abd), np.asarray(Omega),
                      np.asarray(wl), np.asarray(muH)) / \
                      np.asarray(tw)
        elif wl is None:
            raise ValueError(
                "pwind.temp_LTE: insufficient inputs specified")

        # Now that necessary inputs have been derived, compute eta
        eta = self.eta(u, tXtw=tXtw, fj=fj, boltzfac=boltzfac,
                       correlated=correlated, fw=fw, varpi=varpi,
                       varpi_t=varpi_t, thin=thin)
                
        # Compute I_nu / B_nu
        q = tXtw * eta * fj * (1.0 - boltzfac)

        # Convert to TB or TA as requested
        if TA:
            temp = (h*c/(kB*wl)) / \
                   np.log((q - 1.0 + 1.0/boltzfac)/q)
        else:
            temp = (h*c/(kB*wl)) * q / (1.0/boltzfac - 1.0)
            
        # Return output
        return np.squeeze(temp)
    

    def Psi(self, tXtw, fj, boltzfac, correlated=True,
            fw=None, varpi=0.0, varpi_t=0.0, thin=False):
        """
        Return the integrated intensity integral Psi

        Parameters:
           tXtw : float or array
              ratio of timescales tX and tw
           fj : float or array
              fraction of the emitters in the lower state of the
              transition
           boltzfac : float or array
              Boltzmann factor exp(-E_ij/kB T) for the two states of
              the transition
           correlated : bool
              if True, assume correlated winds; if False, uncorrelated
           fw : float or arraylike
              covering factor of wind at launch point; if left as
              None, defaults to zeta_A; only used if correlated is True
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis
           thin : bool or array of bool
              if True, the escape probability is set to 1

        Returns:
           Psi : float or array
              value of the integral Psi
        """
        bcast = np.broadcast(np.asarray(tXtw), np.asarray(fj),
                             np.asarray(boltzfac), np.asarray(correlated),
                             np.asarray(fw), np.asarray(varpi),
                             np.asarray(varpi_t), np.asarray(thin))
        Psi = np.zeros(bcast.shape)
        i = 0
        for (tXtw_, fj_, boltzfac_, correlated_, fw_,
             varpi_, varpi_t_, thin_) in bcast:
            if fw_ is None:
                fw_ = self.zetaA
            print i
            Psi.flat[i] = libpwind.Psi(tXtw_, fj_, boltzfac_,
                                       correlated_, fw_,
                                       varpi_, varpi_t_, thin_,
                                       self.__pw)
            print i, Psi.flat[i]
            i = i+1
        if (not hasattr(tXtw, '__iter__')) and \
           (not hasattr(fj, '__iter__')) and \
           (not hasattr(boltzfac, '__iter__')) and \
           (not hasattr(varpi, '__iter__')) and \
           (not hasattr(varpi_t, '__iter__')) and \
           (not hasattr(thin, '__iter__')):
            Psi = float(Psi)
        return Psi

    def intTA_LTE(self, v0, T,
                  # Emitter data
                  emit=None,
                  # Parameters to determine tXtw
                  tw=None, tXtw=None, abd=None, Omega=None, wl=None,
                  muH=1.4, 
                  # Parameters to determine level populations
                  fj=1.0, boltzfac=0.0, trans=0,
                  # Geometric and radiation transfer parameters
                  varpi=0.0, varpi_t=0.0, correlated=True, fw=None,
                  thin=False):

        """
        Returns the velocity-integrated antenna temperature along 
        a particular line of sight

        Parameters:
           v0 : float or array
              escape speed of system; units are arbitrary, and output
              integrated antenna temperature will be in the same
              velocity units as v0
           T : float or array
              wind kinetic temperature, in K
           emit : emitter or sequence of emitter
              an emitter object for the transition of interest
           tw : float or arraylike
              mass removal timescale, in sec
           tXtw : float or array
              ratio of timescales tX and tw
           abd : float or arraylike
              abundance of absorbers relative to H
           Omega : float or arraylike
              oscillator strength for transition
           wl : float or arraylike
              wavelength for transition, in cm
           muH : float or arraylike
              gas mass per H nucleus, in units of H masses
           fj : float or array
              fraction of the emitters in the lower state of the
              transition
           boltzfac : float or array
              Boltzmann factor exp(-E_ij/kB T) for the two states of
              the transition
           trans : int or array of int
              transition for which to compute the emission
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis
           correlated : bool
              if True, assume correlated winds; if False, uncorrelated
           fw : float or arraylike
              covering factor of wind at launch point; if left as
           thin : bool or array of bool
              if True, the escape probability is set to 1

        Returns:
           vT : float or array
              velocity-integrated antenna temperature, in units of K
              times whatever units v0 is in (i.e., if v0 is in km / s,
              then the returned value will be in units of K km / s)

        Note:
           The user can provide the required input in the following
           ways, which are checked in order.

           1. Set emit and tw, and optionally trans; in this case
              tXtw, abd, Omega, wl, fj, and boltzfac are all ignored;
              the abundance and all properties of the line are taken
              from emit, and the Boltzmann factors are computed
              directly from the emitter object and the input
              temperature.
           2. Set tXtw and wl, and optionally fj and boltzfac; in this
              case emit, tw, abd, and Omega are ignored.
           3. Set abd, Omega, wl, and tw, and optionally fj and
              boltzfac; in this case tXtw is computed from (abd,
              Omega, wl, tw).

           It is an error if not enough parameters are specified to
           satisfy conditions 1, 2, or 3.
        """
        # Check inputs and compute derived quantities if needed
        if emit is not None and tw is not None:
            # Get tXtw, fj, boltzfac from emit and tw
            bcast = np.broadcast(np.asarray(T), np.asarray(emit),
                                 np.asarray(tw), np.asarray(muH),
                                 np.asarray(trans))
            tXtw = np.zeros(bcast.shape)
            fj = np.zeros(bcast.shape)
            boltzfac = np.zeros(bcast.shape)
            wl = np.zeros(bcast.shape)
            i = 0
            for (T_, emit_, tw_, muH_, trans_) in bcast:
                # Get tXtw
                tXtw.flat[i] = emit_.tX(muH_*mH / emit_.abundance,
                                        trans=trans_) / tw_
                # Get wavelength and Boltzmann factor
                up = emit_.data.radUpper[trans_]
                lo = emit_.data.radLower[trans_]
                Eij = (emit_.data.levEnergy[up] -
                       emit_.data.levEnergy[lo])
                wl.flat[i] = h*c / Eij
                boltzfac.flat[i] = np.exp(-Eij / (kB*T_))
                # Get fj
                emit_.setLevPopLTE(T_)
                fj.flat[i] = emit_.levPop[lo]
                # Increment
                i = i+1
        elif tXtw is None:
            # Get tXtw from abd, Omega, wl, tw
            if (abd is None) or (Omega is None) or \
               (wl is None) or (tw is None):
                raise ValueError(
                    "pwind.temp_LTE: insufficient inputs specified")
            tXtw = tX(np.asarray(abd), np.asarray(Omega),
                      np.asarray(wl), np.asarray(muH)) / \
                      np.asarray(tw)
        elif wl is None:
            raise ValueError(
                "pwind.intTA_LTE: insufficient inputs specified")

        # Get Psi
        Psi = self.Psi(tXtw, fj, boltzfac, correlated=correlated,
                       fw=fw, varpi=varpi, varpi_t=varpi_t,
                       thin=thin)

        # Now multiply by prefactor
        intTA = fj * tXtw * Psi * v0/wl

        # Return output
        return np.squeeze(intTA)


    def Xfac(self, T, emit, tw, muH=1.4, trans=0,
             varpi=0.0, varpi_t=0.0, correlated=True, fw=None):
        """
        Return the X factor -- the conversion between mass and
        velocity-integrated antenna temperature.

        Parameters:
           T : float or array
              wind kinetic temperature, in K
           emit : emitter or sequence of emitter
              an emitter object for the transition of interest
           tw : float or arraylike
              mass removal timescale, in sec
           trans : int or array of int
              transition for which to compute the emission
           varpi : float or arraylike
              dimensionless impact parameter along the wind axis
           varpi_t : float or arraylike
              dimensionless impact parameter transverse to the wind axis
           correlated : bool
              if True, assume correlated winds; if False, uncorrelated
           fw : float or arraylike
              covering factor of wind at launch point; if left as

        Returns:
           Xfac : float or array
              X factor, in the same units as Xthin
        """
        # Get tXtw, fj, boltzfac from emit and tw
        bcast = np.broadcast(np.asarray(T), np.asarray(emit),
                             np.asarray(tw), np.asarray(muH),
                             np.asarray(trans))
        tXtw = np.zeros(bcast.shape)
        fi = np.zeros(bcast.shape)
        fj = np.zeros(bcast.shape)
        boltzfac = np.zeros(bcast.shape)
        i = 0
        for (T_, emit_, tw_, muH_, trans_) in bcast:
            # Get tXtw
            tXtw.flat[i] = emit_.tX(muH_*mH / emit_.abundance,
                                    trans=trans_) / tw_
            # Get wavelength and Boltzmann factor
            up = emit_.data.radUpper[trans_]
            lo = emit_.data.radLower[trans_]
            Eij = (emit_.data.levEnergy[up] -
                   emit_.data.levEnergy[lo])
            boltzfac.flat[i] = np.exp(-Eij / (kB*T_))
            # Get fi and fj
            emit_.setLevPopLTE(T_)
            fi.flat[i] = emit_.levPop[up]
            fj.flat[i] = emit_.levPop[lo]
            # Increment
            i = i+1

        # Get Psi and Psi(thin)
        Psi = self.Psi(tXtw, fj, boltzfac, correlated=correlated,
                       fw=fw, varpi=varpi, varpi_t=varpi_t)
        Psi_thin = self.Psi(tXtw, fj, boltzfac, correlated=correlated,
                            fw=fw, varpi=varpi, varpi_t=varpi_t,
                            thin=True)
        
        # Get final result and return
        X = emit.Xthin(trans=trans) * Psi_thin / (fi * Psi)
        return X
