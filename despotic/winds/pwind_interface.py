"""
This module loads the interface to the pwind c++ code
"""

import numpy as np
import os
import os.path as osp
import numpy.ctypeslib as npct
from ctypes import c_double, c_void_p, c_bool, c_ulong, c_int, c_char_p

# Type definition
array_1d_double = npct.ndpointer(dtype=np.double, ndim=1,
                                 flags="CONTIGUOUS")

# Load the library, compiling it if required
libdir = osp.dirname(osp.realpath(__file__))
try:
    libpwind = npct.load_library("libpwind", libdir)
except OSError:
    print("Failed to load libpwind... trying to compile it")
    os.system("cd "+libdir+"; make")
    libpwind = npct.load_library("libpwind", libdir)

# Allocation / de-allocation methods

# Geometries
libpwind.pwind_geom_sphere_new.restype = c_void_p
libpwind.pwind_geom_sphere_new.argtypes = [ c_double ]
libpwind.pwind_geom_cone_new.restype = c_void_p
libpwind.pwind_geom_cone_new.argtypes = [ c_double, c_double ]
libpwind.pwind_geom_cone_sheath_new.restype = c_void_p
libpwind.pwind_geom_cone_sheath_new.argtypes \
    = [ c_double, c_double, c_double ]

# Ideal winds
libpwind.pwind_ideal_pa_new.restype = c_void_p
libpwind.pwind_ideal_pa_new.argtypes \
    = [ c_double, c_double, c_void_p, c_double, c_double ]
libpwind.pwind_ideal_pi_new.restype = c_void_p
libpwind.pwind_ideal_pi_new.argtypes \
    = [ c_double, c_double, c_void_p, c_double, c_double ]
libpwind.pwind_ideal_ps_new.restype = c_void_p
libpwind.pwind_ideal_ps_new.argtypes \
    = [ c_double, c_double, c_void_p, c_double, c_double ]
libpwind.pwind_ideal_ia_new.restype = c_void_p
libpwind.pwind_ideal_ia_new.argtypes \
    = [ c_double, c_double, c_void_p, c_double, c_double ]
libpwind.pwind_ideal_ii_new.restype = c_void_p
libpwind.pwind_ideal_ii_new.argtypes \
    = [ c_double, c_double, c_void_p, c_double, c_double ]
libpwind.pwind_ideal_is_new.restype = c_void_p
libpwind.pwind_ideal_is_new.argtypes \
    = [ c_double, c_double, c_void_p, c_double, c_double ]

# Radiation-driven winds
libpwind.pwind_rad_pa_new.restype = c_void_p
libpwind.pwind_rad_pa_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double ]
libpwind.pwind_rad_pi_new.restype = c_void_p
libpwind.pwind_rad_pi_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double ]
libpwind.pwind_rad_ps_new.restype = c_void_p
libpwind.pwind_rad_ps_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double ]
libpwind.pwind_rad_ia_new.restype = c_void_p
libpwind.pwind_rad_ia_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double ]
libpwind.pwind_rad_ii_new.restype = c_void_p
libpwind.pwind_rad_ii_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double ]
libpwind.pwind_rad_is_new.restype = c_void_p
libpwind.pwind_rad_is_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double ]

# Hot gas-driven winds
libpwind.pwind_hot_pa_new.restype = c_void_p
libpwind.pwind_hot_pa_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double, c_void_p ]
libpwind.pwind_hot_pi_new.restype = c_void_p
libpwind.pwind_hot_pi_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double, c_void_p ]
libpwind.pwind_hot_ps_new.restype = c_void_p
libpwind.pwind_hot_ps_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double, c_void_p ]
libpwind.pwind_hot_ia_new.restype = c_void_p
libpwind.pwind_hot_ia_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double, c_void_p ]
libpwind.pwind_hot_ii_new.restype = c_void_p
libpwind.pwind_hot_ii_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double, c_void_p ]
libpwind.pwind_hot_is_new.restype = c_void_p
libpwind.pwind_hot_is_new.argtypes \
    = [ c_double, c_double, c_double, c_void_p,
        c_double, c_double, c_void_p ]

# Free
libpwind.pwind_free.restype = None
libpwind.pwind_free.argtypes = [ c_void_p ]
libpwind.pwind_geom_free.restype = None
libpwind.pwind_geom_free.argtypes = [ c_void_p ]

# Information return methods
libpwind.Gamma.restype = c_double
libpwind.Gamma.argtypes = [ c_void_p ]
libpwind.mach.restype = c_double
libpwind.mach.argtypes = [ c_void_p ]
libpwind.fcrit.restype = c_double
libpwind.fcrit.argtypes = [ c_void_p ]
libpwind.jsp.restype = c_double
libpwind.jsp.argtypes = [ c_void_p ]
libpwind.xcrit.restype = c_double
libpwind.xcrit.argtypes = [ c_void_p ]
libpwind.sx.restype = c_double
libpwind.sx.argtypes = [ c_void_p ]
libpwind.zetaM.restype = c_double
libpwind.zetaM.argtypes = [ c_void_p ]
libpwind.zetaA.restype = c_double
libpwind.zetaA.argtypes = [ c_void_p ]
libpwind.umax.restype = c_double
libpwind.umax.argtypes = [ c_void_p ]

# Simple setting methods
libpwind.set_mach.restype = None
libpwind.set_mach.argtypes = [ c_double, c_void_p ]
libpwind.set_fcrit.restype = None
libpwind.set_fcrit.argtypes = [ c_double, c_void_p ]
libpwind.set_jsp.restype = None
libpwind.set_jsp.argtypes = [ c_double, c_void_p ]
libpwind.set_geometry.restype = None
libpwind.set_geometry.argtypes = [ c_void_p, c_void_p ]

# Kinematics methods
libpwind.y.restype = c_double
libpwind.y.argtypes = [ c_double, c_void_p ]
libpwind.m.restype = c_double
libpwind.m.argtypes = [ c_double, c_void_p ]
libpwind.dyda.restype = c_double
libpwind.dyda.argtypes = [ c_double, c_void_p ]
libpwind.X.restype = c_double
libpwind.X.argtypes = [ c_double, c_double, c_void_p ]
libpwind.U2.restype = c_double
libpwind.U2.argtypes = [ c_double, c_double, c_void_p ]
libpwind.dU2dx.restype = c_double
libpwind.dU2dx.argtypes = [ c_double, c_double, c_void_p ]
libpwind.dU2da.restype = c_double
libpwind.dU2da.argtypes = [ c_double, c_double, c_void_p ]

# Limit methods
libpwind.alimits.restype = c_ulong
libpwind.alimits.argtypes = [ c_double, c_double, c_double,
                              c_double, c_double, c_void_p,
                              array_1d_double ]
libpwind.xlimits.restype = c_ulong
libpwind.xlimits.argtypes = [ c_double, c_void_p,
                              array_1d_double ]
libpwind.amax.restype = c_double
libpwind.amax.argtypes = [ c_double, c_double, c_double, c_void_p ]
libpwind.amax_abs.restype = c_double
libpwind.amax_abs.argtypes = [ c_void_p ]
libpwind.s_crit.restype = c_ulong
libpwind.s_crit.argtypes = [ c_double, c_double, c_double, c_void_p,
                             array_1d_double ]
libpwind.a_crit.restype = c_ulong
libpwind.a_crit.argtypes = [ c_double, c_double, c_double, c_void_p,
                             array_1d_double ]

# Computation methods
libpwind.drhodx.restype = c_double
libpwind.drhodx.argtypes = [ c_double, c_double, c_void_p ]
libpwind.rho.restype = c_double
libpwind.rho.argtypes = [ c_double, c_double, c_double, c_void_p ]
libpwind.rho_vec.restype = None
libpwind.rho_vec.argtypes = [ c_ulong, array_1d_double, c_double, c_double,
                              c_void_p, array_1d_double ]
libpwind.pdot.restype = c_double
libpwind.pdot.argtypes = [ c_double, c_double, c_double, c_void_p ]
libpwind.pdot_vec.restype = None
libpwind.pdot_vec.argtypes = [ c_ulong, array_1d_double, c_double, c_double,
                              c_void_p, array_1d_double ]
libpwind.Edot.restype = c_double
libpwind.Edot.argtypes = [ c_double, c_double, c_double, c_void_p ]
libpwind.Edot_vec.restype = None
libpwind.Edot_vec.argtypes = [ c_ulong, array_1d_double, c_double, c_double,
                              c_void_p, array_1d_double ]
libpwind.pdotRel_approx.restype = c_double
libpwind.pdotRel_approx.argtypes = [ c_double, c_double, c_double, c_void_p ]
libpwind.pdotRel_exact.restype = c_double
libpwind.pdotRel_exact.argtypes = [ c_double, c_double, c_double,
                                 c_double, c_double, c_void_p ]
libpwind.Phi_uc.restype = c_double
libpwind.Phi_uc.argtypes = [ c_double, c_double, c_double, c_double,
                             c_double, c_double, c_void_p ]
libpwind.Phi_uc_vec.restype = None
libpwind.Phi_uc_vec.argtypes = [ c_ulong, array_1d_double,
                                 c_double, c_double, c_double, c_double,
                                 c_double, c_double, c_void_p,
                                 array_1d_double]
libpwind.Phi_c.restype = c_double
libpwind.Phi_c.argtypes = [ c_double, c_double, c_double, c_double,
                            c_double, c_double, c_double, c_double, c_void_p ]
libpwind.Phi_c_vec.restype = None
libpwind.Phi_c_vec.argtypes = [ c_ulong, array_1d_double, c_double,
                                c_double, c_double, c_double, c_double,
                                c_double, c_double, c_void_p,
                                array_1d_double]
libpwind.tau_uc.restype = c_double
libpwind.tau_uc.argtypes = [ c_double, c_double, c_double, c_double,
                             c_double, c_double, c_double, c_double,
                             c_double, c_double, c_void_p ]
libpwind.tau_uc_vec.restype = None
libpwind.tau_uc_vec.argtypes = [ c_ulong, array_1d_double,
                                 c_double, c_double, c_double,
                                 c_double, c_double, c_double, c_double,
                                 c_double, c_double, c_void_p,
                                 array_1d_double ]
libpwind.tau_uc_multiple.restype = c_double
libpwind.tau_uc_multiple.argtypes = [ c_double, array_1d_double,
                                      array_1d_double,
                                      c_double, c_double,
                                      c_ulong, c_double, c_double, c_double,
                                      c_double, c_double, c_double, c_void_p ]
libpwind.tau_uc_multiple_vec.restype = None
libpwind.tau_uc_multiple_vec.argtypes = [ c_ulong, array_1d_double,
                                          array_1d_double,
                                          array_1d_double,
                                          c_double, c_double,
                                          c_ulong, c_double, c_double,
                                          c_double, c_double, c_double,
                                          c_double, c_void_p,
                                          array_1d_double ]
libpwind.tau_c.restype = c_double
libpwind.tau_c.argtypes = [ c_double, c_double, c_double, c_double, c_double,
                            c_double, c_double, c_double, c_double, c_double,
                            c_double, c_void_p ]
libpwind.tau_c_vec.restype = None
libpwind.tau_c_vec.argtypes = [ c_ulong, array_1d_double,
                                c_double, c_double, c_double, c_double,
                                c_double, c_double, c_double, c_double,
                                c_double, c_double, c_void_p,
                                array_1d_double ]
libpwind.tau_c_multiple.restype = c_double
libpwind.tau_c_multiple.argtypes = [ c_double, array_1d_double,
                                     array_1d_double,
                                     c_double, c_double,
                                     c_ulong, c_double, c_double, c_double,
                                     c_double, c_double,
                                     c_double, c_double, c_void_p ]
libpwind.tau_c_multiple_vec.restype = None
libpwind.tau_c_multiple_vec.argtypes = [ c_ulong, array_1d_double,
                                         array_1d_double,
                                         array_1d_double,
                                         c_double, c_double,
                                         c_ulong, c_double, c_double, c_double,
                                         c_double, c_double,
                                         c_double, c_double, c_void_p,
                                         array_1d_double ]
libpwind.Xi.restype = c_double
libpwind.Xi.argtypes = [ c_double, c_double, c_double,
                         c_double, c_double, c_void_p ]
libpwind.Xi_vec.restype = None
libpwind.Xi_vec.argtypes = [ c_ulong, array_1d_double, c_double, c_double,
                             c_double, c_double, c_void_p, array_1d_double ]
libpwind.xi.restype = c_double
libpwind.xi.argtypes = [ c_double, c_double, c_void_p ]
libpwind.eta.restype = c_double
libpwind.eta.argtypes = [ c_double, c_double, c_double, c_double, c_bool,
                          c_double, c_double, c_double, c_bool,
                          c_double, c_double, c_void_p ]
libpwind.eta_vec.restype = None
libpwind.eta_vec.argtypes = [ c_ulong, array_1d_double, c_double, c_double,
                              c_double, c_bool, c_double, c_double, c_double,
                              c_bool, c_double, c_double, c_void_p,
                              array_1d_double ]
libpwind.Psi.restype = c_double
libpwind.Psi.argtypes = [ c_double, c_double, c_double, c_bool,
                          c_double, c_double, c_double, c_bool,
                          c_double, c_double, c_void_p ]

# Error management methods
libpwind.get_err.restype = c_int
libpwind.get_err.argtypes = [ c_void_p ]
libpwind.set_err.restype = None
libpwind.set_err.argtypes = [ c_int, c_void_p ]
libpwind.clear_err.restype = None
libpwind.clear_err.argtypes = [ c_void_p ]
libpwind.get_err_str.restype = c_char_p
libpwind.get_err_str.argtypes = [ c_void_p ]

# Methods to manage pre-tabulated hot gas data
libpwind.read_hot_wind_table.restype = c_void_p
libpwind.read_hot_wind_table.argtypes = [ c_char_p, c_int, c_int ]
libpwind.free_hot_wind_table.restype = None
libpwind.free_hot_wind_table.argtypes = [ c_void_p ]
libpwind.get_hot_wind_table_limits.restype = None
libpwind.get_hot_wind_table_limits.argtypes = [ c_void_p, array_1d_double ]
