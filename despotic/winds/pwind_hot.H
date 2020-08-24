// A class representing hot gas-driven winds

#ifndef _PWIND_HOT_H_
#define _PWIND_HOT_H_

#include "pwind.H"
#include "pwind_expansion.H"
#include "pwind_geom.H"
#include "pwind_potential.H"
#include <cmath>
#include <limits>
#include <vector>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

////////////////////////////////////////////////////////////////////////
// Little helpful data structure
////////////////////////////////////////////////////////////////////////
struct duda_params {
  double x;
  const void *w;
};

////////////////////////////////////////////////////////////////////////
// The generic class
////////////////////////////////////////////////////////////////////////  
class pwind_hot : public pwind {

public:

  // The constructor
  pwind_hot(const double Gamma_, const double mach_, const double uh_,
	    const double yoverm_,
	    const pwind_potential *potential_,
	    const pwind_expansion *expansion_,
	    const pwind_geom *geom_ =
	    static_cast<const pwind_geom *>(&pw_geom_sphere),
	    const double epsabs_ = 1.0e-6, const double epsrel_ = 1.0e-4,
	    const double interpabs_ = 1.0e-2,
	    const double interprel_ = 1.0e-2,
	    const double fcrit_ = 1.0,
	    const double jsp_ = 0.0,
	    const double amax_grid = 1.0e10);

  // The destructor
  virtual ~pwind_hot();

  // Routines needed to construct the interpolation table
  std::vector<double>
  Ua_solve_tab(const double x, const std::vector<double>& loga,
	       double& loga_stop,
	       const std::vector<double>::size_type istart = 0,
	       const std::vector<double>::size_type iend_ = 0,
	       const double u0 = 0.0) const;
  std::vector<double>
  Ua_solve_auto(const double x, std::vector<double>& loga,
		double& loga_stop) const;
  void add_table_row(const double x);
  void build_table(const std::vector<double>& x_grid_init);
  void build_interpolators();
  void init(const double amax_init);

  // Informational routines
  double getUh() const { return uh; }
  double getXMin() const { return x_grid[0]; }
  double getAMax() const { return exp(loga_grid.back()); }
  double getLogAMax() const { return loga_grid.back(); }

  // Wind kinematics; these are attached to the generic parent class
  // because all these evaluations are interpolations off a table
  virtual double X(const double ur, const double a) const;
  virtual double U2(const double x, const double a) const;
  virtual double dU2dx(const double x, const double a) const;
  virtual double dU2da(const double x, const double a) const;

  // Limits on the range of the wind
  virtual double x_stop(const double a) const;
  virtual double amax(const double x) const;
  virtual std::vector<double> alimits(const double u,
				      const double varpi = 0.0,
				      const double varpi_t = 0.0) const;

protected:

  const double uh;        // Hot gas velocity
  const double yoverm;    // Limit of y/m as a->inf
  double interpabs, interprel; // Error tolerances on interpolation
  double amax_grid; // Absolute upper limit on extent of grid
  double uMin_inf;  // Velocity of gas with x = xcrit as a->infinity
  std::vector<double> loga_grid;  // Grid of points in log radius
  std::vector<double> x_grid;     // Grid of points in x
  std::vector<std::vector<double> > u_grid; // Grid of velocities
  std::vector<double> loga_stop_grid;  // Grid of stopping radii

  // Interpolators and accelerators; note that the accelerators are
  // declared as mutable, since they are modified when used by the
  // const functions
  gsl_spline2d *U2_interp;    // Interpolator for U2(x, a)
  gsl_spline2d *dU2dx_interp; // Interpolator for dU2dx(x, a)
  gsl_spline2d *X_interp ;    // Interpolator for X(u, a)
  gsl_spline *Umin_interp;    // Interpolator for U_min(a) = U(xcrit, a)
  gsl_spline *Umax_interp;    // Interpolator for U_max(a) = U(xmin, a)
  gsl_spline *aMin_interp[2]; // Interpolator for a_min(u) = a(xmin, u)
  gsl_spline *aMax_interp;    // Interpolator for a_min(u) = a(xcrit, u)
  gsl_spline *a_stop_interp;  // Interpolator for a_stop(x)
  gsl_spline *x_stop_interp;  // Interpolator for x_stop(a)
  double aMin_interp_breakpt; // Breakpoint in aMin_interp
  mutable gsl_interp_accel *U2_acc_x, *U2_acc_a;
  mutable gsl_interp_accel *dU2dx_acc_x, *dU2dx_acc_a;
  mutable gsl_interp_accel *X_acc_u, *X_acc_a;
  mutable gsl_interp_accel *Umin_acc, *Umax_acc;
  mutable gsl_interp_accel *aMin_acc[2], *aMax_acc;
  mutable gsl_interp_accel *a_stop_acc, *x_stop_acc;

private:
  gsl_odeiv2_driver *ode_drv;
  gsl_odeiv2_system sys;
  mutable struct duda_params par;
};


////////////////////////////////////////////////////////////////////////
// Specialised by potential and expansion
////////////////////////////////////////////////////////////////////////

// Constant area, point potential
class pwind_hot_pa : public pwind_hot
{
public:
  pwind_hot_pa(const double Gamma_, const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double epsabs_ = 1.0e-6,
	       const double epsrel_ = 1.0e-6,
	       const double interpabs_ = 1.0e-2,
	       const double interprel_ = 1.0e-2,
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       const double amax_grid_ = 1.0e10,
	       const double amax_init = 1000.0);
  virtual ~pwind_hot_pa() { }

  // Limits on the range of the wind
  virtual std::vector<double> xlimits(const double a) const;
};

// Intermediate expansion, point potential
class pwind_hot_pi : public pwind_hot
{
public:
  pwind_hot_pi(const double Gamma_, const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double epsabs_ = 1.0e-6,
	       const double epsrel_ = 1.0e-6,
	       const double interpabs_ = 1.0e-2,
	       const double interprel_ = 1.0e-2,
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       const double amax_grid_ = 1.0e10,
	       const double amax_init = 1000.0);
  virtual ~pwind_hot_pi() { }

  // Limits on the range of the wind
  virtual std::vector<double> xlimits(const double a) const;
};


// Constant solid angle, point potential
class pwind_hot_ps : public pwind_hot
{
public:
  pwind_hot_ps(const double Gamma_, const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double epsabs_ = 1.0e-6,
	       const double epsrel_ = 1.0e-6,
	       const double interpabs_ = 1.0e-2,
	       const double interprel_ = 1.0e-2,
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       const double amax_grid_ = 1.0e10,
	       const double amax_init = 1000.0);
  virtual ~pwind_hot_ps() { }

  // Limits on the range of the wind
  virtual std::vector<double> xlimits(const double a) const;
};


// Constant area, isothermal potential
class pwind_hot_ia : public pwind_hot
{
public:
  pwind_hot_ia(const double Gamma_, const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double epsabs_ = 1.0e-6,
	       const double epsrel_ = 1.0e-6,
	       const double interpabs_ = 1.0e-2,
	       const double interprel_ = 1.0e-2,
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       const double amax_grid_ = 1.0e10,
	       const double amax_init = 1000.0);
  virtual ~pwind_hot_ia() { }

  // Limits on the range of the wind
  virtual std::vector<double> xlimits(const double a) const;
  virtual std::vector<double>
  alimits(const double u,
  	  const double varpi = 0.0,
	  const double varpi_t = 0.0) const; // Overrides default
};

// Intermediate expansion, isothermal potential
class pwind_hot_ii : public pwind_hot
{
public:
  pwind_hot_ii(const double Gamma_, const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double epsabs_ = 1.0e-6,
	       const double epsrel_ = 1.0e-6,
	       const double interpabs_ = 1.0e-2,
	       const double interprel_ = 1.0e-2,
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       const double amax_grid_ = 1.0e10,
	       const double amax_init = 1000.0);
  virtual ~pwind_hot_ii() { }

  // Limits on the range of the wind
  virtual std::vector<double> xlimits(const double a) const;
};


// Constant solid angle, isothermal potential
class pwind_hot_is : public pwind_hot
{
public:
  pwind_hot_is(const double Gamma_, const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double epsabs_ = 1.0e-6,
	       const double epsrel_ = 1.0e-6,
	       const double interpabs_ = 1.0e-2,
	       const double interprel_ = 1.0e-2,
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       const double amax_grid_ = 1.0e10,
	       const double amax_init = 1000.0);
  virtual ~pwind_hot_is() { }

  // Limits on the range of the wind
  virtual std::vector<double> xlimits(const double a) const;
};


#endif
// _PWIND_HOT_H_
