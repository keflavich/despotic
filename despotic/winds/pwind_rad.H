// A class representing radiatively-driven winds

#ifndef _PWIND_RAD_H_
#define _PWIND_RAD_H_

#include "pwind.H"
#include "pwind_expansion.H"
#include "pwind_geom.H"
#include "pwind_potential.H"
#include <cmath>
#include <limits>
#include <vector>

// Radiation-driven winds are somewhat tricky, because depending on
// the potential and expansion factor there are several possible
// behaviors for the trajectories of material with x = xcrit and x =
// -infinity, which bound the velocities occupied by the winds. This
// makes determining limits on a and x for a given u, or on x for a
// given x, complicated. These different classes need to be handled
// individually, since they affect things like the number of roots
// that the equation U2(x, a) = u2 has. For reference, the possible
// classes of behavior are:
//
// Point potential, constant area:
// Velocity u of x = -infinity material is monotonically increasing
// with a, but reaches a finite maximum u_max as a->infinity, while
// velocity of x = xcrit material is always zero. As x->xcrit from
// below, velocity approaches zero, but material is still able to
// travel outward to arbitrary distances a
//
// Point potential, intermediate area:
// x = -infinity material is the same as in the constant area case,
// but material with x = xcrit also rises monotonically to asymptote
// to a finite velocity as a->infinity
//
// Point potential, constant solid angle:
// Same qualitative behavior as the point potential, intermediate area
// case
//
// Isothermal potential, constant area:
// Material with x = -infinity, and at all x < xcrit, has a velocity
// that increases at small a, reaches a finite maximum at finite a,
// and then decreases with a until reaching zero at some finite value
// of a. The radius at which the material reaches zero velocity is a
// decreasing function of x, and reaches 0 at x = xcrit.
//
// Isothermal potential, intermediate area:
// Qualitatively the same as the isothermal potential, constant area
// case
//
// Isothermal potential, constant solid angle:
// Material at all values of x increases in velocity at small a,
// reaches a finite maximum velocity at finite a, then has a velocity
// that declines with a past the maximum. The velocity goes to zero at
// finite a, and the radius a at which this happens is an decreasing
// function of x, but even for x = xcrit the radius a at which the
// velocity reaches 0 is > 1.

////////////////////////////////////////////////////////////////////////
// The generic class
////////////////////////////////////////////////////////////////////////
class pwind_rad : public pwind {

public:

  // The constructor
  pwind_rad(const double Gamma_,
	    const double mach_,
	    const double tau0_,
	    const pwind_potential *potential_,
	    const pwind_expansion *expansion_,
	    const pwind_geom *geom_ =
	    static_cast<const pwind_geom *>(&pw_geom_sphere),
	    const double fcrit_ = 1.0,
	    const double jsp_ = 0.0);

  // The destructor
  virtual ~pwind_rad() { }

  // Method to get xcrit
  virtual double xCrit(const int maxiter = 100,
		       const double epsabs = 1.0e-6,
		       const double epsrel = 1.0e-6) const;

  // Method to get X; since we have to do this numerically, it is the
  // same for all geometries and potentials
  virtual double X(const double ur, const double a) const;

  // The wind acceleration law RHS for a radiatively-driven wind
  virtual double dU2da(const double x, const double a) const;

  // A method that returns the velocity in the limit x -> -infinity;
  // this is generally useful, but must be specified separately for
  // each potential and expansion factor, so we define it as a pure
  // virtual here
  virtual double U2max(const double a) const = 0;

  // Utility functions; these just evaluate certain expressions that
  // require care to avoid numerical problems in some regions of
  // parameter space

  // e^-x [1 - exp(-tau e^x)]
  double fac1(const double x) const;

  // e^-x [Ei(-tau e^x/a) - Ei(-tau) + ln a]
  double fac2(const double x, const double a) const;

  // e^-x a [1 - exp(-tau e^x/a^2)]
  double fac3(const double x, const double a) const;

  // e^-x sqrt(pi tau0 e^x) *
  //    [ erf(sqrt(tau0 e^x)) - erf(sqrt(tau0 e^x/a^2)) ]
  double fac4(const double x, const double a) const;
  
  // e^-x [ exp(-tau e^x/a) - exp(-tau e^x) ]
  double fac5(const double x, const double a) const;

protected:

  double tau0; // Optical depth
};

////////////////////////////////////////////////////////////////////////
// Specialised by potential
////////////////////////////////////////////////////////////////////////

// Point potential; radiatively-driven winds with this potential
// always escape to infinity, reaching asymptotically constant
// velocities as they do so
class pwind_rad_point : public pwind_rad {

public:
  pwind_rad_point(const double Gamma_,
		  const double mach_,
		  const double tau0_,
		  const pwind_expansion *expansion_,
		  const pwind_geom *geom_ =
		  static_cast<const pwind_geom *>(&pw_geom_sphere),
		  const double fcrit_ = 1.0,
		  const double jsp_ = 0.0);
  virtual ~pwind_rad_point() { }

  // Velocity of fluid element with x -> -infinity; same for all point
  // potentials, regardless of expansion factor
  virtual double U2max(const double a) const;

  // Radius at which fluid element with x->-infinity reaches a
  // specified line of sight velocity; same for all point potentials
  double a_from_u_max(const double u,
		      const double varpi = 0.0,
		      const double varpi_t = 0.0) const;

  // Limits on x and maximum value of a at a given x; same for all
  // point potentials
  virtual std::vector<double> xlimits(const double a) const;
  virtual double amax(const double x,
		      const double epsabs = 1.0e-6,
		      const double epsrel = 1.0e-6) const;
};

// Isothermal potential; radiatively-driven winds with this potential
// are always fountains, and always all back; they reach a maximum
// velocity at a finite value of a that depends on x
class pwind_rad_isothermal : public pwind_rad {
public:
  pwind_rad_isothermal(const double Gamma_,
		       const double mach_,
		       const double tau0_,
		       const pwind_expansion *expansion_,
		       const pwind_geom *geom_ =
		       static_cast<const pwind_geom *>(&pw_geom_sphere),
		       const double fcrit_ = 1.0,
		       const double jsp_ = 0.0);
  virtual ~pwind_rad_isothermal() { }

  // Velocity of fluid element with x -> -infinity; same for all 
  // isothermal potentials, regardless of expansion factor
  virtual double U2max(const double a) const;

  // Derivative of U2max wrt a
  virtual double dU2maxda(const double a) const;

  // Function to return the radius at which maximum value of u = given
  // value; analogous to a_from_u_x, but using U2max in place of U2
  double a_from_u_max(const double u,
		      const double varpi = 0.0,
		      const double varpi_t = 0.0,
		      const double alo = 1.0,
		      const double ahi = 0.0,
		      const int maxiter = 100,
		      const double epsabs = 1.0e-6,
		      const double epsrel = 1.0e-6) const;

  // This function returns the radius a at which the line of sight
  // velocity reaches a maximum for material with surface density x,
  // viewed at impact parameter varpi
  double a_max_u_from_varpi(const double x,
			    const double varpi = 0.0,
			    const double varpi_t = 0.0,
			    const int maxiter = 100,
			    const double epsabs = 1.0e-6,
			    const double epsrel = 1.0e-6) const;
  
  // Same as a_max_u_from_varpi, but using u = umax, i.e., in the
  // limit x -> -infinity
  double a_max_umax_from_varpi(const double varpi = 0.0,
			       const double varpi_t = 0.0,
			       const int maxiter = 100,
			       const double epsabs = 1.0e-6,
			       const double epsrel = 1.0e-6) const;


protected:
  double a_maxu;  // Radius at which maximum velocity is reached for
		  // material with x -> -infinity
};

////////////////////////////////////////////////////////////////////////
// Further specialisation by expansion
////////////////////////////////////////////////////////////////////////

// Constant area, point potential
class pwind_rad_pa : public pwind_rad_point
{
public:
  pwind_rad_pa(const double Gamma_,
	       const double mach_,
	       const double tau0_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0) :
    pwind_rad_point(Gamma_, mach_, tau0_,
		    static_cast<const pwind_expansion *>
		    (&pwind_expansion_area), geom_,
		    fcrit_, jsp_) { }
  virtual ~pwind_rad_pa() { }

  // Wind kinematics
  virtual double U2(const double x, const double a) const;
  virtual double dU2dx(const double x, const double a) const;
  virtual std::vector<double> alimits(const double u,
				      const double varpi = 0.0,
				      const double varpi_t = 0.0,
				      const double epsabs = 1.0e-6,
				      const double epsrel = 1.0e-6) const;
  
  // Override Phi_c, since it diverges if we try to integrate to infinity
  virtual double Phi_c(const double u,
		       const double fw = -1.0,
		       const double varpi = 0.0,
		       const double varpi_t = 0.0,
		       const double a0 = 1.0,
		       const double a1 = std::numeric_limits<double>::max(),
		       const double epsabs = 1.0e-6,
		       const double epsrel = 1.0e-6,
		       const std::vector<double>& alimits_ =
		       std::vector<double>()) const;
  virtual std::vector<double> Phi_c(const std::vector<double>& u,
				    const double fw = -1.0,
				    const double varpi = 0.0,
				    const double varpi_t = 0.0,
				    const double a0 = 1.0,
				    const double a1 =
				    std::numeric_limits<double>::max(),
				    const double epsabs = 1.0e-6,
				    const double epsrel = 1.0e-6,
				    const std::vector<double>& alimits_ =
				    std::vector<double>()) const;
};

// Intermediate expansion, point potential
class pwind_rad_pi : public pwind_rad_point
{
public:
  pwind_rad_pi(const double Gamma_,
	       const double mach_,
	       const double tau0_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0);
  virtual ~pwind_rad_pi() { }

  // Wind kinematics
  virtual double U2(const double x, const double a) const;
  virtual double dU2dx(const double x, const double a) const;
  virtual std::vector<double> alimits(const double u,
				      const double varpi = 0.0,
				      const double varpi_t = 0.0,
				      const double epsabs = 1.0e-6,
				      const double epsrel = 1.0e-6) const;

protected:
  double u_min_infinity; // Minimum velocity as a->inifinity
};

// Constant solid angle, point potential
class pwind_rad_ps : public pwind_rad_point
{
public:
  pwind_rad_ps(const double Gamma_,
	       const double mach_,
	       const double tau0_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0);
  virtual ~pwind_rad_ps() { }

  // Wind kinematics
  virtual double U2(const double x, const double a) const;
  virtual double dU2dx(const double x, const double a) const;
  virtual std::vector<double> alimits(const double u,
				      const double varpi = 0.0,
				      const double varpi_t = 0.0,
				      const double epsabs = 1.0e-6,
				      const double epsrel = 1.0e-6) const;

protected:
  double u_min_infinity; // Minimum velocity as a->inifinity
};


// Constant area, isothermal potential
class pwind_rad_ia : public pwind_rad_isothermal
{
public:
  pwind_rad_ia(const double Gamma_,
	       const double mach_,
	       const double tau0_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0) :
    pwind_rad_isothermal(Gamma_, mach_, tau0_,
			 static_cast<const pwind_expansion *>
			 (&pwind_expansion_area), geom_,
			 fcrit_, jsp_) { }
  virtual ~pwind_rad_ia() { }

  // Wind kinematics
  virtual double U2(const double x, const double a) const;
  virtual double dU2dx(const double x, const double a) const;
  virtual std::vector<double> alimits(const double u,
				      const double varpi = 0.0,
				      const double varpi_t = 0.0,
				      const double epsabs = 1.0e-6,
				      const double epsrel = 1.0e-6) const;
  virtual std::vector<double> xlimits(const double a) const;
  virtual double amax(const double x,
		      const double epsabs = 1.0e-6,
		      const double epsrel = 1.0e-6) const;
};


// Intermediate area, isothermal potential
class pwind_rad_ii : public pwind_rad_isothermal
{
public:
  pwind_rad_ii(const double Gamma_,
	       const double mach_,
	       const double tau0_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0) :
    pwind_rad_isothermal(Gamma_, mach_, tau0_,
			 static_cast<const pwind_expansion *>
			 (&pwind_expansion_intermediate), geom_,
			 fcrit_, jsp_) { }
  virtual ~pwind_rad_ii() { }

  // Wind kinematics
  virtual double U2(const double x, const double a) const;
  virtual double dU2dx(const double x, const double a) const;
  virtual std::vector<double> alimits(const double u,
				      const double varpi = 0.0,
				      const double varpi_t = 0.0,
				      const double epsabs = 1.0e-6,
				      const double epsrel = 1.0e-6) const;
  virtual std::vector<double> xlimits(const double a) const;
  virtual double amax(const double x,
		      const double epsabs = 1.0e-6,
		      const double epsrel = 1.0e-6) const;
};


// Constant solid angle, isothermal potential
class pwind_rad_is : public pwind_rad_isothermal
{
public:
  pwind_rad_is(const double Gamma_, const double mach_,
	       const double tau0_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0);
  virtual ~pwind_rad_is() { }

  // Wind kinematics
  virtual double U2(const double x, const double a) const;
  virtual double dU2dx(const double x, const double a) const;
  virtual std::vector<double> alimits(const double u,
				      const double varpi = 0.0,
				      const double varpi_t = 0.0,
				      const double epsabs = 1.0e-6,
				      const double epsrel = 1.0e-6) const;
  virtual std::vector<double> xlimits(const double a) const;
  virtual double amax(const double x,
		      const double epsabs = 1.0e-6,
		      const double epsrel = 1.0e-6) const;

protected:

  // Radius at which maximum velocity occurs for x = xcrit
  double a_maxu_xcrit;

  // Maximum velocity reached for x = xcrit
  double umax_xcrit;

  // Radius at which u = 0 for x = xcrit
  double amax_xcrit;

  // In addition, the following are defined in parent classes:
  //
  // amax_abs = maximum radius a reached by material with x =
  // -infinity; defined in pwind
  //
  // a_maxu = radius at which material with x = -infinity reaches
  // maximum velocity; defined in pwind_rad
  //
  // umax = maximum velocity reached by material with x = -infinity;
  // defined in pwind
};


#endif
// _PWIND_RAD_H_
