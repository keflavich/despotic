// A class representing ideal momentum-driven winds

#ifndef _PWIND_IDEAL_H_
#define _PWIND_IDEAL_H_

#include "pwind.H"
#include "pwind_expansion.H"
#include "pwind_geom.H"
#include "pwind_potential.H"
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>

////////////////////////////////////////////////////////////////////////
// The generic class
////////////////////////////////////////////////////////////////////////
class pwind_ideal : public pwind {

public:

  // The constructor
  pwind_ideal(const double Gamma_,
	      const double mach_,
	      const pwind_potential *potential_,
	      const pwind_expansion *expansion_,
	      const pwind_geom *geom_ =
	      static_cast<const pwind_geom *>(&pw_geom_sphere),
	      const double fcrit_ = 1.0,
	      const double jsp_ = 0.0);

  // The destructor
  virtual ~pwind_ideal() { }
  
  // The wind acceleration law RHS for an ideal wind
  virtual double dU2da(const double x, const double a) const;
};

////////////////////////////////////////////////////////////////////////
// Specalised by expansion and potential
////////////////////////////////////////////////////////////////////////

// Constant area, point potential
class pwind_ideal_pa : public pwind_ideal
{
public:
  pwind_ideal_pa(const double Gamma_,
		 const double mach_,
		 const pwind_geom *geom_ =
		 static_cast<const pwind_geom *>(&pw_geom_sphere),
		 const double fcrit_ = 1.0,
		 const double jsp_ = 0.0) :
    pwind_ideal(Gamma_, mach_,
		static_cast<const pwind_potential *>(&pwind_potential_point),
		static_cast<const pwind_expansion *>(&pwind_expansion_area),
		geom_, fcrit_, jsp_) { }
  virtual ~pwind_ideal_pa() { }

  // Wind kinematics
  virtual double X(const double ur, const double a) const;
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

// Intermediate area, point potential
class pwind_ideal_pi : public pwind_ideal
{
public:
  pwind_ideal_pi(const double Gamma_,
		 const double mach_,
		 const pwind_geom *geom_ =
		 static_cast<const pwind_geom *>(&pw_geom_sphere),
		 const double fcrit_ = 1.0,
		 const double jsp_ = 0.0) :
    pwind_ideal(Gamma_, mach_,
		static_cast<const pwind_potential *>(&pwind_potential_point),
		static_cast<const pwind_expansion *>
		(&pwind_expansion_intermediate),
		geom_, fcrit_, jsp_) { }
  virtual ~pwind_ideal_pi() { }

  // Wind kinematics
  virtual double X(const double ur, const double a) const;
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

// Constant solid angle, point potential
class pwind_ideal_ps : public pwind_ideal
{
public:
  pwind_ideal_ps(const double Gamma_,
		 const double mach_,
		 const pwind_geom *geom_ =
		 static_cast<const pwind_geom *>(&pw_geom_sphere),
		 const double fcrit_ = 1.0,
		 const double jsp_ = 0.0) :
    pwind_ideal(Gamma_, mach_,
		static_cast<const pwind_potential *>(&pwind_potential_point),
		static_cast<const pwind_expansion *>
		(&pwind_expansion_solid_angle),
		geom_, fcrit_, jsp_) { }
  virtual ~pwind_ideal_ps() { }

  // Wind kinematics
  virtual double X(const double ur, const double a) const;
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


// Constant area, isothermal potential
class pwind_ideal_ia : public pwind_ideal
{
public:
  pwind_ideal_ia(const double Gamma_,
		 const double mach_,
		 const pwind_geom *geom_ =
		 static_cast<const pwind_geom *>(&pw_geom_sphere),
		 const double fcrit_ = 1.0,
		 const double jsp_ = 0.0) :
    pwind_ideal(Gamma_, mach_,
		static_cast<const pwind_potential *>
		(&pwind_potential_isothermal),
		static_cast<const pwind_expansion *>(&pwind_expansion_area),
		geom_, fcrit_, jsp_) { }
  virtual ~pwind_ideal_ia() { }

  // Wind kinematics
  virtual double X(const double ur, const double a) const;
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

// Intermediate area, isothermal potential
class pwind_ideal_ii : public pwind_ideal
{
public:
  pwind_ideal_ii(const double Gamma_,
		 const double mach_,
		 const pwind_geom *geom_ =
		 static_cast<const pwind_geom *>(&pw_geom_sphere),
		 const double fcrit_ = 1.0,
		 const double jsp_ = 0.0) :
    pwind_ideal(Gamma_, mach_,
		static_cast<const pwind_potential *>
		(&pwind_potential_isothermal),
		static_cast<const pwind_expansion *>
		(&pwind_expansion_intermediate),
		geom_, fcrit_, jsp_) { }
  virtual ~pwind_ideal_ii() { }

  // Wind kinematics
  virtual double X(const double ur, const double a) const;
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

// Constant solid angle, isothermal potential
class pwind_ideal_is : public pwind_ideal
{
public:
  pwind_ideal_is(const double Gamma_,
		 const double mach_,
		 const pwind_geom *geom_ =
		 static_cast<const pwind_geom *>(&pw_geom_sphere),
		 const double fcrit_ = 1.0,
		 const double jsp_ = 0.0) :
    pwind_ideal(Gamma_, mach_,
		static_cast<const pwind_potential *>
		(&pwind_potential_isothermal),
		static_cast<const pwind_expansion *>
		(&pwind_expansion_solid_angle),
		geom_, fcrit_, jsp_) { }
  virtual ~pwind_ideal_is() { }

  // Wind kinematics
  virtual double X(const double ur, const double a) const;
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

#endif
// _PWIND_IDEAL_H_
