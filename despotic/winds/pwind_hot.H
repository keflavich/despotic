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
// Data structures to hold pre-tabulated solutions
////////////////////////////////////////////////////////////////////////
typedef struct {
  std::vector<double> umax, q_umax, q_stop, q;
  std::vector<double> umax_lo, q_umax_lo, q_stop_lo, q_lo;
  std::vector<double> umin, loggex;
} hot_wind_table_entry;

typedef struct {
  int nu, nq, ngex, ngex_lo;
  double qmin, qmax, loggex_max, dq, du, dlg, dgex_lo;
  std::vector<double> uh;
  std::vector<hot_wind_table_entry> uh_data;
} hot_wind_table;


////////////////////////////////////////////////////////////////////////
// The generic class
////////////////////////////////////////////////////////////////////////  
class pwind_hot : public pwind {

public:
  
  // The constructor
  pwind_hot(const double Gamma_,
	    const double mach_,
	    const double uh_,
	    const pwind_potential *potential_,
	    const pwind_expansion *expansion_,
	    const pwind_geom *geom_ =
	    static_cast<const pwind_geom *>(&pw_geom_sphere),
	    const double fcrit_ = 1.0,
	    const double jsp_ = 0.0,
	    hot_wind_table *full_tab_ = nullptr);

  // The destructor
  virtual ~pwind_hot();

  // Methods to read and delete a pre-computed data table; note that
  // the method to free a table is trivial, but provided so that it
  // can be accessed by the python code
  static hot_wind_table *read_table(const std::string &dirname,
				    const int yidx, const int midx);
  static void free_table(hot_wind_table *t) { delete t; }

  // Informational routines
  double getUh() const { return uh; }

  // Wind kinematics; these are attached to the generic parent class
  // because all these evaluations are interpolations off a table
  virtual double X(const double ur, const double a) const;
  virtual double U2(const double x, const double a) const;
  virtual double dU2dx(const double x, const double a) const;
  virtual double dU2da(const double x, const double a) const;

  // Limits on the range of the wind
  virtual std::vector<double> alimits(const double u,
				      const double varpi = 0.0,
				      const double varpi_t = 0.0,
				      const double epsabs = 1.0e-6,
				      const double epsrel = 1.0e-6) const;
  virtual std::vector<double> xlimits(const double a) const;
  virtual double amax(const double x,
		      const double epsabs = 1.0e-6,
		      const double epsrel = 1.0e-6) const;

  // Utility functions
  bool hasUmin() const { return expansion->yidx() > potential->midx(); }
  bool isFountain() const { return expansion->yidx() < potential->midx(); }
  double interp(const double x,
		const std::vector<double>& xtab,
		const std::vector<double>& ytab =
		std::vector<double>(),
		const int off = 0) const;

protected:
  const double uh;           // Hot gas velocity
  bool manage_table;         // Flag for if we manage the table ourselves
  hot_wind_table *full_tab;  // Tabulated solution data at all uh values
  hot_wind_table_entry tab;  // Table data interpolated to target uh
};


////////////////////////////////////////////////////////////////////////
// Specialised by potential and expansion
////////////////////////////////////////////////////////////////////////

// Constant area, point potential
class pwind_hot_pa : public pwind_hot
{
public:
  pwind_hot_pa(const double Gamma_,
	       const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       hot_wind_table *full_tab_ = nullptr) :
    pwind_hot(Gamma_, mach_, uh_,
	      static_cast<const pwind_potential *>
	      (&pwind_potential_point),
	      static_cast<const pwind_expansion *>
	      (&pwind_expansion_area),
	      geom_, fcrit_, jsp_, full_tab_)
  { }
  virtual ~pwind_hot_pa() { }
};

// Intermediate expansion, point potential
class pwind_hot_pi : public pwind_hot
{
public:
  pwind_hot_pi(const double Gamma_,
	       const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       hot_wind_table *full_tab_ = nullptr) :
    pwind_hot(Gamma_, mach_, uh_,
	      static_cast<const pwind_potential *>
	      (&pwind_potential_point),
	      static_cast<const pwind_expansion *>
	      (&pwind_expansion_intermediate),
	      geom_, fcrit_, jsp_, full_tab_)
  { }
  virtual ~pwind_hot_pi() { }
};


// Constant solid angle, point potential
class pwind_hot_ps : public pwind_hot
{
public:
  pwind_hot_ps(const double Gamma_,
	       const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       hot_wind_table *full_tab_ = nullptr) :
    pwind_hot(Gamma_, mach_, uh_,
	      static_cast<const pwind_potential *>
	      (&pwind_potential_point),
	      static_cast<const pwind_expansion *>
	      (&pwind_expansion_solid_angle),
	      geom_, fcrit_, jsp_, full_tab_)
  { }
  virtual ~pwind_hot_ps() { }
};


// Constant area, isothermal potential
class pwind_hot_ia : public pwind_hot
{
public:
  pwind_hot_ia(const double Gamma_,
	       const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       hot_wind_table *full_tab_ = nullptr) :
    pwind_hot(Gamma_, mach_, uh_,
	      static_cast<const pwind_potential *>
	      (&pwind_potential_isothermal),
	      static_cast<const pwind_expansion *>
	      (&pwind_expansion_area),
	      geom_, fcrit_, jsp_, full_tab_)
  { }
  virtual ~pwind_hot_ia() { }
};

// Intermediate expansion, isothermal potential
class pwind_hot_ii : public pwind_hot
{
public:
  pwind_hot_ii(const double Gamma_,
	       const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       hot_wind_table *full_tab_ = nullptr) :
    pwind_hot(Gamma_, mach_, uh_,
	      static_cast<const pwind_potential *>
	      (&pwind_potential_isothermal),
	      static_cast<const pwind_expansion *>
	      (&pwind_expansion_intermediate),
	      geom_, fcrit_, jsp_, full_tab_)
  { }
  virtual ~pwind_hot_ii() { }
};


// Constant solid angle, isothermal potential
class pwind_hot_is : public pwind_hot
{
public:
  pwind_hot_is(const double Gamma_,
	       const double mach_,
	       const double uh_,
	       const pwind_geom *geom_ =
	       static_cast<const pwind_geom *>(&pw_geom_sphere),
	       const double fcrit_ = 1.0,
	       const double jsp_ = 0.0,
	       hot_wind_table *full_tab_ = nullptr) :
    pwind_hot(Gamma_, mach_, uh_,
	      static_cast<const pwind_potential *>
	      (&pwind_potential_isothermal),
	      static_cast<const pwind_expansion *>
	      (&pwind_expansion_solid_angle),
	      geom_, fcrit_, jsp_, full_tab_)
  { }
  virtual ~pwind_hot_is() { }
};


#endif
// _PWIND_HOT_H_
