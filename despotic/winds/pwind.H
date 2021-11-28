// This defines a generic momentum-driven wind class

#ifndef _PWIND_H_
#define _PWIND_H_

#include <cmath>
#include <limits>
#include <vector>
#include "pwind_expansion.H"
#include "pwind_geom.H"
#include "pwind_potential.H"
#include "pwind_util.H"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

class pwind {

public:

  // Contructor
  pwind(const double Gamma_,
	const double mach_,
	const pwind_potential *potential_,
	const pwind_expansion *expansion_,
	const pwind_geom *geom_ =
	static_cast<const pwind_geom *>(&pw_geom_sphere),
	const double fcrit_ = 1.0,
	const double jsp_ = 0.0);

  // Destructor; virtual so that we force the derived class's delete
  // method to be called when we are deleted; the only thing this does
  // is restore the old GSL error handler
  virtual ~pwind() { gsl_set_error_handler(err_handler); };

  // Shortcuts to expansion and potential functions, for simplicity
  double y(const double a) const {
    return expansion->y(a);
  }
  std::vector<double> y(const std::vector<double> a) const {
    return expansion->y(a);
  }
  double dyda(const double a) const {
    return expansion->dyda(a);
  }
  std::vector<double> dyda(const std::vector<double> a) const {
    return expansion->dyda(a);
  }
  bool const_sa() const {
    return expansion->const_sa();
  }
  double m(const double a) const {
    return potential->m(a);
  }
  std::vector<double> m(const std::vector<double> a) const {
    return potential->m(a);
  }
  double dmda(const double a) const {
    return potential->dmda(a);
  }
  std::vector<double> dmda(const std::vector<double> a) const {
    return potential->dmda(a);
  }

  // Shortcuts to the geometry functions. These return the positions
  // and radii where lines of sight enter and exit the wind
  std::vector<double> s_cone(const double varpi,
			     const double u=0.0) const;
  std::vector<double> a_cone(const double varpi,
			     const double u=0.0) const;
  
  // Conversions between density, velocity, and position
  // X = density from u and a
  // U2 = square of velocity as a function of x and radius
  // dU2dx = d(U2)/dx
  // dU2da = d(U2)/da
  virtual double X(const double ur, const double a) const = 0;
  virtual double U2(const double x, const double a) const = 0;
  virtual double dU2dx(const double x, const double a) const = 0;
  virtual double dU2da(const double x, const double a) const = 0;

  // Information return
  double getGamma() const { return Gamma; }
  double getMach() const { return mach; }
  double getFcrit() const { return fcrit; }
  double getJsp() const { return jsp; }

  // Parameters that can be changed without affecting the internal state
  void setMach(const double mach_) { mach = mach_; }
  void setGeometry(const pwind_geom* geom_) { geom = geom_; }
  void setFcrit(const double fcrit_) {
    fcrit = fcrit_;
    zeta_M = zetaM(xcrit+log(fcrit_), sx);
    zeta_A = zetaA(xcrit+log(fcrit_), sx);
  }
  void setJsp(const double jsp_) { jsp = jsp_; }

  // Utility quantities
  // alimits = return min and max of a at a given velocity u
  // xlimits = return min and max x that exist a given radius
  // amax = return maximum a reached for a given x
  // a_crit = radius at which a line from infinity enters the wind and
  //          encounters material of velocity u
  // s_crit = distance along LOS at which a line from infinity enters the wind
  virtual std::vector<double> alimits(const double u,
				      const double varpi = 0.0,
				      const double varpi_t = 0.0,
				      const double epsabs = 1.0e-6,
				      const double epsrel = 1.0e-6) const = 0;
  virtual std::vector<double> xlimits(const double a) const = 0;
  virtual double amax(const double x,
		      const double epsabs = 1.0e-6,
		      const double epsrel = 1.0e-6) const = 0;
  std::vector<double> a_crit(const double varpi,
			     const double varpi_t,
			     const double u = 0.0) const {
    return geom->a_crit(varpi, varpi_t, u);
  }
  std::vector<double> s_crit(const double varpi,
			     const double varpi_t,
			     const double u = 0.0) const {
    return geom->s_crit(varpi, varpi_t, u);
  }

  // Trivial convenience functions
  double U(const double x, const double a) const {
    return sqrt(U2(x, a));
  }
  double rho(const double x, const double a) const {
    return pM(x, sx) / (pA(x, sx) * U(x, a) * y(a));
  }
  double fc(const double a, const double fw = -1.0) const {
    if (fw < 0) return zeta_A*y(a)/SQR(a);
    else return fw*y(a)/SQR(a);
  }
  double dfcda(const double a, const double fw = -1.0) const {
    if (const_sa()) return 0.0;
    else if (fw < 0) return zeta_A * (a*dyda(a) - 2.0*y(a)) / CUBE(a);
    else return fw * (a*dyda(a) - 2.0*y(a)) / CUBE(a);
  }
  double xcr() const { return xcrit; }
  double s() const { return sx; }
  double zM() const { return zeta_M; }
  double zA() const { return zeta_A; }
  double uMax() const { return umax; }
  double aMaxAbs() const { return amax_abs; }
  double getJproj() const { return jsp * geom->get_mu(); }
  double uRotMax(const double varpi_t) const { return getJproj()/varpi_t; }
  
  // Functions that perform numerical inversions of the wind
  // acceleration law; not all derived classes will use these, because
  // in many cases the inversion can be done analytically;
  // a_from_u_x -- returns the value of a at which material with a
  // given x has velocity u; returns -1 if the input limits do not
  // properly bracket the root; the argument side indicates which side
  // of the wind to look at: side = 1 means far side, side = -1 means
  // near side, though the distinction matters only if jsp != 0
  // x_from_u_a -- same as a_from_u_x, but gets x given u and a
  double a_from_u_x(const double u,
		    const double x,
		    const double varpi = 0.0,
		    const double varpi_t = 0.0,
		    const double alo = 1.0,
		    const double ahi = 1.0e10,
		    const int maxiter = 100,
		    const double epsabs = 1.0e-6,
		    const double epsrel = 1.0e-6) const;
  double x_from_u_a(const double u,
		    const double a,
		    const double varpi = 0.0,
		    const double varpi_t = 0.0,
		    const double xlo = -100.0,
		    const double xhi = 100.0,
		    const int maxiter = 100,
		    const double epsabs = 1.0e-6,
		    const double epsrel = 1.0e-6) const;

  // Functions to return the differential density and mean density of
  // the wind as a function of radius
  // drhodx = pM / a^2 U zeta_M = differential contribution of mean wind
  //                       density by material with surface density x
  // rho = \int drhodx dx = mean density of wind material at
  //                        radius a
  double drhodx(const double x, const double a) const {
    return pM(x, sx) / (SQR(a) * U(x, a) * zeta_M);
  }
  double rho(const double a,
	     const double epsabs,
	     const double epsrel) const;
  std::vector<double> rho(const std::vector<double>& a,
			  const double epsabs,
			  const double epsrel) const;

  // Function to return the total momentum flux of the wind across the
  // shell at radius a; the quantity returned is dimensionless, and
  // can be converted to the dimensional momentum flux by multiplying
  // by Mdot v0
  double pdot(const double a,
	      const double epsabs,
	      const double epsrel) const;
  std::vector<double> pdot(const std::vector<double>& a,
			   const double epsabs,
			   const double epsrel) const;

  // Function to return the total momentum flux of the wind across the
  // shell at radius a, normalised to the driving momentum flux; the
  // first form uses the approximation Mdot ~= fg M_0/(zeta_M tc),
  // while the second form is the completely general one
  double pdotRel(const double a,
		 const double epsabs,
		 const double epsrel) const;
  double pdotRel(const double a,
		 const double fg,
		 const double tctw,
		 const double epsabs,
		 const double epsrel) const;

  // Function to return the total energy flux of the wind across the
  // shell at radius a; the quantity returned is dimensionless, and
  // can be converted to the dimensional energy flux by multiplying by
  // Mdot v0^2
  double Edot(const double a,
	      const double epsabs,
	      const double epsrel) const;
  std::vector<double> Edot(const std::vector<double>& a,
			   const double epsabs,
			   const double epsrel) const;

  // Functions to return absorption quantities; note that Phi_c is
  // declared as virtual because it is not guaranteed to converge for
  // all wind acceleration laws; declaring it as virtual lets us
  // override its definition for laws where it diverges, and avoid
  // making the integrator choke
  double Phi_uc(const double u,
		const double varpi = 0.0,
		const double varpi_t = 0.0,
		const double a0 = 1.0,
		const double a1 = std::numeric_limits<double>::max(),
		const double epsabs = 1.0e-6,
		const double epsrel = 1.0e-6) const;
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
  double tau_uc(const double u,
		const double tXtw,
		const double fj = 1.0,
		const double boltzfac = 0.0,
		const double varpi = 0.0,
		const double varpi_t = 0.0,
		const double a0 = 1.0,
		const double a1 = std::numeric_limits<double>::max(),
		const double epsabs = 1.0e-6,
		const double epsrel = 1.0e-6) const;
  double tau_uc(const double u,
		const std::vector<double>& u_trans,
		const std::vector<double>& tXtw,
		const double fj = 1.0,
		const double boltzfac = 0.0,
		const double varpi = 0.0,
		const double varpi_t = 0.0,
		const double a0 = 1.0,
		const double a1 = std::numeric_limits<double>::max(),
		const double epsabs = 1.0e-6,
		const double epsrel = 1.0e-6) const;
  double tau_c(const double u,
	       const double tXtw,
	       const double fj = 1.0,
	       const double boltzfac = 0.0,	       
	       const double fw = -1.0,
	       const double varpi = 0.0,
	       const double varpi_t = 0.0,
	       const double a0 = 1.0,
	       const double a1 = std::numeric_limits<double>::max(),
	       const double epsabs = 1.0e-6,
	       const double epsrel = 1.0e-6) const;
  double tau_c(const double u,
	       const std::vector<double>& u_trans,
	       const std::vector<double>& tXtw,
	       const double fj = 1.0,
	       const double boltzfac = 0.0,	       
	       const double fw = -1.0,
	       const double varpi = 0.0,
	       const double varpi_t = 0.0,
	       const double a0 = 1.0,
	       const double a1 = std::numeric_limits<double>::max(),
	       const double epsabs = 1.0e-6,
	       const double epsrel = 1.0e-6) const;
  void tau_interp(const double u,
		  const double tXtw,
		  const double fj,
		  const double boltzfac,
		  const double fw,
		  const double varpi,
		  const double varpi_t,
		  const double a0,
		  const double a1,
		  const bool correlated,
		  const bool reverse,
		  const double epsabs,
		  const double epsrel,
		  gsl_spline **spline,
		  gsl_interp_accel **acc) const;
  void tau_interp_var_fc(const double u,
			 const double tXtw,
			 const double fj,
			 const double boltzfac,
			 const double fw,
			 const double varpi,
			 const double varpi_t,
			 const double a0,
			 const double a1,
			 const double epsabs,
			 const double epsrel,
			 gsl_spline **spline,
			 gsl_interp_accel **acc) const;

  // Functions to return thin, subcritical emission quantities
  double Xi(const double u,
	    const double varpi = 0.0,
	    const double varpi_t = 0.0,
	    const double epsabs = 1.0e-6,
	    const double epsrel = 1.0e-6) const;
  virtual double xi(const double varpi = 0.0,
		    const double varpi_t = 0.0,
		    const double epsabs = 1.0e-6,
		    const double epsrel = 1.0e-6) const;

  // Functions to return LTE emission quantities
  double eta(const double u,
	     const double tXtw,
	     const double fj,
	     const double boltzfac,
	     const bool correlated,
	     const double fw = -1.0,
	     const double varpi = 0.0,
	     const double varpi_t = 0.0,
	     const bool thin = false,
	     const double epsabs = 1.0e-6,
	     const double epsrel = 1.0e-6) const;
  double Psi(const double tXtw,
	     const double fj,
	     const double boltzfac,
	     const bool correlated,
	     const double fw = -1.0,
	     const double varpi = 0.0,
	     const double varpi_t = 0.0,
	     const bool thin = false,
	     const double epsabs = 1.0e-6,
	     const double epsrel = 1.0e-6) const;

  // Vectorized versions of the functions defined above; these are
  // identical except that they take vector<double> arguments for u,
  // and return vector<double> values; these routines use openMP
  // threading to parallelise the calculation
  std::vector<double> Phi_uc(const std::vector<double> &u,
			     const double varpi = 0.0,
			     const double varpi_t = 0.0,
			     const double a0 = 1.0,
			     const double a1 =
			     std::numeric_limits<double>::max(),
			     const double epsabs = 1.0e-6,
			     const double epsrel = 1.0e-6) const;
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
  std::vector<double> tau_uc(const std::vector<double>& u,
			     const double tXtw,
			     const double fj = 1.0,
			     const double boltzfac = 0.0,
			     const double varpi = 0.0,
			     const double varpi_t = 0.0,
			     const double a0 = 1.0,
			     const double a1 =
			     std::numeric_limits<double>::max(),
			     const double epsabs = 1.0e-6,
			     const double epsrel = 1.0e-6) const;
  std::vector<double> tau_uc(const std::vector<double>& u,
			     const std::vector<double>& u_trans,
			     const std::vector<double>& tXtw,
			     const double fj = 1.0,
			     const double boltzfac = 0.0,
			     const double varpi = 0.0,
			     const double varpi_t = 0.0,
			     const double a0 = 1.0,
			     const double a1 =
			     std::numeric_limits<double>::max(),
			     const double epsabs = 1.0e-6,
			     const double epsrel = 1.0e-6) const;
  std::vector<double> tau_c(const std::vector<double>& u,
			    const double tXtw,
			    const double fj = 1.0,
			    const double boltzfac = 0.0,	       
			    const double fw = -1.0,
			    const double varpi = 0.0,
			    const double varpi_t = 0.0,
			    const double a0 = 1.0,
			    const double a1 =
			    std::numeric_limits<double>::max(),
			    const double epsabs = 1.0e-6,
			    const double epsrel = 1.0e-6) const;
  std::vector<double> tau_c(const std::vector<double>& u,
			    const std::vector<double>& u_trans,
			    const std::vector<double>& tXtw,
			    const double fj = 1.0,
			    const double boltzfac = 0.0,	       
			    const double fw = -1.0,
			    const double varpi = 0.0,
			    const double varpi_t = 0.0,
			    const double a0 = 1.0,
			    const double a1 =
			    std::numeric_limits<double>::max(),
			    const double epsabs = 1.0e-6,
			    const double epsrel = 1.0e-6) const;
  std::vector<double> Xi(const std::vector<double>& u,
			 const double varpi = 0.0,
			 const double varpi_t = 0.0,
			 const double epsabs = 1.0e-6,
			 const double epsrel = 1.0e-6) const;
  std::vector<double> eta(const std::vector<double>& u,
			  const double tXtw,
			  const double fj,
			  const double boltzfac,
			  const bool correlated,
			  const double fw = -1.0,
			  const double varpi = 0.0,
			  const double varpi_t = 0.0,
			  const bool thin = false,
			  const double epsabs = 1.0e-6,
			  const double epsrel = 1.0e-6) const;

  // Error management functions
  int getErr() const { return gsl_err_stat; }
  void setErr(int err_stat) const { gsl_err_stat = err_stat; }
  void clearErr() const { gsl_err_stat = GSL_SUCCESS; }
  const char *getErrStr() const { return gsl_strerror(gsl_err_stat); }
  
protected:

  // Expansion law, potential, geometry
  const pwind_potential *potential;
  const pwind_expansion *expansion;
  const pwind_geom *geom;

  // Base data
  const double Gamma;
  double mach, theta, phi, fcrit, jsp;

  // GSL error status and handling
  mutable int gsl_err_stat;
  gsl_error_handler_t *err_handler;
  void checkForErr(int gsl_err) const {
    if (gsl_err != GSL_SUCCESS) gsl_err_stat = gsl_err; }
  
  // Derived data; some must be defined by derived classes
  double sx, xcrit, zeta_M, zeta_A, umax, amax_abs;
};

#endif
// _PWIND_H_
