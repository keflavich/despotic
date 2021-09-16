#include "pwind_rad.H"
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>
#include <limits>
#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Functions and data structures required by the GSL
////////////////////////////////////////////////////////////////////////

struct xcrit_params {
  const pwind_rad *w;
};

static inline double xcrit_func(double x, void *params) {
  struct xcrit_params *par = (struct xcrit_params *) params;
  double ret = par->w->getGamma() * par->w->fac1(x) - 1.0;
  return ret;
}

////////////////////////////////////////////////////////////////////////
// Base pwind_rad class
////////////////////////////////////////////////////////////////////////
pwind_rad::pwind_rad(const double Gamma_,
		     const double mach_,
		     const double tau0_,
		     const pwind_potential *potential_,
		     const pwind_expansion *expansion_,
		     const pwind_geom *geom_,
		     const double fcrit_,
		     const double jsp_) :
  pwind(Gamma_, mach_, potential_, expansion_, geom_,
        fcrit_, jsp_),
  tau0(tau0_)
{
  xcrit = xCrit(100, 1.0e-12, 1.0e-12); // Use very high accuracy,
					// since this only needs to be
					// done once
  zeta_M = zetaM(xcrit + log(fcrit_), sx);
  zeta_A = zetaA(xcrit + log(fcrit_), sx);
}

double pwind_rad::xCrit(const int maxiter,
			const double epsabs,
			const double epsrel) const {
  // Prepare data to pass to the GSL
  struct xcrit_params params;
  params.w = this;
  gsl_function F;
  F.function = &xcrit_func;
  F.params = &params;
  
  // Specify that we want to use the Brent solver, and allocate its
  // workspace
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  // Set brackets
  double xlo = -100.0;
  double xhi = log(Gamma)+0.1;
  gsl_root_fsolver_set(s, &F, xlo, xhi);

  // Loop
  double xcr = 0.0;
  for (int iter=0, status=GSL_CONTINUE;
       (iter<maxiter) && (status==GSL_CONTINUE);
       iter++) {
    status = gsl_root_fsolver_iterate(s);
    xcr = gsl_root_fsolver_root(s);
    xlo = gsl_root_fsolver_x_lower(s);
    xhi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(xlo, xhi, epsabs, epsrel);
  }

  // Free workspace and return result
  gsl_root_fsolver_free(s);
  return xcr;
}
  
double
pwind_rad::X(const double ur, const double a) const {
  return x_from_u_a(ur, a);
}

double
pwind_rad::dU2da(const double x, const double a) const {
  double taufac = tau0*exp(x)/y(a);
  if (taufac < 1.0e-10)
    return (Gamma*tau0*(1.0 - 0.5*taufac) - m(a)) / SQR(a);
  else
    return (y(a)*Gamma*exp(-x)*(1.0-exp(-taufac)) - m(a)) / SQR(a);
}

////////////////////////////////////////////////////////////////////////
// Convenience functions
////////////////////////////////////////////////////////////////////////

double
pwind_rad::fac1(const double x) const {
  /* This returns e^-x (1 - exp(-tau e^x)) */
  double taufac = tau0*exp(x);
  if (taufac < 1.0e-10)
    /* Limit as x -> -infinity */
    return tau0*(1.0 - 0.5*taufac);
  else
    /* Finite x */
    return exp(-x)*(1.0 - exp(-taufac));
}

double
pwind_rad::fac2(const double x, const double a) const {
  /* This returns e^-x (Ei(-tau e^x/a) - Ei(-tau) + ln a)  */
  double taufac1 = tau0*exp(x);
  double taufac2 = tau0*exp(x)/a;
  if (taufac1 < 1.0e-10)
    /* Limit as x -> -infinity */
    return tau0*(1.0-1.0/a+exp(x)*tau0/4.0*(1.0/SQR(a)-1.0));
  else if (taufac2 < 1.0e-10) {
    /* Limit as a -> +infinity */
    double ret;
    if (taufac1 < 100) {
      // Case where taufac1 is not too big, so we don't overflow when
      // trying to compute its exponential integral
      ret = exp(-x) * (M_EULER - gsl_sf_expint_Ei(-taufac1) +
		       log(taufac1));
    } else {
      // Case where taufac1 is very large, so we need to use an
      // expansion for Ei
      ret = exp(-taufac1)*(exp(-2.0*x)/tau0) +
	exp(-x)*(M_EULER + x -log(1.0/tau0));
    }
    return ret;
  } else {
    /* Neither limit */
    double ret;
    if (taufac1 < 100) {
      // Case where taufac1 is not too big, so we don't underflow when
      // trying to compute its exponential integral
      ret = exp(-x)*(gsl_sf_expint_Ei(-taufac2) -
		     gsl_sf_expint_Ei(-taufac1) + log(a));
    } else {
      if (taufac2 < 100) {
	// We can compute Ei(taufac2) but not Ei(taufac1), so series
	// expand Ei(taufac1)
	ret = exp(-x) * (gsl_sf_expint_Ei(-taufac2) +
			 exp(-taufac1)*exp(-x)/tau0 + log(a));
      } else {
	// We can't compute Ei(taufac1) or Ei(taufac2)
	ret = exp(-x) * (exp(-taufac1)*exp(-x)/tau0 -
			 exp(-taufac2)*a*exp(-x)/tau0 +
			 log(a));
      }
    }
    return ret;
  }
}

double
pwind_rad::fac3(const double x, const double a) const {
  /* This returns e^-x a (1 - exp(-tau e^x/a^2)) */
  double taufac = tau0*exp(x)/SQR(a);
  if (taufac < 1.0e-10)
    /* Limit as x -> -infinity or a -> +infinity */
    return tau0/a * (1.0 - 0.5*taufac);
  else
    /* Finite x and a */
    return exp(-x)*a*(1.0 - exp(-taufac));
}

double
pwind_rad::fac4(const double x, const double a) const {
  /* This returns e^-x sqrt(pi tau0 e^x) *
     [ erf(sqrt(tau0 e^x)) - erf(sqrt(tau0 e^x/a^2)) ] */
  double taufac1 = sqrt(tau0*exp(x));
  double taufac2 = taufac1/a;
  if (taufac1 < 1.0e-10)
    /* Limit as x -> -infinity */
    return 2.0*(1.0-1.0/a)*tau0 -
      (2.0/3.0)*exp(x)*SQR(tau0)*(1.0-1.0/(CUBE(a)));
  else if (a - 1.0 > 1.0e-10) {
    /* Finite x, a != 1 */
    return exp(-x) * sqrt(M_PI) * taufac1 *
      (gsl_sf_erf(taufac1) - gsl_sf_erf(taufac2));
  } else {
    /* Finite x, a -> 1 */
    return 2.0*tau0*exp(-taufac1*taufac1) *
      (a - 1.0 - fac1(x)*SQR(a-1.0));
  }
}

double
pwind_rad::fac5(const double x, const double a) const {
  /* This returns e^-x [ exp(-tau e^x/a) - exp(-tau e^x) ] */
  double taufac1 = tau0*exp(x);
  double taufac2 = taufac1 / a;
  if (taufac1 < 1.0e-10) {
    // Limit as x -> -infinity
    return tau0*((1.0 - 1.0/a) * 0.5*(1.0-1.0/SQR(a))*taufac1);
  } else {
    // Finite x
    return exp(-x) * (exp(-taufac2) - exp(-taufac1));
  }
}

////////////////////////////////////////////////////////////////////////
// Point potential funcions
////////////////////////////////////////////////////////////////////////

pwind_rad_point::
pwind_rad_point(const double Gamma_,
		const double mach_,
		const double tau0_,
		const pwind_expansion *expansion_,
		const pwind_geom *geom_,
		const double fcrit_,
		const double jsp_) :
  pwind_rad(Gamma_, mach_, tau0_,
	    static_cast<const pwind_potential *>(&pwind_potential_point),
	    expansion_,
	    geom_, fcrit_, jsp_)
{
  umax = sqrt(Gamma_*tau0_-1.0);
  amax_abs = numeric_limits<double>::max();
}
double
pwind_rad_point::U2max(const double a) const {
  return (1.0-1.0/a)*(Gamma*tau0-1.0);
}
double
pwind_rad_point::a_from_u_max(const double u,
			      const double varpi,
			      const double varpi_t) const {
  double vp2 = SQR(varpi) + SQR(varpi_t);
  if (vp2 == 0.0) {
    // varpi = 0 case can be written out explicitly
   return (Gamma*tau0-1.0) / (Gamma*tau0-1.0-SQR(u));
  } else {
    // For varpi != 0, the problem is a 3rd order polynomial; use GSL
    // to solve; there is guaranteed to be one real root > varpi
    double fac = SQR(u)/(Gamma*tau0-1.0) - 1.0;
    double a[3];
    int nroots =
      gsl_poly_solve_cubic(1.0/fac, vp2/fac, -vp2/fac,
			   a, a+1, a+2);
    for (int i=0; i<nroots; i++) {
      if (SQR(a[i]) > vp2 && a[i] > 0.0) {
	return a[i];
      }
    }
  }
  return -1.0; // Sentinel value; should never be encountered
}
vector<double>
pwind_rad_point::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = -numeric_limits<double>::max();
  ret[1] = xcrit;
  return ret;
}
double
pwind_rad_point::amax(const double x,
		      const double epsabs,
		      const double epsrel) const {
  return numeric_limits<double>::max();
}

////////////////////////////////////////////////////////////////////////
// Isothermal potential funcions
////////////////////////////////////////////////////////////////////////

pwind_rad_isothermal::
pwind_rad_isothermal(const double Gamma_,
		     const double mach_,
		     const double tau0_,
		     const pwind_expansion *expansion_,
		     const pwind_geom *geom_,
		     const double fcrit_,
		     const double jsp_) :
  pwind_rad(Gamma_, mach_, tau0_,
	    static_cast<const pwind_potential *>(&pwind_potential_isothermal),
	    expansion_, geom_, fcrit_, jsp_)
{
  umax = sqrt(Gamma_*tau0_ - log(Gamma_*tau0_) - 1.0);
  a_maxu = Gamma_*tau0_;
  // Note: use high tolerance here, because we need to compute this
  // limit numerically, and it only needs to be done once
  amax_abs = a_from_u_max(0.0, 0.0, 0.0, 1.0+1.0e-12, 0.0, 100,
			  1.0e-12, 1.0e-12);
}

double
pwind_rad_isothermal::U2max(const double a) const {
  return (1.0-1.0/a)*Gamma*tau0 - log(a);
}

double
pwind_rad_isothermal::dU2maxda(const double a) const {
  double ret = 1.0/a * (Gamma*tau0/a - 1.0);
  return ret;
}

struct U2max_params {
  double u, a, vp2;
  const pwind_rad *pw;
};

static double U2max_func_loga(double loga, void *params) {
  // This routine returns the residual for U2max(a) - U2_target
  struct U2max_params *par = (struct U2max_params *) params;
  double a = exp(loga);
  double ret = SQR(par->u) - par->pw->U2max(a) * (1.0 - par->vp2/SQR(a));
  // GSL doesn't like infinities
  if (ret > numeric_limits<double>::max())
    ret = numeric_limits<double>::max();
  return ret;
}

double
pwind_rad_isothermal::a_from_u_max(const double u,
				   const double varpi,
				   const double varpi_t,
				   const double alo,
				   const double ahi,
				   const int maxiter,
				   const double epsabs,
				   const double epsrel) const {


  // Prepare data for the GSL
  struct U2max_params params;
  params.u = u;
  params.vp2 = SQR(varpi) + SQR(varpi_t);
  params.pw = this;
  gsl_function F;
  F.function = &U2max_func_loga;
  F.params = &params;

  // Specify that we want to use the Brent solver, and allocate its
  // workspace
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  // Set brakets
  double loga_lo = SQR(alo) < params.vp2 ? 0.5*log(params.vp2) : log(alo);
  double loga_hi;
  if (ahi == 0.0)
    loga_hi = Gamma*tau0;
  else
    loga_hi = log(ahi);
  gsl_root_fsolver_set(s, &F, loga_lo, loga_hi);

  // Ensure that the root is properly bracketed; return -1 if not
  double sgn1 = U2max_func_loga(loga_lo, &params);
  if (loga_lo ==  0.5*log(params.vp2)) sgn1 = 0.0;
  double sgn2 = U2max_func_loga(loga_hi, &params);
  if (sgn1*sgn2 > 0){
    gsl_root_fsolver_free(s);
    return -1;
  }

  // Loop
  double loga = 0.0;
  for (int iter=0, status=GSL_CONTINUE;
       (iter<maxiter) && (status==GSL_CONTINUE);
       iter++) {
    status = gsl_root_fsolver_iterate(s);
    loga = gsl_root_fsolver_root(s);
    loga_lo = gsl_root_fsolver_x_lower(s);
    loga_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(loga_lo, loga_hi, epsabs, epsrel);
  }

  // Free workspace and return result
  gsl_root_fsolver_free(s);
  return exp(loga);
}

struct a_max_u_params {
  double x, vp2;
  const pwind_rad_isothermal *pw;
};

static double a_max_u_loga(double loga, void *params) {
  // This routine returns the residual of equation solved by
  // a_max_umax_from_varpi
  struct a_max_u_params *par = (struct a_max_u_params *) params;
  double a = exp(loga);
  double ret = (1.0 - par->vp2/SQR(a)) * par->pw->dU2da(par->x, a) +
    par->vp2 / (a*SQR(a)) * par->pw->U2(par->x, a);
  // GSL doesn't like infinities
  if (ret > numeric_limits<double>::max())
    ret = numeric_limits<double>::max();
  else if (ret < -numeric_limits<double>::max())
    ret = -numeric_limits<double>::max();
  return ret;
}

double
pwind_rad_isothermal::
a_max_u_from_varpi(const double x,
		   const double varpi,
		   const double varpi_t,
		   const int maxiter,
		   const double epsabs,
		   const double epsrel) const {
  // To understand what this function does, note that
  //
  // u^2 = ur^2 * (1 - varpi^2/a^2),
  //
  // so the maximum of u^2 (or u) occurs at the radius a that
  // satisfies
  //
  // 0 = d/da [ur^2 * (1 - varpi^2/a^2)]
  // 0 = (1 - varpi^2/a^2) * d/da ur^2 + varpi^2 / (2 a^3) * ur^2

  // Prepare data for GSL
  struct a_max_u_params par;
  par.x = x;
  par.vp2 = SQR(varpi) + SQR(varpi_t);
  par.pw = this;
  gsl_function F;
  F.function = &a_max_u_loga;
  F.params = &par;

  // Specify that we want to use the Brent solver, and allocate its
  // workspace
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  // Set brackets
  double loga_lo = par.vp2 > 1.0 ? 0.5*log(par.vp2) : 1.0e-10;
  double loga_hi = log(amax_abs);
  gsl_root_fsolver_set(s, &F, loga_lo, loga_hi);

  // Loop
  double loga = 0.0;
  for (int iter=0, status=GSL_CONTINUE;
       (iter<maxiter) && (status==GSL_CONTINUE);
       iter++) {
    status = gsl_root_fsolver_iterate(s);
    loga = gsl_root_fsolver_root(s);
    loga_lo = gsl_root_fsolver_x_lower(s);
    loga_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(loga_lo, loga_hi, epsabs, epsrel);
  }

  // Free workspace and return result
  gsl_root_fsolver_free(s);
  return exp(loga);  
}

struct a_max_umax_params {
  double vp2;
  const pwind_rad_isothermal *pw;
};

static double a_max_umax_loga(double loga, void *params) {
  // This routine returns the residual of equation solved by
  // a_max_umax_from_varpi
  struct a_max_umax_params *par = (struct a_max_umax_params *) params;
  double a = exp(loga);
  // Handle the case with a == varpi exactly by hand, to avoid
  // producing infinities that will crash the GSL
  if (SQR(a) == par->vp2) return numeric_limits<double>::max();
  //double ret = par->pw->dU2maxda(a) +
  //  2.0*par->pw->U2max(a) / (a * (SQR(a)/par->vp2 - 1.0));
  double ret = (1.0 - par->vp2/SQR(a)) * par->pw->dU2maxda(a) +
    par->vp2 / (a*SQR(a)) * par->pw->U2max(a);
  return ret;
}

double
pwind_rad_isothermal::
a_max_umax_from_varpi(const double varpi,
		      const double varpi_t,
		      const int maxiter,
		      const double epsabs,
		      const double epsrel) const {
  // Same as a_max_u_from_varpi, but with umax in place of u

  // Trivial case
  double vp2 = SQR(varpi) + SQR(varpi_t);
  if (vp2 == 0.0) return(a_maxu);

  // Prepare data for GSL
  struct a_max_umax_params par;
  par.vp2 = vp2;
  par.pw = this;
  gsl_function F;
  F.function = &a_max_umax_loga;
  F.params = &par;

  // Specify that we want to use the Brent solver, and allocate its
  // workspace
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  // Set brackets
  double loga_lo = vp2 > 1.0 ? 0.5*log(vp2) : 0.0;
  double loga_hi = log(amax_abs);
  gsl_root_fsolver_set(s, &F, loga_lo, loga_hi);

  // Loop
  double loga = 0.0;
  for (int iter=0, status=GSL_CONTINUE;
       (iter<maxiter) && (status==GSL_CONTINUE);
       iter++) {
    status = gsl_root_fsolver_iterate(s);
    loga = gsl_root_fsolver_root(s);
    loga_lo = gsl_root_fsolver_x_lower(s);
    loga_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(loga_lo, loga_hi, epsabs, epsrel);
  }

  // Free workspace and return result
  gsl_root_fsolver_free(s);
  return exp(loga);  
}


////////////////////////////////////////////////////////////////////////
// Point potentials specialised by expansion factor
////////////////////////////////////////////////////////////////////////

// Point, area
double
pwind_rad_pa::U2(const double x, const double a) const {
  return (Gamma*fac1(x)-1.0) * (1.0-1.0/a);
}
double
pwind_rad_pa::dU2dx(const double x, const double a) const {
  return Gamma * (tau0*exp(-tau0*exp(x)) - fac1(x)) *
    (1.0-1.0/a);
}
vector<double>
pwind_rad_pa::alimits(const double u,
		      const double varpi,
		      const double varpi_t,
		      const double epsabs,
		      const double epsrel) const {
  vector<double> ret;
  if (fabs(u) > umax) return ret;  // No solutions
  ret.resize(2);    // 2 solutions
  ret[0] = a_from_u_max(u, varpi, varpi_t);
  ret[1] = numeric_limits<double>::max();
  return ret;
}
double
pwind_rad_pa::Phi_c(const double u,
		    const double fw,
		    const double varpi,
		    const double varpi_t,
		    const double a0,
		    const double a1,
		    const double epsabs,
		    const double epsrel,
		    const vector<double>& alimits_) const {
  if (a1 >= numeric_limits<double>::max())
    return numeric_limits<double>::max();
  else
    return pwind::Phi_c(u, fw, varpi, varpi_t, a0, a1, epsabs, epsrel,
			alimits_);
}
vector<double>
pwind_rad_pa::Phi_c(const vector<double>& u,
		    const double fw,
		    const double varpi,
		    const double varpi_t,
		    const double a0,
		    const double a1,
		    const double epsabs,
		    const double epsrel,
		    const vector<double>& alimits_) const {
  // Output holder
  vector<double> result(u.size());

  // Parallel loop over inputs
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (std::vector<double>::size_type i=0; i<u.size(); i++)
    result[i] = Phi_c(u[i], fw, varpi, varpi_t, a0, a1, epsabs, epsrel,
		      alimits_);

  // Return
  return result;
}

// Point, intermediate
pwind_rad_pi::
pwind_rad_pi(const double Gamma_,
	     const double mach_,
	     const double tau0_, 
	     const pwind_geom *geom_,
	     const double fcrit_,
	     const double jsp_) :
    pwind_rad_point(Gamma_, mach_, tau0_,
		    static_cast<const pwind_expansion *>
		    (&pwind_expansion_intermediate), geom_,
		    fcrit_, jsp_) {
  u_min_infinity = sqrt(U2(xcrit, numeric_limits<double>::max()));
}
double
pwind_rad_pi::U2(const double x, const double a) const {
  return Gamma*fac2(x, a)-1.0+1.0/a;
}
double
pwind_rad_pi::dU2dx(const double x, const double a) const {
  return Gamma*(fac5(x, a) - fac2(x, a));
}
vector<double>
pwind_rad_pi::alimits(const double u,
		      const double varpi,
		      const double varpi_t,
		      const double epsabs,
		      const double epsrel) const {
  vector<double> ret;
  if (fabs(u) > umax) return ret;  // No solutions; u too big
  ret.resize(2); // Two solutions
  ret[0] = a_from_u_max(u, varpi, varpi_t);
  if (fabs(u) > u_min_infinity)
    // Velocity is larger than minimum velocity at infinity, so upper
    // limit on a = infinity
    ret[1] = numeric_limits<double>::max();
  else
    // u < minimum velocity at infinity, so there is a solution at
    // finite a
    ret[1] = a_from_u_x(u, xcrit, varpi, varpi_t, epsabs, epsrel);
  return ret;
}

// Point, solid angle
pwind_rad_ps::
pwind_rad_ps(const double Gamma_,
	     const double mach_,
	     const double tau0_, 
	     const pwind_geom *geom_,
	     const double fcrit_,
	     const double jsp_) :
    pwind_rad_point(Gamma_, mach_, tau0_,
		    static_cast<const pwind_expansion *>
		    (&pwind_expansion_solid_angle), geom_,
		    fcrit_, jsp_) {
  u_min_infinity = sqrt(U2(xcrit, numeric_limits<double>::max()));
}
double
pwind_rad_ps::U2(const double x, const double a) const {
  return Gamma*(fac3(x, a) + fac4(x, a) - fac1(x)) - 1.0 + 1.0/a;
}
double
pwind_rad_ps::dU2dx(const double x, const double a) const {
    return Gamma*(fac1(x) - fac3(x, a) - 0.5*fac4(x, a));
}
vector<double>
pwind_rad_ps::alimits(const double u,
		      const double varpi,
		      const double varpi_t,
		      const double epsabs,
		      const double epsrel) const {
  vector<double> ret;
  if (fabs(u) > umax) return ret;  // No solutions; u too big
  ret.resize(2); // Two solutions
  ret[0] = a_from_u_max(u, varpi, varpi_t);
  if (fabs(u) > u_min_infinity)
    // Velocity is larger than minimum velocity at infinity, so upper
    // limit on a = infinity
    ret[1] = numeric_limits<double>::max();
  else
    // u < minimum velocity at infinity, so there is a solution at
    // finite a
    ret[1] = a_from_u_x(u, xcrit, varpi, varpi_t, epsabs, epsrel);
  return ret;
}

////////////////////////////////////////////////////////////////////////
// Isothermal potentials specialised by expansion factors
////////////////////////////////////////////////////////////////////////

// Isothermal, area
double
pwind_rad_ia::U2(const double x, const double a) const {
  return Gamma*fac1(x) * (1.0-1.0/a) - log(a);
}
double
pwind_rad_ia::dU2dx(const double x, const double a) const {
  return (1.0-1.0/a) * Gamma * (tau0*exp(-tau0*exp(x)) - fac1(x));
}
vector<double>
pwind_rad_ia::alimits(const double u,
		      const double varpi,
		      const double varpi_t,
		      const double epsabs,
		      const double epsrel) const {
  vector<double> ret;
  if (fabs(u) > umax) return ret;  // No solutions; u too big
  if (varpi == 0.0 && varpi_t == 0.0) {
    // For varpi = 0, looking down the barrel, two solutions exist
    // if u < umax
    ret.resize(2);
    ret[0] = a_from_u_max(u, varpi, varpi_t, 1.0, a_maxu,
			  100, epsabs, epsrel);
    ret[1] = a_from_u_max(u, varpi, varpi_t, a_maxu, amax_abs,
			  100, epsabs, epsrel);
    return ret;
  } else {
    // For varpi != 0, we must first find the radius at which the line
    // of sight velocity reaches its maximum, and the velocity at that
    // radus, to check if a solution exists
    double vp2 = SQR(varpi) + SQR(varpi_t);
    double apeak = a_max_umax_from_varpi(varpi, varpi_t,
					 100, epsabs, epsrel);
    double u2peak = U2max(apeak) * (1.0 - vp2/SQR(apeak));
    if (SQR(u) > u2peak) return ret; // No solutions; u too big
    ret.resize(2);
    ret[0] = a_from_u_max(u, varpi, varpi_t, 1.0, apeak,
			  100, epsabs, epsrel);
    ret[1] = a_from_u_max(u, varpi, varpi_t, apeak, amax_abs,
			  100, epsabs, epsrel);
    return ret;
  }
}
vector<double>
pwind_rad_ia::xlimits(const double a) const {
  vector<double> ret;
  if (a >= amax_abs) return ret; // No solutions; a too big
  ret.resize(2);
  ret[0] = -numeric_limits<double>::max();
  ret[1] = X(0.0, a);
  if (ret[1] == numeric_limits<double>::max()) {
    // This is a flag value, indicating solver failed, so return no range
    ret.resize(0);
  }
  return ret;
}
double
pwind_rad_ia::amax(const double x,
		   const double epsabs,
		   const double epsrel) const {
  double amx = a_from_u_x(0.0, x, 0.0, 1.0+1.0e-8, amax_abs, 100,
			  epsabs, epsrel);
  if (amx > 0.) return amx;
  else return 1.0+1.0e-8; // Safety value
}

// Intermediate, isothermal
double
pwind_rad_ii::U2(const double x, const double a) const {
  return Gamma*fac2(x, a)-log(a);
}
double
pwind_rad_ii::dU2dx(const double x, const double a) const {
  return Gamma*(fac5(x, a) - fac2(x, a));
}
vector<double>
pwind_rad_ii::alimits(const double u,
		      const double varpi,
		      const double varpi_t,
		      const double epsabs,
		      const double epsrel) const {
  vector<double> ret;
  if (fabs(u) > umax) return ret;  // No solutions; u too big
  if (varpi == 0.0 && varpi_t == 0.0) {
    // For varpi = 0, looking down the barrel, two solutions exist
    // if u < umax
    ret.resize(2);
    ret[0] = a_from_u_max(u, varpi, varpi_t, 1.0, a_maxu,
			  100, epsabs, epsrel);
    ret[1] = a_from_u_max(u, varpi, varpi_t, a_maxu, amax_abs,
			  100, epsabs, epsrel);
    return ret;
  } else {
    // For varpi != 0, we must first find the radius at which the line
    // of sight velocity reaches its maximum, and the velocity at that
    // radus, to check if a solution exists
    double vp2 = SQR(varpi) + SQR(varpi_t);
    double apeak = a_max_umax_from_varpi(varpi, varpi_t, 100,
					 epsabs, epsrel);
    double u2peak = U2max(apeak) / (1.0 - vp2/SQR(apeak));
    if (SQR(u) > u2peak) return ret; // No solutions; u too big
    ret.resize(2);
    ret[0] = a_from_u_max(u, varpi, varpi_t, 1.0, apeak,
			  100, epsabs, epsrel);
    ret[1] = a_from_u_max(u, varpi, varpi_t, apeak, amax_abs,
			  100, epsabs, epsrel);
    return ret;
  }
}
vector<double>
pwind_rad_ii::xlimits(const double a) const {
  vector<double> ret;
  if (a >= amax_abs) return ret; // No solutions; a too big
  ret.resize(2);
  ret[0] = -numeric_limits<double>::max();
  ret[1] = X(0.0, a);
  if (ret[1] == numeric_limits<double>::max()) {
    // This is a flag value, indicating solver failed, so return no range
    ret.resize(0);
  }
  return ret;
}
double
pwind_rad_ii::amax(const double x,
		   const double epsabs,
		   const double epsrel) const {
  double amx = a_from_u_x(0.0, x, 0.0, 1.0+1.0e-8, amax_abs, 100,
			  epsabs, epsrel);
  if (amx > 0.) return amx;
  else return 1.0+1.0e-8; // Safety value
}

// Isothermal, constant solid angle
pwind_rad_is::pwind_rad_is(const double Gamma_,
			   const double mach_,
			   const double tau0_, 
			   const pwind_geom *geom_,
			   const double fcrit_,
			   const double jsp_) :
  pwind_rad_isothermal(Gamma_, mach_, tau0_,
		       static_cast<const pwind_expansion *>
		       (&pwind_expansion_solid_angle), geom_,
		       fcrit_, jsp_)
{
  // Do these computations with very high precision, since they are
  // only done once

  // Solve for radius of maximum velocity along x = xcrit curve
  a_maxu_xcrit = a_max_u_from_varpi(xcrit, 0.0, 0.0, 100,
				    1e-12, 1e-12);
  // Next get velocity maximum at x = xcrit
  umax_xcrit = sqrt(U2(xcrit, a_maxu_xcrit));
  // Get radius at which a = 0 along x = xcrit curve
  amax_xcrit = a_from_u_x(0.0, xcrit, 0.0, 0.0, a_maxu_xcrit,
			  numeric_limits<double>::max(),
			  500, 1e-12, 1e-12);
}
double
pwind_rad_is::U2(const double x, const double a) const {
  double ret = Gamma*(fac3(x, a) + fac4(x, a) - fac1(x)) - log(a);
  return ret;
}
double
pwind_rad_is::dU2dx(const double x, const double a) const {
  return Gamma*(fac1(x) - fac3(x, a) - 0.5*fac4(x, a));
}
vector<double>
pwind_rad_is::alimits(const double u,
		      const double varpi,
		      const double varpi_t,
		      const double epsabs,
		      const double epsrel) const {
  vector<double> ret;
  if (fabs(u) > umax) return ret;  // No solutions; u too big
  if (varpi == 0.0 && varpi_t == 0.0) {
    // For varpi = 0, looking down the barrel, there can be 2 or 4
    // solutions
    if (fabs(u) < umax_xcrit) {
      // Four solutions, which are as follows:
      // 1. Radius at which u = U(-inf, a), in the region of a where U
      // is increasing
      // 2. Radius at which u = U(xcrit, a), in the region of a where
      // U is increasing
      // 3. Radius at which u = U(xcrit, a), in the region of a where
      // U is increasing
      // 4. Radius at which u = U(-inf, a), in the region of a where U
      // is decreasing
      ret.resize(4);
      ret[0] = a_from_u_max(u, varpi, varpi_t, 1.0, a_maxu, 100, epsabs,
			    epsrel);
      ret[1] = a_from_u_x(u, xcrit, varpi, varpi_t, ret[0], a_maxu_xcrit,
			  100, epsabs, epsrel);
      ret[2] = a_from_u_x(u, xcrit, varpi, varpi_t,
			  max(a_maxu_xcrit, ret[1]), amax_xcrit,
			  100, epsabs, epsrel);
      ret[3] = a_from_u_max(u, varpi, varpi_t, ret[2], amax_abs,
			    100, epsabs, epsrel);
      return ret;
    } else {
      // Two solutions; only solutions 1. and 4. of those listed above
      // exist
      ret.resize(2);
      ret[0] = a_from_u_max(u, varpi, varpi_t, 1.0, a_maxu,
			    100, epsabs, epsrel);
      ret[1] = a_from_u_max(u, varpi, varpi_t, a_maxu, amax_abs,
			    100, epsabs, epsrel);
      return ret;
    }
  } else {
    // For varpi != 0, we must first find the radii at which the
    // velocities along x = -inf and x = xcrit reach their maxima, and
    // then find the velocities at those radii. This will determine
    // the number of solutions we have.
    
    // Step 1: find the radius of maximum velocity along the x = -inf
    // curve
    double vp2 = SQR(varpi) + SQR(varpi_t);
    double apeak_xinf = a_max_umax_from_varpi(varpi, varpi_t,
					      100, epsabs, epsrel);
    double u2peak_xinf = U2max(apeak_xinf) /
      (1.0 - vp2/SQR(apeak_xinf));
    if (SQR(u) > u2peak_xinf) return ret; // No solutions; u too big

    // Step 2: find the radius of maximum velocity along the x = xcrit
    // curve
    double apeak_xcrit = a_max_u_from_varpi(xcrit, varpi, varpi_t,
					    100, epsabs, epsrel);
    double u2peak_xcrit = U2(xcrit, apeak_xcrit) /
      (1.0 - vp2/SQR(apeak_xcrit));

    // Step 3: decide how many solutions there are
    if (SQR(u) < u2peak_xcrit) {
      // Step 3a: there are 4 solutions, analogous to the u <
      // umax_xcrit case for varpi = 0
      ret.resize(4);
      ret[0] = a_from_u_max(u, varpi, varpi_t, 1.0, apeak_xinf,
			    100, epsabs, epsrel);
      ret[1] = a_from_u_x(u, xcrit, varpi, varpi_t, ret[0], apeak_xcrit,
			  100, epsabs, epsrel);
      ret[2] = a_from_u_x(u, xcrit, varpi, varpi_t,
			  max(apeak_xcrit, ret[1]), amax_xcrit,
			  100, epsabs, epsrel);
      ret[3] = a_from_u_max(u, varpi, varpi_t, ret[2], amax_abs,
			    100, epsabs, epsrel);
      return ret;
    } else {
      // Step 3b, there are 2 solutions, analogous to the u >=
      // umax_xcrit case for varpi = 0
      ret.resize(2);
      ret[0] = a_from_u_max(u, varpi, varpi_t, 1.0, apeak_xinf,
			    100, epsabs, epsrel);
      ret[1] = a_from_u_max(u, varpi, varpi_t, apeak_xinf, amax_abs,
			    100, epsabs, epsrel);
      return ret;
    }
  }
}
vector<double>
pwind_rad_is::xlimits(const double a) const {
  vector<double> ret;
  if (a >= amax_abs) return ret; // No solutions; a too big
  ret.resize(2);
  ret[0] = -numeric_limits<double>::max();
  if (a < amax_xcrit)
    ret[1] = xcrit;
  else {
    ret[1] = X(0.0, a);
    if (ret[1] == numeric_limits<double>::max()) {
      // This is a flag value, indicating solver failed, so return no range
      ret.resize(0);
    }
  }
  return ret;
}
double
pwind_rad_is::amax(const double x,
		   const double epsabs,
		   const double epsrel) const {
  double amx = a_from_u_x(0.0, x, 0.0, 0.0, 1.0+1.0e-8, amax_abs,
			  100, epsabs, epsrel);
  if (amx > 0.) return amx;
  else return 1.0+1.0e-8; // Safety value
}
