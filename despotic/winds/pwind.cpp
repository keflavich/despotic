// Momentum-driven wind class

#include <cmath>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include "pwind.H"
#if defined(_OPENMP)
#  include "omp.h"
#endif

using namespace std;

////////////////////////////////////////////////////////////////////////
// General note here: we avoid using infinities because we cannot rely
// on c++ compilers to accurately implement and handle IEEE
// infinities, particularly if we turn on optimisation. This results
// in ridiculous outcomes on some compilers. For example, when the
// following code is compiled with g++ version 5 with -Ofast,
//
// double a = std::numeric_limits<double>::infinity()
// std::cout << std::isfinite(a) << std::endl;
//
// the resulting program prints 1, i.e., isfinite of infinity returns
// True! Given this bug, we instead flag infinities by making them
// equal to the maximum finite number. We only use infinities when
// returning things to python through the c inferface.
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Data structures required by the GSL
////////////////////////////////////////////////////////////////////////

struct U2_params {
  double u, a, x, varpi, varpi_t, vp2;
  const pwind *pw;
};

struct Phi_params {
  double u, varpi, varpi_t, vp2, jproj;
  double tXtw, fj, boltzfac;  // This need only be initialised for tau_c
  double tau_target; // Only used by tau_c_func
  double fw; // Only used in the correlated case
  vector<double> alim; // Only used by tau_c
  double epsabs, epsrel; // Only used by tau_c
  const pwind *w;
};

struct Phi_params_vec {
  double u, varpi, varpi_t, vp2, jproj, fj, boltzfac, epsabs, epsrel;
  vector<double> u_trans, tXtw;
  double tau_target;
  double fw;
  vector< vector<double> > alim_k;
  const pwind *w;
};

struct Xi_params {
  double u, vp2;
  const pwind *w;
};

struct xi_params {
  double a, vp2;
  double epsabs, epsrel;
  int gsl_err;
  const pwind *w;
};

struct eta_params {
  double u, vp2, loga_max;
  bool thin;
  gsl_spline *tau_spl;
  gsl_interp_accel *tau_acc;
  const pwind *w;
};

struct tau_params {
  double loga, fw, fc0;
  gsl_spline *tau_spl;
  gsl_interp_accel *tau_acc;
  const pwind *w;
};

struct Psi_params {
  double tXtw, fj, boltzfac, fw, varpi, varpi_t, epsabs, epsrel;
  bool correlated, thin;
  const pwind *w;
};

////////////////////////////////////////////////////////////////////////
// Constructor; note that we turn off GSL error handling here, because
// we will handle errors manually.
////////////////////////////////////////////////////////////////////////
pwind::pwind(const double Gamma_,
	     const double mach_,
	     const pwind_potential *potential_,
	     const pwind_expansion *expansion_,
	     const pwind_geom *geom_,
	     const double fcrit_,
	     const double jsp_) :
    potential(potential_),
    expansion(expansion_),
    geom(geom_),
    Gamma(Gamma_),
    mach(mach_),
    fcrit(fcrit_),
    jsp(jsp_),
    gsl_err_stat(GSL_SUCCESS),
    err_handler(gsl_set_error_handler_off()),
    sx(sxMach(mach))
{ };


////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Numerical inversions of the wind acceleration law
////////////////////////////////////////////////////////////////////////
static double U2_func_loga(double loga, void *params) {
  // This routine returns the residual for U2(x, a) - U2_target
  struct U2_params *par = (struct U2_params *) params;
  double a = exp(loga);
  double jproj = par->pw->getJproj();
  double varpi_a = par->varpi;
  double varpi_t = par->varpi_t;
  double vp2 = par->vp2;
  double u = par->u;
  double U2 = SQR(u - jproj*varpi_t / (SQR(a) - SQR(varpi_a)))
    / (1.0 - vp2/SQR(a));
  // Catch roundoff problems
  if (U2 < 0.0)
    return numeric_limits<double>::max();
  double ret = U2 - par->pw->U2(par->x, a);
  // GSL doesn't like it if we return infinities
  if (ret > numeric_limits<double>::max())
    ret = numeric_limits<double>::max();
  return ret;
}

double pwind::a_from_u_x(const double u,
			 const double x,
			 const double varpi,
			 const double varpi_t,
			 const double alo,
			 const double ahi,
			 const int maxiter,
			 const double epsabs,
			 const double epsrel) const {

  // Prepare data for the GSL
  struct U2_params params;
  params.u = u;
  params.x = x;
  params.varpi = varpi;
  params.varpi_t = varpi_t;
  params.vp2 = SQR(varpi) + SQR(varpi_t);
  params.pw = this;
  gsl_function F;
  F.function = &U2_func_loga;
  F.params = &params;

  // Specify that we want to use the Brent solver, and allocate its
  // workspace
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  // Ensure that the root is properly bracketed; return -1 if not
  double loga_lo = SQR(alo) < params.vp2 ? 0.5*log(params.vp2) : log(alo);
  double loga_hi = log(ahi);
  double sgn1 = U2_func_loga(loga_lo, &params);
  double sgn2 = U2_func_loga(loga_hi, &params);
  if (sgn1*sgn2 > 0){
    gsl_root_fsolver_free(s);
    return -1;
  }

  // Set brackets
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

static double U2_func_x(double x, void *params) {
  // This routine returns the residual for U2(x, a) - U2_target
  struct U2_params *par = (struct U2_params *) params;
  double U2 = SQR(par->u) / (1.0 - par->vp2/SQR(par->a));
  return U2 - par->pw->U2(x, par->a);
}

double pwind::x_from_u_a(const double u,
			 const double a,
			 const double varpi,
			 const double varpi_t,
			 const double xlo,
			 const double xhi,
			 const int maxiter,
			 const double epsabs,
			 const double epsrel) const {
  
  // Prepare data for the GSL
  struct U2_params params;
  params.u = u;
  params.a = a;
  params.vp2 = SQR(varpi) + SQR(varpi_t);
  params.pw = this;
  gsl_function F;
  F.function = &U2_func_x;
  F.params = &params;

  // Ensure that the root is properly bracketed; return sentinel value
  // if not
  double sgn1 = U2_func_x(xlo, &params);
  double sgn2 = U2_func_x(min(xhi, xcrit), &params);
  if (sgn1*sgn2 > 0) {
    if (fabs(sgn1) < epsabs) return xlo;
    else if (fabs(sgn2) < epsabs) return min(xhi, xcrit);
    else return -numeric_limits<double>::max();
  }

  // Specify that we want to use the Brent solver, and allocate its
  // workspace
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  // Set brakets
  gsl_root_fsolver_set(s, &F, xlo, min(xhi, xcrit));

  // Loop
  double x = 0.0;
  for (int iter=0, status=GSL_CONTINUE;
       (iter<maxiter) && (status==GSL_CONTINUE);
       iter++) {
    status = gsl_root_fsolver_iterate(s);
    x = gsl_root_fsolver_root(s);
    double xl = gsl_root_fsolver_x_lower(s);
    double xh = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(xl, xh, epsabs, epsrel);
  }

  // Free workspace and return result
  gsl_root_fsolver_free(s);
  return x;
}


////////////////////////////////////////////////////////////////////////
// Wind differential and total density
////////////////////////////////////////////////////////////////////////

struct rho_integ_params {
  double a;
  const pwind *pw;
};

static double rho_integ(double x, void *params) {
  struct rho_integ_params *par = (struct rho_integ_params *) params;
  double pMval = pM(x, par->pw->s());
  double U = sqrt(par->pw->U2(x, par->a));
  if (pMval == 0) return 0.0;
  else return pMval / (par->a * par->a * U);
}

double pwind::rho(const double a,
		  const double epsabs,
		  const double epsrel) const {
  
  // Load structures for the GSL */
  struct rho_integ_params par;
  par.a = a;
  par.pw = this;
  gsl_function F;
  F.function = &rho_integ;
  F.params = &par;

  // Get integration limits
  vector<double> xlim = xlimits(a);

  // Set up integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);

  // Integrate
  double result, err;
  if (xlim[0] == -numeric_limits<double>::max()) {
    checkForErr(gsl_integration_qagil(&F, xlim[1]+log(fcrit), epsabs, epsrel,
				      MAXINTERVAL, w, &result, &err));
  } else {
    if (xlim[1] - log(fcrit) > xlim[0])
      checkForErr(gsl_integration_qag(&F, xlim[0], xlim[1]+log(fcrit),
				      epsabs, epsrel,
				      MAXINTERVAL, GSL_INTEG_GAUSS61,
				      w, &result, &err));
    else {
      result = 0.0;
    }
  }

  // Take down integration machinery
  gsl_integration_workspace_free(w);

  // Apply prefactor and return
  result /= zeta_M;
  return result;  
}


// Vectorized version of function above
vector<double> pwind::rho(const vector<double>& a,
			  const double epsabs,
			  const double epsrel) const {
  
  vector<double> result(a.size());
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (vector<double>::size_type i=0; i<a.size(); i++)
    result[i] = rho(a[i], epsabs, epsrel);
  return result;
}
 

////////////////////////////////////////////////////////////////////////
// Total momentum flux
////////////////////////////////////////////////////////////////////////
struct pdot_integ_params {
  double a;
  const pwind *pw;
};

static double pdot_integ(double x, void *params) {
  struct pdot_integ_params *par = (struct pdot_integ_params *) params;
  double pMval = pM(x, par->pw->s());
  double U = sqrt(par->pw->U2(x, par->a));
  if (pMval == 0) return 0.0;
  else return pMval * U;
}

double pwind::pdot(const double a,
		   const double epsabs,
		   const double epsrel) const {
  // Load structures for the GSL */
  struct pdot_integ_params par;
  par.a = a;
  par.pw = this;
  gsl_function F;
  F.function = &pdot_integ;
  F.params = &par;

  // Get integration limits
  vector<double> xlim = xlimits(a);

  // Set up integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);

  // Integrate
  double result, err;
  if (xlim[0] == -numeric_limits<double>::max()) {
    checkForErr(gsl_integration_qagil(&F, xlim[1]+log(fcrit), epsabs, epsrel,
				      MAXINTERVAL, w, &result, &err));
  } else {
    if (xlim[1] - log(fcrit) > xlim[0])
      checkForErr(gsl_integration_qag(&F, xlim[0], xlim[1]+log(fcrit),
				      epsabs, epsrel,
				      MAXINTERVAL, GSL_INTEG_GAUSS61,
				      w, &result, &err));
    else {
      result = 0.0;
    }
  }

  // Take down integration machinery
  gsl_integration_workspace_free(w);

  // Multiply by prefactor and return
  result /= zeta_M;
  return result;
}

// Vectorized version of function above
vector<double> pwind::pdot(const vector<double>& a,
			   const double epsabs,
			   const double epsrel) const {
  
  vector<double> result(a.size());
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (vector<double>::size_type i=0; i<a.size(); i++)
    result[i] = pdot(a[i], epsabs, epsrel);
  return result;
}


////////////////////////////////////////////////////////////////////////
// Total energy flux
////////////////////////////////////////////////////////////////////////
struct Edot_integ_params {
  double a;
  const pwind *pw;
};

static double Edot_integ(double x, void *params) {
  struct Edot_integ_params *par = (struct Edot_integ_params *) params;
  double pMval = pM(x, par->pw->s());
  double U2 = par->pw->U2(x, par->a);
  if (pMval == 0) return 0.0;
  else return pMval * U2;
}

double pwind::Edot(const double a,
		   const double epsabs,
		   const double epsrel) const {
  // Load structures for the GSL */
  struct Edot_integ_params par;
  par.a = a;
  par.pw = this;
  gsl_function F;
  F.function = &Edot_integ;
  F.params = &par;

  // Get integration limits
  vector<double> xlim = xlimits(a);

  // Set up integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);

  // Integrate
  double result, err;
  if (xlim[0] == -numeric_limits<double>::max()) {
    checkForErr(gsl_integration_qagil(&F, xlim[1]+log(fcrit), epsabs, epsrel,
				      MAXINTERVAL, w, &result, &err));
  } else {
    if (xlim[1] - log(fcrit) > xlim[0])
      checkForErr(gsl_integration_qag(&F, xlim[0], xlim[1]+log(fcrit),
				      epsabs, epsrel,
				      MAXINTERVAL, GSL_INTEG_GAUSS61,
				      w, &result, &err));
    else {
      result = 0.0;
    }
  }

  // Take down integration machinery
  gsl_integration_workspace_free(w);

  // Multiply by prefactor and return
  result /= 2.0*zeta_M;
  return result;
}

// Vectorized version of function above
vector<double> pwind::Edot(const vector<double>& a,
			   const double epsabs,
			   const double epsrel) const {
  
  vector<double> result(a.size());
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (vector<double>::size_type i=0; i<a.size(); i++)
    result[i] = Edot(a[i], epsabs, epsrel);
  return result;
}

////////////////////////////////////////////////////////////////////////
// Total momentum flux normalised to driving momentum flux
////////////////////////////////////////////////////////////////////////

double pwind::pdotRel(const double a,
		      const double epsabs,
		      const double epsrel) const {

  // Load structures for the GSL */
  struct pdot_integ_params par;
  par.a = a;
  par.pw = this;
  gsl_function F;
  F.function = &pdot_integ;
  F.params = &par;

  // Get integration limits
  vector<double> xlim = xlimits(a);

  // Set up integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);

  // Integrate
  double result, err;
  if (xlim[0] == -numeric_limits<double>::max()) {
    checkForErr(gsl_integration_qagil(&F, xlim[1]+log(fcrit), epsabs, epsrel,
				      MAXINTERVAL, w, &result, &err));
  } else {
    if (xlim[1] - log(fcrit) > xlim[0])
      checkForErr(gsl_integration_qag(&F, xlim[0], xlim[1]+log(fcrit),
				      epsabs, epsrel,
				      MAXINTERVAL, GSL_INTEG_GAUSS61,
				      w, &result, &err));
    else {
      result = 0.0;
    }
  }

  // Take down integration machinery
  gsl_integration_workspace_free(w);

  // Multiply by prefactor and return
  result /= Gamma;
  return result;
}
  
double
pwind::pdotRel(const double a,
	       const double fg,
	       const double tctw,
	       const double epsabs,
	       const double epsrel) const {
  double pd = pdot(a, epsabs, epsrel);
  return 2.0*tctw*pd/(3.0*fg*zeta_M);
}

////////////////////////////////////////////////////////////////////////
// Absorption calculations
////////////////////////////////////////////////////////////////////////

// Static functions that are not part of the class, and are instead
// passed to the GSL integrator

static double Phi_uc_integ(double loga, void *params) {
  struct Phi_params *par = (struct Phi_params *) params;

  // Get a, and check if it is so large or small that we should return
  // zero
  double a = exp(loga);
  if (a >= numeric_limits<double>::max()) return 0.0;
  if (par->vp2 >= SQR(a)) return 0.0;

  // Get x and dU2/dx
  double ur = fabs(par->u -
		   par->jproj * par->varpi_t /
		   (a*a - par->varpi*par->varpi)) /
    sqrt(1.0 - par->vp2/SQR(a));
  if (ur > par->w->uMax()) return 0.0;
  double x = par->w->X(ur, a);

  // Return if x is too large, or if it is -inf -- the latter is
  // possible if we invoke this function with values of u and a near
  // the limit of where solutions exist, and for numerical reasons we
  // fail to find the solution we're looking for
  if (x > par->w->xcr() + log(par->w->getFcrit())) return 0.0;
  if (x == -numeric_limits<double>::max()) return 0.0;

  // Get dU2dx
  double dU2dx = par->w->dU2dx(x, a);

  // Return
  double ret = a*pM(x, par->w->s()) /
    (fabs(dU2dx)*(SQR(a)-par->vp2));
  return ret;
}

static double Phi_c_integ(double loga, void *params) {
  struct Phi_params *par = (struct Phi_params *) params;

  // Get a, and check if it is so large or small that we should return
  // zero
  double a = exp(loga);
  if (!isfinite(a)) return 0.0;
  if (par->vp2 >= SQR(a)) return 0.0;

  // Get x and dU2/dx
  double ur = fabs(par->u) / sqrt(1.0 - par->vp2/SQR(a));
  if (ur > par->w->uMax()) return 0.0;
  double x = par->w->X(ur, a);
  if (x > par->w->xcr() + log(par->w->getFcrit())) return 0.0;
  double dU2dx = par->w->dU2dx(x, a);

  // Return; catch bad values
  if (dU2dx == 0.0) return 0.0;
  double ret = a*pM(x, par->w->s()) /
    (par->w->y(a)*fabs(dU2dx)*(1.0-par->vp2/SQR(a)));
  return ret;
}

static double tau_c_integ(double loga, void *params) {
  struct Phi_params *par = (struct Phi_params *) params;

  // Get a; catch invalid values
  double a = exp(loga);
  if (a > numeric_limits<double>::max()) return 0.0;
  if (SQR(a) < par->vp2) return 0.0;

  // Get optical depth
  double tau = par->tXtw * par->fj * (1.0 - par->boltzfac) *
    par->w->Phi_c(par->u, par->fw, sqrt(par->vp2), 0.0, 1.0, a,
		  par->epsabs, par->epsrel, par->alim);

  // Return integrand
  double ret = a*fabs(par->w->dfcda(a, par->fw))*exp(-tau);
  return ret;
}    

static double tau_c_integ_vec(double loga, void *params) {
  struct Phi_params_vec *par = (struct Phi_params_vec *) params;

  // Get a; catch invalid values
  double a = exp(loga);
  if (a > numeric_limits<double>::max()) return 0.0;
  if (SQR(a) < par->vp2) return 0.0;

  // Get optical depth, summing over transitions
  double tau = 0.0;
  for (vector<double>::size_type k=0; k<par->u_trans.size(); k++) {
    if (par->alim_k[k].size() == 0) continue;
    if (a < par->alim_k[k][0]) continue;
    tau += par->tXtw[k] * par->fj * (1.0 - par->boltzfac) *
      par->w->Phi_c(par->u - par->u_trans[k], par->fw, sqrt(par->vp2), 0.0,
		    1.0, a, par->epsabs, par->epsrel, par->alim_k[k]);
  }

  // Return integrand
  double ret = a*fabs(par->w->dfcda(a, par->fw))*exp(-tau);
  return ret;
}    

double pwind::Phi_uc(const double u,
		     const double varpi,
		     const double varpi_t,
		     const double a0,
		     const double a1,
		     const double epsabs,
		     const double epsrel) const {

  // Declare GSL variables
  gsl_function F;
  struct Phi_params par;
  double result, err;

  // Check if velocity is outside bounds
  if (u == 0) return numeric_limits<double>::max();
  if (fabs(u) >= umax) return 0.0;

  // Get integration limits from velocity and geometry; note that we
  // use a much tighter accuracy limit here than we do for the
  // integral, because small errors in finding the limit can lead to
  // numerical failures in the integration
  vector<double> alim_v = alimits(u, varpi, varpi_t,
				  epsabs/1e4, epsrel/1e4);
  vector<double> alim = geom->a_lim(alim_v, varpi, varpi_t, u);
  if (alim.size() == 0) return 0.0;

  // Apply input limits
  for (vector<double>::size_type i=0; i<alim.size()/2; i++) {
    if (alim[2*i] < a0) alim[2*i] = a0;
    if (alim[2*i+1] > a1) alim[2*i+1] = a1;
  }
  
  // Load the GSL data
  par.u = fabs(u);
  par.varpi = varpi;
  par.varpi_t = varpi_t;
  par.vp2 = SQR(varpi) + SQR(varpi_t);
  par.jproj = getJproj();
  par.epsabs = epsabs/10.0;
  par.epsrel = epsrel/10.0;
  par.w = this;
  F.function = &Phi_uc_integ;
  F.params = &par;

  // Set up integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);
  
  // Integrate
  double integ = 0.0;
  for (vector<double>::size_type i=0; i<alim.size()/2; i++) {
    if (alim[2*i+1] <= alim[2*i]) continue;
    if (alim[2*i+1] < numeric_limits<double>::max()) {
      checkForErr(gsl_integration_qag(&F, log(alim[2*i]), log(alim[2*i+1]),
				      epsabs, epsrel, MAXINTERVAL,
				      GSL_INTEG_GAUSS61, w, &result, &err));
    } else {
      checkForErr(gsl_integration_qagiu(&F, log(alim[2*i]), epsabs, epsrel,
					MAXINTERVAL, w, &result, &err));
    }
    integ += result;
  }
  integ /= zeta_M;
  
  // Take down integration machinery
  gsl_integration_workspace_free(w);
  
  // Return
  return integ;
}

// Vectorized version of previous function
vector<double> pwind::Phi_uc(const vector<double> &u,
			     const double varpi,
			     const double varpi_t,
			     const double a0,
			     const double a1,
			     const double epsabs,
			     const double epsrel) const {
  
  // Output holder
  vector<double> result(u.size());

  // Parallel loop over inputs
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (std::vector<double>::size_type i=0; i<u.size(); i++)
    result[i] = Phi_uc(u[i], varpi, varpi_t, a0, a1, epsabs, epsrel);

  // Return
  return result;
}

  
double pwind::tau_uc(const double u,
		     const double tXtw,
		     const double fj,
		     const double boltzfac,
		     const double varpi,
		     const double varpi_t,
		     const double a0,
		     const double a1,
		     const double epsabs,
		     const double epsrel) const {
  return tXtw * fj * (1.0-boltzfac) *
    Phi_uc(u, varpi, varpi_t, a0, a1, epsabs, epsrel);
}


// Vectorized version of previous function
vector<double> pwind::tau_uc(const vector<double>& u,
			     const double tXtw,
			     const double fj,
			     const double boltzfac,
			     const double varpi,
			     const double varpi_t,
			     const double a0,
			     const double a1,
			     const double epsabs,
			     const double epsrel) const {
  
  // Output holder
  vector<double> result(u.size());
  
  // Parallel loop over inputs
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (std::vector<double>::size_type i=0; i<u.size(); i++)
    result[i] = tau_uc(u[i], tXtw, fj, boltzfac, varpi, varpi_t,
		       a0, a1, epsabs, epsrel);

  // Return
  return result;
}


double pwind::tau_uc(const double u,
		     const vector<double>& u_trans,
		     const vector<double>& tXtw,
		     const double fj,
		     const double boltzfac,
		     const double varpi,
		     const double varpi_t,
		     const double a0,
		     const double a1,
		     const double epsabs,
		     const double epsrel) const {
  double sum = 0.0;
  for (vector<double>::size_type i=0; i<tXtw.size(); i++) {
    sum += tXtw[i] * fj * (1.0 - boltzfac) *
      Phi_uc(u - u_trans[i], varpi, varpi_t, a0, a1, epsabs, epsrel);
  }
  return sum;
}

// Vectorized version of previous function
vector<double> pwind::tau_uc(const vector<double>& u,
			     const vector<double>& u_trans,
			     const vector<double>& tXtw,
			     const double fj,
			     const double boltzfac,
			     const double varpi,
			     const double varpi_t,
			     const double a0,
			     const double a1,
			     const double epsabs,
			     const double epsrel) const {

  // Output holder
  vector<double> result(u.size());

  // Parallel loop over inputs
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (std::vector<double>::size_type i=0; i<u.size(); i++)
    result[i] = tau_uc(u[i], u_trans, tXtw, fj, boltzfac, varpi,
		       varpi_t, a0, a1, epsabs, epsrel);

  // Return
  return result;
}


double pwind::Phi_c(const double u,
		    const double fw,
		    const double varpi,
		    const double varpi_t,
		    const double a0,
		    const double a1,
		    const double epsabs,
		    const double epsrel,
		    const vector<double>& alimits_) const {
  gsl_function F;
  struct Phi_params par;
  double result, err;

  // Check if velocity is outside bounds
  if (u == 0) return numeric_limits<double>::max();
  if (fabs(u) >= umax) return 0.0;

  // Get integration limits; note that we use tighter tolerance on the
  // limits than on the integration itself
  vector<double> alim;
  if (alimits_.size() == 0) {
    vector<double> alim_v = alimits(u, varpi, varpi_t, epsabs/1e4, epsrel/1e4);
    alim = geom->a_lim(alim_v, varpi, varpi_t, u);
  } else {
    alim = alimits_;
  }

  // Apply input limits
  for (vector<double>::size_type i=0; i<alim.size()/2; i++) {
    if (alim[2*i] < a0) alim[2*i] = a0;
    if (alim[2*i+1] > a1) alim[2*i+1] = a1;
  }

  // Load the GSL data
  par.u = fabs(u);
  par.vp2 = SQR(varpi) + SQR(varpi_t);
  par.fw = fw;
  par.w = this;
  F.function = &Phi_c_integ;
  F.params = &par;

  // Set up integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);
  
  // Integrate
  double integ = 0.0;
  for (vector<double>::size_type i=0; i<alim.size()/2; i++) {
    if (alim[2*i+1] <= alim[2*i]) continue;
    if (alim[2*i+1] < numeric_limits<double>::max()) {
      checkForErr(gsl_integration_qag(&F, log(alim[2*i]), log(alim[2*i+1]),
				      epsabs, epsrel, MAXINTERVAL,
				      GSL_INTEG_GAUSS61, w, &result, &err));
    } else {
      checkForErr(gsl_integration_qagiu(&F, log(alim[2*i]), epsabs, epsrel,
					MAXINTERVAL, w, &result, &err));
    }
    integ += result;
  }
  integ /= zeta_M;

  // Take down integration machinery
  gsl_integration_workspace_free(w);

  // Apply covering factor correction
  if (fw <= 0) integ /= zeta_A;
  else integ /= fw;

  // Return
  return integ;
}


// Vectorized version of previous function
vector<double> pwind::Phi_c(const vector<double>& u,
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



// Note that b = log(a - 1); we make this change of variable because
// we need to find roots very close to 1 to high accuracy
static double tau_c_func(double b, void *params) {
  struct Phi_params *par = (struct Phi_params *) params;
  double a = 1.0 + exp(b);
  double tau = par->tXtw * par->fj * (1.0 - par->boltzfac) * 
    par->w->Phi_c(par->u, par->fw, sqrt(par->vp2), 0.0, 1.0, a,
		  par->epsabs, par->epsrel, par->alim);
  return tau - par->tau_target;
}

// Version of tau_c_func for multiple transitions
static double tau_c_func_vec(double b, void *params) {
  struct Phi_params_vec *par = (struct Phi_params_vec *) params;
  double a = 1.0 + exp(b);
  double tau = 0.0;
  for (vector<double>::size_type k=0; k<par->u_trans.size(); k++) {
    if (par->alim_k[k].size() == 0) continue;
    if (a < par->alim_k[k][0]) continue;
    tau += par->tXtw[k] * par->fj * (1.0 - par->boltzfac) * 
      par->w->Phi_c(par->u - par->u_trans[k], par->fw, sqrt(par->vp2), 0.0,
		    1.0, a, par->epsabs, par->epsrel, par->alim_k[k]);
  }
  return tau - par->tau_target;
}

double pwind::tau_c(const double u,
		    const double tXtw,
		    const double fj,
		    const double boltzfac,
		    const double fw,
		    const double varpi,
		    const double varpi_t,
		    const double a0,
		    const double a1,
		    const double epsabs,
		    const double epsrel) const {

  // Get integration limits, using tighter tolerance
  vector<double> alim_v = alimits(u, varpi, varpi_t, epsabs/1e4, epsrel/1e4);
  vector<double> alim = geom->a_lim(alim_v, varpi, varpi_t, u);
  if (alim.size() == 0) return 0.0;
  
  // Apply input limits
  for (vector<double>::size_type i=0; i<alim.size()/2; i++) {
    if (alim[2*i] < a0) alim[2*i] = a0;
    if (alim[2*i+1] > a1) alim[2*i+1] = a1;
  }
  
  // Get tau_c
  double fc0;
  if (fw <= 0.0) fc0 = fc(alim[0]);
  else fc0 = fw * y(alim[0]) / SQR(alim[0]);
  if (const_sa()) {
    // Constant solid angle case
    double tau = tXtw * fj * (1.0 - boltzfac) *
      Phi_c(u, fw, varpi, varpi_t, alim[0], alim.back(),
	    epsabs, epsrel, alim);
    return -log(1.0 - fc0 + fc0*exp(-tau));
  } else {
    // Non-constant solid angle; need to integrate

    // Data for GSL
    struct Phi_params par;
    par.u = u;
    par.fw = fw;
    par.vp2 = SQR(varpi) + SQR(varpi_t);
    par.tXtw = tXtw;
    par.fj = fj;
    par.boltzfac = boltzfac;
    par.alim = alim;
    par.epsabs = epsabs / 10.0; // Tighter tolerance for inner
				// integral
    par.epsrel = epsrel / 10.0;
    par.w = this;
    gsl_function F;
    F.function = &tau_c_integ;
    F.params = &par;

    // Set up integration
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);

    // Integrate; this sometimes runs into trouble because the integrand
    // can have a-near step function discontinuity very close to the
    // end of the integration interval, associated with the optical
    // depth tau changing from < 1 to > 1. Integrators do not handle
    // this gracefully. The approach we take is to break the interval
    // up into intervals that bracket tau = 1.
    double b_lo, b_hi;
    if (alim[0] == 1.0) b_lo = -100.0;
    else b_lo = log(alim[0]-1.0);
    b_hi = log(numeric_limits<double>::max());
    vector<double> loga_int_lim;
    loga_int_lim.push_back(log(alim[0]));
    vector<double> tau_break = { 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0 };
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = nullptr;
    gsl_function G;
    for (vector<double>::size_type i=0; i<tau_break.size(); i++) {
      
      // Search for this value of tau
      par.tau_target = tau_break[i];
      if (tau_c_func(b_lo, &par) * tau_c_func(b_hi, &par) > 0) {
	// tau never reaches this value, so stop
	break;
      }

      // Optical depth reaches target, so find out where
      if (s == nullptr) {
	s = gsl_root_fsolver_alloc(T);
	G.function = &tau_c_func;
	G.params = &par;
      }
      gsl_root_fsolver_set(s, &G, b_lo, b_hi);
      double a_break;
      for (int iter=0, status=GSL_CONTINUE;
	   (iter<100) && (status=GSL_CONTINUE);
	   iter++) {
	status = gsl_root_fsolver_iterate(s);
	a_break = 1.0 + exp(gsl_root_fsolver_root(s));
	b_lo = gsl_root_fsolver_x_lower(s);
	b_hi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(b_lo, b_hi, 1e-2, 1e-2);
      }
      loga_int_lim.push_back(log(a_break));

      // Reset search interval
      b_lo = log(a_break-1.0);
      b_hi = log(numeric_limits<double>::max());
    }
    if (s != nullptr) gsl_root_fsolver_free(s);

    // Now integrate over broken-up interval
    double res, result, err;
    result = 0.0;
    for (vector<double>::size_type i=0; i<loga_int_lim.size()-1; i++) {
      checkForErr(gsl_integration_qag(&F, loga_int_lim[i], loga_int_lim[i+1],
				      epsabs, epsrel,
				      GSL_INTEG_GAUSS61, MAXINTERVAL,
				      w, &res, &err));
      result += res;
    }

    // Final interval, out to infinity
    checkForErr(gsl_integration_qagiu(&F, loga_int_lim.back(),
				      epsabs, epsrel,
				      MAXINTERVAL, w, &res, &err));
    result += res;

    // Take down integration machinery
    gsl_integration_workspace_free(w);
    
    // Return
    return -log(1.0 - fc0 + result);
  }
}

// Vectorized version of previous function
vector<double> pwind::tau_c(const vector<double>& u,
			    const double tXtw,
			    const double fj,
			    const double boltzfac,
			    const double fw,
			    const double varpi,
			    const double varpi_t,
			    const double a0,
			    const double a1,
			    const double epsabs,
			    const double epsrel) const {
  // Output holder
  vector<double> result(u.size());

  // Parallel loop over inputs
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (std::vector<double>::size_type i=0; i<u.size(); i++)
    result[i] = tau_c(u[i], tXtw, fj, boltzfac, fw, varpi, varpi_t,
		      a0, a1, epsabs, epsrel);

  // Return
  return result;
}


double pwind::tau_c(const double u,
		    const vector<double>& u_trans,
		    const vector<double>& tXtw,
		    const double fj,
		    const double boltzfac,
		    const double fw,
		    const double varpi,
		    const double varpi_t,
		    const double a0,
		    const double a1,
		    const double epsabs,
		    const double epsrel) const {

  // Get integration limits using tighter tolerance
  vector< vector<double> > alim(u_trans.size());
  vector<double>::size_type nlim_max = 0;
  for (vector<double>::size_type i=0; i<u_trans.size(); i++) {
    vector<double> alim_v = alimits(u - u_trans[i], varpi, varpi_t,
				    epsabs/1e4, epsrel/1e4);
    alim[i] = geom->a_lim(alim_v, varpi, varpi_t, u - u_trans[i]);
    if (nlim_max < alim[i].size()) nlim_max = alim[i].size();
  }
  // Bail out if all limits have zero size
  if (nlim_max == 0) return 0.0;
  
  // Apply input limits
  for (vector<double>::size_type j=0; j<alim.size(); j++) {
    for (vector<double>::size_type i=0; i<alim[j].size()/2; i++) {
      if (alim[j][2*i] < a0) alim[j][2*i] = a0;
      if (alim[j][2*i+1] > a1) alim[j][2*i+1] = a1;
    }
  }

  // Get global limits on a
  vector<double> alim_glob(2);
  alim_glob[0] = numeric_limits<double>::max();
  alim_glob[1] = 1.0;
  for (vector<double>::size_type i=0; i<alim.size(); i++) {
    if (alim[i].size() == 0) continue;
    if (alim_glob[0] > alim[i][0]) alim_glob[0] = alim[i][0];
    if (alim_glob[1] < alim[i].back()) alim_glob[1] = alim[i].back();
  }
    
  // Get tau_c
  double fc0;
  if (fw <= 0.0) fc0 = fc(alim_glob[0]);
  else fc0 = fw * y(alim_glob[0]) / SQR(alim_glob[0]);
  if (const_sa()) {
    // Constant solid angle case
    double tau = 0.0;
    for (vector<double>::size_type i=0; i<tXtw.size(); i++) {
      if (alim[i].size() == 0) continue;
      tau += tXtw[i] * fj * (1.0 - boltzfac) *
	Phi_c(u - u_trans[i], fw, varpi, varpi_t, alim[i][0],
	      alim[i].back(), epsabs, epsrel, alim[i]);
    }
    return -log(1.0 - fc0 + fc0*exp(-tau));
  } else {
    // Non-constant solid angle; need to integrate

    // Data for GSL
    struct Phi_params_vec par;
    par.u = u;
    par.u_trans = u_trans;
    par.fw = fw;
    par.vp2 = SQR(varpi) + SQR(varpi_t);
    par.tXtw = tXtw;
    par.fj = fj;
    par.boltzfac = boltzfac;
    par.alim_k = alim;
    par.epsabs = epsabs / 10.0; // Tighter tolerance for inner
				// integral
    par.epsrel = epsrel / 10.0;
    par.w = this;
    gsl_function F;
    F.function = &tau_c_integ_vec;
    F.params = &par;

    // Set up integration
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);

    // Integrate; this sometimes runs into trouble because the integrand
    // can have a-near step function discontinuity very close to the
    // end of the integration interval, associated with the optical
    // depth tau changing from < 1 to > 1. Integrators do not handle
    // this gracefully. The approach we take is to break the interval
    // up into intervals that bracket tau = 1.
    double b_lo, b_hi;
    if (alim_glob[0] == 1.0) b_lo = -100.0;
    else b_lo = log(alim_glob[0]-1.0);
    b_hi = log(numeric_limits<double>::max());
    vector<double> loga_int_lim;
    loga_int_lim.push_back(log(alim_glob[0]));
    vector<double> tau_break = { 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0 };
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = nullptr;
    gsl_function G;
    for (vector<double>::size_type i=0; i<tau_break.size(); i++) {
      
      // Search for this value of tau
      par.tau_target = tau_break[i];
      if (tau_c_func_vec(b_lo, &par) * tau_c_func_vec(b_hi, &par) > 0) {
	// tau never reaches this value, so stop
	break;
      }

      // Optical depth reaches target, so find out where
      if (s == nullptr) {
	s = gsl_root_fsolver_alloc(T);
	G.function = &tau_c_func_vec;
	G.params = &par;
      }
      gsl_root_fsolver_set(s, &G, b_lo, b_hi);
      double a_break;
      for (int iter=0, status=GSL_CONTINUE;
	   (iter<100) && (status=GSL_CONTINUE);
	   iter++) {
	status = gsl_root_fsolver_iterate(s);
	a_break = 1.0 + exp(gsl_root_fsolver_root(s));
	b_lo = gsl_root_fsolver_x_lower(s);
	b_hi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(b_lo, b_hi, 1e-2, 1e-2);
      }
      loga_int_lim.push_back(log(a_break));

      // Reset search interval
      b_lo = log(a_break-1.0);
      b_hi = log(numeric_limits<double>::max());
    }
    if (s != nullptr) gsl_root_fsolver_free(s);

    // Now integrate over broken-up interval
    double res, result, err;
    result = 0.0;
    for (vector<double>::size_type i=0; i<loga_int_lim.size()-1; i++) {
      checkForErr(gsl_integration_qag(&F, loga_int_lim[i], loga_int_lim[i+1],
				      epsabs, epsrel,
				      GSL_INTEG_GAUSS61, MAXINTERVAL,
				      w, &res, &err));
      result += res;
    }

    // Final interval, out to infinity
    checkForErr(gsl_integration_qagiu(&F, loga_int_lim.back(), epsabs,
				      epsrel, MAXINTERVAL, w,
				      &res, &err));
    result += res;

    // Take down integration machinery
    gsl_integration_workspace_free(w);
    
    // Return
    return -log(1.0 - fc0 + result);
  }
}

// Vectorized version of previous function
vector<double> pwind::tau_c(const vector<double>& u,
			    const vector<double>& u_trans,
			    const vector<double>& tXtw,
			    const double fj,
			    const double boltzfac,
			    const double fw,
			    const double varpi,
			    const double varpi_t,
			    const double a0,
			    const double a1,
			    const double epsabs,
			    const double epsrel) const {

  // Output holder
  vector<double> result(u.size());

  // Parallel loop over inputs
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (std::vector<double>::size_type i=0; i<u.size(); i++)
    result[i] = tau_c(u[i], u_trans, tXtw, fj, boltzfac, fw, varpi, varpi_t,
		      a0, a1, epsabs, epsrel);

  // Return
  return result;
}


// This routine produces an interpolating function that can be used to
// evaluate tau over the range in radius from a0 to a1;
// optical depths are computed starting from a0 unless reverse is
// true, in which case they are computed from a1. The accuracy of the
// interpolated exp(-tau) should be substantially higher than the
// overall error tolerance, so that it is safe to use the
// interpolating function in another numerical integration.
void pwind::tau_interp(const double u,
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
		       gsl_interp_accel **acc) const {

  // Set limits
  double loga0 = log(a0);
  double loga1 = log(a1);

  // Initialize a coarse grid
  vector<double> loga(gsl_interp_type_min_size(gsl_interp_cspline));
  for (vector<double>::size_type i=0; i<loga.size(); i++)
    loga[i] = loga0 + i * (loga1 - loga0) / ((double) loga.size()-1);

  // Compute e^-tau on the coarse grid; use higher accuracy in this
  // step than target final accuracy
  vector<double> Phi(loga.size(), 0.0), etau(loga.size(), 0.0);
  etau[0] = 1.0;
  for (vector<double>::size_type i=1; i<loga.size(); i++) {
    if (correlated)
      Phi[i] = Phi[i-1] +
	Phi_c(u, fw, varpi, varpi_t, exp(loga[i-1]), exp(loga[i]),
	      epsabs/10., epsrel/10.);
    else
      Phi[i] = Phi[i-1] +
	Phi_uc(u, varpi, varpi_t, exp(loga[i-1]), exp(loga[i]),
	       epsabs/10., epsrel/10.);
    etau[i] = exp(-tXtw * fj * (1.0-boltzfac) * Phi[i]);
  }

  // Now loop to improve accuracy
  while (1) {

    // Build interpolator for current grid
    gsl_interp *interp_lin = gsl_interp_alloc(gsl_interp_linear, loga.size());
    gsl_interp *interp_cube = gsl_interp_alloc(gsl_interp_cspline, loga.size());
    gsl_interp_accel *acc_lin = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_cube = gsl_interp_accel_alloc();
    gsl_interp_init(interp_lin, loga.data(), etau.data(), loga.size());
    gsl_interp_init(interp_cube, loga.data(), etau.data(), loga.size());

    // Go through grid and find all intervals that require subdivision
    // to hit accuracy target; note that we use 10x normal accuracy
    // for the integrations here, but sqrt(10)x normal accuracy for testing
    // the interpolation
    vector<vector<double>::size_type> intervals;
    for (vector<double>::size_type i=0; i<loga.size()-1; i++) {
      unsigned int ncheck = 4;
      for (unsigned int j=1; j<ncheck; j++) {
	double loga_pt = (ncheck-j) * loga[i]/((double) ncheck) +
	  j * loga[i+1]/((double) ncheck);
	double etau1 = gsl_interp_eval(interp_lin, loga.data(),
				       etau.data(), loga_pt, acc_lin);
	double etau2 = gsl_interp_eval(interp_cube, loga.data(),
				       etau.data(), loga_pt, acc_lin);
	double abserr = fabs(etau1-etau2);
	double relerr = abserr/(etau1 + numeric_limits<double>::epsilon());
	if (abserr > epsabs/sqrt(10.) && relerr > epsrel/sqrt(10.)) {
	  intervals.push_back(i);
	  break;
	}
      }
    }

    // Free interpolation functions
    gsl_interp_free(interp_lin);
    gsl_interp_free(interp_cube);
    gsl_interp_accel_free(acc_lin);
    gsl_interp_accel_free(acc_cube);

    // If we found no intervals that require further division, we're
    // done
    if (intervals.size() == 0) break;
    
    // If we're here, we need to subdivide some intervals by placing
    // points in the middle of them; solve for the new values of Phi
    // at these points
    for (vector<double>::size_type i=0; i<intervals.size(); i++) {
      double loga_tmp = 0.5*(loga[i+intervals[i]] +
			     loga[i+intervals[i]+1]);
      double Phi_tmp;
      if (correlated)
	Phi_tmp = Phi[i+intervals[i]] +
	  Phi_c(u, fw, varpi, varpi_t,
		exp(loga[i+intervals[i]]),
		exp(loga_tmp),
		epsabs/10.,
		epsrel/10.);
      else
	Phi_tmp = Phi[i+intervals[i]] +
	  Phi_uc(u, varpi, varpi_t,
		 exp(loga[i+intervals[i]]),
		 exp(loga_tmp),
		 epsabs/10.,
		 epsrel/10.);
      loga.insert(loga.begin()+i+intervals[i]+1, loga_tmp);
      Phi.insert(Phi.begin()+i+intervals[i]+1, Phi_tmp);
      etau.insert(etau.begin()+i+intervals[i]+1,
		  exp(-tXtw*fj*(1.0-boltzfac)*Phi_tmp));
    }
  }

  // Construct a spline interpolation to tau
  vector<double> tau(Phi.size());
  if (!reverse)
    for (vector<double>::size_type i=0; i<Phi.size(); i++)
      tau[i] = tXtw * fj * (1.0-boltzfac) * Phi[i];
  else
    for (vector<double>::size_type i=0; i<Phi.size(); i++)
      tau[i] = tXtw * fj * (1.0-boltzfac) * (Phi.back()-Phi[i]);    
  *spline = gsl_spline_alloc(gsl_interp_linear, loga.size());
  *acc = gsl_interp_accel_alloc();
  gsl_spline_init(*spline, loga.data(), tau.data(), loga.size());
}

static double tau_interp_var_fc_integ(double loga_pr, void *params) {
  // This routine returns the integrand of the integral that defines
  // tau_c for self-absorption
  struct tau_params *par = (struct tau_params *) params;
  
  // Get a; catch invalid values
  double a_pr = exp(loga_pr);
  if (a_pr > numeric_limits<double>::max()) return 0.0;

  // Get values
  double tau_a_pr = gsl_spline_eval(par->tau_spl, loga_pr, par->tau_acc);
  double tau_a = gsl_spline_eval(par->tau_spl, par->loga, par->tau_acc);
  double dfcda = fabs(par->w->dfcda(a_pr, par->fw))/par->fc0;
  
  // Get integrand
  double integ = a_pr * dfcda * exp(-(tau_a_pr - tau_a));
  return integ;
}

// This routine produces the same output as tau_interp, but for the
// case of correlated winds, a variable covering fraction, and
// emission at velocities u < 0; this case requires special handling
void pwind::tau_interp_var_fc(const double u,
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
			      gsl_interp_accel **acc) const {

  // As a first step, we build an interpolating function for tau
  // exactly as we would without the covering factor weighting, just
  // assuming full coverage. This needs very high error tolerance.
  gsl_spline *spl_cover;
  gsl_interp_accel *acc_cover;
  tau_interp(u, tXtw, fj, boltzfac, fw, varpi, varpi_t, a0, a1,
	     true, false,
	     epsabs/100., epsrel/100.,
	     &spl_cover, &acc_cover);

  // Prepare for integration
  struct tau_params par;
  par.fw = fw;
  par.tau_spl = spl_cover;
  par.tau_acc = acc_cover;
  par.w = this;
  gsl_function F;
  F.function = &tau_interp_var_fc_integ;
  F.params = &par;

  // Set up integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);

  // Now we need to build a grid of tau

  // Set limits
  double loga0 = log(a0);
  double loga1 = log(a1);

  // Initialize a coarse grid
  vector<double> loga(gsl_interp_type_min_size(gsl_interp_cspline));
  for (vector<double>::size_type i=0; i<loga.size(); i++)
    loga[i] = loga0 + i * (loga1 - loga0) / ((double) loga.size()-1);
  vector<double>::size_type ngood = 0;

  // Compute attenuation on a coarse grid; use higher than target
  // accuracy, since we will need to interpolate on these data to get
  // the final answer
  vector<double> tau(loga.size());
  for (vector<double>::size_type i=0; i<loga.size()-1; i++) {
    par.loga = loga[i];
    par.fc0 = fc(exp(par.loga), fw);
    if (par.fc0 == 0.0) tau[i] = nan("");
    else {

      // Figure out where the optical depth relative to the current
      // position goes through various critical values
      double tau_start = gsl_spline_eval(spl_cover, par.loga, acc_cover);
      vector<double> loga_int_lim = { par.loga };
      vector<double> tau_break = { 0.1, 0.316, 1.0, 3.16, 10.0 };
      size_t ptr = 1;
      for (vector<double>::size_type k=0; k<tau_break.size(); k++) {
	for ( ; ptr<spl_cover->size; ptr++) {
	  if (spl_cover->y[ptr] > tau_break[i]+tau_start) {
	    double loga_val = spl_cover->x[ptr-1] +
	      (1.0 - (spl_cover->y[ptr] - tau_break[i] - tau_start) /
	       (spl_cover->y[ptr] - spl_cover->y[ptr-1])) *
	      (spl_cover->x[ptr] - spl_cover->x[ptr-1]);
	    loga_int_lim.push_back(loga_val);
	    break;
	  }
	}
      }
      loga_int_lim.push_back(numeric_limits<double>::max());

      // Now integrate
      double integ = 0.0;
      for (vector<double>::size_type j=0; j<loga_int_lim.size()-1; j++) {
	double res, err;
	double loga_lo = max(loga[i], loga_int_lim[j]);
	double loga_hi = min(loga1, loga_int_lim[j+1]);
	if (loga_lo >= loga_hi) continue;
	if (loga_hi < log(numeric_limits<double>::max()))
	  checkForErr(gsl_integration_qag(&F, loga_lo, loga_hi, epsabs/10.,
					  epsrel/10., MAXINTERVAL,
					  GSL_INTEG_GAUSS61, w, &res, &err));
	else
	  checkForErr(gsl_integration_qagiu(&F, loga_lo, epsabs/10.,
					    epsrel/10., MAXINTERVAL, w,
					    &res, &err));
	integ += res;
      }
      if (integ == 0) {
	// Catch values where we can't evaluate tau sensibly
	tau[i] = nan("");
      } else {
	tau[i] = -log(integ);
	ngood++;
      }
    }
  }

  // Now loop to improve accuracy
  while (1) {

    // Figure out where we need to refine
    vector<vector<double>::size_type> intervals;
    if (ngood >= gsl_interp_type_min_size(gsl_interp_cspline)) {
	
      // Build interpolator for current grid
      gsl_interp *interp_lin = gsl_interp_alloc(gsl_interp_linear, ngood);
      gsl_interp *interp_cube = gsl_interp_alloc(gsl_interp_cspline, ngood);
      gsl_interp_accel *acc_lin = gsl_interp_accel_alloc();
      gsl_interp_accel *acc_cube = gsl_interp_accel_alloc();
      gsl_interp_init(interp_lin, loga.data(), tau.data(), ngood);
      gsl_interp_init(interp_cube, loga.data(), tau.data(), ngood);

      // Go through grid and find all intervals that require subdivision
      // to hit accuracy target; note that we use 10x normal accuracy
      // for the integrations here, but sqrt(10)x normal accuracy for testing
      // the interpolation
      for (vector<double>::size_type i=0; i<ngood-1; i++) {
	if (loga[i+1] - loga[i] <= epsabs/10.) continue;
	if (loga[i+1]/loga[i] - 1.0 <= epsrel/10.) continue;
	unsigned int ncheck = 8;
	for (unsigned int j=1; j<ncheck; j++) {
	  double loga_pt = (ncheck-j) * loga[i]/((double) ncheck) +
	    j * loga[i+1]/((double) ncheck);
	  double tau1 = gsl_interp_eval(interp_lin, loga.data(),
					tau.data(), loga_pt, acc_lin);
	  double tau2 = gsl_interp_eval(interp_cube, loga.data(),
					tau.data(), loga_pt, acc_lin);
	  double abserr = fabs(exp(-tau1)-exp(-tau2));
	  double relerr = abserr /
	    (exp(-tau1) + numeric_limits<double>::epsilon());
	  if (abserr > epsabs && relerr > epsrel) {
	    intervals.push_back(i);
	    break;
	  }
	}
      }

      // Refine final interval between valid point and NaN if needed
      if (ngood < loga.size()) {
	if (exp(-tau[ngood-1]) > epsabs &&
	    loga[ngood] - loga[ngood-1] > epsabs &&
	    loga[ngood]/loga[ngood-1] - 1.0 > epsrel)
	  intervals.push_back(ngood-1);
      }
      
      // Free interpolation functions
      gsl_interp_free(interp_lin);
      gsl_interp_free(interp_cube);
      gsl_interp_accel_free(acc_lin);
      gsl_interp_accel_free(acc_cube);
      
    } else {
      // Handle special case where number of points where we
      // successfully evaluated the integral is fewer than required for
      // cubic interpolation; if this happens, we automatically flag all
      // intervals between good points for refinement; if there are no
      // good intervals, we flag the first one for refinement
      for (vector<double>::size_type i=0; i<ngood; i++)
	intervals.push_back(i);

    }  

    // If we found no intervals that require further division, we're
    // done
    if (intervals.size() == 0) break;
    
    // If we're here, we need to subdivide some intervals by placing
    // points in the middle of them; solve for the new values of tau
    // at these points
    for (vector<double>::size_type i=0; i<intervals.size(); i++) {
      par.loga = 0.5*(loga[i+intervals[i]] +
		      loga[i+intervals[i]+1]);
      par.fc0 = fc(exp(par.loga), fw);
      // Catch values where we overflow
      double tau_tmp;
      if (par.fc0 == 0.0) tau_tmp = nan("");
      else {

	// Figure out where the optical depth relative to the current
	// position goes through various critical values
	double tau_start = gsl_spline_eval(spl_cover, par.loga, acc_cover);
	vector<double> loga_int_lim = { par.loga };
	vector<double> tau_break = { 0.1, 0.316, 1.0, 3.16, 10.0, 31.6 };
	size_t ptr = 1;
	for (vector<double>::size_type i=0; i<tau_break.size(); i++) {
	  for ( ; ptr<spl_cover->size; ptr++) {
	    if (spl_cover->y[ptr] > tau_break[i]+tau_start) {
	      double loga_val = spl_cover->x[ptr-1] +
		(1.0 - (spl_cover->y[ptr] - tau_break[i] - tau_start) /
		 (spl_cover->y[ptr] - spl_cover->y[ptr-1])) *
		(spl_cover->x[ptr] - spl_cover->x[ptr-1]);
	      loga_int_lim.push_back(loga_val);
	      break;
	    }
	  }
	}
	loga_int_lim.push_back(numeric_limits<double>::max());

	// Now integrate
	double integ = 0.0;
	for (vector<double>::size_type j=0; j<loga_int_lim.size()-1; j++) {
	  double res, err;
	  double loga_lo = max(par.loga, loga_int_lim[j]);
	  double loga_hi = min(loga1, loga_int_lim[j+1]);
	  if (loga_lo >= loga_hi) continue;
	  if (loga_hi < log(numeric_limits<double>::max())) {
	    checkForErr(gsl_integration_qag(&F, loga_lo, loga_hi, epsabs/10.,
					    epsrel/10., MAXINTERVAL,
					    GSL_INTEG_GAUSS61, w, &res, &err));
	    integ += res;
	  } else {
	    checkForErr(gsl_integration_qagiu(&F, loga_lo, epsabs/10.,
					      epsrel/10., MAXINTERVAL, w,
					      &res, &err));
	    integ += res;
	  }
	}

	// Store result
	if (integ > 0) {
	  tau_tmp = -log(integ);
	  ngood++;
	} else {
	  tau_tmp = nan("");
	}
      }
      loga.insert(loga.begin()+i+intervals[i]+1, par.loga);
      tau.insert(tau.begin()+i+intervals[i]+1, tau_tmp);
    }
  }

  // Replace the NaN's with copies of the last valid value
  for (vector<double>::size_type i=ngood; i<loga.size(); i++)
    tau[i] = tau[ngood-1];

  // Construct a spline interpolation to the final data set  
  *spline = gsl_spline_alloc(gsl_interp_linear, loga.size());
  *acc = gsl_interp_accel_alloc();
  gsl_spline_init(*spline, loga.data(), tau.data(), loga.size());

  // Free memory
  gsl_integration_workspace_free(w);
  gsl_spline_free(spl_cover);
  gsl_interp_accel_free(acc_cover);
}

////////////////////////////////////////////////////////////////////////
// Subcritical, thin emission calculations
////////////////////////////////////////////////////////////////////////

static double Xi_integ(double loga, void *params) {
  struct Xi_params *par = (struct Xi_params *) params;

  // Get a; catch invalid values
  double a = exp(loga);
  if (a > numeric_limits<double>::max()) return 0.0;
  if (SQR(a) < par->vp2) return 0.0;

  // Get x and dU2/dx
  double ur = fabs(par->u) / sqrt(1.0 - par->vp2/SQR(a));
  if (ur > par->w->uMax()) return 0.0;
  double x = par->w->X(ur, a);
  if (x > par->w->xcr() + log(par->w->getFcrit())) return 0.0;
  double dU2dx = par->w->dU2dx(x, a);
  
  // Compute integrand
  double pMval = pM(x, par->w->s());
  double pAval = pA(x, par->w->s());
  double ret = 1.0 / (par->w->y(a) * (SQR(a) - par->vp2)) *
    SQR(pMval) / (pAval * ur * fabs(dU2dx));
  return ret;
}

double pwind::Xi(const double u,
		 const double varpi,
		 const double varpi_t,
		 const double epsabs,
		 const double epsrel) const {

  // Check if velocity is outside bounds
  if (u == 0) return numeric_limits<double>::max();
  if (fabs(u) >= umax) return 0.0;

  // Get integration limits; use very high accuracy for this step
  vector<double> alim_v = alimits(u, varpi, varpi_t,
				  epsabs/1e4, epsrel/1e4);
  vector<double> alim = geom->a_lim(alim_v, varpi, varpi_t, u);

  // Load the GSL data
  struct Xi_params par;
  par.u = fabs(u);
  par.vp2 = SQR(varpi) + SQR(varpi_t);
  par.w = this;
  gsl_function F;
  F.function = &Xi_integ;
  F.params = &par;

  // Set up integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);

  // Integrate
  double integ = 0.0;
  for (vector<double>::size_type i=0; i<alim.size()/2; i++) {
    double result, err;
    if (alim[2*i+1] < numeric_limits<double>::max()) {
      checkForErr(gsl_integration_qag(&F, log(alim[2*i]), log(alim[2*i+1]),
				      epsabs, epsrel, MAXINTERVAL,
				      GSL_INTEG_GAUSS61, w, &result, &err));
    } else {
      checkForErr(gsl_integration_qagiu(&F, log(alim[2*i]), epsabs, epsrel,
					MAXINTERVAL, w, &result, &err));
    }
    integ += result;
  }
  integ /= SQR(zeta_M);

  // Free and return
  gsl_integration_workspace_free(w);
  return integ;
}

// Vectorized version of previous function
vector<double> pwind::Xi(const vector<double>& u,
			 const double varpi,
			 const double varpi_t,
			 const double epsabs,
			 const double epsrel) const {
  // Output holder
  vector<double> result(u.size());

  // Parallel loop over inputs
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (std::vector<double>::size_type i=0; i<u.size(); i++)
    result[i] = Xi(u[i], varpi, varpi_t, epsabs, epsrel);

  // Return
  return result;
}


static double xi_integ_x(double x, void *params) {
  struct xi_params *par = (struct xi_params *) params;

  // Get U2
  double U2 = par->w->U2(x, par->a);

  // Return
  double ret = SQR(pM(x, par->w->s())) /
    (pA(x, par->w->s()) * U2);
  return ret;
}

static double xi_integ_s(double s, void *params) {
  struct xi_params *par = (struct xi_params *) params;

  // Get a; catch invalid values
  par->a = sqrt(SQR(s) + par->vp2);

  // Get x range
  vector<double> xlim = par->w->xlimits(par->a);
  if (xlim.size() == 0) return 0.0;

  // Prepare GSL integration
  gsl_function F;
  F.function = &xi_integ_x;
  F.params = par;
  gsl_integration_workspace *w =
    gsl_integration_workspace_alloc(MAXINTERVAL);

  // Integrate
  double result, err;
  int success;
  if (xlim[0] == -numeric_limits<double>::max()) {
    success = gsl_integration_qagil(&F, xlim[1]+log(par->w->getFcrit()),
				    par->epsabs, par->epsrel,
				    MAXINTERVAL, w, &result, &err);
  } else {
    if (xlim[0] < xlim[1]+log(par->w->getFcrit()))
      success = gsl_integration_qag(&F, xlim[0],
				    xlim[1]+log(par->w->getFcrit()),
				    par->epsabs, par->epsrel,
				    MAXINTERVAL, GSL_INTEG_GAUSS61, w,
				    &result, &err);
    else {
      result = 0.0;
      success = GSL_SUCCESS;
    }
  }

  // Store error status
  if (success != GSL_SUCCESS) par->gsl_err = success;

  // Multiply by prefactor
  result /= (par->w->y(par->a) * SQR(par->a));

  // Free memory and return
  gsl_integration_workspace_free(w);
  return result;
}

double pwind::xi(const double varpi,
		 const double varpi_t,
		 const double epsabs,
		 const double epsrel) const {
  
  // Get integration limits
  vector<double> slim = geom->s_crit(varpi, varpi_t, 0.0);
  if (slim.size() == 0) return 0.0;

  // Catch if we've been asked to integrate over a radius < 1; if so,
  // the result is undefined, so return infinity
  if (SQR(varpi) + SQR(varpi_t) <= 1.0)
    return numeric_limits<double>::max();

  // Load the GSL data; set tolerance higher for inner integral than
  // final tolerance
  struct xi_params par;
  par.vp2 = SQR(varpi) + SQR(varpi_t);
  par.epsrel = epsrel/10.;
  par.epsabs = epsabs/10.;
  par.w = this;
  par.gsl_err = GSL_SUCCESS;
  gsl_function F;
  F.function = &xi_integ_s;
  F.params = &par;
  
  // Allocate workspace
  gsl_integration_workspace *w =
    gsl_integration_workspace_alloc(MAXINTERVAL);

  // Integrate
  double integ = 0.0;
  for (vector<double>::size_type i=0; i<slim.size()/2; i++) {
    double result, err;
    if (slim[2*i+1] > sqrt(SQR(amax_abs)+par.vp2))
      slim[2*i+1] = sqrt(SQR(amax_abs)+par.vp2);
    if (slim[2*i] < -sqrt(SQR(amax_abs)+par.vp2))
      slim[2*i] = -sqrt(SQR(amax_abs)+par.vp2);
    if (slim[2*i] > -numeric_limits<double>::max() &&
	slim[2*i+1] < numeric_limits<double>::max()) {
      checkForErr(gsl_integration_qag(&F, slim[2*i], slim[2*i+1],
				      epsabs, epsrel, MAXINTERVAL,
				      GSL_INTEG_GAUSS61, w, &result, &err));
    } else if (slim[2*i] > -numeric_limits<double>::max()) {
      checkForErr(gsl_integration_qagiu(&F, slim[2*i], epsabs, epsrel,
					MAXINTERVAL, w, &result, &err));
    } else if (slim[2*i+1] < numeric_limits<double>::max()) {
      checkForErr(gsl_integration_qagil(&F, slim[2*i+1], epsabs, epsrel,
					MAXINTERVAL, w, &result, &err));
    } else {
      checkForErr(gsl_integration_qagi(&F, epsabs, epsrel,
				       MAXINTERVAL, w, &result, &err));
    }
    integ += result;
  }
  integ /= SQR(zeta_M);

  // If one of the inner integrals threw an error, be sure to record
  // that
  if (par.gsl_err != GSL_SUCCESS) gsl_err_stat = par.gsl_err;

  // Free and return
  gsl_integration_workspace_free(w);
  return integ;
}


////////////////////////////////////////////////////////////////////////
// LTE emission calculations
////////////////////////////////////////////////////////////////////////

static double eta_integ(double loga, void *params) {
  struct eta_params *par = (struct eta_params *) params;

  // Get a; catch invalid values
  double a = exp(loga);
  if (a > numeric_limits<double>::max()) return 0.0;
  if (SQR(a) < par->vp2) return 0.0;
  
  // Get x and dU2/dx
  double ur = fabs(par->u) / sqrt(1.0 - par->vp2/SQR(a));
  if (ur > par->w->uMax()) return 0.0;
  double x = par->w->X(ur, a);
  if (x > par->w->xcr() + log(par->w->getFcrit())) return 0.0;
  double dU2dx = par->w->dU2dx(x, a);

  // Get beta
  double beta;
  if (par->thin) beta = 1.0;
  else beta = exp(-gsl_spline_eval(par->tau_spl, loga, par->tau_acc));

  // Return integrand; note that we add a numerical safety to dU2dx to
  // handle a specific case: if both X(ur, a) and the limits on the
  // integration have to be determined numerically, then when one is
  // at the limit of the integral, it is possible that the inversion
  // required to get X(ur, a) will fail, and X will wind up being set
  // to xlo. When this happens, evaluation of pM(x) / dU2dx can
  // produce 0 / 0 = NaN. We avoid this by adding a safety in the
  // denominator of the intregrand.
  double ret = 1.0 / (a - par->vp2/a) * beta *
    pM(x, par->w->s()) /
    (fabs(dU2dx) + numeric_limits<double>::epsilon());
  return ret;
}


double pwind::eta(const double u,
		  const double tXtw,
		  const double fj,
		  const double boltzfac,
		  const bool correlated,
		  const double fw,
		  const double varpi,
		  const double varpi_t,
		  const bool thin,
		  const double epsabs,
		  const double epsrel) const {

  // Check if velocity is outside bounds
  if (u == 0) return numeric_limits<double>::max();
  if (fabs(u) >= umax) return 0.0;

  // Get integration limits; note that we tighten the tolerance
  // considerably before this step, because getting the limits off can
  // have unexpectedly large consequences
  vector<double> alim_v = alimits(u, varpi, varpi_t,
				  epsabs/1e4, epsrel/1e4);
  vector<double> alim = geom->a_lim(alim_v, varpi, varpi_t, u);
  if (alim.size() == 0) return 0.0;
  double a0 = numeric_limits<double>::max(), a1 = 0.0;
  for (vector<double>::size_type i=0; i<alim.size(); i++) {
    if (alim[i] < a0) a0 = alim[i];
    if (alim[i] > a1) a1 = alim[i];
  }

  // Load the GSL data
  struct eta_params par;
  par.u = u;
  par.vp2 = SQR(varpi) + SQR(varpi_t);
  par.thin = thin;
  par.loga_max = log(a1);
  par.w = this;
  gsl_function F;
  F.function = &eta_integ;
  F.params = &par;

  // Build an interpolation function that gives tau as a
  // function of a
  if (!thin) {
    if (!correlated || const_sa() || u > 0) {
      tau_interp(u, tXtw, fj, boltzfac, fw, varpi, varpi_t, a0,
		 a1, correlated, u<0, epsabs/10., epsrel/10.,
		 &par.tau_spl, &par.tau_acc);
    } else {
      tau_interp_var_fc(u, tXtw, fj, boltzfac, fw, varpi, varpi_t,
			a0, a1, epsabs/10., epsrel/10., &par.tau_spl,
			&par.tau_acc);
    }
  } else {
    par.tau_spl = NULL;
    par.tau_acc = NULL;
  }
  
  // Allocate workspace
  gsl_integration_workspace *w =
    gsl_integration_workspace_alloc(MAXINTERVAL);

  // Integration needs a bit of care, because there's an e^(-tau)
  // factor inside our integral that can cause a near-step function
  // discontinuity that causes havoc in the integrator. We therefore
  // want to break the integral up into different segments based on
  // the value of tau, then add them. Here we find the list of
  // critical points at which we want to break things up.
  vector<double> loga_int_lim = { 0.0 };
  if (!thin) {
    vector<double> tau_break = { 0.1, 0.316, 1.0, 3.16, 10.0 };
    // Find points in the table where tau crosses a critical point
    for (size_t i=0; i<par.tau_spl->size-1; i++) {
      for (vector<double>::size_type j=0; j<tau_break.size(); j++) {
	if ((par.tau_spl->y[i]-tau_break[j]) *
	    (par.tau_spl->y[i+1]-tau_break[j]) < 0.0) {
	  loga_int_lim.push_back(par.tau_spl->x[i+1]);
	  break;
	}
      }
    }
  }
  loga_int_lim.push_back(numeric_limits<double>::max());

  // Now itegrate; we have a set of geometric intervals and a set of
  // critical points. Our approach is that we loop over the geometric
  // intervals, and within each one we calculate an integral by
  // breaking up the interval over the critical points.
  double integ = 0.0;
  for (vector<double>::size_type i=0; i<alim.size()/2; i++) {
    double res, err;
    for (vector<double>::size_type j=0; j<loga_int_lim.size()-1; j++) {
      double loga_lo = max(loga_int_lim[j], log(alim[2*i]));
      double loga_hi = min(loga_int_lim[j+1], log(alim[2*i+1]));
      if (loga_lo >= loga_hi) continue;
      if (loga_hi < log(numeric_limits<double>::max())) {
	checkForErr(gsl_integration_qag(&F, loga_lo, loga_hi, epsabs,
					epsrel, MAXINTERVAL,
					GSL_INTEG_GAUSS61, w,
					&res, &err));
      } else {
	if ((fabs(2*loga_lo - log(par.vp2)) > 0.1) &&
	    (loga_lo > 0.1)) {
	  gsl_integration_qagiu(&F, loga_lo, epsabs, epsrel,
				MAXINTERVAL, w, &res, &err);
	} else {
	  // The integral out to infinity can have issues if there is a
	  // singularity near the lower limit because a_min = varpi or
	  // a_min = 1. We check for this case and avoid it by
	  // breaking up the integral.
	  double res1; 
	  checkForErr(gsl_integration_qag(&F, loga_lo, loga_lo+0.1, epsabs,
					  epsrel, MAXINTERVAL,
					  GSL_INTEG_GAUSS61, w,
					  &res1, &err));
	  checkForErr(gsl_integration_qagiu(&F, loga_lo+0.1, epsabs, epsrel,
					    MAXINTERVAL, w, &res, &err));
	  res += res1;
	}
      }
      integ += res;
    }
  }  
  integ /= zeta_M;
  
  // Free and return
  gsl_integration_workspace_free(w);
  gsl_spline_free(par.tau_spl);
  gsl_interp_accel_free(par.tau_acc);
  return integ;
}

// Vectorized version of previous function
vector<double> pwind::eta(const vector<double>& u,
			  const double tXtw,
			  const double fj,
			  const double boltzfac,
			  const bool correlated,
			  const double fw,
			  const double varpi,
			  const double varpi_t, 
			  const bool thin,
			  const double epsabs,
			  const double epsrel) const {
  // Output holder
  vector<double> result(u.size());

  // Parallel loop over inputs
#ifdef _OPENMP
#  pragma omp parallel for schedule(dynamic, 4)
#endif
  for (std::vector<double>::size_type i=0; i<u.size(); i++)
    result[i] = eta(u[i], tXtw, fj, boltzfac, correlated, fw, varpi,
		    varpi_t, thin, epsabs, epsrel);

  // Return
  return result;
}


static double Psi_integ(double u, void *params) {
  struct Psi_params *par = (struct Psi_params *) params;
  double eta = par->w->eta(u, par->tXtw, par->fj, par->boltzfac,
			   par->correlated, par->fw, par->varpi,
			   par->varpi_t, par->thin, par->epsabs,
			   par->epsrel);
  return eta;
}


double pwind::Psi(const double tXtw,
		  const double fj,
		  const double boltzfac,
		  const bool correlated,
		  const double fw,
		  const double varpi,
		  const double varpi_t, 
		  const bool thin,
		  const double epsabs,
		  const double epsrel) const {

  // Psi diverges at varpi = 1 exactly
  if (SQR(varpi) + SQR(varpi_t) == 1.0)
    return numeric_limits<double>::max();

  // Load the GSL data
  struct Psi_params par;
  par.tXtw = tXtw;
  par.fj = fj;
  par.boltzfac = boltzfac;
  par.correlated = correlated;
  par.fw = fw;
  par.varpi = varpi;
  par.varpi_t = varpi_t;
  par.epsabs = epsabs/10.0;   // Tighter tolerance on inner integral
  par.epsrel = epsrel/10.0;
  par.thin = thin;
  par.w = this;
  gsl_function F;
  F.function = &Psi_integ;
  F.params = &par;

  // Set up integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(MAXINTERVAL);

  // Integrate
  double res, err;
  checkForErr(gsl_integration_qagi(&F, epsabs, epsrel, MAXINTERVAL,
				   w, &res, &err));

  // Take down integration machinery
  gsl_integration_workspace_free(w);

  // Return
  return res;
}
  
