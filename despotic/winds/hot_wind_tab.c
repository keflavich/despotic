// This little bit of code generates a table of solutions for the wind
// acceleration law for hot winds

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#define SQR(a) ((a)*(a))

// Structure used to pass data to the ODE derivatives function
typedef struct ode_params {
  double gex, uh;
  int yidx, midx;
} ode_params;

typedef struct root_params {
  double u_target, a;
  gsl_odeiv2_system *sys;
  gsl_odeiv2_driver *drv;
} root_params;

// Global variables
double eps_ode = 1.0e-10;  // accuracy goal for ODE integrations
double eps_solve = 1.0e-6; // accuracy goal of implicit equation solutions
bool verbose = false;      // verbosity level


/////////////////////////////////////////////////////////////////////////////
// Derivatives functions
/////////////////////////////////////////////////////////////////////////////

// Function to return dU/da =
// 1/sqrt(2 u a^2) [Gamma exp(-x) y (1 - u/uh)^2 - m]
double duda(double a, double u, double gex, double uh, double y, double m) {
  return 1.0 / (2.0*SQR(a)*u) *
    (gex * y * SQR(1.0 - u/uh) - m);
}

// Wrapper for the previous function, in the format required by GSL;
// also handled bounds checking
int duda_gsl(double a, const double u[], double uprime[],
	     void *params) {

  // Extact parameters
  ode_params *par = (ode_params *) params;
  double gex = par->gex;
  double uh = par->uh;
  double yidx = par->yidx;
  double midx = par->midx;

  // We need to be a little careful about what happens near u = 0,
  // because that is singular. In reality we should only ever have u >
  // 0, but the ODE driver will sometimes give us u < 0, and we need
  // to handle this properly. There are two distinct cases, which we
  // have to separate: the first is that we're near a = 1, u is small,
  // and the ODE integrator is exploring different values of u, not
  // knowing that there is a singularity at u = 0. In this case, we
  // should let integration continue, and just return du/da evaluated
  // as if u were positive instead of negative. The other case is that
  // we're approaching u = 0 from above because we're at a > 1, and
  // the flow has turned around. In this case we want to signal to
  // halt integration. The way we differentiate between these two
  // possibilities is based on the sign of du/da: in the former case
  // it is positive (under the replacement u -> |u|), while in the
  // latter case it is negative.
  if (u[0] < uh) {
    double y = pow(a, yidx);
    double m = pow(a, midx);
    uprime[0] = duda(a, u[0], gex, uh, y, m);
  } else {
    uprime[0] = 0.0;
  }
  if (u[0] < 0.0 && uprime[0] <= 0.0)
    return GSL_EBADFUNC;
  else
    return GSL_SUCCESS;
}


/////////////////////////////////////////////////////////////////////////////
// Utility function that computes the velocity near a = 1 using a
// series expansion. In this function, eps0 is the distance from a = 1
// at which the velocity is computed, i.e., it is computed at a = 1 + eps0
/////////////////////////////////////////////////////////////////////////////

double u0(double eps0, double gex, double uh,
	  double dyda, double dmda) {
  if (gex == 1.0) {
    return eps0 / (2.0*uh) * (-1.0 + sqrt(1.0 + 2.0*SQR(uh) * (dyda-dmda)));
  } else {
    return sqrt((gex-1.0) * eps0);
  }
}


/////////////////////////////////////////////////////////////////////////////
// Function to integrate the ODE describing the wind acceleration law
// to some specified radius a_stop, then return the velocity at that
// radius; if the wind reaches zero velocity before reaching a, the
// code returns 0. This routine assumes that the GSL setup has already
// been done, and takes pointers to the GSL system and driver as
// arguments. The value ustart specifies how small the initial value
// of u needs to be
/////////////////////////////////////////////////////////////////////////////

double ua(double a_stop,
	  double ustart,
	  gsl_odeiv2_system *sys,
	  gsl_odeiv2_driver *drv) {

  // Extract the quantities we need to initialize the integration
  ode_params *par = (ode_params *) sys->params;
  double gex = par->gex;
  double uh = par->uh;
  double dyda = par->yidx;
  double dmda = par->midx;

  // Set starting point for integration
  double a, u, eps0;
  eps0 = 1.0e-3;  // Starting guess
  while (1) {
    u = u0(eps0, gex, uh, dyda, dmda);
    if (u <= ustart/10.0 && eps0 < 0.1*(a_stop-1)) break;
    else eps0 /= 10.0;
  }
  a = 1.0 + eps0;
  
  // Integrate
  gsl_odeiv2_driver_reset_hstart(drv, eps0);
  int status = gsl_odeiv2_driver_apply(drv, &a, a_stop, &u);

  // If we reached target, return the value of u; otherwise return -1
  // as a flag value
  if (status == GSL_SUCCESS)
    return u;
  else
    return -1.0;
}


/////////////////////////////////////////////////////////////////////////////
// Function that takes as input a radius a and velocity u and finds
// the value of Gamma e^-x that produces the target velocity at that
// radius. This routine assumes that the GSL setup has already
// been done, and takes pointers to the GSL system and driver as
// arguments.
/////////////////////////////////////////////////////////////////////////////

// Helper function that returns residual
double gex_solve_residual(double gex, void *params) {
  root_params *par = (root_params *) params;
  ode_params *p = par->sys->params;
  p->gex = gex;
  double resid =
    par->u_target - ua(par->a, par->u_target, par->sys, par->drv);
  return resid;
}

double gex_solve(double a, double u,
		 double eps,
		 gsl_odeiv2_system *sys,
		 gsl_odeiv2_driver *drv,
		 gsl_root_fsolver *solver) {

  // If given a velocity of 0 and yidx >= midx, automatically return
  // gex = 1; for yidx < midx, set u to a small but non-zero value to
  // avoid numerical problems
  ode_params *p = (ode_params *) sys->params;
  if (u == 0) {
    if (p->yidx >= p->midx) return 1.0;
    else u = eps/10.0;
  }

  // Set parameters to pass to root solver
  root_params par = { u, a, sys, drv };

  // Set up GSL function to pass to solver
  gsl_function F = { &gex_solve_residual, &par };

  // Set lower bound on gex by using the fact that the acceleration
  // for the hot wind case is always slower than the corresponding
  // ideal case, so the gex we need for the hot gas case will always
  // be larger.
  double gex_min = 1.0;
  if (p->midx == 0) {
    if (p->yidx == 0)
      gex_min = SQR(u) * a/(a-1.0) + 1.0;
    else if (p->yidx == 1)
      gex_min = (SQR(u) + a/(a-1.0)) / log(a);
    else if (p->yidx == 2)
      gex_min = SQR(u) / (a-1.0) + 1.0/a;
  } else if (p->midx == 1) {
    if (p->yidx == 0)
      gex_min = (SQR(u) + log(a)) * a/(a-1.0);
    else if (p->yidx == 1)
      gex_min = SQR(u) / log(a) + 1.0;
    else if (p->yidx == 2)
      gex_min = (SQR(u) + log(a)) / (a-1.0);
  }
  if (gex_min < 1.0) gex_min = 1.0;

  // Evaluate velocity using gex_min; if it exceeds target, there is
  // no solution, return exactly 1 as a flag value
  p->gex = gex_min;
  if (ua(a, u, sys, drv) > u) return 1.0;

  // For the upper limit, we need to bracket manually; we start with
  // an initial guess for the upper limit on gex, integrate the
  // acceleration law, and then increase our guess if u is still too
  // small
  double gex_max = 100.0;
  while (1) {
    p->gex = gex_max;
    double u_gex_max = ua(a, u, sys, drv);
    if (u_gex_max > u) break;
    else gex_max *= 2.0;
  }
  
  // Initialize solver
  gsl_root_fsolver_set(solver, &F, gex_min, gex_max);

  // Solve
  int status = GSL_SUCCESS;
  double gex_lo, gex_hi;
  while (status == GSL_SUCCESS) {
    status = gsl_root_fsolver_iterate(solver);
    gex_lo = gsl_root_fsolver_x_lower(solver);
    gex_hi = gsl_root_fsolver_x_upper(solver);
    if (gsl_root_test_interval(gex_lo, gex_hi, eps, eps) == GSL_SUCCESS)
      break;
  }

  // Return solution if found; otherwise return 1
  if (status == GSL_SUCCESS)
    return 0.5*(gex_lo + gex_hi);
  else
    return 1.0;
}

/////////////////////////////////////////////////////////////////////////////
// Function to find the critical points along an acceleration curve
/////////////////////////////////////////////////////////////////////////////

void compute_critical_points(double qmax, double eps, ode_params par,
			     double *umax, double *q_umax,
			     double *q_stop,
			     gsl_odeiv2_driver *drv) {

  // Handle the special case where gex = 1 and yidx <= midx
  if (par.gex == 1.0 && par.yidx <= par.midx) {
    *umax = 0.0;
    *q_umax = *q_stop = -DBL_MAX;
    return;
  }

  // Reset ODE solver
  gsl_odeiv2_driver_reset(drv);
  
  // Set starting point for integration
  double dyda = par.yidx;
  double dmda = par.midx;
  double a0 = 1.0 + eps;
  double u0;
  if (par.gex == 1.0) {
    u0 = eps / (2.0*par.uh) * (-1.0 +
			       sqrt(1.0 + 2.0*SQR(par.uh) * (dyda-dmda)));
  } else {
    u0 = sqrt((par.gex-1.0) * eps);
  }
  
  // Integrate to find maximum velocity and stopping point. These occur in
  // different locations and for different reasons depending on
  // whether yidx < midx or yidx >= midx. For yidx < midx, the maximum
  // velocity occurs at finite a, and u -> 0 at finite a as well; we
  // therefore look for the maximum velocity within the interval, and
  // stopping point as the point where u -> 0. For yidx >= midx, the
  // wind asymptotes to a finite u as a -> inifinty, so we define the
  // stopping point as the point at which u stops increasing
  // significantly (relative to our accuracy goal). We therefore have
  // two different procedures here, depending on the relative values
  // of yidx and midx
  double amax = 1.0 + pow(10., qmax);
  if (par.yidx < par.midx) {

    // Case where velocity turns around; first integrate region where
    // velocity is increasing
    double a = a0;             // Starting position
    double u = u0;             // Starting velocity
    double da = 1.0e-3;        // Starting step size
    bool bracketed = false;    // Flag for if we have bracketed maximum
    bool reached_amax = false; // Flag for if we hit ceiling
    while (!reached_amax) {

      // Store previous values so we can back up when needed
      double a_last = a;
      double u_last = u;

      // Limit time step
      double a_next = a + da;
      if (a_next > amax) {
	a_next = amax;
	reached_amax = true;
      }

      // Take a step
      int status = gsl_odeiv2_driver_apply(drv, &a, a_next, &u);

      // If step was unsuccesful, reduce time step and try again
      if (status != GSL_SUCCESS) {
	a = a_last;
	u = u_last;
	da /= 2.0;
	reached_amax = false;
	gsl_odeiv2_driver_reset(drv);
	continue;
      }

      // Get derivative of velocity at current position
      double dlogu_dloga = (a/u) *
	duda(a, u, par.gex, par.uh, pow(a, par.yidx), pow(a, par.midx));

      // Action depends on value of dlogu_dloga
      if (fabs(dlogu_dloga) < eps ||
	  fabs(a - a_last) < DBL_MIN) {

	// Value is zero to within requested tolerance, so store
	// result and stop iterating
	*q_umax = log10(a-1);
	*umax = u;
	break;

      } else if (dlogu_dloga > 0.0) {

	// Derivative is positive, so wind is still accelerating;
	// increase step size if we have not yet bracketed the maximum
	if (!bracketed) da *= 2.0;

      } else {

	// Derivative is negative; flag that we have now bracketed the
	// root, back up to previous position, and decrease step size
	bracketed = true;
	a = a_last;
	u = u_last;
	reached_amax = false;
	da /= 2.0;

      }

    } // end integration loop

    // We have now found the maximum; we now integrate from there
    // until we hit the stopping point
    da = a/2.0;
    bool reached_stop = false;
    while (!reached_amax) {

      // Save last value
      double a_last = a;
      double u_last = u;

      // Limit time step
      double a_next = a + da;
      if (a_next > amax) {
	a_next = amax;
	reached_amax = true;
      }

      // Take a step
      int status = gsl_odeiv2_driver_apply(drv, &a, a_next, &u);

      // Check if this step succeeded
      if (status == GSL_SUCCESS) {

	// Success; now check if u has dropped below our tolerance
	if (u < eps) {
	  *q_stop = log10(a-1);
	  break;
	}

	// u is still above our tolerance; if we have not yet hit the
	// stopping point, increase step size
	if (!reached_stop) da *= 2.0;
	
      } else {

	// If we're here, we need to reduce the time step and try
	// again
	reached_stop = true;
	da /= 2.0;
	a = a_last;
	u = u_last;
	reached_amax = false;
	gsl_odeiv2_driver_reset(drv);

      }
      
    } // End integration loop

    // If we hit amax, store values at amax
    if (reached_amax) {
      *q_stop = qmax;
      if (*umax == 0.0) {
	*umax = u;
	*q_umax = qmax;
      }
    }

  } else {
    
    // Case for midx < yidx; in this case we integrate until we reach
    // the maximum radius to which we're willing to go, or until d(log
    // u)/d(log a) drops below our tolerance
    double a = a0;             // Starting position
    double u = u0;             // Starting velocity
    double da = 1.0e-3;        // Starting step size
    bool reached_amax = false; // Flag for if we hit ceiling
    while (!reached_amax) {

      // Store previous values so we can back up when needed
      double a_last = a;
      double u_last = u;

      // Limit time step
      double a_next = a + da;
      if (a_next > amax) {
	a_next = amax;
	reached_amax = true;
      }

      // Take a step
      int status = gsl_odeiv2_driver_apply(drv, &a, a_next, &u);

      // If step was unsuccesful, reduce time step and try again
      if (status != GSL_SUCCESS) {
	a = a_last;
	u = u_last;
	da /= 2.0;
	reached_amax = false;
	gsl_odeiv2_driver_reset(drv);
	continue;
      }

      // Get derivative of velocity at current position
      double dlogu_dloga = (a/u) *
	duda(a, u, par.gex, par.uh, pow(a, par.yidx), pow(a, par.midx));

      // Action depends on value of dlogu_dloga
      if (dlogu_dloga < eps) {

	// d(log u)/d(log a) is below our tolerance, so record result
	// and stop
	*umax = u;
	*q_umax = log10(a-1);
	*q_stop = log10(a-1);
	break;

      } else {

	// d(log u)/d(log a) is still above our tolerance, so increase
	// step size and continue integrating
	da *= 2.0;

      }
      
    } // end integration loop

    // If we terminated because we hit amax, record final values
    if (reached_amax) {
      *q_stop = qmax;
      *q_umax = qmax;
      *umax = u;
    } 
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function to integrate the acceleration curve from a given starting
// value of u and a until the velocity reaches some target value u_stop
/////////////////////////////////////////////////////////////////////////////

bool a_from_u(double *a, double *u, ode_params par,
	      double da, double u_stop,
	      double a_stop, double eps,
	      gsl_odeiv2_driver *drv) {

  // If we have been given a = 1, u = 0 as a starting point, compute a
  // starting point using the series expansion solution near a = 1
  if (*u == 0.0) {
    double gex = par.gex;
    double uh = par.uh;
    double dyda = par.yidx;
    double dmda = par.midx;
    double eps0 = 1.0e-3;
    while (1) {
      *u = u0(eps0, gex, uh, dyda, dmda);
      if (*u <= u_stop/10.0 && eps0 < 0.1*(a_stop-1)) break;
      else eps0 /= 10.0;
    }
    *a = 1.0 + eps0;
    gsl_odeiv2_driver_reset_hstart(drv, eps0);
  }

  // Integrate until we hit target u
  bool overshot = false;
  double u_start = *u;
  while (true) {
  
    // Store previous values so we can back up when needed
    double a_last = *a;
    double u_last = *u;

    // Avoid overshoot
    double a_next = *a + da;
    if (a_next > a_stop) a_next = a_stop;

    // Take a step
    int status = gsl_odeiv2_driver_apply(drv, a, a_next, u);

    // If step was unsuccesful, reduce time step and try again
    if (status != GSL_SUCCESS) {
      *a = a_last;
      *u = u_last;
      da /= 2.0;
      gsl_odeiv2_driver_reset(drv);
      continue;
    }

    // Compare value of u to target value
    if (fabs(*u/u_stop - 1.0) < eps ||
	da/(*a) < eps/1.0e6) {
	
      // We hit target, so we're done
      return true;

    } else if ((*u > u_stop && u_start < u_stop) ||
	       (*u < u_stop && u_start > u_stop)) {

      // We overshot, so back up
      gsl_odeiv2_driver_reset_hstart(drv, 0.01 * (*a - a_last));
      overshot = true;
      *a = a_last;
      *u = u_last;
      da /= 2.0;

    } else if (a_next == a_stop) {

      // We have reached the largest allowed value of a without
      // reaching u_stop, so stop
      return false;

    } else {

      // We have not reached target, do double step size unless we
      // already overshot
      if (!overshot) da *= 2.0;
      
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function to compute a table of Gamma e^-x vs. u and q, where q =
// log10 (a - 1). The arguments are nu and umax, which specify the
// grid in u, qmin, qmax, and nq, which specify the grid in q. The
// returned values are loggex, which gives the values of
// log10(Gamma e^-x) at each point of the (u, q) grid, and umin, which
// gives the minimum value of u for which a solution exists at each q;
// value in the loggex grid for values of u < umin are filled with
// 0.0 as a flag value. The function returns estimates of the mean and
// maximum interpolation errors.
/////////////////////////////////////////////////////////////////////////////

void build_gex_u_q_grid(double uh, double yidx, double midx,
			int nu, double umax,
			int nq, double qmin, double qmax,
			double *loggex, double *umin,
			double *mean_err, double *max_err) {
  
  // Build GSL machinery we will be using
  ode_params par = { 0.0, uh, yidx, midx };
  gsl_odeiv2_system sys = { &duda_gsl, NULL, 1, &par };
  gsl_odeiv2_driver *drv =
    gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
				  0.1*(1.0+pow(10., qmin)),
				  eps_ode, eps_ode);
  gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

  // Fill values on grid
  if (verbose) printf("Bulding loggex grid\n");
  for (int j=0; j<nq; j++) {
    double q = qmin + j * (qmax-qmin) / (nq-1.0);
    double a = 1.0 + pow(10., q);
    if (verbose)
      printf("   ...row %d / %d at q = %f\n", j+1, nq, q);
    for (int i=0; i<nu; i++) {
      double u = umax * i / (nu-1.0);
      loggex[i + nu*j] =
	log10(gex_solve(a, u, eps_solve, &sys, drv, solver));
    }

    // Get umin if needed
    if (yidx > midx) {
      par.gex = 1.0;
      umin[j] = ua(a, fmin(1.0e-3, 0.1*umax/nu), &sys, drv);
    } else {
      umin[j] = 0.0;
    }
  }

  // Estimate interpolation error by compute the difference between
  // the value at every grid point and the results of a linear
  // interpolation from its neighbours
  *mean_err = *max_err = 0.0;
  for (int j=1; j<nq-1; j++) {
    for (int i=1; i<nu-1; i++) {
      double loggex_interp = 0.25 *
	(loggex[i-1 + nu*(j-1)] +
	 loggex[i+1 + nu*(j-1)] +
	 loggex[i-1 + nu*(j+1)] +
	 loggex[i+1 + nu*(j+1)]);
      double err = fabs(loggex_interp - loggex[i + nu*j]);
      *max_err = err > *max_err ? err : *max_err;
      *mean_err += err;
    }
  }
  *mean_err /= ((nu-2) * (nq-2));

  // Free GSL machinery
  gsl_root_fsolver_free(solver);
  gsl_odeiv2_driver_free(drv);
}


/////////////////////////////////////////////////////////////////////////////
// Function to generate a grid of q(u/umax, log10(gex)) values. The
// inputs are the number of grid points in u, nu, the minimum values,
// maximum values, the number of grid points and maximum value in in
// loggex, ngex and loggex_max, and the maximum q value at which we will stop
// integrating, qmax. The outputs are a grid of q values, together
// with auxiliary arrays of size ngex giving the value of the critical
// points on the acceleration curve for each gex value. The function
// returns estimates of the mean and maximum interpolation errors.
/////////////////////////////////////////////////////////////////////////////

void build_q_gex_u_grid(double uh, int yidx, int midx,
			int nu,
			int ngex, double loggex_max,
			double qmax,
			double *q, double *umax, double *q_umax,
			double *q_stop,
			double *mean_err, double *max_err,
			double *max_err_u, double *max_err_loggex,
			bool gex_lo) {
  
  // Build GSL machinery we will be using
  ode_params par = { 0.0, uh, yidx, midx };
  gsl_odeiv2_system sys = { &duda_gsl, NULL, 1, &par };
  gsl_odeiv2_driver *drv =
    gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
				  eps_solve, eps_ode, eps_ode);

  // Fill values on grid
  if (verbose) printf("Bulding q grid\n");
  for (int j=0; j<ngex; j++) {
    double loggex;
    if (!gex_lo) loggex = j * loggex_max / (ngex-1.0);
    else loggex = log10(1.0 + j * (pow(10., loggex_max) - 1.0) /
			(ngex - 1));
    if (verbose)
      printf("   ...row %d / %d at loggex = %f\n", j+1, ngex, loggex);

    // Handle special cases gex = 1 and yidx <= midx, in which case
    // the wind fails to launch at all
    if (j == 0 && yidx <= midx) {
      umax[j] = 0;
      q_umax[j] = -DBL_MAX;
      q_stop[j] = -DBL_MAX;
      for (int i=0; i<nu; i++) q[i + nu*j] = -DBL_MAX;
      if (yidx < midx) 
	for (int i=0; i<nu; i++) q[ngex*nu + i + nu*j] = -DBL_MAX;
      continue;
    } 

    // Get critical points at this loggex value
    par.gex = pow(10., loggex);
    compute_critical_points(qmax, eps_solve, par,
			    umax+j, q_umax+j, q_stop+j,
			    drv);

    // Set values at minimum and maximum u
    q[nu*j] = -DBL_MAX;
    q[nu*j + nu-1] = q_umax[j];

    // Loop over u values
    double amax = 1.0 + pow(10., q_umax[j]);
    double a = 1.0;
    double u = 0.0;
    for (int i=1; i<nu-1; i++) {
      double u_target = umax[j] * i / (nu-1.0);
      double da;
      if (i != 1) da = pow(10., q[nu*j + i-1]) - pow(10., q[nu*j + i-2]);
      else da = 1.0e-3;
      bool success = a_from_u(&a, &u, par, da, u_target, amax,
			      eps_solve, drv);
      if (success) q[nu*j + i] = log10(a-1);
      else {
	// If we hit the stopping point, stop and fill remaining
	// points with max value	
	for (int k=i; k<nu-1; k++) q[nu*j + k] = q_umax[j];
	break;
      }
    }

    // If we have yidx < midx, we now need to repeat this process with
    // decreasing velocity
    if (yidx < midx) {
      int off = nu * (ngex + j);

      // Set values at end points
      q[off] = q_stop[j];
      q[off + nu-1] = q_umax[j];

      // Now integrate the decelerating part of the velocity curve
      a = 1.0 + pow(10., q_umax[j]);
      u = umax[j];
      amax = 1.0 + pow(10., q_stop[j]);
      for (int i=nu-2; i>0; i--) {
	double u_target = umax[j] * i / (nu-1.0);
	double da;
	if (i != nu-2) da = pow(10., q[off + i+1]) - pow(10., q[off + i+2]);
	else da =  pow(10., q[nu*j + nu-1]) - pow(10., q[nu*j + nu-2]);
	if (da == 0) da = 1.0e-3;
	bool success = a_from_u(&a, &u, par, da, u_target, amax,
				eps_solve, drv);
	if (success) {
	  q[off + i] = log10(a-1);
	  u = u_target;
	} else {
	  // If we hit the stopping point, stop and fill remaining
	  // points with max value	
	  for (int k=i; k>0; k--) q[off + k] = q_stop[j];
	  break;
	}
      }
    }
  }

  // Estimate interpolation error by compute the difference between
  // the value at every grid point and the results of a linear
  // interpolation from its neighbours; note that we have a special
  // case here: if yidx = midx, then the function q(u/umax, gex)
  // becomes undefined as gex -> 1, because umax -> 0. Evaluating the
  // error at the gex -> 1 boundary therefore results in an error that
  // never converges no matter what resolution is used; we therefore
  // exclude j = 1 from the error analysis for this case only
  *mean_err = *max_err = 0.0;
  *max_err_u = *max_err_loggex = -1.0;
  int jstart = 1;
  if (yidx == midx) jstart = 2;
  for (int j=jstart; j<ngex-1; j++) {
    for (int i=1; i<nu-1; i++) {
      double q_interp = 0.25 *
	(q[i-1 + nu*(j-1)] +
	 q[i+1 + nu*(j-1)] +
	 q[i-1 + nu*(j+1)] +
	 q[i+1 + nu*(j+1)]);
      double a_interp = 1.0 + pow(10., q_interp);
      double a = 1.0 + pow(10., q[i + nu*j]);
      double err = fabs(a - a_interp) / fmin(a, a_interp);
      *mean_err += err;
      if (err > *max_err) {
	*max_err = err;
	*max_err_u = i / (nu-1.0);
	*max_err_loggex = j * loggex_max / (ngex-1.0);
      }
      if (yidx < midx) {
	int off = nu * (ngex+j) + i;
	q_interp = 0.25 *
	  (q[off - 1 - nu] +
	   q[off + 1 - nu] +
	   q[off - 1 + nu] +
	   q[off + 1 + nu]);
	a_interp = 1.0 + pow(10., q_interp);
	a = 1.0 + pow(10., q[off]);
	err = fabs(a - a_interp) / fmin(a, a_interp);
	*mean_err += err;
	if (err > *max_err) {
	  *max_err = err;
	  *max_err_u = i / (nu-1.0);
	  *max_err_loggex = j * loggex_max / (ngex-1.0);
	}
      }
    }
  }
  *mean_err /= ((nu-2) * (ngex-2));
  if (yidx < midx) *mean_err /= 2.0;

  // Free GSL machinery
  gsl_odeiv2_driver_free(drv);
}


/////////////////////////////////////////////////////////////////////////////
// Main function
/////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

  // Parse inputs
  char usage[] = "Usage: hot_wind_tab uh yidx midx nu nq qmin qmax ngex loggex_max ngex_lo [-v, --verbose] [-ovr, --overwrite] [-d, --dir dir] [-a, --ascii]\n   uh = hot gas velocity\n   yidx = expansion parameter for clouds, y = a^yidx\n   midx = gravitational potential parameter, m = a^idx\n   nu = number of grid points in velocity\n   nq = number of grid points in q, where q = log_10(a - 1)\n   qmin = minimum value of q in grid\n   qmax = maximum value of q in grid\n   ngex = number of grid points in log_10(Gamma e^-x)\n   loggex_max = maximum value of log_10(Gamma e^-x) in grid\n   ngex_lo = number of grid points in the grid of log_10(Gamma e^-x) near 0\n   --verbose = produce verbose output\n   --ovewrite = overwrite existing output files (default behavior is to skip if existing output is found)\n   --dir = write output to directory dir (default is run directory)\n   --ascii = write ASCII output (default is raw binary output)\n";
  if (argc < 11 || argc > 16) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  double uh = atof(argv[1]);
  int yidx = atoi(argv[2]);
  int midx = atoi(argv[3]);
  int nu = atoi(argv[4]);
  int nq = atoi(argv[5]);
  double qmin = atof(argv[6]);
  double qmax = atof(argv[7]);
  int ngex = atoi(argv[8]);
  double loggex_max = atof(argv[9]);
  int ngex_lo = atoi(argv[10]);
  bool overwrite = false;
  bool ascii = false;
  char *dir = NULL;
  verbose = false;
  for (int n=11; n<argc; n++) {
    if (strcmp(argv[n], "-v") == 0 ||
	strcmp(argv[n], "--verbose") == 0) {
      verbose = true;
    } else if (strcmp(argv[n], "-ovr") == 0 ||
	       strcmp(argv[n], "--overwrite") == 0) {
      overwrite = true;
    } else if (strcmp(argv[n], "-a") == 0 ||
	       strcmp(argv[n], "--ascii") == 0) {
      ascii = true;
    } else if (strcmp(argv[n], "-d") == 0 ||
	       strcmp(argv[n], "--dir") == 0) {
      if (n == argc-1) {
	fprintf(stderr, "%s", usage);
	exit(1);
      } else {
	n++;
	// Ensure dir is properly terminated by /
	int l = strlen(argv[n]);
	if (argv[n][l-1] == '/') {
	  dir = malloc(l+1);
	  sprintf(dir, "%s", argv[n]);
	} else {
	  dir = malloc(l+2);
	  sprintf(dir, "%s/", argv[n]);
	}
      }
    } else {
      fprintf(stderr, "%s", usage);
      exit(1);
    }
  }

  // Check that parameters are valid
  if (uh <= 0.0 || yidx < 0 || yidx > 2 || midx < 0 || midx > 1 ||
      qmin >= qmax || loggex_max <= 0.0 || nu < 4 || nq < 4 || ngex < 4 ||
      ngex_lo < 4) {
    fprintf(stderr, "Parameter out of range; valid ranges are:\n");
    fprintf(stderr, "   uh > 0\n");
    fprintf(stderr, "   yidx = 0, 1, or 2\n");
    fprintf(stderr, "   midx = 0 or 1\n");
    fprintf(stderr, "   nu > 3\n");
    fprintf(stderr, "   nq > 3\n");
    fprintf(stderr, "   qmin < qmax\n");
    fprintf(stderr, "   ngex > 3\n");
    fprintf(stderr, "   loggex_max > 0\n");
    fprintf(stderr, "   ngex_lo > 3\n");
    exit(1);
  }

  // q(u, gex) calculation

  // Construct output filenames
  char *fname1, *fname2, ext[4];
  if (ascii) sprintf(ext, "txt"); else sprintf(ext, "bin");
  if (dir) {
    fname1 = malloc(strlen(dir) + 100);
    fname2 = malloc(strlen(dir) + 100);
    sprintf(fname1, "%sqtab_u_uh%f_y%d_m%d.%s", dir, uh, yidx, midx, ext);
    sprintf(fname2, "%sqtab_q_uh%f_y%d_m%d.%s", dir, uh, yidx, midx, ext);
  } else {
    fname1 = malloc(100);
    fname2 = malloc(100);
    sprintf(fname1, "qtab_u_uh%f_y%d_m%d.%s", uh, yidx, midx, ext);
    sprintf(fname2, "qtab_q_uh%f_y%d_m%d.%s", uh, yidx, midx, ext);
  }

  // Check if files exist; if it does and overwriting is not set, skip
  FILE *fp1, *fp2, *fp3;
  fp1 = fopen(fname1, "r");
  fp2 = fopen(fname2, "r");
  if (fp1 && fp2 && !overwrite) {

    // Just close
    fclose(fp1);
    fclose(fp2);
    
  } else {
    
    // Files do not exist, or we are overwriting existing ones
    if (fp1) fclose(fp1);
    if (fp2) fclose(fp2);

    // Allocate memory to hold q(u, gex) results
    int nc = yidx < midx ? 2 : 1;
    double *q = calloc(nc*nu*ngex, sizeof(double));
    double *umax = calloc(ngex, sizeof(double));
    double *q_umax = calloc(ngex, sizeof(double));
    double *q_stop = calloc(ngex, sizeof(double));
    if (!q || !umax || !q_umax || !q_stop) {
      fprintf(stderr, "Error: unable to allocate memory to hold result\n");
      exit(1);
    }
  
    // Build q(u, gex) grid
    double mean_err, max_err, max_err_u, max_err_loggex;
    build_q_gex_u_grid(uh, yidx, midx, nu, ngex, loggex_max, qmax,
		       q, umax, q_umax, q_stop,
		       &mean_err, &max_err,
		       &max_err_u, &max_err_loggex, false);
    if (verbose)
      printf("Mean error = %f, max error = %f at u/umax = %f, loggex = %f\n",
	     mean_err, max_err, max_err_u, max_err_loggex);

    // Write output
    fp1 = fopen(fname1, "w");
    fp2 = fopen(fname2, "w");
    if (ascii) {

      // ASCII mode
      for (int j=0; j<ngex; j++) {
	double loggex = j * loggex_max / (ngex-1);
	fprintf(fp1, "%e   %e   %e   %e\n",
		loggex, umax[j], q_umax[j], q_stop[j]);
      }
      for (int j=0; j<ngex; j++) {
	for (int i=0; i<nu; i++) {
	  fprintf(fp2, "   %e", q[nu*j+i]);
	}
	fprintf(fp2, "\n");
      }
      if (yidx < midx) {
	for (int j=0; j<ngex; j++) {
	  for (int i=0; i<nu; i++) {
	    fprintf(fp2, "   %e", q[ngex*nu + nu*j+i]);
	  }
	  fprintf(fp2, "\n");
	}
      }
    } else {

      // Binary mode
      for (int j=0; j<ngex; j++) {
	double loggex = j * loggex_max / (ngex-1);
	fwrite(&loggex, sizeof(double), 1, fp1);
      }
      fwrite(umax, sizeof(double), ngex, fp1);
      fwrite(q_umax, sizeof(double), ngex, fp1);
      fwrite(q_stop, sizeof(double), ngex, fp1);
      if (yidx >= midx)
	fwrite(q, sizeof(double), ngex*nu, fp2);
      else
	fwrite(q, sizeof(double), 2*ngex*nu, fp2);
    }
    fclose(fp1);
    fclose(fp2);

    // Free memory
    free(q);
    free(umax);
    free(q_umax);
    free(q_stop);
  }

  // Handle q(u, gex) grid for values near gex = 1 -- only required
  // for yidx <= midx
  if (yidx <= midx) {
  
    // Construct output filenames
    char *fname1, *fname2, ext[4];
    if (ascii) sprintf(ext, "txt"); else sprintf(ext, "bin");
    if (dir) {
      fname1 = malloc(strlen(dir) + 100);
      fname2 = malloc(strlen(dir) + 100);
      sprintf(fname1, "%sqtab_gexlo_u_uh%f_y%d_m%d.%s", dir, uh, yidx, midx, ext);
      sprintf(fname2, "%sqtab_gexlo_q_uh%f_y%d_m%d.%s", dir, uh, yidx, midx, ext);
    } else {
      fname1 = malloc(100);
      fname2 = malloc(100);
      sprintf(fname1, "qtab_gexlo_u_uh%f_y%d_m%d.%s", uh, yidx, midx, ext);
      sprintf(fname2, "qtab_gexlo_q_uh%f_y%d_m%d.%s", uh, yidx, midx, ext);
    }

    // Check if files exist; if it does and overwriting is not set, skip
    FILE *fp1, *fp2, *fp3;
    fp1 = fopen(fname1, "r");
    fp2 = fopen(fname2, "r");
    if (fp1 && fp2 && !overwrite) {
      
      // Just close
      fclose(fp1);
      fclose(fp2);
    
    } else {
    
      // Files do not exist, or we are overwriting existing ones
      if (fp1) fclose(fp1);
      if (fp2) fclose(fp2);

      // Allocate memory to hold q(u, gex) results
      int nc = yidx < midx ? 2 : 1;
      double *q = calloc(nc*nu*ngex_lo, sizeof(double));
      double *umax = calloc(ngex_lo, sizeof(double));
      double *q_umax = calloc(ngex_lo, sizeof(double));
      double *q_stop = calloc(ngex_lo, sizeof(double));
      if (!q || !umax || !q_umax || !q_stop) {
	fprintf(stderr, "Error: unable to allocate memory to hold result\n");
	exit(1);
      }
  
      // Build q(u, gex) grid
      double mean_err, max_err, max_err_u, max_err_loggex;
      double loggex_lo_max = loggex_max / (ngex - 1);
      build_q_gex_u_grid(uh, yidx, midx, nu, ngex_lo, loggex_lo_max, qmax,
			 q, umax, q_umax, q_stop,
			 &mean_err, &max_err,
			 &max_err_u, &max_err_loggex, true);
      if (verbose)
	printf("Mean error = %f, max error = %f at u/umax = %f, loggex = %f\n",
	       mean_err, max_err, max_err_u, max_err_loggex);

      // Write output
      fp1 = fopen(fname1, "w");
      fp2 = fopen(fname2, "w");
      if (ascii) {

	// ASCII mode
	for (int j=0; j<ngex_lo; j++) {
	  double loggex = log10(1.0 + j * (pow(10., loggex_lo_max) - 1.0) /
				(ngex_lo - 1));
	fprintf(fp1, "%e   %e   %e   %e\n",
		loggex, umax[j], q_umax[j], q_stop[j]);
	}
	for (int j=0; j<ngex_lo; j++) {
	  for (int i=0; i<nu; i++) {
	    fprintf(fp2, "   %e", q[nu*j+i]);
	  }
	  fprintf(fp2, "\n");
	}
	if (yidx < midx) {
	  for (int j=0; j<ngex; j++) {
	    for (int i=0; i<nu; i++) {
	      fprintf(fp2, "   %e", q[ngex_lo*nu + nu*j+i]);
	    }
	    fprintf(fp2, "\n");
	  }
	}
      } else {

	// Binary mode
	for (int j=0; j<ngex_lo; j++) {
	  double loggex = log10(1.0 + j * (pow(10., loggex_lo_max) - 1.0) /
				(ngex_lo - 1));
	  fwrite(&loggex, sizeof(double), 1, fp1);
	}
	fwrite(umax, sizeof(double), ngex_lo, fp1);
	fwrite(q_umax, sizeof(double), ngex_lo, fp1);
	fwrite(q_stop, sizeof(double), ngex_lo, fp1);
	if (yidx >= midx)
	  fwrite(q, sizeof(double), ngex_lo*nu, fp2);
	else
	  fwrite(q, sizeof(double), 2*ngex_lo*nu, fp2);
      }
      fclose(fp1);
      fclose(fp2);

      // Free memory
      free(q);
      free(umax);
      free(q_umax);
      free(q_stop);
    }
  }

  // gex(u, q) calculation

  // Set file names
  char *fname3;
  if (dir) {
    fname3 = malloc(strlen(dir) + 100);
    sprintf(fname1, "%sgextab_q_uh%f_y%d_m%d.%s", dir, uh, yidx, midx, ext);
    sprintf(fname2, "%sgextab_u_uh%f_y%d_m%d.%s", dir, uh, yidx, midx, ext);
    sprintf(fname3, "%sgextab_gex_uh%f_y%d_m%d.%s", dir, uh, yidx, midx, ext);
  } else {
    fname3 = malloc(100);
    sprintf(fname1, "gextab_q_uh%f_y%d_m%d.%s", uh, yidx, midx, ext);
    sprintf(fname2, "gextab_u_uh%f_y%d_m%d.%s", uh, yidx, midx, ext);
    sprintf(fname3, "gextab_gex_uh%f_y%d_m%d.%s", uh, yidx, midx, ext);
  }
  
  // Check if files exist; if it does and overwriting is not set, skip
  fp1 = fopen(fname1, "r");
  fp2 = fopen(fname2, "r");
  fp3 = fopen(fname3, "r");
  if (fp1 && fp2 && fp3 && !overwrite) {

    // Just close
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    
  } else {
    
    // Files do not exist, or we are overwriting existing ones
    if (fp1) fclose(fp1);
    if (fp2) fclose(fp2);
    if (fp3) fclose(fp3);

    // Allocate memory to hold gex(u, q) results
    double *loggex = calloc(nu*nq, sizeof(double));
    double *umin = calloc(nq, sizeof(double));
    if (!loggex || !umin) {
      fprintf(stderr, "Error: unable to allocate memory to hold result\n");
      exit(1);
    }

    // Build gex(u, q) grid
    double ufac = 1.0 - 100.*eps_solve;
    double mean_err, max_err;
    build_gex_u_q_grid(uh, yidx, midx, nu, ufac*uh, nq, qmin, qmax,
		       loggex, umin,
		       &mean_err, &max_err);
    if (verbose)
      printf("Mean error = %f, max error = %f\n",
	     mean_err, max_err);

    // Write results
    fp1 = fopen(fname1, "w");
    fp2 = fopen(fname2, "w");
    fp3 = fopen(fname3, "w");
    if (ascii) {

      // ASCII mode
      for (int j=0; j<nq; j++) {
	double q = qmin + j * (qmax-qmin) / (nq-1.0);
	fprintf(fp1, "%e   %f\n", q, umin[j]);
      }
      for (int i=0; i<nu; i++) {
	double u = i * ufac*uh / (nu-1.0);
	fprintf(fp2, "%e\n", u);
      }
      for (int j=0; j<nq; j++) {
	for (int i=0; i<nu; i++) {
	  fprintf(fp3, "   %e", loggex[j*nu + i]);
	}
	fprintf(fp3, "\n");
      }
    } else {

      // Binary mode
      for (int j=0; j<nq; j++) {
	double q = qmin + j * (qmax-qmin) / (nq-1.0);
	fwrite(&q, sizeof(double), 1, fp1);
      }
      fwrite(umin, sizeof(double), nq, fp1);
      for (int i=0; i<nu; i++) {
	double u = i * ufac*uh / (nu-1.0);
	fwrite(&u, sizeof(double), 1, fp2);
      }
      fwrite(loggex, sizeof(double), nu*nq, fp3);
    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);

    // Free memory
    free(umin);
    free(loggex);
  }

  // Free memory
  free(fname1);
  free(fname2);
  free(fname3);
}
