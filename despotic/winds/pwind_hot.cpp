#include "pwind_hot.H"
#include "pwind_util.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <algorithm>
#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Functions and data structures used by the GSL
////////////////////////////////////////////////////////////////////////

// This returns
// du/da = 1/sqrt(2 u a^2) [Gamma exp(-x) y (1 - u/uh)^2 - m]
static int duda(double a, const double u[], double duda[],
		void *params) {
  struct duda_params *par = (struct duda_params *) params;
  const pwind_hot *w = static_cast<const pwind_hot *>(par->w);
  double x = par->x;
  double uh = w->getUh();
  double Gamma = w->getGamma();
  double y = w->y(a);
  double m = w->m(a);
  double gex = max(Gamma*exp(-x), 1.0); // Prevent problems due to
					// roundoff making this
					// slightly less than 1

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
  if (u[0] < uh)
    duda[0] = 1.0 / (2.0*SQR(a)*fabs(u[0])) *
      (gex * y * SQR(1.0 - u[0]/uh) - m);
  else
    duda[0] = 0.0;
  if (u[0] < 0.0 && duda[0] <= 0.0)
    return GSL_EBADFUNC;
  else
    return GSL_SUCCESS;
}

struct xmin_params {
  double abstol, sx, zeta_M;
};

double xmin_resid(double x, void *params) {
  struct xmin_params *par = (struct xmin_params *) params;
  double ret=zetaM(x, par->sx) / (par->zeta_M*par->abstol) - 1.0;
  return ret;
}

struct alim_params {
  double x, u2, vp2;
  const pwind_hot *w;
};

static double alim_func(double loga, void *params) {
  struct alim_params *par = (struct alim_params *) params;
  double a = exp(loga);
  double ret = par->u2 -
    par->w->U2(par->x, a) * (1.0 - par->vp2/SQR(a));
  return ret;
}

struct alim_quart_params {
  double u2, vp2, f2;
};

static double alim_quart_func(double a, void *params) {
  struct alim_quart_params *par = (struct alim_quart_params *) params;
  return SQR(a-1.0)*par->f2 - par->u2 / (1.0 - par->vp2/SQR(a));
}

struct u_max_params {
  double vp2, x;
  const pwind_hot *w;
};

static double u_max_func(double loga, void *params) {
  struct u_max_params *par = (struct u_max_params *) params;
  double a = exp(loga);
  double ret = -(1.0-par->vp2/SQR(a)) * par->w->U2(par->x, a);
  return ret;
}


////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////

pwind_hot::pwind_hot(const double Gamma_, const double mach_,
		     const double uh_, const double yoverm_,
		     const pwind_potential *potential_,
		     const pwind_expansion *expansion_,
		     const pwind_geom *geom_,
		     const double epsabs_,
		     const double epsrel_, const double interpabs_,
		     const double interprel_, const double fcrit_,
		     const double amax_grid_) :
  pwind(Gamma_, mach_, potential_, expansion_, geom_,
	epsabs_, epsrel_, fcrit_),
  uh(uh_), yoverm(yoverm_), interpabs(interpabs_), interprel(interprel_),
  amax_grid(amax_grid_)
{
  xcrit = log(Gamma_);
  zeta_M = zetaM(xcrit + log(fcrit_), sx);
  zeta_A = zetaA(xcrit + log(fcrit_), sx);
}

////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////
pwind_hot::~pwind_hot() {
  // Free GSL interpolation apparatus
  gsl_spline2d_free(U2_interp);
  gsl_spline2d_free(X_interp);
  gsl_spline_free(Umin_interp);
  gsl_spline_free(Umax_interp);
  gsl_interp_accel_free(U2_acc_x);
  gsl_interp_accel_free(U2_acc_a);
  gsl_interp_accel_free(X_acc_u);
  gsl_interp_accel_free(X_acc_a);
  gsl_interp_accel_free(Umin_acc);
  gsl_interp_accel_free(Umax_acc);
  gsl_spline_free(aMin_interp[0]);
  gsl_interp_accel_free(aMin_acc[0]);
  if (aMin_interp[1] != nullptr) {
    gsl_spline_free(aMin_interp[1]);
    gsl_interp_accel_free(aMin_acc[1]);
  }
  if (aMax_interp != nullptr) {
    gsl_spline_free(aMax_interp);
    gsl_interp_accel_free(aMax_acc);
  }
  if (a_stop_interp != nullptr) {
    gsl_spline_free(a_stop_interp);
    gsl_spline_free(x_stop_interp);
    gsl_interp_accel_free(a_stop_acc);
    gsl_interp_accel_free(x_stop_acc);
  }
}


////////////////////////////////////////////////////////////////////////
// Routines to evaluate the kinematics using the interpolation
// functions; note that we handle the case of log(a) < log(a_grid[0])
// or u < u_grid[0] using our series expansions:
//
// for x != xcrit, near a = 1 we have
// u = sqrt((Gamma e^-x - 1) (a-1))
// for x = xcrit, near a = 1 we have
// u = (a-1)/(2 uh) * (-1 + sqrt(1 + 2*uh**2 (y' - m')))
//
// We also impose safety cuts on x. We should never be calling these
// routines with x outside the allowed range, but roundoff error can
// cause the value of x to wander outside the allowed range by an
// amount of order machine precision. To avoid calculations dying
// when this happens, we check x and correct it if necessary.
////////////////////////////////////////////////////////////////////////
inline double
pwind_hot::X(const double ur, const double a) const {
  double loga = log(a);
  if (loga < loga_grid.front() || ur < u_grid.back()[0]) {
    // Analytic result for a -> 0
    return -log((SQR(ur)/(a-1.0) + 1.0)/Gamma);
  } else {
    return gsl_spline2d_eval(X_interp, ur, min(loga, loga_grid.back()),
			     X_acc_u, X_acc_a);
  }
}
inline double
pwind_hot::U2(const double x, const double a) const {
  double loga = log(a);
  if (loga > loga_grid[0]) {
    double xtmp = x;
    xtmp = max(xtmp, x_grid[0]+1.0e-10);
    // The 1.0e-10 prevents catastrophic failures due to roundoff
    // problems in the gsl inerpolation routine
    xtmp = min(xtmp, x_grid.back());
    double u2 = gsl_spline2d_eval(U2_interp, xtmp,
				  min(loga, loga_grid.back()),
				  U2_acc_x, U2_acc_a);
    return u2;
  } else {
    if (x != xcrit) {
      return (Gamma*exp(-x)-1.0)*(a-1.0);
    } else {
      double u = (a-1.0) / (2.0*uh) *
	(-1.0 + sqrt(1.0+2.0*SQR(uh)*(dyda(1.0) - dmda(1.0))));
      return SQR(u);
    }
  }
}
inline double
pwind_hot::dU2dx(const double x, const double a) const {
  double loga = log(a);
  if (loga > loga_grid[0]) {
    double xtmp = x;
    xtmp = max(xtmp, x_grid[0]);
    xtmp = min(xtmp, x_grid.back());
    double val = gsl_spline2d_eval_deriv_x(U2_interp, xtmp,
					   min(loga, loga_grid.back()),
					   U2_acc_x, U2_acc_a);
    // Safety check here: we can get junk here because we've
    // interpolated to get x, and that interpolation has given us a
    // combination (x, a) that is not actually allowed, or right up
    // against the boundary of the disallowed region. This happens
    // most frequently at small velocities, because for some
    // potentials and expansions, u(a) at fixed x can go to zero at
    // finite a with a derivative that diverges; when we try to
    // interpolate back our numerical approximation to this behaviour
    // to get x(u, a), the resulting x is not very accurate. The way
    // we check for this is we see if any of the four points that
    // define our interpolation box have a value of -1. If so, this
    // indicates that we are too close to the edge of the allowed (x,
    // a) region to give a meaningful value of dU2dx. The safest
    // course of action in this case is to return dU2/dx = -infinity,
    // which becomes true as we approach the boundary. Since dU2dx
    // generally appears in the denominator in the integrands we want
    // to evaluate, this choice forces the integrand to zero, and
    // gives the best approximation we're going to get in these odd
    // corners of parameter space.
    for (int i=0; i<=1; i++) {
      for (int j=0; j<=1; j++) {
	if (u_grid[U2_acc_x->cache+i][U2_acc_a->cache+j] == -1.0)
	  return -numeric_limits<double>::max();
      }
    }
    return val;
  } else {
    return -Gamma*exp(-x)*(a-1.0);
  }
}
inline double
pwind_hot::dU2da(const double x, const double a) const {
  return (Gamma*y(a)*exp(-x)*SQR(1.0-U(x,a)/uh) - m(a)) / SQR(a);
}

////////////////////////////////////////////////////////////////////////
// Routines to construct the interpolation table
////////////////////////////////////////////////////////////////////////

// This routine solves for U(x, a) at a specified value of x and a
// specified vector of log radii loga; the first element of loga is
// required to be 0. If istart > 0, then the integration is
// started from element istart of the vector loga, and the returned
// array is appropriately shortened; the velocity u0 at the first grid
// point must be supplied in this case. Integration goes to the end of
// loga or to iend, whichever is smaller. The values of the velocity
// are returned by the function, and the final radius reached is
// returned in loga_stop. If the integration encounters a singularity
// before reaching the end of the integration, then all grid points
// after the singularity will be set to -1.
vector<double>
pwind_hot::Ua_solve_tab(const double x,
			const vector<double>& loga,
			double& loga_stop,
			const vector<double>::size_type istart,
			const vector<double>::size_type iend_,
			const double u0) const {

  // See where to end
  vector<double>::size_type iend =
    iend_ == 0 ? loga.size()-1 : min(iend_, loga.size()-1);

  // Allocate array to hold results
  vector<double> u(iend-istart+1);

  // Handle the special case x = xcrit and dy/da <= dm/da
  if (x == xcrit) {
    if (dyda(1.0) < dmda(1.0)) {
      // In this case the gas has u = 0 at a = 1, and the velocity is
      // undefined elsewhere, which we flag by returning u = -1
      for (vector<double>::size_type i=0; i<iend-istart; i++) {
	if (i==0 && istart==0) u[istart+i] = -1.0;
	else u[i] = 0.0;
      }
      loga_stop = 0.0;
      return u;
    } else if (dyda(1.0) == dmda(1.0)) {
      // In this case the gas has u = 0 at all radii, so return that
      for (vector<double>::size_type i=0; i<iend-istart; i++)
	u[i] = 0.0;
      loga_stop = loga[iend];
      return u;
    }
  }

  // Set up initial point; if istart = 0, then find it by series
  // expansion
  double Ua[1];
  if (istart == 0) {
    if (x == xcrit)
      Ua[0] = (exp(loga[0])-1.0) / (2.0*uh) *
	(-1.0 + sqrt(1.0+2.0*SQR(uh)*(dyda(1.0)-dmda(1.0))));
    else
      Ua[0] = sqrt((Gamma*exp(-x)-1.0)*(exp(loga[0])-1.0));
  } else {
    Ua[0] = u0;
  }
  u[0] = Ua[0];

  // Catch case where u[0] < -1, indicating that we're starting past
  // the turnaround point
  if (u[0] < 0.0) {
    for (vector<double>::size_type i=1; i<u.size(); i++)
	u[i] = -1.0;
    return u;
  }
  
  // Integrate
  par.x = x;
  double a = exp(loga[istart]);
  gsl_odeiv2_driver_reset(ode_drv);
  for (vector<double>::size_type i=istart+1; i<=iend; i++) {
    int status = gsl_odeiv2_driver_apply(ode_drv, &a, exp(loga[i]), Ua);
    if (status == GSL_EBADFUNC) {
      // u is approaching 0, so the equation is singular. Inch the
      // solution forward until the final value of u is less than the
      // interpolation tolerance
      double da = (exp(loga[i]) - a)/2.0;
      while (1) {
	// Make sure not to overshoot target
	if (a + da > exp(loga[i])) da = exp(loga[i]) - a;
	// Try to advance with this step
	status = gsl_odeiv2_driver_apply(ode_drv, &a, a+da, Ua);
	// Check if we managed to get to the endpoint of this
	// interval, of if we reached a velocity small enough to stop
	if (a >= exp(loga[i]) || Ua[0] < interpabs) break;
	// If neither stopping condition is met, continue, reducing
	// time step if this advanced failed
	if (status == GSL_EBADFUNC) {
	  da /= 2.0;
	  gsl_odeiv2_driver_reset(ode_drv);
	}
	// Stop if a is no longer changing
	if (da/a < epsabs || da < epsrel) break;
      }
      // Two possibilities: we managed to reach the endpoint, or we
      // did not; in the former case, we just continue as if nothing
      // happened, while in the latter case we halt
      if (a < exp(loga[i])) {
	// Store final radius, flag the remaining grid points, and
	// return
	loga_stop = log(a);
	for (vector<double>::size_type j=i; j<=iend; j++)
	  u[j-istart] = -1.0;
	return u;
      }
    }
    u[i-istart] = Ua[0];
  }
  loga_stop = loga[iend];

  // Return
  return u;
}

// This routine solves for U(x, a) at a specified value of x. It
// accepts an initial set of radial points at which to solve, but
// implements automatic error control and refines the grid
// automatically if necessary to achieve the required interpolation
// accuracy on the final grid. The function returns the tabulated U
// values. The input a values form the initial grid, and on return the
// grid will have been expanded if necessary to achieve the required
// error tolerance; the new grid is returned in loga.
vector<double>
pwind_hot::Ua_solve_auto(const double x,
			 vector<double>& loga,
			 double& loga_stop) const {
  
  // Start by solving for the acceleration law on the grid we've been
  // given
  vector<double> u = Ua_solve_tab(x, loga, loga_stop);

  // If the wind acceleration law for this x is all zeros, we're done
  vector<double>::size_type ctr;
  for (ctr=0; ctr<u.size(); ctr++) if (u[ctr] > 0.0) break;
  if (ctr == u.size()) return u;

  // First check if we need to extend the grid. This check takes two
  // forms, depending on the value of y/m as a -> infinity:
  //
  // Case 1: for Gamma exp(-x) y/m >= 1 as a->inf, the velocity
  // should asymptote to a constant (which is equal to uh if the
  // limiting value is > 1, and is between 0 and uh if the limiting
  // value is exactly 1). In this case, we check for convergence by
  // fitting the last two data points to a functional form u = u0 +
  // u1/a, and comparing the results of that to a zeroth-order
  // extraplation that u = constant past the edge of our grid.
  //
  // Case 2: for Gamma exp(-x) y/m < 1 as a->inf, the velocity reaches
  // 0 at finite a. In this case we want to integrate far enough that
  // we reach the singular point where u -> 0, or until we get to the
  // maximum radius to which we're willing to integrate
  bool case1 = Gamma*exp(-x)*yoverm >= 1;
  while (1) {

    // Grab the end of the grid
    vector<double>::size_type i = u.size()-1;
    double a1 = exp(loga[i]);

    if (case1) {
      // Case 1
    
      // Find estimate of u as a->infinity fitting to u = u0 + u1/a
      if (a1 >= amax_grid) break; // Stop if we hit the absolute limit
      double a0 = exp(loga[i-1]);
      double uinf = (a1*u[i] - a0*u[i-1]) / (a1-a0);

      // Compare to constant extraplation to infinity
      double abserr = fabs(u[i] - uinf);
      double relerr = abserr / uinf;

      // If error is small enough, we're done
      if (abserr <= interpabs || relerr <= interprel) break;

    } else {
      // Case 2
      if (loga_stop <= loga.back()) break;
    }

    // If we're here, we need to extend the grid; grow it by a factor
    // of 2, and extend the solution
    double loga_new = loga.back() + log(2.0);
    if (loga_new > log(amax_grid)) loga_new = log(amax_grid);
    loga.push_back(loga_new);
    vector<double> u_new =
      Ua_solve_tab(x, loga, loga_stop, i, i+1, u[i]);
    u.push_back(u_new[1]);
  }
  
  // March through the acceleration law we have obtained, checking the
  // error in each interval, and adding points where necessary to
  // reduce it; be sure to handle properly cases where the velocity is
  // going to zero
  while (1) {

    // Grab the data we want to check
    vector<double> loga_tmp, u_tmp;
    if (loga_stop == loga.back()) {
      // Use all the data, because u never reaches 0
      loga_tmp = loga;
      u_tmp = u;
    } else {
      // The velocity does go to zero, so construct a temporary vector
      // of loga and u values, including the stopping point
      vector<double>::size_type i;
      for (i=0; loga[i] < loga_stop; i++) {
	loga_tmp.push_back(loga[i]);
	u_tmp.push_back(u[i]);
      }
      loga_tmp.push_back(loga_stop);
      u_tmp.push_back(u[i]);
    }
    // If the accelerate law doesn't contain enough points to form a
    // spline, because we hit the stopping point too soon, bail out
    // here
    if (loga_tmp.size() < gsl_interp_type_min_size(gsl_interp_cspline))
      break;
    // Set up the linear and cubic interpolators on the data; we will
    // use these to estimate the error    
    gsl_interp *interp_lin =
      gsl_interp_alloc(gsl_interp_linear, loga_tmp.size());
    gsl_interp *interp_cube =
      gsl_interp_alloc(gsl_interp_akima, loga_tmp.size());
    gsl_interp_accel *acc_lin = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_cube = gsl_interp_accel_alloc();
    gsl_interp_init(interp_lin, loga_tmp.data(), u_tmp.data(),
		    loga_tmp.size());
    gsl_interp_init(interp_cube, loga_tmp.data(), u_tmp.data(),
		    loga_tmp.size());

    // Check all intervals to find any that require subdivision
    vector<vector<double>::size_type> intervals;
    for (vector<double>::size_type i=0; i<loga_tmp.size()-1; i++) {
      unsigned int ncheck = 8;
      for (unsigned int j=1; j<ncheck; j++) {
	double loga_pt = (ncheck-j)*loga_tmp[i]/((float) ncheck) +
	  j*loga_tmp[i+1]/((float) ncheck);
	double u1 = gsl_interp_eval(interp_lin, loga_tmp.data(),
				    u_tmp.data(),
				    loga_pt, acc_lin);
	double u2 = gsl_interp_eval(interp_cube, loga_tmp.data(),
				    u_tmp.data(),
				    loga_pt, acc_cube);
	double abserr = fabs(u1-u2);
	double relerr = abserr/fabs(u1);
	if (abserr > interpabs && relerr > interprel) {
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
    // points in the middle of them; solve for the new velocities at
    // these points by integrating from the previous ones
    for (vector<double>::size_type i=0; i<intervals.size(); i++) {
      vector<double> u_copy = u;
      vector<double> loga_copy = loga;
      double loga_stop_dummy;
      loga.insert(loga.begin()+i+intervals[i]+1, 
		  0.5*(loga[i+intervals[i]] +
		       loga[i+intervals[i]+1]));
      vector<double> u_new =
	Ua_solve_tab(x, loga, loga_stop_dummy, i+intervals[i],
		     i+intervals[i]+1, u[i+intervals[i]]);
      u.insert(u.begin()+i+intervals[i]+1, u_new[1]);
    }
  }

  // Return
  return u;
}

// This routine adds a row to the interpolation table at a specified
// x; if this new x requird an expanded grid in log a, it recomputes
// the remainder of the table rows on the expanded grid as well
void pwind_hot::add_table_row(const double x) {

  // First compute the table row for this x on the existing grid
  vector<double> loga_new(loga_grid);
  double loga_stop;
  vector<double> u_row = Ua_solve_auto(x, loga_new, loga_stop);

  // Find where in the table this row should be inserted
  vector<double>::size_type idx;
  for (idx = 0; idx<x_grid.size(); idx++)
    if (x < x_grid[idx]) break;

  // Insert the new x and the new velocities into the table
  x_grid.insert(x_grid.begin()+idx, x);
  u_grid.insert(u_grid.begin()+idx, u_row);

  // Also insert the stopping radius if we are in a case where we have
  // to keep track of it
  if (yoverm < 1.0)
    loga_stop_grid.insert(loga_stop_grid.begin()+idx, loga_stop);

  // Check we had to extend the grid further in a; if so, go back and
  // extend any previously-computed grids using the new grid points
  if (loga_new.back() > loga_grid.back()) {

    // Find the last point in the new grid that is also in the old
    // grid
    vector<double>::size_type newptr = 0;
    while (loga_new[newptr] < loga_grid.back()) newptr++;

    // Loop over existing values of U(a) at each x
    for (vector<double>::size_type i=0; i<x_grid.size(); i++) {
      
      // Skip the row we just added
      if (i == idx) continue;

      // Compute u values for this x on the new part of the grid
      double loga_stop;
      vector<double> u_new
	= Ua_solve_tab(x_grid[i], loga_new, loga_stop, newptr,
		       loga_new.size(), u_grid[i].back());

      // Add these new points to the grid
      for (vector<double>::size_type j=1; j<u_new.size(); j++) {
	u_grid[i].push_back(u_new[j]);
      }
    }

    // Extend the grid
    for (newptr++; newptr < loga_new.size(); newptr++)
      loga_grid.push_back(loga_new[newptr]);
  }

  // See if this x required adding any points to the acceleration
  // law grid; if so, go back and re-compute for all previous values
  // of x on the finer grid
  if (loga_new.size() != loga_grid.size()) {
    for (vector<double>::size_type i=0; i<x_grid.size(); i++) {

      // Skip the row we just added
      if (i == idx) continue;

      // March through the old and new grids, seeing where points
      // were added
      vector<double>::size_type oldptr = 0, newptr = 0, uptr = 0;
      while (newptr < loga_new.size()-1) {
	if (loga_new[newptr+1] != loga_grid[oldptr+1]) {
	  // Next grid points don't match, so advance through the new
	  // grid until we find a point that lines up with the old
	  // grid
	  vector<double>::size_type istart = newptr;
	  while (loga_new[newptr+1] < loga_grid[oldptr+1]) newptr++;
	  // Now generate the velocity solution on the new grid
	  // points
	  double loga_stop;
	  vector<double> u_new
	    = Ua_solve_tab(x_grid[i], loga_new, loga_stop, istart,
			   newptr, u_grid[i][uptr]);
	  // Insert new velocities into data; note that the first
	  // point in u_new is a duplicate, so don't insert it
	  u_grid[i].insert(u_grid[i].begin()+uptr+1, u_new.begin()+1,
			   u_new.end());
	  uptr += u_new.size()-1;
	}
	// Advance pointers
	oldptr++;
	newptr++;
	uptr++;
      }
    }
    
    // Store new grid in a
    loga_grid = loga_new;
  }
}

// This is the master routine that builds the interpolation table; it
// features error control in both the x and a dimensions
void pwind_hot::build_table(const vector<double>& x_grid_init) {

  // Generate solutions on starting x grid
  for (vector<double>::size_type i=0; i<x_grid_init.size(); i++)
    add_table_row(x_grid_init[i]);
  
  // Now iterate until accuracy in x direction is adequate; this
  // requires a little care to properly handle the case where u goes
  // to zero at finite a, and that a depends on x
  while (1) {

    // Set up the linear and cubic interpolators on the data; we will
    // use these to estimate the error
    vector<gsl_interp *> interp_lin(loga_grid.size()),
      interp_cube(loga_grid.size());
    vector<gsl_interp_accel *> acc_lin(loga_grid.size()),
      acc_cube(loga_grid.size());
    vector<vector<double> > u(loga_grid.size());   // Transpose of u_grid
    for (vector<double>::size_type i=0; i<loga_grid.size(); i++) {
      u[i].resize(x_grid.size());
      for (vector<double>::size_type j=0; j<x_grid.size(); j++)
	u[i][j] = u_grid[j][i];
      acc_lin[i] = gsl_interp_accel_alloc();
      acc_cube[i] = gsl_interp_accel_alloc();
      if (loga_stop_grid.size() == 0) {
	// Case where the velocity is always finite
	interp_lin[i] = gsl_interp_alloc(gsl_interp_linear, x_grid.size());
	interp_cube[i] = gsl_interp_alloc(gsl_interp_cspline, x_grid.size());
	gsl_interp_init(interp_lin[i], x_grid.data(), u[i].data(),
			x_grid.size());
	gsl_interp_init(interp_cube[i], x_grid.data(), u[i].data(),
			x_grid.size());
      } else {
	// Case where the velocity goes to zero at finite radius
	vector<double>::size_type j;
	for (j=0; j<x_grid.size(); j++)
	  if (u[i][j] < 0.0) break;
	if (j >= gsl_interp_type_min_size(gsl_interp_cspline)) {
	  interp_lin[i] = gsl_interp_alloc(gsl_interp_linear, j);
	  interp_cube[i] = gsl_interp_alloc(gsl_interp_cspline, j);
	  gsl_interp_init(interp_lin[i], x_grid.data(), u[i].data(),
			  j);
	  gsl_interp_init(interp_cube[i], x_grid.data(), u[i].data(),
			  j);
	} else {
	  interp_lin[i] = nullptr;
	  interp_cube[i] = nullptr;
	}
      }
    }
    
    // Loop over the x direction
    vector<vector<double>::size_type> intervals;
    for (vector<double>::size_type j=0; j<x_grid.size()-1; j++) {
      unsigned int ncheck=8;
      for (unsigned int k=1; k<ncheck; k++) {
	double x_pt = (ncheck-k)*x_grid[j]/((float) ncheck) +
	  k * x_grid[j+1]/((float) ncheck);
	// Loop over the log a direction
	for (vector<double>::size_type i=0; i<loga_grid.size(); i++) {
	  // Stop if we've gone beyond the grid
	  if (loga_stop_grid.size() > 0)
	    if (loga_grid[i] > loga_stop_grid[j+1]) break;
	  // Stop if we've gone out so far in a that there were too
	  // few x's that got this far for us to build interpolators
	  if (interp_lin[i] == nullptr) break;
	  // Stop if we've run off the edge for for this interpolator
	  if (j+1 >= interp_lin[i]->size) break;
	  // Check this interval in x and log a
	  double u1 = gsl_interp_eval(interp_lin[i], x_grid.data(),
				      u[i].data(), x_pt, acc_lin[i]);
	  double u2 = gsl_interp_eval(interp_cube[i], x_grid.data(),
				      u[i].data(), x_pt, acc_cube[i]);
	  double abserr = fabs(u1-u2);
	  double relerr = abserr/fabs(u1);
	  if (abserr > interpabs && relerr > interprel) {
	    intervals.push_back(j);
	    break;
	  }
	}
	if (k < ncheck) break;
      }
    }

    // Free the existing interpolators
    for (vector<double>::size_type i=0; i<loga_grid.size(); i++) {
      gsl_interp_free(interp_lin[i]);
      gsl_interp_free(interp_cube[i]);
      gsl_interp_accel_free(acc_lin[i]);
      gsl_interp_accel_free(acc_cube[i]);
    }

    // If we have stopping radii, we also need to ensure that those
    // are well-measured; form linear and spline interpolators of
    // loga_stop versus x, and check them for accuracy
    if (loga_stop_grid.size() > 0) {

      // Set up interpolators for stopping radius
      gsl_interp *interp_lin =
	gsl_interp_alloc(gsl_interp_linear, x_grid.size());
      gsl_interp *interp_cube =
	gsl_interp_alloc(gsl_interp_cspline, x_grid.size());
      gsl_interp_accel *acc_lin = gsl_interp_accel_alloc();
      gsl_interp_accel *acc_cube = gsl_interp_accel_alloc();
      gsl_interp_init(interp_lin, x_grid.data(), loga_stop_grid.data(),
		      x_grid.size());
      gsl_interp_init(interp_cube, x_grid.data(), loga_stop_grid.data(),
		      x_grid.size());

      // Loop and check accuracy
      for (vector<double>::size_type i=0; i<x_grid.size()-1; i++) {
	unsigned int ncheck = 8;
	for (unsigned int j=1; j<ncheck; j++) {
	  double x_pt = (ncheck-j)*x_grid[i]/((float) ncheck) +
	    j*x_grid[i+1]/((float) ncheck);
	  double a_stop_1
	    = exp(gsl_interp_eval(interp_lin, x_grid.data(),
				  loga_stop_grid.data(),
				  x_pt, acc_lin));
	  double a_stop_2
	    = exp(gsl_interp_eval(interp_cube, x_grid.data(),
				  loga_stop_grid.data(),
				  x_pt, acc_lin));
	  double abserr = fabs(a_stop_1 - a_stop_2);
	  // Note: divide by 2 for safety margin -- since our
	  // interpolation as u vs a are only good to a relative
	  // accuracy of relerr, our errors the value of a at which u
	  // goes to zero may not be able to achieve an accuracy
	  // better than 2 * relerr, and demanding a relative error of
	  // relerr may result in an infinite loop
	  double relerr = abserr / a_stop_1 / 2.0;
	  if (abserr > interpabs && relerr > interprel) {
	    intervals.push_back(i);
	    break;
	  }
	}
      }

      // Remove duplicates
      sort(intervals.begin(), intervals.end());
      vector<vector<double>::size_type>::iterator it
	= unique(intervals.begin(), intervals.end());
      intervals.resize(distance(intervals.begin(), it));

      // Free memory
      gsl_interp_free(interp_lin);
      gsl_interp_free(interp_cube);
      gsl_interp_accel_free(acc_lin);
      gsl_interp_accel_free(acc_cube);
    }
    
    // We have now found all the intervals in the x direction where
    // the resolution is inadequate; if there are none, we are done
    if (intervals.size() == 0) break;

    // If we're here, we need to add some table rows by subdividing
    // some of our existing intervals in x
    for (vector<double>::size_type i=0; i<intervals.size(); i++) {
      double x_new = 0.5*(x_grid[i+intervals[i]] +
			  x_grid[i+intervals[i]+1]);
      add_table_row(x_new);
    }
  }
}

// This builds all the interpolators we require from the tabulated
// data
void pwind_hot::build_interpolators() {

  // Build the 1D interpolators for Umin and Umax
  Umin_interp = gsl_spline_alloc(gsl_interp_cspline, loga_grid.size());
  Umax_interp = gsl_spline_alloc(gsl_interp_cspline, loga_grid.size());
  gsl_spline_init(Umin_interp, loga_grid.data(), u_grid.back().data(),
		  loga_grid.size());
  gsl_spline_init(Umax_interp, loga_grid.data(), u_grid.front().data(),
		  loga_grid.size());
  Umin_acc = gsl_interp_accel_alloc();
  Umax_acc = gsl_interp_accel_alloc();

  // Build a 1D interpolator that returns a(u) along the x = xmin
  // trajectory; note that we may need two of these if we are in a
  // fountain solution, since u is not a monotonic function of a in
  // that case
  vector<double> u_grid_tmp, loga_grid_tmp;
  if (loga_stop_grid.size() == 0) {
    // Wind case
    u_grid_tmp.push_back(u_grid.front()[0]);
    loga_grid_tmp.push_back(loga_grid[0]);
    for (vector<double>::size_type i=1; i<u_grid.front().size(); i++) {
      // Note: this if statement is here to prevent crashes in case a
      // user asks for an interpolation accuracy that is too high to be
      // achieved given the accuracy of our ODE solver, and this results
      // in an acceleration law that does not increase over some interval
      if (u_grid.front()[i] > u_grid_tmp.back()) {
	u_grid_tmp.push_back(u_grid.front()[i]);
	loga_grid_tmp.push_back(loga_grid[i]);
      }
    } 
    aMin_interp[0] = gsl_spline_alloc(gsl_interp_cspline,
				      loga_grid_tmp.size());
    gsl_spline_init(aMin_interp[0], u_grid_tmp.data(),
		    loga_grid_tmp.data(), loga_grid_tmp.size());
    aMin_acc[0] = gsl_interp_accel_alloc();
    aMin_interp[1] = nullptr;
    aMin_acc[1] = nullptr;
    aMin_interp_breakpt = numeric_limits<double>::max();
  } else {
    // Fountain case. In this case we build two interpolators, one for
    // the rising and one for the falling part of the acceleration
    // law, and we record the breakpoint between them.
    vector<double>::size_type i=0;
    // First find the global maximum
    umax = 0.0;
    for (vector<double>::size_type j=0; j<loga_grid.size(); j++) {
      if (u_grid[0][j] > umax) {
	umax = u_grid[0][j];
	i = j;
	aMin_interp_breakpt = loga_grid[j];
      }
    }
    // Build interpolator for region where u is rising with a
    u_grid_tmp.push_back(u_grid.front()[0]);
    loga_grid_tmp.push_back(loga_grid[0]);
    for (vector<double>::size_type j=1; j<=i; j++) {
      if (u_grid.front()[j] > u_grid_tmp.back()) {
	u_grid_tmp.push_back(u_grid.front()[j]);
	loga_grid_tmp.push_back(loga_grid[j]);
      }
    }
    aMin_interp[0] = gsl_spline_alloc(gsl_interp_cspline,
				      u_grid_tmp.size());
    gsl_spline_init(aMin_interp[0], u_grid_tmp.data(),
		    loga_grid_tmp.data(), u_grid_tmp.size());
    aMin_acc[0] = gsl_interp_accel_alloc();
    // Now get the second part; note that we have to reverse the order
    // because the GSL wants things in strictly increasing order
    u_grid_tmp.resize(0);
    loga_grid_tmp.resize(0);
    for (vector<double>::size_type j=loga_grid.size()-1; j>=i; j--) {
      if (u_grid.front()[j] < 0.0) continue;
      if (u_grid_tmp.size() > 0)
	if (u_grid.front()[j] <= u_grid_tmp.back()) continue;
      u_grid_tmp.push_back(u_grid.front()[j]);
      loga_grid_tmp.push_back(loga_grid[j]);
    }
    aMin_interp[1] = gsl_spline_alloc(gsl_interp_cspline,
				      u_grid_tmp.size());
    gsl_spline_init(aMin_interp[1], u_grid_tmp.data(),
		    loga_grid_tmp.data(), u_grid_tmp.size());
    aMin_acc[1] = gsl_interp_accel_alloc();
  }

  // If neded, build a 1D interpolator that returns a(u) along the x =
  // xcrit trajectory
  if (dyda(1.0) <= dmda(1.0))
    aMax_interp = nullptr;
  else {
    u_grid_tmp.resize(0);
    loga_grid_tmp.resize(0);
    u_grid_tmp.push_back(u_grid.back()[0]);
    loga_grid_tmp.push_back(loga_grid[0]);
    for (vector<double>::size_type i=1; i<u_grid.back().size(); i++) {
      if (u_grid.back()[i] > u_grid_tmp.back()) {
	u_grid_tmp.push_back(u_grid.back()[i]);
	loga_grid_tmp.push_back(loga_grid[i]);
      }
    } 
    aMax_interp = gsl_spline_alloc(gsl_interp_cspline, loga_grid_tmp.size());
    gsl_spline_init(aMax_interp, u_grid_tmp.data(), loga_grid_tmp.data(),
		    loga_grid_tmp.size());
    aMax_acc = gsl_interp_accel_alloc();    
  }

  // If needed, build 1d interpolators to return the stopping radius
  // as a function of x, and the stopping x as a function of radius
  if (loga_stop_grid.size() > 0) {
    a_stop_interp =
      gsl_spline_alloc(gsl_interp_cspline, loga_stop_grid.size());
    gsl_spline_init(a_stop_interp, x_grid.data(), loga_stop_grid.data(),
		    loga_stop_grid.size());
    a_stop_acc = gsl_interp_accel_alloc();
    vector<double> x_tmp(x_grid.size()), loga_stop_tmp(x_grid.size());
    vector<double>::size_type idx;
    for (idx=0; idx<x_grid.size(); idx++) {
      if (idx < x_grid.size()-1) {
	// Stop if we are going to hit the grid edge, so the maximum
	// log a is not increasing
	if (loga_stop_grid[x_grid.size()-(idx+1)-1] <=
	    loga_stop_grid[x_grid.size()-idx-1]) {
	  break;
	}
      }
      x_tmp[idx] = x_grid[x_grid.size()-idx-1];
      loga_stop_tmp[idx] = loga_stop_grid[x_grid.size()-idx-1];
    }
    // If we stopped early, fill in the last data point, then truncate
    if (idx < x_grid.size()) {
      x_tmp[idx] = x_grid[0];
      loga_stop_tmp[idx] = loga_stop_grid[0];
      x_tmp.resize(idx+1);
      loga_stop_tmp.resize(idx+1);
    }
    x_stop_interp = 
      gsl_spline_alloc(gsl_interp_cspline, loga_stop_tmp.size());
    gsl_spline_init(x_stop_interp, loga_stop_tmp.data(), x_tmp.data(),
		    loga_stop_tmp.size());
    x_stop_acc = gsl_interp_accel_alloc();
  } else {
    a_stop_interp = x_stop_interp = nullptr;
    a_stop_acc = x_stop_acc = nullptr;
  }

  // Build the 2D interpolator required to compute U(x, a); use 0 at
  // (x, a) values that are forbidden, since u -> 0 as we approach
  // these regions
  double *U2flat = new double[x_grid.size()*loga_grid.size()];
  for (vector<double>::size_type i=0; i<x_grid.size(); i++) {
    for (vector<double>::size_type j=0; j<loga_grid.size(); j++) {
      if (u_grid[i][j] > 0)
	U2flat[j*x_grid.size()+i] = SQR(u_grid[i][j]);
      else
	U2flat[j*x_grid.size()+i] = 0.0;
    }
  }
  U2_interp = gsl_spline2d_alloc(gsl_interp2d_bilinear, x_grid.size(),
				 loga_grid.size());
  gsl_spline2d_init(U2_interp, x_grid.data(), loga_grid.data(),
		    U2flat, x_grid.size(), loga_grid.size());
  U2_acc_x = gsl_interp_accel_alloc();
  U2_acc_a = gsl_interp_accel_alloc();
  double u2m = 0.0;
  for (unsigned int i=0; i<x_grid.size()*loga_grid.size(); i++) {
    if (U2flat[i] > u2m) u2m = U2flat[i];
  }
  delete [] U2flat;  // Spline interface makes an internal copy

  // Now build the interpolator to get X(u, a). This is a bit tricky,
  // because our data are not uniformly sampled in u, and not all
  // values of u are allowed at all a. We follow this procedure:
  // (1) Set up a grid in u by taking all the u values in our table
  //     (at all x and a) and ordering them.
  // (2) At each log a, make a 1D interpolation function that returns
  //     X(u).
  // (3) At each (u, a) point:
  //   (3a) If U_min(a) < u < U_max(a), find X(u, a) using the 1D
  //        interpolator build in step 2.
  //   (3b) If u < U_min(a), set X(u, a) = xcrit
  //   (3c) If u > U_max(a), set X(u, a) = xmin
  // (5) Construct 2D interpolation function for X(u, a)

  // Step (1): set up a grid in u, consisting of sorted, unique u values
  vector<double> ux_grid;
  for (vector<double>::size_type i=0; i<x_grid.size(); i++) {
    for (vector<double>::size_type j=0; j<loga_grid.size(); j++) {
      ux_grid.push_back(u_grid[i][j]);
    }
  }
  sort(ux_grid.begin(), ux_grid.end());
  vector<double>::iterator it = unique(ux_grid.begin(), ux_grid.end());
  ux_grid.resize(distance(ux_grid.begin(), it));

  // Steps (2) - (3): loop over log a
  double *Xflat = new double[ux_grid.size()*loga_grid.size()];
  for (vector<double>::size_type j=0; j<loga_grid.size(); j++) {

    // Step (2): build 1D interpolator to return X(u, a) at fixed a;
    // note that u is decreasing with x, so that we need to traverse
    // the set of (u, x) values backwards to get something in
    // increasing order that will make the GSL happy
    vector<double> utmp, xtmp;
    for (long i=x_grid.size()-1; i>=0; i--) {
      if (u_grid[i][j] < 0.0) continue;
      if (utmp.size() > 0)
	if (u_grid[i][j] <= utmp.back()) continue;
      utmp.push_back(u_grid[i][j]);
      xtmp.push_back(x_grid[i]);
    }
    // If at this a there is a finite x such that x = X(0, a) (i.e.,
    // if the velocity goes to zero at finite a), then we need to add
    // an extra grid point giving that solution
    if (u_grid.back()[j] < 0.0) {
      utmp.insert(utmp.begin(), 0.0);
      xtmp.insert(xtmp.begin(), x_stop(exp(loga_grid[j])));
    }
    gsl_interp *xu_interp;
    gsl_interp_accel *xu_acc;
    // Need to check how many points we have; for a fountain flow, at
    // high a may only have one point
    if (utmp.size() >= gsl_interp_type_min_size(gsl_interp_akima))
      xu_interp = gsl_interp_alloc(gsl_interp_akima, utmp.size());
    else if (utmp.size() >= gsl_interp_type_min_size(gsl_interp_linear))
      xu_interp = gsl_interp_alloc(gsl_interp_linear, utmp.size());
    else
      xu_interp = nullptr;
    if (xu_interp != nullptr) {
      gsl_interp_init(xu_interp, utmp.data(), xtmp.data(), utmp.size());
      xu_acc = gsl_interp_accel_alloc();
    } else {
      xu_acc = nullptr;
    }

    // Step (3): construct a grid of X(u, a) values
    for (vector<double>::size_type i=0; i<ux_grid.size(); i++) {
      double u = ux_grid[i];

      // Check which case we're in
      if (u <= utmp.front()) {
	// Case (3b)
	Xflat[j*ux_grid.size()+i] = xcrit;
      } else if (u >= utmp.back() || xu_interp == nullptr) {
	// Case (3c)
	Xflat[j*ux_grid.size()+i] = x_grid[0];
      } else {
	// Case (3a)
	Xflat[j*ux_grid.size()+i] =
	  gsl_interp_eval(xu_interp, utmp.data(), xtmp.data(),
			  u, xu_acc);
      }
    }

    // Free
    if (xu_interp != nullptr) {
      gsl_interp_free(xu_interp);
      gsl_interp_accel_free(xu_acc);
    }
  }

  // Step (4): build a 2D interpolator for X(u, a) from the data; use
  // bilinear interpolation for safety; free unneeded memory
  X_interp = gsl_spline2d_alloc(gsl_interp2d_bilinear, ux_grid.size(),
				loga_grid.size());
  gsl_spline2d_init(X_interp, ux_grid.data(), loga_grid.data(),
		    Xflat, ux_grid.size(), loga_grid.size());    
  X_acc_u = gsl_interp_accel_alloc();
  X_acc_a = gsl_interp_accel_alloc();
  delete [] Xflat;
}

// This initialization routine computes the wind acceleration law,
// tabulates it, and makes the interpolators required to use it. We do
// this in a separate routine and not in the constructor because we
// require access to the wind expansion factor and potential, which
// are not defined until the derived classes are instantiatied, which
// happens after pwind_hot's constructor has run. Thus the calling
// order must be pwind_hot::pwind_hot, followed by
// pwind_hot_XY::pwind_hot_XY, and only then can we call
// pwind_hot::init().
void pwind_hot::init(const double amax_init) {

  // Establish the initial grid in radius; we start with a grid spaced
  // by factors of 2 in radius; set the first grid point to not
  // exaclty 1
  loga_grid.push_back(0.0);
  while (loga_grid.back() < log(amax_init))
    loga_grid.push_back(loga_grid.back()+log(2.0));
  loga_grid[0] = 1.0e-8;

  // Figure out how low in x we need to go to capture enough of the mass
  struct xmin_params xpar;
  xpar.abstol = interpabs;
  xpar.sx = sx;
  xpar.zeta_M = zeta_M;
  gsl_function F;
  F.function = &xmin_resid;
  F.params = &xpar;

  // Specify that we want to use the Brent solver, and allocate its
  // workspace
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  // Set brackets
  double xlo = -100, xhi = xcrit;
  gsl_root_fsolver_set(s, &F, xlo, xhi);

  // Iterate
  double xmin = 0.0;
  for (int iter=0, status=GSL_CONTINUE;
       (iter<100) && (status=GSL_CONTINUE);
       iter++) {
    status = gsl_root_fsolver_iterate(s);
    xmin = gsl_root_fsolver_root(s);
    xlo = gsl_root_fsolver_x_lower(s);
    xhi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(xlo, xhi, epsabs, epsrel);
  }
  gsl_root_fsolver_free(s);

  // Establish the initial grid in x; start with factors of 2 in
  // Sigma, corresponding to steps of log(2) in x
  vector<double> x_grid_init;
  x_grid_init.push_back(xcrit);
  while (x_grid_init[0] > xmin) {
    x_grid_init.insert(x_grid_init.begin(), 1, x_grid_init[0]-log(2.0));
  }

  // Initialise the ODE integrator
  par.w = this;
  sys = { &duda, nullptr, 1, &par };
  ode_drv =
    gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
  				  1.0e-6, epsrel, epsabs);

  // Build the interpolation table
  build_table(x_grid_init);

  // Free the integrator
  gsl_odeiv2_driver_free(ode_drv);

  // Build the interpolators
  build_interpolators();

  // Set limits on u; note that we only handle the case of a
  // monotonically-increasing velocity here; the case where the
  // velocity peaks at some finite radius is handled inside
  // build_interpolators
  if (a_stop_interp == nullptr)
    umax = u_grid[0].back();
}

////////////////////////////////////////////////////////////////////////
// Limits on the range of the wind
////////////////////////////////////////////////////////////////////////

inline double
pwind_hot::amax(const double x) const {
  if (a_stop_interp == nullptr) {
    // Wind reaches infinity
    return numeric_limits<double>::max();
  } else {
    // Wind turns around at finite radius
    return exp(gsl_spline_eval(a_stop_interp, x, a_stop_acc));
  }
}

inline double
pwind_hot::x_stop(const double a) const {
  if (x_stop_interp == nullptr) {
    // Any material launched can reach any a
    return xcrit;
  } else {
    double loga = log(a);
    // Avoid roundoff issues
    if (loga > loga_stop_grid[0]) loga = loga_stop_grid[0];
    return gsl_spline_eval(x_stop_interp, loga, x_stop_acc);
  }
}

// Get limits on a at a given u, varpi. This version of the routine
// works for all cases except isothermal potental, constant area,
// because that case is a fountain rather than a wind. We override
// this function for that case below.
vector<double>
pwind_hot::alimits(const double u, const double varpi,
		   const double varpi_t) const {
  vector<double> alim;
  double uabs = fabs(u);

  // Check if the wind ever reaches the required velocity; if not,
  // return an empty vector
  if (uabs >= umax) return alim;
  
  // If we're here, the wind does reach the requested line of sight
  // velocity; we need to find out where numerically
  alim.resize(2);

  // Decide what to to based on varpi^2 + varpi_t^2
  double vp2 = SQR(varpi) + SQR(varpi_t);
  if (vp2 == 0.0) {

    // In this case, we can just use our tabulated aMin and aMax
    // functions directly, without needing to solve numerically. The
    // only trick is that we have to catch the case where the input
    // value of u is so small that it is below the smallest velocity
    // in our table, and thus outside the range of validity of the
    // interpolator. In this case, we use the analytic series
    // expansions for the behavior of the solution near a = 1.

    // Limit along x = xmin trajectory
    if (uabs > u_grid[0][0])
      // Interpolation
      alim[0] = gsl_spline_eval(aMin_interp[0], uabs, aMin_acc[0]);
    else
      // Analytic limit
      alim[0] = 1.0 + SQR(u)/((Gamma*exp(-x_grid[0])-1.0));

    // Limit along x = xcrit trajectory
    if (uabs > uMin_inf)
      // Velocity exceeds velocity at infinity, so the upper limit is
      // infinite
      alim[1] = numeric_limits<double>::max();
    else {
      // Velocity is lower than velocity along x = xcrit at infinity,
      // so there is a finite solution
      if (uabs > u_grid.back()[0])
	// Interpolation
	alim[1] = gsl_spline_eval(aMax_interp, uabs, aMax_acc);
      else
	// Analytic limit
	alim[1] = 1.0 + 2.0*uabs*uh *
	  (-1.0 + sqrt(1.0+2.0*SQR(uh)*(dyda(1.0)-dmda(1.0))));
    }

    // Return
    return alim;
    
  } else {

    // For vp2 != 0, we need to find the required velocity
    // numerically; as for vp2 == 0, be careful to handle the case
    // of very small u using the analytic solution. Note that if vp2
    // > a_grid[0]^2, we are guaranteed to find a solution inside the
    // grid, since u = 0 at a = vp2

    // Check if input velocity is within the range achieved on our
    // interpolation grid
    double amn = exp(loga_grid[0]), amx = exp(loga_grid.back());
    double umn, umx;
    if (vp2 >= SQR(amn)) umn = 0.0;
    else umn = u_grid[0][0] * sqrt(1.0-vp2/SQR(amn));
    if (vp2 >= SQR(amx)) umx = 0.0;
    else umx = umax * sqrt(1.0-vp2/SQR(amx));

    // Handle cases
    if (uabs <= umn) {
      
      // Input velocity is below the lowest achieved on our
      // interpolation grid, so solve for a using the analytic series
      // expansion for u(a) as a -> 1; equation in this case is a cubic
      double gex = Gamma*exp(-x_grid[0]) - 1.0;
      double roots[3];
      int nroots =
	gsl_poly_solve_cubic(-(SQR(u)+gex)/gex, -vp2, vp2,
			     roots, roots+1, roots+2);
      for (int i=0; i<nroots; i++) {
	if (roots[i] > 0.0) {
	  alim[0] = roots[i];
	  break;
	}
      }
    } else if (uabs >= umx) {

      // Input velocity is above the highest achieved on our
      // interpolation grid (but below the highest radial velocity
      // possible, or we would have caught it above). In this case
      // solve analytically assuming that ur = const
      alim[0] = sqrt(vp2) / sqrt(1.0 - SQR(u/umax));
      
    } else {

      // Desired velocity is achieved within the range of our
      // interpolation function, so solve numerically

      // Set pointers for GSL
      alim_params aparams;
      aparams.u2 = SQR(u);
      aparams.vp2 = vp2;
      aparams.x = x_grid[0];
      aparams.w = this;
      gsl_function F;
      F.function = &alim_func;
      F.params = &aparams;

      // Specify that we want to use the Brent solver, and allocate its
      // workspace
      const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
      gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

      // Set brackets
      double loga_lo = log(max(amn, sqrt(vp2)));
      double loga_hi = loga_grid.back();
      gsl_root_fsolver_set(s, &F, loga_lo, loga_hi);

      // Loop
      for (int iter=0, status=GSL_CONTINUE;
	   (iter<100) && (status==GSL_CONTINUE);
	   iter++) {
	status = gsl_root_fsolver_iterate(s);
	alim[0] = exp(gsl_root_fsolver_root(s));
	loga_lo = gsl_root_fsolver_x_lower(s);
	loga_hi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(loga_lo, loga_hi, epsabs, epsrel);
      }

      // Free workspace
      gsl_root_fsolver_free(s);
    }

    // Now we need to handle the upper limit on a, which is reached
    // along the x = xcrit trajectory

    // First check if there is a solution at all; one exists only if
    // the velocity is less than the velocity at infinity
    if (uabs > uMin_inf)
      alim[1] = numeric_limits<double>::max();
    else {

      // Solution does exist. Next check if it is found within our
      // interpolation grid
      if (vp2 >= SQR(amn)) umn = 0.0;
      else umn = u_grid.back()[0] * sqrt(1.0-vp2/SQR(amn));
      if (vp2 >= SQR(amx)) umx = 0.0;
      else umx = uMin_inf * sqrt(1.0-vp2/SQR(amx));
      
      // Handle cases
      if (uabs < umn) {

	// Solution is below interpolation grid, so solve using the
	// analytic solution along x = xcrit. The resulting equation is
	// a quartic, for which we find the root numerically using the
	// GSL root solver.
	alim_quart_params aparams;
	aparams.u2 = SQR(u);
	aparams.vp2 = vp2;
	aparams.f2 =
	  SQR( (-1.0+sqrt(1.0+2.0*SQR(uh)*(dyda(1.0)-dmda(1.0)))) /
	       (2.0*uh) );
	gsl_function F;
	F.function = &alim_quart_func;
	F.params = &aparams;
	const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
	double a_lo = 1.0;
	double a_hi = amn;
	gsl_root_fsolver_set(s, &F, a_lo, a_hi);
	for (int iter=0, status=GSL_CONTINUE;
	     (iter<100) && (status==GSL_CONTINUE);
	     iter++) {
	  status = gsl_root_fsolver_iterate(s);
	  alim[1] = exp(gsl_root_fsolver_root(s));
	  a_lo = gsl_root_fsolver_x_lower(s);
	  a_hi = gsl_root_fsolver_x_upper(s);
	  status = gsl_root_test_interval(a_lo, a_hi, epsabs, epsrel);
	}
	gsl_root_fsolver_free(s);
	
      } else if (uabs > umx) {

	// Velocity is above highest found in our interpolation table;
	// solve analytically, assuming that radial velocity =
	// constant above the upper limit of our table
	alim[1] = sqrt(vp2) / sqrt(1.0 - SQR(u/uMin_inf));

      } else {

	// Input velocity can be matched within our interpolation
	// grid, so find the solution numerically

	// Set pointers for GSL
	alim_params aparams;
	aparams.u2 = SQR(u);
	aparams.vp2 = vp2;
	aparams.x = xcrit;
	aparams.w = this;
	gsl_function F;
	F.function = &alim_func;
	F.params = &aparams;

	// Specify that we want to use the Brent solver, and allocate its
	// workspace
	const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

	// Set brackets
	double loga_lo = log(max(amn, sqrt(vp2)));
	double loga_hi = loga_grid.back();
	gsl_root_fsolver_set(s, &F, loga_lo, loga_hi);

	// Loop
	for (int iter=0, status=GSL_CONTINUE;
	     (iter<100) && (status==GSL_CONTINUE);
	     iter++) {
	  status = gsl_root_fsolver_iterate(s);
	  alim[1] = exp(gsl_root_fsolver_root(s));
	  loga_lo = gsl_root_fsolver_x_lower(s);
	  loga_hi = gsl_root_fsolver_x_upper(s);
	  status = gsl_root_test_interval(loga_lo, loga_hi, epsabs, epsrel);
	}

	// Free workspace
	gsl_root_fsolver_free(s);
      }
    }
      
    // Return
    return alim;
  }
}

////////////////////////////////////////////////////////////////////////
// Winds specialised by expansion and potential
////////////////////////////////////////////////////////////////////////

// Point potential, constant area
pwind_hot_pa::pwind_hot_pa(const double Gamma_, const double mach_,
			   const double uh_,
			   const pwind_geom* geom_,
			   const double epsabs_,
			   const double epsrel_,
			   const double interpabs_,
			   const double interprel_,
			   const double fcrit_,
			   const double amax_grid_,
			   const double amax_init) :
  pwind_hot(Gamma_, mach_, uh_, 1.0,
	    static_cast<const pwind_potential *>(&pwind_potential_point),
	    static_cast<const pwind_expansion *>(&pwind_expansion_area),
	    geom_, epsabs_,
	    epsrel_, interpabs_, interprel_, fcrit_, amax_grid_) {
  init(amax_init);
  amax_abs = numeric_limits<double>::max(); // Gas reaches infinity
  uMin_inf = 0.0; // Material with x = xcrit has 0 velocity as a->infinity
}
inline vector<double>
pwind_hot_pa::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = getXMin();
  ret[1] = xcrit;
  return ret;
}

// Point potential, intermediate expansion
pwind_hot_pi::pwind_hot_pi(const double Gamma_, const double mach_,
			   const double uh_,
			   const pwind_geom* geom_,
			   const double epsabs_,
			   const double epsrel_,
			   const double interpabs_,
			   const double interprel_,
			   const double fcrit_,
			   const double amax_grid_,
			   const double amax_init) :
  pwind_hot(Gamma_, mach_, uh_, numeric_limits<double>::max(),
	    static_cast<const pwind_potential *>(&pwind_potential_point),
	    static_cast<const pwind_expansion *>
	    (&pwind_expansion_intermediate),
	    geom_, epsabs_,
	    epsrel_, interpabs_, interprel_, fcrit_, amax_grid_) {
  init(amax_init);
  amax_abs = numeric_limits<double>::max();  // Gas reaches infinity
  uMin_inf = u_grid.back().back();
}
inline vector<double>
pwind_hot_pi::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = getXMin();
  ret[1] = xcrit;
  return ret;
}

// Point potential, constant solid angle
pwind_hot_ps::pwind_hot_ps(const double Gamma_, const double mach_,
			   const double uh_,
			   const pwind_geom* geom_,
			   const double epsabs_,
			   const double epsrel_,
			   const double interpabs_,
			   const double interprel_,
			   const double fcrit_,
			   const double amax_grid_,
			   const double amax_init) :
  pwind_hot(Gamma_, mach_, uh_, numeric_limits<double>::max(),
	    static_cast<const pwind_potential *>(&pwind_potential_point),
	    static_cast<const pwind_expansion *>(&pwind_expansion_solid_angle),
	    geom_, epsabs_,
	    epsrel_, interpabs_, interprel_, fcrit_, amax_grid_) {
  init(amax_init);
  amax_abs = numeric_limits<double>::max();  // Gas reaches infinity
  uMin_inf = u_grid.back().back();
}
inline vector<double>
pwind_hot_ps::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = getXMin();
  ret[1] = xcrit;
  return ret;
}


// Isothermal potential, constant area
pwind_hot_ia::pwind_hot_ia(const double Gamma_, const double mach_,
			   const double uh_,
			   const pwind_geom* geom_,
			   const double epsabs_,
			   const double epsrel_,
			   const double interpabs_,
			   const double interprel_,
			   const double fcrit_,
			   const double amax_grid_,
			   const double amax_init) :
  pwind_hot(Gamma_, mach_, uh_, 0.0,
	    static_cast<const pwind_potential *>(&pwind_potential_isothermal),
	    static_cast<const pwind_expansion *>(&pwind_expansion_area),
	    geom_, epsabs_,
	    epsrel_, interpabs_, interprel_, fcrit_, amax_grid_) {
  init(amax_init);
  amax_abs = exp(loga_grid.back());  // Maximum radius reached
}
inline vector<double>
pwind_hot_ia::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = getXMin();
  ret[1] = x_stop(a);
  return ret;
}

// Version of alimits specialised to this case, which is different
// from all the others because it is a fountain and not a wind. The
// limits on a occur along the maximum velocity trajectory.
vector<double>
pwind_hot_ia::alimits(const double u, const double varpi,
		      const double varpi_t) const {
  vector<double> alim;
  double uabs = fabs(u);

  // Decide what to to based on varpi
  double vp2 = SQR(varpi) + SQR(varpi_t);
  if (vp2 == 0.0) {

    // See if solution exists at all; return if not
    if (uabs >= umax) return alim;

    // In this case we can just read the answers directly from our
    // interpolators; the only trick is what to do if we're outside
    // their limits, because the input velocity is below the smallest
    // we've recorded. If this happens, we handle the small a side
    // using the series expansion, and the large a side just by
    // returning the radius of the smallest u we have
    alim.resize(2);
    if (uabs < aMin_interp[0]->interp->xmin)
      alim[0] = 1.0 + SQR(u)/((Gamma*exp(-x_grid[0])-1.0));
    else
      alim[0] = gsl_spline_eval(aMin_interp[0], uabs, aMin_acc[0]);
    alim[1] = gsl_spline_eval(aMin_interp[1],
			      max(uabs, aMin_interp[1]->interp->xmin),
			      aMin_acc[1]);
    return alim;

  } else {

    // In this case we must solve numerically. This has to be done in
    // two steps. First, we find the maximum line of sight velocity
    // achieved for the input value of vp2; this determines the
    // number of solutions, and, if they exist, the brackets in a
    // where they are confined. Second, if solutions exist, we find
    // them.

    // Step 0: make sure vp2 isn't so large as to be off our grid
    // entirely
    if (0.5*log(vp2) > loga_grid.back()) return alim;
    else alim.resize(2);

    // Step 1: find the maximum
    struct u_max_params params;
    params.vp2 = vp2;
    params.x = x_grid[0];
    params.w = this;
    gsl_function F;
    F.function = &u_max_func;
    F.params = &params;
    const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);
    double loga_lo = max(loga_grid[0], 0.5*log(vp2));
    double loga_hi = loga_grid.back();
    double loga_max;
    gsl_min_fminimizer_set(s, &F, 0.5*(loga_lo+loga_hi),
			   loga_lo, loga_hi);
    for (int iter=0; iter<100; iter++) {
      int status = gsl_min_fminimizer_iterate(s);
      loga_max = gsl_min_fminimizer_x_minimum(s);
      loga_lo = gsl_min_fminimizer_x_lower(s);
      loga_hi = gsl_min_fminimizer_x_upper(s);
      status = gsl_min_test_interval(loga_lo, loga_hi,
				     epsabs, epsrel);
      if (status != GSL_CONTINUE) break;
    }
    gsl_min_fminimizer_free(s);
    double ulos_max = U(x_grid[0], exp(loga_max)) *
      sqrt(1.0-vp2/SQR(exp(loga_max)));

    // Step 2: we have found the maximum velocity; now check if a
    // solution exists, and, if so, find it numerically
    if (uabs >= ulos_max) return alim;
    struct alim_params apar;
    apar.u2 = SQR(u);
    apar.x = x_grid[0];
    apar.vp2 = vp2;
    apar.w = this;
    gsl_function F1;
    F1.function = &alim_func;
    F1.params = &apar;
    const gsl_root_fsolver_type *T1 = gsl_root_fsolver_brent;
    gsl_root_fsolver *s1 = gsl_root_fsolver_alloc(T1);

    // Step 2a: find the root below the maximum
    loga_lo = max(loga_grid[0], 0.5*log(vp2));
    loga_hi = loga_max;
    gsl_root_fsolver_set(s1, &F1, loga_lo, loga_hi);
    for (int iter=0, status=GSL_CONTINUE;
	 (iter<100) && (status==GSL_CONTINUE);
	 iter++) {
      status = gsl_root_fsolver_iterate(s1);
      alim[0] = exp(gsl_root_fsolver_root(s1));
      loga_lo = gsl_root_fsolver_x_lower(s1);
      loga_hi = gsl_root_fsolver_x_upper(s1);
      status = gsl_root_test_interval(loga_lo, loga_hi, epsabs, epsrel);
    }

    // Step 2b: find the root above the maximum; note that this root
    // might be beyond the edge of our grid, in which case we will set
    // the limit to where our grid is truncated
    loga_lo = loga_max;
    loga_hi = loga_grid.back();
    if (F1.function(loga_lo, F1.params) *
	F1.function(loga_hi, F1.params) < 0) {
      gsl_root_fsolver_set(s1, &F1, loga_lo, loga_hi);
      for (int iter=0, status=GSL_CONTINUE;
	   (iter<100) && (status==GSL_CONTINUE);
	   iter++) {
	status = gsl_root_fsolver_iterate(s1);
	alim[1] = exp(gsl_root_fsolver_root(s1));
	loga_lo = gsl_root_fsolver_x_lower(s1);
	loga_hi = gsl_root_fsolver_x_upper(s1);
	status = gsl_root_test_interval(loga_lo, loga_hi, epsabs, epsrel);
      }
    } else {
      alim[1] = exp(loga_hi);
    }
    gsl_root_fsolver_free(s1);

    // Return
    return alim;
  }
}

// Isothermal potential, intermediate expansion
pwind_hot_ii::pwind_hot_ii(const double Gamma_, const double mach_,
			   const double uh_,
			   const pwind_geom* geom_,
			   const double epsabs_,
			   const double epsrel_,
			   const double interpabs_,
			   const double interprel_,
			   const double fcrit_,
			   const double amax_grid_,
			   const double amax_init) :
  pwind_hot(Gamma_, mach_, uh_, 1.0,
	    static_cast<const pwind_potential *>(&pwind_potential_isothermal),
	    static_cast<const pwind_expansion *>(&pwind_expansion_intermediate),
	    geom_, epsabs_,
	    epsrel_, interpabs_, interprel_, fcrit_, amax_grid_) {
  init(amax_init);
  amax_abs = numeric_limits<double>::max();  // Gas reaches infinity
  uMin_inf = u_grid.back().back();
}
inline vector<double>
pwind_hot_ii::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = getXMin();
  ret[1] = xcrit;
  return ret;
}

// Isothermal potential, constant solid angle
pwind_hot_is::pwind_hot_is(const double Gamma_, const double mach_,
			   const double uh_,
			   const pwind_geom* geom_,
			   const double epsabs_,
			   const double epsrel_,
			   const double interpabs_,
			   const double interprel_,
			   const double fcrit_,
			   const double amax_grid_,
			   const double amax_init) :
  pwind_hot(Gamma_, mach_, uh_, numeric_limits<double>::max(),
	    static_cast<const pwind_potential *>(&pwind_potential_isothermal),
	    static_cast<const pwind_expansion *>(&pwind_expansion_solid_angle),
	    geom_, epsabs_,
	    epsrel_, interpabs_, interprel_, fcrit_, amax_grid_) {
  init(amax_init);
  amax_abs = numeric_limits<double>::max();  // Gas reaches infinity
  uMin_inf = u_grid.back().back();
}
inline vector<double>
pwind_hot_is::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = getXMin();
  ret[1] = xcrit;
  return ret;
}
