#include "pwind_hot.H"
#include "pwind_util.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <unistd.h>
#if __cplusplus >= 201703L
#   include <filesystem>
#endif

using namespace std;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////

pwind_hot::pwind_hot(const double Gamma_,
		     const double mach_,
		     const double uh_,
		     const pwind_potential *potential_,
		     const pwind_expansion *expansion_,
		     const pwind_geom *geom_,
		     const double epsabs_,
		     const double epsrel_,
		     const double fcrit_,
		     const double jsp_,
		     hot_wind_table *full_tab_) :
  pwind(Gamma_, mach_, potential_, expansion_, geom_,
	epsabs_, epsrel_, fcrit_, jsp_),
  uh(uh_),
  full_tab(full_tab_)
{
  // Get xcrit, zeta_M, zeta_A
  xcrit = log(Gamma_);
  zeta_M = zetaM(xcrit + log(fcrit_), sx);
  zeta_A = zetaA(xcrit + log(fcrit_), sx);

  // If we were not given a pointer to the data, try to read off disk
  if (!full_tab) {
  
    // Get the name of the directory where we should go looking for data
    char *desp_home = getenv("DESPOTIC_HOME");
    string despotic_home;
    if (desp_home) {
      despotic_home = string(desp_home);
    } else {
      char cwd[1000];
      despotic_home = string(getcwd(cwd, 1000));
    }

    // Read
    full_tab = read_table(despotic_home, expansion->yidx(), potential->midx());

    // Flag that we own the table
    manage_table = true;
    
  } else {

    // Flag that we do not manage the data table
    manage_table = false;
    
  }

  // Make sure input uh is in bounds, and bail out if not
  if (uh < full_tab->uh.front() || uh > full_tab->uh.back()) {
    cerr << "pwind::hot: input uh out of bounds for tabulated data; "
	 << "table range is uh = " << full_tab->uh.front()
	 << " - " << full_tab->uh.back()
	 << ", requested uh = " << uh
	 << endl;
    exit(1);
  }

  // Allocate space in the interpolated data
  int nuh = full_tab->uh.size();
  int ngex = full_tab->uh_data[0].umax.size();
  int nu = full_tab->uh_data[0].q.size() / ngex;
  int nq = full_tab->uh_data[0].umin.size();
  tab.umax.resize(ngex);
  tab.q_umax.resize(ngex);
  tab.q_stop.resize(ngex);
  tab.q.resize(ngex * nu);
  tab.umin.resize(nq);
  tab.loggex.resize(nu * nq);

  // Interpolate in uh to our target value
  int idx;
  for (idx = 0; idx < full_tab->uh.size()-1; idx++) {
    if (uh <= full_tab->uh[idx+1]) break;
  }
  double w = (uh - full_tab->uh[idx]) /
    (full_tab->uh[idx+1] - full_tab->uh[idx]);
  for (int i=0; i<ngex; i++) {
    tab.umax[i] = (1.0 - w) * full_tab->uh_data[idx].umax[i] +
      w * full_tab->uh_data[idx+1].umax[i];
    tab.q_stop[i] = (1.0 - w) * full_tab->uh_data[idx].q_stop[i] +
      w * full_tab->uh_data[idx+1].q_stop[i];
    if (isFountain()) {
      tab.q_umax[i] = (1.0 - w) * full_tab->uh_data[idx].q_umax[i] +
	w * full_tab->uh_data[idx+1].q_umax[i];
    }
  }
  for (int i=0; i<ngex*nu; i++) {
    tab.q[i] = (1.0 - w) * full_tab->uh_data[idx].q[i] +
      w * full_tab->uh_data[idx+1].q[i];
  }
  if (hasUmin()) {
    for (int i=0; i<nq; i++) {
      tab.umin[i] = (1.0 - w) * full_tab->uh_data[idx].umin[i] +
	w * full_tab->uh_data[idx+1].umin[i];
    }
  }
  for (int i=0; i<nu*nq; i++) {
    tab.loggex[i] = (1.0 - w) * full_tab->uh_data[idx].loggex[i] +
      w * full_tab->uh_data[idx+1].loggex[i];
  }

  // Set the maximum values of a and u
  amax_abs = pow(10., tab.q_stop[ngex-1]) + 1.0;
  umax = tab.umax[ngex-1];  
}

////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////
pwind_hot::~pwind_hot() {
  if (manage_table) free_table(full_tab);
}


////////////////////////////////////////////////////////////////////////
// Method to read the hot wind table
////////////////////////////////////////////////////////////////////////
hot_wind_table *
pwind_hot::read_table(const string &dirname,
		      const int yidx,
		      const int midx) {

  // Allocate the object to hold output
  hot_wind_table *dat = new hot_wind_table;

  // Build name of index file
  stringstream ss;
#if __cplusplus >= 201703L
  filesystem::path p(dirname);
  ss << "table_params_y" << yidx << "_m" << midx << ".txt";
  string idx_file = (p / ss.str()).string();
#else
  ss << dirname << "/table_params_y" << yidx << "_m" << midx << ".txt";
  string idx_file = ss.str();
#endif

  // Open the index file; if we encounter an error, just return null
  ifstream idx_if(idx_file, ios::in);
  if (!idx_if) {
    delete dat;
    return nullptr;
  }
  
  // Read the index file
  int nuh;
  idx_if >> nuh >> dat->nu >> dat->nq >> dat->ngex
	 >> dat->qmin >> dat->qmax >> dat->loggex_max;
  dat->uh.resize(nuh);
  dat->dq = (dat->qmax - dat->qmin) / (dat->nq-1);
  dat->du = 1.0 / (dat->nu-1);
  dat->dlg = dat->loggex_max / (dat->ngex-1);
  for (int i=0; i<nuh; i++) idx_if >> dat->uh[i];
  idx_if.close();

  // Allocate memory for table entries
  dat->uh_data.resize(nuh);
  for (int i=0; i<nuh; i++) {
    dat->uh_data[i].umax.resize(dat->ngex);
    dat->uh_data[i].q_umax.resize(dat->ngex);
    dat->uh_data[i].q_stop.resize(dat->ngex);
    dat->uh_data[i].umin.resize(dat->nq);
    if (yidx != 0 || midx != 1) 
      dat->uh_data[i].q.resize(dat->ngex * dat->nu);
    else
      dat->uh_data[i].q.resize(2 * dat->ngex * dat->nu);
    dat->uh_data[i].loggex.resize(dat->nu * dat->nq);
  }

  // Now begin to read the data files; first file: qtab_u
  ss.str("");
#if __cplusplus >= 201703L
  ss << "qtab_u_y" << yidx << "_m" << midx << ".bin";
  p /= ss.str();
  string fname = p.string();
#else
  ss << dirname << "/qtab_u_y" << yidx << "_m" << midx << ".bin";
  string fname = ss.str();
#endif
  ifstream qtab_u_if(fname, ios::in | ios::binary);
  if (!qtab_u_if) {
    delete dat;
    return nullptr;
  }

  // Read the qtab_u file
  for (int i=0; i<nuh; i++) {
    // Burn first ngex entries, because these list gex values; we
    // can construct these ourselves later as needed
    qtab_u_if.seekg(dat->ngex * sizeof(double), qtab_u_if.cur);
    qtab_u_if.read((char *) (dat->uh_data[i].umax.data()),
		   dat->ngex * sizeof(double));
    qtab_u_if.read((char *) (dat->uh_data[i].q_umax.data()),
		   dat->ngex * sizeof(double));
    qtab_u_if.read((char *) (dat->uh_data[i].q_stop.data()),
		   dat->ngex * sizeof(double));
  }
  qtab_u_if.close();

  // Next data file: qtab_q
  ss.str("");
#if __cplusplus >= 201703L
  ss << "qtab_q_y" << yidx << "_m" << midx << ".bin";
  fname = (p / ss.str()).string();
#else
  ss << dirname << "/qtab_q_y" << yidx << "_m" << midx << ".bin";
  fname = ss.str();
#endif
  ifstream qtab_q_if(fname, ios::in | ios::binary);
  if (!qtab_q_if) {
    delete dat;
    return nullptr;
  }
  for (int i=0; i<nuh; i++) {
    qtab_q_if.read((char *) (dat->uh_data[i].q.data()),
		   dat->uh_data[i].q.size() * sizeof(double));
  }
  qtab_q_if.close();

  // Next data file: gextab_q
  ss.str("");
#if __cplusplus >= 201703L
  ss << "gextab_q_y" << yidx << "_m" << midx << ".bin";
  fname = (p / ss.str()).string();
#else
  ss << dirname << "/gextab_q_y" << yidx << "_m" << midx << ".bin";
  fname = ss.str();
#endif
  ifstream gextab_q_if(fname, ios::in | ios::binary);
  if (!gextab_q_if) {
    delete dat;
    return nullptr;
  }
  for (int i=0; i<nuh; i++) {
    // Burn first nq entries, because these list q values; we
    // can construct these ourselves later as needed
    gextab_q_if.seekg(dat->nq * sizeof(double), gextab_q_if.cur);
    gextab_q_if.read((char *) (dat->uh_data[i].umin.data()),
		     dat->nq * sizeof(double));
  }
  gextab_q_if.close();

  // Last data file: gextab_gex
  ss.str("");
#if __cplusplus >= 201703L
  ss << "gextab_u_y" << yidx << "_m" << midx << ".bin";
  fname = (p / ss.str()).string();
#else
  ss << dirname << "/gextab_gex_y" << yidx << "_m" << midx << ".bin";
  fname = ss.str();
#endif
  ifstream gextab_gex_if(fname, ios::in | ios::binary);
  if (!gextab_gex_if) {
    delete dat;
    return nullptr;
  }
  for (int i=0; i<nuh; i++) {
    gextab_gex_if.read((char *) (dat->uh_data[i].loggex.data()),
		       dat->uh_data[i].loggex.size() * sizeof(double));
  }

  // Return the data
  return dat;
}

////////////////////////////////////////////////////////////////////////
// This is a little utility to do interpolation on a monatonic table
////////////////////////////////////////////////////////////////////////
inline double
pwind_hot::interp(const double x,
		  const vector<double>& xtab,
		  const vector<double>& ytab,
		  const int off) const {

  // Index variables
  int xl = 0;
  int xr = xtab.size()-1;

  // Table can be monotonically increasing or decreasing; handle both
  // cases
  if (xtab.front() < xtab.back()) {

    // Increasing
    while (xr - xl > 1) {
      int xi = (xl + xr) / 2;
      if (x < xtab[xi]) xr = xi;
      else xl = xi;
    }

  } else {

    // Decreasing
    while (xr - xl > 1) {
      int xi = (xl + xr) / 2;
      if (x > xtab[xi]) xr = xi;
      else xl = xi;
    }
  }

  // Get weight
  double w = (x - xtab[xl]) / (xtab[xr] - xtab[xl]);

  // Compute interpolated value, using ytab if provided or just from 0
  // to 1 if ytab is empty
  double y;
  if (ytab.size() > 0) {
    y = (1.0-w) * ytab[xl + off] + w * ytab[xr + off];
  } else {
    y = ((1.0-w) * xl + w * xr) / xtab.size();
  }

  // Return
  return y;
}


////////////////////////////////////////////////////////////////////////
// Routines to evaluate the kinematics using the interpolation
// functions
////////////////////////////////////////////////////////////////////////
inline double
pwind_hot::X(const double ur, const double a) const {
  
  // Get indices and weights for bilinear interpolation
  double q = log10(a - 1.0);
  int qidx = (int) ((q - full_tab->qmin) / full_tab->dq);
  double wq = (q - full_tab->qmin) / full_tab->dq - qidx;
  int uidx = (int) ((ur / uh) / full_tab->du);
  double wu = (ur / uh) / full_tab->du - uidx;

  // Handle values of q that are out of bounds
  if (qidx < 0) {
    // Analytic result for a -> 1
    return -log((SQR(ur)/(a-1.0) + 1.0)/Gamma);
  } else if (qidx >= full_tab->nq-1) {
    // For a -> infinity, return -max_double as a flag value
    return -numeric_limits<double>::max();
  }

  // Handle values of u that are out of bounds
  if (uidx >= full_tab->nu-1) return -numeric_limits<double>::max();
  if (hasUmin()) {
    double umin = (1.0 - wq) * tab.umin[qidx] + wq * tab.umin[qidx+1];
    if (ur < umin) return -numeric_limits<double>::max();
  }

  // Do binlinear interpolation in to get log (Gamma exp^-x)
  double loggex =
    (1.0 - wq) * (1.0 - wu) * tab.loggex[uidx + full_tab->nu * qidx] +
    (1.0 - wq) * wu * tab.loggex[uidx + 1 + full_tab->nu * qidx] +
    wq * (1.0 - wu) * tab.loggex[uidx + full_tab->nu * (qidx + 1)] +
    wq * wu * tab.loggex[uidx + 1 + full_tab->nu * (qidx + 1)];

  // Get x and return
  double x = -log(pow(10., loggex) / Gamma);
  return x;
}

inline double
pwind_hot::U2(const double x, const double a) const {

  // Get index and weight for this x in our data table
  double loggex = log10(Gamma * exp(-x));
  int lgidx = (int) (loggex / full_tab->dlg);
  double wlg = loggex / full_tab->dlg - lgidx;
  double q = log10(a - 1.0);

  // Handle special case of q below minimum value in our grid using
  // series solution
  if (q < full_tab->qmin) {
    if (x != xcrit) {
      return (Gamma*exp(-x)-1.0)*(a-1.0);
    } else {
      double u = (a-1.0) / (2.0*uh) *
	(-1.0 + sqrt(1.0+2.0*SQR(uh)*(dyda(1.0) - dmda(1.0))));
      return SQR(u);
    }
  }

  // Handle special case of x off the grid
  if (lgidx >= full_tab->ngex-1) return 0.0;
  
  // Interpolate to get stopping point and maximum u value at this x
  double q_stop = (1.0 - wlg) * tab.q_stop[lgidx] +
    wlg * tab.q_stop[lgidx+1];
  double umax = (1.0 - wlg) * tab.umax[lgidx] +
    wlg * tab.umax[lgidx+1];

  // Fountain and non-fountain cases are handled differently
  if (!isFountain()) {

    // Handle special case where q is off table
    if (q > q_stop) return umax;

    // If we're here, this q is between the min and max values, and is
    // therefore in our table. We therefore need to interpolate the
    // table to this x, then locate this point in it.
    vector<double> qtab(full_tab->nu);
    for (int i=0; i<full_tab->nu; i++)
      qtab[i] = (1.0-wlg) * tab.q[lgidx * full_tab->nu + i] +
	wlg * tab.q[(lgidx+1) * full_tab->nu + i];

    // Now use table to get interpolated value
    return umax * interp(q, qtab);

  } else {

    // Handle special case where q is off table
    if (q > q_stop) return 0.0;

    // Get turnaround point for fountain
    double q_umax = (1.0 - wlg) * tab.q_umax[lgidx] +
      wlg * tab.q_umax[lgidx+1];

    // Figure out whether we are in the accelerating or decelerating
    // part of the fountain, and interpolate the correct part of the q
    // vs. u curve
    vector<double> qtab(full_tab->nu);
    int off = q < q_umax ? 0 : full_tab->ngex * full_tab->nu;
    for (int i=0; i<full_tab->nu; i++)
      qtab[i] = (1.0-wlg) * tab.q[off + lgidx * full_tab->nu + i] +
	wlg * tab.q[off + (lgidx+1) * full_tab->nu + i];

    // Use interpolated table to get u
    return umax * interp(q, qtab);

  }
}

inline double
pwind_hot::dU2dx(const double x, const double a) const {

  // Here we compute the derivative numerically by computing the
  // interpolated values at x values displaced by a fraction of a cell
  // spacing
  double loggex = log10(Gamma * exp(-x));
  double x1 = log(pow(10., -(loggex + 0.01*full_tab->dlg)) * Gamma);
  double x2 = log(pow(10., -(loggex - 0.01*full_tab->dlg)) * Gamma);
  double u1 = U(x1, a);
  double u2 = U(x2, a);
  return (u2*u2 - u1*u1) / (x2 - x1);
}

inline double
pwind_hot::dU2da(const double x, const double a) const {
  return (Gamma*y(a)*exp(-x)*SQR(1.0-U(x,a)/uh) - m(a)) / SQR(a);
}

////////////////////////////////////////////////////////////////////////
// Limits on the range of the wind
////////////////////////////////////////////////////////////////////////

inline double
pwind_hot::amax(const double x) const {
  if (!isFountain()) {
    
    // Wind reaches infinity
    return numeric_limits<double>::max();
    
  } else {
    
    // Get index and weight for this x in our data table
    double loggex = log10(Gamma * exp(-x));
    int lgidx = (int) (loggex / full_tab->dlg);
    double wlg = loggex / full_tab->dlg - lgidx;

    // Interpolate to get stopping point
    double q_stop = (1.0 - wlg) * tab.q_stop[lgidx] +
      wlg * tab.q_stop[lgidx+1];
    return q_stop;

  }
}

// Get limits on x at a given a
vector<double>
pwind_hot::xlimits(const double a) const {

  // Fountain and non-fountain cases handled differently
  if (!isFountain()) {
    
    // Non-fountain case: just return full range covered by table
    vector<double> x(2);    
    x[0] = -log(pow(10., full_tab->loggex_max) / Gamma);
    x[1] = xcrit;
    return x;
    
  } else {

    // Fountain case; in this case we use our table of q_stop to
    // interpolate to find the limits on x at a given a
    double q = log10(a - 1.0);

    // Handle special case where q is bigger than largest q_stop; in
    // this case just return an empty vector
    if (q > tab.q_stop[full_tab->ngex-1]) {
      vector<double> x;
      return x;
    }

    // Lower limit on x is set to limit of table
    vector<double> x(2);
    x[0] = -log(pow(10., full_tab->loggex_max) / Gamma);

    // Interpolate to get upper limit
    x[1] = -log(pow(10., full_tab->loggex_max * interp(q, tab.q_stop))
		/ Gamma);

    // Return interpolated value
    return x;

  }

}

// Get limits on a at a given u, varpi.
vector<double>
pwind_hot::alimits(const double u, const double varpi,
		   const double varpi_t) const {

  // Setup
  double vp2 = SQR(varpi) + SQR(varpi_t);
  double uabs = fabs(u);
  vector<double> alim;

  // To get the lower limit on a, we construct the line of sight
  // velocity along the velocity curve for our lowest column density
  // case
  vector<double> aLOSmax, uLOSmax;
  int off = full_tab->nu * (full_tab->ngex-1);
  if (vp2 >= 1.0) {
    // For vp2 > 1, add a point to the table at a = vp2 exactly
    aLOSmax.push_back( sqrt(vp2) );
    uLOSmax.push_back( 0.0 );
  }
  for (int i=0; i<full_tab->nu; i++) {
    double a = pow(10., tab.q[i + off]) + 1.0;
    double a2 = a*a;
    if (a2 >= vp2) {
      aLOSmax.push_back(a);
      uLOSmax.push_back(((double) i) / (full_tab->nu - 1) *
			tab.umax[full_tab->ngex-1] * sqrt(1.0 - vp2/a2));
    }
  }

  // Handle special case of input velocity > largest value in our table
  if (uabs > uLOSmax.back()) return alim;

  // Find radius at which the uLOSmax curve hits the input value of
  // u; this will be the minimum radius
  alim.push_back( interp(uabs, uLOSmax, aLOSmax) );

  // How we proceed next depends on the type of wind we have. We
  // distinguish three cases
  if (!isFountain() && !hasUmin()) {

    // Case of a wind that is not a fountain, and also has material
    // going down to zero velocity at all radii. In this case the
    // solution is trivial, and the maximum radius is just wherever
    // our table ends
    alim.push_back( pow(10., tab.q_stop[full_tab->ngex-1]) + 1 );

  } else if (!isFountain()) {

    // Case of a wind that is not a fountain, but where even material
    // at x = x_crit is at fininte velocity > 0 for a > 1. In this
    // case, there may be a finite maximum a, which we must find by
    // tracing out the velocity of material at x = xcrit
    vector<double> aLOSmin, uLOSmin;
    if (vp2 >= 1.0) {
      aLOSmin.push_back( sqrt(vp2) );
      uLOSmin.push_back( 0.0 );
    }
    for (int i=0; i<full_tab->nq; i++) {
      double q = full_tab->qmin + i * full_tab->dq;
      double a = pow(10., q) + 1.0;
      double a2 = a*a;
      if (a2 >= vp2) {
	aLOSmin.push_back(a);
	uLOSmin.push_back(tab.umin[i] * sqrt(1.0 - vp2/a2));
      }
    }

    // Handle special case where the velocity of the material at x =
    // xcrit is sill below our target LOS velocity
    if (uLOSmin.back() < uabs) {
      alim.push_back( pow(10., tab.q_stop[full_tab->ngex-1]) + 1 );
    } else {
      alim.push_back( interp(uabs, uLOSmin, aLOSmin) );
    }
    
  } else {

    // Fountain case; in this case our highest velocity value of
    // loggex has both an accelerating and decelerating part of its
    // U(a) curve. We have already found the solution on the
    // accelerating part of the velocity curve, and we now need to
    // find the solution on the decelerating part
    uLOSmax.resize(0);
    aLOSmax.resize(0);
    off = full_tab->nu * (2 * full_tab->ngex-1);
    for (int i=0; i<full_tab->nu; i++) {
      double a = pow(10., tab.q[i + off]) + 1.0;
      double a2 = a*a;
      if (a2 >= vp2) {
	aLOSmax.push_back(a);
	uLOSmax.push_back(((double) i) / (full_tab->nu - 1) *
			  tab.umax[full_tab->ngex-1] * sqrt(1.0 - vp2/a2) );
      }
    }

    // For decelerating part, there are two possibilities: either the
    // largest value of loggex that we follow reaches velocity < u at
    // on our grid, or it does not; in the latter case, we just return
    // the grid edge as the upper limit on a, and in the former case
    // we find the solution as for the accelerating part
    if (uabs < uLOSmax.back()) {
      alim.push_back( pow(10., full_tab->qmax) + 1.0 );
    } else {
      alim.push_back ( interp(uabs, uLOSmax, aLOSmax) );
    }

  }

  // Return result
  return alim;
}


