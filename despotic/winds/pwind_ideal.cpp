#include "pwind_ideal.H"
#include <limits>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Base pwind_ideal class
////////////////////////////////////////////////////////////////////////
pwind_ideal::pwind_ideal(const double Gamma_, const double mach_,
			 const pwind_potential *potential_,
			 const pwind_expansion *expansion_,
			 const pwind_geom *geom_,
			 const double fcrit_,
			 const double jsp_) :
  pwind(Gamma_, mach_, potential_, expansion_, geom_, 
	fcrit_, jsp_)
{
  xcrit = log(Gamma_);
  umax = numeric_limits<double>::max();
  zeta_M = zetaM(xcrit + log(fcrit_), sx);
  zeta_A = zetaA(xcrit + log(fcrit_), sx);
  amax_abs = numeric_limits<double>::max();
}

double
pwind_ideal::dU2da(const double x, const double a) const {
  return (y(a)*Gamma*exp(-x) - m(a))/SQR(a);
}

////////////////////////////////////////////////////////////////////////
// Classes specialised by geometry and potential
////////////////////////////////////////////////////////////////////////

// Point, area
double 
pwind_ideal_pa::X(const double ur, const double a) const {
  return -log((SQR(ur)/(1.0-1.0/a)+1.0)/Gamma);
}
double 
pwind_ideal_pa::U2(const double x, const double a) const {
  return (Gamma*exp(-x)-1.0)*(1.0-1.0/a);
}
double 
pwind_ideal_pa::dU2dx(const double x, const double a) const {
  return -(1.0-1.0/a)*Gamma*exp(-x);
}
vector<double>
pwind_ideal_pa::alimits(const double u,
			const double varpi,
			const double varpi_t,
			const double epsabs,
			const double epsrel) const {
  (void) epsabs;   // Not used; this can be done analytically
  (void) epsrel;   // Not used; this can be done analytically
  vector<double> ret(2);
  ret[0] = max(1.0, sqrt(SQR(varpi)+SQR(varpi_t)));
  ret[1] = numeric_limits<double>::max();
  return ret;
}
vector<double>
pwind_ideal_pa::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = -numeric_limits<double>::max();
  ret[1] = xcrit;
  return ret;
}
double
pwind_ideal_pa::amax(const double x,
		     const double epsabs,
		     const double epsrel) const {
  return numeric_limits<double>::max();
}
double
pwind_ideal_pa::Phi_c(const double u,
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
pwind_ideal_pa::Phi_c(const vector<double>& u,
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
    result[i] = Phi_c(u[i], fw, varpi, varpi_t, a0, a1, epsabs, epsrel);

  // Return
  return result;
}

// Point, intermediate
double 
pwind_ideal_pi::X(const double ur, const double a) const {
  return -log((SQR(ur)+(1.0-1.0/a)) / (Gamma*log(a)));
}
double 
pwind_ideal_pi::U2(const double x, const double a) const {
  return Gamma*exp(-x)*log(a) - (1.0-1.0/a);
}
double 
pwind_ideal_pi::dU2dx(const double x, const double a) const {
  return -Gamma*log(a)*exp(-x);
}
vector<double>
pwind_ideal_pi::alimits(const double u,
			const double varpi,
			const double varpi_t,
			const double epsabs,
			const double epsrel) const {
  // Lower limit is easy, but upper limit must be obtained numerically  
  vector<double> ret(2);
  ret[0] = max(1.0, sqrt(SQR(varpi)+SQR(varpi_t)));
  ret[1] = a_from_u_x(u, xcrit, varpi, varpi_t, ret[0], epsabs, epsrel);
  if (ret[1] < 0) {
    // Since u->infinity with a but only logarithmically, for large u
    // we may not be able to find the limit; handle this by just
    // setting the upper limit to infinity
    ret[1] = numeric_limits<double>::max();
  }
  return ret;
}
vector<double>
pwind_ideal_pi::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = -numeric_limits<double>::max();
  ret[1] = xcrit;
  return ret;
}
double
pwind_ideal_pi::amax(const double x,
		     const double epsabs,
		     const double epsrel) const {
  return numeric_limits<double>::max();
}

// Point, solid angle
double 
pwind_ideal_ps::X(const double ur, const double a) const {
  return -log((SQR(ur)/(a-1.0)+1.0/a)/Gamma);
}
double 
pwind_ideal_ps::U2(const double x, const double a) const {
  return (Gamma*exp(-x)-1.0/a)*(a-1.0);
}
double 
pwind_ideal_ps::dU2dx(const double x, const double a) const {
  return -(a-1.0)*Gamma*exp(-x);
}
vector<double>
pwind_ideal_ps::alimits(const double u,
			const double varpi,
			const double varpi_t,
			const double epsabs,
			const double epsrel) const {
  vector<double> ret(2);
  ret[0] = max(1.0, sqrt(SQR(varpi)+SQR(varpi_t)));
  ret[1] = a_from_u_x(u, xcrit, varpi, varpi_t, ret[0], epsabs, epsrel);  
  if (ret[1] < 0) {
    // Since u->infinity with a but only logarithmically, for large u
    // we may not be able to find the limit; handle this by just
    // setting the upper limit to infinity
    ret[1] = numeric_limits<double>::max();
  }
  return ret;
}
vector<double>
pwind_ideal_ps::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = -numeric_limits<double>::max();
  ret[1] = xcrit;
  return ret;
}
double
pwind_ideal_ps::amax(const double x,
		     const double epsabs,
		     const double epsrel) const {
  return numeric_limits<double>::max();
}

// Isothermal, area
double 
pwind_ideal_ia::X(const double ur, const double a) const {
  return log((1-1.0/a)*Gamma/(SQR(ur)+log(a)));
}
double 
pwind_ideal_ia::U2(const double x, const double a) const {
  return Gamma*exp(-x)*(1.0-1.0/a) - log(a);
}
double 
pwind_ideal_ia::dU2dx(const double x, const double a) const {
  return -(1.0-1.0/a)*Gamma*exp(-x);
}
vector<double>
pwind_ideal_ia::alimits(const double u,
			const double varpi,
			const double varpi_t,
			const double epsabs,
			const double epsrel) const {
  (void) epsabs;   // Not used; this case can be done analytically
  (void) epsrel;   // Not used; this case can be done analytically
  vector<double> ret(2);
  ret[0] = max(1.0, sqrt(SQR(varpi)+SQR(varpi_t)));
  ret[1] = numeric_limits<double>::max();
  return ret;
}
vector<double>
pwind_ideal_ia::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = -numeric_limits<double>::max();
  ret[1] = -log(a*log(a) / (Gamma*(a-1.0)));
  return ret;
}
double
pwind_ideal_ia::amax(const double x,
		     const double epsabs,
		     const double epsrel) const {
  double amx = a_from_u_x(0.0, x, 0.0, 1.0+1.0e-8,
			  exp(Gamma*exp(-x)),
			  100, epsabs, epsrel);
  if (amx > 0.0) return amx;
  else return 1.0+1.0e-8; // Safety value
}
double
pwind_ideal_ia::Phi_c(const double u,
		      const double fw,
		      const double varpi,
		      const double varpi_t,
		      const double a0,
		      const double a1,
		      const double epsabs,
		      const double epsrel,
		      const vector<double>& alimits_) const {
  if (a1 >= numeric_limits<double>::max() && u < umax)
    return numeric_limits<double>::max();
  else
    return pwind::Phi_c(u, fw, varpi, varpi_t, a0, a1,
			epsabs, epsrel, alimits_);
}
vector<double>
pwind_ideal_ia::Phi_c(const vector<double>& u,
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
    result[i] = Phi_c(u[i], fw, varpi, varpi_t, a0, a1, epsabs, epsrel);

  // Return
  return result;
}

// Isothermal, intermediate
double 
pwind_ideal_ii::X(const double ur, const double a) const {
  return -log((SQR(ur)/log(a)+1.0)/Gamma);
}
double 
pwind_ideal_ii::U2(const double x, const double a) const {
  return (Gamma*exp(-x)-1.0) * log(a);
}
double 
pwind_ideal_ii::dU2dx(const double x, const double a) const {
  return -Gamma*exp(-x)*log(a);
}
vector<double>
pwind_ideal_ii::alimits(const double u,
			const double varpi,
			const double varpi_t,
			const double epsabs,
			const double epsrel) const {
  vector<double> ret(2);
  ret[0] = max(1.0, sqrt(SQR(varpi)+SQR(varpi_t)));
  ret[1] = numeric_limits<double>::max();
  return ret;
}
vector<double>
pwind_ideal_ii::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = -numeric_limits<double>::max();
  ret[1] = xcrit;
  return ret;
}
double
pwind_ideal_ii::amax(const double x,
		     const double epsabs,
		     const double epsrel) const {
  return numeric_limits<double>::max();
}
double
pwind_ideal_ii::Phi_c(const double u,
		      const double fw,
		      const double varpi,
		      const double varpi_t,
		      const double a0,
		      const double a1,
		      const double epsabs,
		      const double epsrel,
		      const vector<double>& alimits_) const {
  if (a1 >= numeric_limits<double>::max() && u < umax)
    return numeric_limits<double>::max();
  else
    return pwind::Phi_c(u, fw, varpi, varpi_t, a0, a1,
			epsabs, epsrel, alimits_);
}
vector<double>
pwind_ideal_ii::Phi_c(const vector<double>& u,
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
    result[i] = Phi_c(u[i], fw, varpi, varpi_t, a0, a1, epsabs, epsrel);

  // Return
  return result;
}

// Isothermal, solid angle
double 
pwind_ideal_is::X(const double ur, const double a) const {
  return log((a-1.0)*Gamma/(SQR(ur)+log(a)));
}
double 
pwind_ideal_is::U2(const double x, const double a) const {
  return Gamma*exp(-x)*(a-1.0) - log(a);
}
double 
pwind_ideal_is::dU2dx(const double x, const double a) const {
  return -Gamma*exp(-x)*(a-1.0);
}
vector<double>
pwind_ideal_is::alimits(const double u,
			const double varpi,
			const double varpi_t,
			const double epsabs,
			const double epsrel) const {
  vector<double> ret(2);
  ret[0] = max(1.0, sqrt(SQR(varpi)+SQR(varpi_t)));
  ret[1] = a_from_u_x(u, xcrit, varpi, varpi_t, epsabs, epsrel);
  return ret;
}
vector<double>
pwind_ideal_is::xlimits(const double a) const {
  vector<double> ret(2);
  ret[0] = -numeric_limits<double>::max();
  ret[1] = xcrit;
  return ret;
}
double
pwind_ideal_is::amax(const double x,
		     const double epsabs,
		     const double epsrel) const {
  return numeric_limits<double>::max();
}

