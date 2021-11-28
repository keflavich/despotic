#include "pwind_interface.H"
#include <iostream>
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Allocation and de-allocation functions
////////////////////////////////////////////////////////////////////////

// Free -- same for all pwinds
void pwind_free(pwind *pw) { delete pw; }
void pwind_geom_free(pwind_geom *geom) { delete geom; }

// Geometries
pwind_geom_sphere *
pwind_geom_sphere_new(const double phi) {
  return new pwind_geom_sphere(phi);
}

pwind_geom_cone *
pwind_geom_cone_new(const double theta, const double phi) {
  return new pwind_geom_cone(theta, phi);
}

pwind_geom_cone_sheath *
pwind_geom_cone_sheath_new(const double theta_out,
			   const double theta_in,
			   const double phi) {
  return new pwind_geom_cone_sheath(theta_out, theta_in, phi);
}

// Ideal
pwind_ideal_pa *
pwind_ideal_pa_new(const double Gamma,
		   const double mach,
		   const pwind_geom* geom,
		   const double fcrit,
		   const double jsp) {
  return new pwind_ideal_pa(Gamma, mach, geom, fcrit, jsp);
}
pwind_ideal_pi *
pwind_ideal_pi_new(const double Gamma,
		   const double mach,
		   const pwind_geom* geom,
		   const double fcrit,
		   const double jsp) {
  return new pwind_ideal_pi(Gamma, mach, geom, fcrit, jsp);
}
pwind_ideal_ps *
pwind_ideal_ps_new(const double Gamma,
		   const double mach,
		   const pwind_geom* geom,
		   const double fcrit,
		   const double jsp) {
  return new pwind_ideal_ps(Gamma, mach, geom, fcrit, jsp);
}
pwind_ideal_ia *
pwind_ideal_ia_new(const double Gamma,
		   const double mach,
		   const pwind_geom* geom,
		   const double fcrit,
		   const double jsp) {
  return new pwind_ideal_ia(Gamma, mach, geom, fcrit, jsp);
}
pwind_ideal_ii *
pwind_ideal_ii_new(const double Gamma,
		   const double mach,
		   const pwind_geom* geom,
		   const double fcrit,
		   const double jsp) {
  return new pwind_ideal_ii(Gamma, mach, geom, fcrit, jsp);
}
pwind_ideal_is *
pwind_ideal_is_new(const double Gamma,
		   const double mach,
		   const pwind_geom* geom,
		   const double fcrit,
		   const double jsp) {
  return new pwind_ideal_is(Gamma, mach, geom, fcrit, jsp);
}

// Radiation-driven
pwind_rad_pa *
pwind_rad_pa_new(const double Gamma,
		 const double mach,
		 const double tau0,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp) {
  return new pwind_rad_pa(Gamma, mach, tau0, geom, fcrit, jsp);
}
pwind_rad_pi *
pwind_rad_pi_new(const double Gamma,
		 const double mach,
		 const double tau0,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp) {
  return new pwind_rad_pi(Gamma, mach, tau0, geom, fcrit, jsp);
}
pwind_rad_ps *
pwind_rad_ps_new(const double Gamma,
		 const double mach,
		 const double tau0,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp) {
  return new pwind_rad_ps(Gamma, mach, tau0, geom, fcrit, jsp);
}
pwind_rad_ia *
pwind_rad_ia_new(const double Gamma,
		 const double mach,
		 const double tau0,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp) {
  return new pwind_rad_ia(Gamma, mach, tau0, geom, fcrit, jsp);
}
pwind_rad_ii *
pwind_rad_ii_new(const double Gamma,
		 const double mach,
		 const double tau0,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp) {
  return new pwind_rad_ii(Gamma, mach, tau0, geom, fcrit, jsp);
}
pwind_rad_is *
pwind_rad_is_new(const double Gamma,
		 const double mach,
		 const double tau0,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp) {
  return new pwind_rad_is(Gamma, mach, tau0, geom, fcrit, jsp);
}

// Hot gas-driven
pwind_hot_pa *
pwind_hot_pa_new(const double Gamma,
		 const double mach,
		 const double uh,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp,
		 hot_wind_table *tab) {
  return new pwind_hot_pa(Gamma, mach, uh, geom,
			  fcrit, jsp, tab);
}
pwind_hot_pi *
pwind_hot_pi_new(const double Gamma,
		 const double mach,
		 const double uh,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp,
		 hot_wind_table *tab) {
  return new pwind_hot_pi(Gamma, mach, uh, geom, 
			  fcrit, jsp, tab);
}
pwind_hot_ps *
pwind_hot_ps_new(const double Gamma,
		 const double mach,
		 const double uh,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp,
		 hot_wind_table *tab) {
  return new pwind_hot_ps(Gamma, mach, uh, geom, 
			  fcrit, jsp, tab);
}
pwind_hot_ia *
pwind_hot_ia_new(const double Gamma,
		 const double mach,
		 const double uh,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp,
		 hot_wind_table *tab) {
  return new pwind_hot_ia(Gamma, mach, uh, geom, 
			  fcrit, jsp, tab);
}
pwind_hot_ii *
pwind_hot_ii_new(const double Gamma,
		 const double mach,
		 const double uh,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp,
		 hot_wind_table *tab) {
  return new pwind_hot_ii(Gamma, mach, uh, geom, 
			  fcrit, jsp, tab);
}
pwind_hot_is *
pwind_hot_is_new(const double Gamma,
		 const double mach,
		 const double uh,
		 const pwind_geom* geom,
		 const double fcrit,
		 const double jsp,
		 hot_wind_table *tab) {
  return new pwind_hot_is(Gamma, mach, uh, geom, 
			  fcrit, jsp, tab);
}


////////////////////////////////////////////////////////////////////////
// Computation functions
////////////////////////////////////////////////////////////////////////

// Limits
unsigned long alimits(const double u,
		      const double varpi,
		      const double varpi_t,
		      const double epsabs,
		      const double epsrel,
		      const pwind *pw,
		      double *alim) {
  std::vector<double> alim_ = pw->alimits(u, varpi, varpi_t,
					  epsabs, epsrel);
  for (std::vector<double>::size_type i=0; i<alim_.size(); i++)
    alim[i] = alim_[i];
  return alim_.size();
}
unsigned long xlimits(const double a, const pwind *pw, double *xlim) {
  std::vector<double> xlim_ = pw->xlimits(a);
  for (std::vector<double>::size_type i=0; i<xlim_.size(); i++)
    xlim[i] = xlim_[i];
  return xlim_.size();
}
double amax(const double x,
	    const double epsabs,
	    const double epsrel,
	    const pwind *pw) {
  return pw->amax(x, epsabs, epsrel);
}
unsigned long s_crit(const double varpi,
		     const double varpi_t,
		     const double u,
		     const pwind *pw, double *s_crit_) {
  std::vector<double> sc = pw->s_crit(varpi, varpi_t, u);
  for (std::vector<double>::size_type i=0; i<sc.size(); i++)
    s_crit_[i] = sc[i];
  return sc.size();
}
unsigned long a_crit(const double varpi,
		     const double varpi_t,
		     const double u,
		     const pwind *pw,
		     double *a_crit_) {
  std::vector<double> ac = pw->a_crit(varpi, varpi_t, u);
  for (std::vector<double>::size_type i=0; i<ac.size(); i++)
    a_crit_[i] = ac[i];
  return ac.size();
}

// Density
double drhodx(const double x,
	      const double a,
	      const pwind *pw) {
  return pw->drhodx(x, a);
}
double rho(const double a,
	   const double epsabs,
	   const double epsrel,
	   const pwind *pw) {
  return pw->rho(a, epsabs, epsrel);
}
void rho_vec(const unsigned long na,
	     const double *a,
	     const double epsabs,
	     const double epsrel,
	     const pwind *pw,
	     double *result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<na; i++)
    result[i] = pw->rho(a[i], epsabs, epsrel);
}

// Momentum flux
double pdot(const double a,
	    const double epsabs,
	    const double epsrel,
	    const pwind *pw) {
  return pw->pdot(a, epsabs, epsrel);
}
void pdot_vec(const unsigned long na,
	      const double *a,
	      const double epsabs,
	      const double epsrel,
	      const pwind *pw,
	      double *result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<na; i++)
    result[i] = pw->pdot(a[i], epsabs, epsrel);
}


// Momentum flux normalised to driving momentum flux
double pdotRel_approx(const double a,
		      const double epsabs,
		      const double epsrel,
		      const pwind *pw) {
  return pw->pdotRel(a, epsabs, epsrel);
}
double pdotRel_exact(const double a,
		     const double fg,
		     const double tctw,
		     const double epsabs,
		     const double epsrel,
		     const pwind *pw) {
  return pw->pdotRel(a, fg, tctw, epsabs, epsrel);
}

// Energy flux
double Edot(const double a,
	    const double epsabs,
	    const double epsrel,
	    const pwind *pw) {
  return pw->Edot(a, epsabs, epsrel);
}
void Edot_vec(const unsigned long na,
	      const double *a,
	      const double epsabs,
	      const double epsrel,
	      const pwind *pw,
	      double *result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<na; i++)
    result[i] = pw->Edot(a[i], epsabs, epsrel);
}


// Absorption
double Phi_uc(const double u,
	      const double varpi,
	      const double varpi_t,
	      const double a0,
	      const double a1,
	      const double epsabs,
	      const double epsrel,
	      const pwind *pw) {
  return pw->Phi_uc(u, varpi, varpi_t, a0, a1, epsabs, epsrel);
}
void Phi_uc_vec(const unsigned long nu,
		const double *u,
		const double varpi,
		const double varpi_t,
		const double a0,
		const double a1,
		const double epsabs,
		const double epsrel,
		const pwind *pw,
		double *result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<nu; i++)
    result[i] = pw->Phi_uc(u[i], varpi, varpi_t, a0, a1, epsabs, epsrel);
}
double tau_uc(const double u,
	      const double tXtw,
	      const double fj,
	      const double boltzfac,
	      const double varpi,
	      const double varpi_t,
	      const double a0,
	      const double a1,
	      const double epsabs,
	      const double epsrel,
	      const pwind *pw) {
  return pw->tau_uc(u, tXtw, fj, boltzfac, varpi, varpi_t, a0, a1,
		    epsabs, epsrel);
}
void tau_uc_vec(const unsigned long nu,
		const double *u,
		const double tXtw,
		const double fj,
		const double boltzfac,
		const double varpi,
		const double varpi_t,
		const double a0,
		const double a1,
		const double epsabs,
		const double epsrel,
		const pwind *pw,
		double *result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<nu; i++)
    result[i] = pw->tau_uc(u[i], tXtw, fj, boltzfac, varpi, varpi_t, a0, a1,
			   epsabs, epsrel);
}
double tau_uc_multiple(const double u,
		       const double *u_trans,
		       const double *tXtw,
		       const double fj,
		       const double boltzfac,
		       const unsigned long ntrans,
		       const double varpi,
		       const double varpi_t,
		       const double a0,
		       const double a1,
		       const double epsabs,
		       const double epsrel,
		       const pwind *pw) {
  vector<double> u_trans_, tXtw_;
  u_trans_.assign(u_trans, u_trans+ntrans);
  tXtw_.assign(tXtw, tXtw+ntrans);
  return pw->tau_uc(u, u_trans_, tXtw_, fj, boltzfac, varpi, varpi_t,
		    a0, a1, epsrel, epsabs);
}
void tau_uc_multiple_vec(const unsigned long nu,
			 const double *u,
			 const double *u_trans,
			 const double *tXtw,
			 const double fj,
			 const double boltzfac,
			 const unsigned long ntrans,
			 const double varpi,
			 const double varpi_t,
			 const double a0,
			 const double a1,
			 const double epsabs,
			 const double epsrel,
			 const pwind *pw,
			 double *result) {
  vector<double> u_trans_, tXtw_;
  u_trans_.assign(u_trans, u_trans+ntrans);
  tXtw_.assign(tXtw, tXtw+ntrans);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<nu; i++)
    result[i] = pw->tau_uc(u[i], u_trans_, tXtw_, fj, boltzfac,
			   varpi, varpi_t,
			   a0, a1, epsrel, epsabs);
}
double Phi_c(const double u,
	     const double fw,
	     const double varpi,
	     const double varpi_t,
	     const double a0,
	     const double a1,
	     const double epsabs,
	     const double epsrel,
	     const pwind *pw) {
  return pw->Phi_c(u, fw, varpi, varpi_t, a0, a1, epsabs, epsrel);
}
void Phi_c_vec(const unsigned long nu,
	       const double *u,
	       const double fw,
	       const double varpi,
	       const double varpi_t,
	       const double a0,
	       const double a1,
	       const double epsabs,
	       const double epsrel,
	       const pwind *pw,
	       double *result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<nu; i++)
    result[i] = pw->Phi_c(u[i], fw, varpi, varpi_t, a0, a1, epsabs, epsrel);
}
double tau_c(const double u,
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
	     const pwind *pw) {
  return pw->tau_c(u, tXtw, fj, boltzfac, fw, varpi, varpi_t, a0, a1,
		   epsabs, epsrel);
}
void tau_c_vec(const unsigned long nu,
	       const double *u,
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
	       const pwind *pw,
	       double *result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<nu; i++)
    result[i] = pw->tau_c(u[i], tXtw, fj, boltzfac, fw, varpi, varpi_t, a0, a1,
			  epsabs, epsrel);
}
double tau_c_multiple(const double u,
		      const double *u_trans,
		      const double *tXtw,
		      const double fj,
		      const double boltzfac,
		      const unsigned long ntrans,
		      const double fw,
		      const double varpi,
		      const double varpi_t,
		      const double a0,
		      const double a1,
		      const double epsabs,
		      const double epsrel,
		      const pwind *pw) {
  vector<double> u_trans_, tXtw_;
  u_trans_.assign(u_trans, u_trans+ntrans);
  tXtw_.assign(tXtw, tXtw+ntrans);
  return pw->tau_c(u, u_trans_, tXtw_, fj, boltzfac, fw, varpi, varpi_t,
		   a0, a1, epsabs, epsrel);
}
void tau_c_multiple_vec(const unsigned long nu,
			const double *u,
			const double *u_trans,
			const double *tXtw,
			const double fj,
			const double boltzfac,
			const unsigned long ntrans,
			const double fw,
			const double varpi,
			const double varpi_t,
			const double a0,
			const double a1,
			const double epsabs,
			const double epsrel,
			const pwind *pw,
			double *result) {
  vector<double> u_trans_, tXtw_;
  u_trans_.assign(u_trans, u_trans+ntrans);
  tXtw_.assign(tXtw, tXtw+ntrans);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<nu; i++)
    result[i] = pw->tau_c(u[i], u_trans_, tXtw_, fj, boltzfac, fw,
			  varpi, varpi_t,
			  a0, a1, epsabs, epsrel);
}

// Thin, subcritical emission
double Xi(const double u,
	  const double varpi,
	  const double varpi_t,
	  const double epsabs,
	  const double epsrel,
	  const pwind *pw) {
  return pw->Xi(u, varpi, varpi_t, epsabs, epsrel);
}
void Xi_vec(const unsigned long nu,
	    const double *u,
	    const double varpi,
	    const double varpi_t,
	    const double epsabs,
	    const double epsrel,
	    const pwind *pw,
	    double *result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<nu; i++)
    result[i] = pw->Xi(u[i], varpi, varpi_t, epsabs, epsrel);
}
double xi(const double varpi,
	  const double varpi_t,
	  const double epsabs,
	  const double epsrel,
	  const pwind *pw) {
  return pw->xi(varpi, varpi_t, epsabs, epsrel);
}

// LTE emission
double eta(const double u,
	   const double tXtw,
	   const double fj,
	   const double boltzfac,
	   const bool correlated,
	   const double fw,
	   const double varpi,
	   const double varpi_t,
	   const bool thin,
	   const double epsabs,
	   const double epsrel,
	   const pwind *pw) {
  return pw->eta(u, tXtw, fj, boltzfac, correlated, fw, varpi, varpi_t,
		 thin, epsabs, epsrel);
}
void eta_vec(const unsigned long nu,
	     const double *u,
	     const double tXtw,
	     const double fj,
	     const double boltzfac,
	     const bool correlated,
	     const double fw,
	     const double varpi,
	     const double varpi_t,
	     const bool thin,
	     const double epsabs,
	     const double epsrel,
	     const pwind *pw,
	     double *result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
  for (unsigned long i=0; i<nu; i++)
    result[i] = pw->eta(u[i], tXtw, fj, boltzfac, correlated, fw,
			varpi, varpi_t,
			thin, epsabs, epsrel);
}
double Psi(const double tXtw,
	   const double fj,
	   const double boltzfac,
	   const bool correlated,
	   const double fw,
	   const double varpi,
	   const double varpi_t,
	   const bool thin,
	   const double epsabs,
	   const double epsrel,
	   const pwind *pw) {
  return pw->Psi(tXtw, fj, boltzfac, correlated, fw, varpi, varpi_t, thin,
		 epsabs, epsrel);
}

// Expose methods that set and get simple information from pwind objects
double y(const double a, const pwind *pw) {
  return pw->y(a);
}
double dyda(const double a, const pwind *pw) {
  return pw->dyda(a);
}
double m(const double a, const pwind *pw) {
  return pw->m(a);
}
double X(const double ur, const double a, const pwind *pw) {
  return pw->X(ur, a);
}
double U2(const double x, const double a, const pwind *pw) {
  return pw->U2(x, a);
}
double dU2dx(const double x, const double a, const pwind *pw) {
  return pw->dU2dx(x, a);
}
double dU2da(const double x, const double a, const pwind *pw) {
  return pw->dU2da(x, a);
}
double Gamma(const pwind *pw) { return pw->getGamma(); }
double mach(const pwind *pw) { return pw->getMach(); }
double fcrit(const pwind *pw) { return pw->getFcrit(); }
double jsp(const pwind *pw) { return pw->getJsp(); }
double xcrit(const pwind *pw) { return pw->xcr(); }
double sx(const pwind *pw) { return pw->s(); }
double zetaM(const pwind *pw) { return pw->zM(); }
double zetaA(const pwind *pw) { return pw->zA(); }
double umax(const pwind *pw) { return pw->uMax(); }
double amax_abs(const pwind *pw) { return pw->aMaxAbs(); }
void set_mach(const double mach, pwind *pw) { pw->setMach(mach); }
void set_geometry(const pwind_geom *geom, pwind *pw)
{ pw->setGeometry(geom); }
void set_fcrit(const double fcrit, pwind *pw)
{ pw->setFcrit(fcrit); }
void set_jsp(const double jsp, pwind *pw) { pw->setJsp(jsp); }

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

// Expose methods for handling errors
int get_err(const pwind *pw) { return pw->getErr(); }
void set_err(const int err, pwind *pw) { pw->setErr(err); }
void clear_err(pwind *pw) { pw->clearErr(); }
const char *get_err_str(const pwind *pw) { return pw->getErrStr(); }

// Functions to manage the hot gas tabulated data
hot_wind_table *read_hot_wind_table(const char *dirname,
				    const int yidx,
				    const int midx) {
  return pwind_hot::read_table(std::string(dirname), yidx, midx);
}
void free_hot_wind_table(hot_wind_table *tab) {
  pwind_hot::free_table(tab);
}
void get_hot_wind_table_limits(hot_wind_table *tab,
			       double *uh_limits) {
  uh_limits[0] = tab->uh.front();
  uh_limits[1] = tab->uh.back();
}
