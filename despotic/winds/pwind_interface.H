// Expose methods to c so we can call them from python

#include "pwind_geom.H"
#include "pwind_ideal.H"
#include "pwind_rad.H"
#include "pwind_hot.H"
#include <iostream>
#ifdef _OPENMP
#include "omp.h"
#endif

extern "C" {
  ////////////////////////////////////////////////////////////////////////
  // Allocation and de-allocation functions
  ////////////////////////////////////////////////////////////////////////

  // Free
  void pwind_free(pwind *pw);
  void pwind_geom_free(pwind_geom *geom);

  // Geometric allocation
  pwind_geom_sphere *pwind_geom_sphere_new(const double phi);

  pwind_geom_cone *
  pwind_geom_cone_new(const double theta, const double phi);

  pwind_geom_cone_sheath *
  pwind_geom_cone_sheath_new(const double theta_out,
			     const double theta_in,
			     const double phi);
  
  // Ideal
  pwind_ideal_pa *
  pwind_ideal_pa_new(const double Gamma, const double mach,
		     const pwind_geom* geom,
		     const double fcrit,
		     const double jsp);
  pwind_ideal_pi *
  pwind_ideal_pi_new(const double Gamma, const double mach,
		     const pwind_geom* geom,
		     const double fcrit,
		     const double jsp);
  pwind_ideal_ps *
  pwind_ideal_ps_new(const double Gamma, const double mach,
		     const pwind_geom* geom,
		     const double fcrit,
		     const double jsp);
  pwind_ideal_ia *
  pwind_ideal_ia_new(const double Gamma, const double mach,
		     const pwind_geom* geom,
		     const double fcrit,
		     const double jsp);
  pwind_ideal_ii *
  pwind_ideal_ii_new(const double Gamma, const double mach,
		     const pwind_geom* geom,
		     const double fcrit,
		     const double jsp);
  pwind_ideal_is *
  pwind_ideal_is_new(const double Gamma, const double mach,
		     const pwind_geom* geom,
		     const double fcrit,
		     const double jsp);

  // Radiation-driven
  pwind_rad_pa *
  pwind_rad_pa_new(const double Gamma, const double mach,
		   const double tau0, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp);
  pwind_rad_pi *
  pwind_rad_pi_new(const double Gamma, const double mach,
		   const double tau0, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp);
  pwind_rad_ps *
  pwind_rad_ps_new(const double Gamma, const double mach,
		   const double tau0, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp);
  pwind_rad_ia *
  pwind_rad_ia_new(const double Gamma, const double mach,
		   const double tau0, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp);
  pwind_rad_ii *
  pwind_rad_ii_new(const double Gamma, const double mach,
		   const double tau0, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp);
  pwind_rad_is *
  pwind_rad_is_new(const double Gamma, const double mach,
		   const double tau0, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp);

  // Hot gas-driven
  pwind_hot_pa *
  pwind_hot_pa_new(const double Gamma, const double mach,
		   const double uh, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp,
		   hot_wind_table *tab);
  pwind_hot_pi *
  pwind_hot_pi_new(const double Gamma, const double mach,
		   const double uh, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp,
		   hot_wind_table *tab);
  pwind_hot_ps *
  pwind_hot_ps_new(const double Gamma, const double mach,
		   const double uh, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp,
		   hot_wind_table *tab);
  pwind_hot_ia *
  pwind_hot_ia_new(const double Gamma, const double mach,
		   const double uh, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp,
		   hot_wind_table *tab);
  pwind_hot_ii *
  pwind_hot_ii_new(const double Gamma, const double mach,
		   const double uh, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp,
		   hot_wind_table *tab);
  pwind_hot_is *
  pwind_hot_is_new(const double Gamma, const double mach,
		   const double uh, const pwind_geom* geom,
		   const double fcrit,
		   const double jsp,
		   hot_wind_table *tab);
 
  ////////////////////////////////////////////////////////////////////////
  // Expose methods that set and return simple information from pwind
  // objects
  ////////////////////////////////////////////////////////////////////////
  double y(const double a, const pwind *pw);
  double dyda(const double a, const pwind *pw);
  double m(const double a, const pwind *pw);
  double X(const double ur, const double a, const pwind *pw);
  double U2(const double x, const double a, const pwind *pw);
  double dU2dx(const double x, const double a, const pwind *pw);
  double dU2da(const double x, const double a, const pwind *pw);
  double Gamma(const pwind *pw);
  double mach(const pwind *pw);
  double fcrit(const pwind *pw);
  double jsp(const pwind *pw);
  double xcrit(const pwind *pw);
  double sx(const pwind *pw);
  double zetaM(const pwind *pw);
  double zetaA(const pwind *pw);
  double umax(const pwind *pw);
  double amax_abs(const pwind *pw);
  void set_mach(const double mach, pwind *pw);
  void set_geometry(const pwind_geom *geom, pwind *pw);
  void set_fcrit(const double fcrit, pwind *pw);
  void set_jsp(const double jsp, pwind *pw);
  
  ////////////////////////////////////////////////////////////////////////
  // Computation functions
  ////////////////////////////////////////////////////////////////////////

  // Range limits
  unsigned long alimits(const double u,
			const double varpi,
			const double varpi_t,
			const double epsabs,
			const double epsrel,
			const pwind *pw,
			double *alim);
  unsigned long xlimits(const double a,
			const pwind *pw,
			double *xlim);
  double amax(const double x,
	      const double epsabs,
	      const double epsrel,
	      const pwind *pw);
  unsigned long s_crit(const double varpi,
		       const double varpi_t,
		       const double u,
		       const pwind *pw,
		       double *s_crit_);
  unsigned long a_crit(const double varpi,
		       const double varpi_t,
		       const double u,
		       const pwind *pw,
		       double *a_crit_);

  // Wind density
  double drhodx(const double x,
		const double a,
		const pwind *pw);
  double rho(const double a,
	     const double epsabs,
	     const double epsrel,
	     const pwind *pw);
  void rho_vec(const unsigned long na,
	       const double *a,
	       const double epsabs,
	       const double epsrel,
	       const pwind *pw,
	       double *result);

  // Momentum flux
  double pdot(const double a,
	      const double epsabs,
	      const double epsrel,
	      const pwind *pw);
  void pdot_vec(const unsigned long na,
		const double *a,
		const double epsabs,
		const double epsrel,
		const pwind *pw,
		double *result);

  // Energy flux
  double Edot(const double a,
	      const double epsabs,
	      const double epsrel,
	      const pwind *pw);
  void Edot_vec(const unsigned long na,
		const double *a,
		const double epsabs,
		const double epsrel,
		const pwind *pw,
		double *result);

  // Momentum flux normalised to driving momentum flux
  double pdotRel_approx(const double a,
			const double epsabs,
			const double epsrel,
			const pwind *pw);
  double pdotRel_exact(const double a,
		       const double fg,
		       const double tctw,
		       const double epsabs,
		       const double epsrel,
		       const pwind *pw);
  
  // Absorption
  double Phi_uc(const double u,
		const double varpi,
		const double varpi_t,
		const double a0,
		const double a1,
		const double epsabs,
		const double epsrel,
		const pwind *pw);
  void Phi_uc_vec(const unsigned long nu,
		  const double *u,
		  const double varpi,
		  const double varpi_t,
		  const double a0,
		  const double a1,
		  const double epsabs,
		  const double epsrel,
		  const pwind *pw,
		  double *result);
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
		const pwind *pw);
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
		  double *result);
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
			 const pwind *pw);
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
			   double *result);
  double Phi_c(const double u,
	       const double fw,
	       const double varpi,
	       const double varpi_t,
	       const double a0,
	       const double a1,
	       const double epsabs,
	       const double epsrel,
	       const pwind *pw);
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
		 double *result);
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
	       const pwind *pw);
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
		 double *result);
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
			const pwind *pw);
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
			  double *result);

  // Thin, subcritical emission
  double Xi(const double u,
	    const double varpi,
	    const double varpi_t,
	    const double epsabs,
	    const double epsrel,
	    const pwind *pw);
  void Xi_vec(const unsigned long nu,
	      const double *u,
	      const double varpi,
	      const double varpi_t,
	      const double epsabs,
	      const double epsrel,
	      const pwind *pw,
	      double *result);
  double xi(const double varpi,
	    const double varpi_t,
	    const double epsabs,
	    const double epsrel,
	    const pwind *pw);

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
	     const pwind *pw);
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
	       double *result);
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
	     const pwind *pw);

  // Error management
  int get_err(const pwind *pw);
  void set_err(const int err, pwind *pw);
  void clear_err(pwind *pw);
  const char *get_err_str(const pwind *pw);

  // Functions to manage the hot data tabulated data
  hot_wind_table *read_hot_wind_table(const char *dirname,
				      const int yidx,
				      const int midx);
  void free_hot_wind_table(hot_wind_table *tab);
  void get_hot_wind_table_limits(hot_wind_table *tab,
				 double *uh_limits);
}
