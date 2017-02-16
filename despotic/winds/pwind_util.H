/* This module defines some generic utilities that are used by all
   the other c routines. */

#ifndef _PWIND_UTIL_H_
#define _PWIND_UTIL_H_

#include <cmath>
#include <iostream>
#include <gsl/gsl_sf_erf.h>

#define SQR(a) ((a)*(a))
#define CUBE(a) ((a)*(a)*(a))
#define MAXINTERVAL 65536

static inline double pM(const double x, const double sx) {
  double sx2 = SQR(sx);
  return 1.0/(sqrt(2.0*M_PI)*sx) * exp(-SQR(x-sx2/2.0) / (2.0*sx2));
}

static inline double pA(const double x, const double sx) {
  double sx2 = SQR(sx);
  return 1.0/(sqrt(2.0*M_PI)*sx) * exp(-SQR(x+sx2/2.0) / (2.0*sx2));
}

static inline double zetaM(const double xcrit, const double sx) {
  return 0.5*(1.0-gsl_sf_erf((-2.0*xcrit+SQR(sx))/(2.*sqrt(2)*sx)));
}

static inline double zetaA(const double xcrit, const double sx) {
  return 0.5*(1.0+gsl_sf_erf((2.0*xcrit+SQR(sx))/(2.*sqrt(2)*sx)));
}

static inline double sxMach(const double mach) {
  double alpha=2.5;
  double rfac = 0.5*(3.0-alpha)/(2.0-alpha) *
    (1.0-pow(mach, 2*(2-alpha))) / (1.0-pow(mach, 2*(3-alpha)));
  return sqrt(log(1.0+rfac*SQR(mach)/4.));
}

#endif
// _PWIND_UTIL_H_
