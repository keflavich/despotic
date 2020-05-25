#include "pwind_geom.H"
#include "pwind_util.h"
#include <algorithm>
#include <limits>

using namespace std;

////////////////////////////////////////////////////////////////////////
// General for all geometries
////////////////////////////////////////////////////////////////////////

std::vector<double>
pwind_geom::a_crit(const double varpi, const double varpi_t,
		   const double u) const {

  // Get s values
  vector<double> s = s_crit(varpi, varpi_t, u);
  if (s.size() == 0) return s;

  // Get a
  vector<double> a(s.size());
  for (vector<double>::size_type i=0; i<s.size(); i++)
    a[i]= sqrt(SQR(s[i]) + SQR(varpi) + SQR(varpi_t));

  // Swap pairs of elements if necessary to ensure radii are increasing
  for (vector<double>::size_type i=0; i<a.size()/2; i++) {
    if (a[2*i+1] < a[2*i]) {
      double atmp = a[2*i+1];
      a[2*i+1] = a[2*i];
      a[2*i] = atmp;
    }
  }

  // Return a
  return a;
}

vector<double>
pwind_geom::a_lim(const vector<double> a_lim_v, const double varpi,
		  const double varpi_t, const double u) const {
  
  // Get geometric limits along line of sight
  vector<double> a_geom = a_crit(varpi, varpi_t, u);

  // Construct intersection of geometric and velocity limits
  vector<double> a_lim;
  for (vector<double>::size_type i=0; i<a_geom.size(); i+=2) {
    for (vector<double>::size_type j=0; j<a_lim_v.size(); j+=2) {
      double a_lo = max(a_geom[i], a_lim_v[j]);
      double a_hi = min(a_geom[i+1], a_lim_v[j+1]);
      if (a_lo < a_hi) {
	a_lim.push_back(a_lo);
	a_lim.push_back(a_hi);
      }
    }
  }

  // Return
  return a_lim;
}

////////////////////////////////////////////////////////////////////////
// Spherical geometry
////////////////////////////////////////////////////////////////////////

vector<double>
pwind_geom_sphere::s_crit(const double varpi,
			  const double varpi_t,
			  const double u) const {
  vector<double> s;
  if (u <= 0.0) {
    s.push_back(-numeric_limits<double>::max());
    s.push_back(0.0);
  }
  if (u >= 0.0) {
    s.push_back(0.0);
    s.push_back(-numeric_limits<double>::max());
  }
  return s;
}

////////////////////////////////////////////////////////////////////////
// Conical geometry
////////////////////////////////////////////////////////////////////////

// This routine finds the intersection of a line with a cone. The cone
// is tilted at an angle phi relative to the z axis, with the cone's
// central axis constrained to lie in the xz-plane. The line is
// horizontal, described by y = constant and z = constant. To find the
// intersection, we proceed as follows:
//
// 1. We rotate the coordinate system to align the cone's central axis
// with the z axis. To do this, we rotate by an angle -phi about a unit
// vector in the \hat{y} direction. Using Rodrigues' rotation formula,
// this causes an arbitrary point (x, y, z) to be rotated to
// (x cos phi - z sin phi, y, z cos phi - x sin phi).
//
// 2. The equation of the line in this rotated coordinate system is
// then given by
// x = s cos phi - z0 sin phi
// y = y0
// z = z0 cos phi + s sin phi
// where y = y0, z = z0 were the equations of the line in the original
// coordinate system, and s measures position along the line.
//
// 3. In the rotated coordiante system where the cone's central axis
// is \hat{z}, the cone is described by the equation
// x^2 + y^2 - z^2 tan^2 theta = 0,
// which we must solve simultaneously with the equation of the line
// derived in step 2.
//
// 4. Inserting the equation of the line into the equation of the cone
// and solving for the position along the line s gives
// s = +- [sqrt(2) cos(theta) / (cos(2 theta) + cos(2 phi))] *
//    { sqrt[z^2 - (y^2+z^2) cos (2 theta) - y^2 cos(2 phi)] +
//      z sin(2 phi) }
// Thus a solution exists only if the term under the square root,
// det = z^2 - (y^2+z^2) cos (2 theta) - y^2 cos(2 phi) >= 0.
//
// 5. If a solution does exist, the intersection points (x_r, y_r,
// z_r) in the rotated coordinate system can be found by inserting the
// values of s obtaind in step 4 into the equation of the line found
// in step 2. The points in the original coordinate system can then be
// found by performing the inverse rotation, giving
// x = x_r cos(phi) + z_r sin(phi)
// y = y0
// z = z_r cos(phi) - x_r sin(phi)
// However, this simplifies to simply the statement that x = s, so the
// solution found in step 4 gives the x coordinate of the
// intersection, which is what we are after.

vector<double>
pwind_geom_cone::s_intersect(const double varpi,
			     const double varpi_t) const {
  const double vp2 = SQR(varpi) + SQR(varpi_t);
  vector<double> s;

  // Compute the determinant, and if it is <= 0 (indicating no
  // solution), return an empty vector
  const double det = SQR(varpi) - vp2*cos(2*theta)
    - SQR(varpi_t)*cos(2*phi);
  if (det <= 0) return s;

  // Find the points where the line of sight intersects the wind cone;
  // if one or both of them fall within the sphere x^2 + y^2 + z^2 < 1
  // where the wind is absent, displace the intersection point onto
  // the surface of the sphere by moving along the line of sight
  const double denom = cos(2*theta) + cos(2*phi);
  s.push_back((-cos(theta) * sqrt(2*det) + varpi*sin(2*phi)) / denom);
  s.push_back(( cos(theta) * sqrt(2*det) + varpi*sin(2*phi)) / denom);

  // If s = -+infinity is within the wind cone, add those points
  if (phi - theta <= -M_PI/2 || phi + theta >= M_PI/2) {
    s.push_back(-numeric_limits<double>::max());
    s.push_back(numeric_limits<double>::max());
  }

  // If s = 0 is within the wind cone, add two copies of it. We check
  // if s = 0 is within the wind cone by computing the dot product of
  // (0, y0, z0) with the unit vector along with wind central axis,
  // (sin phi, 0, cos phi). This gives the angle phi' between the wind
  // central axis and the s = 0 point; s = 0 is within the wind if
  // this angle is <= theta.
  const double phi_prime = acos(varpi*cos(phi) / sqrt(vp2));
  if (phi_prime < theta) {
    s.push_back(0.0);
    s.push_back(0.0);
  }

  // If the line of sight hits the sphere where the wind ends, add
  // those intersection points too
  if (vp2 < 1.0) {
    s.push_back(-sqrt(1.0-vp2));
    s.push_back(sqrt(1.0-vp2));
  }

  // Sort the various intersection points we've identified
  sort(s.begin(), s.end());

  // Return intersection points
  return s;
}

void
pwind_geom_cone::s_clean(vector<double>& s, const double vp2,
			 const double u) const {

  // Walk through the intersection points, getting rid of any
  // intervals that lie within the unit sphere
  for (long i=s.size()-2; i>=0; i-=2) {
    if (SQR(s[i]) + vp2 <= 1.0 && SQR(s[i+1]) + vp2 <= 1.0) {
      for (vector<double>::size_type j=i; j<s.size()-2; j++) {
	s[j] = s[j+2];
	s[j+1] = s[j+3];
      }
      s.resize(s.size()-2);
    }
  }

  // If we are looking only for intersection points at u < 0 or u > 0,
  // further truncate the list
  vector<double> s_return;
  for (vector<double>::size_type i=0; i<s.size(); i+=2) {
    if (u <= 0 && s[i] < 0.0) {
      s_return.push_back(s[i]);
      s_return.push_back(s[i+1]);
    } else if (u >= 0 && s[i+1] > 0.0) {
      s_return.push_back(s[i]);
      s_return.push_back(s[i+1]);
    }
  }

  // Store final result back to original s
  s = s_return;
}  

vector<double>
pwind_geom_cone::s_crit(const double varpi, const double varpi_t,
			const double u) const {
  vector<double> s = s_intersect(varpi, varpi_t);
  s_clean(s, SQR(varpi)+SQR(varpi_t), u);
  return s;
}

////////////////////////////////////////////////////////////////////////
// Conical sheath geometry
////////////////////////////////////////////////////////////////////////
vector<double>
pwind_geom_cone_sheath::
s_crit(const double varpi, const double varpi_t, const double u) const {
  
  // Get points of impact on inner and outer cones
  vector<double> s_outer = s_intersect(varpi, varpi_t);
  vector<double> s_inner = inner_cone.s_intersect(varpi, varpi_t);

  // Merge the two lists; to do this, we walk through the two lists,
  // keeping track of when we enter each cone
  vector<double> s;
  vector<double>::size_type outptr=0, inptr=0;
  while (outptr < s_outer.size()) {
    if (inptr < s_inner.size()) {
      // We have not yet passed the last point where we hit the inner
      // cone, so check if we hit the inner cone or outer cone next
      if (s_inner[inptr] <= s_outer[outptr]) {
	// We hit the inner cone next
	s.push_back(s_inner[inptr]);
	inptr++;
	continue;
      }
    }
    // If we are here, we hit the outer cone next; record the hit
    s.push_back(s_outer[outptr]);
    outptr++;
  }
  
  // Clean
  s_clean(s, SQR(varpi)+SQR(varpi_t), u);

  // Return
  return s;
}

