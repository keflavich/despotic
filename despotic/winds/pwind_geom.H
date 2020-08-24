// These classes define geometric structures in the wind

#ifndef _PWIND_GEOM_H_
#define _PWIND_GEOM_H_

#include <cmath>
#include <vector>

class pwind_geom {
public:
  pwind_geom(const double phi_ = 0.0) : phi(phi_), mu(cos(phi_)) {}
  virtual ~pwind_geom() {}

  // Return angle phi or its cosine
  double get_phi() const { return phi; }
  double get_mu() const { return mu; }

  // Radii where line from infinity crosses wind boundary
  virtual std::vector<double>
  a_crit(const double varpi, const double varpi_t,
	 const double u=0.0) const;

  // Radial integration limits coming from geometric and velocity
  // considerations
  virtual std::vector<double>
  a_lim(const std::vector<double> a_lim_v, const double varpi,
	const double varpi_t, const double u=0.0) const;

  // Positions where line from infinity crosses wind boundary --
  // defined in derived class
  virtual std::vector<double>
  s_crit(const double varpi, const double varpi_t,
	 const double u=0.0) const = 0;

  // Geometric covering fraction -- defined in derived class
  virtual double fcov() const = 0;

protected:
  // Orientation of wind relative to line of sight; only used /
  // meaningful if geometry is cone or cone_sheath, or if wind
  // specific angular momentum jsp != 0
  double phi;
  double mu; // mu = cos(phi) -- cached for speed
};

// Geometries we know about

// A trivial case of a spherical wind
class pwind_geom_sphere : public pwind_geom {
public:
  pwind_geom_sphere(const double phi_ = 0.0) : pwind_geom(phi_) {}
  virtual ~pwind_geom_sphere() {}
  virtual std::vector<double> s_crit(const double varpi,
				     const double varpi_t,
				     const double u=0.0) const;
  virtual double fcov() const { return 1.0; }
};

// A biconical wind
class pwind_geom_cone : public pwind_geom {
public:
  pwind_geom_cone(const double theta_, const double phi_) :
    pwind_geom(phi_), theta(theta_) {}
  virtual ~pwind_geom_cone() {}
  virtual std::vector<double> s_crit(const double varpi,
				     const double varpi_t,
				     const double u=0.0) const;
  std::vector<double> s_intersect(const double varpi,
				  const double varpi_t) const;
  void s_clean(std::vector<double>& s, const double vp2,
	       const double u=0.0) const;
  virtual double fcov() const { return sin(theta); }
protected:
  double theta;
};

// A biconical sheath
class pwind_geom_cone_sheath : public pwind_geom_cone {
public:
  pwind_geom_cone_sheath(const double theta_out, const double theta_in,
			 const double phi_) :
    pwind_geom_cone(theta_out, phi_),
    inner_cone(theta_in, phi_) {}
  virtual ~pwind_geom_cone_sheath() {}
  virtual std::vector<double> s_crit(const double varpi,
				     const double varpi_t,
				     const double u=0.0) const;
  virtual double fcov() const { return sin(theta) - inner_cone.fcov(); }
protected:
  pwind_geom_cone inner_cone;
};

// Create a sphere object that we can use as a default
static const pwind_geom_sphere pw_geom_sphere;

#endif
// _PWIND_GEOM_H_
