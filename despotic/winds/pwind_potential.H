// Classes that describe potentials for winds

#ifndef _PWIND_POTENTIAL_H_
#define _PWIND_POTENTIAL_H_

#include <vector>
#include <iostream>

class pwind_potential {
public:
  pwind_potential() {};
  ~pwind_potential() {};
  virtual double m(const double a) const = 0;
  virtual std::vector<double> m(const std::vector<double> a) const = 0;
  virtual double dmda(const double a) const = 0;
  virtual std::vector<double> dmda(const std::vector<double> a) const = 0;
  virtual int midx() const = 0;
};

// Child classes
class pwind_point : public pwind_potential {
public :
  pwind_point() {};
  ~pwind_point() {};
  virtual double m(const double a) const { return 1.0; }
  virtual std::vector<double> m(const std::vector<double> a) const {
    std::vector<double> m_(a.size());
    for (std::vector<double>::size_type i=0; i<a.size(); i++) m_[i] = 1.0;
    return m_;
  }
  virtual double dmda(const double a) const { return 0.0; }
  virtual std::vector<double> dmda(const std::vector<double> a) const {
    std::vector<double> dmda_(a.size());
    for (std::vector<double>::size_type i=0; i<a.size(); i++) dmda_[i] = 0.0;
    return dmda_;
  }
  virtual int midx() const { return 0; }
};

class pwind_isothermal : public pwind_potential {
public :
  pwind_isothermal() {};
  ~pwind_isothermal() {};
  virtual double m(const double a) const { return a; }
  virtual std::vector<double> m(const std::vector<double> a) const {
    std::vector<double> m_ = a;
    return m_;
  }
  virtual double dmda(const double a) const { return 1.0; }
  virtual std::vector<double> dmda(const std::vector<double> a) const {
    std::vector<double> dmda_(a.size());
    for (std::vector<double>::size_type i=0; i<a.size(); i++) dmda_[i] = 1.0;
    return dmda_;
  }
  virtual int midx() const { return 1; }
};

// Create static instances of these classes to that they can be used
// elsewhere
static const pwind_point pwind_potential_point;
static const pwind_isothermal pwind_potential_isothermal;

#endif
// _PWIND_POTENTIAL_H_
