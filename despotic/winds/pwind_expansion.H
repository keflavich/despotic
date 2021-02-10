// Classes that describe how clouds expand in a wind

#ifndef _PWIND_EXPANSION_H_
#define _PWIND_EXPANSION_H_

#include <vector>

// Pure virtual base class that defines cloud expansion
class pwind_expansion {
public:
  pwind_expansion() {}
  ~pwind_expansion() {}
  virtual double y(const double a) const = 0;
  virtual std::vector<double> y(const std::vector<double> a) const = 0;
  virtual double dyda(const double a) const = 0;
  virtual std::vector<double> dyda(const std::vector<double> a) const = 0;
  virtual bool const_sa() const = 0;
  virtual int yidx() const = 0;
};

// Derived classes
class pwind_area : public pwind_expansion {
public:
  pwind_area() {}
  ~pwind_area() {}
  virtual double y(const double a) const { return 1.0; }
  virtual std::vector<double> y(const std::vector<double> a) const {
    std::vector<double> y_(a.size());
    for (std::vector<double>::size_type i=0; i<a.size(); i++) y_[i] = 1.0;
    return y_;
  }
  virtual double dyda(const double a) const { return 0.0; }
  virtual std::vector<double> dyda(const std::vector<double> a) const {
    std::vector<double> dyda_(a.size());
    for (std::vector<double>::size_type i=0; i<a.size(); i++) dyda_[i] = 0.0;
    return dyda_;
  }
  virtual bool const_sa() const { return false; }
  virtual int yidx() const { return 0; }
};

class pwind_intermediate : public pwind_expansion {
public:
  pwind_intermediate() {}
  ~pwind_intermediate() {}
  virtual double y(const double a) const { return a; }
  virtual std::vector<double> y(const std::vector<double> a) const {
    std::vector<double> y_ = a;
    return y_;
  }
  virtual double dyda(const double a) const { return 1.0; }
  virtual std::vector<double> dyda(const std::vector<double> a) const {
    std::vector<double> dyda_(a.size());
    for (std::vector<double>::size_type i=0; i<a.size(); i++) dyda_[i] = 1.0;
    return dyda_;
  }
  virtual bool const_sa() const { return false; }
  virtual int yidx() const { return 1; }
};

class pwind_solid_angle : public pwind_expansion {
public:
  pwind_solid_angle() {}
  ~pwind_solid_angle() {}
  virtual double y(const double a) const { return a*a; }
  virtual std::vector<double> y(const std::vector<double> a) const {
    std::vector<double> y_(a.size());
    for (std::vector<double>::size_type i=0; i<a.size(); i++)
      y_[i] = a[i]*a[i];
    return y_;
  }
  virtual double dyda(const double a) const { return 2.0*a; }
  virtual std::vector<double> dyda(const std::vector<double> a) const {
    std::vector<double> dyda_(a.size());
    for (std::vector<double>::size_type i=0; i<a.size(); i++)
      dyda_[i] = 2.0*a[i];
    return dyda_;
  }
  virtual bool const_sa() const { return true; }
  virtual int yidx() const { return 2; }
};

// Create static instances of these classes to that they can be used
// elsewhere
static const pwind_area pwind_expansion_area;
static const pwind_intermediate pwind_expansion_intermediate;
static const pwind_solid_angle pwind_expansion_solid_angle;

#endif
// _PWIND_EXPANSION_H_
