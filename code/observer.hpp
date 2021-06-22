#pragma once

struct Observer
{
  Observer(double dt_write = 0.1, double tend = 1.0, double t0 = 0.0)
    : dt_write(dt_write), tend(tend), t0(t0) {}

  bool operator()(double t) const
  {
    return allowWrite && (t >= (t0 + write_idx * dt_write) - tol || t > tend - tol);
  }

  Observer& operator++()
  {
    write_idx++;
    return *this;
  }

  Observer operator++(int)
  {
    Observer cpy(*this);
    write_idx++;
    return cpy;
  }

  void enable() { allowWrite = true; }
  void disable() { allowWrite = false; }

  void setEndTime(double tend_) { tend = tend_; }
  void setStartTime(double t0_) { t0 = t0_; }
  void setTimestep(double dt_write_) { dt_write = dt_write_; }

  size_t get_write_idx() const  { return write_idx; }
  double getTimestep() const { return dt_write; }

private:
  double dt_write;
  double tend, t0;
  size_t write_idx = 1;
  double tol = 1.e-10;
  bool allowWrite = true;
};

template <class Real>
struct ScalarContainer
{
  ScalarContainer(Real value_) : value(value_) {}

  Real const& operator[](size_t) const { return value; }

private:
  Real value;
};
