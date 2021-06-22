#pragma once

#include "tools.hpp"

template <class Real>
struct InitialPfcTriangular
{
  InitialPfcTriangular(Real A_, Real rho0_, std::array<Real, 3>& x0_, Real radius_ = 10.0, Real q_ = 1.0)
  : A(A_),
    rho0(rho0_),
    x0(x0_),
    radius(radius_),
    q(q_)
  {
    idx0 = 1;
    idx1 = 2;
    d = 4.0*M_PI/std::sqrt(3.0);
  }

  Real operator()(const std::array<Real, 3>& x) const
  {
    Real dist = std::sqrt(sqr(x[0]-x0[0]) + sqr(x[1]-x0[1]) + sqr(x[2]-x0[2]));

    return dist < radius ? rho0 + A*(2.0*cos(q*(x[idx0] - x0[idx0]))
				- 4.0*cos(q*sqrt(3.0)/2.0 * (x[idx1] - x0[idx1] - d/2.0))*cos(q/2.0 * (x[idx0] - x0[idx0])))
			 : rho0;
  }

protected:
  Real A;
  Real rho0;
  std::array<Real, 3> x0;
  Real radius;
  Real q;

  size_t idx0, idx1;
  Real d;
};


struct InitialRandom
{
  InitialRandom(long int seed = std::time(0))
  {
    std::srand(seed);
  }

  template <class Real>
  Real operator()(const Real& amplitude, const Real& mean) const
  {
    return amplitude * 2.0 * (std::rand() / (Real)(RAND_MAX) - 0.5) + mean;
  }
};
