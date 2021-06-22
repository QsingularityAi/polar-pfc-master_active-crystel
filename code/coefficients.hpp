#pragma once

#include "flags.hpp"

const Flag INIT_LAPLACE    = 1 << 0;
const Flag INIT_DIVERGENCE = 1 << 1;
const Flag INIT_GRADIENT   = 1 << 2;
const Flag INIT_CURL       = 1 << 3;
const Flag INIT_CURLN      = 1 << 4;
const Flag INIT_CURLINV    = 1 << 5;

const Flag INIT_ALL = INIT_LAPLACE | INIT_DIVERGENCE | 
		      INIT_GRADIENT | INIT_CURL |
		      INIT_CURLN | INIT_CURLINV ;
		      
template <class Real>
struct ConstantCoeff
{
  ConstantCoeff(Real value) : value(value) {}
  inline Real operator[](size_t l) { return value; }
private:
  Real value;
};
