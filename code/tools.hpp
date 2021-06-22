#pragma once

#include <cmath>

#include "Math.hpp"
#include "Shtns.hpp"


namespace Shtns
{

  template <class Real, class Container>
  Real integrate(shtns_cfg& shtns, Real R, Container const& psi, Real* weights, size_t& n_weights)
  {
    if (weights[0] == 0.0 || n_weights == 0) {
      n_weights = shtns_gauss_wts(shtns, weights);
      std::cout << "n_weights = " << n_weights << ", nlat = " << shtns->nlat << "\n";
    }

    Real dphi = 2.0*M_PI / shtns->nphi;

    Real result = 0.0;
    for (size_t i = 0; i < shtns->nlat; ++i) {
      Real sum = 0.0;
      for (size_t j = 0; j < shtns->nphi; ++j)
        sum += psi[get_idx(shtns->nlat, j, i)];

      sum *= weights[sym_idx(i, n_weights)]*dphi;
      result += sum;
    }

    return R*R*result;
  }

  template <class Real, class Container, class Functor>
  Real integrate(shtns_cfg& shtns, Real R, Container const& psi, Functor f, Real* weights, size_t& n_weights)
  {
    if (weights[0] == 0.0 || n_weights == 0) {
      n_weights = shtns_gauss_wts(shtns, weights);
      std::cout << "n_weights = " << n_weights << ", nlat = " << shtns->nlat << "\n";
    }

    Real dphi = 2.0*M_PI / shtns->nphi;

    Real result = 0.0;
    for (size_t i = 0; i < shtns->nlat; ++i) {
      Real sum = 0.0;
      for (size_t j = 0; j < shtns->nphi; ++j)
        sum += f(psi[get_idx(shtns->nlat, j, i)]);

      sum *= weights[sym_idx(i, n_weights)]*dphi;
      result += sum;
    }

    return R*R*result;
  }

  // template <class Real, class Complex, class Container2, class Functor>
  // Real integrate(shtns_cfg& shtns, Real R, VectorValues<Real,Complex> const& grad_psi, Container2 const& psi, Functor f, Real* weights, size_t& n_weights)
  // {
  //   if (weights[0] == 0.0 || n_weights == 0) {
  //     n_weights = shtns_gauss_wts(shtns, weights);
  //     std::cout << "n_weights = " << n_weights << ", nlat = " << shtns->nlat << "\n";
  //   }

  //   Real dphi = 2.0*M_PI / shtns->nphi;
  //   std::array<Real, 3> vec;

  //   Real result = 0.0;
  //   for (size_t i = 0; i < shtns->nlat; ++i) {
  //     Real sum = 0.0;
  //     for (size_t j = 0; j < shtns->nphi; ++j) {
  //       sphere_to_cartesian(shtns, grad_psi.values_t(), grad_psi.values_p(), j, i, vec);
  //       sum += f(vec, psi[get_idx(shtns->nlat, j, i)]);
  //     }

  //     sum *= weights[sym_idx(i, n_weights)]*dphi;
  //     result += sum;
  //   }

  //   return R*R*result;
  // }

} // end namespace Shtns
