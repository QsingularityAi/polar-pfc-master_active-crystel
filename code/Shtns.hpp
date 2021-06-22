#pragma once

#include <complex>
#include <fftw3.h>
#include <shtns.h>

namespace Shtns
{

  inline size_t sym_idx(size_t i, size_t N)
  {
    return i < N ? i : (2*N - i - 1);
  }


  inline size_t get_idx(size_t nlat, size_t i/*phi*/, size_t j/*theta*/)
  {
    return i*nlat + j;
  }

  template <class T>
  T _sqr(T x) { return x*x; }

  template <class Real>
  inline Real sphere_dot(shtns_cfg const& shtns, Real Vt, Real Vp, Real Wt, Real Wp, size_t j/*theta*/, Real radius)
  {
    return _sqr(radius)*Vt*Wt + _sqr(radius*(shtns->st[j]))*Vp*Wp;
  }

  template <class Real>
  inline Real sphere_norm2(shtns_cfg const& shtns, Real Vt, Real Vp, size_t j/*theta*/, Real radius)
  {
    return _sqr(radius*Vt) + _sqr(radius*(shtns->st[j])*Vp);
  }

  template <class Real>
  inline Real sphere_norm(shtns_cfg const& shtns, Real Vt, Real Vp, size_t j/*theta*/, Real radius)
  {
    return std::sqrt(sphere_norm2(shtns, Vt, Vp, j, radius));
  }


  template <class Real, class Array>
  inline size_t sphere_to_cartesian_orth(shtns_cfg const& shtns, Real const* Vt, Real const* Vp, size_t i/*phi*/, size_t j/*theta*/, Array& out)
  {
    Real phi = i * 2.0*M_PI/shtns->mres / shtns->nphi;
    Real cos_theta = shtns->ct[j];
    Real sin_theta = shtns->st[j];
    size_t idx = get_idx(shtns->nlat, i, j);

    auto sin_theta_relaxes = (sin_theta > 0 ? 1 : -1) * std::max(1.e-8, std::abs(sin_theta));

    auto V1 = -Vp[idx]/sin_theta_relaxes;
    auto V2 = sin_theta * Vt[idx];

    out[0] = -V2*std::sin(phi) + V1*cos_theta*std::cos(phi);
    out[1] = V2*std::cos(phi) + V1*std::sin(phi)*cos_theta;
    out[2] = -V1*sin_theta;

    return idx;
  }

  template <class Real, class Array>
  inline size_t sphere_to_cartesian(shtns_cfg const& shtns, Real const* Vt, Real const* Vp, size_t i/*phi*/, size_t j/*theta*/, Array& out)
  {
    Real phi = i * 2.0*M_PI/shtns->mres / shtns->nphi;
    Real cos_theta = shtns->ct[j];
    Real sin_theta = shtns->st[j];
    size_t idx = get_idx(shtns->nlat, i, j);

    out[0] = -Vp[idx]*std::sin(phi) + Vt[idx]*cos_theta*std::cos(phi);
    out[1] = Vp[idx]*std::cos(phi) + Vt[idx]*std::sin(phi)*cos_theta;
    out[2] = -Vt[idx]*sin_theta;

    return idx;
  }

  template <class Real, class Array>
  inline void cartesian_to_sphere(shtns_cfg const& shtns, Array const& in, Real* Vt, Real* Vp, size_t i/*phi*/, size_t j/*theta*/)
  {
    Real phi = i * 2.0*M_PI/shtns->mres / shtns->nphi;
    Real cos_theta = shtns->ct[j];
    Real sin_theta = shtns->st[j];
    size_t idx = get_idx(shtns->nlat, i, j);

    Vt[idx] = cos_theta*(std::cos(phi)*in[0] + std::sin(phi)*in[1]) - sin_theta*in[2];
    Vp[idx] = -std::sin(phi)*in[0] + std::cos(phi)*in[1];
  }


  template <class Real, class Array, class Array2>
  inline void cartesian_to_sphere(shtns_cfg const& shtns, Array const& in, Array2& out, size_t i/*phi*/, size_t j/*theta*/)
  {
    Real phi = i * 2.0*M_PI/shtns->mres / shtns->nphi;
    Real cos_theta = shtns->ct[j];
    Real sin_theta = shtns->st[j];

    out[0] = cos_theta*(std::cos(phi)*in[0] + std::sin(phi)*in[1]) - sin_theta*in[2];
    out[1] = -std::sin(phi)*in[0] + std::cos(phi)*in[1];
  }


  template <class Real, class Array>
  inline void get_coord(shtns_cfg const& shtns, size_t i/*phi*/, size_t j/*theta*/, Real R, Array& x)
  {
    Real phi = i * 2.0*M_PI/shtns->mres / shtns->nphi;
    x[0] = R * shtns->st[j] * std::cos(phi);
    x[1] = R * shtns->st[j] * std::sin(phi);
    x[2] = R * shtns->ct[j];
  }

} // end namespace Shtns
