#pragma once

#include <array>
#include <boost/mpl/int.hpp>

namespace Shtns
{
  using boost::mpl::int_;

  template <class T>
  T sqr(T x) { return x*x; }

  template <class T>
  T pow3(T x) { return x*x*x; }


  template <class T>
  T pow_aux(T x, int_<1>) { return x; }

  template <class T>
  T pow_aux(T x, int_<0>) { return T(1); }

  template <int p, class T>
  T pow_aux(T x, int_<p>) { return x * pow_aux(x, int_<p-1>()); }

  template <int p, class T>
  T pow(T x) { return pow_aux(x, int_<p>()); }

  //template <class T, std::size_t N>
  //T norm(std::array<T, N> const& x)
  //{
  //  return std::sqrt(norm2(x));
  //}

  template <class T>
  T norm2(std::array<T, 3> const& x)
  {
    return sqr(x[0]) + sqr(x[1]) + sqr(x[2]);
  }

  template <class T>
  T norm2(std::array<T, 2> const& x)
  {
    return sqr(x[0]) + sqr(x[1]);
  }

  template <class T>
  T norm2(std::array<T, 1> const& x)
  {
    return sqr(x[0]);
  }
  template <class T, std::size_t N>
  T norm(std::array<T, N> const& x)
  {
    return std::sqrt(norm2(x));
  }

  template <class T>
  std::array<T, 3>& operator-=(std::array<T, 3>& lhs, std::array<T, 3> const& rhs)
  {
    lhs[0] -= rhs[0]; lhs[1] -= rhs[1]; lhs[2] -= rhs[2];
    return lhs;
  }

  template <class T>
  std::array<T, 2>& operator-=(std::array<T, 2>& lhs, std::array<T, 2> const& rhs)
  {
    lhs[0] -= rhs[0]; lhs[1] -= rhs[1];
    return lhs;
  }

  template <class T>
  std::array<T, 1>& operator-=(std::array<T, 1>& lhs, std::array<T, 1> const& rhs)
  {
    lhs[0] -= rhs[0];
    return lhs;
  }

  template <class T, size_t N>
  std::array<T, N> operator-(std::array<T, N> lhs, std::array<T, N> const& rhs)
  {
    lhs-= rhs;
    return lhs;
  }

  template <class T>
  std::array<T, 3>& operator+=(std::array<T, 3>& lhs, std::array<T, 3> const& rhs)
  {
    lhs[0] += rhs[0]; lhs[1] += rhs[1]; lhs[2] += rhs[2];
    return lhs;
  }

  template <class T>
  std::array<T, 2>& operator+=(std::array<T, 2>& lhs, std::array<T, 2> const& rhs)
  {
    lhs[0] += rhs[0]; lhs[1] += rhs[1];
    return lhs;
  }

  template <class T>
  std::array<T, 1>& operator+=(std::array<T, 1>& lhs, std::array<T, 1> const& rhs)
  {
    lhs[0] += rhs[0];
    return lhs;
  }

  template <class T, size_t N>
  std::array<T, N> operator+(std::array<T, N> lhs, std::array<T, N> const& rhs)
  {
    lhs+= rhs;
    return lhs;
  }

  template <class T>
  std::array<T, 3>& operator*=(std::array<T, 3>& lhs, double factor)
  {
    lhs[0] *= factor; lhs[1] *= factor; lhs[2] *= factor;
    return lhs;
  }

  template <class T>
  std::array<T, 2>& operator*=(std::array<T, 2>& lhs, double factor)
  {
    lhs[0] *= factor; lhs[1] *= factor;
    return lhs;
  }

  template <class T>
  std::array<T, 1>& operator*=(std::array<T, 1>& lhs, double factor)
  {
    lhs[0] *= factor;
    return lhs;
  }

  template <class T, size_t N>
  std::array<T, N> operator*(std::array<T, N> lhs, double factor)
  {
    lhs*= factor;
    return lhs;
  }

  template <class T, size_t N>
  std::array<T, N> operator*(double factor, std::array<T, N> lhs)
  {
    lhs*= factor;
    return lhs;
  }

  template <class T>
  std::array<T, 3> cross(std::array<T, 3> const& a, std::array<T, 3> const& b)
  {
    return {{
      a[1]*b[2] - a[2]*b[1],
      a[2]*b[0] - a[0]*b[2],
      a[0]*b[1] - a[1]*b[0]}};
  }



  template <class T>
  T distance(std::array<T, 3> const& lhs, std::array<T, 3> const& rhs)
  {
    return std::sqrt(sqr(lhs[0] - rhs[0]) + sqr(lhs[1] - rhs[1]) + sqr(lhs[2] - rhs[2]));
  }

  template <class T>
  T distance(std::array<T, 2> const& lhs, std::array<T, 2> const& rhs)
  {
    return std::sqrt(sqr(lhs[0] - rhs[0]) + sqr(lhs[1] - rhs[1]));
  }


  template <class T>
  T dot(std::array<T, 3> const& lhs, std::array<T, 3> const& rhs)
  {
    return (lhs[0] * rhs[0]) + (lhs[1] * rhs[1]) + (lhs[2] * rhs[2]);
  }

  template <class T>
  T dot(std::array<T, 2> const& lhs, std::array<T, 2> const& rhs)
  {
    return (lhs[0] * rhs[0]) + (lhs[1] * rhs[1]);
  }

  template <class T>
  T sign(T x) { return x > 0 ? T(1) : (x < 0 ? T(-1) : T(0)); }

} // end namespace Shtns
