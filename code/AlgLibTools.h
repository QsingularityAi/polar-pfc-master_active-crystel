#pragma once

#include "interpolation.h" // from alglib
#include "kdtree_nanoflann.h"
#include "Math.hpp"

namespace Shtns {

namespace pfc_tools {

  #define DOW 3
  using WorldVector = std::array<double,DOW>;

  /// ALGLIB routines
  ///_______________________________________________________________________________________________

  inline void vector2real_1d_array(std::vector<double> const& dofvector, alglib::real_1d_array& array)
  {
    array.setlength(dofvector.size());
    for(size_t i = 0; i < dofvector.size(); ++i)
    {
      array(i) = dofvector[i];
    }
  }

  inline void vector2real_2d_array(std::vector<WorldVector> const& dofvector, alglib::real_2d_array& array)
  {
    array.setlength(dofvector.size(), DOW);
    for(size_t i = 0; i < dofvector.size(); ++i)
    {
      for (int j = 0; j < DOW; ++j)
	      array(i,j) = dofvector[i][j];
    }
  }

  inline void vector2real_2d_array(std::vector<double> const& dofvector, alglib::real_2d_array& array)
  {
    array.setlength(dofvector.size(), 1);
    for(size_t i = 0; i < dofvector.size(); ++i)
	    array(i,0) = dofvector[i];
  }

  // ---------------------------------------------------------------------------


  inline void function_pfc(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, void *ptr)
  {
    auto x_trans = std::sqrt(sqr(x[0]-c[3]) + sqr(x[1]-c[4]) + sqr(x[2]-c[5]));
    func = c[2]*x_trans <=-M_PI ? -c[1]+c[0] :
           c[2]*x_trans >= M_PI ?  c[1]+c[0] :
           c[0] + c[1]*std::cos(c[2]*x_trans);
  }

  inline void gradient_pfc(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr)
  {
    auto x_trans = std::sqrt(sqr(x[0]-c[3]) + sqr(x[1]-c[4]) + sqr(x[2]-c[5]));
    func = c[2]*x_trans <=-M_PI ? -c[1]+c[0] :
           c[2]*x_trans >= M_PI ?  c[1]+c[0] :
           c[0] + c[1]*std::cos(c[2]*std::sqrt(sqr(x[0]-c[3]) + sqr(x[1]-c[4]) +sqr(x[2]-c[5])));

    grad[0] = 1.0;
    grad[1] = c[2]*x_trans <=-M_PI ? -1.0 :
              c[2]*x_trans >= M_PI ?  1.0 : std::cos(x_trans);

    if (c[2]*x_trans <=-M_PI || c[2]*x_trans >= M_PI) {
      grad[2] = 0.0;
      grad[3] = 0.0;
      grad[4] = 0.0;
      grad[5] = 0.0;
    } else {
      auto tmp = -c[1]*std::sin(c[2]*x_trans);
      grad[2] = tmp*x_trans;
      grad[3] = tmp*c[2]*(x[0]-c[3])/x_trans;
      grad[4] = tmp*c[2]*(x[1]-c[4])/x_trans;
      grad[5] = tmp*c[2]*(x[2]-c[5])/x_trans;
    }
  }

  // ---------------------------------------------------------------------------

  inline void function_pfc2(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, void *ptr)
  {
    auto x_nrm2 = sqr(x[0]-c[3]) + sqr(x[1]-c[4]) + sqr(x[2]-c[5]);
    func = c[0] + c[1] * std::exp( -x_nrm2 / sqr(c[2]) );
  }

  inline void gradient_pfc2(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr)
  {
    auto x_nrm2 = sqr(x[0]-c[3]) + sqr(x[1]-c[4]) + sqr(x[2]-c[5]);
    auto f0 = std::exp( -x_nrm2 / sqr(c[2]) );

    func = c[0] + c[1] * f0;

    grad[0] = 1.0;
    grad[1] = f0;
    grad[2] = 2*c[1]/sqr(c[2]) * f0 * x_nrm2/c[2];
    grad[3] = 2*c[1]/sqr(c[2]) * f0 * (x[0]-c[3]);
    grad[4] = 2*c[1]/sqr(c[2]) * f0 * (x[1]-c[4]);
    grad[5] = 2*c[1]/sqr(c[2]) * f0 * (x[2]-c[5]);
  }

  // ------------------------------------------------------------------------------

  inline void nonlinearPtcFit(std::vector<WorldVector> const& points, std::vector<double> const& values, WorldVector &p, double min_psi, double max_psi)
  {
    using namespace alglib;

    real_2d_array x;
    vector2real_2d_array(points, x);

    real_1d_array y;
    vector2real_1d_array(values, y);

    real_1d_array c; c.setlength(6);
    c(0) = 0.5 * (max_psi + min_psi);
    c(1) = 0.5 * (max_psi - min_psi);
    c(2) = 1;
    c(3) = p[0];
    c(4) = p[1];
    c(5) = p[2];

    double d = 4*M_PI/std::sqrt(3.0);

    real_1d_array lb; lb.setlength(6);
    lb(0) = c(0)-1.0;
    lb(1) = c(1)/2.0;
    lb(2) = 0.0;
    lb(3) = p[0] - d/4.0;
    lb(4) = p[1] - d/4.0;
    lb(5) = p[2] - d/4.0;

    real_1d_array ub; ub.setlength(6);
    ub(0) = c(0)+1.0;
    ub(1) = c(1)*2.0;
    ub(2) = std::numeric_limits<double>::infinity();
    ub(3) = p[0] + d/4.0;
    ub(4) = p[1] + d/4.0;
    ub(5) = p[2] + d/4.0;

    double epsf = 0.0;
    double epsx = 0.0;
    ae_int_t maxits = 0;
    ae_int_t info;
    lsfitstate state;
    lsfitreport rep;
    double diffstep = 1.e-6;

    lsfitcreatefg(x, y, c, points.size(), DOW, 6, true, state);
    lsfitsetcond(state, epsx, maxits);
    lsfitsetbc(state, lb, ub);
    alglib::lsfitfit(state, pfc_tools::function_pfc2, pfc_tools::gradient_pfc2);
    lsfitresults(state, info, c, rep);

    WorldVector p0 = p;
    p[0] = c(3); p[1] = c(4); p[2] = c(5);
  }
} // end namespace

}
