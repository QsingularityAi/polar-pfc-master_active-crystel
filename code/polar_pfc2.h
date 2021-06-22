#pragma once

#include <cmath>
#include <iostream>
#include <complex>
#include <memory>

#include "ProblemSHTNS.hpp"
#include "initial.hpp"
#include "ParticleWriter.hh"
#include "utility/Filesystem.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

// #define USE_C4
// #define USE_EXPLICIT_TRANSPORT

typedef double Real;
typedef std::complex<double> Complex;

typedef Eigen::Matrix<Complex, 3, 3>  EigenMatrix;
typedef Eigen::Matrix<Real, 3, 3>     EigenMatrixReal;
typedef Eigen::Matrix<Complex, 3, 1>  EigenVector;



template <class R>
std::array<R,3> sp2cart(std::array<R,2> const& p, R radius)
{
  return {{
    radius * std::sin(p[1]) * std::cos(p[0]),
    radius * std::sin(p[1]) * std::sin(p[0]),
    radius * std::cos(p[1])
  }};
}

using namespace Shtns;

struct PolarPfc : ProblemSHTNS<Real, Complex>
{
  typedef ProblemSHTNS<Real, Complex> super;

  PolarPfc(std::string name = "polar_pfc_shtns")
      : super(name, 2), psi3(NULL), P3(NULL)
  {
    eps =     pt.get("pfc.r", -0.98);
    density = pt.get("pfc.rho", -0.4);
    C1 =      pt.get("pfc.C1", 0.2);
    C4 =      pt.get("pfc.C4", 0.0);
    v0 =      pt.get("pfc.v0", 0.5);
    Dr =      pt.get("pfc.Dr", 0.5);

    read_from_file = pt.get(name + ".read_from_file", false);
    read_fn = pt.get<std::string>(name + ".read_filename", "");
    write_fn = pt.get<std::string>(name + ".write_filename", "");

    #ifndef USE_C4
    if (std::abs(C4) > 1.e-10)
      std::cerr << "C4 > 0 but !USE_C4\n";
    #endif

    // write statistics only after timestep dt2
    double dt_write3 = pt.get("output.dt3", 1.0);
    observer3.setTimestep(dt_write3);
    observer3.setEndTime(end_time);
    observer3.setStartTime(start_time);

    bool enableObserver3 = pt.get("output.enable3", true);
    if (!enableObserver3)
      observer3.disable();

    particle_writer.reset(new ParticleWriter<Real, Complex>(this->workspace, this->radius));

    std::string out_dir = "polar_pfc/output";
    std::string filename2 = pt.get("output.dir", out_dir) + "/order_parameters.csv";
    if (filesystem::exists(filename2))
      std::remove(filename2.c_str());
  }

  /// Destructor: free allocated data
  ~PolarPfc()
  {
    if (psi3) {
      delete psi3; psi3 = NULL;

      fftw_free(psi3_old);
    }
    if (P3) {
      delete P3; P3 = NULL;

      fftw_free(P3_old_x);
      fftw_free(P3_old_y);
    }
  }

  /// allocate data for nonlinear terms and corresponding old solutions
  virtual void initData(Flag init_flag = INIT_ALL) override
  {
    super::initData(init_flag);

    psi3 = new ScalarValues<Real, Complex>("psi3", workspace);
    psi3_old = (Complex*) fftw_malloc(super::size_dual * sizeof(Complex));

#ifdef USE_C4
      P3 = new VectorValues<Real, Complex>("P3", workspace);
      P3_old_x = (Complex*) fftw_malloc(super::size_dual * sizeof(Complex));
      P3_old_y = (Complex*) fftw_malloc(super::size_dual * sizeof(Complex));
#endif

    weights = std::make_shared<ScalarValues<Real, Complex>>("weights", workspace, false);
    energy = std::make_shared<ScalarValues<Real, Complex>>("energy", workspace);
    energyPsi = std::make_shared<ScalarValues<Real, Complex>>("energyPsi", workspace, false);
    energyP = std::make_shared<ScalarValues<Real, Complex>>("energyP", workspace, false);

    createMarix();
  }

  /// precalculate coefficients
  void createMarix()
  {
    mat.resize(workspace->lmax + 1);
#ifdef USE_EXPLICIT_TRANSPORT
    mat2.resize(workspace->lmax + 1);
#else
    mat3.resize(workspace->lmax + 1);
#endif

    double tau = timestep;
    for (long int l = 0; l <= workspace->lmax; ++l) {
      double alpha  = 1.0 -     tau*laplace[l]*(eps + sqr(1.0 + laplace[l]));
      double alpha2 = 3.0 - 2.0*tau*laplace[l]*(eps + sqr(1.0 + laplace[l]));
      double beta   = 1.0 -     tau*C1*(laplace[l] - Dr);
      double beta2  = 3.0 - 2.0*tau*C1*(laplace[l] - Dr);

      double t0 = tau*v0*div[l];
      double t1 = tau*v0*grad[l];

      EigenMatrixReal B;
      B(0,0) = alpha; B(0,1) = t0;   B(0,2) = 0.0;
      B(1,0) = t1;    B(1,1) = beta; B(1,2) = 0.0;
      B(2,0) = 0.0;   B(2,1) = 0.0;  B(2,2) = beta;
      mat[l].compute(B);

#ifdef USE_EXPLICIT_TRANSPORT
      B(0.0) = alpha2;  B(0,1) = 0.0;   B(0,2) = 0.0;
      B(1,0) = 0.0;     B(1,1) = beta2; B(1,2) = 0.0;
      B(2,0) = 0.0;     B(2,1) = 0.0;   B(2,2) = beta2;
      mat2[l].compute(B);
#else
      B(0.0) = alpha2; B(0,1) = 2.0*t0; B(0,2) = 0.0;
      B(1,0) = 2.0*t1; B(1,1) = beta2;  B(1,2) = 0.0;
      B(2,0) = 0.0;    B(2,1) = 0.0;    B(2,2) = beta2;
      mat3[l].compute(B);
#endif
    }
  }

  /// set the initial solution for psi and P,
  /// and set the nonlinear terms to zero
  virtual void solveInitialProblem() override
  {

    Real* psi = solution[0]->values();
    Real* Px = solution_vec[0]->values_t();
    Real* Py = solution_vec[0]->values_p();

    double A = 0.5;
    double rho0 = density;
    double d = 4.0*M_PI/std::sqrt(3.0);


    int initial = pt.get("pfc.initial", 0);
    double cut = pt.get("pfc.cut", 1.5);

    if (initial == 0) { // random noise
      // fill u0 and f
      InitialRandom fct;
      for (size_t i = 0; i < super::size; ++i) {
        psi[i] = fct(A, rho0);
      }
    }
    else if (initial == 1) { // triangular phase at poles
      double q = 1.0;
      auto fct = [A,rho0,q](std::array<Real, 3> const& x) {
        return rho0 + A*(2.0*std::cos(q*x[0]) - 4.0*std::cos(q*(std::sqrt(3.0)/2.0 * x[1] - M_PI))*std::cos(q/2.0 * x[0]));
      };

      std::array<Real, 3> X;
      for (size_t i = 0; i < workspace->nphi; ++i) {
        for (size_t j = 0; j < workspace->nlat; ++j) {
          get_coord(workspace, i, j, radius, X);
          size_t idx = get_idx(workspace->nlat, i, j);

          double nrm = std::sqrt(sqr(X[0]) + sqr(X[1]));

          psi[idx] = nrm < cut*d ? fct(X) : rho0;
        }
      }
    }
    else if (initial == 2) { // triangular phase at poles with twist
      double q = 1.0;
      auto fct = [A,rho0,q](std::array<Real, 3> const& x) {
        return rho0 + A*(2.0*std::cos(q*x[0]) - 4.0*std::cos(q*(std::sqrt(3.0)/2.0 * x[1] - M_PI))*std::cos(q/2.0 * x[0]));
      };

      double twist = pt.get("pfc.twist", 0.0);

      std::array<Real, 3> X, X2;
      for (size_t i = 0; i < workspace->nphi; ++i) {
        for (size_t j = 0; j < workspace->nlat; ++j) {
          get_coord(workspace, i, j, radius, X);
          size_t idx = get_idx(workspace->nlat, i, j);

          double angle = twist * (X[2] > 0.0 ? 1.0 : -1.0);

          X2[0] = std::cos(angle) * X[0] - std::sin(angle) * X[1];
          X2[1] = std::sin(angle) * X[0] + std::cos(angle) * X[1];
          X2[2] = X[2];

          double nrm = std::sqrt(sqr(X2[0]) + sqr(X2[1]));

          psi[idx] = nrm < cut*d ? fct(X2) : rho0;
        }
      }
    }
    else if (initial == 3) { // triangular phase at one pole
      double q = 1.0;
      auto fct = [A,rho0,q](std::array<Real, 3> const& x) {
        return rho0 + A*(2.0*std::cos(q*x[0]) - 4.0*std::cos(q*(std::sqrt(3.0)/2.0 * x[1] - M_PI))*std::cos(q/2.0 * x[0]));
      };

      std::array<Real, 3> X;
      for (size_t i = 0; i < workspace->nphi; ++i) {
        for (size_t j = 0; j < workspace->nlat; ++j) {
          get_coord(workspace, i, j, radius, X);
          size_t idx = get_idx(workspace->nlat, i, j);

          double nrm = std::sqrt(sqr(X[0]) + sqr(X[1]));

          psi[idx] = X[2] > 0.0 && nrm < cut*d ? fct(X) : rho0;
        }
      }
    }

    solution[0]->space_filled(true);


    for (size_t i = 0; i < super::size; ++i) {
      (*psi3)[i] = Px[i] = Py[i] = 0.0;
#ifdef USE_C4
      P3->t(i) = P3->p(i) = 0.0;
#endif
    }
    solution_vec[0]->space_filled(true);

    transform(0);
    back_transform(0);
    transform_vec(0);
  }

  /// copy old solutions
  virtual void initTimestep() override
  {
    super::initTimestep();

    std::memcpy(psi3_old, psi3->dual(), super::size_dual*sizeof(Complex));
#ifdef USE_C4
    std::memcpy(P3_old_x, P3->dual_t(), super::size_dual*sizeof(Complex));
    std::memcpy(P3_old_y, P3->dual_p(), super::size_dual*sizeof(Complex));
#endif
  }

  /// calculate nonlinear terms
  virtual void beginIteration(size_t iter) override
  {
    Real* psi = solution[0]->values();
    Real* Px = solution_vec[0]->values_t();
    Real* Py = solution_vec[0]->values_p();

    for (size_t i = 0; i < super::size; ++i)
      (*psi3)[i] = pow3(psi[i]);
    psi3->space_filled(true);

#ifdef USE_C4
    std::array<Real, 3> X;
    for (size_t i = 0; i < workspace->nphi; ++i) {
      for (size_t j = 0; j < workspace->nlat; ++j) {
        sphere_to_cartesian(workspace, Px, Py, i, j, X);

        Real nrm2 = sqr(X[0]) + sqr(X[1]) + sqr(X[2]);
        X[0] *= nrm2; X[1] *= nrm2; X[2] *= nrm2;
        cartesian_to_sphere(workspace, X, P3->values_t(), P3->values_p(), i, j);
      }
    }
    P3->space_filled(true);
#endif

    transform(*psi3);
#ifdef USE_C4
      transform_vec(*P3);
#endif
  }

  /// use a simple Euler discretization at the first timestep
  virtual void firstIteration(size_t iter) override
  {
    EigenVector rhs(3), erg(3);
    size_t it = iteration;
    double tau = timestep;

    for (size_t lm = 0; lm < super::size_dual; ++lm) {
      long int k = workspace->li[lm];

      rhs(0) = getOldSolution(0,it)[lm] + tau*laplace[k] * psi3->dual(lm);
      rhs(1) = getOldSolutionVec<0>(0,it)[lm];
      rhs(2) = getOldSolutionVec<1>(0,it)[lm];

#ifdef USE_C4
      rhs(1)+= C4*tau*(laplace[k] - Dr) * P3->dual_t(lm);
      rhs(2)+= C4*tau*(laplace[k] - Dr) * P3->dual_p(lm);
#endif

      erg = mat[k].solve(rhs);

      solution[0]->dual(lm) = erg(0);
      solution_vec[0]->dual_t(lm) = erg(1);
      solution_vec[0]->dual_p(lm) = erg(2);
    }
  }

  /// use a splitting method for the time-discretization
  virtual void oneIteration(size_t iter) override
  {
    EigenVector rhs(3), erg(3);
    size_t it = iteration;
    double tau = timestep;

    for (size_t lm = 0; lm < super::size_dual; ++lm) {
      long int k = workspace->li[lm];

      rhs(0) = 4.0*getOldSolution(0,it)[lm] - getOldSolution(0,it-1)[lm]
                    + 2.0*tau*laplace[k]*(2.0*psi3->dual(lm) - psi3_old[lm]);
      rhs(1) = 4.0*getOldSolutionVec<0>(0,it)[lm] - getOldSolutionVec<0>(0,it-1)[lm];
      rhs(2) = 4.0*getOldSolutionVec<1>(0,it)[lm] - getOldSolutionVec<1>(0,it-1)[lm];

#ifdef USE_C4
      rhs(1)+= 2.0*C4*tau*(laplace[k] - Dr)*(2.0*P3->dual_t(lm) - P3_old_x[lm]);
      rhs(2)+= 2.0*C4*tau*(laplace[k] - Dr)*(2.0*P3->dual_p(lm) - P3_old_y[lm]);
#endif

#ifdef USE_EXPLICIT_TRANSPORT
      rhs(0)-= 2.0*tau*v0*div[k]*getOldSolutionVec<0>(0,it)[lm];
      rhs(1)-= 2.0*tau*v0*grad[k]*getOldSolution(0,it)[lm];
      erg = mat2[k].solve(rhs);
#else
      erg = mat3[k].solve(rhs);
#endif

      solution[0]->dual(lm) = erg(0);
      solution_vec[0]->dual_t(lm) = erg(1);
      solution_vec[0]->dual_p(lm) = erg(2);
    }
  }

  /// apply a filter to the solution and transform some of the solutions back to real space
  virtual void endIteration(size_t iter) override
  {
    back_transform(0);
  }

  virtual void closeTimestep() override
  {
    back_transform_vec(0);
    auto energies = calc_energy(workspace, *solution[0], *solution_vec[0], weights->values(), n_weights);

    if (observer3(time)) {
      char buf[31];
      int my_val = write_iteration++;
      sprintf(buf, "%05d", my_val);

      std::string out_dir = "polar_pfc/output";
      std::string filename = pt.get("output.dir", out_dir) + "/particles_" + std::string(buf) + ".csv";
      particle_writer->write(filename, getSolutionVec(), getSolution(), observer3.getTimestep(), v0);

      std::string filename2 = pt.get("output.dir", out_dir) + "/order_parameters.csv";
      write_order_parameters(filename2, energies);

      ++observer3;
    }

    super::closeTimestep();
  }

  std::string get_read_filename() const { return read_fn; }
  std::string get_write_filename() const { return write_fn; }


  void write_order_parameters(std::string filename, std::array<Real,2> const& energy_)
  {
    auto P = global_polarization(workspace, *solution_vec[0], weights->values(), n_weights);
//     auto w = angular_momentum(workspace, *solution_vec[0], weights->values(), n_weights);

    auto P_ = global_polarization_vec(workspace, *solution_vec[0], weights->values(), n_weights);
    auto w_ = angular_momentum_vec(workspace, *solution_vec[0], weights->values(), n_weights);
    auto Ppsi_ = global_polarization2_vec(workspace, *solution[0], *solution_vec[0], weights->values(), n_weights);

    std::ofstream out(filename, std::ios_base::app);
    out << time;
    out << "," << P/(4*M_PI*sqr(radius));
    out << "," << P/(4*M_PI*sqr(radius));

    for (auto p_i : P_)
      out << "," << p_i;

    for (auto w_i : w_)
      out << "," << w_i;

    for (auto p_i : Ppsi_)
      out << "," << p_i;

    for (auto e : energy_)
      out << "," << e;

    out << "\n";
  }



  template <class Real, class Complex>
  Real global_polarization(shtns_cfg& shtns, VectorValues<Real,Complex> const& P, Real* weights, size_t& n_weights)
  {
    if (weights[0] == 0.0 || n_weights == 0) {
      n_weights = shtns_gauss_wts(shtns, weights);
    }

    Real dphi = 2.0*M_PI / shtns->nphi;
    std::array<Real, 3> p;

    Real result = 0.0;
    for (size_t i = 0; i < shtns->nlat; ++i) {
      Real sum = 0.0;
      for (size_t j = 0; j < shtns->nphi; ++j) {
        sphere_to_cartesian(shtns, P.values_t(), P.values_p(), j, i, p);
        sum += norm(p);
      }

      sum *= weights[sym_idx(i, n_weights)]*dphi;
      result += sum;
    }

    return sqr(radius)*result;
  }



  template <class Real, class Complex>
  std::array<Real, 3> global_polarization_vec(shtns_cfg& shtns, VectorValues<Real,Complex> const& P, Real* weights, size_t& n_weights)
  {
    if (weights[0] == 0.0 || n_weights == 0) {
      n_weights = shtns_gauss_wts(shtns, weights);
    }

    Real dphi = 2.0*M_PI / shtns->nphi;
    std::array<Real, 3> p;

    std::array<Real, 3> result = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < shtns->nlat; ++i) {
      std::array<Real, 3> sum = {0.0, 0.0, 0.0};
      for (size_t j = 0; j < shtns->nphi; ++j) {
        sphere_to_cartesian(shtns, P.values_t(), P.values_p(), j, i, p);
        sum += p;
      }

      sum *= weights[sym_idx(i, n_weights)]*dphi;
      result += sum;
    }

    return sqr(radius)*result;
  }



  template <class Real, class Complex>
  std::array<Real, 3> global_polarization2_vec(shtns_cfg& shtns, ScalarValues<Real,Complex> const& psi, VectorValues<Real,Complex> const& P, Real* weights, size_t& n_weights)
  {
    if (weights[0] == 0.0 || n_weights == 0) {
      n_weights = shtns_gauss_wts(shtns, weights);
    }

    Real dphi = 2.0*M_PI / shtns->nphi;
    std::array<Real, 3> p;


    double min_psi =  1.0e50;
    for (size_t i = 0; i < shtns->nphi; ++i) {
      for (size_t j = 0; j < shtns->nlat; ++j) {
        size_t idx = get_idx(shtns->nlat, i, j);
        min_psi = std::min( min_psi, psi[idx] );
      }
    }


    std::array<Real, 3> result = {0.0, 0.0, 0.0};
    Real int_psi = 0.0;
    for (size_t i = 0; i < shtns->nlat; ++i) {
      std::array<Real, 3> sum = {0.0, 0.0, 0.0};
      Real sum2 = 0.0;
      for (size_t j = 0; j < shtns->nphi; ++j) {
        size_t idx = sphere_to_cartesian(shtns, P.values_t(), P.values_p(), j, i, p);
        sum += (psi[idx] - min_psi) * p;
        sum2 += (psi[idx] - min_psi);
      }

      sum *= weights[sym_idx(i, n_weights)]*dphi;
      sum2 *= weights[sym_idx(i, n_weights)]*dphi;
      result += sum;
      int_psi += sum2;
    }

    return result * (1.0/int_psi);
  }


//   template <class Real, class Complex>
//   Real angular_momentum(shtns_cfg& shtns, VectorValues<Real,Complex> const& P, Real* weights, size_t& n_weights)
//   {
//     if (weights[0] == 0.0 || n_weights == 0) {
//       n_weights = shtns_gauss_wts(shtns, weights);
//     }
//
//     Real dphi = 2.0*M_PI / shtns->nphi;
//     std::array<Real, 3> p, x;
//
//     Real result = 0.0;
//     for (size_t i = 0; i < shtns->nlat; ++i) {
//       Real sum = 0.0;
//       for (size_t j = 0; j < shtns->nphi; ++j) {
//         get_coord(shtns, j, i, radius, x);
//         auto normal = x*(1.0/norm(x));
//         sphere_to_cartesian(shtns, P.values_t(), P.values_p(), j, i, p);
//         sum += norm(cross(normal, p));
//       }
//
//       sum *= weights[sym_idx(i, n_weights)]*dphi;
//       result += sum;
//     }
//
//     return sqr(radius)*result;
//   }


  template <class Real, class Complex>
  std::array<Real, 3> angular_momentum_vec(shtns_cfg& shtns, VectorValues<Real,Complex> const& P, Real* weights, size_t& n_weights)
  {
    if (weights[0] == 0.0 || n_weights == 0) {
      n_weights = shtns_gauss_wts(shtns, weights);
    }

    Real dphi = 2.0*M_PI / shtns->nphi;
    std::array<Real, 3> p, x;

    std::array<Real, 3> result = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < shtns->nlat; ++i) {
      std::array<Real, 3> sum = {0.0, 0.0, 0.0};
      for (size_t j = 0; j < shtns->nphi; ++j) {
        get_coord(shtns, j, i, radius, x);
        x *= (1.0/norm(x));
        sphere_to_cartesian(shtns, P.values_t(), P.values_p(), j, i, p);
        sum += cross(x, p);
      }

      sum *= weights[sym_idx(i, n_weights)]*dphi;
      result += sum;
    }

    return sqr(radius)*result;
  }

  // compute the positions of all orientational defects and return it as vector of
  // spherical coordinates
  template <class Real, class Complex>
  std::vector<std::array<Real,3>>
  defect_positions (shtns_cfg& shtns,
                    ScalarValues<Real,Complex> const& psi,
                    VectorValues<Real,Complex> const& P,
                    Real* weights, size_t& n_weights)
  {
    ScalarValues<Real,Complex> Q0{"Q0",shtns}, Q1{"Q1",shtns}, Q2{"Q2",shtns},
                               Q3{"Q3",shtns}, Q4{"Q4",shtns};

    // 1. compute the Qi q-tensor coefficient functions

    for (size_t i = 0; i < shtns->nphi; ++i) {
      for (size_t j = 0; j < shtns->nlat; ++j) {
        // (i,j) is the index pair of a spherical point

        // x={sin(theta_j) * cos(phi_i), sin(theta_j) * cos(phi_i), cos(theta_j)}
        std::array<Real,3> x;
        get_coord(shtns, i, j, radius, x);

        // compute the normal vector
        std::array<Real,3> n = x;
        n *= 1.0/norm(n);

        // p = P(x) is the evaluation of the vector-valued function P at x
        std::array<Real,3> p;
        size_t idx = sphere_to_cartesian(shtns, P.values_t(), P.values_p(), i, j, p);
        // idx = i*nlat + j is the global index to address scalar or vector values on the sphere

        // normalize the polar vector
        p *= 1.0/norm(p);

        // the director
        std::array<Real,3> d = cross(n, p);

        // symmetric part of q-tensor
        std::array<Real,5> q = {
          d[0]*d[0] - 0.5*(1.0 - n[0]*n[0]),
          d[0]*d[1] - 0.5*(0.0 - n[0]*n[1]),
          d[0]*d[2] - 0.5*(0.0 - n[0]*n[2]),
          d[1]*d[1] - 0.5*(1.0 - n[1]*n[1]),
          d[1]*d[2] - 0.5*(0.0 - n[1]*n[2])
        };

        // TODO: smoothing of the qi

        // compute the norm of the tensor from the symmetric components only
        double nrm_q = std::sqrt(2.0 * (q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4] + q[0]*q[3]));

        // normalize q
        q *= 1.0/nrm_q;

        // Assign computed q-tensor components to scalar functions Qi
        Q0[idx] = q[0];
        Q1[idx] = q[1];
        Q2[idx] = q[2];
        Q3[idx] = q[3];
        Q4[idx] = q[4];
      }
    }
    // explicitly state that we have filled in the values
    Q0.space_filled(true);
    Q1.space_filled(true);
    Q2.space_filled(true);
    Q3.space_filled(true);
    Q4.space_filled(true);

    // make the spherical fourier transform of the scalar fields Qi
    transform(Q0);
    transform(Q1);
    transform(Q2);
    transform(Q3);
    transform(Q4);

    // 2. compute the gradients of the scalar fields: grad(Qi)

    VectorValues<Real,Complex> GradQ0{"GradQ0",shtns}, GradQ1{"GradQ1",shtns},
                               GradQ2{"GradQ2",shtns}, GradQ3{"GradQ3",shtns},
                               GradQ4{"GradQ4",shtns};

    for (size_t lm = 0; lm < super::size_dual; ++lm) {
      // lm=(l,m) is a global contiguous index of fourier coefficients

      // Assign coefficient to first component of fourier coefficients
      GradQ0.dual_t(lm) = Q0.dual(lm)/radius;
      GradQ1.dual_t(lm) = Q1.dual(lm)/radius;
      GradQ2.dual_t(lm) = Q2.dual(lm)/radius;
      GradQ3.dual_t(lm) = Q3.dual(lm)/radius;
      GradQ4.dual_t(lm) = Q4.dual(lm)/radius;
    }

    // explicitly state that we have filled in the fourier coefficients
    GradQ0.dual_filled(true);
    GradQ1.dual_filled(true);
    GradQ2.dual_filled(true);
    GradQ3.dual_filled(true);
    GradQ4.dual_filled(true);

    // make the spherical fourier transform of the gradient fields GradQi back to real space
    back_transform(GradQ0);
    back_transform(GradQ1);
    back_transform(GradQ2);
    back_transform(GradQ3);
    back_transform(GradQ4);

    // 3. detect the defects

    ScalarValues<Real,Complex> FD{"FD",shtns};
    for (size_t i = 0; i < shtns->nphi; ++i) {
      for (size_t j = 0; j < shtns->nlat; ++j) {

        // evaluate the gradient functions in spherical coordinates
        std::array<Real,3> grad_q0, grad_q1, grad_q2, grad_q3, grad_q4;
        sphere_to_cartesian(shtns, GradQ0.values_t(), GradQ0.values_p(), i, j, grad_q0);
        sphere_to_cartesian(shtns, GradQ1.values_t(), GradQ1.values_p(), i, j, grad_q1);
        sphere_to_cartesian(shtns, GradQ2.values_t(), GradQ2.values_p(), i, j, grad_q2);
        sphere_to_cartesian(shtns, GradQ3.values_t(), GradQ3.values_p(), i, j, grad_q3);
        sphere_to_cartesian(shtns, GradQ4.values_t(), GradQ4.values_p(), i, j, grad_q4);

        // compute the distorsion indicator
        double fD = norm2(grad_q0) + norm2(grad_q1) + norm2(grad_q2)
                  + norm2(grad_q3) + norm2(grad_q4);

        // assign it to a global vector
        size_t idx = get_idx(shtns->nlat,i,j);
        FD[idx] = fD;
      }
    }

    std::vector<std::array<Real,3>> positions;
    for (size_t i = 0; i < shtns->nphi; ++i) {
      for (size_t j = 0; j < shtns->nlat; ++j) {

        std::array<Real,3> x;
        get_coord(shtns, i, j, radius, x);

        size_t idx = get_idx(shtns->nlat,i,j);

        // compute global index of right, left, top, bottom vertex
        size_t idx_r = get_idx(shtns->nlat, (i + shtns->nphi + 1)%(shtns->nphi), j);
        size_t idx_l = get_idx(shtns->nlat, (i + shtns->nphi - 1)%(shtns->nphi), j);
        size_t idx_t = get_idx(shtns->nlat, i, std::max(0,std::min(int(shtns->nlat),int(j)+1)));
        size_t idx_b = get_idx(shtns->nlat, i, std::max(0,std::min(int(shtns->nlat),int(j)-1)));

        if (FD[idx]<FD[idx_l] && FD[idx]<FD[idx_r] && FD[idx]<FD[idx_b] && FD[idx]<FD[idx_t])
          positions.push_back(x);
      }
    }

    return positions;
  }

  template <class Real, class Complex>
  std::array<Real, 2> calc_energy(shtns_cfg& shtns, ScalarValues<Real,Complex> const& psi, VectorValues<Real,Complex> const& P, Real* weights, size_t& n_weights)
  {
    for (size_t lm = 0; lm < super::size_dual; ++lm) {
      long int k = workspace->li[lm];

      energy->dual(lm) = (eps + sqr(1.0 + laplace[k])) * psi.dual(lm);
    }
    energy->dual_filled(true);
    back_transform(*energy);


    if (weights[0] == 0.0 || n_weights == 0) {
      n_weights = shtns_gauss_wts(shtns, weights);
    }

    Real dphi = 2.0*M_PI / shtns->nphi;
    std::array<Real, 3> p, x;

    std::array<Real, 2> result = {0.0, 0.0};
    for (size_t i = 0; i < shtns->nlat; ++i) {
      std::array<Real, 2> sum = {0.0, 0.0};
      for (size_t j = 0; j < shtns->nphi; ++j) {
        size_t idx = sphere_to_cartesian(shtns, P.values_t(), P.values_p(), j, i, p);

        (*energyPsi)[idx] = 0.5*psi[idx]*(*energy)[idx] + 0.25*sqr(sqr(psi[idx]));
        sum[0] += (*energyPsi)[idx];

        auto nrm2 = dot(p,p);
        (*energyP)[idx] = 0.5*C1*nrm2 + 0.25*C4*sqr(nrm2);
        sum[1] += (*energyP)[idx];
      }

      sum *= weights[sym_idx(i, n_weights)]*dphi;
      result += sum;
    }
    energyPsi->space_filled(true);
    energyP->space_filled(true);

    return sqr(radius)*result;
  }


  void setEnergyPsi(std::shared_ptr<ScalarValues<Real, Complex>> vec)
  {
    energyPsi = std::move(vec);
  }

  void setEnergyP(std::shared_ptr<ScalarValues<Real, Complex>> vec)
  {
    energyP = std::move(vec);
  }

  std::shared_ptr<ScalarValues<Real, Complex>> const& getEnergyPsi() const
  {
    return energyPsi;
  }

  std::shared_ptr<ScalarValues<Real, Complex>> const& getEnergyP() const
  {
    return energyP;
  }

private:
  double eps;
  double density;
  double C1, C4;
  double v0;
  double Dr;

  int write_iteration = 0;

  bool read_from_file;
  std::string read_fn, write_fn;

  std::vector<Eigen::FullPivLU<EigenMatrixReal> > mat;
#ifdef USE_EXPLICIT_TRANSPORT
  std::vector<Eigen::FullPivLU<EigenMatrixReal> > mat2;
#else
  std::vector<Eigen::FullPivLU<EigenMatrixReal> > mat3;
#endif

  ScalarValues<Real, Complex>* psi3;
  Complex *psi3_old;
  VectorValues<Real, Complex>* P3;
  Complex *P3_old_x, *P3_old_y;

  std::shared_ptr<ScalarValues<Real, Complex>> weights;
  std::shared_ptr<ScalarValues<Real, Complex>> energy;
  std::shared_ptr<ScalarValues<Real, Complex>> energyPsi, energyP;
  size_t n_weights = 0;

  std::unique_ptr<ParticleWriter<Real, Complex>> particle_writer;
  Observer observer3;
};
