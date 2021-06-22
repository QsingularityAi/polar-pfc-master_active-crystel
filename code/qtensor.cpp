#include <cmath>
#include <iostream>
#include <complex>
#include <memory>

#include "DOFVectorBase.hpp"
#include "ProblemSHTNS.hpp"
#include "initial.hpp"
#include "utility/Filesystem.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

typedef double Real;
typedef std::complex<double> Complex;

using EigenMatrix = Eigen::Matrix<Real, 4, 4>;
using EigenVector = Eigen::Matrix<Complex, 4, 1>;
using WorldVector = std::array<Real, 3>;


using namespace Shtns;

struct QTensorPotential : ProblemSHTNS<Real, Complex>
{
  typedef ProblemSHTNS<Real, Complex> super;

  QTensorPotential(std::string name = "qtensor")
      : super(name, 1)
  {
    L = pt.get("Q.L", 1.0);
    a = pt.get("Q.a", 2.0/3.0);
    c = pt.get("Q.c", 1.0/3.0);
  }


  /// allocate data for nonlinear terms and corresponding old solutions
  virtual void initData(Flag init_flag = INIT_ALL) override
  {
    super::initData(init_flag);

    f1.reset(new ScalarValues<Real, Complex>("f1", workspace));
    f2.reset(new ScalarValues<Real, Complex>("f2", workspace));
    f3.reset(new ScalarValues<Real, Complex>("f3", workspace));
//     f4.reset(new ScalarValues<Real, Complex>("f4", workspace));
    f5.reset(new ScalarValues<Real, Complex>("f5", workspace));
    f6.reset(new ScalarValues<Real, Complex>("f6", workspace));
    f7.reset(new ScalarValues<Real, Complex>("f7", workspace));
    f8.reset(new ScalarValues<Real, Complex>("f8", workspace));

//     v1.reset(new VectorValues<Real, Complex>("v1", workspace));
    v2.reset(new VectorValues<Real, Complex>("v2", workspace));
    v3.reset(new VectorValues<Real, Complex>("v3", workspace));
    v4.reset(new VectorValues<Real, Complex>("v4", workspace));

    createMarix();
  }

  /// precalculate coefficients
  void createMarix()
  {
    mat.resize(workspace->lmax + 1);
    laplaceQ.resize(workspace->lmax + 1);

    double tau = timestep;
    for (long int l = 0; l <= workspace->lmax; ++l)
    {
      laplaceQ[l] = (l*(l+1.0)/sqr(radius)) + 2.0/sqr(radius);
      double alpha = -laplaceQ[l];
      double beta   = 1.0 - tau*L*(laplaceQ[l]) - tau*a;

      alpha = (alpha == 0 ? 1.0 : alpha);
      beta = (beta == 0 ? 1.0 : beta);


      laplaceQ[l] = (laplaceQ[l] == 0 ? 1.0 : laplaceQ[l]);

      EigenMatrix B;
      B(0,0) = alpha; B(0,1) = 0.0;   B(0,2) = 1.0;  B(0,3) = 0.0;
      B(1,0) = 0.0;   B(1,1) = alpha; B(1,2) = 0.0;  B(1,3) = 1.0;
      B(2,0) = 0.0;   B(2,1) = 0.0;   B(2,2) = beta; B(2,3) = 0.0;
      B(3,0) = 0.0;   B(3,1) = 0.0;   B(3,2) = 0.0;  B(3,3) = beta;
      mat[l].compute(B);

      if (!mat[l].isInvertible())
        throw std::runtime_error("Matrix not invertible!");
    }
  }

  /// set the initial solution for psi and P,
  /// and set the nonlinear terms to zero
  virtual void solveInitialProblem() override
  {
    Real* Qt = solution_vec[1]->values_t();
    Real* Qp = solution_vec[1]->values_p();

    // fill u0 and f
    InitialRandom fct;
    for (size_t i = 0; i < super::size; ++i) {
      Qt[i] = fct(M_PI/10.0, 0.0);
      Qp[i] = fct(M_PI/10.0, 0.0);
    }
    solution_vec[1]->space_filled(true);
    transform_vec(1);
    back_transform_vec(1); // -> q

    transform_vec(1);
    for (size_t lm = 0; lm < super::size_dual; ++lm) {
      long int l = workspace->li[lm];
      solution_vec[0]->dual_t(lm) = solution_vec[1]->dual_t(lm)/laplaceQ[l];
      solution_vec[0]->dual_p(lm) = solution_vec[1]->dual_p(lm)/laplaceQ[l];

      f2->dual(lm) = div[l] * solution_vec[0]->dual_t(lm);
      f3->dual(lm) = rot[l] * solution_vec[0]->dual_p(lm);
    }
    solution_vec[0]->dual_filled(true);
    back_transform_vec(0); // -> p = [laplaceDR + 2/R^2]^(-1) * q

    f2->dual_filled(true);
    f3->dual_filled(true);
    back_transform(*f2); // div(p)
    back_transform(*f3); // rot(p)

#if 0
    Real nrm_f2 = 0.0, nrm_f3 = 0.0, nrm_q = 0.0;
    for (size_t i = 0; i < workspace->nphi; ++i) {
      for (size_t j = 0; j < workspace->nlat; ++j) {
        const size_t idx = get_idx(workspace->nlat, i, j);

        nrm_f2 += sqr((*f2)[idx]);
        nrm_f3 += sqr((*f3)[idx]);

        nrm_q += sqr(solution_vec[1]->t(idx)) + sqr(solution_vec[1]->p(idx));
      }
    }
    nrm_f2 /= super::size;
    nrm_f3 /= super::size;
    nrm_q /= super::size;
    std::cout << "nrm_f2 = " << nrm_f2 << "\n";
    std::cout << "nrm_f3 = " << nrm_f3 << "\n";
    std::cout << "nrm_q = " << nrm_q << "\n";
#endif
  }

  // <v,w>
  template <class T>
  T dot(T vt, T vp, T wt, T wp, size_t j) const
  {
    return sqr(radius)*(vt*wt + sqr(workspace->st[j])*vp*wp);
  }

  // <*v, w>
  template <class T>
  T dot_orth1(T vt, T vp, T wt, T wp, size_t j) const
  {
    return sqr(radius)*workspace->st[j]*(-vp*wt + vt*wp);
  }

  template <class T>
  T norm2(T vt, T vp, size_t j) const
  {
    return sqr(radius*vt) + sqr(radius*workspace->st[j]*vp);
  }

  /// calculate nonlinear terms
  virtual void beginIteration(size_t iter) override
  {
    Real* Pt = solution_vec[0]->values_t();
    Real* Pp = solution_vec[0]->values_p();

    Real* Qt = solution_vec[1]->values_t();
    Real* Qp = solution_vec[1]->values_p();

    Real kappa = 1.0/sqr(radius);

    WorldVector P, Q;
    Real nrm_f1 = 0.0;
    for (size_t i = 0; i < workspace->nphi; ++i) {
      for (size_t j = 0; j < workspace->nlat; ++j) {
        const size_t idx = get_idx(workspace->nlat, i, j);

        (*f1)[idx] = norm2(Pt[idx], Pp[idx], j);
        (*f5)[idx] = 2.0*(2.0*kappa*(*f1)[idx] - 2.0*dot(Qt[idx], Qp[idx], Pt[idx], Pp[idx], j) - sqr((*f2)[idx]) - sqr((*f3)[idx]));

        nrm_f1 += sqr((*f1)[idx]);
      }
    }
    nrm_f1 /= super::size;
    std::cout << "nrm_f1 = " << nrm_f1 << "\n";

    f1->space_filled(true);
    f5->space_filled(true);
    transform(*f1);  // |p|^2
    transform(*f5);

    for (size_t lm = 0; lm < super::size_dual; ++lm) {
      long int l = workspace->li[lm];

      auto f4 = laplace[l] * f1->dual(lm);
      f5->dual(lm) += 2.0*f4;
      f6->dual(lm) = laplace[l] * f5->dual(lm);

      v2->dual_t(lm) = grad[l] * f5->dual(lm);
      v2->dual_p(lm) = 0.0;
    }
    f5->dual_filled(true);
    f6->dual_filled(true);
    v2->dual_filled(true);
    back_transform(*f5); // |DQ(p)|^2
    back_transform(*f6); // laplace(|DQ(p)|^2)
    back_transform(*v2); // grad(|DQ(p)|^2)


    WorldVector V2, Porth;
    Real nrm_f7 = 0.0, nrm_f8 = 0.0, nrm_f5 = 0.0, nrm_v3 = 0.0, nrm_v4 = 0.0;
    for (size_t i = 0; i < workspace->nphi; ++i) {
      for (size_t j = 0; j < workspace->nlat; ++j) {
        const size_t idx = get_idx(workspace->nlat, i, j);

        v3->t(idx) = (*f5)[idx] * Qt[idx];
        v3->p(idx) = (*f5)[idx] * Qp[idx];

        v4->t(idx) = (*f6)[idx] * Pt[idx];
        v4->p(idx) = (*f6)[idx] * Pp[idx];

        (*f7)[idx] = dot(Pt[idx], Pp[idx], v2->t(idx), v2->p(idx), j);
        (*f8)[idx] = dot_orth1(Pt[idx], Pp[idx], v2->t(idx), v2->p(idx), j);

        nrm_f5 += sqr((*f5)[idx]);
        nrm_f7 += sqr((*f7)[idx]);
        nrm_f8 += sqr((*f8)[idx]);

        nrm_v3 += sqr(v3->t(idx)) + sqr(v3->p(idx));
        nrm_v4 += sqr(v4->t(idx)) + sqr(v4->p(idx));
      }
    }
    v3->space_filled(true);
    v4->space_filled(true);
    f7->space_filled(true);
    f8->space_filled(true);

    transform(*v3); // |DQ(p)|^2*q
    transform(*v4); // laplace(|DQ(p)|^2)*p
    transform(*f7); // <p, grad(|DQ(p)|^2)>
    transform(*f8); // <*p, grad(|DQ(p)|^2)>

    nrm_f5 /= super::size;
    nrm_f7 /= super::size;
    nrm_f8 /= super::size;
    nrm_v3 /= super::size;
    nrm_v4 /= super::size;
    std::cout << "nrm_f5 = " << nrm_f5 << "\n";
    std::cout << "nrm_f7 = " << nrm_f7 << "\n";
    std::cout << "nrm_f8 = " << nrm_f8 << "\n";
    std::cout << "nrm_v3 = " << nrm_v3 << "\n";
    std::cout << "nrm_v4 = " << nrm_v4 << "\n";
  }

  /// use a splitting method for the time-discretization
  virtual void oneIteration(size_t iter) override
  {
    EigenVector rhs(4), erg(4);
    size_t it = iteration;
    double tau = timestep;

    rhs(0) = 0;
    rhs(1) = 0;

    for (size_t lm = 0; lm < super::size_dual; ++lm) {
      long int l = workspace->li[lm];

      rhs(2) = getOldSolutionVec<0>(1,it)[lm] - c * tau * (v3->dual_t(lm) - v4->dual_t(lm) + grad[l]  * f7->dual(lm));
      rhs(3) = getOldSolutionVec<1>(1,it)[lm] - c * tau * (v3->dual_p(lm) - v4->dual_p(lm) - rot_n[l] * f8->dual(lm));

      erg = mat[l].solve(rhs);

      solution_vec[0]->dual_t(lm) = erg(0);
      solution_vec[0]->dual_p(lm) = erg(1);
      solution_vec[1]->dual_t(lm) = erg(2);
      solution_vec[1]->dual_p(lm) = erg(3);
    }
  }

  /// apply a filter to the solution and transform some of the solutions back to real space
  virtual void endIteration(size_t iter) override
  {
    for (size_t lm = 0; lm < super::size_dual; ++lm) {
      long int l = workspace->li[lm];
      f2->dual(lm) = div[l] * solution_vec[0]->dual_t(lm); // div(p)
      f3->dual(lm) = rot[l] * solution_vec[0]->dual_p(lm); // rot(p)
    }
    f2->dual_filled(true);
    f3->dual_filled(true);

    back_transform_vec(0);
    back_transform_vec(1);

    back_transform(*f2);
    back_transform(*f3);
  }

private:
  double L;
  double a;
  double c;

  int write_iteration = 0;

  std::vector<Eigen::FullPivLU<EigenMatrix> > mat;
  std::vector<Real> laplaceQ;

  std::unique_ptr<ScalarValues<Real, Complex>> f1, f2, f3, /*f4,*/ f5, f6, f7, f8;
  std::unique_ptr<VectorValues<Real, Complex>> /*v1,*/ v2, v3, v4;
};



int main(int argc, char** argv)
{
  std::unique_ptr<QTensorPotential> prob;

  if (argc > 1)
    prob.reset(new QTensorPotential(std::string(argv[1])));
  else
    prob.reset(new QTensorPotential);

  // specify the solution components
  VectorValues<Real, Complex> P("P", prob->getWorkspace());
  VectorValues<Real, Complex> Q("Q", prob->getWorkspace());
  prob->add_data(P).add_data(Q);

  DOFVector P_(prob->getMesh(), P);
  DOFVector Q_(prob->getMesh(), Q);

  // initialize the data
  prob->initData(INIT_LAPLACE | INIT_DIVERGENCE | INIT_GRADIENT | INIT_CURL | INIT_CURLN);

  // configure the file-writer
  if (prob->get_pt().get("output.write_vtu", true) ) {
    prob->getFilewriter().attach(P_).attach(Q_);
  }

  // solve the problem
  prob->solve();

  return 0;
}

