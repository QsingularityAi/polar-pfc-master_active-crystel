#pragma once

#include <cassert>

#include "Output.hpp"
#include "Shtns.hpp"

namespace Shtns
{
  ///----------------------------------------------------------
  struct Values
  {
    Values(std::string name, shtns_cfg shtns)
      : name_(name)
      , shtns_(shtns)
    {}

    virtual ~Values() {};

    std::string name() const { return name_; }

    shtns_cfg shtns() { return shtns_; }
    shtns_cfg shtns() const { return shtns_; }

    void space_filled(bool value) { space_filled_ = value; }
    bool space_filled() { return space_filled_; }

    void dual_filled(bool value) { dual_filled_ = value; }
    bool dual_filled() { return dual_filled_; }

    virtual int num_components() const = 0;

  protected:
    shtns_cfg shtns_;

  private:
    std::string name_;

    bool space_filled_ = false;
    bool dual_filled_ = false;
  };


  ///----------------------------------------------------------
  template <class Real, class Complex = std::complex<Real> >
  struct ScalarValues
      : public Values
  {
    ScalarValues(std::string name, shtns_cfg shtns, bool allocate_dual = true, bool create_plan_ = false, bool create_plan_dual_ = false)
      : Values(name, shtns)
      , allocate_dual_(allocate_dual)
      , create_plan(create_plan_)
      , create_plan_dual(create_plan_dual_)
    {
      values_ = (Real*) fftw_malloc( NSPAT_ALLOC(shtns_) * sizeof(Real));

      if (allocate_dual_)
        dual_ = (Complex*) fftw_malloc( shtns_->nlm * sizeof(Complex));
    }

    ~ScalarValues()
    {
      if (create_plan)
        fftw_destroy_plan(plan_);
      if (create_plan_dual)
        fftw_destroy_plan(plan_dual_);

      if (values_) {
        fftw_free(values_);
        values_ = NULL;
      }

      if (allocate_dual_ && dual_) {
        fftw_free(dual_);
        dual_ = NULL;
      }
    }

    Real* values() { return values_; }
    Real const* values() const { return values_; }

    Complex* dual() { return dual_; }
    Complex const* dual() const { return dual_; }

    Real& operator[](size_t i) { return values_[i]; }
    Real const& operator[](size_t i) const { return values_[i]; }

    Complex& dual(size_t lm) { return dual_[lm]; }
    Complex const& dual(size_t lm) const { return dual_[lm]; }

    virtual int num_components() const override { return 1; }

    fftw_plan& plan() { return plan_; }
    fftw_plan& plan_dual() { return plan_dual_; }

  private:
    Real* values_ = nullptr;
    Complex* dual_ = nullptr;

    bool allocate_dual_;

    bool create_plan;
    bool create_plan_dual;
    fftw_plan plan_;
    fftw_plan plan_dual_;
  };


  ///----------------------------------------------------------
  template <class Real, class Complex = std::complex<Real> >
  struct VectorValues
      : public Values
  {
    VectorValues(std::string name, shtns_cfg shtns, bool allocate_dual_ = true)
      : Values(name, shtns)
      , allocate_primal(true)
      , allocate_dual(allocate_dual_)
    {
      values_t_ = (Real*) fftw_malloc( NSPAT_ALLOC(shtns_) * sizeof(Real));
      values_p_ = (Real*) fftw_malloc( NSPAT_ALLOC(shtns_) * sizeof(Real));

      if (allocate_dual) {
        dual_t_ = (Complex*) fftw_malloc( shtns_->nlm * sizeof(Complex));
        dual_p_ = (Complex*) fftw_malloc( shtns_->nlm * sizeof(Complex));
      }
    }

    VectorValues(std::string name, ScalarValues<Real, Complex> &theta_values, ScalarValues<Real, Complex> &phi_values)
      : Values(name, theta_values.shtns())
      , values_t_(theta_values.values())
      , values_p_(phi_values.values())
      , dual_t_(theta_values.dual())
      , dual_p_(phi_values.dual())
      , allocate_primal(false)
      , allocate_dual(false)
    {}

    ~VectorValues()
    {
      if (allocate_primal) {
        if (values_t_) {
          fftw_free(values_t_);
          values_t_ = NULL;
        }

        if (values_p_) {
          fftw_free(values_p_);
          values_p_ = NULL;
        }
      }

      if (allocate_dual) {
        if (dual_t_) {
          fftw_free(dual_t_);
          dual_t_ = NULL;
        }
        if (dual_p_) {
          fftw_free(dual_p_);
          dual_p_ = NULL;
        }
      }
    }

    Real* values_t() { return values_t_; }
    Real* values_p() { return values_p_; }

    Real const* values_t() const { return values_t_; }
    Real const* values_p() const { return values_p_; }

    Complex* dual_t() { return dual_t_; }
    Complex* dual_p() { return dual_p_; }

    Complex const* dual_t() const { return dual_t_; }
    Complex const* dual_p() const { return dual_p_; }


    Real const& t(size_t i) const { return values_t_[i]; }
    Real const& p(size_t i) const { return values_p_[i]; }

    Real& t(size_t i) { return values_t_[i]; }
    Real& p(size_t i) { return values_p_[i]; }

    Complex const& dual_t(size_t i) const { return dual_t_[i]; }
    Complex const& dual_p(size_t i) const { return dual_p_[i]; }

    Complex& dual_t(size_t i) { return dual_t_[i]; }
    Complex& dual_p(size_t i) { return dual_p_[i]; }

    virtual int num_components() const override { return 3; }

  private:
    Real* values_t_ = nullptr;
    Real* values_p_ = nullptr;

    Complex* dual_t_ = nullptr;
    Complex* dual_p_ = nullptr;

    bool allocate_primal;
    bool allocate_dual;
  };

  // ______________________________________________________________________________
  // some shortcuts

  template <class R, class C>
  inline void spat_to_SH(ScalarValues<R,C>& data)
  {
    warn_msg( data.space_filled(), "W0: Space-data of ", data.name(), " probably not set!" );
    spat_to_SH(data.shtns(), data.values(), data.dual());
    data.dual_filled(true);
    data.space_filled(false);
  }

  template <class R, class C>
  inline void SH_to_spat(ScalarValues<R,C>& data)
  {
    warn_msg( data.dual_filled(), "W1: Dual-data of ", data.name(), " probably not set!" );
    SH_to_spat(data.shtns(), data.dual(), data.values());
    data.space_filled(true);
  }

  template <class R, class C>
  inline void spat_to_VSH(VectorValues<R,C>& data)
  {
    warn_msg( data.space_filled(), "W2: Space-data of ", data.name(), " probably not set!" );
    spat_to_SHsphtor(data.shtns(), data.values_t(), data.values_p(), data.dual_t(), data.dual_p());
    data.dual_filled(true);
    data.space_filled(false);
  }

  template <class R, class C>
  inline void VSH_to_spat(VectorValues<R,C>& data)
  {
    warn_msg( data.dual_filled(), "W3: Dual-data of ", data.name(), " probably not set!" );
    SHsphtor_to_spat(data.shtns(), data.dual_t(), data.dual_p(), data.values_t(), data.values_p());
    data.space_filled(true);
  }

} // end namespace Shtns
