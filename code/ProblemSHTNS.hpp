#pragma once

#include "tools.hpp"
#include "observer.hpp"
#include "flags.hpp"
#include "coefficients.hpp"
#include "container.hpp"
#include "io/VtkWriter.hpp"
#include "io/VtkTimeseriesWriter.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace Shtns
{

  /// Abstract Problem description, as in AMDiS
  template <class Real, class Complex>
  struct ProblemSHTNS
  {
    typedef ProblemSHTNS self;

    ProblemSHTNS(std::string name_, size_t n_old = 1)
        : name(name_),
          use_old_solution(n_old > 0),
          n_old_solutions(n_old),
          grad(1.0), rot_n(1.0), rot_inv(1.0)
    {
      // init-file is alwas named as follows
      read_json(name + ".json", pt);

      // create workspace
      size_t lmax = pt.get("shtns.lmax", 196), nlat = pt.get("shtns.nlat", 200);
      size_t mmax = pt.get("shtns.mmax", 196), nphi = pt.get("shtns.nphi", 400);
      size_t mres = 1;
      radius = pt.get("shtns.R", 1.0);
      size_t n_threads = pt.get("shtns.threads", 1);

      shtns_verbose(1);
      shtns_use_threads(n_threads);
      workspace = shtns_init( sht_gauss, lmax, mmax, mres, nlat, nphi );
      size = workspace->nphi * workspace->nlat;
      size_dual = workspace->nlm;

      // init time parameters
      timestep =   pt.get("time.dt", 1.e-3);
      start_time = pt.get("time.t0", 0.0);
      end_time =   pt.get("time.tend", 100.0);
      n_iter =     pt.get("time.n_iter", 1); // number of fixed point iteration loops

      tol = 1.e-10;
      time = start_time;

      // init output parameters
      double dt_write = pt.get("output.dt", 0.1);
      observer.setTimestep(dt_write);
      observer.setEndTime(end_time);
      observer.setStartTime(start_time);

      // write statistics only after timestep dt2
      double dt_write2 = pt.get("output.dt2", dt_write);
      observer2.setTimestep(dt_write2);
      observer2.setEndTime(end_time);
      observer2.setStartTime(start_time);

      std::string dir = pt.get("output.dir", std::string("./data"));
      std::string fname = pt.get("output.filename", name);

      write_vtu = pt.get("output.write_vtu", write_vtu);

      // create a file-writer that must be configured by the user
      using namespace Shtns::filesystem;
      auto output = path(dir) /= path(fname);
      filewriter.reset(new VtkWriter(output.string(), VtkWriter::COMPRESSED, VtkWriter::FLOAT32));
      pvdwriter.reset(new VtkTimeseriesWriter(*filewriter));

      mesh.reset(new Mesh(workspace, radius));
    }

    /// Destructor, frees all old-solutions
    ~ProblemSHTNS()
    {
      if (use_old_solution) {
        for (size_t i = 0; i < old_solutions.size(); ++i) {
          for (size_t j = 0; j < old_solutions[i]->size(); ++j) {
            fftw_free((*old_solutions[i])[j]); (*old_solutions[i])[j] = NULL;
          }
        }
        for (size_t i = 0; i < old_solutions_vec.size(); ++i) {
          for (size_t j = 0; j < old_solutions_vec[i]->size(); ++j) {
            fftw_free((*old_solutions_vec[i])[j].first); (*old_solutions_vec[i])[j].first = NULL;
            fftw_free((*old_solutions_vec[i])[j].second); (*old_solutions_vec[i])[j].second = NULL;
          }
        }
      }
    }


    /// initialize the old_solution
    virtual void initData(Flag init_flag = INIT_ALL)
    {
      use_old_solution = pt.get(name + ".use_old_solution", true);

      size_t n_old = std::max((size_t)1,(size_t)n_old_solutions);
      for (size_t i = 0; i < n_old; ++i) {
        old_solutions.push_back(new std::vector<Complex*>(solution.size()));
        old_solutions_vec.push_back(new std::vector<std::pair<Complex*,Complex*> >(solution_vec.size()));
      }

      if (use_old_solution) {
        for (size_t i = 0; i < n_old; ++i) {
          for (size_t j = 0; j < solution.size(); ++j)
            (*old_solutions[i])[j] = (Complex*) fftw_malloc(size_dual * sizeof(Complex));
          for (size_t j = 0; j < solution_vec.size(); ++j) {
            (*old_solutions_vec[i])[j].first = (Complex*) fftw_malloc(size_dual * sizeof(Complex));
            (*old_solutions_vec[i])[j].second = (Complex*) fftw_malloc(size_dual * sizeof(Complex));
          }
        }
      } else {
        for (size_t i = 0; i < solution.size(); ++i)
          (*old_solutions[0])[i] = solution[i]->dual();
        for (size_t i = 0; i < solution_vec.size(); ++i) {
          (*old_solutions_vec[0])[i].first = solution_vec[i]->dual_t();
          (*old_solutions_vec[0])[i].second = solution_vec[i]->dual_p();
        }
      }

      // provide operators for some derivative terms
      size_t max_idx = workspace->lmax + 1;
      if (init_flag.isSet(INIT_LAPLACE))    laplace.resize(max_idx);
      if (init_flag.isSet(INIT_DIVERGENCE)) div.resize(max_idx);
      if (init_flag.isSet(INIT_CURL))       rot.resize(max_idx);

      if (init_flag.isSet(INIT_GRADIENT))
        grad = ConstantCoeff<Real>(1.0/radius);

      if (init_flag.isSet(INIT_CURLN))
        rot_n = ConstantCoeff<Real>(-1.0/radius);

      if (init_flag.isSet(INIT_CURLINV))
        rot_inv = ConstantCoeff<Real>(radius);

      for (long int l = 0; l < max_idx; ++l)
      {
        if (init_flag.isSet(INIT_LAPLACE))
          laplace[l] = -(l*(l+1.0)/(radius*radius));

        if (init_flag.isSet(INIT_DIVERGENCE))
          div[l] = -(l*(l+1.0)/radius);

        if (init_flag.isSet(INIT_CURL))
          rot[l] = -(l*(l+1.0)/radius);
      }
    }


    /// create an initial solution
    virtual void solveInitialProblem() {}


    /// call \ref transform for all plans and write initial solution
    virtual void transferInitialSolution()
    {
      std::cout << "timestep t = " << time << "\n";
      if (write_vtu) {
        pvdwriter->write(time);
      }

  //     for (size_t i = 0; i < solution.size(); ++i)
  //       transform(i);
  //     for (size_t i = 0; i < solution_vec.size(); ++i)
  //       transform_vec(i);

      first_seq = boost::posix_time::microsec_clock::local_time();
    }


    /// called at the beginning of each timestep
    virtual void initTimestep()
    {
      if (use_old_solution) {
        for (size_t i = 0; i < solution.size(); ++i) {
          std::memcpy(getOldSolution(i, iteration), solution[i]->dual(), size_dual*sizeof(Complex));
        }
        for (size_t i = 0; i < solution_vec.size(); ++i) {
          std::memcpy(getOldSolutionVec<0>(i, iteration), solution_vec[i]->dual_t(), size_dual*sizeof(Complex));
          std::memcpy(getOldSolutionVec<1>(i, iteration), solution_vec[i]->dual_p(), size_dual*sizeof(Complex));
        }
      }
    }


    /// write solution vector to file
    virtual void closeTimestep()
    {
      using namespace boost::posix_time;

      if (observer(time)) {
        ++observer;
        if (write_vtu) {
          pvdwriter->write(time);
        }
        std::cout << "timestep t = " << time << "\n";
      }

      if (observer2(time)) {
        if (write_vtu) {
          // pvdwriter->write_timeseries();
          pvdwriter->write_pvdfile();
        }
        time_duration td = microsec_clock::local_time()-first_seq;
        double measured_time = static_cast<double>(td.total_microseconds())*1.e-6;
        std::cout << "# Statistics: time/timestep = " << measured_time/double(iteration) << " sec\n";
        ++observer2;
      }
    }


    ///
    virtual void beginIteration(size_t iter) {}
    virtual void oneIteration(size_t iter) = 0;
    virtual void firstIteration(size_t iter) { oneIteration(iter); }
    virtual void endIteration(size_t iter) {}


    /// implementation of the timestep loop, including a fixed point iteration
    void solve()
    {
      solveInitialProblem();
      transferInitialSolution();

      time = start_time + timestep;
      iteration = 0;

      // first iteration
      initTimestep();
      for (size_t iter = 0; iter < n_iter; ++iter) {
        beginIteration(iter);
        firstIteration(iter);
        endIteration(iter);
      }
      closeTimestep();
      ++iteration;

      // iteration loop
      while (time < end_time) {
        time += timestep;

        initTimestep();

        // fixed point iteration
        for (size_t iter = 0; iter < n_iter; ++iter) {
          beginIteration(iter);
          oneIteration(iter);
          endIteration(iter);
        }

        closeTimestep();
        ++iteration;
      }
    }


    /// add a solution component and create a plan (and dual plan) for this component
    self& add_data(ScalarValues<Real, Complex>& data)
    {
      solution.push_back(&data);
      return *this;
    }

    self& add_data(VectorValues<Real, Complex>& data)
    {
      solution_vec.push_back(&data);
      return *this;
    }


    /// return \ref filewriter
    VtkWriter& getFilewriter() { return *filewriter; }

    /// return \ref workspace
    shtns_cfg& getWorkspace() { return workspace; }
    shtns_cfg const& getWorkspace() const { return workspace; }

    Complex* getOldSolution(size_t i = 0, int iter = 0)
    {
      size_t j = (iter + n_old_solutions) % n_old_solutions;
      return (*old_solutions[j])[i];
    }

    Complex const* getOldSolution(size_t i = 0, int iter = 0) const
    {
      size_t j = (iter + n_old_solutions) % n_old_solutions;
      return (*old_solutions[j])[i];
    }

    template <int Comp>
    Complex* getOldSolutionVec(size_t i = 0, int iter = 0)
    {
      size_t j = (iter + n_old_solutions) % n_old_solutions;
      return Comp == 0 ? (*old_solutions_vec[j])[i].first : (*old_solutions_vec[j])[i].second;
    }

    template <int Comp>
    Complex const* getOldSolutionVec(size_t i = 0, int iter = 0) const
    {
      size_t j = (iter + n_old_solutions) % n_old_solutions;
      return Comp == 0 ? (*old_solutions_vec[j])[i].first : (*old_solutions_vec[j])[i].second;
    }

    ScalarValues<Real, Complex>& getSolution(size_t i = 0) { return (*solution[i]); }
    ScalarValues<Real, Complex> const& getSolution(size_t i = 0) const { return (*solution[i]); }

    VectorValues<Real, Complex>& getSolutionVec(size_t i = 0) { return (*solution_vec[i]); }
    VectorValues<Real, Complex> const& getSolutionVec(size_t i = 0) const { return (*solution_vec[i]); }

    boost::property_tree::ptree& get_pt() { return pt; }

    Mesh const& getMesh()
    {
      return *mesh;
    }

    VtkTimeseriesWriter& getTimeseriesWriter()
    {
      return *pvdwriter;
    }

  protected:

    void transform(size_t idx) { spat_to_SH(*solution[idx]); }
    void transform_vec(size_t idx) { spat_to_VSH(*solution_vec[idx]); }

    void back_transform(size_t idx) { SH_to_spat(*solution[idx]); }
    void back_transform_vec(size_t idx) { VSH_to_spat(*solution_vec[idx]); }

    void transform(ScalarValues<Real, Complex>& U) const { spat_to_SH(U); }
    void transform(VectorValues<Real, Complex>& U) const { spat_to_VSH(U); }

    void back_transform(ScalarValues<Real, Complex>& U) const { SH_to_spat(U); }
    void back_transform(VectorValues<Real, Complex>& U) const { VSH_to_spat(U); }

  protected:
    std::string name;

    shtns_cfg workspace;
    Observer observer, observer2;
    std::unique_ptr<VtkWriter> filewriter;
    std::unique_ptr<VtkTimeseriesWriter> pvdwriter;

    double radius;
    double time;
    double start_time;
    double end_time;
    double timestep;
    double tol;

    size_t size;
    size_t size_dual;
    size_t iteration;
    size_t n_iter;
    size_t n_threads;
    unsigned int default_flag, detailed_flag;

    bool use_old_solution;
    bool write_vtu = true;

    mutable std::vector<ScalarValues<Real, Complex>*> solution;
    mutable std::vector<VectorValues<Real, Complex>*> solution_vec;
    mutable std::vector<std::vector<Complex*>*> old_solutions;
    mutable std::vector<std::vector<std::pair<Complex*,Complex*> >*> old_solutions_vec;
    size_t n_old_solutions;

    std::vector<Real> laplace, div, rot;
    ConstantCoeff<Real> grad, rot_n, rot_inv;

    boost::property_tree::ptree pt;
    boost::posix_time::ptime first_seq;

    std::unique_ptr<Mesh> mesh;
  };

} // end namespace Shtns
