#include "polar_pfc2.h"

#include "DOFVectorBase.hpp"

using namespace Shtns;

int main(int argc, char** argv)
{
  std::unique_ptr<PolarPfc> prob;

  if (argc > 1)
    prob.reset(new PolarPfc(std::string(argv[1])));
  else
    prob.reset(new PolarPfc);

  // specify the solution components
  ScalarValues<Real, Complex> psi("psi", prob->getWorkspace());
  VectorValues<Real, Complex> P("P", prob->getWorkspace());
  prob->add_data(psi).add_data(P);

  // initialize the data
  prob->initData(INIT_LAPLACE | INIT_DIVERGENCE | INIT_GRADIENT);

  // configure the file-writer
  DOFVector psi_(prob->getMesh(), psi);
  DOFVector P_(prob->getMesh(), P);
  // DOFVector energyPsi_(prob->getMesh(), *prob->getEnergyPsi());
  // DOFVector energyP_(prob->getMesh(), *prob->getEnergyP());

  prob->getFilewriter().attach(psi_).attach(P_); //.attach(energyPsi_).attach(energyP_);

  // solve the problem
  prob->solve();

  prob->getTimeseriesWriter().write_pvdfile();
  // prob->getTimeseriesWriter().write_timeseries();

  // write file solution to file, in order to restart simulation
//   prob.getFilewriter().write_dual(prob.get_write_filename());
//   prob.getFilewriter().write_dual_ascii(prob.get_write_filename() + "_ascii");

  return 0;
}

