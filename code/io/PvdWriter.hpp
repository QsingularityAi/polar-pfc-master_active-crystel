#pragma once

#include <map>
#include <fstream>
#include <string>

#include "VtkWriter.hpp"
#include "utility/Filesystem.hpp"
#include "utility/Timer.hpp"

namespace Shtns
{
  class PvdWriter
  {
  public:
    PvdWriter(VtkWriter& filewriter)
      : filewriter(filewriter)
    {
      using namespace Shtns::filesystem;

      auto p = path(filewriter.filename());
      filename = p.stem().string();
      directory = p.remove_filename();

      auto d = path(directory) /= path("data");
      create_directories(d);
    }

    void write(double t)
    {
      using namespace Shtns::filesystem;

      Timer timer;
      auto fn_time = path("data") /= path(filename + "_" + std::to_string(t) + ".vtu");

      auto fn_vtu = path(directory) /= fn_time;
      filewriter.filename(fn_vtu.string());
      filewriter.write();

      auto fn_pvd = path(directory) /= path(filename + ".pvd");
      write_pvd_file(fn_pvd.string(), fn_time.string(), t);

      msg("Write file ",fn_time," needed ",timer.elapsed()," sec");
    }

    template <class Observer>
    void write(Observer& observer, double t)
    {
      if (observer(t)) {
        write(t);
        ++observer;
      }
    }

  protected:
    void write_pvd_file(std::string pvd_filename, std::string vtu_filename, double t)
    {
      std::ofstream out(pvd_filename);
      out << "<?xml version=\"1.0\"?>\n";
      out << "<VTKFile type=\"Collection\" version=\"0.1\" >\n<Collection>\n";

      for (auto const& timestep : timestep_list)
        out << "<DataSet timestep=\"" << timestep.first << "\" part=\"0\" file=\"" << timestep.second << "\"/>\n";

      out << "<DataSet timestep=\"" << t << "\" part=\"0\" file=\"" << vtu_filename << "\"/>\n";
      out << "</Collection>\n</VTKFile>";

      timestep_list[t] = vtu_filename;
    }

  private:
    VtkWriter& filewriter;

    std::string filename;
    Shtns::filesystem::path directory;
    std::map<double, std::string> timestep_list;
  };

} // end namespace Shtns

