#pragma once

#include <map>
#include <fstream>
#include <string>

#include "VtkWriter.hpp"
#include "utility/Filesystem.hpp"

namespace Shtns
{
  class VtkTimeseriesWriter
  {
    using pos_type = VtkWriter::pos_type;
    using size_type = VtkWriter::size_type; // used in block headers
    using index_type = VtkWriter::index_type; // used in connectivity and offsets

  public:
    VtkTimeseriesWriter(VtkWriter& filewriter)
      : filewriter(filewriter)
      , _datatype(filewriter.datatype())
      , _format(filewriter.format())
    {
      using namespace Shtns::filesystem;

      auto p = path(filewriter.filename());
      filename = p.stem().string();
      directory = p.remove_filename();

      auto d = path(directory) /= path("data");
      create_directories(d);
    }

    // write in the first call the mesh as binary file and in all other calls
    // the data as binary files
    void write(double t);

    template <class Observer>
    void write(Observer& observer, double t)
    {
      if (observer(t)) {
        write(t);
        ++observer;
      }
    }

    // create a timeseries file, i.e. a vtu file containing all the timesteps
    void write_timeseries() const;

    // create pvd file, i.e. a sequence of vtu files and a container file that redirects to the timestep files
    void write_pvdfile(bool force = false) const;


  protected:
    size_type write_data(std::ofstream& out, DOFVector const& values) const
    {
      if (_datatype == VtkWriter::FLOAT32)
        return filewriter.write_data_appended<float>(filewriter.mesh(), out, values);
      else
        return filewriter.write_data_appended<double>(filewriter.mesh(), out, values);
    }

    size_type write_points(std::ofstream& out) const
    {
      if (_datatype == VtkWriter::FLOAT32)
        return filewriter.write_points_appended<float>(out, filewriter.mesh().nodes());
      else
        return filewriter.write_points_appended<double>(out, filewriter.mesh().nodes());
    }

    std::vector<size_type> write_cells(std::ofstream& out) const
    {
      return filewriter.write_cells_appended(out, filewriter.mesh().cells());
    }

  private:
    VtkWriter& filewriter;

    VtkWriter::DataTypes _datatype;
    VtkWriter::FormatTypes _format;

    std::string filename;
    Shtns::filesystem::path directory;
    std::map<double, std::string> timestep_list;

    // stores block_sizes
    size_type bs_points;
    std::vector<size_type> bs_cells;
    std::vector<size_type> bs_data;

    std::string filename_mesh;

    bool initialized = false;
  };

} // end namespace Shtns

