#include "VtkTimeseriesWriter.hpp"

#include "utility/Timer.hpp"

namespace Shtns {


void VtkTimeseriesWriter::write(double t)
{
  using namespace Shtns::filesystem;
  Timer timer;

  if (!initialized) {
    auto fn_mesh = path("data") /= path(filename + "_mesh.data");
    auto path_mesh = path(directory) /= fn_mesh;
    std::ofstream file_mesh(path_mesh.string(), std::ios_base::out);

    bs_points = write_points(file_mesh);
    bs_cells = write_cells(file_mesh);

    filename_mesh = path_mesh.string();
    initialized = true;
  }

  auto fn_time = path("data") /= path(filename + "_" + std::to_string(t) + ".data");
  auto path_data = path(directory) /= fn_time;

  std::ofstream file_data(path_data.string(), std::ios_base::out);

  for (auto const& v : filewriter.values()) {
    if (v->type() == DOFVector::NODE_DATA) {
      bs_data.push_back( write_data(file_data, *v) );
    }
  }

  timestep_list[t] = path_data.string();

  msg("Write file ",fn_time," needed ",timer.elapsed()," sec");
}


void VtkTimeseriesWriter::write_timeseries() const
{
  using namespace Shtns::filesystem;
  auto fn_pvd = path(directory) /= path(filename + ".timeseries.vtu");


  Mesh const& mesh = filewriter.mesh();

  std::map<double, std::vector<pos_type>> offsets; // time -> (pos => offset)

  std::ofstream out(fn_pvd.string(), std::ios_base::out);
  out << "<?xml version=\"1.0\"?>\n";

  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
      << "byte_order=\"" << VtkWriter::get_endian() << "\" header_type=\"" << VtkWriter::to_string<size_type>() << "\""
      << (_format == VtkWriter::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\">\n" : ">\n");
  out << "<UnstructuredGrid TimeValues=\"\n";
  for (auto const& timestep : timestep_list)
    out << timestep.first << "\n";
  out << "\">\n";

  out << "<Piece NumberOfPoints=\"" << mesh.nodes().size() << "\" "
      << "NumberOfCells=\"" << mesh.cells().size() << "\">\n";

  // write vertices
  out << "<Points>\n";
  size_t i = 0;
  for (auto const& timestep : timestep_list) {
    out << "<DataArray type=\"" << VtkWriter::to_string(_datatype) << "\""
        << " NumberOfComponents=\"3\" TimeStep=\"" << i << "\" format=\"appended\" offset=";
    offsets[timestep.first].push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";
    ++i;
  }
  out << "</Points>\n";

  // write cells
  out << "<Cells>\n";
  i = 0;
  for (auto const& timestep : timestep_list) {
    out << "<DataArray type=\"" << VtkWriter::to_string<index_type>() << "\""
        << " Name=\"connectivity\" TimeStep=\"" << i << "\" format=\"appended\" offset=";
    offsets[timestep.first].push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"" << VtkWriter::to_string<index_type>() << "\""
        << " Name=\"offsets\" TimeStep=\"" << i << "\" format=\"appended\" offset=";
    offsets[timestep.first].push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"" << VtkWriter::to_string<uint8_t>() << "\" Name=\"types\" TimeStep=\"" << i << "\" format=\"appended\" offset=";
    offsets[timestep.first].push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";
    ++i;
  }
  out << "</Cells>\n";

  // write data vectors
  auto itScalar = std::find_if(filewriter.values().begin(), filewriter.values().end(),
                        [](DOFVector const* v){ return v->type() == DOFVector::NODE_DATA && v->shape() == DOFVector::SCALAR; });
  auto itVector = std::find_if(filewriter.values().begin(), filewriter.values().end(),
                        [](DOFVector const* v){ return v->type() == DOFVector::NODE_DATA && v->shape() == DOFVector::VECTOR; });

  out << "<PointData" << (itScalar != filewriter.values().end() ? " Scalars=\"" + (*itScalar)->name() + "\"" : "")
                      << (itVector != filewriter.values().end() ? " Vectors=\"" + (*itVector)->name() + "\"" : "")
                      << ">\n";
  i = 0;
  for (auto const& timestep : timestep_list) {
    for (auto const& v : filewriter.values()) {
      if (v->type() == DOFVector::NODE_DATA) {

        out << "<DataArray Name=\"" << v->name() << "\" TimeStep=\"" << i << "\" type=\"" << VtkWriter::to_string(_datatype) << "\""
            << " NumberOfComponents=\"" << v->numComponents() << "\" format=\"appended\" offset=";
        offsets[timestep.first].push_back(out.tellp());
        out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
        out << "/>\n";
      }
    }
    ++i;
  }
  out << "</PointData>\n";

  out << "</Piece>\n";
  out << "</UnstructuredGrid>\n";


  pos_type appended_pos = 0;
  out << "<AppendedData encoding=\"raw\">\n_";
  appended_pos = out.tellp();

  std::ifstream file_mesh(filename_mesh, std::ios_base::in | std::ios_base::binary);
  out << file_mesh.rdbuf();
  file_mesh.close();
  assert( out.tellp() == appended_pos + pos_type(bs_points) + pos_type(bs_cells[0]) + pos_type(bs_cells[1]) + pos_type(bs_cells[2]) );

  for (auto const& timestep : timestep_list) {
    std::ifstream file(timestep.second, std::ios_base::in | std::ios_base::binary);
    out << file.rdbuf();
  }
  out << "</AppendedData>\n";
  out << "</VTKFile>";

  // write correct offsets in file.
  pos_type offset = 0;
  for (auto const& timestep : timestep_list) {
    offset = 0;
    auto const& off = offsets[timestep.first];

    out.seekp(off[0]);
    out << '"' << 0 << '"';
    offset += pos_type(bs_points);
    for (size_t i = 0; i < bs_cells.size(); ++i) { // write offsets of cells only
      out.seekp(off[1+i]);
      out << '"' << offset << '"';
      offset += pos_type(bs_cells[i]);
    }
  }

  size_t shift = bs_cells.size() + 1;
  size_t j = 0;
  for (auto const& timestep : timestep_list) {
    auto const& off = offsets[timestep.first];

    for (size_t i = shift; i < off.size(); ++i) {
      out.seekp(off[i]);
      out << '"' << offset << '"';
      offset += pos_type(bs_data[j++]);
    }
  }
}


void VtkTimeseriesWriter::write_pvdfile(bool force) const
{
  using namespace Shtns::filesystem;

  Mesh const& mesh = filewriter.mesh();

  size_t data_idx = 0;
  for (auto const& timestep : timestep_list) {

    auto p = path(timestep.second);
    auto fn_vtu = path("data") /= path(p.stem().string() + ".vtu");
    auto path_vtu = path(directory) /= fn_vtu;

    if (!force && exists(path_vtu)) {
      for (auto const& v : filewriter.values())
        if (v->type() == DOFVector::NODE_DATA)
          ++data_idx;
      continue;
    }

    std::vector<pos_type> offsets; // pos => offset

    std::ofstream out(path_vtu.string(), std::ios_base::out);
    out << "<?xml version=\"1.0\"?>\n";

    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
        << "byte_order=\"" << VtkWriter::get_endian() << "\" header_type=\"" << VtkWriter::to_string<size_type>() << "\""
        << (_format == VtkWriter::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\">\n" : ">\n");
    out << "<UnstructuredGrid>\n";

    out << "<Piece NumberOfPoints=\"" << mesh.nodes().size() << "\" "
        << "NumberOfCells=\"" << mesh.cells().size() << "\">\n";

    // write vertices
    out << "<Points>\n";
    out << "<DataArray type=\"" << VtkWriter::to_string(_datatype) << "\""
        << " NumberOfComponents=\"3\" format=\"appended\" offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";
    out << "</Points>\n";

    // write cells
    out << "<Cells>\n";
    out << "<DataArray type=\"" << VtkWriter::to_string<index_type>() << "\""
        << " Name=\"connectivity\" format=\"appended\" offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"" << VtkWriter::to_string<index_type>() << "\""
        << " Name=\"offsets\" format=\"appended\" offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"" << VtkWriter::to_string<uint8_t>() << "\" Name=\"types\" format=\"appended\" offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";
    out << "</Cells>\n";

    // write data vectors
    auto itScalar = std::find_if(filewriter.values().begin(), filewriter.values().end(),
                          [](DOFVector const* v){ return v->type() == DOFVector::NODE_DATA && v->shape() == DOFVector::SCALAR; });
    auto itVector = std::find_if(filewriter.values().begin(), filewriter.values().end(),
                          [](DOFVector const* v){ return v->type() == DOFVector::NODE_DATA && v->shape() == DOFVector::VECTOR; });

    out << "<PointData" << (itScalar != filewriter.values().end() ? " Scalars=\"" + (*itScalar)->name() + "\"" : "")
                        << (itVector != filewriter.values().end() ? " Vectors=\"" + (*itVector)->name() + "\"" : "")
                        << ">\n";
    for (auto const& v : filewriter.values()) {
      if (v->type() == DOFVector::NODE_DATA) {

        out << "<DataArray Name=\"" << v->name() << "\" type=\"" << VtkWriter::to_string(_datatype) << "\""
            << " NumberOfComponents=\"" << v->numComponents() << "\" format=\"appended\" offset=";
        offsets.push_back(out.tellp());
        out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
        out << "/>\n";
      }
    }
    out << "</PointData>\n";
    out << "</Piece>\n";
    out << "</UnstructuredGrid>\n";

    out << "<AppendedData encoding=\"raw\">\n_";
    pos_type appended_pos = out.tellp();

    std::ifstream file_mesh(filename_mesh, std::ios_base::in | std::ios_base::binary);
    out << file_mesh.rdbuf();
    file_mesh.close();
    assert( out.tellp() == appended_pos + pos_type(bs_points) + pos_type(bs_cells[0]) + pos_type(bs_cells[1]) + pos_type(bs_cells[2]) );

    std::ifstream file_data(timestep.second, std::ios_base::in | std::ios_base::binary);
    out << file_data.rdbuf();
    file_data.close();
    out << "</AppendedData>\n";
    out << "</VTKFile>";

    // write correct offsets in file.
    pos_type offset = 0;
    out.seekp(offsets[0]);
    out << '"' << offset << '"';
    offset += pos_type(bs_points);
    for (size_t i = 0; i < bs_cells.size(); ++i) { // write offsets of cells only
      out.seekp(offsets[1+i]);
      out << '"' << offset << '"';
      offset += pos_type(bs_cells[i]);
    }

    size_t shift = bs_cells.size() + 1;
    for (size_t i = shift; i < offsets.size(); ++i) {
      out.seekp(offsets[i]);
      out << '"' << offset << '"';
      offset += pos_type(bs_data[data_idx++]);
    }

    out.close();
    std::remove(timestep.second.c_str());
  }

  // write pvd file
  auto fn_pvd = path(directory) /= path(filename + ".pvd");

  std::ofstream out(fn_pvd.string(), std::ios_base::out);
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"Collection\" version=\"0.1\" >\n<Collection>\n";

  for (auto const& timestep : timestep_list) {
    auto p = path(timestep.second);
    auto fn_vtu = path("data") /= path(p.stem().string() + ".vtu");
    out << "<DataSet timestep=\"" << timestep.first << "\" part=\"0\" file=\"" << fn_vtu.string() << "\"/>\n";
  }
  out << "</Collection>\n</VTKFile>";
}


} // end namespace Shtns

