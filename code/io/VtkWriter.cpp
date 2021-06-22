#include "VtkWriter.hpp"

#include <sstream>
#include <fstream>
#include <iterator>
#include <string>
#include <regex>

#ifdef SHTNS_HAS_COMPRESSION
#include <zlib.h>
#endif

#include "Mesh.hpp"
#include "Output.hpp"
#include "Shtns.hpp"
#include "utility/Filesystem.hpp"
#include "utility/String.hpp"

namespace Shtns {

VtkWriter::VtkWriter(std::string const& filename, FormatTypes format, DataTypes datatype)
  : _filename(filename)
  , _format(format)
  , _datatype(datatype)
{
#ifndef SHTNS_HAS_COMPRESSION
  if (_format == COMPRESSED) {
    msg("dec-pde is compiled without compression. Falling back to BINARY VTK output!");
    _format = BINARY;
  }
#endif

#ifdef SHTNS_HAS_MPI
  int rank = -1;
  int num_ranks = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  if (num_ranks > 1) {
    auto p = filesystem::path(_filename);
    auto fn = p.stem();
    p.remove_filename();
    p /= fn.string() + "_p" + std::to_string(rank) + ".vtu";
    _filename = p.string();
  }
#endif
}


VtkWriter& VtkWriter::attach(DOFVector const& data)
{
  assert_msg(_mesh == nullptr || _mesh == &data.mesh(), "Can not write multiple meshes!");
  _mesh = &data.mesh();
  _values.push_back(&data);
  return *this;
}


void VtkWriter::write(Mesh const& mesh,
                      std::vector<DOFVector const*> const& values) const
{
  std::ofstream out(_filename, std::ios_base::ate | std::ios::binary);
  if (_format == ASCII) {
    if (_datatype == FLOAT32)
      out << std::setprecision(std::numeric_limits<float>::digits10);
    else
      out << std::setprecision(std::numeric_limits<double>::digits10);
  }

  std::vector<pos_type> offsets; // pos => offset
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
      << "byte_order=\"" << get_endian() << "\" header_type=\"" << to_string<size_type>() << "\""
      << (_format == COMPRESSED ? " compressor=\"vtkZLibDataCompressor\">\n" : ">\n");
  out << "<UnstructuredGrid>\n";
  out << "<Piece NumberOfPoints=\"" << mesh.nodes().size() << "\" "
      << "NumberOfCells=\"" << mesh.cells().size() << "\">\n";

  auto itScalar = std::find_if(values.begin(), values.end(),
                         [](DOFVector const* v){ return v->type() == DOFVector::NODE_DATA && v->shape() == DOFVector::SCALAR; });
  auto itVector = std::find_if(values.begin(), values.end(),
                         [](DOFVector const* v){ return v->type() == DOFVector::NODE_DATA && v->shape() == DOFVector::VECTOR; });

  out << "<PointData" << (itScalar != values.end() ? " Scalars=\"" + (*itScalar)->name() + "\"" : "")
                      << (itVector != values.end() ? " Vectors=\"" + (*itVector)->name() + "\"" : "")
                      << ">\n";
  for (auto const& v : values)
    if (v->type() == DOFVector::NODE_DATA)
      write_data(mesh, out, offsets, *v);
  out << "</PointData>\n";

  auto it2Scalar = std::find_if(values.begin(), values.end(),
                         [](DOFVector const* v){ return v->type() != DOFVector::NODE_DATA && v->shape() == DOFVector::SCALAR; });
  auto it2Vector = std::find_if(values.begin(), values.end(),
                         [](DOFVector const* v){ return v->type() != DOFVector::NODE_DATA && v->shape() == DOFVector::VECTOR; });
  out << "<CellData" << (it2Scalar != values.end() ? " Scalars=\"" + (*it2Scalar)->name() + "\"" : "")
                     << (it2Vector != values.end() ? " Vectors=\"" + (*it2Vector)->name() + "\"" : "")
                     << ">\n";
  for (auto const& v : values)
    if (v->type() != DOFVector::NODE_DATA)
      write_data(mesh, out, offsets, *v);
  out << "</CellData>\n";

  out << "<Points>\n";
  write_points(out, offsets, mesh.nodes());
  out << "</Points>\n";

  out << "<Cells>\n";
  write_cells(out, offsets, mesh.cells());
  out << "</Cells>\n";
  out << "</Piece>\n";
  out << "</UnstructuredGrid>\n";

  std::vector<size_type> blocks; // size of i'th appended block
  pos_type appended_pos = 0;
  if (is_a(_format, APPENDED)) {
    out << "<AppendedData encoding=\"raw\">\n_";
    appended_pos = out.tellp();
    for (auto const& v : values) {
      if (v->type() == DOFVector::NODE_DATA) {
        if (_datatype == FLOAT32)
          blocks.push_back( write_data_appended<float>(mesh, out, *v) );
        else
          blocks.push_back( write_data_appended<double>(mesh, out, *v) );
      }
    }
    for (auto const& v : values) {
      if (v->type() != DOFVector::NODE_DATA) {
        if (_datatype == FLOAT32)
          blocks.push_back( write_data_appended<float>(mesh, out, *v) );
        else
          blocks.push_back( write_data_appended<double>(mesh, out, *v) );
      }
    }
    if (_datatype == FLOAT32)
      blocks.push_back( write_points_appended<float>(out, mesh.nodes()) );
    else
      blocks.push_back( write_points_appended<double>(out, mesh.nodes()) );
    auto bs = write_cells_appended(out, mesh.cells());
    blocks.insert(blocks.end(), bs.begin(), bs.end());
    out << "</AppendedData>\n";
  }

  out << "</VTKFile>";

  // fillin offset values and block sizes
  if (is_a(_format, APPENDED)) {
    pos_type offset = 0;
    for (size_t i = 0; i < offsets.size(); ++i) {
      out.seekp(offsets[i]);
      out << '"' << offset << '"';
      offset += pos_type(blocks[i]);
    }
  }
}


void VtkWriter::write_data(Mesh const& mesh, std::ofstream& out, std::vector<pos_type>& offsets,
                           DOFVector const& vec) const
{
  out << "<DataArray Name=\"" << vec.name() << "\" type=\"" << to_string(_datatype) << "\""
      << " NumberOfComponents=\"" << vec.numComponents() << "\" format=\"" << (_format == ASCII ? "ascii\">\n" : "appended\"");

  if (_format == ASCII) {
    size_t i = 0;
    for (auto const& v : vec) {
      for (auto const& vi : v)
        out << vi << (++i % 6 != 0 ? ' ' : '\n');
    }
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";
  } else {
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


void VtkWriter::write_points(std::ofstream& out, std::vector<pos_type>& offsets,
                             Mesh::NodeRange const& nodes) const
{
  out << "<DataArray type=\"" << to_string(_datatype) << "\""
      << " NumberOfComponents=\"3\" format=\"" << (_format == ASCII ? "ascii\">\n" : "appended\"");

  if (_format == ASCII) {
    size_t i = 0;
    for (auto const& v : nodes) {
      for (auto const& v_i : v)
        out << v_i << (++i % 6 != 0 ? ' ' : '\n');
    }
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";
  } else {
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


void VtkWriter::write_cells(std::ofstream& out, std::vector<pos_type>& offsets,
                            Mesh::CellRange const& cells) const
{
  if (_format == ASCII) {
    out << "<DataArray type=\"" << to_string<index_type>() << "\" Name=\"connectivity\" format=\"ascii\">\n";
    size_t i = 0;
    for (auto const& c : cells) {
      for (auto const& c_i : c)
        out << c_i << (++i % 6 != 0 ? ' ' : '\n');
    }
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";

    out << "<DataArray type=\"" << to_string<index_type>() << "\" Name=\"offsets\" format=\"ascii\">\n";
    i = 0;
    size_t old_o = 0;
    for (auto const& c : cells)
      out << (old_o += Element::nodes(cells.type)) << (++i % 6 != 0 ? ' ' : '\n');
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";

    out << "<DataArray type=\"" << to_string<uint8_t>() << "\" Name=\"types\" format=\"ascii\">\n";
    i = 0;
    for (auto const& c : cells)
      out << uint16_t(type_map(cells.type)) << (++i % 6 != 0 ? ' ' : '\n');
    out << (i % 6 != 0 ? "\n" : "") << "</DataArray>\n";
  } else {
    out << "<DataArray type=\"" << to_string<index_type>() << "\" Name=\"connectivity\" format=\"appended\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"" << to_string<index_type>() << "\" Name=\"offsets\" format=\"appended\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";

    out << "<DataArray type=\"" << to_string<uint8_t>() << "\" Name=\"types\" format=\"appended\"";
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<size_type>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class T>
VtkWriter::size_type VtkWriter::write_values_to_buffer(size_t max_num_values, unsigned char* buffer,
                                                       std::vector<T> const& vec, size_t shift) const
{
  size_t num_values = std::min(max_num_values, vec.size()-shift);
  size_type bs = num_values*sizeof(T);
  std::memcpy(buffer, (unsigned char*)(vec.data()+shift), size_t(bs));
  return bs;
}


template <class OStream>
VtkWriter::size_type VtkWriter::write_compressed(unsigned char const* buffer, unsigned char* buffer_out,
                                                 size_type bs, size_type cbs, int level, OStream& outb) const
{
#ifdef SHTNS_HAS_COMPRESSION
  uLongf uncompressed_space = uLongf(bs);
  uLongf compressed_space = uLongf(cbs);

  Bytef* out = reinterpret_cast<Bytef*>(buffer_out);
  Bytef const* in = reinterpret_cast<Bytef const*>(buffer);

  if (compress2(out, &compressed_space, in, uncompressed_space, level) != Z_OK)
    error_exit("Zlib error while compressing data.");
  else
    outb.write((char*)out, compressed_space);

  return compressed_space;
#else
  error_exit("Can not call write_compressed without compression enabled!");
  return 0;
#endif
}


template <class T>
VtkWriter::size_type VtkWriter::write_appended(std::ofstream& out, std::vector<T> const& values) const
{
  assert_msg(is_a(_format, APPENDED), "Function should by called only in appended mode!");
  pos_type begin_pos = out.tellp();

  size_type size = values.size() * sizeof(T);

  size_type num_full_blocks = size / block_size;
  size_type last_block_size = size % block_size;
  size_type num_blocks = num_full_blocks + (last_block_size > 0 ? 1 : 0);

  // write block-size(s)
  size_type zero = 0;
  if (_format == COMPRESSED) {
    out.write((char*)&num_blocks, sizeof(size_type));
    out.write((char*)&block_size, sizeof(size_type));
    out.write((char*)&last_block_size, sizeof(size_type));
    for (size_type i = 0; i < num_blocks; ++i)
      out.write((char*)&zero, sizeof(size_type));
  } else {
    out.write((char*)&size, sizeof(size_type));
  }

  size_type compressed_block_size = block_size + (block_size + 999)/1000 + 12;
  std::vector<unsigned char> buffer(block_size);
  std::vector<unsigned char> buffer_out;
  size_t num_values = block_size / sizeof(T);

  std::vector<size_type> cbs(size_t(num_blocks), 0); // compressed block sizes
  for (size_t i = 0; i < size_t(num_blocks); ++i) {
    size_type bs = write_values_to_buffer<T>(num_values, buffer.data(), values, i*num_values);

    if (_format == COMPRESSED) {
      buffer_out.resize(size_t(compressed_block_size));
      cbs[i] = write_compressed(buffer.data(), buffer_out.data(), bs,
                                compressed_block_size, compression_level, out);
    } else
      out.write((char*)buffer.data(), bs);
  }

  pos_type end_pos = out.tellp();
  if (_format == COMPRESSED) {
    out.seekp(begin_pos + int64_t(3*sizeof(size_type)));
    out.write((char*)cbs.data(), num_blocks*sizeof(size_type));
    out.seekp(end_pos);
  }

  return size_type(end_pos - begin_pos);
}


std::vector<VtkWriter::size_type> VtkWriter::write_cells_appended(std::ofstream& out, Mesh::CellRange const& cells) const
{
  assert_msg(is_a(_format, APPENDED), "Function should by called only in appended mode!");

  // write conncetivity
  std::vector<index_type> connectivity;
  connectivity.reserve(cells.size() * (*cells.begin()).size());
  for (auto const& c : cells) {
    for (auto const& ci : c)
      connectivity.push_back( index_type(ci) );
  }
  size_type bs0 = write_appended(out, connectivity);

  // write offsets
  std::vector<index_type> cell_offsets;
  cell_offsets.reserve(cells.size());
  index_type old_o = 0;
  for (auto const& c : cells)
    cell_offsets.push_back(old_o += Element::nodes(cells.type));
  size_type bs1 = write_appended(out, cell_offsets);

  // write cell types
  std::vector<uint8_t> types;
  types.reserve(cells.size());
  for (auto const& c : cells)
    types.push_back(type_map(cells.type));
  size_type bs2 = write_appended(out, types);

  return {bs0, bs1, bs2};
}

} // end namespace Dec
