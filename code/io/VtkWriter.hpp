#pragma once

#include <iosfwd>
#include <map>

#include "Element.hpp"
#include "FileWriter.hpp"

namespace Shtns
{
  /// File-Writer for Vtk .vtu files
  class VtkWriter : public FileWriter
  {
    using pos_type = typename std::ostream::pos_type;
    using size_type = uint32_t; // uint64_t
    using index_type = int32_t; // int64_t

    friend class PvdWriter;
    friend class VtkTimeseriesWriter;

  public:
    enum FormatTypes {
      ASCII      = 1<<0,
      BINARY     = 1<<1,
      COMPRESSED = 1<<2,
      APPENDED = BINARY | COMPRESSED
    };

    enum DataTypes {
      FLOAT32 = 32,
      FLOAT64 = 64
    };

  public:
    /// Constructor, takes the \p filename of a .vtu file as argument
    VtkWriter(std::string const& filename, FormatTypes format = BINARY, DataTypes datatype = FLOAT32);

    /// Implementation of \ref FileWriter::write
    virtual void write(Mesh const& mesh,
                       std::vector<DOFVector const*> const& values = {}) const override;

    /// Write the attached data to the set filename
    void write() const
    {
      assert_msg(_mesh != nullptr, "A mesh must be attached to the VtkWriter!");
      write(*_mesh, _values);
    }

    /// Returns the currently set filename
    std::string filename() const
    {
      return _filename;
    }

    /// Update the filename
    void filename(std::string const& fn)
    {
      _filename = fn;
    }

    /// Attach a mesh to the writer
    VtkWriter& attach(Mesh const& mesh)
    {
      _mesh = &mesh;
      return *this;
    }

    /// Attach node data to the writer
    VtkWriter& attach(DOFVector const& data);

    Mesh const& mesh() const
    {
      assert( _mesh != nullptr );
      return *_mesh;
    }

    std::vector<DOFVector const*> const& values() const
    {
      return _values;
    }

    DataTypes datatype() const
    {
      return _datatype;
    }

    FormatTypes format() const
    {
      return _format;
    }

  private:
    void write_data(Mesh const& mesh, std::ofstream& out, std::vector<pos_type>& offsets,
                    DOFVector const& values) const;

    void write_points(std::ofstream& out, std::vector<pos_type>& offsets,
                      Mesh::NodeRange const& nodes) const;

    void write_cells(std::ofstream& out, std::vector<pos_type>& offsets,
                     Mesh::CellRange const& cells) const;

    template <class T>
    size_type write_appended(std::ofstream& out, std::vector<T> const& values) const;

    template <class T>
    size_type write_data_appended(Mesh const& mesh, std::ofstream& out, DOFVector const& values) const;

    template <class T>
    size_type write_points_appended(std::ofstream& out, Mesh::NodeRange const& nodes) const;

    std::vector<size_type> write_cells_appended(std::ofstream& out, Mesh::CellRange const& cells) const;

    // Returns endianness
    static std::string get_endian()
    {
      short i = 1;
      return (reinterpret_cast<char*>(&i)[1] == 1 ? "BigEndian" : "LittleEndian");
    }

    static uint8_t type_map(Element::Types T)
    {
      return T == Element::VERTEX    ? 1 :
             T == Element::LINE2     ? 3 :
             T == Element::TRIANGLE3 ? 5 :
             T == Element::QUAD4     ? 9 :
             T == Element::LINE3     ? 21 :
             T == Element::TRIANGLE6 ? 22 :
             T == Element::QUAD8     ? 23 : 0;
    }

    template <class T>
    static std::string to_string()
    {
      return std::is_same<T, int8_t>::value   ? "Int8"   :
             std::is_same<T, int16_t>::value  ? "Int16"  :
             std::is_same<T, int32_t>::value  ? "Int32"  :
             std::is_same<T, int64_t>::value  ? "Int64"  :
             std::is_same<T, uint8_t>::value  ? "UInt8"  :
             std::is_same<T, uint16_t>::value ? "UInt16" :
             std::is_same<T, uint32_t>::value ? "UInt32" :
             std::is_same<T, uint64_t>::value ? "UInt64" : "Unknown";
    }

    static std::string to_string(DataTypes t)
    {
      return t == FLOAT32 ? "Float32" :
             t == FLOAT64 ? "Float64" : "Unknown";
    }

  private:
    template <class T>
    size_type write_values_to_buffer(size_t max_num_values, unsigned char* buffer,
                                      std::vector<T> const& vec, size_t shift) const;
    template <class OStream>
    size_type write_compressed(unsigned char const* buffer, unsigned char* buffer_out,
                                size_type bs, size_type cbs, int level, OStream& outb) const;

  private:
    mutable std::string _filename;
    FormatTypes         _format;
    DataTypes           _datatype;

    // attached data
    Mesh const*                   _mesh = nullptr;
    std::vector<DOFVector const*> _values;

    size_t const block_size = 1024*32;
    int compression_level = -1; // in [0,9], -1 ... use default value
  };


  template <class T>
  VtkWriter::size_type VtkWriter::write_data_appended(Mesh const& mesh, std::ofstream& out, DOFVector const& vec) const
  {
    assert_msg(is_a(_format, APPENDED), "Function should by called only in appended mode!");

    std::vector<T> data;
    data.reserve(vec.size() * vec.numComponents());
    for (auto const& v : vec) {
      for (auto const& vi : v)
        data.push_back(T(vi));
    }
    return write_appended(out, data);
  }


  template <class T>
  VtkWriter::size_type VtkWriter::write_points_appended(std::ofstream& out, Mesh::NodeRange const& nodes) const
  {
    assert_msg(is_a(_format, APPENDED), "Function should by called only in appended mode!");

    std::vector<T> data;
    data.reserve(nodes.size() * (*nodes.begin()).size());
    for (auto const& v : nodes) {
      for (auto const& vi : v)
        data.push_back( T(vi) );
    }
    return write_appended(out, data);
  }

} // end namespace Shtns
