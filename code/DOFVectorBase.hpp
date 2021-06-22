#pragma once

#include <complex>
#include <string>
#include <vector>

#include "container.hpp"
#include "Mesh.hpp"

namespace Shtns
{
  class DOFVectorBase
  {
  protected:
    using Real = double;
    using Complex = std::complex<double>;

  public:
    enum DataType {
      CELL_DATA, EDGE_DATA, NODE_DATA
    };

    enum Shape {
      SCALAR, VECTOR
    };

  public:
    DOFVectorBase(Mesh const& mesh, std::string name, DataType type, Shape shape)
      : _mesh(&mesh)
      , _type(type)
      , _name(name)
      , _shape(shape)
    {}

    DOFVectorBase(DOFVectorBase const&) = default;
    DOFVectorBase(DOFVectorBase&&) = default;

    DOFVectorBase& operator=(DOFVectorBase const&) = default;
    DOFVectorBase& operator=(DOFVectorBase&&) = default;

    size_t size() const
    {
      return mesh().vertices();
    }

    size_t num_rows() const
    {
      return size();
    }

    size_t num_cols() const
    {
      return 1;
    }

    Mesh const& mesh() const
    {
      return *_mesh;
    }

    std::string name() const
    {
      return _name;
    }

    DataType type() const
    {
      return _type;
    }

    Shape shape() const
    {
      return _shape;
    }

  private:
    Mesh const* _mesh;
    std::string _name;
    DataType    _type;
    Shape      _shape;
  };



  class DOFVector
      : public DOFVectorBase
  {
    struct iterator
    {
    private:
      shtns_cfg workspace;
      size_t i = 0, max_i = 0;
      std::vector<double> x;

      ScalarValues<Real, Complex> const* scalar = nullptr;
      VectorValues<Real, Complex> const* vector = nullptr;

    public:
      // begin-Constructor for scalar containers
      iterator(shtns_cfg w, ScalarValues<Real, Complex> const* s)
        : workspace(w)
        , max_i(w->nphi * w->nlat)
        , x(1)
        , scalar(s)
      {
        update();
      }

      // begin-Constructor for vector-valued containers
      iterator(shtns_cfg w, VectorValues<Real, Complex> const* v)
        : workspace(w)
        , max_i(w->nphi * w->nlat)
        , x(3)
        , vector(v)
      {
        update();
      }

      // end-Constructor for scalar containers
      iterator(shtns_cfg w, ScalarValues<Real, Complex> const* s, bool /*end*/)
        : workspace(w)
        , i(w->nphi * w->nlat)
        , max_i(w->nphi * w->nlat)
        , scalar(s)
      {}

      // end-Constructor for vector-valued containers
      iterator(shtns_cfg w, VectorValues<Real, Complex> const* v, bool /*end*/)
        : workspace(w)
        , i(w->nphi * w->nlat)
        , max_i(w->nphi * w->nlat)
        , vector(v)
      {}

      iterator(iterator const&) = default;
      iterator(iterator&&) = default;

      iterator& operator++()
      {
        ++i;
        update();
        return *this;
      }

      iterator  operator++(int) { iterator tmp(*this); operator++(); return tmp; }

      bool operator==(const iterator& rhs) { return i == rhs.i; }
      bool operator!=(const iterator& rhs) { return !(this->operator==(rhs)); }

      std::vector<double> const& operator*() const { return x; }

    private:

      void update()
      {
        if (i < max_i) {
          if (scalar)
            x[0] = (*scalar)[i];
          else if (vector) {
            size_t _i = i / workspace->nlat;
            size_t _j = i % workspace->nlat;
            sphere_to_cartesian(workspace, vector->values_t(), vector->values_p(), _i, _j, x);
          }
        }
      }
    };

  public:
    DOFVector(Mesh const& mesh, ScalarValues<Real, Complex> const& s)
      : DOFVectorBase(mesh, s.name(), NODE_DATA, SCALAR)
      , scalar(&s)
    {}

    DOFVector(Mesh const& mesh, VectorValues<Real, Complex> const& v)
      : DOFVectorBase(mesh, v.name(), NODE_DATA, VECTOR)
      , vector(&v)
      , _numComponents(v.num_components())
    {}

    size_t numComponents() const
    {
      return _numComponents;
    }

    iterator begin() const
    {
      if (scalar)
        return {scalar->shtns(), scalar};
      else if (vector)
        return {vector->shtns(), vector};
      std::abort();
    }

    iterator end() const
    {
      if (scalar)
        return {scalar->shtns(), scalar, true};
      else if (vector)
        return {vector->shtns(), vector, true};
      std::abort();
    }

  private:
    ScalarValues<Real, Complex> const* scalar = nullptr;
    VectorValues<Real, Complex> const* vector = nullptr;

    size_t _numComponents = 1;
  };

} // end namespace Shtns
