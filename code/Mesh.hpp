#pragma once

#include "Shtns.hpp"
#include "Element.hpp"

namespace Shtns
{

  class Mesh
  {
    friend class VtkWriter;

    /// A range over the nodes of a mesh
    struct NodeRange
    {
      struct iterator
      {
      private:
        shtns_cfg workspace;
        double R;
        size_t max_i;
        size_t idx;
        std::array<double, 3> x;

      public:
        // begin-Constructor
        iterator(shtns_cfg w, double radius)
          : workspace(w)
          , R(radius)
          , max_i(w->nphi * w->nlat)
          , idx(0)
        {
          update();
        }

        // end-Constructor
        iterator(shtns_cfg w, double radius, bool /*end*/)
          : workspace(w)
          , R(radius)
          , max_i(w->nphi * w->nlat)
          , idx(max_i)
        {}

        iterator(iterator const&) = default;
        iterator(iterator&&) = default;

        iterator& operator++()
        {
          ++idx;
          update();
          return *this;
        }

        iterator operator++(int) { iterator tmp(*this); operator++(); return tmp; }

        iterator& operator+=(size_t shift)
        {
          idx += shift;
          update();
          return *this;
        }

        iterator operator+(size_t shift) { iterator tmp(*this); operator+=(shift); return tmp; }

        size_t operator-(iterator const& that) { return idx - that.idx; }

        bool operator==(const iterator& rhs) { return idx == rhs.idx; }
        bool operator!=(const iterator& rhs) { return !(this->operator==(rhs)); }

        std::array<double, 3> const& operator*() const { return x; }

      private:
        void update()
        {
          if (idx < max_i) {
            size_t _i = idx / workspace->nlat;
            size_t _j = idx % workspace->nlat;
            get_coord(workspace, _i, _j, R, x);
          }
        }
      };

    private:
      shtns_cfg workspace;
      double R;

    public:
      /// Constructor
      NodeRange(shtns_cfg w, double radius) : workspace(w), R(radius) {}

      iterator begin() const { return {workspace, R}; }
      iterator end() const { return {workspace, R, true}; }

      size_t size() const
      {
        return (workspace->nphi)*(workspace->nlat);
      }
    };


    /// A range over the cells created during traversal
    struct CellRange
    {
      struct iterator
      {
      private:
        shtns_cfg workspace;
        size_t max_j;
        size_t max_i, idx;
        std::array<size_t, 4> indices;

      public:
        // begin-Constructor
        iterator(shtns_cfg w)
          : workspace(w)
          , max_j(w->nlat-1)
          , max_i(w->nphi * max_j)
          , idx(0)
        {
          update();
        }

        // end-Constructor
        iterator(shtns_cfg w, bool /*end*/)
          : workspace(w)
          , max_j(w->nlat-1)
          , max_i(w->nphi * max_j)
          , idx(max_i)
        {}

        iterator(iterator const&) = default;
        iterator(iterator&&) = default;

        iterator& operator++()
        {
          ++idx;
          update();
          return *this;
        }

        iterator  operator++(int) { iterator tmp(*this); operator++(); return tmp; }

        iterator& operator+=(size_t shift)
        {
          idx += shift;
          update();
          return *this;
        }

        iterator operator+(size_t shift) { iterator tmp(*this); operator+=(shift); return tmp; }

        size_t operator-(iterator const& that) { return idx - that.idx; }

        bool operator==(const iterator& rhs) { return idx == rhs.idx; }
        bool operator!=(const iterator& rhs) { return !(this->operator==(rhs)); }

        std::array<size_t, 4> const& operator*() const { return indices; }

      private:

        void update()
        {
          if (idx < max_i) {
            size_t _i = idx / max_j;
            size_t _j = idx % max_j;
            size_t _ip1 = (_i+1) % (workspace->nphi);

            indices = {get_idx(workspace->nlat, _i,   _j),
                       get_idx(workspace->nlat, _ip1, _j),
                       get_idx(workspace->nlat, _ip1, _j+1),
                       get_idx(workspace->nlat, _i,   _j+1)};
          }
        }
      };

    public:
      shtns_cfg workspace;
      static const Element::Types type = Element::QUAD4;

      CellRange(shtns_cfg w) : workspace(w) {}

      iterator begin() const { return {workspace}; }
      iterator end() const { return {workspace, true}; }

      size_t size() const
      {
        return (workspace->nphi)*(workspace->nlat-1);
      }
    };


  public:
    /// Constructor
    Mesh(shtns_cfg w, double radius)
      : _workspace(w)
      , R(radius)
    {}


    /// Return a proxy to the mesh nodes
    NodeRange nodes() const
    {
      return {_workspace, R};
    }

    /// Return a proxy to the mesh cells
    CellRange cells() const
    {
      return {_workspace};
    }

    /// iterator over cells
    CellRange::iterator begin() const
    {
      return {_workspace};
    }

    /// end-iterator over cells
    CellRange::iterator end() const
    {
      return {_workspace, true};
    }

    size_t size() const
    {
      return CellRange(_workspace).size();
    }

    size_t vertices() const
    {
      return NodeRange(_workspace,R).size();
    }

    shtns_cfg workspace() const
    {
      return _workspace;
    }

  private:
    double R;
    shtns_cfg _workspace;
  };

} // end namespace Shtns
