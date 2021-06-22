#pragma once

#include "nanoflann.hpp"

#include <cstdlib>
#include <iostream>

using namespace nanoflann;


/** A simple vector-of-vectors adaptor for nanoflann, without duplicating the storage.
  *  The i'th vector represents a point in the state space.
  *
  *  \tparam DIM If set to >0, it specifies a compile-time fixed dimensionality for the points in the data set, allowing more compiler optimizations.
  *  \tparam value_type The type of the point coordinates (typically, double or float).
  *  \tparam Distance The distance metric to use: nanoflann::metric_L1, nanoflann::metric_L2, nanoflann::metric_L2_Simple, etc.
  *  \tparam index_type The type for indices in the KD-tree index (typically, size_t of int)
  */
template <typename value_type,
          int DIM,
          class Distance = nanoflann::metric_L2_Simple,
          typename index_type = size_t>
struct KDTreeVectorOfArrayAdaptor
{
  using PointType = std::array<value_type, DIM>;
  using VectorOfVectorsType = std::vector<PointType>;


  using self_t = KDTreeVectorOfArrayAdaptor;

  using metric_t = typename Distance::template traits<value_type, self_t>::distance_t;
  using index_t = KDTreeSingleIndexAdaptor< metric_t, self_t, DIM, index_type>;

  index_t* index; //! The kd-tree index for the user to call its methods as usual with any other FLANN index.

  /// Constructor: takes a const ref to the vector of vectors object with the data points
  KDTreeVectorOfArrayAdaptor(const VectorOfVectorsType &mat, const int leaf_max_size = 10)
    : m_data(mat)
  {
    assert( !mat.empty() );

    index = new index_t( DIM, *this /* adaptor */, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size) );
    index->buildIndex();
  }

  ~KDTreeVectorOfArrayAdaptor()
  {
    delete index;
  }

  const VectorOfVectorsType &m_data;

  void reinit(const int leaf_max_size = 10)
  {
    delete index;
    index = new index_t( DIM, *this /* adaptor */, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size) );
    index->buildIndex();
  }

  /** Query for the \a num_closest closest points to a given point (entered as query_point[0:dim-1]).
    *  Note that this is a short-cut method for index->findNeighbors().
    *  The user can also call index->... methods as desired.
    * \note nChecks_IGNORED is ignored but kept for compatibility with the original FLANN interface.
    */
  void query(const value_type *query_point, const size_t num_closest, index_type *out_indices, value_type *out_distances_sq, const int /*nChecks_IGNORED*/ = 10) const
  {
    nanoflann::KNNResultSet<double,index_type> resultSet(num_closest);
    resultSet.init(out_indices, out_distances_sq);
    index->findNeighbors(resultSet, query_point, nanoflann::SearchParams());
  }

  /** @name Interface expected by KDTreeSingleIndexAdaptor
    * @{ */

  const self_t& derived() const { return *this; }
  self_t&       derived()       { return *this; }

  // Must return the number of data points
  size_t kdtree_get_point_count() const {
    return m_data.size();
  }

  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
  value_type kdtree_distance(const value_type *p1, const size_t idx_p2, size_t size) const
  {
    value_type s=0;
    for (int i = 0; i < DIM; ++i)
    {
      const value_type d = p1[i] - m_data[idx_p2][i];
      s += d*d;
    }
    return s;
  }

  // Returns the dim'th component of the idx'th point in the class:
  value_type kdtree_get_pt(const size_t idx, int dim) const {
    return m_data[idx][dim];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const {
          return false;
  }

  /** @} */

}; // end of KDTreeVectorOfVectorsAdaptor

using KD_Tree = KDTreeVectorOfArrayAdaptor<double, 3>;
