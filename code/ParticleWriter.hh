#pragma once

#include <cmath>
#include <iostream>
#include <fstream>

#include "AlgLibTools.h"
#include "tools.hpp"

#include "kdtree_nanoflann.h"

namespace Shtns
{

  template <class R>
  std::array<R,3> project(std::array<R,3> const& p, std::array<R, 3> normal)
  {
    R s = dot(normal, p);
    normal *= s;
    return p - normal;
  }


  template <class R, class C>
  class ParticleWriter
  {
  public:
    using WorldVector = std::array<R, 3>;

    ParticleWriter(shtns_cfg workspace, R radius = 1.0)
      : workspace(workspace)
      , radius(radius)
      , old_maxima(1)
      , old_tau(1)
    {
      nphi = workspace->nphi;
      nlat = workspace->nlat;

      init_mesh();
    }


    void init_mesh()
    {
      points.reserve(nphi * nlat);

      WorldVector X;
      for (size_t i = 0; i < nphi; ++i) {
        for (size_t j = 0; j < nlat; ++j) {
          get_coord(workspace, i, j, radius, X);
          points.push_back(X);
        }
      }

      mesh.reset(new KD_Tree(points));
    }


    std::tuple<double, double, double> // min, max, threshold
    get_threshold(ScalarValues<R,C> const& psi, double rel_threshold = 0.3)
    {
      double min_psi =  1.0e50;
      double max_psi = -1.0e50;
      for (size_t i = 0; i < nphi; ++i) {
        for (size_t j = 0; j < nlat; ++j) {
          size_t idx = get_idx(nlat, i, j);
          min_psi = std::min( min_psi, psi[idx] );
          max_psi = std::max( max_psi, psi[idx] );
        }
      }

      double threshold = max_psi - rel_threshold * (max_psi - min_psi);

      return std::make_tuple(min_psi, max_psi, threshold);
    }


    void find_maxima(ScalarValues<R,C> const& psi, double threshold, std::vector<WorldVector>& maxima)
    {
      WorldVector X, X0;
      std::vector<std::pair<size_t, size_t>> points_idx;
      for (size_t i = 0; i < nphi; ++i) {
        for (size_t j = 1; j < nlat-1; ++j) {

          size_t idx = get_idx(nlat, i, j); // center
          if (psi[idx] > threshold) {
            size_t idx_r = get_idx(nlat, (i+1) % nphi, j);
            size_t idx_l = get_idx(nlat, (i-1+nphi) % nphi, j);
            size_t idx_t = get_idx(nlat, i, j+1);
            size_t idx_b = get_idx(nlat, i, j-1);

            if (psi[idx] > psi[idx_r] && psi[idx] > psi[idx_l] &&
                psi[idx] > psi[idx_t] && psi[idx] > psi[idx_b])
            {
              get_coord(workspace, i, j, radius, X0);

              bool already_found = false;

              auto found_it = points_idx.begin();
              for (auto it = points_idx.begin(); it != points_idx.end(); ++it) {
                get_coord(workspace, it->first, it->second, radius, X);
                R dist = std::sqrt(sqr(X[0] - X0[0]) + sqr(X[1] - X0[1]) + sqr(X[2] - X0[2]));
                if (dist < 0.5*d) {
                  already_found = true;
                  found_it = it;
                  break;
                }
              }

              if (already_found) {
                size_t idx_k = get_idx(nlat, found_it->first, found_it->second);
                if (psi[idx] > psi[idx_k])
                  *found_it = std::make_pair(i,j);
              } else {
                points_idx.push_back(std::make_pair(i,j));
              }
            }
          }
        }
      }

      maxima.reserve(points_idx.size());
      for (auto const& p : points_idx) {
        get_coord(workspace, p.first, p.second, radius, X);
        maxima.push_back(X);
      }
    }


    void refine_maxima(ScalarValues<R,C> const& psi, std::vector<WorldVector>& maxima, double min_psi, double max_psi)
    {
      std::vector<R> values;
      values.reserve(nphi * nlat);

      for (size_t i = 0; i < nphi; ++i) {
        for (size_t j = 0; j < nlat; ++j) {
          size_t idx = get_idx(nlat, i, j);
          values.push_back(psi[idx]);
        }
      }

      const size_t npoints = 10;
      std::vector<size_t> out_indices(npoints);
      std::vector<double> out_distances_sq(npoints);

      for (auto& p : maxima) {
        // get 10 nearest vertices
        mesh->query(p.data(), npoints, out_indices.data(), out_distances_sq.data());

        std::vector<WorldVector> p_; p_.reserve(npoints);
        std::vector<double> v_; v_.reserve(npoints);
        for (auto i : out_indices) {
          p_.push_back(points[i]);
          v_.push_back(values[i]);
        }

        pfc_tools::nonlinearPtcFit(p_, v_, p, min_psi, max_psi);
        p *= radius / norm(p);
      }
    }


    // calculate mean direction
    void get_direction(std::shared_ptr<KD_Tree> maxima_tree,
                      VectorValues<R,C> const& P, ScalarValues<R,C> const& psi,
                      std::vector<WorldVector> const& maxima, std::vector<WorldVector>& direction,
                      double min_psi, double max_psi)
    {
      direction.resize(maxima.size());
      std::vector<size_t> num_direction(maxima.size(),0);

      WorldVector P0, P1;
      size_t l = 0;
      for (size_t i = 0; i < nphi; ++i) {
        for (size_t j = 0; j < nlat; ++j, ++l) {

          auto const& X = points[l];

          std::vector<std::pair<size_t,R> > indicesDists; indicesDists.reserve(1);

          // find closest maxima
          maxima_tree->index->radiusSearch(
            X.data(),
            sqr(0.4*d),
            indicesDists,
            nanoflann::SearchParams()
          );

          if (indicesDists.size() > 0) {
            const size_t idx = get_idx(nlat, i, j);
            sphere_to_cartesian(workspace, P.values_t(), P.values_p(), i, j, P0);

            for (auto const& indexDist : indicesDists) {
              const size_t idx_store = indexDist.first;

              P1 = P0 * ((psi[idx] - min_psi) / (max_psi - min_psi));
              direction[idx_store] += P1;
              num_direction[idx_store]++;
            }
          }
        }
      }

      for (size_t i = 0; i < direction.size(); ++i) {
        if (num_direction[i] > 0) {
          direction[i] *= 1.0/R(num_direction[i]);
        } else {
          std::cout << "no directions added!\n";
        }
        auto n = maxima[i]; n *= 1.0/std::max(1.e-5, norm(maxima[i]));
        auto q = project(direction[i], n);
        direction[i] = q * (1.0/(std::max(1.e-5,norm(q))));
      }
    }



    // calculate order parameter
    void get_order_parameter(std::shared_ptr<KD_Tree> maxima_tree,
                            std::vector<WorldVector> const& maxima, std::vector<WorldVector> const& direction,
                            std::vector<double>& polar, std::vector<double>& vop)
    {
      polar.resize(maxima.size(), 0.0);
      vop.resize(maxima.size(), 0.0);
      for (size_t i = 0; i < maxima.size(); ++i) {

        std::vector<std::pair<size_t,R> > indicesDists; indicesDists.reserve(14);
        // find closest maxima
        maxima_tree->index->radiusSearch(
          maxima[i].data(),
          sqr(2.5*d),
          indicesDists,
          nanoflann::SearchParams()
        );

  //       maxima_tree->query(maxima[i].data(), npoints, out_indices.data(), out_distances_sq.data());

        auto n = maxima[i]; n *= 1.0/std::max(1.e-5, norm(n));

        double phi = std::atan2(maxima[i][1], maxima[i][0]);

        WorldVector ex; ex[0] = 1.0;
        WorldVector ey; ey[1] = 1.0;
        WorldVector t_phi = -std::sin(phi) * ex + std::cos(phi) * ey;
        t_phi *= 1.0/std::max(1.e-10, norm(t_phi));

        // TODO: implement VOP-Formel
        double weight = 0;

        double nominator = 0.0;
        double denominator = 0.0;
        for (auto const& neighbor : indicesDists) {
          if (neighbor.first == i)
            continue;

          auto q = project(direction[neighbor.first], n); q *= 1.0/std::max(1.e-5, norm(q));

          double x = dot(direction[i], q);

          polar[i] += x / std::sqrt(neighbor.second);
          nominator += dot(direction[i], t_phi);
          denominator += std::sqrt(dot(direction[i],direction[i]));
          weight += 1.0/std::sqrt(neighbor.second);
        }

        vop[i] = 1.0/(1.0 - 2.0/M_PI) * (nominator/denominator - 2.0/M_PI);

        if (weight > 0) {
          polar[i] /= R(weight);
        }
      }
    }


    void get_velocity(std::shared_ptr<KD_Tree> maxima_tree,
                      std::vector<WorldVector> const& maxima,
                      std::vector< std::vector<WorldVector> > const& old_maxima,
                      double tau, std::vector<double> const& old_tau, double v0,
                      std::vector<WorldVector>& velocity)
    {

      size_t idx = (last_idx-1)%old_maxima.size();

      velocity.resize(maxima.size());
      WorldVector zero{0.0, 0.0, 0.0};
      for (auto& v : velocity)
        v = zero;

      std::vector<char> found(maxima.size(), 0);
      for (auto const& m_old : old_maxima[idx]) {

        size_t index = 0;
        double distance = -1.0;
        maxima_tree->query(m_old.data(), 1, &index, &distance);
        if (distance < sqr(d/2.0)) {
          velocity[index] = (maxima[index] - m_old) * (1.0/tau);
          if (norm(velocity[index]) > 2*v0)
            velocity[index] *= 0.0;
          else
            found[index] = 1;
        }
      }

      const size_t npoints = 7;
      std::vector<size_t> out_indices(npoints);
      std::vector<double> out_distances(npoints);
      for (size_t i = 0; i < found.size(); ++i) {
        if (found[i] == 0) {
          maxima_tree->query(maxima[i].data(), npoints, out_indices.data(), out_distances.data());

          WorldVector vel{0.0, 0.0, 0.0};
          double w = 0.0;
          for (size_t j = 0; j < npoints; ++j) {
            size_t const index = out_indices[j];
            if (i != index && found[index] != 0) {
              vel += velocity[index] * (1.0/std::sqrt(out_distances[j]));
              w += 1.0/std::sqrt(out_distances[j]);
            }
          }
          if (w > 1.e-10)
            vel *= 1.0/w;
          velocity[i] = vel;
        }
      }

    }


    void write(std::string filename,
              VectorValues<R,C> const& P, ScalarValues<R,C> const& psi,
              double tau = 1.0,
              double v0 = 0.5)
    {
      std::vector<WorldVector> maxima;
      std::vector<WorldVector> direction;
      std::vector<WorldVector> velocity;
      std::vector<double> polar;
      std::vector<double> vop;

      double min_psi, max_psi, threshold;
      std::tie(min_psi, max_psi, threshold) = get_threshold(psi);

      find_maxima(psi, threshold, maxima);
      refine_maxima(psi, maxima, min_psi, max_psi);

      auto tree = std::make_shared<KD_Tree>(maxima);

      get_direction(tree, P, psi, maxima, direction, min_psi, max_psi);
      get_order_parameter(tree, maxima, direction, polar, vop);

      velocity.resize(maxima.size());
      if (last_idx >= old_maxima.size())
        get_velocity(tree, maxima, old_maxima, tau, old_tau, v0, velocity);

      std::ofstream out(filename, std::ios_base::out);

      //out << points.size() << "\n";
      out << "x,y,z,d0,d1,d2,v0,v1,v2,polar,vop\n";
      for (size_t i = 0; i < maxima.size(); ++i) {
        for (size_t j = 0; j < 3; ++j)
          out << maxima[i][j] << ',';

        for (size_t j = 0; j < 3; ++j)
          out << direction[i][j] << ',';

        for (size_t j = 0; j < 3; ++j)
          out << velocity[i][j] << ',';

        out << polar[i] << ',' << vop[i] << '\n';
      }
      out.close();


      size_t new_idx = last_idx%old_maxima.size();
      std::swap(old_maxima[new_idx], maxima);
      old_tau[new_idx] = tau;

      last_idx++;
    }


  private:
    std::string filename;
    shtns_cfg workspace;
    double radius;

    size_t nphi;
    size_t nlat;

    std::vector<WorldVector> points;
    std::unique_ptr<KD_Tree> mesh;

    std::vector< std::vector<WorldVector> > old_maxima;
    std::vector< double > old_tau;
    size_t last_idx = 0;

    const double d = 4.0*M_PI/std::sqrt(3.0);
  };

} // end namespace Shtns
