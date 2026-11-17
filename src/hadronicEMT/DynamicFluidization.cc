/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "DynamicFluidization.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>


namespace Jetscape {

DynamicFluidization::DynamicFluidization(double sigma, double e_threshold,
                                         double cutoff)
    : sigma_(sigma), e_threshold_(e_threshold), cutoff_(cutoff) {}

DynamicFluidization::~DynamicFluidization() = default;

namespace {
struct CellKey {
  int ix, iy, iz;
  bool operator==(const CellKey &o) const {
    return ix == o.ix && iy == o.iy && iz == o.iz;
  }
};
struct CellHash {
  std::size_t operator()(const CellKey &k) const {
    return (static_cast<std::size_t>(k.ix) * 73856093u) ^
           (static_cast<std::size_t>(k.iy) * 19349663u) ^
           (static_cast<std::size_t>(k.iz) * 83492791u);
  }
};
}  // namespace

std::vector<bool> DynamicFluidization::Evaluate(
    const std::vector<Candidate> &hadrons) {
  const std::size_t n = hadrons.size();
  std::vector<bool> fluidize(n, false);
  last_e_.assign(n, 0.);
  if (n == 0) return fluidize;

  const double two_s2 = 2. * sigma_ * sigma_;
  const double norm = 1. / std::pow(2. * M_PI * sigma_ * sigma_, 1.5);
  const double rcut = cutoff_ * sigma_;
  const double rcut2 = rcut * rcut;

  // cell list
  auto key_of = [&](const Candidate &h) {
    return CellKey{static_cast<int>(std::floor(h.x / rcut)),
                   static_cast<int>(std::floor(h.y / rcut)),
                   static_cast<int>(std::floor(h.z / rcut))};
  };
  std::unordered_map<CellKey, std::vector<std::size_t>, CellHash> cells;
  for (std::size_t i = 0; i < n; i++) {
    if (hadrons[i].mass <= 0. || hadrons[i].E <= 0.) continue;
    cells[key_of(hadrons[i])].push_back(i);
  }

  for (std::size_t j = 0; j < n; j++) {
    const Candidate &hj = hadrons[j];
    if (!hj.eligible || hj.mass <= 0.) continue;
    std::array<std::array<double, 4>, 4> T = {};
    const CellKey kj = key_of(hj);
    for (int dx = -1; dx <= 1; dx++)
      for (int dy = -1; dy <= 1; dy++)
        for (int dz = -1; dz <= 1; dz++) {
          auto it = cells.find(CellKey{kj.ix + dx, kj.iy + dy, kj.iz + dz});
          if (it == cells.end()) continue;
          for (const std::size_t i : it->second) {
            const Candidate &hi = hadrons[i];
            const double rx = hj.x - hi.x, ry = hj.y - hi.y, rz = hj.z - hi.z;
            const double r2 = rx * rx + ry * ry + rz * rz;
            if (r2 > rcut2) continue;
            const double gamma = hi.E / hi.mass;
            const double rdotu =
                (rx * hi.px + ry * hi.py + rz * hi.pz) / hi.mass;
            const double K =
                gamma * norm * std::exp(-(r2 + rdotu * rdotu) / two_s2);
            const double p[4] = {hi.E, hi.px, hi.py, hi.pz};
            const double w = K / hi.E;
            for (int a = 0; a < 4; a++)
              for (int b = a; b < 4; b++) T[a][b] += w * p[a] * p[b];
          }
        }
    for (int a = 0; a < 4; a++)
      for (int b = 0; b < a; b++) T[a][b] = T[b][a];

    double e = 0.;
    std::array<double, 4> u = {1., 0., 0., 0.};
    if (T[0][0] > 0. && !LandauMatch(T, e, u)) e = 0.;
    last_e_[j] = e;
    fluidize[j] = (e - hj.mass * norm >= e_threshold_);
  }
  return fluidize;
}

bool DynamicFluidization::LandauMatch(
    const std::array<std::array<double, 4>, 4> &T, double &e,
    std::array<double, 4> &u) {
  // M = T^mu_nu = T^{mu a} g_{a nu}; its eigenvalues are e (timelike
  // eigenvector u) and -P_i (spacelike), with |P_i| <= e for matter.  Power
  // iteration on M + c 1, c = T^{00} >= e: e + c dominates |c - P_i|.
  double M[4][4];
  for (int a = 0; a < 4; a++)
    for (int b = 0; b < 4; b++) M[a][b] = (b == 0 ? T[a][b] : -T[a][b]);
  const double c = T[0][0];
  double v[4] = {1., 0., 0., 0.};
  double lambda = 0.;
  for (int it = 0; it < 500; it++) {
    double w[4];
    for (int a = 0; a < 4; a++) {
      w[a] = c * v[a];
      for (int b = 0; b < 4; b++) w[a] += M[a][b] * v[b];
    }
    // normalize to u.u = 1 (timelike) if possible, else by the time component
    const double norm2 = w[0] * w[0] - w[1] * w[1] - w[2] * w[2] - w[3] * w[3];
    if (!(norm2 > 0.) || !(w[0] > 0.)) return false;
    const double nrm = std::sqrt(norm2);
    double change = 0.;
    for (int a = 0; a < 4; a++) {
      const double vn = w[a] / nrm;
      change = std::max(change, std::abs(vn - v[a]));
      v[a] = vn;
    }
    lambda = nrm - c;   // |(M + c) v| for the normalized eigenvector
    if (change < 1e-12) break;
  }
  // e = u_mu T^{mu nu} u_nu (exact for the eigenvector, robust otherwise)
  double Tu[4] = {0., 0., 0., 0.};
  const double ul[4] = {v[0], -v[1], -v[2], -v[3]};
  for (int a = 0; a < 4; a++)
    for (int b = 0; b < 4; b++) Tu[a] += T[a][b] * ul[b];
  e = 0.;
  for (int a = 0; a < 4; a++) e += ul[a] * Tu[a];
  (void)lambda;
  for (int a = 0; a < 4; a++) u[a] = v[a];
  return e > 0.;
}

}  // namespace Jetscape
