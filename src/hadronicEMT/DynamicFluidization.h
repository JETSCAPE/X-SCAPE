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

#ifndef DYNAMICFLUIDIZATION_H
#define DYNAMICFLUIDIZATION_H

#include <array>
#include <memory>
#include <vector>

namespace Jetscape {

/**
 * The fluidization criterion of dynamic (concurrent) hydro + transport
 * running, after arXiv:2507.19389 (SMASH + vHLLE, "dynamic fluidization").
 *
 * At one time t, for hadrons at the same t (a transport step): the local
 * energy-momentum tensor at hadron j,
 *
 *   T^{mu nu}(x_j) = sum_i p_i^mu p_i^nu / p_i^0  K(x_j - x_i, p_i/m_i),
 *   K(r, u)        = gamma / (2 pi sigma^2)^{3/2}
 *                    exp(-(|r|^2 + (r.u)^2) / (2 sigma^2)),
 *
 * Landau-matched (T^mu_nu u^nu = e u^mu), and hadron j fluidizes when it
 * is eligible and e - m_j/(2 pi sigma^2)^{3/2} >= e_th; the subtracted term
 * is its own contribution in its rest frame, so a heavy hadron does not
 * fluidize by itself.  Eligibility (string formation time, unscattered
 * nucleons, ...) is decided by the caller, which knows the transport's
 * history.  Every hadron contributes to T^{mu nu}.
 *
 * Neighbours are found with a cell list of side cutoff*sigma, so the cost
 * is ~ N x (hadrons within the cutoff), not N^2.
 */
class DynamicFluidization {
 public:
  struct Candidate {
    double t, x, y, z;           ///< position [fm]
    double E, px, py, pz;        ///< momentum [GeV]
    double mass;                 ///< [GeV]
    bool eligible;               ///< may fluidize (caller's history checks)
  };

  /// sigma [fm], e_threshold [GeV/fm^3], cutoff in units of sigma
  DynamicFluidization(double sigma, double e_threshold, double cutoff = 4.);
  ~DynamicFluidization();

  /// one flag per candidate: true = fluidize now
  std::vector<bool> Evaluate(const std::vector<Candidate> &hadrons);

  /// the rest-frame energy density at each candidate of the last Evaluate,
  /// before subtracting the own contribution [GeV/fm^3] (diagnostics)
  const std::vector<double> &LastEnergyDensity() const { return last_e_; }

  /// Landau matching of a T^{mu nu} (upper indices, metric +---): the
  /// timelike eigenvector of T^mu_nu and its eigenvalue; false if none found
  static bool LandauMatch(const std::array<std::array<double, 4>, 4> &T,
                          double &e, std::array<double, 4> &u);

  double sigma() const { return sigma_; }
  double e_threshold() const { return e_threshold_; }

 private:
  double sigma_;
  double e_threshold_;
  double cutoff_;
  std::vector<double> last_e_;
};

}  // namespace Jetscape

#endif  // DYNAMICFLUIDIZATION_H
