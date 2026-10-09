/*************************************************************************************
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
 *
 *************************************************************************************/

// The HTL jet quenching parameter used by MATTER (QhatParametrizationType 1,
// and 5 below Q0 where the virtuality factor is 1), so that MyEloss can run
// with exactly the qhat(T, E) of the JETSCAPE MATTER+LBT tune:
//   qhat = (C_A 50.4864 / pi) alpha_s^run(2ET) alpha_s T^3 ln(2ET / m_D^2),
//   m_D^2 = alpha_s 4 pi T^2 (6 + n_f) / 6,   2ET >= 1 GeV^2,
// with alpha_s^run from Matter::RunningAlphaS (Lambda_QCD fixed by alpha_s at
// 1 GeV^2). Gluon qhat in GeV^3; T and E in GeV.

#ifndef HTLQHAT_H
#define HTLQHAT_H

#include <cmath>

namespace HTLQhat {

inline double RunningAlphaS(double muSquare, double alphas, int nf = 3) {
  if (muSquare < 1.0)
    return alphas;
  const double lambda2 = std::exp(-12.0 * M_PI / ((33 - 2 * nf) * alphas));
  return 12.0 * M_PI / ((33.0 - 2.0 * nf) * std::log(muSquare / lambda2));
}

inline double GluonQhat(double T, double E, double alphas, int nf = 3) {
  if (T <= 0.)
    return 0.;
  const double CA = 3.0;
  double scale = std::max(2.0 * E * T, 1.0);
  double mD2 = alphas * 4.0 * M_PI * T * T * (6.0 + nf) / 6.0;
  double q = (CA * 50.4864 / M_PI) * RunningAlphaS(scale, alphas, nf) * alphas *
             T * T * T * std::log(scale / mD2);
  return q > 0. ? q : 0.;
}

}  // namespace HTLQhat

#endif  // HTLQHAT_H
