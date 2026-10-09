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

// Holographic (Einstein-Maxwell-Dilaton) jet quenching parameter at mu_B = 0.
//
// The EMD calculation gives Q(T) = qhat / (sqrt(lambda_t) T^3), fitted by
//   Q(T) = A/2 [1 - tanh((T-T0)/dT)]
//        + 1/2 [1 + tanh((T-T0)/dT)] [B + C exp(-(T-T0)/tau)]
// with T in MeV. The fit reproduces the tabulated EMD result (20-800 MeV) to
// 0.1% on average and 2.5% at worst. Q(T) peaks at T ~ 213 MeV.
//
// The holographic qhat is extracted from an adjoint (light-like) Wilson loop,
// i.e. it is the gluon qhat; quarks are obtained by scaling with C_F/C_A.

#ifndef EMDQHAT_H
#define EMDQHAT_H

#include <cmath>

namespace EMDQhat {

// qhat / (sqrt(lambda_t) T^3); T in GeV.
inline double QhatOverSqrtLambdaT3(double T) {
  const double A = 3.98564, B = 8.10975, C = 4.10263;
  const double T0 = 133.6481, dT = 50.3029, tau = 254.864;  // MeV
  const double TMeV = 1000. * T;
  const double th = std::tanh((TMeV - T0) / dT);
  return 0.5 * A * (1. - th) +
         0.5 * (1. + th) * (B + C * std::exp(-(TMeV - T0) / tau));
}

// Gluon qhat in GeV^3 for temperature T [GeV] and 't Hooft coupling lambda_t.
// conformal_Tref > 0 freezes Q(T) at Q(conformal_Tref): a conformal qhat ~ T^3
// with the same value at T = conformal_Tref (control without the EMD shape).
inline double GluonQhat(double T, double lambda_t, double conformal_Tref = 0.) {
  if (T <= 0. || lambda_t <= 0.)
    return 0.;
  double Q = QhatOverSqrtLambdaT3(conformal_Tref > 0. ? conformal_Tref : T);
  return Q * std::sqrt(lambda_t) * T * T * T;
}

}  // namespace EMDQhat

#endif  // EMDQHAT_H
