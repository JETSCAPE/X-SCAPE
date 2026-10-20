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

// JEWEL 2.6.0 (C++ port, external_packages/jewel) as final-state shower.
//
// JEWEL does not evolve in time steps: on the first call for a hard parton the
// complete medium-modified shower of that parton is generated (vacuum
// splittings, elastic and radiative scatterings with LPM interference) in the
// medium given by the JETSCAPE fluid (GetHydroCellSignal, as in LBT). The
// module returns the final shower partons (status 0), the medium recoils
// (status 1) and the thermal momenta of the scattering centres before the
// scattering (status -1, to be subtracted), the status convention of LBT. All
// returned partons carry a JewelUserInfo and are left untouched afterwards.
// JEWEL replaces both MATTER and LBT and cannot be combined with other energy
// loss modules.

#ifndef JEWEL_H
#define JEWEL_H

#include "JetEnergyLossModule.h"

using namespace Jetscape;

// marks partons that have been produced by JEWEL
class JewelUserInfo : public fjcore::PseudoJet::UserInfoBase {
 public:
  JewelUserInfo() = default;
};

class Jewel : public JetEnergyLossModule<Jewel> {
 public:
  Jewel();
  virtual ~Jewel();

  void InitTask();
  void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton> &pIn,
                    vector<Parton> &pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);
  // summary of the JEWEL counters at the end of the run (once per process)
  void FinishTask();

 private:
  double qmax_factor_ = 1.0;  // maximum virtuality = qmax_factor * pT

  // JEWEL has global state: it is initialised once per process
  static bool jewel_initialised_;

  // Allows the registration of the module so that it is available to be used
  // by the Jetscape framework.
  static RegisterJetScapeModule<Jewel> reg;
};

#endif  // JEWEL_H
