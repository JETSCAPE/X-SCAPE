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

#ifndef MEDIUMCONTOURDUMP_H
#define MEDIUMCONTOURDUMP_H

#include "JetEnergyLossModule.h"

using namespace Jetscape;

// Writes the fluid of the current event at mid-rapidity (z = 0) in the contour
// format of the JEWEL (2+1)D medium interface (I. Kolbe, arXiv:2303.14166,
// github.com/isobelkolbe/jewel-2.4.0-2D): one file Tcontour<tau> ("x y T") and
// one file Vcontour<tau> ("x y v_x v_y") per proper time, and a histogram of
// the binary-collision positions of the initial state ("x y N_coll") for her
// NCOLLHISTO. Used as the only Eloss module of a one-event run; the partons
// are left unchanged. The files can be read back by CustomModuleHydroContours.
class MediumContourDump : public JetEnergyLossModule<MediumContourDump> {
 public:
  MediumContourDump();
  virtual ~MediumContourDump();

  void InitTask();
  void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton> &pIn,
                    vector<Parton> &pOut);

 private:
  void Dump();

  std::string dir_;
  double xmax_ = 20.0, dx_ = 0.2;   // grid -xmax..xmax in x and y [fm]
  double tau0_ = -1.0, dtau_ = 0.2, taumax_ = 20.0;  // [fm]; tau0 < 0: fluid start
  long ncoll_samples_ = 200000;     // samples of the binary-collision positions
  bool contours_ = true;            // false: only the N_coll histogram

  static bool done_;
  static RegisterJetScapeModule<MediumContourDump> reg;
};

#endif  // MEDIUMCONTOURDUMP_H
