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

#ifndef NCOLLHISTO_H
#define NCOLLHISTO_H

#include <random>
#include <string>
#include <vector>

#include "InitialState.h"

using namespace Jetscape;

// Hard-scattering vertices sampled from a histogram of binary-collision
// positions ("x y N_coll", bin centres on a uniform grid), as NCOLLHISTO in the
// JEWEL (2+1)D medium interface (I. Kolbe, github.com/isobelkolbe/jewel-2.4.0-2D):
// a bin with probability proportional to N_coll, then a uniform position
// within the bin.
class NcollHisto : public Jetscape::InitialState {
 public:
  NcollHisto();
  ~NcollHisto();
  void InitTask();
  void ExecuteTask() {}
  void ClearTask() {}
  double GetNcoll() { return ncoll_; }
  void SampleABinaryCollisionPoint(double &t, double &x, double &y, double &z);

 private:
  std::vector<double> x_, y_, cumulative_;
  double dx_ = 0.0, ncoll_ = 0.0;
  std::uniform_real_distribution<double> uniform_{0.0, 1.0};

  static RegisterJetScapeModule<NcollHisto> reg;
};

#endif  // NCOLLHISTO_H
