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
// Parton Gun

#ifndef PGUN_H
#define PGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

#include <utility>
#include <vector>

using namespace Jetscape;

//! Set when PGun understands <rapidity>, <phi> and <points>.
#define PGUN_HAS_POINTS 1

/** Parton gun: one parton per event from the origin (x = 0, t = 0).

    <pT>, <parID>: as before.  Optional, defaults keep the old gun:
    <rapidity> (0) and <phi> (0) of the parton.
    <points>: "pT y pT y ...", several (pT, rapidity) points in one run;
    event i uses point i mod K (as PythiaGun's <pTHatBins>), and with
    <setReuseHydro> nReuseHydro must be a multiple of K so that every hydro
    event gets the same events per point.  <points> overrides <pT> and
    <rapidity>.
*/
class PGun : public HardProcess {
  static Pythia8::Pythia InternalHelperPythia;

 private:
  double fixed_pT;
  double parID;
  double rapidity_ = 0.0;
  double phi_ = 0.0;
  //! (pT, rapidity) per point; one entry without <points>
  std::vector<std::pair<double, double>> points_;
  int activePoint_ = 0;

  // Allows the registration of the module so that it is available to be used by
  // the Jetscape framework.
  static RegisterJetScapeModule<PGun> reg;

 public:
  PGun();
  ~PGun();

  void InitTask();
  void ExecuteTask();

  int GetNPoints() const { return static_cast<int>(points_.size()); }
  double GetPointPt(int k) const { return points_.at(k).first; }
  double GetPointRapidity(int k) const { return points_.at(k).second; }
  //! the current event's point (i mod K)
  int GetActivePoint() const { return activePoint_; }
  double GetPhi() const { return phi_; }
  double GetParID() const { return parID; }
};

#endif  // PGUN_H
