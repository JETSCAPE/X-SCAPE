/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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
//Parton Gun

#ifndef PYTHIATGUN_H
#define PYTHIATGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "TransverseEPS09s.h"
//#include "Pythia8/Pythia.h"

using namespace Jetscape;

class PythiaTGun : public HardProcess {

private:
  double pTHatMin;
  double pTHatMax;
  double ImpParMean;
  double s_1x;
  double s_1y;
  int n_1x;
  int n_1y;
  int BeamId;
  double eCM;
  double vir_factor;
  bool initial_virtuality_pT;
  double softMomentumCutoff;
  bool FSR_on;
  bool softQCD;
  std::vector<std::unique_ptr<Pythia8::Pythia>> pythia_vec;
  void ConfigurePythia(Pythia8::Pythia& py, unsigned int seed);
  int GetPythiaGridIndex(double x, double y) const;
  double  CrossSection;
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<PythiaTGun> reg;

public:
  PythiaTGun();
  ~PythiaTGun();
  void InitTask();
  void ExecuteTask();
  void GetCrossSec();
  // Getters
  double GetpTHatMin() const { return pTHatMin; }
  double GetpTHatMax() const { return pTHatMax; }

  // Cross-section information in mb and event weight.
  double GetSigmaGen() { return CrossSection; }
  double GetSigmaErr() { return 0;}//infoPtr.sigmaErr(); };
  double GetPtHat() { return 0;}//infoPtr.pTHat(); };
  double GetEventWeight() { return 0;}//infoPtr.weight(); };
  //void Test(int ipy);
};

#endif // PYTHIATGUN_H
