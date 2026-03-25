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

// Create a pythia collision at a specified point and return the two inital hard partons

#ifndef PYTHIAISRGUN_H
#define PYTHIAISRGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;
using std::uniform_real_distribution;

class PythiaIsrGun : public HardProcess, public Pythia8::Pythia {

private:
  double pTHatMin;
  double pTHatMax;
  double eCM;
  double cross_section; //User provided cross section for event weighting
  bool FSR_on;
  bool multi_scatter;
  int proj_A;
  int targ_A;
  //First hard collision information when multiple hard collisions
  double first_sigmaGen = -1.0;
  double first_sigmaErr = -1.0;
  double first_ptHat = -1.0;
  double first_weight = -1.0;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<PythiaIsrGun> reg;

public:
  int NumberOfHardPartons;
  /** standard ctor
      @param xmlDir: Note that the environment variable PYTHIA8DATA takes precedence! So don't use it.
      @param printBanner: Suppress starting blurb. Should be set to true in production, credit where it's due
  */
  PythiaIsrGun(string xmlDir = "DONTUSETHIS", bool printBanner = false)
      : Pythia8::Pythia(xmlDir, printBanner), HardProcess() {
    SetId("UninitializedPythiaIsrGun");
  }

  ~PythiaIsrGun();

  void InitTask();
  void ExecuteTask();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

  // Getters
  double GetpTHatMin() const { return pTHatMin; }
  double GetpTHatMax() const { return pTHatMax; }

  // Cross-section information in mb and event weight.
  double GetSigmaGen() { return first_sigmaGen; };
  double GetSigmaErr() { return first_sigmaErr; };
  double GetPtHat() { return first_ptHat; };
  double GetEventWeight() { return first_weight; };

protected:
  uniform_real_distribution<double> ZeroOneDistribution;
  // uniform_int_distribution<int> RandIntDistribution;
};

#endif // PYTHIAISRGUN_H
