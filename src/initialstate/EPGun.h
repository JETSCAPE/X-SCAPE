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

#ifndef EPGUN_H
#define EPGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class EPGun : public HardProcess, public Pythia8::Pythia {

private:
  double eProton   = 920.;
  double eElectron = 27.5;
  double Q2min     = 25.;
  double vir_factor;
  double softMomentumCutoff;
  bool FSR_on = false;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<EPGun> reg;

public:
  /** standard ctor
      @param xmlDir: Note that the environment variable PYTHIA8DATA takes precedence! So don't use it.
      @param printBanner: Suppress starting blurb. Should be set to true in production, credit where it's due
  */
  EPGun(string xmlDir = "DONTUSETHIS", bool printBanner = false)
      : Pythia8::Pythia(xmlDir, printBanner), HardProcess() {
    SetId("UninitializedEPGun");
  }

  ~EPGun();

  void InitTask();
  void ExecuteTask();

  // Cross-section information in mb and event weight.
  double GetSigmaGen() { return info.sigmaGen(); };
  double GetSigmaErr() { return info.sigmaErr(); };
  double GetEventWeight() { return info.weight(); };

  std::shared_ptr<Hadron> PythiaToJSHadron(Pythia8::Particle &particle){

    std::shared_ptr<Hadron> jshadron = std::make_shared<Hadron>
      (Hadron(1,particle.id(),801,particle.pT(),particle.eta(),particle.phi(),particle.e(),0));
    return jshadron;
  }
};

#endif // EPGun_H
