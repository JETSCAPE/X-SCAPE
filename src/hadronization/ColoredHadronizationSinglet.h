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

#ifndef COLOREDHADRONIZATIONSINGLET_H
#define COLOREDHADRONIZATIONSINGLET_H

#include "HadronizationModule.h"
#include "Pythia8/Pythia.h"
#include "JetScapeSignalManager.h"

using namespace Jetscape;

class ColoredHadronizationSinglet : public HadronizationModule<ColoredHadronizationSinglet> {
public:
  ColoredHadronizationSinglet();
  virtual ~ColoredHadronizationSinglet();

  void InitTask();
  void DoHadronization(vector<vector<shared_ptr<Parton>>> &shower,
                       vector<shared_ptr<Hadron>> &hOut,
                       vector<shared_ptr<Parton>> &pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

  /** @brief Pointer to the InitialState module. */
  std::shared_ptr<InitialState> ini;

private:
  double p_fake;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<ColoredHadronizationSinglet> reg;

protected:
  static Pythia8::Pythia pythia;
};

#endif // COLOREDHADRONIZATIONSINGLET_H
