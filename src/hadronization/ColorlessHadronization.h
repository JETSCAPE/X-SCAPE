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

#ifndef COLORLESSHADRONIZATION_H
#define COLORLESSHADRONIZATION_H

#include "HadronizationModule.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class ColorlessHadronization
    : public HadronizationModule<ColorlessHadronization> {
 public:
  ColorlessHadronization();
  virtual ~ColorlessHadronization();

  void InitTask();
  void DoHadronization(vector<vector<shared_ptr<Parton>>> &shower,
                       vector<shared_ptr<Hadron>> &hOut,
                       vector<shared_ptr<Parton>> &pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

  /**
   * Seed the next DoHadronization() call with `seed` (one-shot): Pythia's
   * generator and the module's own one (remnant direction) are both reseeded
   * from it. This is how stored final partons are re-hadronized exactly.
   */
  void SetNextRandomSeed(unsigned int seed) {
    next_seed_ = seed;
    has_next_seed_ = true;
  }

  /// The seed the last DoHadronization() used; 0 if it was not reseeded.
  unsigned int GetLastRandomSeed() const { return last_seed_; }

 private:
  double p_fake;
  bool take_recoil;
  double Lambda_QCD;

  // <JetHadronization><reseed_per_event> 1: every event draws a seed from the
  // module's generator and reseeds with it, so each event can be reproduced
  // on its own (GetLastRandomSeed). 0 (default): Pythia runs on unbroken.
  bool reseed_per_event_ = false;
  bool has_next_seed_ = false;
  unsigned int next_seed_ = 0;
  unsigned int last_seed_ = 0;
  void Reseed(unsigned int seed);

  // Allows the registration of the module so that it is available to be used by
  // the Jetscape framework.
  static RegisterJetScapeModule<ColorlessHadronization> reg;

 protected:
  static Pythia8::Pythia pythia;
  std::uniform_real_distribution<double> ZeroOneDistribution;
};

#endif  // COLORLESSHADRONIZATION_H
