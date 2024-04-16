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
// -----------------------------------------------------------------------------
// This is a wrapper for SMASH hadronic afterburner with the JETSCAPE framework
// for the collider modus of SMASH to generate initial conditions
// -----------------------------------------------------------------------------

#ifndef SMASHINITIALSTATEWRAPPER_H
#define SMASHINITIALSTATEWRAPPER_H

#include "smash/configuration.h"
#include "smash/experiment.h"
#include "smash/collidermodus.h"

#include "Transport.h"
#include "JetScapeWriter.h"

using namespace Jetscape;

class SmashInitialConditionWrapper : public Transport {
private:

  double end_time_ = -1.0;
  shared_ptr<smash::Experiment<smash::ColliderModus>> smash_collider_experiment_;

  /// Convert Jetscape (JS) hadron list to smash particle list
  smash::ParticleList get_smash_plist_from_JS_hadrons(const std::vector<shared_ptr<Hadron>>& JS_hadrons);

  std::vector<shared_ptr<Hadron>> TestHadronList();

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<SmashInitialConditionWrapper> reg;

public:
  /// Fill the provided Jetscape (JS) hadron list from the SMASH particles
  void fill_JS_hadrons_from_smash_particles(const smash::Particles &smash_particles,
                                            std::vector<shared_ptr<Hadron>> &JS_hadrons);
  SmashInitialConditionWrapper();

  void InitTask();
  void ExecuteTask();
  void WriteTask(weak_ptr<JetScapeWriter> w);

  void InitPerEvent() override;
  void CalculateTimeTask() override;
  void FinishPerEvent() override;

  std::vector<Hadron> GetCurrentHadronList() const override;

  virtual any GetHistory() {return GetCurrentHadronList();}

  void reset_event_numbering() { event_number_ = 0; }
  int current_event_number() {return event_number_;}

  std::vector<std::vector<shared_ptr<Hadron>>> jetscape_hadrons_;
  int event_number_ = 0;
};

#endif // SMASHINITIALSTATEWRAPPER_H