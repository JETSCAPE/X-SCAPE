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
// This is a general basic class for hadronic afterburner

#include "./Afterburner.h"
#include "./JetScapeSignalManager.h"

using namespace std;

namespace Jetscape {
void Afterburner::Init() {
  // Makes sure that XML file with options and parameters is loaded
  JetScapeModuleBase::Init();
  JSINFO << "Initializing Afterburner : " << GetId() << " ...";

  // Get the pointer to soft sampler
  soft_particlization_sampler_ =
      JetScapeSignalManager::Instance()->GetSoftParticlizationPointer().lock();
  if (!soft_particlization_sampler_) {
    JSWARN << "No soft particlization module found. It is necessary to provide"
           << " hadrons to afterburner.";
    // exit(1);
  }

  if (GetXMLElementInt({"Afterburner", "include_fragmentation_hadrons"})) {
    // Get the pointer to the hard sampler
    hard_particlization_module_ = JetScapeSignalManager::Instance()
                                      ->GetHadronizationManagerPointer()
                                      .lock();
  }
  InitTask();
}

// TODO(stdnmr) Make intermediate class more useful for other Afterburner
// e.g. * add afterburner function to get hadrons ftom soft particilizing
// TODO(stdnmr) Add also functionality to get fragmentation hadrons

void Afterburner::Exec() {
  VERBOSE(2) << "Afterburner running: " << GetId() << " ...";
  ExecuteTask();
}

void Afterburner::CalculateTime() {
  VERBOSE(2) << "Afterburner running for time: " << GetId() << " ...";
  CalculateTimeTask();
}

std::vector<std::vector<std::shared_ptr<Hadron>>>
Afterburner::GetSoftParticlizationHadrons() {
  return soft_particlization_sampler_->Hadron_list_;
}

std::vector<shared_ptr<Hadron>> Afterburner::GetFragmentationHadrons() {
  JSINFO << "Get fragmentation hadrons in Afterburner";
  std::vector<shared_ptr<Hadron>> h_list;
  hard_particlization_module_->GetHadrons(h_list);
  JSINFO << "Got " << h_list.size()
         << " fragmentation hadrons from HadronizationManager.";
  for (auto h : h_list) {
    if (h->has_no_position()) {
      // No position info set in hadronization module
      JSWARN << "Found fragmentation hadron without properly set position in "
                "Afterburner.\nInclusion of fragmentation hadrons only "
                "possible for HybridHadronization.\nExiting.";
      exit(1);
    }
  }
  return h_list;
}

std::vector<std::vector<std::shared_ptr<Hadron>>>
Afterburner::GatherAfterburnerHadrons() {
  std::vector<std::vector<shared_ptr<Hadron>>> afterburner_had_events;
  afterburner_had_events = GetSoftParticlizationHadrons();
  if (GetXMLElementInt({"Afterburner", "include_fragmentation_hadrons"})) {
    if (afterburner_had_events.size() != 1) {
      JSWARN << "Fragmentation hadrons in Afterburner are only possible without "
                "repeated sampling from SoftParticlization. Exiting.";
      exit(1);
    }
    std::vector<shared_ptr<Hadron>> frag_hadrons = GetFragmentationHadrons();
    afterburner_had_events[0].insert(afterburner_had_events[0].end(),
                                     frag_hadrons.begin(), frag_hadrons.end());
  }
  return afterburner_had_events;
}

} // end namespace Jetscape
