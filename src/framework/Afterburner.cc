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
// This is a general basic class for hadronic afterburner

#include "./Afterburner.h"
#include "./JetScapeSignalManager.h"

using namespace std;

namespace Jetscape {
/**
 * @brief Initialize the Afterburner module.
 *
 * Ensures the XML configuration is loaded via the base class, initializes
 * the random number distribution used for smearing fragmentation hadron
 * positions, and calls module-specific initialization routines.
 */
void Afterburner::Init() {
  // Makes sure that XML file with options and parameters is loaded
  JetScapeModuleBase::InitTask();
  JSINFO << "Initializing Afterburner : " << GetId() << " ...";
  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double>{0.0, 1.0};
  InitTask();
  InitTasks();
}

/**
 * @brief Execute the Afterburner task for the current event.
 *
 * This method is invoked by the framework during the task execution phase.
 * It currently emits a verbose log message and can be extended to perform
 * per-event afterburner operations.
 */
void Afterburner::ExecuteTask() {
  VERBOSE(2) << "Afterburner running: " << GetId() << " ...";
}

/**
 * @brief Perform time-stepping calculations for the Afterburner.
 *
 * Emits a verbose log message and delegates to the time-calculation task
 * implementation via `CalculateTimeTask()`.
 */
void Afterburner::CalculateTime() {
  VERBOSE(2) << "Afterburner running for time: " << GetId() << " ...";
  CalculateTimeTask();
}

/**
 * @brief Retrieve hadrons produced by the soft particlization module.
 *
 * Queries the `JetScapeSignalManager` for the soft particlization module
 * pointer. If the module is not present, a warning is emitted and a
 * placeholder (dummy) vector containing an empty hadron vector is returned.
 *
 * @return A vector of hadron-event vectors produced by soft particlization,
 * or a dummy empty event list if the module is not available.
 */
std::vector<std::vector<std::shared_ptr<Hadron>>>
Afterburner::GetSoftParticlizationHadrons() {
  auto soft_particlization =
      JetScapeSignalManager::Instance()->GetSoftParticlizationPointer().lock();
  if (!soft_particlization) {
    JSWARN << "No soft particlization module found. Check if fragmentation"
           << " hadrons are handed to afterburner.";
    std::vector<std::shared_ptr<Hadron>> hadrons;
    dummy.push_back(hadrons);
    return dummy;
  } else {
    return soft_particlization->Hadron_list_;
  }
}

/**
 * @brief Retrieve fragmentation (hard) hadrons for inclusion in the
 * afterburner.
 *
 * Obtains the `HadronizationManager` from the `JetScapeSignalManager`,
 * requests the fragmentation hadron list, applies a small random spatial
 * smearing to avoid exact position overlaps, converts certain kaon
 * identifiers to K0/anti-K0 where appropriate, and filters out partonic
 * entries. If the hadronization manager is missing, the program will exit
 * with an error.
 *
 * @return A vector of shared pointers to fragmentation `Hadron` objects
 * that are suitable for handing to the afterburner.
 */
std::vector<shared_ptr<Hadron>> Afterburner::GetFragmentationHadrons() {
  JSINFO << "Get fragmentation hadrons in Afterburner";
  auto hadronization_mgr = JetScapeSignalManager::Instance()
                               ->GetHadronizationManagerPointer()
                               .lock();
  if (!hadronization_mgr) {
    JSWARN << "No hardronization module found. It is necessary to include"
           << " fragmentation hadrons to afterburner as requested.";
    exit(1);
  }
  std::vector<shared_ptr<Hadron>> h_list;
  hadronization_mgr->GetHadrons(h_list);
  JSINFO << "Got " << h_list.size()
         << " fragmentation hadrons from HadronizationManager.";

  std::vector<shared_ptr<Hadron>> h_list_new;
  rand_int_ptr_ = (std::make_shared<std::uniform_int_distribution<int>>(0, 1));
  for (auto h : h_list) {
    if (h->has_no_position()) {
      JSDEBUG << "Found fragmentation hadron without properly set position in "
                 "Afterburner.\nInclusion of fragmentation hadrons only "
                 "possible for HybridHadronization.";
    }

    // move all the fragmentation hadrons a little bit around to avoid having
    // multiple hadrons at the same position if they are at the same position
    const FourVector r = h->x_in();
    const double rand_x =
        ZeroOneDistribution(*GetMt19937Generator()) * 2e-4 - 1e-4;
    const double rand_y =
        ZeroOneDistribution(*GetMt19937Generator()) * 2e-4 - 1e-4;
    const double rand_z =
        ZeroOneDistribution(*GetMt19937Generator()) * 2e-4 - 1e-4;
    double position_smeared[4] = {r.t(), r.x() + rand_x, r.y() + rand_y,
                                  r.z() + rand_z};
    h->set_x(position_smeared);

    if ((std::abs(h->pid()) > 10) && (h->pid() != 21)) {
      if (h->pstat() > 0) {
        // convert Kaon-L or Kaon-S into K0 or Anti-K0
        if (h->pid() == 310 || h->pid() == 130) {
          const int rand_int = (*rand_int_ptr_)(*GetMt19937Generator());
          const int id = (rand_int == 0) ? 311 : -311;
          h->set_id(id);
        }
        h_list_new.push_back(h);
      } else if (h->pstat() < 0) {
        // convert Kaon-L or Kaon-S into K0 or Anti-K0
        // change id of negative Kaons to make them consistent with the SMASH
        // output
        if (h->pid() == 310 || h->pid() == 130) {
          const int rand_int = (*rand_int_ptr_)(*GetMt19937Generator());
          const int id = (rand_int == 0) ? 311 : -311;
          h->set_id(id);
        }
      }
    } else if ((std::abs(h->pid()) < 10) || (h->pid() == 21)) {
      JSWARN << "Found a free quark or gluon! This can not be handed over to "
                "SMASH.\n"
                "Check confinement in hadronization module!";
    }
  }
  return h_list_new;
}

/**
 * @brief Gather all hadrons that should be processed by the afterburner.
 *
 * Collects hadrons from soft particlization and, optionally depending on
 * configuration flags, includes fragmentation hadrons. Handles configuration
 * options that control whether only final-state hadrons are output and
 * ensures that hadron lists are cleared where necessary to avoid duplicate
 * outputs.
 *
 * @return A vector of event-wise hadron lists prepared for afterburner
 * processing (and eventual hand-off to transport afterburners like SMASH).
 */
std::vector<std::vector<std::shared_ptr<Hadron>>>
Afterburner::GatherAfterburnerHadrons() {
  std::vector<std::vector<shared_ptr<Hadron>>> afterburner_had_events;
  afterburner_had_events = GetSoftParticlizationHadrons();

  if (GetXMLElementInt({"Afterburner", "output_only_final_state_hadrons"})) {
    // clear Hadron_list_ in soft_particlization, otherwise the final hadron
    // output of the writer contains also the soft hadrons which were used as
    // input for SMASH
    auto soft_particlization = JetScapeSignalManager::Instance()
                                   ->GetSoftParticlizationPointer()
                                   .lock();
    if (soft_particlization) {
      soft_particlization->Hadron_list_.clear();
    }
  }

  if (GetXMLElementInt({"Afterburner", "include_fragmentation_hadrons"})) {
    if (afterburner_had_events.size() > 1) {
      JSWARN
          << "Fragmentation hadrons in Afterburner are only possible without "
             "repeated sampling from SoftParticlization. Exiting.";
      exit(1);
    }
    std::vector<shared_ptr<Hadron>> frag_hadrons = GetFragmentationHadrons();

    if (GetXMLElementInt({"Afterburner", "output_only_final_state_hadrons"})) {
      // empty the hadron vector in the hadronization manager to circumvent the
      // output of these hadrons if they are implemented in the SMASH
      // afterburner
      auto hadronization_mgr = JetScapeSignalManager::Instance()
                                   ->GetHadronizationManagerPointer()
                                   .lock();
      hadronization_mgr->DeleteRealHadrons();
    }

    afterburner_had_events[0].insert(afterburner_had_events[0].end(),
                                     frag_hadrons.begin(), frag_hadrons.end());
    dummy.clear();
  }
  return afterburner_had_events;
}

/**
 * @brief Get hadrons produced during the current timestep by bulk dynamics.
 *
 * Queries the `BulkDynamicsManager` for newly produced hadrons and clears
 * the manager's internal buffer of those hadrons. If no bulk manager is
 * available, a warning is logged and an empty list is returned.
 *
 * @return A vector of shared pointers to `Hadron` objects produced this
 * timestep (may be empty).
 */
std::vector<std::shared_ptr<Hadron>>
Afterburner::GetTimestepParticlizationHadrons() {
  auto bdm = JetScapeSignalManager::Instance()->GetBulkPointer().lock();
  if (!bdm) {
    JSWARN
        << "No BulkDynamicsManager module found. Returning empty hadron list.";
    return {};
  }
  return bdm->GetNewHadronsAndClear();
}

/**
 * @brief Get hadrons that should be removed this timestep.
 *
 * Requests the list of hadrons to remove from the `BulkDynamicsManager`.
 * If the bulk manager is not present, logs a warning and returns an empty
 * list.
 *
 * @return A vector of shared pointers to `Hadron` objects that should be
 * removed this timestep (may be empty).
 */
std::vector<std::shared_ptr<Hadron>> Afterburner::GetTimestepHadronsToRemove() {
  auto bdm = JetScapeSignalManager::Instance()->GetBulkPointer().lock();
  if (!bdm) {
    JSWARN
        << "No BulkDynamicsManager module found. Returning empty hadron list.";
    return {};
  }
  return bdm->GetHadronsToRemoveAndClear();
}

}  // end namespace Jetscape
