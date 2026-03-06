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

#ifndef AFTERBURNER_H
#define AFTERBURNER_H

#include "JetScapeModuleBase.h"
#include "SoftParticlization.h"
#include "HadronizationManager.h"
#include "RealType.h"
#include "BulkMediaInfo.h"
#include "sigslot.h"

namespace Jetscape {

/**
 * @class Afterburner
 * @brief Interface to hadronic afterburner.
 *
 * The Afterburner class gathers hadrons from soft particlization and
 * fragmentation, and performs additional processing before passing them to
 * the final event output.
 */
class Afterburner : public JetScapeModuleBase {
 public:
  /**
   * @brief Default constructor.
   */
  Afterburner() {
    VERBOSE(8);
    SetId("Afterburner");
  }

  /**
   * @brief Destructor.
   */
  ~Afterburner() {
    VERBOSE(8);
    disconnect_all();
  }

  /**
   * @brief Initialize the Afterburner module.
   *
   * This function ensures that required configuration parameters are loaded
   * and initializes necessary random number distributions.
   */
  virtual void Init() override;

  /**
   * @brief Perform the main afterburner work for a timestep.
   *
   * Implementations must define this to carry out the afterburner-specific
   * processing (rescattering, decays, etc.) for the current timestep.
   */
  virtual void ExecuteTask();

  /**
   * @brief Calculate time advancement for the afterburner.
   *
   * Used to compute or update any timing-related internal state between
   * timesteps (e.g. when advancing microscopic afterburner evolution).
   */
  virtual void CalculateTime();

  /**
   * @brief Get the current list of hadrons in the afterburner as Jetscape
   * Hadrons.
   *
   * This must be provided by all concrete afterburner implementations and
   * returns the current snapshot of hadrons managed by the afterburner.
   *
   * @return Vector of `Hadron` objects representing the current state.
   */
  virtual std::vector<Hadron> GetCurrentHadronList() const = 0;

  /**
   * @brief Provide an abstract history object for diagnostics or recording.
   *
   * Default implementation returns the current hadron list; concrete
   * implementations may override to provide richer history structures.
   *
   * @return An `any` containing the history representation.
   */
  virtual any GetHistory() { return GetCurrentHadronList(); }

 protected:
  /**
   * @brief Gather all hadrons from soft particlization and fragmentation.
   *
   * Collects hadrons produced by different sources and merges them into a
   * common structure for further afterburner processing.
   *
   * @return Nested vector of shared pointers to `Hadron` objects (events).
   */
  std::vector<std::vector<std::shared_ptr<Hadron>>> GatherAfterburnerHadrons();

  /**
   * @brief Get the events of soft particlization hadrons.
   *
   * Fetches hadrons generated through the soft particlization module.
   *
   * @return Nested vector of shared pointers to `Hadron` objects (soft events).
   */
  std::vector<std::vector<std::shared_ptr<Hadron>>>
  GetSoftParticlizationHadrons();

  /**
   * @brief Get the list of fragmentation hadrons.
   *
   * Returns hadrons produced by fragmentation/hard processes managed by the
   * `HadronizationManager`.
   *
   * @return Vector of shared pointers to `Hadron` objects.
   */
  std::vector<std::shared_ptr<Hadron>> GetFragmentationHadrons();

  /**
   * @brief Get the list of hadrons for the upcoming timestep from
   * `BulkDynamicsManager`.
   *
   * This call will clear the internal list in the manager and transfer
   * ownership of the hadrons for the timestep processing.
   *
   * @return Vector of shared pointers to `Hadron` objects for the timestep.
   */
  std::vector<std::shared_ptr<Hadron>> GetTimestepParticlizationHadrons();

  /**
   * @brief Get the list of hadrons to be removed for the upcoming timestep
   * from `BulkDynamicsManager`.
   *
   * This will clear the removal list in the manager and return the entries
   * that should be removed from the afterburner.
   *
   * @return Vector of shared pointers to `Hadron` objects to remove.
   */
  std::vector<std::shared_ptr<Hadron>> GetTimestepHadronsToRemove();

  /**< Placeholder storage for hadrons. */
  std::vector<std::vector<std::shared_ptr<Hadron>>> dummy;

  /**< Uniform random number distribution in [0,1] for position smearing. */
  std::uniform_real_distribution<double> ZeroOneDistribution;

  /**< RNG for the Kaon-L / Kaon-S switch to K0 / Anti-K0. */
  std::shared_ptr<std::uniform_int_distribution<int>> rand_int_ptr_;
};

}  // end namespace Jetscape

#endif  // AFTERBURNER_H
