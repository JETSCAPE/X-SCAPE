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

#ifndef BULKDYNAMICSMANAGER_H
#define BULKDYNAMICSMANAGER_H

//#include "JetScapeTask.h"
#include "JetScapeModuleBase.h"
#include "BulkMediaBase.h"
#include "JetClass.h"
#include "FluidCellInfo.h"
#include "BulkMediaInfo.h"
#include "LiquefierBase.h"
#include "HadronicLiquefier.h"
#include "HadronicEMT.h"
#include "sigslot.h"

#include <vector>

namespace Jetscape {
/**
 * @class BulkDynamicsManager
 * @brief Orchestrates the bulk-medium stages and their hand-off to transport.
 *
 * The manager coordinates initial-state hadronic transport, hydrodynamics,
 * soft particlization, and afterburner evolution. It also mediates access to
 * bulk-medium properties and tracks hadron lists exchanged across timesteps.
 */
class BulkDynamicsManager
    : public JetScapeModuleBase,
      public std::enable_shared_from_this<BulkDynamicsManager> {
 public:
  /**
   * @brief Construct a bulk dynamics manager.
   *
   * Initializes the module id to "BulkDynamicsManager".
   */
  BulkDynamicsManager();

  /**
   * @brief Destroy the bulk dynamics manager.
   */
  virtual ~BulkDynamicsManager();

  /**
   * @brief Initialize attached bulk-dynamics tasks and configuration.
   */
  virtual void InitTask();

  /**
   * @brief Execute the manager task-level entry point.
   */
  virtual void ExecuteTask();

  /**
   * @brief Clear attached tasks and signal connections.
   *
   * Can be overridden by derived tasks.
   */
  virtual void ClearTask();

  /**
   * @brief Perform per-timestep calculation for all active child tasks.
   */
  virtual void CalculateTime();

  /**
   * @brief Execute one manager-controlled timestep.
   */
  virtual void ExecTime();

  /**
   * @brief Perform event-level initialization.
   */
  virtual void InitPerEvent();

  /**
   * @brief Perform event-level finalization and cleanup.
   */
  virtual void FinishPerEvent();

  /**
   * @brief Write manager-controlled output via a JetScape writer.
   * @param w Weak pointer to the output writer.
   */
  void WriteTask(weak_ptr<JetScapeWriter> w);

  /**
   * @brief Forward an energy-deposit update to attached media modules.
   * @param t Discrete time index.
   * @param edop Energy deposit to inject.
   */
  void UpdateEnergyDeposit(int t, double edop) {
    UpdateEnergyDepositFromModules(t, edop);
  }

  /**
   * @brief Query energy density from attached media modules.
   * @param t Discrete time index.
   * @param edensity Output reference for energy density.
   */
  void GetEnergyDensity(int t, double &edensity) {
    GetEnergyDensityFromModules(t, edensity);
  }

  /**
   * @brief Query hydro cell information at spacetime point $(t,x,y,z)$.
   * @param t Time coordinate.
   * @param x Spatial x coordinate.
   * @param y Spatial y coordinate.
   * @param z Spatial z coordinate.
   * @param fCell Output fluid-cell container.
   */
  void GetHydroCell(double t, double x, double y, double z,
                    std::unique_ptr<FluidCellInfo> &fCell) {
    GetHydroInfoFromModules(t, x, y, z, fCell);
  }

  /**
   * @brief Query the hydro start time from attached media modules.
   * @param tau0 Output reference to the hydro start proper time.
   */
  void GetHydroStartTime(double &tau0) { GetHydroStartTimeFromModules(tau0); }

  /**
   * @brief Forward energy-deposit update to all compatible child modules.
   * @param t Discrete time index.
   * @param edop Energy deposit to inject.
   */
  void UpdateEnergyDepositFromModules(int t, double edop);

  /**
   * @brief Retrieve energy density from all compatible child modules.
   * @param t Discrete time index.
   * @param edensity Output reference for energy density.
   */
  void GetEnergyDensityFromModules(int t, double &edensity);

  /**
   * @brief Retrieve hydro information from attached modules.
   * @param t Time coordinate.
   * @param x Spatial x coordinate.
   * @param y Spatial y coordinate.
   * @param z Spatial z coordinate.
   * @param fluid_cell_info_ptr Output fluid-cell information.
   */
  void GetHydroInfoFromModules(
      Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
      std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);

  /**
   * @brief Query hydro start time from attached fluid-dynamics modules.
   * @param tau0 Output reference to hydro start proper time.
   */
  void GetHydroStartTimeFromModules(double &tau0);

  /**
   * @brief Query active bulk information with hydro/hadronic fallback logic.
   * @param t Time coordinate.
   * @param x Spatial x coordinate.
   * @param y Spatial y coordinate.
   * @param z Spatial z coordinate.
   * @param fluid_cell_info_ptr Output fluid-cell information.
   */
  void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                   Jetscape::real z,
                   std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);

  /**
   * @brief Convert generic bulk-medium info into `FluidCellInfo`.
   * @param fluid_cell_info_ptr Output fluid-cell information object.
   * @param bulk_info_ptr Input bulk-medium information object.
   */
  void InfoWrapper(std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr,
                   std::unique_ptr<BulkMediaInfo> &bulk_info_ptr);

  /**
   * @brief Get and clear hadrons scheduled to be injected next timestep.
   * @return Vector of hadrons to add to transport.
   */
  std::vector<shared_ptr<Hadron>> GetNewHadronsAndClear();

  /**
   * @brief Get and clear hadrons scheduled for removal next timestep.
   * @return Vector of hadrons to remove from transport.
   */
  std::vector<shared_ptr<Hadron>> GetHadronsToRemoveAndClear();

  /**
   * @brief Mark hadrons that crossed the extraction iso-$\tau$ hypersurface.
   * @param current_hadrons Current hadron list to inspect.
   */
  void DetermineHadronsCrossingIsoTau(
      std::vector<shared_ptr<Hadron>> &current_hadrons);

  /**
   * @brief Extract hadrons from transport IC at the extraction iso-$\tau$.
   *
   * Participant hadrons are stored in
   * `store_source_term_hadrons_iso_tau_`, spectator hadrons in
   * `store_spectator_hadrons_iso_tau_`.
   *
   * @param AllHadronsCrossedIsoTau Output flag indicating completion.
   */
  void ExtractHadronsFromTransportInitialConditionIsoTau(
      bool &AllHadronsCrossedIsoTau);

  /**
   * @brief Free-stream a hadron to a target proper time.
   * @param tau Target proper time.
   * @param hadron Hadron to propagate.
   */
  void PropagateHadronFreeStreamingToTau(double tau,
                                         std::shared_ptr<Hadron> &hadron);

  /**
   * @brief Build hydro source terms from hadrons at extraction iso-$\tau$.
   */
  void CreateHadronicSourceTermsForHydroInitializationIsoTau();

  /**
   * @brief Store newly produced hadrons from soft particlization.
   */
  void StoreHadronsFromSoftParticlization();

  /**
   * @brief Queue one hadron to be injected into transport next timestep.
   * @param new_hadron Hadron to add.
   */
  void AddNewHadron(const shared_ptr<Hadron> &new_hadron) {
    new_hadrons_for_timestep_.push_back(new_hadron);
  }

  /**
   * @brief Queue multiple hadrons to be injected next timestep.
   * @param new_hadrons Hadrons to add.
   */
  void AddNewHadrons(const std::vector<shared_ptr<Hadron>> &new_hadrons) {
    for (const auto &had : new_hadrons) {
      AddNewHadron(had);
    }
  }

  /**
   * @brief Queue one hadron for removal from transport next timestep.
   * @param hadron Hadron to remove.
   */
  void RemoveHadron(const shared_ptr<Hadron> &hadron) {
    remove_hadrons_for_timestep_.push_back(hadron);
  }

  /**
   * @brief Queue multiple hadrons for removal next timestep.
   * @param hadrons Hadrons to remove.
   */
  void RemoveHadrons(const std::vector<shared_ptr<Hadron>> &hadrons) {
    for (const auto &had : hadrons) {
      RemoveHadron(had);
    }
  }

 private:
  /** New hadrons for upcoming timestep of transport evolution,
   * to be filled at end of timestep by particlization routine.
   */
  std::vector<shared_ptr<Hadron>> new_hadrons_for_timestep_;

  /** Hadrons to be removed from the transport evolution in the
   * upcoming timestep, to be determined at the end of timestep
   * by some criterion (iso-tau surface, energy density, ...)
   */
  std::vector<shared_ptr<Hadron>> remove_hadrons_for_timestep_;

  /** Store the hadrons extracted when the SMASH initial condition is used and
   * the hydro does not run in Cartesian coordinates.
   * In this case, SMASH has to run first (up to some large time), and the
   * particles have to be extracted at an iso-tau surface. Afterwards they can
   * be used as source terms in the hydro.
   */
  std::vector<shared_ptr<Hadron>> store_source_term_hadrons_iso_tau_;
  std::vector<shared_ptr<Hadron>> store_spectator_hadrons_iso_tau_;

  /** Store the hadrons from the soft particlization, when the hydro runs in
   * Milne coordinates. Then they are fed into SMASH after the hydro has run.
   * Add also hadrons that have a too large gamma factor to be added to the
   * hydro evolution.
   */
  std::vector<shared_ptr<Hadron>> store_hadrons_soft_particlization_;

  /** Store BDM final state hadrons for output
   */
  std::vector<std::vector<shared_ptr<Hadron>>> BDM_final_state_hadrons_;

  /** Flag to decide if the energy density criterion is used to distinguish
   * between different media.
   */
  bool energy_density_criterion_;

  /** Switching temperature between media.
   */
  float Tc_;

  /** Switching energy density between media
   */
  double ec_;

  /** Switching proper time, when particles from transport initial condition are
   * fed into the hydro
   */
  double IC_particle_extraction_tau_;

  /** Flag and pointer to create a file output of the hadronic time evolution.
   * This can be used to create a video of the hadronic evolution.
   */
  bool hadronic_time_evolution_to_file_;
  std::unique_ptr<ofstream> hadronic_time_evolution_file_;

  /**
   * @brief Create hadronic time-evolution output file when enabled.
   */
  void CreateHadronicTimeEvolutionFileIfNecessary() {
    if (hadronic_time_evolution_to_file_) {
      hadronic_time_evolution_file_ =
          std::make_unique<ofstream>("hadronic_time_evolution_BDM.dat");
    }
  }

  /**
   * @brief Close hadronic time-evolution output file when enabled.
   */
  void CloseHadronicTimeEvolutionFileIfNecessary() {
    if (hadronic_time_evolution_to_file_) {
      hadronic_time_evolution_file_->close();
    }
  }

  /**
   * @brief Write hadronic state snapshot to the optional evolution file.
   */
  void PrintHadronicTimeEvolutionToFileIfNecessary();

  /**
   * Is the hydro in cartesian or not? Needed to decide whether SMASH IC has to
   * run first, or if it can run concurrently with hydro.
   */
  bool hydro_Cartesian_;
  bool SMASH_IC_attached_;

  bool SMASH_IC_in_progress_;
  bool reset_time_hydro_Milne_;
  bool hydro_in_progress_;
  bool afterburner_in_progress_;

  double deltaT_main_clock_;

  /**
   * HadronicEMT object to obtain bulk quantities of hadronic media.
   */
  HadronicEMT hadronic_emt_;

  /**
   * Kinematic cuts for the hadrons to be added to the hydrodynamic evolution.
   */
  bool enforce_pT_cut_;
  bool enforce_rapidity_cut_;
  bool ignore_spectator_hadrons_;
  double pT_cut_;
  double rapidity_cut_;
  int counter_hadrons_already_added_;

 protected:
  std::weak_ptr<LiquefierBase> liquefier_ptr_;
  std::weak_ptr<HadronicLiquefier> hadronic_liquefier_ptr_;

  std::uniform_real_distribution<double> ZeroOneDistribution;
};

}  // end namespace Jetscape

#endif
