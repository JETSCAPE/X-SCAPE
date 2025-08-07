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
/** @class Bulk dynamics manager manager.
   */
class BulkDynamicsManager
    : public JetScapeModuleBase,
      public std::enable_shared_from_this<BulkDynamicsManager> {

public:
  /** Default constructor to create a bulk dynamics manager. Sets task ID as "BulkDynamicsManager".
   */
  BulkDynamicsManager();

  /** Destructor for the bulk dynamics manager.
   */
  virtual ~BulkDynamicsManager();

  /** It initializes the tasks attached to the bulk dynamics manager.
   */
  virtual void InitTask();

  /**
  */
  virtual void ExecuteTask();

  /** It erases the tasks attached with the bulk dynamics manager. It can be overridden by other tasks.
   */
  virtual void ClearTask();

  virtual void CalculateTime();

  virtual void ExecTime();

  virtual void InitPerEvent();

  virtual void FinishPerEvent();

  void WriteTask(weak_ptr<JetScapeWriter> w);

  void UpdateEnergyDeposit(int t, double edop){ UpdateEnergyDepositFromModules(t, edop); }

  void GetEnergyDensity(int t, double &edensity){ GetEnergyDensityFromModules(t, edensity); }

  void GetHydroCell(double t, double x, double y, double z,
                            std::unique_ptr<FluidCellInfo> &fCell) { GetHydroInfoFromModules(t, x, y, z, fCell); }

  void GetHydroStartTime(double &tau0){GetHydroStartTimeFromModules(tau0); }

  void UpdateEnergyDepositFromModules(int t, double edop);

  void GetEnergyDensityFromModules(int t, double &edensity);

  void GetHydroInfoFromModules(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
			    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);

  void GetHydroStartTimeFromModules(double &tau0);

  void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
                            std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);

  void InfoWrapper(std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr, std::unique_ptr<BulkMediaInfo> &bulk_info_ptr);

  /** Get the new hadrons for the upcoming timesteps and clear the vector for the next timestep
   */
  std::vector<shared_ptr<Hadron>> GetNewHadronsAndClear();

  /** Get the hadrons to be removed in the upcoming timestep and clear the vector for the next timestep
   */
  std::vector<shared_ptr<Hadron>> GetHadronsToRemoveAndClear();

  /** Determine if a hadron has crossed the iso-tau hypersurface
   */
  void DetermineHadronsCrossingIsoTau(std::vector<shared_ptr<Hadron>> &current_hadrons);

  /** Extract hadrons from the transport initial condition at the iso-tau surface.
   * This adds participant hadrons to store_source_term_hadrons_iso_tau_
   * and spectator hadrons to store_spectator_hadrons_iso_tau_.
  */
  void ExtractHadronsFromTransportInitialConditionIsoTau(bool &AllHadronsCrossedIsoTau);

  /** Create hadronic source terms for hydro initialization from the hadrons
   * at an iso-tau surface.
  */
  void CreateHadronicSourceTermsForHydroInitializationIsoTau();

  /** Store hadrons from the soft particlization in store_hadrons_soft_particlization_
  */
  void StoreHadronsFromSoftParticlization();

  /** Add one new hadron for transport to hadron list for new timestep
   */
  void AddNewHadron(const shared_ptr<Hadron>& new_hadron) {
    new_hadrons_for_timestep_.push_back(new_hadron);
  }

  /** Add list of new hadrons for transport to hadron list for new timestep
   */
  void AddNewHadrons(const std::vector<shared_ptr<Hadron>>& new_hadrons) {
    for (const auto& had : new_hadrons) {
      AddNewHadron(had);
    }
  }

  /** Remove one hadron from transport, add it to remove hadron list for new timestep
   */
  void RemoveHadron(const shared_ptr<Hadron>& hadron) {
    remove_hadrons_for_timestep_.push_back(hadron);
  }

  /** Remove list of hadrons from transport, add it to remove hadron list for new timestep
   */
  void RemoveHadrons(const std::vector<shared_ptr<Hadron>>& hadrons) {
    for (const auto& had : hadrons) {
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

  /** Function to create the hadronic_time_evolution_file_ 
   * if the flag hadronic_time_evolution_to_file_ is set to true.
  */
  void CreateHadronicTimeEvolutionFileIfNecessary() {
    if (hadronic_time_evolution_to_file_) {
      hadronic_time_evolution_file_ = 
        std::make_unique<ofstream>("hadronic_time_evolution_BDM.dat");
    }
  }

  /** Function to close the hadronic_time_evolution_file_ 
   * if the flag hadronic_time_evolution_to_file_ is set to true.
  */
  void CloseHadronicTimeEvolutionFileIfNecessary() {
    if (hadronic_time_evolution_to_file_) {
      hadronic_time_evolution_file_->close();
    }
  }

  /** Function to print the hadronic content of the time evolution to the file
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

} // end namespace Jetscape

#endif
