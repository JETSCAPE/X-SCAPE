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

  /** Extract particles at iso-tau hypersurface
   */
  void ExtractParticlesIsoTau(double tau_surface, std::vector<shared_ptr<Hadron>> &current_hadrons);

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

  /** Switching temperature between media.
   */
  float Tc_;

  /** Switching proper time, when particles from transport initial condition are
   * fed into the hydro 
   */
  double IC_particle_extraction_tau_;

  /**
   * Is the hydro in cartesian or not? Needed to decide whether SMASH IC has to 
   * run first, or if it can run concurrently with hydro.
   */
  bool hydro_Cartesian_;
  bool SMASH_IC_attached_;
  
  bool SMASH_IC_in_progress_;
  bool hydro_in_progress_;

  protected:
    std::weak_ptr<LiquefierBase> liquefier_ptr_;
    std::weak_ptr<HadronicLiquefier> hadronic_liquefier_ptr_;

};

} // end namespace Jetscape

#endif
