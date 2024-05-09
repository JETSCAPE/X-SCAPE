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

#include "BulkDynamicsManager.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include "MakeUniqueHelper.h"
#include "QueryHistory.h"
#include <string>

#include <iostream>
#include <vector>
#include <thread>

using namespace std;

namespace Jetscape {

BulkDynamicsManager::BulkDynamicsManager() : JetScapeModuleBase() {
  SetId("BulkDynamicsManager");
  VERBOSE(8);
}

BulkDynamicsManager::~BulkDynamicsManager() {
  // Check if this is all really needed with shared_ptr ...
  JSDEBUG;
  ClearTask();

  if (GetNumberOfTasks() > 0)
    EraseTaskLast();
}

void BulkDynamicsManager::ClearTask() {
  JSDEBUG << "BulkDynamicsManager ClearTask() ...";

  int n = GetNumberOfTasks();
  for (int i = 1; i < n; i++)
    EraseTaskLast();

  // Clean Up not really working with iterators (see also above!!!) Some logic not clear for me.
  JetScapeSignalManager::Instance()->CleanUp();
}

void BulkDynamicsManager::InitTask() {
  JSINFO << "Intialize BulkDynamicsManager ...";

  //Critical temperature to switch from hydro to something else
  Tc_ = GetXMLElementDouble({"BDM", "Tc"});
  IC_particle_extraction_tau_ = GetXMLElementDouble({"BDM", "IC_particle_extraction_tau"});
  hydro_Cartesian_ = false;
  std::string strCartesianHydro = GetXMLElementText({"Hydro", "CartesianHydro"});
  if ((int)strCartesianHydro.find("true") >= 0) {
    hydro_Cartesian_ = true;
    JSINFO << "BulkDynamicsManager set up for run with Cartesian hydro ...";
  } else {
    JSINFO << "BulkDynamicsManager set up for run with Milne hydro ...";
  }

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid bulk dynamics Manager modules found ...";
    exit(-1);
  }

  for (auto task : GetTaskList()) {
    // Check if the task is an instance of FluidDynamics, then set liquefier_ptr
    if (auto fluidDynamics = std::dynamic_pointer_cast<FluidDynamics>(task)) {
        liquefier_ptr_ = fluidDynamics->get_liquefier();
        hadronic_liquefier_ptr_ = fluidDynamics->get_hadronic_liquefier();
    }
  }

  JSINFO << "Found " << GetNumberOfTasks()
         << " Bulk Dynamics Manager Tasks/Modules Initialize them ... ";

  for(auto it : GetTaskList()) {
        auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
        JSWARN << "Module in task list: " << module->GetId();
  }
}

void BulkDynamicsManager::ExecuteTask() {
  VERBOSE(1) << "Run BulkDynamicsManager Manager ...";
  JSDEBUG << "Task Id = " << this_thread::get_id();

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid Bulk Dynamics Manager modules found ...";
    exit(-1);
  }
}

void BulkDynamicsManager::CalculateTime()
{
  VERBOSE(3) << "Calculate Bulk Dynamics Manager per timestep ... Current Time = " 
            << GetModuleCurrentTime();
  VERBOSE(3) << "Task Id = " << this_thread::get_id();

  VERBOSE(3) << "Size of new hadron list at beginning of CalculateTime in BDM (should be something) = " 
            << new_hadrons_for_timestep_.size();

  JetScapeModuleBase::CalculateTimeTasks();

  VERBOSE(3) << "Size of new hadron list at end of CalculateTime in BDM (should be empty) = " 
            << new_hadrons_for_timestep_.size();

}

void BulkDynamicsManager::ExecTime()
{
  VERBOSE(3) << "Execute Bulk Dynamics Manager at timestep (end) ... Current Time = "
            << GetModuleCurrentTime() << " Thread Id = " << this_thread::get_id();
  VERBOSE(3) << "Task Id = " << this_thread::get_id();
 
  if(SMASH_IC_attached_) {

    if (SMASH_IC_in_progress_) {
      VERBOSE(3) << "Size of new hadron list at beginning of ExecTime (should be empty) = " 
                  << new_hadrons_for_timestep_.size();
      
      linb::any current_hadrons_IC = QueryHistory::Instance()->GetHistoryFromModule("SMASHInitialState");
      std::vector<Hadron> hadrons = any_cast<std::vector<Hadron>>(current_hadrons_IC);

      // Convert to vector of shared pointers using std::transform
      std::vector<std::shared_ptr<Hadron>> shared_hadrons;
      std::transform(hadrons.begin(), hadrons.end(), std::back_inserter(shared_hadrons),
                  [](const Hadron& h) { return std::make_shared<Hadron>(h); });

      // Determine which particles should be removed from the SMASH initial condition
      ExtractParticlesIsoTau(IC_particle_extraction_tau_, shared_hadrons);

      if(!hydro_Cartesian_) {
        for(const auto& hadron : remove_hadrons_for_timestep_) {
          JSINFO << "Remove hadron with participant status " << hadron->participant();
          if (hadron->participant()) {
            store_source_term_hadrons_iso_tau_.push_back(hadron);
          } else {
            store_spectator_hadrons_iso_tau_.push_back(hadron);
          }
        }
      }
      JSINFO << "Currently " << store_source_term_hadrons_iso_tau_.size() << " source term hadrons in storage.";
      JSINFO << "Currently " << store_spectator_hadrons_iso_tau_.size() << " spectator hadrons in storage.";

      if(shared_hadrons.empty()) {
        SMASH_IC_in_progress_ = false;
        hydro_in_progress_ = true;
      }

      if(!hydro_Cartesian_ && !SMASH_IC_in_progress_) {
        JSINFO << "SMASH IC is empty, resetting time to " << IC_particle_extraction_tau_-GetMainClock()->GetDeltaT();
        GetMainClock()->ResetToTime(IC_particle_extraction_tau_-GetMainClock()->GetDeltaT());
        JSINFO << "Time reset to " << GetMainClock()->GetCurrentTime();

        // Set SMASH IC to inactive and the other modules to active
        for(auto it : GetTaskList()) {
          auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
          if(dynamic_pointer_cast<SmashInitialConditionWrapper>(module)) {
            JSWARN << "SetActive(false) = " << module->GetId();
            module->SetActive(false);
          } else {
            JSWARN << "SetActive(true) = " << module->GetId();
            module->SetActive(true);
          }
        }

        // Check if the hydro is activated
        // ONLY FOR TESTING
        bool hydro_activated = false;
        for(auto it : GetTaskList()) {
          auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
          if(dynamic_pointer_cast<FluidDynamics>(module)) {
            if(module->GetActive()) {
              hydro_activated = true;
              JSWARN << "Hydro is active " << hydro_activated;
              break;
            }
          }
        }
        // ONLY FOR TESTING END
      }
    }
    
    // if the hydro is in Milne coordinates, and still running
    if(!hydro_Cartesian_ && hydro_in_progress_) {
      if(!weak_ptr_is_uninitialized(hadronic_liquefier_ptr_)) {
        hadronic_liquefier_ptr_.lock()->clear_hadron_droplet_list();
      }

      // If the hydro is activated, then add the hadrons to the hydro source terms
      // Create hydro source terms from the removed hadrons
      bool add_sources_at_current_tau = false;
      if (std::abs(GetMainClock()->GetCurrentTime() - IC_particle_extraction_tau_) < (GetMainClock()->GetDeltaT() - rounding_error)) {
        add_sources_at_current_tau = true;
      }
      JSWARN << "Add sources at current tau = " << add_sources_at_current_tau;
      if (!weak_ptr_is_uninitialized(hadronic_liquefier_ptr_) && add_sources_at_current_tau) {
        JSWARN << "Create a hydro source with " 
                  << store_source_term_hadrons_iso_tau_.size() << " hadrons";

        // Convert std::shared_ptr<Jetscape::Hadron> to raw pointers and store in a vector
        std::vector<Jetscape::Hadron> hadron_objects;
        for (const auto& ptr : store_source_term_hadrons_iso_tau_) {
            hadron_objects.push_back(*ptr);
        }
        hadronic_liquefier_ptr_.lock()->add_hydro_sources_hadrons(IC_particle_extraction_tau_,hadron_objects);

        JSWARN << "Number of hadronic droplets: " << hadronic_liquefier_ptr_.lock()->get_dropletlist_size();
      }

    }

  }

  // Get the soft particlization hadrons from iSS for the further evolution in SMASH
  // This function fills new_hadrons_for_timestep_, which is handed to SMASH




  JetScapeModuleBase::ExecTimeTasks();

  VERBOSE(3) << "Size of new hadron list at end of ExecTime (should be something) = " 
            << new_hadrons_for_timestep_.size();
}

void BulkDynamicsManager::InitPerEvent()
{
  VERBOSE(3) << "InitPerEvent Bulk Dynamics Manager when used per timestep ...";
  VERBOSE(3) << "Task Id = " << this_thread::get_id();

  JetScapeModuleBase::InitPerEventTasks();

  /**
   * If SMASH IC is attached to BDM and the hydro runs in Milne coordinates, 
   * then we have to set all other modules to inactive and run SMASH first until 
   * it is empty. Then the time is reset to 0 and the other modules can run.
  */
  SMASH_IC_attached_ = false;
  if (!hydro_Cartesian_) {
    for(auto it : GetTaskList()) {
      auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
      if(dynamic_pointer_cast<SmashInitialConditionWrapper>(module)) {
        JSWARN << "SetActive(true) = " << module->GetId();
        SMASH_IC_attached_ = true;
        SMASH_IC_in_progress_ = true;
        hydro_in_progress_ = false;
      } else {
        JSWARN << "SetActive(false) = " << module->GetId();
        module->SetActive(false);
      }
    }
  }

}

void BulkDynamicsManager::FinishPerEvent()
{
  VERBOSE(3) << "FinishPerEvent Bulk Dynamics Manager when used per timestep ...";
  VERBOSE(3) << "Task Id = " << this_thread::get_id();

  JetScapeModuleBase::FinishPerEventTasks();

  //JP: Quick fix, to be discussed, similar to writer, clear is only called for active tasks, so call here directly ...
  ClearTask();
}

void BulkDynamicsManager::UpdateEnergyDepositFromModules(int t, double edop){

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid bulk manager modules found ...";
    exit(-1);
  }
  for (auto it : GetTaskList()) {
    if(dynamic_pointer_cast<FluidDynamics>(it))dynamic_pointer_cast<FluidDynamics>(it)->UpdateEnergyDeposit(t,edop);
  }
}

void BulkDynamicsManager::GetEnergyDensityFromModules(int t, double &edensity){
  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid bulk manager modules found ...";
    exit(-1);
  }
  for (auto it : GetTaskList()) {
    if(dynamic_pointer_cast<FluidDynamics>(it))dynamic_pointer_cast<FluidDynamics>(it)->GetEnergyDensity(t,edensity);
  }
}

void BulkDynamicsManager::GetHydroInfoFromModules(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
						    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr){
  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid bulk manager modules found ...";
    exit(-1);
  }
  //If and only if there is one media module and it is hydro do this like JETSCAPE
  if(GetNumberOfTasks() == 1){
    for (auto it : GetTaskList()) {
      if(dynamic_pointer_cast<FluidDynamics>(it)) {
        dynamic_pointer_cast<FluidDynamics>(it)->GetHydroInfo(t,x,y,z,fluid_cell_info_ptr);
      } else {
        GetBulkInfo(t,x,y,z,fluid_cell_info_ptr);
      }
    }
  }
  else
    GetBulkInfo(t,x,y,z,fluid_cell_info_ptr);
}

void BulkDynamicsManager::GetHydroStartTimeFromModules(double &tau0){
  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid bulk manager modules found ...";
    exit(-1);
  }
  for (auto it : GetTaskList()) {
    if(dynamic_pointer_cast<FluidDynamics>(it))dynamic_pointer_cast<FluidDynamics>(it)->GetHydroStartTime(tau0);
  }
}

void BulkDynamicsManager::GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
                                                    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr){

  bool validHydro = false;

  //Need a cleaner way of getting media info
  //Would be great place to implement std::variant
  //variant<std::unique_ptr<FluidCellInfo>,std::unique_ptr<BulkMediaInfo>> info;

  for (auto it : GetTaskList()) {
    if(dynamic_pointer_cast<FluidDynamics>(it)){
      dynamic_pointer_cast<FluidDynamics>(it)->GetHydroInfo(t,x,y,z,fluid_cell_info_ptr);
      if(fluid_cell_info_ptr->temperature > Tc_) validHydro = true;
    }
  }
  //if validHydro = true, we are done; if not get info from other modules
  if(validHydro == false){
    std::unique_ptr<BulkMediaInfo> bulk_info_ptr;
    // this has to be fixed with the new HadronicEMT class
    /*for (auto it : GetTaskList()) {
    for (auto it : GetTaskList()) {
      if(dynamic_pointer_cast<Afterburner>(it)){
	      dynamic_pointer_cast<Afterburner>(it)->GetBulkInfo(t,x,y,z,bulk_info_ptr);
      }
    }*/
    InfoWrapper(fluid_cell_info_ptr,bulk_info_ptr);
  }
}

void BulkDynamicsManager::InfoWrapper(std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr,std::unique_ptr<BulkMediaInfo> &bulk_info_ptr){
  fluid_cell_info_ptr = make_unique<FluidCellInfo>();
  fluid_cell_info_ptr->temperature = bulk_info_ptr->temperature;
  fluid_cell_info_ptr->pressure = bulk_info_ptr->pressure;
  fluid_cell_info_ptr->entropy_density = bulk_info_ptr->entropy_density;
  fluid_cell_info_ptr->energy_density = bulk_info_ptr->energy_density;
  fluid_cell_info_ptr->qgp_fraction = bulk_info_ptr->qgp_fraction;
  fluid_cell_info_ptr->mu_B = bulk_info_ptr->mu_B;
  fluid_cell_info_ptr->mu_C = bulk_info_ptr->mu_C;
  fluid_cell_info_ptr->mu_S = bulk_info_ptr->mu_S;
  fluid_cell_info_ptr->vx = bulk_info_ptr->vx;
  fluid_cell_info_ptr->vy = bulk_info_ptr->vy;
  fluid_cell_info_ptr->vz = bulk_info_ptr->vz;
  fluid_cell_info_ptr->bulk_Pi = bulk_info_ptr->bulk_Pi;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      fluid_cell_info_ptr->pi[i][j] = bulk_info_ptr->pi[i][j];
    }
  }
  // tmn from bulk info not converted as not present in fluid cell info
}

std::vector<shared_ptr<Hadron>> BulkDynamicsManager::GetNewHadronsAndClear() {
  std::vector<shared_ptr<Hadron>> new_h_to_return;
  // The swap puts the empty vector for new_hadrons_for_timestep_
  // and therefore clears the vector (to be filled again at next timestep)
  new_h_to_return.swap(new_hadrons_for_timestep_);
  return new_h_to_return;
}

std::vector<shared_ptr<Hadron>> BulkDynamicsManager::GetHadronsToRemoveAndClear() {
  std::vector<shared_ptr<Hadron>> new_h_to_remove;
  // The swap puts the empty vector for remove_hadrons_for_timestep_
  // and therefore clears the vector (to be filled again at the next timestep)
  new_h_to_remove.swap(remove_hadrons_for_timestep_);
  return new_h_to_remove;
}

void BulkDynamicsManager::ExtractParticlesIsoTau(double tau_surface, 
                          std::vector<shared_ptr<Hadron>> &current_hadrons) {
  int i = 0;
  for (const auto& had : current_hadrons) {
    const FourVector r = had->x_in();
    const double t = r.t();
    const double z = r.z();
    const double tau = sqrt(t*t - z*z);

    if(tau >= tau_surface) {
      i++;
      RemoveHadron(had);
    }
  }
  JSINFO << "Found " << i << " hadrons to remove from SMASH";
  JSINFO << "Size remove_hadrons_for_timestep_ = " << remove_hadrons_for_timestep_.size();
}

} // end namespace Jetscape
