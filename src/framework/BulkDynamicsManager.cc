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
  JSINFO << "Initialize BulkDynamicsManager ...";

  ZeroOneDistribution = uniform_real_distribution<double>{0.0, 1.0};

  //Critical temperature to switch from hydro to something else
  Tc_ = GetXMLElementDouble({"BDM", "Tc"});
  pT_cut_ = GetXMLElementDouble({"BDM", "pT_cut"});
  enforce_pT_cut_ = false;
  if (pT_cut_ > rounding_error) {
    enforce_pT_cut_ = true;
    JSINFO << "BulkDynamicsManager set up with pT cut = " << pT_cut_ << " ...";
  } else {
    JSINFO << "BulkDynamicsManager set up without pT cut ...";
  }
  rapidity_cut_ = GetXMLElementDouble({"BDM", "rapidity_cut"});
  enforce_rapidity_cut_ = false;
  if (rapidity_cut_ > rounding_error) {
    enforce_rapidity_cut_ = true;
    JSINFO << "BulkDynamicsManager set up with rapidity cut = " << rapidity_cut_ << " ...";
  } else {
    JSINFO << "BulkDynamicsManager set up without rapidity cut ...";
  }
  IC_particle_extraction_tau_ = GetXMLElementDouble({"BDM", "IC_particle_extraction_tau"});
  hydro_Cartesian_ = false;
  std::string strCartesianHydro = GetXMLElementText({"Hydro", "CartesianHydro"});
  if ((int)strCartesianHydro.find("true") != std::string::npos) {
    hydro_Cartesian_ = true;
    JSINFO << "BulkDynamicsManager set up for run with Cartesian hydro ...";
  } else {
    JSINFO << "BulkDynamicsManager set up for run with Milne hydro ...";
  }
  hadronic_time_evolution_to_file_ = false;
  std::string strHadronicTimeEvolutionToFile = GetXMLElementText({"BDM", "CreateHadronicTimeEvolutionOutput"});
  if ((int)strHadronicTimeEvolutionToFile.find("true") != std::string::npos) {
    hadronic_time_evolution_to_file_ = true;
    JSINFO << "BulkDynamicsManager set up to write hadronic time evolution to file ...";
  } else {
    JSINFO << "BulkDynamicsManager set up to not write hadronic time evolution to file ...";
  }

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid bulk dynamics Manager modules found ...";
    exit(-1);
  }

  bool hydro_module_attached = false;
  for (auto task : GetTaskList()) {
    // Check if the task is an instance of FluidDynamics, then set liquefier_ptr
    if (auto fluidDynamics = std::dynamic_pointer_cast<FluidDynamics>(task)) {
      liquefier_ptr_ = fluidDynamics->get_liquefier();
      hadronic_liquefier_ptr_ = fluidDynamics->get_hadronic_liquefier();

      JetScapeSignalManager::Instance()->SetHydroPointer(
          dynamic_pointer_cast<FluidDynamics>(fluidDynamics));
      hydro_module_attached = true;
    }
    if (auto particlization = std::dynamic_pointer_cast<SoftParticlization>(task)) {
      JSWARN << "Connect signals for SoftParticlization";
      JetScapeSignalManager::Instance()->ConnectGetHydroHyperSurfaceSignal(particlization);
      JetScapeSignalManager::Instance()->ConnectClearHydroHyperSurfaceSignal(particlization);
    }
  }

  for (auto task : GetTaskList()) {
    if (auto hadronization = std::dynamic_pointer_cast<SoftParticlization>(task)) {
      JSWARN << "Connection check HydroHyperSurfaceConnected = " << hadronization->GetGetHydroHyperSurfaceConnected();
      JSWARN << "Connection check ClearHydroHyperSurfaceConnected = " << hadronization->GetClearHydroHyperSurfaceConnected();
    }
  }

  JSINFO << "Found " << GetNumberOfTasks()
         << " Bulk Dynamics Manager Tasks/Modules Initialize them ... ";

  for(auto it : GetTaskList()) {
        auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
        JSWARN << "Module in task list: " << module->GetId();
  }

  // Create the hadronic_emt_ object
  if (hydro_module_attached) {
    HadronicEMT hadronic_emt_;
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
  //JSWARN << "Execute Bulk Dynamics Manager at timestep (end) ... Current Time = "
  //          << GetModuleCurrentTime() << " Thread Id = " << this_thread::get_id();
  VERBOSE(3) << "Task Id = " << this_thread::get_id();
  VERBOSE(3) << "Size of new hadron list at beginning of ExecTime (should be empty) = " 
                  << new_hadrons_for_timestep_.size();

  PrintHadronicTimeEvolutionToFileIfNecessary();

  // Check if the SMASH IC is attached and the hydro runs in Milne coordinates
  if (SMASH_IC_attached_ && !hydro_Cartesian_) {
    // If the SMASH IC is still running, check if hadrons have to be removed
    // according to the iso-tau surface criterion
    bool AllHadronsCrossedIsoTau = false;
    if (SMASH_IC_in_progress_) {
      ExtractHadronsFromTransportInitialConditionIsoTau(AllHadronsCrossedIsoTau);
      //JSINFO << "Currently " << store_source_term_hadrons_iso_tau_.size() << " source term hadrons in storage.";
      //JSINFO << "Currently " << store_spectator_hadrons_iso_tau_.size() << " spectator hadrons in storage.";
    }

    // If all hadrons have crossed the iso-tau surface, then the SMASH IC is
    // done and the hydro can start.
    // Therefore, set SMASH IC to inactive and the hydro and soft particlization
    // to active. Also reset the time to close before the iso-tau surface.
    if (AllHadronsCrossedIsoTau && !reset_time_hydro_Milne_) {
      SMASH_IC_in_progress_ = false;

      // Set the SMASH IC to inactive and hydro and soft particlization to active
      for(auto it : GetTaskList()) {
        auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
        if (dynamic_pointer_cast<FluidDynamics>(module)) {
          JSWARN << "SetActive(true) = " << module->GetId();
          module->SetActive(true);
          reset_time_hydro_Milne_ = true;
          hydro_in_progress_ = true;
        } else if (dynamic_pointer_cast<SoftParticlization>(module)) {
          JSWARN << "SetActive(true) = " << module->GetId();
          module->SetActive(true);
        } else {
          JSWARN << "SetActive(false) = " << module->GetId();
          module->SetActive(false);
        }
      }

      // Reset the deltaT of the main clock to the original value after the 
      // SMASH IC run
      GetMainClock()->SetDeltaT(deltaT_main_clock_);

      JSINFO << "SMASH IC is empty, resetting time to " << IC_particle_extraction_tau_-GetMainClock()->GetDeltaT();
      GetMainClock()->ResetToTime(IC_particle_extraction_tau_-GetMainClock()->GetDeltaT());
      JSINFO << "Time reset to " << GetMainClock()->GetCurrentTime();

      // Create the hadronic source terms for the hydro
      CreateHadronicSourceTermsForHydroInitializationIsoTau();
    }

    // If the hydro is running, check if it is finished
    if (hydro_in_progress_) {
      for(auto it : GetTaskList()) {
        auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
        if(auto fluid_dynamics = dynamic_pointer_cast<FluidDynamics>(module)) {
          if(fluid_dynamics->GetHydroStatus() == FINISHED) {
            hydro_in_progress_ = false;
            afterburner_in_progress_ = true;
            JSWARN << "Hydro finished";
            break;
          }
        }
      }

      // In case the hydro is still running, check for produced hadrons and store them
      StoreHadronsFromSoftParticlization();
      JSWARN << "Size of store_hadrons_soft_particlization_ = " << store_hadrons_soft_particlization_.size();

      // If the hydro is finished, then set SMASH as active again
      afterburner_in_progress_ = false;
      if (!SMASH_IC_in_progress_ && !hydro_in_progress_ && !afterburner_in_progress_) {
        for(auto it : GetTaskList()) {
          auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
          if(dynamic_pointer_cast<Afterburner>(module) && (module->GetId() == "SMASH")) {
            JSWARN << "SetActive(true) = " << module->GetId();
            module->SetActive(true);
            afterburner_in_progress_ = true;
          } else {
            JSWARN << "SetActive(false) = " << module->GetId();
            module->SetActive(false);
          }
        }

        JSINFO << "Hydro is done, resetting time to " << IC_particle_extraction_tau_;
        GetMainClock()->ResetToTime(IC_particle_extraction_tau_);
        JSWARN << "Time reset to " << GetMainClock()->GetCurrentTime();
        JSWARN << "End of the next timestep is " << GetMainClock()->GetCurrentTime()+GetMainClock()->GetDeltaT();
      }
    }

    if (afterburner_in_progress_) {
      // check the hadrons in store_spectator_hadrons_iso_tau_ and store_hadrons_soft_particlization_ if they have times larger than the current time
      // and smaller than the current time + deltaT, then add them to the new_hadrons_for_timestep_ list
      for (const auto& had : store_spectator_hadrons_iso_tau_) {
        const FourVector r = had->x_in();
        const double t = r.t();

        if((t >= GetMainClock()->GetCurrentTime()) 
          && (t < GetMainClock()->GetCurrentTime()+GetMainClock()->GetDeltaT())) {
          new_hadrons_for_timestep_.push_back(had);
        }
      }

      for (const auto& had : store_hadrons_soft_particlization_) {
        const FourVector r = had->x_in();
        const double t = r.t();

        if((t >= GetMainClock()->GetCurrentTime()) 
          && (t < GetMainClock()->GetCurrentTime()+GetMainClock()->GetDeltaT())) {
          new_hadrons_for_timestep_.push_back(had);
        }
      }
    }
  }

  JetScapeModuleBase::ExecTimeTasks();

  VERBOSE(3) << "Size of new hadron list at end of ExecTime (should be something) = " 
            << new_hadrons_for_timestep_.size();
}

void BulkDynamicsManager::InitPerEvent()
{
  VERBOSE(3) << "InitPerEvent Bulk Dynamics Manager when used per timestep ...";
  VERBOSE(3) << "Task Id = " << this_thread::get_id();

  JSWARN << "InitPerEvent Bulk Dynamics Manager when used per timestep ...";
  JetScapeModuleBase::InitPerEventTasks();

  /**
   * If SMASH IC is attached to BDM and the hydro runs in Milne coordinates, 
   * then we have to set all other modules to inactive and run SMASH first until 
   * it is empty. Then the time is reset and the other modules can run.
  */
  SMASH_IC_attached_ = false;
  SMASH_IC_in_progress_ = false;
  reset_time_hydro_Milne_ = false;
  if (!hydro_Cartesian_) {
    // store the deltaT of the main clock
    deltaT_main_clock_ = GetMainClock()->GetDeltaT();
    // reset the deltaT of the main clock to small value for high iso-tau
    // extraction accuracy
    GetMainClock()->SetDeltaT(0.01);

    for(auto it : GetTaskList()) {
      auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
      if(dynamic_pointer_cast<SmashInitialConditionWrapper>(module)) {
        JSWARN << "SetActive(true) = " << module->GetId();
        SMASH_IC_attached_ = true;
        SMASH_IC_in_progress_ = true;
        hydro_in_progress_ = false;
        afterburner_in_progress_ = false;
      } else {
        JSWARN << "SetActive(false) = " << module->GetId();
        module->SetActive(false);
      }
    }
    // Clear the hadron droplet list in the hadronic liquefier
    if(!weak_ptr_is_uninitialized(hadronic_liquefier_ptr_)) {
      hadronic_liquefier_ptr_.lock()->clear_hadron_droplet_list();
    }
  }

  // JUST FOR CHECKING
  for (auto task : GetTaskList()) {
    if (auto hadronization = std::dynamic_pointer_cast<SoftParticlization>(task)) {
      JSWARN << "Connection check HydroHyperSurfaceConnected InitPerEvent = " << hadronization->GetGetHydroHyperSurfaceConnected();
      JSWARN << "Connection check ClearHydroHyperSurfaceConnected InitPerEvent = " << hadronization->GetClearHydroHyperSurfaceConnected();
    }
  }

  CreateHadronicTimeEvolutionFileIfNecessary();
}

void BulkDynamicsManager::FinishPerEvent()
{
  VERBOSE(3) << "FinishPerEvent Bulk Dynamics Manager when used per timestep ...";
  VERBOSE(3) << "Task Id = " << this_thread::get_id();

  JSWARN << "FinishPerEvent Bulk Dynamics Manager when used per timestep ...";
  JSWARN << "Task Id = " << this_thread::get_id();

  JetScapeModuleBase::FinishPerEventTasks();

  // If the SMASH IC is attached and the hydro runs in Milne coordinates, then
  // fill the BDM_final_state_hadrons_ vector with the final state hadrons for the event
  if (SMASH_IC_attached_ && !hydro_Cartesian_) {
    // get the hadrons from the SMASH afterburner
    linb::any current_hadrons_afterburner = QueryHistory::Instance()->GetHistoryFromModule("SMASH");
    if (!current_hadrons_afterburner.empty()) {
      try {
        std::vector<Hadron> hadrons = any_cast<std::vector<Hadron>>(current_hadrons_afterburner);
        // convert to vector of shared pointers using std::transform
        std::vector<std::shared_ptr<Hadron>> shared_hadrons;
        std::transform(hadrons.begin(), hadrons.end(), std::back_inserter(shared_hadrons),
                    [](const Hadron& h) { return std::make_shared<Hadron>(h); });
        // print the number of hadrons in the event
        JSWARN << "Number of hadrons in event in FinishPerEventTasks = " << shared_hadrons.size();
        
        // store the hadrons in the BDM_final_state_hadrons_ vector
        BDM_final_state_hadrons_.push_back(shared_hadrons);
      } catch (const linb::bad_any_cast& e) {
        JSWARN << "Failed to retrieve hadrons from module SMASH: " << e.what();
      }
    } else {
      JSWARN << "Failed to retrieve hadrons from module SMASH, invalid QueryHistory::GetHistoryFromModule() call";
    }
  }

  CloseHadronicTimeEvolutionFileIfNecessary();

  // Clear some vectors
  new_hadrons_for_timestep_.clear();
  remove_hadrons_for_timestep_.clear();
  store_source_term_hadrons_iso_tau_.clear();
  store_spectator_hadrons_iso_tau_.clear();
  store_hadrons_soft_particlization_.clear();

  //JP: Quick fix, to be discussed, similar to writer, clear is only called for active tasks, so call here directly ...
  ClearTask();
}

void BulkDynamicsManager::WriteTask(weak_ptr<JetScapeWriter> w) {
  JSWARN << "BDM hadron printout";
  auto f = w.lock();
  if (!f) {
    return;
  }
  if (SMASH_IC_attached_ && !hydro_Cartesian_) {
    // write BDM_final_state_hadrons_ hadrons for all events to file
    f->WriteComment("BDM final state hadrons");
    for (const auto &event : BDM_final_state_hadrons_) {
      int i = -1;
      for (const auto hadron : event) {
        //f->WriteWhiteSpace("[" + to_string(++i) + "] H");
        f->Write(hadron);
      }
    }
  }
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

    // Get the current hadrons from the initial condition or afterburner module
    // Check if SMASH IC or afterburner is active (it should be only one of them)
    bool SMASH_IC_active = false;
    bool SMASH_Afterburner_active = false;
    for(auto it : GetTaskList()) {
      auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
      if(dynamic_pointer_cast<SmashInitialConditionWrapper>(module)) {
        SMASH_IC_active = true;
      }
      if(dynamic_pointer_cast<Afterburner>(module) && (module->GetId() == "SMASH")) {
        SMASH_Afterburner_active = true;
      }
    }
    // If SMASH IC or afterburner is active, get the current hadrons
    // If both are active, we assume the afterburner is the one to use
    // If neither is active, we cannot get the bulk info, so return
    std::vector<Hadron> current_hadrons;
    if (SMASH_IC_active) {
      linb::any current_hadrons_IC = QueryHistory::Instance()->GetHistoryFromModule("SMASHInitialState");
      if (!current_hadrons_IC.empty()) {
        try {
          current_hadrons = any_cast<std::vector<Hadron>>(current_hadrons_IC);
        } catch (const linb::bad_any_cast& e) {
          JSWARN << "Failed to retrieve hadrons from module SMASHInitialState: " << e.what();
        }
      } else {
        JSWARN << "Failed to retrieve hadrons from module SMASHInitialState, invalid QueryHistory::GetHistoryFromModule() call";
      }
    } else if (SMASH_Afterburner_active) {
      linb::any current_hadrons_afterburner = QueryHistory::Instance()->GetHistoryFromModule("SMASH");
      if (!current_hadrons_afterburner.empty()) {
        try {
          current_hadrons = any_cast<std::vector<Hadron>>(current_hadrons_afterburner);
        } catch (const linb::bad_any_cast& e) {
          JSWARN << "Failed to retrieve hadrons from module SMASH: " << e.what();
        }
      } else {
        JSWARN << "Failed to retrieve hadrons from module SMASH, invalid QueryHistory::GetHistoryFromModule() call";
      }
    }

    // Use hadronic_emt_.GetBulkInfo(t,x,y,z,bulk_info_ptr,current_hadrons) to get the bulk info from the hadronic medium
    if (SMASH_IC_active) {
      hadronic_emt_.GetBulkInfo(t,x,y,z,bulk_info_ptr,current_hadrons);
    } else if (SMASH_Afterburner_active) {
      hadronic_emt_.GetBulkInfo(t,x,y,z,bulk_info_ptr,current_hadrons);
    } else {
      JSWARN << "No SMASH IC or Afterburner active, cannot get bulk info!";
      return;
    }
    InfoWrapper(fluid_cell_info_ptr,bulk_info_ptr);
  }
}

void BulkDynamicsManager::InfoWrapper(
  std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr,
  std::unique_ptr<BulkMediaInfo> &bulk_info_ptr){
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
  // T^{\mu\nu} from bulk info not converted as not present in fluid cell info
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

void BulkDynamicsManager::DetermineHadronsCrossingIsoTau( 
                          std::vector<shared_ptr<Hadron>> &current_hadrons) {
  int i = 0;
  for (const auto& had : current_hadrons) {
    const FourVector r = had->x_in();
    const double t = r.t();
    const double z = r.z();
    const double tau = sqrt(t*t - z*z);

    if(tau >= IC_particle_extraction_tau_) {
      i++;
      RemoveHadron(had);
    }
  }
  //JSINFO << "Found " << i << " hadrons to remove from SMASH";
  //JSINFO << "Size remove_hadrons_for_timestep_ = " << remove_hadrons_for_timestep_.size();
}

void BulkDynamicsManager::ExtractHadronsFromTransportInitialConditionIsoTau(bool &AllHadronsCrossedIsoTau) {
  linb::any current_hadrons_IC = QueryHistory::Instance()->GetHistoryFromModule("SMASHInitialState");
  std::vector<Hadron> hadrons = any_cast<std::vector<Hadron>>(current_hadrons_IC);

  // Convert to vector of shared pointers using std::transform
  std::vector<std::shared_ptr<Hadron>> shared_hadrons;
  std::transform(hadrons.begin(), hadrons.end(), std::back_inserter(shared_hadrons),
              [](const Hadron& h) { return std::make_shared<Hadron>(h); });

  // Determine which particles should be removed from the SMASH initial condition
  DetermineHadronsCrossingIsoTau(shared_hadrons);

  for(const auto& hadron : remove_hadrons_for_timestep_) {
    bool participant = hadron->participant();
    bool add_hadron_to_source_term = false;
    // Check if the participant is in the kinematic cuts (if applied)
    if (participant) {
      bool hadron_above_pT_cut_threshold = false;
      if (enforce_pT_cut_) {
        // Check if the hadron has a pT larger than the cut
        const FourVector p = hadron->p_in();
        const double pT = sqrt(p.x()*p.x() + p.y()*p.y());
        if (pT > pT_cut_) {
          hadron_above_pT_cut_threshold = true;
        }
      }
      bool hadron_above_rapidity_cut_threshold = false;
      if (enforce_rapidity_cut_) {
        // Check if the hadron has a rapidity larger than the cut
        const FourVector p = hadron->p_in();
        const double rapidity = 0.5*log((p.t()+p.z())/(p.t()-p.z()));
        if (abs(rapidity) > rapidity_cut_) {
          hadron_above_rapidity_cut_threshold = true;
        }
      }
      // create a bool if the hadron is outside one of the kinematic cuts
      bool hadron_outside_cut_threshold = false;
      if (enforce_pT_cut_ && hadron_above_pT_cut_threshold) {
        hadron_outside_cut_threshold = true;
      }
      if (enforce_rapidity_cut_ && hadron_above_rapidity_cut_threshold) {
        hadron_outside_cut_threshold = true;
      }
      // If the hadron is outside the cuts, then it is a 'spectator'
      // and should not be added to the source term
      if (hadron_outside_cut_threshold) {
        add_hadron_to_source_term = false;
      } else {
        add_hadron_to_source_term = true;
      }
    } else {
      // Not a participant, so it is a spectator
      add_hadron_to_source_term = false;
    }

    if (add_hadron_to_source_term) {
      store_source_term_hadrons_iso_tau_.push_back(hadron);
    } else {
      store_spectator_hadrons_iso_tau_.push_back(hadron);
    }
  }

  if(shared_hadrons.empty()) {
    AllHadronsCrossedIsoTau = true;
  }
}

void BulkDynamicsManager::CreateHadronicSourceTermsForHydroInitializationIsoTau() {
  if (!weak_ptr_is_uninitialized(hadronic_liquefier_ptr_)) {
    //JSWARN << "Create a hydro source with " << store_source_term_hadrons_iso_tau_.size() << " hadrons";
    // Convert std::shared_ptr<Jetscape::Hadron> to raw pointers and store in a vector
    std::vector<Jetscape::Hadron> hadron_objects;
    for (const auto& ptr : store_source_term_hadrons_iso_tau_) {
      hadron_objects.push_back(*ptr);
    }
    hadronic_liquefier_ptr_.lock()->add_hydro_sources_hadrons(hadron_objects);
  } else {
    JSWARN << "Hadronic liquefier in BDM not initialized";
    exit(1);
  }
}

void BulkDynamicsManager::StoreHadronsFromSoftParticlization() {
  for(auto it : GetTaskList()) {
    auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
    if (auto hadronization = std::dynamic_pointer_cast<SoftParticlization>(module)) {
      std::vector<std::vector<shared_ptr<Hadron>>> hadron_list = hadronization->Hadron_list_;
      // Get the hadrons from the last event in the hadron list add them to 
      // store_hadrons_soft_particlization_, clean the Hadron_list_
      if (!hadron_list.empty()) {
        for (const auto& hadrons : hadron_list.back()) {
          store_hadrons_soft_particlization_.push_back(hadrons);
        }
        hadronization->ClearHadronList();
      }
    }
  }
}

void BulkDynamicsManager::PrintHadronicTimeEvolutionToFileIfNecessary() {
  // Get all the hadrons from SMASH (initial condition of afterburner) and print them into the file
  if (hadronic_time_evolution_to_file_) {
    double current_time = GetMainClock()->GetCurrentTime();
    // check if the SMASH IC is attached and get the hadrons from there
    // print the initial state hadrons only every 0.1 fm/c
    bool print_IC_hadrons = false;
    if (fmod(abs(current_time), 0.1) < 1e-3) {
      print_IC_hadrons = true;
    }
    if (SMASH_IC_in_progress_ && print_IC_hadrons) {
      linb::any current_hadrons_IC = QueryHistory::Instance()->GetHistoryFromModule("SMASHInitialState");
      std::vector<Hadron> hadrons = any_cast<std::vector<Hadron>>(current_hadrons_IC);

      // Convert to vector of shared pointers using std::transform
      std::vector<std::shared_ptr<Hadron>> shared_hadrons;
      std::transform(hadrons.begin(), hadrons.end(), std::back_inserter(shared_hadrons),
                  [](const Hadron& h) { return std::make_shared<Hadron>(h); });

      // Print the hadrons to the file
      for (const auto& had : shared_hadrons) {
        const FourVector r = had->x_in();
        const double t = r.t();
        const double x = r.x();
        const double y = r.y();
        const double z = r.z();
        const FourVector p = had->p_in();
        const double E = p.t();
        const double px = p.x();
        const double py = p.y();
        const double pz = p.z();
        *hadronic_time_evolution_file_ << current_time << " " << x << " " << y 
          << " " << z << " " << had->pid() << " " << E << " " << px << " " 
          << py << " " << pz << " " << had->participant() << endl;
      }
    }
    if (afterburner_in_progress_) {
      linb::any current_hadrons_afterburner = QueryHistory::Instance()->GetHistoryFromModule("SMASH");
      std::vector<Hadron> hadrons = any_cast<std::vector<Hadron>>(current_hadrons_afterburner);

      // Convert to vector of shared pointers using std::transform
      std::vector<std::shared_ptr<Hadron>> shared_hadrons;
      std::transform(hadrons.begin(), hadrons.end(), std::back_inserter(shared_hadrons),
                  [](const Hadron& h) { return std::make_shared<Hadron>(h); });

      // Print the hadrons to the file
      for (const auto& had : shared_hadrons) {
        const FourVector r = had->x_in();
        const double t = r.t();
        const double x = r.x();
        const double y = r.y();
        const double z = r.z();
        const FourVector p = had->p_in();
        const double E = p.t();
        const double px = p.x();
        const double py = p.y();
        const double pz = p.z();
        // participant status -1 indicates here that this hadron is from the
        // afterburner phase (last column in the file)
        *hadronic_time_evolution_file_ << current_time << " " << x << " " << y 
          << " " << z << " " << had->pid() << " " << E << " " << px << " " 
          << py << " " << pz << " " << -1 << endl;
      }
    }
  }
}

} // end namespace Jetscape
