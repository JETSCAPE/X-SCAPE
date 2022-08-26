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
  Clear();

  if (GetNumberOfTasks() > 0)
    EraseTaskLast();
}

void BulkDynamicsManager::Clear() {
  JSDEBUG << "BulkDynamicsManager clear() ...";

  int n = GetNumberOfTasks();
  for (int i = 1; i < n; i++)
    EraseTaskLast();

  // Clean Up not really working with iterators (see also above!!!) Some logic not clear for me.
  JetScapeSignalManager::Instance()->CleanUp();
  JetScapeTask::ClearTasks();

}

void BulkDynamicsManager::Init() {
  JSINFO << "Intialize BulkDynamicsManager ...";

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid bulk dynamics Manager modules found ...";
    exit(-1);
  }

  JSINFO << "Found " << GetNumberOfTasks()
         << " Bulk Dynamics Manager Tasks/Modules Initialize them ... ";

  //Critical temperature to switch from hydro to something else
  Tc = GetXMLElementDouble({"BDM", "Tc"});

  BulkDynamicsManager::InitTasks();

}
void BulkDynamicsManager::Exec() {
  VERBOSE(1) << "Run BulkDynamicsManager Manager ...";
  JSDEBUG << "Task Id = " << this_thread::get_id();

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid Bulk Dynamics Manager modules found ...";
    exit(-1);
  }

  JetScapeTask::ExecuteTasks();

  //
  VERBOSE(3) << " " << GetNumberOfTasks()
             << " Bulk Dynamics Manager Tasks/Modules finished.";
}

void BulkDynamicsManager::CalculateTime()
{
  VERBOSE(3) << "Calculate Bulk Dynamics Manager per timestep ... Current Time = "<<GetModuleCurrentTime();
  VERBOSE(3) << "Task Id = " << this_thread::get_id();

  JetScapeModuleBase::CalculateTimeTasks();
}

void BulkDynamicsManager::ExecTime()
{
  VERBOSE(3) << "Execute Bulk Dynamics Manager per timestep ... Current Time = "<<GetModuleCurrentTime()<<" Thread Id = "<<this_thread::get_id();
  VERBOSE(3) << "Task Id = " << this_thread::get_id();

  JetScapeModuleBase::ExecTimeTasks();
}

void BulkDynamicsManager::InitPerEvent()
{
  VERBOSE(3) << "InitPerEvent Bulk Dynamics Manager when used per timestep ...";
  VERBOSE(3) << "Task Id = " << this_thread::get_id();

  JetScapeModuleBase::InitPerEventTasks();
}

void BulkDynamicsManager::FinishPerEvent()
{
  VERBOSE(3) << "FinishPerEvent Bulk Dynamics Manager when used per timestep ...";
  VERBOSE(3) << "Task Id = " << this_thread::get_id();

  JetScapeModuleBase::FinishPerEventTasks();

  //JP: Quick fix, to be discussed, similar to writer, clear is only called for active tasks, so call here directly ...
  Clear();
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
	if(dynamic_pointer_cast<FluidDynamics>(it))
	  dynamic_pointer_cast<FluidDynamics>(it)->GetHydroInfo(t,x,y,z,fluid_cell_info_ptr);
      	else
	  GetBulkInfo(t,x,y,z,fluid_cell_info_ptr);
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
      if(fluid_cell_info_ptr->temperature > Tc) validHydro = true;
    }
  }
  //if validHydro = true, we are done; if not get info from other modules
  if(validHydro == false){
    std::unique_ptr<BulkMediaInfo> bulk_info_ptr;
    for (auto it : GetTaskList()) {
      if(dynamic_pointer_cast<BulkMediaBase>(it)){
	dynamic_pointer_cast<BulkMediaBase>(it)->GetBulkInfo(t,x,y,z,bulk_info_ptr);
      }
    }
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
}

std::vector<shared_ptr<Hadron>> BulkDynamicsManager::GetNewHadronsAndClear() {
  std::vector<shared_ptr<Hadron>> new_h_to_return;
  // The swap puts the empty vector for new_hadrons_for_timestep
  // and therefore clears the vector to be filled again at next timestep
  new_h_to_return.swap(new_hadrons_for_timestep);
  return new_h_to_return;
}

} // end namespace Jetscape
