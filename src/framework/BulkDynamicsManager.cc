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

} // end namespace Jetscape
