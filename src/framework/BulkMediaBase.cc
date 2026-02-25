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
// This is a general basic class for bulk media

#include "./BulkMediaBase.h"
#include "./JetScapeSignalManager.h"

using namespace std;

namespace Jetscape {
/**
 * @brief Initialize the BulkMediaBase module.
 *
 * Ensures that the XML configuration is loaded by calling the base-class
 * initialization and then initializes this task and any registered
 * sub-tasks.
 */
void BulkMediaBase::Init() {
  // Makes sure that XML file with options and parameters is loaded
  JetScapeModuleBase::InitTask();
  JSINFO << "Initializing BulkMediaBase : " << GetId() << " ...";

  InitTask();
  InitTasks();
}

/**
 * @brief Execute the module's task.
 *
 * Default implementation emits a verbose log entry. Subclasses should
 * override to provide actual execution logic.
 */
void BulkMediaBase::ExecuteTask() {
  VERBOSE(2) << "BulkMediaBase running: " << GetId() << " ...";
}

/**
 * @brief Perform time-step calculations.
 *
 * Default implementation logs a verbose message and calls
 * CalculateTimeTask() to perform the concrete time-dependent work.
 */
void BulkMediaBase::CalculateTime() {
  VERBOSE(2) << "BulkMediaBase running for time: " << GetId() << " ...";
  CalculateTimeTask();
}

}  // end namespace Jetscape
