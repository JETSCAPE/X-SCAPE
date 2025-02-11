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
// This is a general basic class for hadronic transport

#include "./Transport.h"
#include "./JetScapeSignalManager.h"

using namespace std;

namespace Jetscape {

void Transport::Init() {
  // Makes sure that XML file with options and parameters is loaded
  JetScapeModuleBase::InitTask();
  JSINFO << "Initializing Hadronic Transport : " << GetId() << " ...";
  InitTask();
  InitTasks();
}

void Transport::ExecuteTask() {
  VERBOSE(2) << "Hadronic transport running: " << GetId() << " ...";
}

void Transport::CalculateTime() {
  VERBOSE(2) << "Hadronic transport running for time: " << GetId() << " ...";
  CalculateTimeTask();
}

}