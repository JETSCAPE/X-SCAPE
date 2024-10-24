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

#include "JetScapeModuleBase.h"
#include "JetScapeXML.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeLogger.h"

#include <iostream>

namespace Jetscape {

// Create an instance of the static map to register modules
JetScapeModuleFactory::map_type *JetScapeModuleFactory::moduleMap =
    new JetScapeModuleFactory::map_type;

int JetScapeModuleBase::current_event = 0;

// ---------------------------------------------------------------------------
/** Default constructor to create a JetScapeModuleBase. It sets the XML file name to a default string value.
   */
JetScapeModuleBase::JetScapeModuleBase()
    : JetScapeTask(), xml_main_file_name(""), xml_user_file_name(""), time_stepped(false),
      mt19937_generator_(nullptr), TimeModule(0.,100.) {}

// ---------------------------------------------------------------------------
/** This is a destructor for the JetScapeModuleBase.
   */
JetScapeModuleBase::~JetScapeModuleBase() { disconnect_all(); }

// ---------------------------------------------------------------------------
/** A virtual function for a default initialization of a JetScapeModuleBase. It also checks whether a XML file is loaded or not.
   */
void JetScapeModuleBase::InitTask() {
  if (!JetScapeXML::Instance()->GetXMLRootMain()) {
    JSWARN << "Not a valid JetScape Main XML file or no XML file loaded!";
    exit(-1);
  }
  if (!JetScapeXML::Instance()->GetXMLRootUser()) {
    JSWARN << "Not a valid JetScape XML file or no XML file loaded!";
    exit(-1);
  }
}

// ---------------------------------------------------------------------------
/** This function returns a random number based on Mersenne-Twister algorithm.
   */
shared_ptr<std::mt19937> JetScapeModuleBase::GetMt19937Generator() {
  // Instantiate if it isn't there yet
  if (!mt19937_generator_) {
    mt19937_generator_ =
        JetScapeTaskSupport::Instance()->GetMt19937Generator(GetMyTaskNumber());
  }
  return mt19937_generator_;
}


void JetScapeModuleBase::ExecuteTasks()
{
  auto tasks =  GetTaskList();
  VERBOSE(7) << " : # Subtasks = " << tasks.size();
  for (auto it : tasks) {
    auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
    if (module && module->GetActive() && !module->IsTimeStepped()) {
      JSDEBUG << "Executing " << it->GetId();
      it->Exec();
    }
  }
}

void JetScapeModuleBase::ClearTasks() {
  auto tasks =  GetTaskList();
  VERBOSE(7) << " : # Subtasks = " << tasks.size();
  for (auto it : tasks) {
    auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
    if (module && module->GetActive() && !module->IsTimeStepped()) {
      it->Clear();
    }
  }
}

void JetScapeModuleBase::CheckExec()
{
  VERBOSE(7) << "JetScapeModuleBase::CheckExec()";
  auto tasks =  GetTaskList();
  for (auto it : tasks) {
    auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
    if (module && module->GetActive()) {
      if (IsTimeStepped() != std::dynamic_pointer_cast<JetScapeModuleBase>(it)->IsTimeStepped()) {
        //if (!std::dynamic_pointer_cast<JetEnergyLoss>(it) && !std::dynamic_pointer_cast<JetEnergyLoss>()) {
          JSWARN<<" ERROR: "<<GetId() <<" and "<< it->GetId()<<" are not set consistently in time step mode. Proper execution can not be ensured. EXIT!"; exit(-1);
        //}
      }
    }
  }

  CheckExecs();
}

void JetScapeModuleBase::CheckExecs()
{
  auto tasks =  GetTaskList();
  VERBOSE(7) << " : # Subtasks = " << tasks.size();
  for (auto it : tasks) {
    auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
    if (module && module->GetActive()) {
      std::dynamic_pointer_cast<JetScapeModuleBase>(it)->CheckExec();
    }
  }
}

void JetScapeModuleBase::CalculateTimeTasks()
{
  if (ClockUsed()) {
    auto tasks =  GetTaskList();
    for (auto it : tasks) {
      auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
      if (module && module->IsTimeStepped() && module->IsValidModuleTime()) {
        VERBOSE(3) << "Calculate Time Step = " << it->GetId();
        module->CalculateTime();
      }
    }
 }
}

void JetScapeModuleBase::ExecTimeTasks()
{
  if (ClockUsed()) {
    auto tasks =  GetTaskList();
    for (auto it : tasks) {
      auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
      if (module && module->IsTimeStepped() && module->IsValidModuleTime()) {
  	     VERBOSE(3) << "Execute Time Step = " << it->GetId();
  	     module->ExecTime();
      }
    }
  }
}

void JetScapeModuleBase::InitPerEventTasks()
{
  if (ClockUsed()) {
    auto tasks =  GetTaskList();
    for (auto it : tasks) {
      auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
      if (module && module->IsTimeStepped()) {
         VERBOSE(3) << "InitPerEventTasks " << it->GetId();
  	     module->InitPerEvent();
      }
    }
  }
}

void JetScapeModuleBase::FinishPerEventTasks()
{
  if (ClockUsed()) {
    auto tasks =  GetTaskList();
    for (auto it : tasks) {
      auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
      if (module && module->IsTimeStepped()) {
         VERBOSE(3) << "FinishPerEventTasks " << it->GetId();
  	     module->FinishPerEvent();
      }
    }
  }
}

} // end namespace Jetscape
