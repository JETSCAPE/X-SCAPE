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

#include "JetScapePerEvent.h"
#include "JetScapeWriter.h"
#include "FluidDynamics.h"
#include "InitialState.h"
#include "PreequilibriumDynamics.h"
#include "SoftParticlization.h"

#include "QueryHistory.h"

#include <cstdlib>

using namespace std;

namespace Jetscape {

JetScapePerEvent::JetScapePerEvent() : JetScape() { VERBOSE(8); }

JetScapePerEvent::~JetScapePerEvent() { VERBOSE(8); }

//________________________________________________________________
void JetScapePerEvent::ExecInit() {
  if (exec_initialized_)
    return;

  JSINFO << BOLDRED << "Run JetScape (per event) ...";
  JSINFO << BOLDRED << "Number of Events = " << GetNumberOfEvents();

  // Collect active writers and snapshot the original active flags, mirroring the
  // pre-loop setup of JetScape::Exec().
  for (auto it : GetTaskList()) {
    if (dynamic_pointer_cast<JetScapeWriter>(it)) {
      if (it->GetActive()) {
        vWriter.push_back(dynamic_pointer_cast<JetScapeWriter>(it));
      }
    }
    if (dynamic_pointer_cast<JetScapeModuleBase>(it) && it->GetActive()) {
      dynamic_pointer_cast<JetScapeModuleBase>(it)->CheckExec();
    }

    taskOrgActiveMap.emplace(it->GetTaskNumber(), it->GetActive());
  }

  // Shift the framework-global event counter to the configured starting event.
  // Default is 0, leaving behavior unchanged. Only the public IncrementCurrentEvent()
  // mutator is available, so advance up to the requested start.
  if (start_event_ < 0) {
    JSWARN << "Negative start event (" << start_event_ << ") requested; using 0.";
    start_event_ = 0;
  }
  if (GetCurrentEvent() > start_event_) {
    JSWARN << "Event counter (" << GetCurrentEvent()
           << ") is already past the requested start event (" << start_event_
           << "); cannot move it backwards.";
  }
  while (GetCurrentEvent() < start_event_)
    IncrementCurrentEvent();

  exec_initialized_ = true;
}

//________________________________________________________________
void JetScapePerEvent::ExecPerEvent() {
  if (!exec_initialized_)
    ExecInit();

  int i = GetCurrentEvent();

  if (i % n_events_printout == 0) {
    JSINFO << BOLDRED << "Run Event # = " << i;
  }
  VERBOSE(1) << BOLDRED << "Run Event # = " << i;
  JSDEBUG << "Found " << GetNumberOfTasks() << " Modules Execute them ... ";

  // Execute and run per time step for modules if implemented ...
  if (ClockUsed()) {
    VERBOSE(3) << "Main Clock Reset ...";

    // Do per event execution except for Hadronization and Afterburner (if not
    // timestepped) ... Set the proper pre per event active etc flags first ...
    SetPerEventExecFlags(true);
    ExecuteTasks();

    GetMainClock()->Reset();

    JetScapeModuleBase::InitPerEventTasks();

    // Quick and dirty to see all tasks ... make recursive if needed
    QueryHistory::Instance()->UpdateTaskMap();

    do {
      VERBOSE(3) << BOLDRED
                 << "Current Main Clock Time = "
                 << GetMainClock()->GetCurrentTime()
                 << " dT = " << GetMainClock()->GetDeltaT();

      JetScapeModuleBase::CalculateTimeTasks();
      JetScapeModuleBase::ExecTimeTasks();

    } while (GetMainClock()->Tick());

    // Follow up with per event execution of Hadronization and Afterburner (if
    // not timestepped) ... Set the proper pre per event active etc flags first.
    SetPerEventExecFlags(false);
    ExecuteTasks();

    // Reset per event flags to original state to allow ClearTasks etc to be
    // executed properly and as expected ...
    ResetPerEventExecFlags();
  } else {
    ExecuteTasks();
    // Quick and dirty to see all tasks ... make recursive if needed
    QueryHistory::Instance()->UpdateTaskMap();
  }

  // Reusal of hydro events: deactivate task after it has finished but before it
  // gets cleaned up.
  if (reuse_hydro_) {
    if (n_reuse_hydro_ <= 0) {
      JSWARN << " reuse_hydro is set, but n_reuse_hydro = " << n_reuse_hydro_;
      throw std::runtime_error("Incompatible reusal settings.");
    }
    // Check if iMatter/ISR is used
    bool imatter_is_used = false;
    for (auto it : GetTaskList()) {
      if (it->GetId() == "PythiaGun") {
        for (auto itt : it->GetTaskList()) {
          if (itt->GetId() == "IsrManager") {
            VERBOSE(1) << " iMatter is used with reuse_hydro,"
                       << " so initial state is rerun for each event.";
            imatter_is_used = true;
            break;
          }
        }
      }
    }
    bool hydro_pointer_is_set = false;
    bool iss_pointer_is_set = false;
    for (auto it : GetTaskList()) {
      if (!dynamic_pointer_cast<FluidDynamics>(it) &&
          !dynamic_pointer_cast<PreequilibriumDynamics>(it) &&
          !dynamic_pointer_cast<InitialState>(it) &&
          !dynamic_pointer_cast<SoftParticlization>(it)) {
        continue;
      }

      // IS: For ISR+3DGlauber, the initial state 3D Glauber is not used
      // if imatter is not used, then initial state is rerun
      // This behavior must be rethaught for Au+Au
      // where we would expect the InitialState to be run per hydro event
      if (imatter_is_used && dynamic_pointer_cast<InitialState>(it)) {
        continue;
      }

      if (dynamic_pointer_cast<FluidDynamics>(it))
        if (dynamic_pointer_cast<FluidDynamics>(it)->IsTimeStepped()) {
          JSWARN << " Reusing hydro with per time stepped = true not allowed!";
          throw std::runtime_error(
              "Reusing hydro with per time stepped = true not allowed.");
        }

      // only deactivate the first hydro
      if (dynamic_pointer_cast<FluidDynamics>(it) && hydro_pointer_is_set) {
        continue;
      }

      if (i % n_reuse_hydro_ == n_reuse_hydro_ - 1) {
        JSDEBUG << " i was " << i
                << " i%n_reuse_hydro_ = " << i % n_reuse_hydro_
                << " --> ACTIVATING";
        it->SetActive(true);
        if (dynamic_pointer_cast<FluidDynamics>(it)) {
          hydro_pointer_is_set = true;
        }
      } else {
        JSDEBUG << " i was " << i
                << " i%n_reuse_hydro_ = " << i % n_reuse_hydro_
                << " --> DE-ACTIVATING";
        it->SetActive(false);
        if (dynamic_pointer_cast<FluidDynamics>(it)) {
          hydro_pointer_is_set = true;
        }
      }
      // Do the soft hadronization only at once
      if (dynamic_pointer_cast<SoftParticlization>(it))
        if (dynamic_pointer_cast<SoftParticlization>(it)->IsTimeStepped()) {
          JSWARN << " Reusing hydro with per time stepped = true not allowed!";
          throw std::runtime_error(
              "Reusing hydro with per time stepped = true not allowed.");
        }

      // only deactivate the first iSS
      if (dynamic_pointer_cast<SoftParticlization>(it) && iss_pointer_is_set) {
        continue;
      }

      if (i % n_reuse_hydro_ == n_reuse_hydro_ - 1) {
        JSDEBUG << " i was " << i
                << " i%n_reuse_hydro_ = " << i % n_reuse_hydro_
                << " --> ACTIVATING";
        it->SetActive(true);
        if (dynamic_pointer_cast<SoftParticlization>(it)) {
          iss_pointer_is_set = true;
        }
      } else {
        JSDEBUG << " i was " << i
                << " i%n_reuse_hydro_ = " << i % n_reuse_hydro_
                << " --> DE-ACTIVATING";
        it->SetActive(false);
        if (dynamic_pointer_cast<SoftParticlization>(it)) {
          iss_pointer_is_set = true;
        }
      }
    }
  }

  if (ClockUsed()) {
    VERBOSE(3) << "Clock is used and FinishPerEventTasks is called!";
    JetScapeModuleBase::FinishPerEventTasks();
  }

  // print all tasks and if they are active or not
  for (auto it : GetTaskList()) {
    auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);
    if (module) {
      if (module->GetActive()) {
        VERBOSE(3) << "IsActive(true) = " << module->GetId();
      } else {
        VERBOSE(3) << "IsActive(false) = " << module->GetId();
      }
    }
  }

  VERBOSE(3) << "Start the writing process ...";
  // collect module header data
  for (auto w : vWriter) {
    auto f = w.lock();
    if (f) {
      JetScapeTask::CollectHeaders(w);
    }
  }
  // official header
  for (auto w : vWriter) {
    auto f = w.lock();
    if (f) {
      f->WriteHeaderToFile();
    }
  }

  // event data
  for (auto w : vWriter) {
    auto f = w.lock();
    if (f) {
      JetScapeTask::WriteTasks(w);
    }
  }

  // Finalize
  for (auto w : vWriter) {
    auto f = w.lock();
    if (f) {
      f->WriteEvent();
    }
  }

  // NOTE: ClearTasks() and IncrementCurrentEvent() are intentionally NOT called
  // here. The module results stay in memory so an external program can read them.
  // Call ClearPerEvent() once that is done.
  VERBOSE(3) << "End of Event " << i;
}

//________________________________________________________________
void JetScapePerEvent::ClearPerEvent() {
  // Now clean up, only affects active tasks
  VERBOSE(3) << "Clearing tasks ...";
  JetScapeModuleBase::ClearTasks();

  IncrementCurrentEvent();
}

//________________________________________________________________
void JetScapePerEvent::Exec() {
  JSWARN << "JetScapePerEvent::Exec() is disabled. Drive this class from an "
            "external loop using ExecInit() / ExecPerEvent(i) / ClearPerEvent() "
            "instead. EXIT!";
  exit(-1);
}

}  // end namespace Jetscape
