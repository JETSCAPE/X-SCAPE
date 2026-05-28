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
// ------------------------------------------------------------
// JetScape per-event driver example: the event loop lives here in the
// macro (not inside Exec), so the module results of each event can be
// inspected before the per-event memory is released.
// -------------------------------------------------------------

#include <iostream>
#include <time.h>

// JetScape Framework includes ...
#include "JetScapePerEvent.h"
#include "JetScapeWriterStream.h"
#include "QueryHistory.h"
#include "Version.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

#include <chrono>
#include <thread>

using namespace Jetscape;

// -------------------------------------

int main(int argc, char** argv) {
  clock_t t;
  t = clock();
  time_t start, end;
  time(&start);

  // Logger settings (can be also set also via XML file, although note in that
  // case they will apply only after they are initialized)
  JetScapeLogger::Instance()->SetInfo(true);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  JetScapeLogger::Instance()->SetVerboseLevel(0);

  // Create the per-event Jetscape task, and assign XML configuration files from
  // command line arguments. The user can supply 0, 1, 2 arguments, where the
  // first (second) corresponds to the user (main) XML path.
  auto jetscape = make_shared<JetScapePerEvent>();
  const char* mainXMLName = "../config/jetscape_main.xml";
  const char* userXMLName = "../config/jetscape_user.xml";
  if (argc == 2) {
    if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
      std::cout << "Command line options:" << std::endl;
      std::cout << "    First (optional) argument: path to user XML file       "
                   "  ./JetScapePerEventTest /path/to/user.xml"
                << std::endl;
      std::cout << "    Second (optional) argument: path to main XML file      "
                   "./JetScapePerEventTest /path/to/user.xml /path/to/main.xml"
                << std::endl;
      std::cout << "    If no command line options are given, defaults are "
                   "used: config/jetscape_user.xml config/jetscape_main.xml"
                << std::endl;
      return -1;
    } else if (strcmp(argv[1], "--version") == 0 ||
               strcmp(argv[1], "-v") == 0) {
      std::cout << " XSCAPE version = " << XscapeVersion
                << " (includes JETSCAPE version = " << JetScapeVersion << ")"
                << std::endl;
      return -1;
    } else {
      userXMLName = argv[1];
    }
  } else if (argc == 3) {
    userXMLName = argv[1];
    mainXMLName = argv[2];
  }
  jetscape->SetXMLMainFileName(mainXMLName);
  jetscape->SetXMLUserFileName(userXMLName);

  // Initialize all modules tasks
  jetscape->Init();

  // Optionally start the event counter at a non-zero value (default 0).
  //jetscape->SetStartEvent(100);

  // One-time pre-loop setup (writers, CheckExec, active-flag snapshot).
  // Optional: ExecPerEvent() also calls this lazily on the first event.
  jetscape->ExecInit();

  // External event loop. The execution of one event is steered from here,
  // rather than from the internal loop in JetScape::Exec(). The event number
  // is tracked by the global counter (GetCurrentEvent()), advanced by
  // ClearPerEvent().
  int nEvents = jetscape->GetNumberOfEvents();
  for (int i = 0; i < nEvents; i++) {
    // Run a single event. All module results stay in memory afterwards.
    jetscape->ExecPerEvent();

    // --- External access point ---
    // At this point the per-event results are still live and can be read by
    // this program before the memory is released, e.g. via the QueryHistory
    // singleton or via the module task list. As a minimal demonstration we
    // just report the current event number here; replace this with your own
    // analysis of the module histories.
    INFO_NICE << "External access for event "
              << JetScapeModuleBase::GetCurrentEvent()
              << ": module data is available here.";

    // Release the per-event memory and advance the event counter.
    jetscape->ClearPerEvent();
  }

  // For the future, cleanup is mostly already done in write and clear
  jetscape->Finish();

  INFO_NICE << "Finished!";
  cout << endl;

  t = clock() - t;
  time(&end);
  printf("CPU time: %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
  printf("Real time: %f seconds.\n", difftime(end, start));
  return 0;
}
