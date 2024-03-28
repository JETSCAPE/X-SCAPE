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

// -----------------------------------------------------------------------------
// Test the generation of SMASH initial conditions in the collider modus
// -----------------------------------------------------------------------------

#include <iostream>
#include <time.h>
#include <chrono>
#include <thread>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetScapeLogger.h"
#include "JetScapeWriterStream.h"
#include "JetScapeXML.h"

//#include "HadronicLiquefier.h"
#include "BulkDynamicsManager.h"
#include "SMASHInitialStateWrapper.h"
#include "SmashWrapper.h"
#include "MusicWrapper.h"
#include "iSpectraSamplerWrapper.h"

using namespace Jetscape;

// Forward declaration
void Show();

int main(int argc, char** argv)
{
  clock_t t; t = clock();
  time_t start, end; time(&start);

  cout<<endl;

  // DEBUG=true by default and REMARK=false
  // can be also set also via XML file (at least partially)
  JetScapeLogger::Instance()->SetInfo(true);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  //SetVerboseLevel (9 a lot of additional debug output ...)
  //If you want to suppress it: use SetVerboseLevel(0) or max  SetVerboseLevel(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);

  Show();

  // clocks here are defaulted for testing, clocks can customized via inheriting from the MainClock/ModuleClock base classes ...
  auto mClock = make_shared<MainClock>("SpaceTime",-2,3,0.1); // JP: make consistent with reading from XML in init phase ...
  mClock->Info();

  auto jetscape = make_shared<JetScape>();
  const char* mainXMLName = "../config/jetscape_main.xml";
  const char* userXMLName = "../config/jetscape_user_SMASHInitialConditionTest.xml";

  jetscape->SetXMLMainFileName(mainXMLName);
  jetscape->SetXMLUserFileName(userXMLName);
  JetScapeXML::Instance()->OpenXMLMainFile(mainXMLName);
  JetScapeXML::Instance()->OpenXMLUserFile(userXMLName);

  jetscape->AddMainClock(mClock);
  jetscape->ClockInfo();
  jetscape->SetTimeStepped(true);

  // Initial conditions
  auto initial_state = make_shared<InitialState>();
  jetscape->Add(initial_state);

  auto smash_ic = make_shared<SmashInitialConditionWrapper>();
  // per time step for the IC
  smash_ic->SetTimeStepped(true);

  // Liquefier
  //auto hadronic_liquefier = make_shared<HadronicLiquefier>();

  // Hydro evolution
  auto hydro = make_shared<MpiMusic>();
  // per time step for the hydro
  hydro->SetTimeStepped(true);
  //hydro->add_a_hadronic_liquefier(hadronic_liquefier);

  // Soft particlization
  //auto iSS = make_shared<iSpectraSamplerWrapper>();
  // per time step soft particlization
  //iSS->SetTimeStepped(true);

  // Bulk Dynamics Manager (BDM)
  auto bdm = make_shared<BulkDynamicsManager>();
  bdm->SetTimeStepped(true);
  bdm->Add(smash_ic);
  //bdm->Add(hydro);
  //bdm->Add(iSS);

  // Add BDM to X-SCAPE
  jetscape->Add(bdm);

  // Output
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  writer->SetId("Writer");
  jetscape->Add(writer);

  // Initialize all modules tasks
  jetscape->Init();

  // Run JetScape with all task/modules as specified ...
  jetscape->Exec();

  // "dummy" so far ...
  // Most things done in write and clear ...
  jetscape->Finish();

  INFO_NICE<<"Finished!";
  cout<<endl;

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"-------------------------------------------------------";
  INFO_NICE<<"| SMASH Initial Condition Test X-SCAPE Framework ...  |";
  INFO_NICE<<"-------------------------------------------------------";
  INFO_NICE;
}