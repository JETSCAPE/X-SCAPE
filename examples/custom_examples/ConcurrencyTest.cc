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
// ------------------------------------------------------------
// JetScape Framework hydro from file Test Program
// (use either shared library (need to add paths; see setup.csh)
// (or create static library and link in)
// -------------------------------------------------------------

#include <iostream>
#include <time.h>

// Move it here to avoid conflicts:
// 1. Conflicts with MUSIC macros: hbarc, theta, limit, etc
// 2. Conflict of make_unique from smash and from JetScape
#include "SmashWrapper.h"

// JetScape Framework includes ...
#include "BulkDynamicsManager.h"
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

// User modules derived from jetscape framework clasess
#include "AdSCFT.h"
#include "Matter.h"
#include "Martini.h"
#include "MusicWrapper.h"
// Make sure that nasty MUSIC macros are neutralized
#undef PI
#undef hbarc
#undef default_tol
#undef absol
#undef maxi
#undef mini
#undef sgn
#undef theta
#undef gmn
#undef limit

#include "FreestreamMilneWrapper.h"
#include "iSpectraSamplerWrapper.h"
#include "TrentoInitial.h"
#include "PGun.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColorlessHadronization.h"
#include "ColoredHadronization.h"
#include "HybridHadronization.h"

#include <chrono>
#include <thread>

using namespace std;

using namespace Jetscape;

// Forward declaration
void Show();

// -------------------------------------

int main(int argc, char** argv)
{
  clock_t t; t = clock();
  time_t start, end; time(&start);

  cout<<endl;

  // DEBUG=true by default and REMARK=false
  // can be also set also via XML file (at least partially)
  JetScapeLogger::Instance()->SetInfo(true);
  JetScapeLogger::Instance()->SetDebug(true);
  JetScapeLogger::Instance()->SetRemark(false);
  //SetVerboseLevel (9 a lot of additional debug output ...)
  //If you want to suppress it: use SetVerboseLevel(0) or max  SetVerboseLevel(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(8);

  Show();

  // clocks here are defaulted for testing, clocks can costumized via inhererting from the MainClock/ModuleClock base classes ...
  auto mClock = make_shared<MainClock>("SpaceTime",0,100,10.0); // JP: make consistent with reading from XML in init phase ...
  mClock->Info();

  auto jetscape = make_shared<JetScape>();

  const char* mainXMLName = "../config/jetscape_main.xml";
  const char* userXMLName = "../config/jetscape_user.xml";
  jetscape->SetXMLMainFileName(mainXMLName);
  jetscape->SetXMLUserFileName(userXMLName);  // Needs to be mostly empty

  jetscape->AddMainClock(mClock);
  jetscape->ClockInfo();

  jetscape->SetReuseHydro (false);
  jetscape->SetNReuseHydro (0);
  // jetscape->SetReuseHydro (true);
  // jetscape->SetNReuseHydro (5);


  // NOT YET CONCURRENT ////////////////////////////////////////////////////////

  // Initial conditions and hydro
  // auto trento = make_shared<TrentoInitial>();
  // auto freestream = make_shared<FreestreamMilneWrapper> ();
  // auto hydro = make_shared<MpiMusic> ();

  // auto pGun= make_shared<PGun> ();
  // jetscape->Add(pGun);

  // jetscape->Add(trento);
  // jetscape->Add(freestream);
  // jetscape->Add(hydro);

  // Energy loss
  // auto jlossmanager = make_shared<JetEnergyLossManager> ();
  // auto jloss = make_shared<JetEnergyLoss> ();
  //
  // auto matter = make_shared<Matter> ();
  // // auto lbt = make_shared<LBT> ();
  // // auto martini = make_shared<Martini> ();
  // // auto adscft = make_shared<AdSCFT> ();
  //
  // // Note: if you use Matter, it MUST come first (to set virtuality)
  // jloss->Add(matter);
  // // jloss->Add(lbt);  // go to 3rd party and ./get_lbtTab before adding this module
  // // jloss->Add(martini);
  // // jloss->Add(adscft);
  // jlossmanager->Add(jloss);
  // jetscape->Add(jlossmanager);

  // Hadronization
  // auto hadroMgr = make_shared<HadronizationManager> ();
  // auto hadro = make_shared<Hadronization> ();
  // auto hybrid = make_shared<HybridHadronization> ();
  // hadro->Add(hybrid);
  // auto colorless = make_shared<ColorlessHadronization> ();
  // // hadro->Add(colorless);
  // hadroMgr->Add(hadro);
  // jetscape->Add(hadroMgr);

  // surface sampler
  // auto iSS = make_shared<iSpectraSamplerWrapper>();
  // jetscape->Add(iSS);

  //////////////////////////////////////////////////////////////////////////////


  // afterburner through BDM
  auto bulkmanager = make_shared<BulkDynamicsManager> ();
  auto smash = make_shared<SmashWrapper>();
  bulkmanager->SetTimeStepped(true);
  smash->SetTimeStepped(true); // if attached modules are not set to timestepped there will be not executed
  bulkmanager->Add(smash);
  // bulkmanager->Add(cartestianhydro) // todo: attach hydro in cartesian coordinates here as well
  jetscape->Add(bulkmanager);

  // afterburner (w/o BDM)
  // auto smash = make_shared<SmashWrapper>();
  // smash->SetTimeStepped(true);
  // jetscape->Add(smash);

  // Output
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  writer->SetId("Writer");
  // same as JetScapeWriterAscii but gzipped
  // auto writer= make_shared<JetScapeWriterAsciiGZ> ("test_out.dat.gz");
  // HEPMC3
#ifdef USE_HEPMC
  // auto writer= make_shared<JetScapeWriterHepMC> ("test_out.hepmc");
#endif
  jetscape->Add(writer);

  // Intialize all modules tasks
  jetscape->Init();

  // Run JetScape with all task/modules as specified ...
  jetscape->Exec();

  // "dummy" so far ...
  // Most thinkgs done in write and clear ...
  jetscape->Finish();

  INFO_NICE<<"Finished!";
  cout<<endl;

  // wait for 5s
  //std::this_thread::sleep_for(std::chrono::milliseconds(500000));

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));
  //printf ("Real time: %f seconds.\n",(start-end));
  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"-----------------------------------------------";
  INFO_NICE<<"| SMASH Test JetScape Framework ... |";
  INFO_NICE<<"-----------------------------------------------";
  INFO_NICE;
}