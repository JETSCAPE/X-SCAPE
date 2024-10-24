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
// -------------------------------------------------
// XSCAPE Framework Clock Pythia Brick Test Program
// -------------------------------------------------

#include <iostream>
#include <time.h>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif


// User modules derived from jetscape framework clasess
#include "TrentoInitial.h"
#include "AdSCFT.h"
#include "Matter.h"
#include "LBT.h"
#include "Martini.h"
#include "NullPreDynamics.h"
#include "MusicWrapper.h"
#include "iSpectraSamplerWrapper.h"
#include "PythiaGun.h"
#include "InitialStateRadiationTest.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"
#include "ColorlessHadronization.h"
#include "CascadeTest.h"

#include "MainClock.h"
#include "ModuleClock.h"
#include "MilneClock.h"

#include "QueryHistory.h"

#include <chrono>
#include <thread>

using namespace std;

using namespace Jetscape;

// ------------------------------------------
//JP: Quick test module for access to history ...
//To be done currently hardcoded for a quick test in JetScape.cc ...
class HistTest : public JetScapeModuleBase
{
  public:

  HistTest() : JetScapeModuleBase() {SetId("HistTest");}

  //virtual void InitPerEvent() {QueryHistory::Instance()->PrintTaskMap();}

  virtual void ExecTime()
  {

    if (GetMainClock()->GetCurrentTime()<GetMainClock()->GetDeltaT()) QueryHistory::Instance()->PrintTaskMap();

    vector<any> eLossHistories = QueryHistory::Instance()->GetHistoryFromModules("JetEnergyLoss");

    if (GetMainClock()->GetCurrentTime()<2) {

    cout<< "HistTest::ExecTime(): Current Main Clock Time = "<<GetMainClock()->GetCurrentTime() << endl;
    cout<< "HistTest::ExecTime(): Print Histories via vector<any> eLossHistories = QueryHistory::Instance()->GetHistoryFromModules(\"JetEnergyLoss\")" <<endl;

    for (auto mHist : eLossHistories)
    {
      any_cast<std::shared_ptr<PartonShower>>(mHist)->PrintEdges(false);
    }
   }
  }

  private:
};
// ------------------------------------------
class ClockTest : public JetScapeModuleBase
{
public:

  ClockTest() : JetScapeModuleBase() {SetId("ClockTest");}

  virtual void ExecTime(){
    cout<< "ClockTest::ExecTime():============================ " <<     endl;
    cout<< "ClockTest::ExecTime(): Current Main Clock Time = "<<GetMainClock()->GetCurrentTime() << endl;
    cout<< "ClockTest::ExecTime(): Current Main Clock Detla T = "<<GetMainClock()->GetDeltaT() << endl;
    if(UseModuleClock())cout<< "ClockTest::ExecTime(): Current Module Clock Time = "<<GetModuleClock()->GetCurrentTime() << endl;
    if(UseModuleClock())cout<< "ClockTest::ExecTime(): Current Module Clock Detla T = "<<GetModuleClock()->GetDeltaT() << endl;
    cout<< "ClockTest::ExecTime(): Current Module Start Time = "<< GetTStart() << endl;
    cout<< "ClockTest::ExecTime(): Current Module End Time = "<< GetTEnd() << endl;
    cout<< "ClockTest::ExecTime(): IsValidModuleTime?  "<< IsValidModuleTime() << endl;
    cout<< "ClockTest::ExecTime():============================ " <<     endl;
  }

private:
};
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
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  //SetVerboseLevel (9 a lot of additional debug output ...)
  //If you want to suppress it: use SetVerboseLevle(0) or max  SetVerboseLevle(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);


  Show();

  // -------------
  //Test clock ...

  //auto mClock = make_shared<MainClock>();
  //mClock->SetTimeRefFrameId("SpaceTime");

  // clocks here are defaulted for testing, clocks can costumized via inhererting from the MainClock/ModuleClock base classes ...

  //auto mClock = make_shared<MainClock>("SpaceTime",-1,5,0.1); // JP: make consistent with reading from XML in init phase ...
  //auto mClock = make_shared<MainClock>("SpaceTime",-0.1,0.1,0.1);
  auto mClock = make_shared<MainClock>("SpaceTime", 0.5,5.0,0.1);

  auto mModuleClock = make_shared<ModuleClock>();
  mModuleClock->SetTimeRefFrameId("SpaceTime * 2");
  auto mMilneClock = make_shared<MilneClock>();
  mMilneClock->setEtaMax(5.0);

  mClock->Info();
  mModuleClock->Info();
  mMilneClock->Info();

  /*
  mClock->Info();

  //while (mClock->Next()) {

  mClock->Tick();
  mClock->Info();

  mModuleClock->Transform(mClock);
  mModuleClock->Info();

  //};
  */
  // -------------

  auto jetscape = make_shared<JetScape>();
  jetscape->SetXMLMainFileName("../config/jetscape_main.xml");
  jetscape->SetXMLUserFileName("../config/jetscape_user_MUSICMainClockTest.xml");
  jetscape->SetId("primary");
  jetscape->AddMainClock(mClock);
  jetscape->ClockInfo();

  auto clockTest1 = make_shared<ClockTest>();
  clockTest1->SetTimeStepped(true);
  clockTest1->SetTimeRange(-1,0);

  auto clockTest2 = make_shared<ClockTest>();
  clockTest2->SetTimeStepped(true);
  //clockTest2->SetTimeRange(0,3.5);

  auto clockTest3 = make_shared<ClockTest>();
  clockTest3->SetTimeStepped(true);
  //clockTest3->SetTimeRange(0,4.);
  clockTest3->AddModuleClock(mModuleClock);

  jetscape->Add(clockTest1);
  jetscape->Add(clockTest2);
  jetscape->Add(clockTest3);

  // Initial conditions and hydro
  auto trento = make_shared<TrentoInitial>();
  auto pythiaGun= make_shared<PythiaGun> ();
  auto isr = make_shared<InitialStateRadiationTest> ();
  auto null_predynamics = make_shared<NullPreDynamics> ();
  auto hydro = make_shared<MpiMusic> ();
  //hydro->AddModuleClock(mMilneClock);
  hydro->SetTimeStepped(true);

  // surface sampler
  auto iSS = make_shared<iSpectraSamplerWrapper> ();
  iSS->SetTimeStepped(true);

  //auto hydroTest = make_shared<BrickTest> ();
  //hydroTest->SetMultiThread(true);
  //hydroTest->SetTimeStepped(true);

  jetscape->Add(trento);
  jetscape->Add(pythiaGun);
  //jetscape->Add(isr);
  jetscape->Add(null_predynamics);
  jetscape->Add(hydro);
  //jetscape->Add(hydroTest);
  jetscape->Add(iSS);

  // Energy loss
  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();

  //Do per time step for these modules with main clock attached ...
  //Needed to overwrite functions: CalculateTime() and ExecTime(), in these functions get
  //time, either main clock time or if module clock attached the tranformed time via: GetModuleCurrentTime();

  jlossmanager->SetTimeStepped(true);
  jloss->SetTimeStepped(true);

  // To test for time step consistency execution settings uncomment next line ...
  //jloss->SetTimeStepped(false);

  //quick and dirty to check if module clock transformation is working conceptually ...
  //jloss->AddModuleClock(mModuleClock);

  //Matter is added but not executed, need to implement the per time step execution in JetEnergyLoss::DoShower()...
  auto matter = make_shared<Matter> ();
  // auto lbt = make_shared<LBT> ();
  //auto martini = make_shared<Martini> ();
  //auto adscft = make_shared<AdSCFT> ();

  // Note: if you use Matter, it MUST come first (to set virtuality)
  //jloss->Add(matter);
  // jloss->Add(lbt);  // go to 3rd party and ./get_lbtTab before adding this module
  // jloss->Add(martini);
  //jloss->Add(adscft);
  //jlossmanager->Add(jloss);
  //jetscape->Add(jlossmanager);

  auto cascadeTest = make_shared<CascadeTest> ();
  cascadeTest->SetMultiThread(true);
  cascadeTest->SetTimeStepped(true);
  jetscape->Add(cascadeTest);

  //Test task for access of History via QueryHistory instance and use any data-type for generic access via JetScapeModuleBase::GetHistory()
  auto histTest = make_shared<HistTest>();
  histTest->SetTimeStepped(true);// to be executed per time step
  //jetscape->Add(histTest);

  // JP: Leave out for now for testing clock(s) ... has to be updated accordingly ...
  // (Hadronization not yet modified for per-timestep evolution see JetEnergyLossManager as an example ...)
  // Hadronization
  /*
  auto hadroMgr = make_shared<HadronizationManager> ();
  auto hadro = make_shared<Hadronization> ();
  //auto hadroModule = make_shared<ColoredHadronization> ();
  //hadro->Add(hadroModule);
  auto colorless = make_shared<ColorlessHadronization> ();
  hadro->Add(colorless);
  hadroMgr->Add(hadro);
  jetscape->Add(hadroMgr);
  */

  // Output
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  writer->SetId("AsciiWriter"); //for task search test ...
  jetscape->Add(writer);

  /*
#ifdef USE_GZIP
  // same as JetScapeWriterAscii but gzipped
  auto writergz= make_shared<JetScapeWriterAsciiGZ> ("test_out.dat.gz");
  jetscape->Add(writergz);
#endif
  // HEPMC3
#ifdef USE_HEPMC
  auto hepmcwriter= make_shared<JetScapeWriterHepMC> ("test_out.hepmc");
  jetscape->Add(hepmcwriter);
#endif
  */

  //test ...
  //QueryHistory::Instance()->AddMainTask(jetscape);
  //QueryHistory::Instance()->PrintTasks();
  //QueryHistory::Instance()->PrintTaskMap();

  //check with quick and dirty ... make recursive ...
  /*
  cout<<jetscape->GetNumberOfTasks()<<endl;
  auto taskList = jetscape->GetTaskList();
  for (auto it : taskList) {
    cout<<it->GetId()<<endl;
    for (auto it2 : it->GetTaskList()) {
      cout<<" "<<it2->GetId()<<endl;
      for (auto it3 : it2->GetTaskList())
        cout<<"  "<<it3->GetId()<<endl;}
  }
  */

  //printAllTasks(taskList);

  // Initialize all modules tasks
  jetscape->Init();

  // Run JetScape with all task/modules as specified
  jetscape->Exec();

  // For the future, cleanup is mostly already done in write and clear
  jetscape->Finish();

  INFO_NICE<<"Finished!";
  cout<<endl;

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));

  // Print pythia statistics
  // pythiaGun->stat();

  // // Demonstrate how to work with pythia statistics
  // //Pythia8::Info& info = pythiaGun->info;
  // cout << " nTried    = " << info.nTried() << endl;
  // cout << " nSelected = " << info.nSelected()  << endl;
  // cout << " nAccepted = " << info.nAccepted()  << endl;
  // cout << " sigmaGen  = " <<   info.sigmaGen()  << endl;
  // cout << " sigmaErr  = " <<   info.sigmaErr()  << endl;

  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"-------------------------------------------";
  INFO_NICE<<"| Clock Brick Test XSCAPE Framework ...   |";
  INFO_NICE<<"-------------------------------------------";
  INFO_NICE;
}
