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
#include "JetScapeWriterRootHepMC.h"
#endif

// User modules derived from jetscape framework clasess
#include "TrentoInitial.h"
#include "AdSCFT.h"
#include "Matter.h"
#include "LBT.h"
#include "Martini.h"
#include "Brick.h"
#include "BrickTest.h"
#include "GubserHydro.h"
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

// Forward declaration
void Show();

// -------------------------------------
// QA Task Demo with ROOT Hisograms etc
#include <Riostream.h>
#include "TRandom.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

class DemoQA : public JetScapeModuleBase
{
  public:

  DemoQA() : JetScapeModuleBase() {SetId("DemoQA"); oName = "";}
  DemoQA(string m_oName) : JetScapeModuleBase() {SetId("DemoQA"); oName = m_oName;}
  virtual ~DemoQA() {f->cd();hPt->Write("hPt");f->ls();f->Write();f->Close();}

  void Init() {f=new TFile(oName.c_str(),"RECREATE"); hPt = new TH1D("hPt","",100,0,100);}
  void Exec() {
    QueryHistory::Instance()->UpdateTaskMap();
    vector<any> eLossHistories = QueryHistory::Instance()->GetHistoryFromModules("JetEnergyLoss");
    for (auto mHist : eLossHistories)
    {
      auto ps = any_cast<std::shared_ptr<PartonShower>>(mHist);
      //For demonstation purpose only ...
      cout<<"Shower iniating parton pT = "<<ps->GetPartonAt(0)->pt()<<endl;

      hPt->Fill(ps->GetPartonAt(0)->pt());
    }
    //cout<<hPt->GetEntries()<<endl;
  }
  //Bummer, finish actually not yet recursively implemented !!!
  void Finish() {cout<<"Finish called ..."<<endl; f->cd();hPt->Write("hPt");f->ls();f->Write();f->Close();}

  private:

  string oName;
  TFile *f;
  TH1D *hPt;

};

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

  auto jetscape = make_shared<JetScape>();
  jetscape->SetXMLMainFileName("../config/jetscape_main.xml");
  jetscape->SetXMLUserFileName("../config/jetscape_user_root_test.xml");
  jetscape->SetId("primary");

  // Initial conditions and hydro
  //auto trento = make_shared<TrentoInitial>();
  auto trento = make_shared<InitialState>();
  auto pythiaGun= make_shared<PythiaGun> ();
  auto isr = make_shared<InitialStateRadiationTest> ();
  auto hydro = make_shared<Brick> ();
  jetscape->Add(trento);
  jetscape->Add(pythiaGun);
  jetscape->Add(hydro);

  // Energy loss
  auto jlossmanager = make_shared<JetEnergyLossManager> ();
  auto jloss = make_shared<JetEnergyLoss> ();

  //Matter is added but not executed, need to implement the per time step execution in JetEnergyLoss::DoShower()...
  auto matter = make_shared<Matter> ();
  //auto lbt = make_shared<LBT> ();
  //auto martini = make_shared<Martini> ();
  //auto adscft = make_shared<AdSCFT> ();

  // Note: if you use Matter, it MUST come first (to set virtuality)
  jloss->Add(matter);
  jlossmanager->Add(jloss);
  jetscape->Add(jlossmanager);

  auto hadroMgr = make_shared<HadronizationManager> ();
  auto hadro = make_shared<Hadronization> ();
  auto colorless = make_shared<ColorlessHadronization> ();
  hadro->Add(colorless);
  hadroMgr->Add(hadro);
  jetscape->Add(hadroMgr);

  // Output
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  writer->SetId("AsciiWriter"); //for task search test ...
  jetscape->Add(writer);
  //auto writer= make_shared<JetScapeWriterRoot> ("test_out.root");
  //writer->SetId("RootWriter");
  //auto writergz= make_shared<JetScapeWriterAsciiGZ> ("test_out.dat.gz");
  //jetscape->Add(writergz);
  #ifdef USE_HEPMC
  //auto hepmcwriter= make_shared<JetScapeWriterHepMC> ("test_out.hepmc");
  //jetscape->Add(hepmcwriter);
  auto hepmcwriterRoot= make_shared<JetScapeWriterRootHepMC> ("test_out_hepmc.root");
  jetscape->Add(hepmcwriterRoot);

  auto writer2= make_shared<JetScapeWriterHepMC> ("test_out_hepmc.hepmc");
  writer2->SetId("hepMCWriter");
  jetscape->Add(writer2);
  #endif

  auto QATest = make_shared<DemoQA>("qa_test.root");
  jetscape->Add(QATest);

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


  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"-------------------------------------------";
  INFO_NICE<<"| ROOT Brick Test XSCAPE Framework ...   |";
  INFO_NICE<<"-------------------------------------------";
  INFO_NICE;
}
