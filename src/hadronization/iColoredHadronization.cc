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

#include "JetScapeSignalManager.h"
#include "iColoredHadronization.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include "tinyxml2.h"
#include "MCGlauberGenStringWrapper.h"
#include <memory>

#define DEBUG_ISMAIL_4

using namespace Jetscape;
using namespace Pythia8;

// Register the module with the base class
RegisterJetScapeModule<iColoredHadronization>
    iColoredHadronization::reg("iColoredHadronization");

Pythia8::Pythia iColoredHadronization::pythia("IntentionallyEmpty", false);

iColoredHadronization::iColoredHadronization() {
  SetId("ISRMyHadroTest");
  VERBOSE(8);
}

iColoredHadronization::~iColoredHadronization() { VERBOSE(8); }

void iColoredHadronization::Init() {

  std::string s = GetXMLElementText({"JetHadronization", "name"});
  JSDEBUG << s << " to be initializied ...";

  double p_read_xml =
      GetXMLElementDouble({"JetHadronization", "eCMforHadronization"});
  p_fake = p_read_xml;

  std::string weak_decays =
      GetXMLElementText({"JetHadronization", "weak_decays"});

  VERBOSE(2) << "Start Hadronizing using the PYTHIA module...";

  // Show initialization at DEBUG or high verbose level
  pythia.readString("Init:showProcesses = off");
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showMultipartonInteractions = off");
  pythia.readString("Init:showChangedParticleData = off");
  if (JetScapeLogger::Instance()->GetDebug() ||
      JetScapeLogger::Instance()->GetVerboseLevel() > 2) {
    pythia.readString("Init:showProcesses = on");
    pythia.readString("Init:showChangedSettings = on");
    pythia.readString("Init:showMultipartonInteractions = on");
    pythia.readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  if (JetScapeLogger::Instance()->GetDebug() ||
      JetScapeLogger::Instance()->GetVerboseLevel() > 2) {
    pythia.readString("Next:numberShowInfo = 1");
    pythia.readString("Next:numberShowProcess = 1");
    pythia.readString("Next:numberShowEvent = 1");
  }

  pythia.readString("ProcessLevel:all = off");
  pythia.readString("PartonLevel:FSR=off");
  if (weak_decays == "off") {
    JSINFO << "Weak decays are turned off";
    pythia.readString("HadronLevel:Decay = off");
  } else {
    JSINFO << "Weak decays are turned on";
    pythia.readString("HadronLevel:Decay = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = 10.0");
  }
  pythia.init();
}

void iColoredHadronization::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  auto f = w.lock();
  if (!f)
    return;
  f->WriteComment("Hadronization Module : " + GetId());
  f->WriteComment("Hadronization to be implemented accordingly ...");
}

void iColoredHadronization::DoHadronization(
    vector<vector<shared_ptr<Parton>>> &shower,
    vector<shared_ptr<Hadron>> &hOut, vector<shared_ptr<Parton>> &pOut) {
  
  // JSINFO << "Starting "
  Event &event = pythia.event;
  event.reset();
  double pz = p_fake;
  #ifdef DEBUG_ISMAIL_4
    std::ofstream File3;
    File3.open("ISR-FinalPartons.dat", std::ofstream::out);
    File3 << "## &&&&&&&&&&&&&&&&&&& the number of showers are: " << shower.size()  << " time " << GetModuleCurrentTime() << std::endl;
    File3 << "# event status label pid col acol max_col px py pz E" << std::endl;
  #endif
  // auto Particles = shower->GetFinalPartons();
  for (unsigned int ishower = 0; ishower < shower.size(); ++ishower) {
    JSDEBUG << "&&&&&&&&&&&&&&&&&&& there are " << shower.at(ishower).size()
            << " partons in the shower number " << ishower;
    for (unsigned int ipart = 0; ipart < shower.at(ishower).size(); ++ipart) {
      
      // if(shower.at(ishower).at(ipart)->pstat() < 0 ) continue;

      double onshellE = pow(pow(shower.at(ishower).at(ipart)->px(), 2) +
                                pow(shower.at(ishower).at(ipart)->py(), 2) +
                                pow(shower.at(ishower).at(ipart)->pz(), 2),
                            0.5);

      if (shower.at(ishower).at(ipart)->pid() == 22) {

        VERBOSE(1) << BOLDYELLOW
                   << " photon found in colored hadronization with ";
        VERBOSE(1) << BOLDYELLOW
                   << "px = " << shower.at(ishower).at(ipart)->px();
        //cin >> blurb;
      }

      event.append(shower.at(ishower).at(ipart)->pid(), 23,
                   shower.at(ishower).at(ipart)->color(),
                   shower.at(ishower).at(ipart)->anti_color(),
                   shower.at(ishower).at(ipart)->px(),
                   shower.at(ishower).at(ipart)->py(),
                   shower.at(ishower).at(ipart)->pz(), onshellE);
      #ifdef DEBUG_ISMAIL_4
        File3 << GetCurrentEvent() << " " << shower.at(ishower).at(ipart)->pstat() << " "<< shower.at(ishower).at(ipart)->plabel() << " " << shower.at(ishower).at(ipart)->pid() << " " <<
                    shower.at(ishower).at(ipart)->color()<< " " <<
                    shower.at(ishower).at(ipart)->anti_color()<< " " <<
                    shower.at(ishower).at(ipart)->max_color()<< " " <<
                    shower.at(ishower).at(ipart)->px()<< " " <<
                    shower.at(ishower).at(ipart)->py()<< " " <<
                    shower.at(ishower).at(ipart)->pz()<< " " << onshellE << std::endl;
      #endif
    }
  }
  #ifdef DEBUG_ISMAIL_4
    File3 << "# Remnants now "<< std::endl;
  #endif

  auto ini  = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  auto Hard = JetScapeSignalManager::Instance()->GetHardProcessPointer().lock();

  auto MCGsecond = std::dynamic_pointer_cast<MCGlauberGenStringWrapper> (Hard->GetTaskList()[1]);
  auto Remnants = Hard->Remnants;
  if(2 * ini->pTHat.size() != Hard->Remnants.size()){
    throw std::runtime_error("Not enough remnants = " + std::to_string(Hard->Remnants.size()) + " Scattering = " + std::to_string(ini->pTHat.size()));
  }
  double NHardScatterings = double(ini->pTHat.size());

  for (unsigned int ipart = 0; ipart < Remnants.size(); ++ipart) {
      auto Rem = Remnants[ipart];
      double Pz, Px, Py, En;
      if(Rem.pz() >=0){
        En = MCGsecond->Get_Proj_Remnant()[0] / double(NHardScatterings);
        Px = MCGsecond->Get_Proj_Remnant()[1] / double(NHardScatterings);
        Py = MCGsecond->Get_Proj_Remnant()[2] / double(NHardScatterings);
        Pz = MCGsecond->Get_Proj_Remnant()[3] / double(NHardScatterings);
      } else {
        En = MCGsecond->Get_Targ_Remnant()[0] / double(NHardScatterings);
        Px = MCGsecond->Get_Targ_Remnant()[1] / double(NHardScatterings);
        Py = MCGsecond->Get_Targ_Remnant()[2] / double(NHardScatterings);
        Pz = MCGsecond->Get_Targ_Remnant()[3] / double(NHardScatterings);
      }

      std::cout << "Px = " << Px << " Py = " << Py << " Pz = " << Pz << " En = "<< En << " " << NHardScatterings <<  std::endl;

      double onshellE = pow(pow(Rem.px() + Px, 2) + pow(Rem.py() + Py, 2) + pow(Pz, 2),0.5);
      event.append(Rem.pid(), 23,
                   Rem.color(),
                   Rem.anti_color(),
                   Rem.px() + Px,
                   Rem.py() + Py,
                   Pz, onshellE);
      #ifdef DEBUG_ISMAIL_4
      File3 << GetCurrentEvent() << " "  << Rem.pstat() << " "<< Rem.plabel() << " " << Rem.pid() << " " <<
                   Rem.color()<< " " <<
                   Rem.anti_color()<< " " <<
                   Rem.max_color()<< " " <<
                   Rem.px()<< " " <<
                   Rem.py()<< " " <<
                   Pz << " " << onshellE << std::endl;
      #endif
  }

    //first, find unpaired color and anticolor tags.
    std::vector<int> cols;
    std::vector<int> acols;
    for (unsigned int ipart = 0; ipart < event.size(); ++ipart) {
      if (event[ipart].id() == 22) {
        continue;
      }
      if (event[ipart].col() != 0) {
        cols.push_back(event[ipart].col());
      }
      if (event[ipart].acol() != 0) {
        acols.push_back(event[ipart].acol());
      }
    }


    //the outcomes are: 1-unpaired color tag, 2-unpaired anticolor tag, 3-both an unpaired color & anticolor tag, 4-no unpaired tags
    //1-add an antiquark, 2-add a quark, 3-add a gluon, 4-add nothing (possibly photon only event)
    int icol = 0;
    while (icol < cols.size()) {
      bool foundpair = false;
      for (int iacol = 0; iacol < acols.size(); ++iacol) {
        if (cols[icol] == acols[iacol]) {
          cols.erase(cols.begin() + icol);
          acols.erase(acols.begin() + iacol);
          foundpair = true;
          continue;
        }
      }
      if (!foundpair) {
        ++icol;
      }
    }



    #ifdef DEBUG_ISMAIL_4
    File3 << "##########\n### unpaired Colors " << std::endl;
    for(int i =0; i<cols.size();i++ ){
      File3 << cols[i] << " ";
    }
    File3 << "\n### unpaired aColors " << std::endl;
    for(int i =0; i<acols.size();i++ ){
      File3 << acols[i] << " ";
    }
    File3 << "\n##########" << std::endl;
    #endif
    if(cols.size() > 0 || acols.size() > 0){
      std::cerr << "Unpaired colors sent to Pythia" << std::endl;
      exit(1);
    }

    int pid = 0;
    int color = 0;
    int anti_color = 0;
    if ((cols.size() > 0) && (acols.size() > 0)) {
      pid = 21;
      color = cols[0];
      anti_color = acols[0];
    } else if ((cols.size() > 0) && (acols.size() == 0)) {
      pid = -1;
      color = cols[0];
      anti_color = 0;
    } else if ((cols.size() == 0) && (acols.size() > 0)) {
      pid = 1;
      color = 0;
      anti_color = acols[0];
    }

    if (pid != 0) {
      pz = -1 * pz;
      event.append(pid, 23, anti_color, color, 0.2, 0.2, pz,
                   sqrt(pz * pz + 0.08));
    }

    VERBOSE(2) << "There are " << hOut.size() << " Hadrons and " << pOut.size()
               << " partons after Hadronization";
  
  
  pythia.next();
  // event.list();

  unsigned int ip = hOut.size();
  for (unsigned int i = 0; i < event.size(); ++i) {
    if (!event[i].isFinal())
      continue;
    //if ( !event[i].isHadron() )  continue;
    if (fabs(event[i].eta()) > 20)
      continue; //To prevent "nan" from propagating, very rare though

    double x[4] = {0, 0, 0, 0};
    hOut.push_back(make_shared<Hadron>(ip, event[i].id(), event[i].status(),
                                       event[i].pT(), event[i].eta(),
                                       event[i].phi(), event[i].e(), x));
    ++ip;
  }

  shower.clear();
  #ifdef DEBUG_ISMAIL_4
    File3.close();
  #endif
}
