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

#include "ColoredHadronizationSinglet.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include "tinyxml2.h"

using namespace Jetscape;
using namespace Pythia8;

// Register the module with the base class
RegisterJetScapeModule<ColoredHadronizationSinglet>
    ColoredHadronizationSinglet::reg("ColoredHadronizationSinglet");

Pythia8::Pythia ColoredHadronizationSinglet::pythia("IntentionallyEmpty", false);

ColoredHadronizationSinglet::ColoredHadronizationSinglet() {
  SetId("MyHadroTest");
  VERBOSE(8);
}

ColoredHadronizationSinglet::~ColoredHadronizationSinglet() { VERBOSE(8); }

void ColoredHadronizationSinglet::InitTask() {

  std::string s = GetXMLElementText({"JetHadronization", "name"});
  JSDEBUG << s << " to be initializied ...";

  double p_read_xml =
      GetXMLElementDouble({"JetHadronization", "eCMforHadronization"});
  p_fake = p_read_xml / 6.;

  /*std::string weak_decays =
      GetXMLElementText({"JetHadronization", "weak_decays"});*/

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

  // General settings for hadron decays
  std::string pythia_decays = GetXMLElementText({"JetHadronization", "pythia_decays"});
  double tau0Max = 10.0;
  double tau0Max_xml = GetXMLElementDouble({"JetHadronization", "tau0Max"});
	if(tau0Max_xml >= 0){tau0Max = tau0Max_xml;}
  else{JSWARN << "tau0Max should be larger than 0. Set it to 10.";}
  if(pythia_decays == "on"){
    JSINFO << "Pythia decays are turned on for tau0Max < " << tau0Max;
    pythia.readString("HadronLevel:Decay = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = " + std::to_string(tau0Max));
  } else {
    JSINFO << "Pythia decays are turned off";
    pythia.readString("HadronLevel:Decay = off");
  }

  // Settings for decays (old flag, will be depracted at some point)
  // This overwrites the previous settings if the user xml file contains the flag
  std::string weak_decays =
    GetXMLElementText({"JetHadronization", "weak_decays"});
  if (weak_decays == "off") {
    JSINFO << "Hadron decays are turned off.";
    JSWARN << "This parameter will be depracted at some point. Use 'pythia_decays' instead.\nOverwriting 'pythia_decays'.";
    pythia.readString("HadronLevel:Decay = off");
  } else if(weak_decays == "on") {
    JSINFO << "Hadron decays inside a range of 10 mm/c are turned on.";
    JSWARN << "This parameter will be depracted at some point. Use 'pythia_decays' and 'tau0Max' for more control on decays.\nOverwriting 'pythia_decays' and fix 'tau0Max' to 10.";
    pythia.readString("HadronLevel:Decay = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = 10.0");
  }

  std::stringstream lines;
  lines << GetXMLElementText({"JetHadronization", "LinesToRead"}, false);
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue; // skip empty lines
    JSINFO << "Also reading in: " << s;
    pythia.readString(s);
  }

  pythia.init();
}

void ColoredHadronizationSinglet::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  auto f = w.lock();
  if (!f)
    return;
  f->WriteComment("Hadronization Module : " + GetId());
  f->WriteComment("Hadronization to be implemented accordingly ...");
}

void ColoredHadronizationSinglet::DoHadronization(
    vector<vector<shared_ptr<Parton>>> &shower,
    vector<shared_ptr<Hadron>> &hOut, vector<shared_ptr<Parton>> &pOut) {

  // cout << "HADRONIZING\n\n\n" << endl;
  Event &event = pythia.event;
  event.reset();
  double pz = p_fake;

  //initial hadron list
  //pre existing hadrons
  if(hOut.size() > 0){
    for(vector<shared_ptr<Hadron>>::iterator hadIter = hOut.begin(); hadIter<hOut.end();){
      if(hadIter->get()->pid() > 23){ //skipping leptons
        double massnow = hadIter->get()->e()*hadIter->get()->e() -
                        (hadIter->get()->px()*hadIter->get()->px() + hadIter->get()->py()*hadIter->get()->py() + hadIter->get()->pz()*hadIter->get()->pz());
        massnow = (massnow >= 0.) ? sqrt(massnow) : -sqrt(-massnow);
        event.append(hadIter->get()->pid(),0,0,0,hadIter->get()->px(),hadIter->get()->py(),hadIter->get()->pz(),hadIter->get()->e(),massnow);
        event[event.size()-1].vProd(0., 0., 0., 0.);
      }
      hadIter++;
    }
  }

  int pythiastat;
  JSDEBUG << "&&&&&&&&&&&&&&&&&&& the number of showers are: " << shower.size();
  for (unsigned int ishower = 0; ishower < shower.size(); ++ishower) {
    JSDEBUG << "&&&&&&&&&&&&&&&&&&& there are " << shower.at(ishower).size()
            << " partons in the shower number " << ishower;
    for (unsigned int ipart = 0; ipart < shower.at(ishower).size(); ++ipart) {
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

      if (shower.at(ishower).at(ipart)->plabel() < 0) { pythiastat = 63; } //beam remnant
      else { pythiastat = 23; } //outgoing particle

      event.append(shower.at(ishower).at(ipart)->pid(), pythiastat,
                   shower.at(ishower).at(ipart)->color(),
                   shower.at(ishower).at(ipart)->anti_color(),
                   shower.at(ishower).at(ipart)->px(),
                   shower.at(ishower).at(ipart)->py(),
                   shower.at(ishower).at(ipart)->pz(), onshellE);
      event[event.size() - 1].vProd(shower.at(ishower).at(ipart)->x_in().comp(1)*Pythia8::FM2MM,
                                    shower.at(ishower).at(ipart)->x_in().comp(2)*Pythia8::FM2MM,
                                    shower.at(ishower).at(ipart)->x_in().comp(3)*Pythia8::FM2MM,
                                    shower.at(ishower).at(ipart)->x_in().comp(0)*Pythia8::FM2MM);
      //JSINFO << "parton " << shower.at(ishower).at(ipart)->pid() << " location from jetscape: " << event[event.size() - 1].vProd()*Pythia8::MM2FM;
    }
  }

  pythia.next();
  // event.list();

  unsigned int ip = hOut.size();
  for (unsigned int i = 0; i < event.size(); ++i) {
    if (!event[i].isFinal())
      continue;
    //if ( !event[i].isHadron() )  continue;
    if (fabs(event[i].eta()) > 20)
      continue; //To prevent "nan" from propagating, very rare though

    double x[4] = {event[i].vProd().e()*Pythia8::MM2FM, event[i].vProd().px()*Pythia8::MM2FM, event[i].vProd().py()*Pythia8::MM2FM, event[i].vProd().pz()*Pythia8::MM2FM};
    //JSINFO << "hadron " << event[i].id() << " location from pythia: " << event[i].vProd()*Pythia8::MM2FM;
    hOut.push_back(make_shared<Hadron>(ip, event[i].id(), event[i].status(),
                                       event[i].pT(), event[i].eta(),
                                       event[i].phi(), event[i].e(), x));
    ++ip;
  }

  shower.clear();
}
