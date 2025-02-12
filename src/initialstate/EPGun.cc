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

// Create a pythia collision at a specified point and return the two inital hard partons

#include "EPGun.h"
#include "Matter.h"
#include <sstream>
#include <iostream>
#include <fstream>
#define MAGENTA "\033[35m"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<EPGun> EPGun::reg("EPGun");

EPGun::~EPGun() { VERBOSE(8); }

void EPGun::InitTask() {

  JSDEBUG << "Initialize EPGun";
  VERBOSE(8);

  // Show initialization at INFO level
  readString("Init:showProcesses = off");
  readString("Init:showChangedSettings = off");
  readString("Init:showMultipartonInteractions = off");
  readString("Init:showChangedParticleData = off");
  if (JetScapeLogger::Instance()->GetInfo()) {
    readString("Init:showProcesses = on");
    readString("Init:showChangedSettings = on");
    readString("Init:showMultipartonInteractions = on");
    readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  readString("Next:numberShowInfo = 0");
  readString("Next:numberShowProcess = 0");
  readString("Next:numberShowEvent = 0");

  // For parsing text
  stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);
  numbf.setf(ios::showpoint);
  numbf.precision(1);
  stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

  std::string s = GetXMLElementText({"Hard", "EPGun", "name"});
  SetId(s);
  // cout << s << endl;

  // initial kinematics
  eElectron = GetXMLElementDouble({"Hard", "EPGun", "electron_energy"});
  eProton = GetXMLElementDouble({"Hard", "EPGun", "proton_energy"});
  use_positron = GetXMLElementInt({"Hard", "EPGun", "use_positron"});
  photoproduction = GetXMLElementInt({"Hard", "EPGun", "photoproduction"});
  breitVir = GetXMLElementInt({"Hard", "EPGun", "breit_vir"});
  Q2pow = GetXMLElementDouble({"Hard", "EPGun", "Q2_pow"});
  Q2factor = GetXMLElementDouble({"Hard", "EPGun", "Q2_factor"});
  initial_virtuality_pT = GetXMLElementInt({"Eloss", "Matter", "initial_virtuality_pT"});

  //kinematic cuts
  Q2min = GetXMLElementDouble({"Hard", "EPGun", "Q2min"});
  Q2max = GetXMLElementDouble({"Hard", "EPGun", "Q2max"});
  W2min = GetXMLElementDouble({"Hard", "EPGun", "W2min"});
  W2max = GetXMLElementDouble({"Hard", "EPGun", "W2max"});
  xmin = GetXMLElementDouble({"Hard", "EPGun", "xmin"});
  xmax = GetXMLElementDouble({"Hard", "EPGun", "xmax"});
  ymin = GetXMLElementDouble({"Hard", "EPGun", "ymin"});
  ymax = GetXMLElementDouble({"Hard", "EPGun", "ymax"});

  //other Pythia settings
  readString("HadronLevel:Decay = off");
  readString("HadronLevel:all = off");
  
  // EP Gun stuff
  readString("Beams:frameType = 2");
  // BeamA = proton.
  readString("Beams:idA = 2212");
  settings.parm("Beams:eA", eProton);
  // BeamB = electron.
  if(use_positron){
    readString("Beams:idB = -11");
    JSINFO << "Running with positron beam.";
  }else
    readString("Beams:idB = 11");
  settings.parm("Beams:eB", eElectron);

  if(photoproduction){
    readString("PDF:lepton2gamma = on");
    readString("PhotonParton:all = on");
    readString("Photon:Q2max = 1.0");
    readString("Photon:ProcessType = 0");
    readString("SoftQCD:nonDiffractive = on");
  }
  else{
    // Set up DIS process within some phase space.
    // Neutral current (with gamma/Z interference).
    readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
    // Uncomment to allow charged current.
    //readString("WeakBosonExchange:ff2ff(t:W) = on");
    // Phase-space cut: minimal Q2 of process.
    settings.parm("PhaseSpace:Q2Min", Q2min);

    // Set dipole recoil on. Necessary for DIS + shower.
    readString("SpaceShower:dipoleRecoil = on");

    // Allow emissions up to the kinematical limit,
    // since rate known to match well to matrix elements everywhere.
    readString("SpaceShower:pTmaxMatch = 2");

    // QED radiation off lepton not handled yet by the new procedure.
    readString("TimeShower:QEDshowerByL = off");
    readString("PartonShowers:model = 1");
    readString("TimeShower:pTmaxMatch = 1");

    //special PDF
    readString("PDF:lepton = off");
    readString("PDF:useHard = on");
    //readString("PDF:pHardSet = LHAPDF6:PDF4LHC21_40"); //for special PDF setting
  }

  // SC: read flag for FSR
  FSR_on = GetXMLElementInt({"Hard", "EPGun", "FSR_on"});
  if (FSR_on)
    readString("PartonLevel:FSR = on");
  else
    readString("PartonLevel:FSR = off");

  JSINFO << MAGENTA << "EP Gun with FSR_on: " << FSR_on;

  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
  readString("Random:setSeed = on");
  numbi.str("Random:seed = ");
  unsigned int seed = 0;
  if (RandomXmlDescription) {
    tinyxml2::XMLElement *xmle =
        RandomXmlDescription->FirstChildElement("seed");
    if (!xmle)
      throw std::runtime_error("Cannot parse xml");
    xmle->QueryUnsignedText(&seed);
  } else {
    JSWARN << "No <Random> element found in xml, seeding to 0";
  }
  VERBOSE(7) << "Seeding pythia to " << seed;
  numbi << seed;
  readString(numbi.str());

  //Reading vir_factor from xml for MATTER
  vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});
  softMomentumCutoff = GetXMLElementDouble({"Hard", "EPGun", "softMomentumCutoff"});
  initial_virtuality_pT = GetXMLElementInt({"Eloss", "Matter", "initial_virtuality_pT"});
  if(vir_factor < rounding_error) {
    JSWARN << "vir_factor should not be zero or negative";
    exit(1);
  }

  std::stringstream lines;
  lines << GetXMLElementText({"Hard", "EPGun", "LinesToRead"}, false);
  int i = 0;
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue; // skip empty lines
    VERBOSE(7) << "Also reading in: " << s;
    readString(s);
  }

  // And initialize
  if (!init()) { // Pythia>8.1
    throw std::runtime_error("Pythia init() failed.");
  }

  std::ofstream sigma_printer;
  sigma_printer.open(printer, std::ios::trunc);

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double>{0.0, 1.0};
}

void EPGun::ExecuteTask() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();

  bool flag62 = false;
  vector<Pythia8::Particle> p62;

  // sort by pt
  struct greater_than_pt {
    inline bool operator()(const Pythia8::Particle &p1,
                           const Pythia8::Particle &p2) {
      return (p1.pT() > p2.pT());
    }
  };

  do {
    bool check = next();
    if(check ==  false) continue;

    //getting scattered electron index
    int elecID = 6;
    if(photoproduction){
      for(int iElec=0; iElec<event.size(); iElec++){
        if(abs(event[iElec].id()) == 11 and event[iElec].status() == 23)
          elecID = iElec;
      }
    }

    //kinematic cuts
    Pythia8::Vec4 pProton = event[1].p();
    Pythia8::Vec4 peIn    = event[2].p();
    Pythia8::Vec4 peOut   = event[6].p();
    Pythia8::Vec4 pPhoton = peIn - peOut;

    // Q2, W2, Bjorken x, y.
    double Q2    = - pPhoton.m2Calc();
    double W2    = (pProton + pPhoton).m2Calc();
    double x     = Q2 / (2. * pProton * pPhoton);
    double y     = (pProton * pPhoton) / (pProton * peIn);

    if(x < xmin or x > xmax) continue;
    if(y < ymin or y > ymax) continue;
    if(Q2 < Q2min or Q2 > Q2max) continue;
    if(W2 < W2min or W2 > W2max) continue;

    //JSINFO << "Q2 = " << Q2 << "; W2 = " << W2 << "; x = " << x << "; y = " << y;

    p62.clear();
      if (!printer.empty()){
            std::ofstream sigma_printer;
            sigma_printer.open(printer, std::ios::out | std::ios::app);

            sigma_printer << "sigma = " << GetSigmaGen() << " Err =  " << GetSigmaErr() << endl ;
            //sigma_printer.close();


//      JSINFO << BOLDYELLOW << " sigma = " << GetSigmaGen() << " sigma err = " << GetSigmaErr() << " printer = " << printer << " is " << sigma_printer.is_open() ;
    };

    // pTarr[0]=0.0; pTarr[1]=0.0;
    // pindexarr[0]=0; pindexarr[1]=0;

    for (int parid = 0; parid < event.size(); parid++) {
      if (parid < 3)
        continue; // 0, 1, 2: total event and beams
      Pythia8::Particle &particle = event[parid];

      //skipping everything decayed
      if (!particle.isFinal())
        continue;

      //replacing diquarks with antiquarks (and anti-dq's with quarks)
      //the id is set to the heaviest quark in the diquark (except down quark)
      //this technically violates baryon number conservation over the entire event
      //also can violate electric charge conservation
      if( (std::abs(particle.id()) > 1100) && (std::abs(particle.id()) < 6000) && ((std::abs(particle.id())/10)%10 == 0) ){
        if(particle.id() > 0){particle.id( -1*particle.id()/1000 );}
        else{particle.id( particle.id()/1000 );}
      }

      //catching scattered electron and beam remenants
      if(particle.isHadron() or particle.isLepton()){
        AddHadron(EPGun::PythiaToJSHadron(particle));
        continue;
      }

      if (!FSR_on) {
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        if (fabs(particle.id()) > 5 &&
            (particle.id() != 21 && particle.id() != 22))
          continue;

        // reject rare cases of very soft particles that don't have enough e to get
        // reasonable virtuality
        if (initial_virtuality_pT && (particle.pT() < softMomentumCutoff)) {
          // this cutoff was 1.0/sqrt(vir_factor) in versions < 3.6
          continue;
        } else if(!initial_virtuality_pT && (particle.pAbs() < softMomentumCutoff)) {
          continue;
        }

        //if(particle.id()==22) cout<<"########this is a photon!######" <<endl;
        // accept
      } else { // FSR_on true: use Pythia vacuum shower instead of MATTER
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        if (fabs(particle.id()) > 5 &&
            (particle.id() != 21 && particle.id() != 22))
          continue;
      }
      p62.push_back(particle);
    }

    // if you want at least 2
    //if (p62.size() < 2) continue;
    if ( p62.size() < 1 ) continue;

    // Now have all candidates, sort them
    // sort by pt
    std::sort(p62.begin(), p62.end(), greater_than_pt());
    // // check...
    // for (auto& p : p62 ) cout << p.pT() << endl;

    flag62 = true;

  } while (!flag62);

  double p[4], xLoc[4];

  // This location should come from an initial state
  for (int i = 0; i <= 3; i++) {
    xLoc[i] = 0.0;
  };

  // // Roll for a starting point
  // // See: https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  // std::random_device device;
  // std::mt19937 engine(device()); // Seed the random number engine

  if (!ini) {
    VERBOSE(1) << "No initial state module, setting the starting location to "
                  "0. Make sure to add e.g. trento before EPGun.";
  } else {
    double t,x, y,z;
    ini->SampleABinaryCollisionPoint(t,x, y,z);
    xLoc[1] = x;
    xLoc[2] = y;
  }

  // Loop through particles

  // Only top two
  //for(int np = 0; np<2; ++np){

  // Accept them all

  //getting Breit frame for the event to set virtualities in
  Pythia8::Vec4 pProton = event[1].p();
  Pythia8::Vec4 peIn    = event[4].p();
  Pythia8::Vec4 peOut   = event[6].p();
  Pythia8::Vec4 pPhoton = peIn - peOut;
  double Q2    = - pPhoton.m2Calc();
  double W2    = (pProton + pPhoton).m2Calc();
  double x     = Q2 / (2. * pProton * pPhoton);
  double y     = (pProton * pPhoton) / (pProton * peIn);
  //Pythia8::Vec4 pBreit  = 2*x*pProton + pPhoton;
  Pythia8::Vec4 pQuark  = 2*x*pProton;
  Pythia8::Vec4 pBreit2  = 2*x*pProton + pPhoton;
  Pythia8::RotBstMatrix breitBoost = Pythia8::toCMframe(pQuark,pPhoton);
  Pythia8::Vec4 pQuarkI= event[3].p();
  Pythia8::Vec4 pQuarkF= event[5].p();
  pQuarkI.rotbst(breitBoost);
  pQuarkF.rotbst(breitBoost);

  //test statements
  pBreit2.rotbst(breitBoost);
  pQuark.rotbst(breitBoost);
  pPhoton.rotbst(breitBoost);
  JSINFO << "pQuark Initial: " << pQuarkI.px() << " " << pQuarkI.py() << " " << pQuarkI.pz() << " ";
  JSINFO << "pQuark Final: " << pQuarkF.px() << " " << pQuarkF.py() << " " << pQuarkF.pz() << " ";
  JSINFO << "pPhoton: " << pPhoton.px() << " " << pPhoton.py() << " " << pPhoton.pz() << " ";
  JSINFO << "Breit frame test: " << pBreit2.px() << " " << pBreit2.py() << " " << pBreit2.pz() << " ";
  const double QS = 0.9;

  int hCounter = 0;
  for (int np = 0; np < p62.size(); ++np) {
    Pythia8::Particle &particle = p62.at(np);
    Pythia8::Vec4 partp = particle.p();
    partp.rotbst(breitBoost);
    double mass = particle.m();
    double eCM = info.eCM();

    //only doing vir setting for DIS
    if(breitVir and particle.status() == 62){
      //setting up max vir
      double max_vir = (partp.pAbs() * partp.pAbs() - mass*mass) * vir_factor;
      double min_vir = (QS * QS / 2.0) * (1.0 + std::sqrt(1.0 + 4.0 * particle.m() * particle.m() / QS / QS));
      double tQ2 = 0.;

      //using z axis for pT since thats the axis the photon quark collision happens on
      if(initial_virtuality_pT){
        max_vir = (partp.pz() * partp.pz() - mass*mass) * vir_factor;
      }

      //JSINFO << Q2factor;
      max_vir *= pow(Q2/(info.s()),Q2pow) * Q2factor/sqrt(x);
      //JSINFO << max_vir;

      int iSplit = 0; // quark
      if (particle.id() == gid) {
        JSDEBUG << " parton is a gluon ";
        iSplit = 1; // gluon
      } else {
        JSDEBUG << " parton is a quark ";
      }

      //evaluating virtuality for different cases
      if (max_vir <= QS * QS){
        tQ2 = 0.0;
      }else{
        double nu = (partp.e() + partp.pAbs())/sqrt(2.0);

        if (abs(particle.id()) == 4 || abs(particle.id()) == 5) {
          if (max_vir > min_vir) {
              tQ2 =
                  matterHelper.generate_vac_t_w_M(particle.id(), particle.m(), nu,
                                    QS * QS / 2.0, max_vir, 0, iSplit);
            } else {
              tQ2 = QS * QS;
            }
            //  std::ofstream tdist;
            //  tdist.open("tdist_heavy.dat", std::ios::app);
            //  tdist << tQ2 << endl;
            //  tdist.close();

            VERBOSE(8) << BOLDYELLOW << " virtuality calculated as = " << tQ2;
        } else if (particle.id() == gid) {
            tQ2 = matterHelper.generate_vac_t_w_M(particle.id(), particle.m(), nu,
                                    QS * QS / 2.0, max_vir, 0, iSplit);
        } else {
            tQ2 = matterHelper.generate_vac_t(particle.id(), nu, QS * QS / 2.0,
                                max_vir, 0, iSplit);

        }

        //tQ2 = test_vir;
        //catching virtualitiies that are too high
        if(sqrt(tQ2) > particle.pAbs() /*or sqrt(tQ2) > partp.pAbs()*/) tQ2 = min_vir;
      }

      //applying virtuality change to the parton
      //JSINFO << BOLDYELLOW << "Particle with ID: " << particle.id();
      //JSINFO << BOLDYELLOW << "initial momentum: " << particle.px() << " " << particle.py() << " " << particle.pz();
      //JSINFO << BOLDYELLOW << "breit momentum: " << partp.px() << " " << partp.py() << " " << partp.pz();
      //JSINFO << BOLDYELLOW << "Virtuality: " << sqrt(tQ2) << " ";

      double scale = sqrt(particle.e()*particle.e() - tQ2 - particle.m2())/particle.pAbs();
      particle.px(particle.px()*scale);
      particle.py(particle.py()*scale);
      particle.pz(particle.pz()*scale);
    }

    VERBOSE(7) << "Adding particle with pid = " << particle.id()
               << " at x=" << xLoc[1] << ", y=" << xLoc[2] << ", z=" << xLoc[3];

    VERBOSE(7) << "Adding particle with pid = " << particle.id()
               << ", pT = " << particle.pT() << ", y = " << particle.y()
               << ", phi = " << particle.phi() << ", e = " << particle.e();

    VERBOSE(7) << " at x=" << xLoc[1] << ", y=" << xLoc[2] << ", z=" << xLoc[3];

    auto ptn = make_shared<Parton>(0, particle.id(), 0, particle.pT(), particle.eta(),particle.phi(), particle.e(), xLoc);
    ptn->set_color(particle.col());
    ptn->set_anti_color(particle.acol());
    ptn->set_max_color(1000 * (np + 1));

    //adding mean formtime for partons that need it set
    if(breitVir and particle.status() == 62){
      double mean_form_time = (2.*ptn->e()) / (ptn->e()*ptn->e()
                            - ptn->px()*ptn->px() - ptn->py()*ptn->py()
                            - ptn->pz()*ptn->pz() - ptn->restmass()*ptn->restmass()
                            + rounding_error) / fmToGeVinv;
      ptn->set_form_time(mean_form_time);
      ptn->set_mean_form_time();

      double velocity[4];
      velocity[0] = 1.0;
      for (int j = 1; j <= 3; j++) {
        velocity[j] = ptn->p(j) / ptn->e();
      }
      ptn->set_jet_v(velocity);
    }

    AddParton(ptn);
  }

  VERBOSE(8) << GetNHardPartons();

}
