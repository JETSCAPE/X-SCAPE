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

#include "PythiaTGun.h"
#include <sstream>
#include <cstdlib>
//#include <iostream>
//#include <fstream>
#define MAGENTA "\033[35m"
 
using namespace std;

// Register the module with the base class
RegisterJetScapeModule<PythiaTGun> PythiaTGun::reg("CustomModulePythiaTGun");

PythiaTGun::PythiaTGun() : HardProcess() {
  SetId("PythiaTGun");
  VERBOSE(8);
}


PythiaTGun::~PythiaTGun() { VERBOSE(8); }

void PythiaTGun::InitTask() {

  JSINFO <<MAGENTA<< "Initialize PythiaTGun";

  std::string s = GetXMLElementText({"Hard", "CustomModulePythiaTGun", "name"});
  SetId(s);

  pTHatMin = GetXMLElementDouble({"Hard", "CustomModulePythiaTGun", "pTHatMin"});
  pTHatMax = GetXMLElementDouble({"Hard", "CustomModulePythiaTGun", "pTHatMax"});
  ImpParMean = GetXMLElementDouble({"Hard", "CustomModulePythiaTGun", "ImpParMean"});
  s_1x = GetXMLElementDouble({"Hard", "CustomModulePythiaTGun", "s_1x"});
  s_1y = GetXMLElementDouble({"Hard", "CustomModulePythiaTGun", "s_1y"});
  n_1x = GetXMLElementInt({"Hard", "CustomModulePythiaTGun", "n_1x"});
  n_1y = GetXMLElementInt({"Hard", "CustomModulePythiaTGun", "n_1y"});
  BeamId = GetXMLElementInt({"Hard", "CustomModulePythiaTGun", "BeamId"});

  FSR_on = GetXMLElementInt({"Hard", "CustomModulePythiaTGun", "FSR_on"});

  eCM = GetXMLElementDouble({"Hard", "CustomModulePythiaTGun", "eCM"});
  vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});
  initial_virtuality_pT =
      GetXMLElementInt({"Eloss", "Matter", "initial_virtuality_pT"});
  softMomentumCutoff =
      GetXMLElementDouble({"Hard", "CustomModulePythiaTGun", "softMomentumCutoff"});

  unsigned int baseSeed = 0;
  tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
  if (RandomXmlDescription) {
    tinyxml2::XMLElement *xmle =
        RandomXmlDescription->FirstChildElement("seed");
    if (!xmle)
      throw std::runtime_error("Cannot parse xml");
    xmle->QueryUnsignedText(&baseSeed);
  }

  if (vir_factor < rounding_error) {
    JSWARN << "vir_factor should not be zero or negative";
    exit(1);
  }

  JSINFO << MAGENTA << "Pythia TGun with FSR_on: " << FSR_on;
  JSINFO << MAGENTA << "Pythia TGun with " << pTHatMin
         << " < pTHat < " << pTHatMax;

  double dx = 2.0 * s_1x / double(n_1x);
  //double dx = s_1x / double(n_1x-1);
  double dy = s_1y / double(n_1y);
  //double dy = (s_1y) / double(n_1y-1);

  for (int ix = 0; ix < n_1x; ix++) {
    for (int iy = 0; iy < n_1y; iy++) {
      //JSINFO << MAGENTA << "Initializing Pythia for grid point (" << ix << ", "
      //       << iy << ") with x = " << -s_1x + dx / 2.0 + ix * dx
      //       << " and y = " << 0.0 + dy / 2.0 + iy * dy;
      JSINFO << MAGENTA << "Initializing Pythia for grid point (" << ix << ", "
             << iy << ") with x = " << -s_1x + dx / 2.0 + ix * dx
             << " and y = " << 0.0 + dy / 2.0 + iy * dy;
      double x = -s_1x + dx / 2.0 + ix * dx;
      //double x = ix * dx;
      double y = 0.0 + dy / 2.0 + iy * dy;
      //double y = iy * dy;

      double s1 = std::sqrt(std::pow(x - ImpParMean / 2.0, 2)
                            + std::pow(y, 2));
      double s2 = std::sqrt(std::pow(x + ImpParMean / 2.0, 2)
                            + std::pow(y, 2));
      //system("ps -o rss= -p $$");
      auto py = std::make_unique<Pythia8::Pythia>("", false);
      //system("ps -o rss= -p $$");

      ConfigurePythia(*py, baseSeed + ix * n_1y + iy);

      auto protonA = py->getPDFPtr(2212, 1, "A", true);
      auto protonB = py->getPDFPtr(2212, 1, "B", true);

      auto pdfA = std::make_shared<EPS09s>(BeamId, 1, 1, protonA);
      auto pdfB = std::make_shared<EPS09s>(BeamId, 1, 1, protonB);

      pdfA->setTpos(s1);
      pdfB->setTpos(s2);

      pdfA->setSideLabel("A");
      pdfB->setSideLabel("B");

      py->setPDFPtr(pdfA, pdfB);

      if (!py->init()) {
        throw std::runtime_error("Pythia init() failed.");
      }

      pythia_vec.push_back(std::move(py));
    }
  }

  std::ofstream sigma_printer;
  sigma_printer.open(printer, std::ios::trunc);
}

void PythiaTGun::Test(int ipy){
  const double xmin  = 1.0e-6;
  const double xmax  = 1.0;
  const double q2min = 1.69;
  const double q2max = 1000000.0;
  int nQ2 = 80;
  int nX = 50;
  auto logGridValue = [](double minVal, double maxVal, int i, int n) {
    if (n <= 1) return minVal;
    const double t = double(i) / double(n - 1);
    return minVal * std::exp(t * std::log(maxVal / minVal));
  };
  const std::vector<int> partonIds = {
      21,   // gluon
      2, -2,
      1, -1,
      3, -3,
      4, -4,
      5, -5
  };
  std::string outFile = "xf_table_" + std::to_string(ipy) + ".txt";
  std::ofstream fout(outFile.c_str(), std::ios::out);
  fout << std::scientific << std::setprecision(10);
  double xfA, xfB;
  for (int id : partonIds) {
    for (int iq = 0; iq < nQ2; ++iq) {
      const double Q2 = logGridValue(q2min, q2max, iq, nQ2);

      for (int ix = 0; ix < nX; ++ix) {
        const double x = logGridValue(xmin, xmax, ix, nX);
        Pythia8::Pythia& py = *pythia_vec[ipy];
        auto pdfptrA = py.getInUsePDFPtr("A");
        auto pdfptrB = py.getInUsePDFPtr("B");
        //pdfptrA->xfUpdate(id, x, Q2);
        //pdfptrB->xfUpdate(id, x, Q2);
        xfA = pdfptrA->xf(id, x, Q2);
        xfB = pdfptrB->xf(id, x, Q2);
         fout << id << "  "
             << x << "  "
             << Q2 << "  "
             << xfA << "  "
             << xfB << "\n";
      }

      fout << "\n";
    }

    fout << "\n\n";
  }
  fout.close();

  JSINFO << MAGENTA
         << "Wrote average xf nuclear modification table to: "
         << outFile;
}

void PythiaTGun::ExecuteTask() {
  Test(1);
  exit(1);
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();
  double p[4], xLoc[4];
  for (int i = 0; i <= 3; i++) {
    xLoc[i] = 0.0;
  };

  if (!ini) {
    VERBOSE(1) << "No initial state module, setting the starting location to "
                  "0. Make sure to add e.g. trento before PythiaGun.";
  } else {
    double t,x, y,z;
    ini->SampleABinaryCollisionPoint(t,x, y,z);
    xLoc[1] = x;
    xLoc[2] = y;
  }
  // Now map sampled position to trained PYTHIA grid
  int ipy = GetPythiaGridIndex(xLoc[1], xLoc[2]);

  Pythia8::Pythia& py = *pythia_vec[ipy];

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
    py.next();
    p62.clear();
      if (!printer.empty()){
            std::ofstream sigma_printer;
            sigma_printer.open(printer, std::ios::out | std::ios::app);

            sigma_printer << "sigma = " << py.info.sigmaGen() << " Err =  " << py.info.sigmaErr() << endl ;
            //sigma_printer.close();

//      JSINFO << BOLDYELLOW << " sigma = " << GetSigmaGen() << " sigma err = " << GetSigmaErr() << " printer = " << printer << " is " << sigma_printer.is_open() ;
    };

    // pTarr[0]=0.0; pTarr[1]=0.0;
    // pindexarr[0]=0; pindexarr[1]=0;

    for (int parid = 0; parid < py.event.size(); parid++) {
      if (parid < 3)
        continue; // 0, 1, 2: total event and beams
      Pythia8::Particle &particle = py.event[parid];

      //replacing diquarks with antiquarks (and anti-dq's with quarks)
      //the id is set to the heaviest quark in the diquark (except down quark)
      //this technically violates baryon number conservation over the entire event
      //also can violate electric charge conservation
      if( (std::abs(particle.id()) > 1100) && (std::abs(particle.id()) < 6000) && ((std::abs(particle.id())/10)%10 == 0) ){
        if(particle.id() > 0){particle.id( -1*particle.id()/1000 );}
        else{particle.id( particle.id()/1000 );}
      }

      if (!FSR_on) {
        // only accept particles after MPI
        if (particle.status() != 62)
          continue;
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
        if (!particle.isFinal())
          continue;
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        if (fabs(particle.id()) > 5 &&
            (particle.id() != 21 && particle.id() != 22))
          continue;
      }
      p62.push_back(particle);
    }

    // if you want at least 2
    if (p62.size() < 2)
      continue;
    //if ( p62.size() < 1 ) continue;

    // Now have all candidates, sort them
    // sort by pt
    std::sort(p62.begin(), p62.end(), greater_than_pt());
    // // check...
    // for (auto& p : p62 ) cout << p.pT() << endl;

    //skipping event if softQCD is on & pThat exceeds max (where next bin is HardQCD with this as pThatmin)
    if(softQCD && (py.info.pTHat() >= pTHatMax)){continue;}

    flag62 = true;

  } while (!flag62);
  CrossSection += py.info.sigmaGen();

  int hCounter = 0;
  for (int np = 0; np < p62.size(); ++np) {
    Pythia8::Particle &particle = p62.at(np);

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
    AddParton(ptn);
  }

  VERBOSE(8) << GetNHardPartons();

}

void PythiaTGun::ConfigurePythia(Pythia8::Pythia& py, unsigned int seed) {

  py.readString("Init:showProcesses = off");
  py.readString("Init:showChangedSettings = off");
  py.readString("Init:showMultipartonInteractions = off");
  py.readString("Init:showChangedParticleData = off");

  if (JetScapeLogger::Instance()->GetInfo()) {
    py.readString("Init:showProcesses = on");
    py.readString("Init:showChangedSettings = on");
    py.readString("Init:showMultipartonInteractions = on");
    py.readString("Init:showChangedParticleData = on");
  }

  py.readString("Print:quiet = on");
  py.readString("Next:numberShowInfo = 0");
  py.readString("Next:numberShowProcess = 0");
  py.readString("Next:numberShowEvent = 0");

  py.readString("HadronLevel:Decay = off");
  py.readString("HadronLevel:all = off");
  py.readString("PartonLevel:ISR = on");
  py.readString("PartonLevel:MPI = on");
  py.readString("PromptPhoton:all=on");
  py.readString("WeakSingleBoson:all=off");
  py.readString("WeakDoubleBoson:all=off");

  std::stringstream numbf;
  numbf.setf(std::ios::fixed, std::ios::floatfield);
  numbf.setf(std::ios::showpoint);
  numbf.precision(1);

  if (pTHatMin < 0.01) {
    py.readString("HardQCD:all = off");
    py.readString("SoftQCD:nonDiffractive = on");
  } else {
    py.readString("HardQCD:all = on");

    numbf.str("");
    numbf.clear();
    numbf << "PhaseSpace:pTHatMin = " << pTHatMin;
    py.readString(numbf.str());

    numbf.str("");
    numbf.clear();
    numbf << "PhaseSpace:pTHatMax = " << pTHatMax;
    py.readString(numbf.str());
  }

  if (FSR_on)
    py.readString("PartonLevel:FSR = on");
  else
    py.readString("PartonLevel:FSR = off");

  py.readString("Random:setSeed = on");

  std::stringstream numbi;
  numbi << "Random:seed = " << seed;
  py.readString(numbi.str());

  py.readString("Beams:idA = 2212");
  py.readString("Beams:idB = 2212");

  numbf.str("");
  numbf.clear();
  numbf << "Beams:eCM = " << eCM;
  py.readString(numbf.str());

  std::stringstream lines;
  std::string s;
  lines << GetXMLElementText({"Hard", "CustomModulePythiaTGun", "LinesToRead"}, false);

  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue;
    VERBOSE(7) << "Also reading in: " << s;
    py.readString(s);
  }
}

int PythiaTGun::GetPythiaGridIndex(double x, double y) const {

  y = std::abs(y);  // top-bottom symmetry

  int ix = static_cast<int>(std::floor((x - (-s_1x)) / (2.0 * s_1x / double(n_1x))));
  //int ix = static_cast<int>(std::floor(x));
  int iy = static_cast<int>(std::floor((y - 0.0) / (s_1y / double(n_1y))));
  //int iy = static_cast<int>(std::floor(y));

  if (ix < 0) ix = 0;
  if (ix >= n_1x) ix = n_1x - 1;

  if (iy < 0) iy = 0;
  if (iy >= n_1y) iy = n_1y - 1;
  return ix * n_1y + iy;
}