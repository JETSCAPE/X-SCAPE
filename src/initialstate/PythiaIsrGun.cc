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

// Create a pythia collision at a specified point and return the two inital hard
// partons
#include "PythiaIsrGun.h"
#include <cmath>
#include <cstdint>
#include <random>
#include <sstream>
#include <iostream>
#include <fstream>
#define MAGENTA "\033[35m"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<PythiaIsrGun> PythiaIsrGun::reg("PythiaIsrGun");

PythiaIsrGun::~PythiaIsrGun() { VERBOSE(8); }

void PythiaIsrGun::InitTask() {
  JSDEBUG << "Initialize PythiaIsrGun";
  VERBOSE(8);

  std::string s = GetXMLElementText({"Hard", "PythiaGun", "name"});
  SetId(s);

  // SC: read flag for FSR
  FSR_on = GetXMLElementInt({"Hard", "PythiaGun", "FSR_on"});
  eCM = GetXMLElementDouble({"Hard", "PythiaGun", "eCM"});
  ReadPtHatBins();
  const int nBins = GetNPtHatBins();

  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
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

  // Window 0 keeps the seed; the others need their own streams (as in
  // PythiaGun).  Pythia's seeds run from 1 to 900,000,000; seed 0 (Pythia:
  // from the clock) gives the others fresh ones.
  seedBin_.assign(nBins, seed);
  std::random_device entropy;
  for (int k = 1; k < nBins; ++k) {
    if (seed == 0) {
      seedBin_[k] = 1 + entropy() % 900000000u;
    } else {
      uint64_t z = (static_cast<uint64_t>(seed) << 16) + k;  // splitmix64
      z += 0x9e3779b97f4a7c15ULL;
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      z ^= z >> 31;
      seedBin_[k] = 1 + static_cast<unsigned int>(z % 900000000ULL);
    }
  }

  // The extra instances copy this one's settings and particle data while they
  // are still the defaults (nothing has been read into them yet).
  extraPythia_.clear();
  for (int k = 1; k < nBins; ++k)
    extraPythia_.emplace_back(
        std::make_unique<Pythia8::Pythia>(settings, particleData, false));
  for (int k = 0; k < nBins; ++k)
    ConfigurePythia(PythiaOf(k), k);

  activeBin_ = 0;
  pTHatMin = pTHatBins_[0].first;
  pTHatMax = pTHatBins_[0].second;

  std::ofstream sigma_printer;
  sigma_printer.open(printer, std::ios::trunc);
}

//! <pTHatBins>: pairs of numbers, "min max min max ...".  Absent or empty: the
//! one window pTHatMin .. pTHatMax.  As PythiaGun::ReadPtHatBins, except that
//! the hydro is never reused (every event is complete).
void PythiaIsrGun::ReadPtHatBins() {
  pTHatMin = GetXMLElementDouble({"Hard", "PythiaGun", "pTHatMin"});
  pTHatMax = GetXMLElementDouble({"Hard", "PythiaGun", "pTHatMax"});
  pTHatBins_.clear();

  std::stringstream text(
      GetXMLElementText({"Hard", "PythiaGun", "pTHatBins"}, false));
  std::vector<double> edges;
  std::string token;
  while (text >> token) {
    try {
      size_t used = 0;
      edges.push_back(std::stod(token, &used));
      if (used != token.size())
        throw std::invalid_argument(token);
    } catch (const std::exception &) {
      throw std::runtime_error("PythiaIsrGun: <pTHatBins> takes numbers, got '" +
                               token + "'");
    }
  }
  if (edges.empty()) {
    pTHatBins_.emplace_back(pTHatMin, pTHatMax);
    return;
  }
  if (edges.size() % 2)
    throw std::runtime_error("PythiaIsrGun: <pTHatBins> needs pairs 'min "
                             "max', got an odd number of values");
  for (size_t j = 0; j < edges.size(); j += 2) {
    if (!(edges[j] > 0 && edges[j + 1] > edges[j]))
      throw std::runtime_error(
          "PythiaIsrGun: every <pTHatBins> window needs 0 < min < max");
    // ConfigurePythia hands the edges over with one decimal, as for pTHatMin.
    for (size_t e = j; e < j + 2; ++e)
      if (std::fabs(edges[e] * 10 - std::round(edges[e] * 10)) > 1e-6)
        throw std::runtime_error(
            "PythiaIsrGun: <pTHatBins> edges take at most one decimal "
            "(Pythia gets them rounded to one), got " +
            std::to_string(edges[e]));
    pTHatBins_.emplace_back(edges[j], edges[j + 1]);
  }
  const int nBins = GetNPtHatBins();
  if (nBins == 1)
    return;

  std::string lines =
      GetXMLElementText({"Hard", "PythiaGun", "LinesToRead"}, false);
  if (lines.find("PhaseSpace:pTHat") != std::string::npos)
    throw std::runtime_error(
        "PythiaIsrGun: <LinesToRead> sets PhaseSpace:pTHat..., which would "
        "override every <pTHatBins> window");

  std::ostringstream msg;
  for (const auto &b : pTHatBins_)
    msg << " [" << b.first << ", " << b.second << "]";
  JSINFO << MAGENTA << "Pythia ISR Gun with " << nBins
         << " pTHat windows (event i uses window i mod " << nBins
         << ", every event complete, no hydro reuse):" << msg.str();
}

Pythia8::Pythia &PythiaIsrGun::PythiaOf(int bin) {
  if (bin == 0)
    return *this;
  return *extraPythia_.at(bin - 1);
}

const Pythia8::Info &PythiaIsrGun::InfoOf(int bin) {
  if (bin == 0)
    return info;
  return extraPythia_.at(bin - 1)->info;
}

//! One window's Pythia instance: the settings of pTHatMin/pTHatMax with that
//! window's range and seed.  For window 0 this is the former InitTask body.
void PythiaIsrGun::ConfigurePythia(Pythia8::Pythia &py, int bin) {
  // Show initialization at INFO level
  py.readString("Init:showProcesses = off");
  py.readString("Init:showChangedSettings = off");
  py.readString("Init:showMultipartonInteractions = on");
  py.readString("Init:showChangedParticleData = off");
  if (JetScapeLogger::Instance()->GetInfo()) {
    py.readString("Init:showProcesses = on");
    py.readString("Init:showChangedSettings = on");
    py.readString("Init:showMultipartonInteractions = on");
    py.readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  py.readString("Next:numberShowInfo = 0");
  py.readString("Next:numberShowProcess = 0");
  py.readString("Next:numberShowEvent = 0");

  // Standard settings
  py.readString("HardQCD:all = on");  // will repeat this line in the xml for
                                      // demonstration
  py.readString("HardQCD:nQuarkNew = 3");  // Number Of Quark flavours
  py.readString("MultipartonInteractions:processLevel = 0");
  py.readString(
      "MultipartonInteractions:nQuarkIn = 3");  // Number Of Quark flavours

  py.readString("HadronLevel:Decay = off");
  py.readString("HadronLevel:all = on");
  py.readString("PartonLevel:ISR = off");
  py.readString("PartonLevel:MPI = on");
  py.readString("PromptPhoton:all=on");
  py.readString("WeakSingleBoson:all=off");
  py.readString("WeakDoubleBoson:all=off");

  // For parsing text
  stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);
  numbf.setf(ios::showpoint);
  numbf.precision(1);
  stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

  if (FSR_on)
    py.readString("PartonLevel:FSR = on");
  else
    py.readString("PartonLevel:FSR = off");

  const double binMin = pTHatBins_[bin].first;
  const double binMax = pTHatBins_[bin].second;
  JSINFO << MAGENTA << "Pythia Gun with FSR_on: " << FSR_on;
  JSINFO << MAGENTA << "Pythia Gun with " << binMin << " < pTHat < " << binMax
         << (GetNPtHatBins() > 1 ? " (window " + std::to_string(bin) + ")"
                                 : std::string());

  numbf.str("PhaseSpace:pTHatMin = ");
  numbf << binMin;
  py.readString(numbf.str());
  numbf.str("PhaseSpace:pTHatMax = ");
  numbf << binMax;
  py.readString(numbf.str());

  py.readString("Random:setSeed = on");
  numbi.str("Random:seed = ");
  VERBOSE(7) << "Seeding pythia to " << seedBin_[bin];
  numbi << seedBin_[bin];
  py.readString(numbi.str());

  // Species
  py.readString("Beams:idA = 2212");
  py.readString("Beams:idB = 2212");

  // Energy
  numbf.str("Beams:eCM = ");
  numbf << eCM;
  py.readString(numbf.str());

  std::stringstream lines;
  lines << GetXMLElementText({"Hard", "PythiaGun", "LinesToRead"}, false);
  std::string s;
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue;  // skip empty lines
    VERBOSE(7) << "Also reading in: " << s;
    py.readString(s);
  }

  // And initialize
  if (!py.init()) {  // Pythia>8.1
    throw std::runtime_error("Pythia init() failed.");
  }
}

void PythiaIsrGun::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  JetScapeTask::WriteTasks(w);
}

void PythiaIsrGun::ExecuteTask() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();

  // This event's pTHat window: event i uses window i mod K.
  const int nBins = GetNPtHatBins();
  activeBin_ = nBins > 1 ? GetCurrentEvent() % nBins : 0;
  pTHatMin = pTHatBins_[activeBin_].first;
  pTHatMax = pTHatBins_[activeBin_].second;
  Pythia8::Pythia &py = PythiaOf(activeBin_);
  if (nBins > 1)
    VERBOSE(1) << "pTHat window " << activeBin_ << ": " << pTHatMin
               << " < pTHat < " << pTHatMax;

  // Reading vir_factor from xml for MATTER
  double vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});

  bool flag62 = false;
  vector<Pythia8::Particle> p62;

  // sort by pt
  struct greater_than_pt {
    inline bool operator()(const Pythia8::Particle &p1,
                           const Pythia8::Particle &p2) {
      return (p1.pT() > p2.pT());
    }
  };

  FourVector p_p;

  int NSamplings = 0;
ReDoSampling:
  do {
    NSamplings++;
    py.next();
    p62.clear();

    std::vector<int> IndexToSkip;

    for (int parid = 0; parid < py.event.size(); parid++) {
      if (parid < 3)
        continue;  // 0, 1, 2: total event and beams
      Pythia8::Particle &particle = py.event[parid];

      if (!(particle.isGluon() ||
            particle.isQuark())) {  // Getting rid of diquark
        if (particle.status() == -31) {
          IndexToSkip.push_back(particle.daughter1());
          IndexToSkip.push_back(particle.daughter2());
          IndexToSkip.push_back(py.event[particle.daughter1()].mother1());
          IndexToSkip.push_back(py.event[particle.daughter1()].mother2());
        } else if (particle.status() == -33) {
          IndexToSkip.push_back(particle.mother1());
          IndexToSkip.push_back(particle.mother2());
          IndexToSkip.push_back(py.event[particle.mother1()].daughter1());
          IndexToSkip.push_back(py.event[particle.mother1()].daughter2());
        }
      }
    }

    for (auto &ToSkip : IndexToSkip) {
      JSWARN << " Skipping non-parton index " << ToSkip;
    }
    if (!printer.empty()) {
      std::ofstream sigma_printer;
      sigma_printer.open(printer, std::ios::out | std::ios::app);

      sigma_printer << "sigma = " << GetSigmaGen()
                    << " Err =  " << GetSigmaErr() << endl;
      // sigma_printer.close();

      //      JSINFO << BOLDYELLOW << " sigma = " << GetSigmaGen() << " sigma
      //      err = " << GetSigmaErr() << " printer = " << printer << " is " <<
      //      sigma_printer.is_open() ;
    };

    // pTarr[0]=0.0; pTarr[1]=0.0;
    // pindexarr[0]=0; pindexarr[1]=0;

    for (int parid = 0; parid < py.event.size(); parid++) {
      if (parid < 3)
        continue;  // 0, 1, 2: total event and beams

      Pythia8::Particle &particle = py.event[parid];
      if (!FSR_on) {
        if (!((particle.status() == -21) || (particle.status() == -23) ||
              (particle.status() == -31) || (particle.status() == -33)))
          continue;
        for (auto &ToSkip : IndexToSkip) {
          if (ToSkip == parid) {
            goto SkipParton;
          }
        }
        // if ( (particle.status() > -21)||(particle.status()<-23) ) continue ;

        // only accept particles after MPI
        // if (particle.status() != 62)
        // continue;
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        // if (fabs(particle.id()) > 5 &&
        //  (particle.id() != 21 && particle.id() != 22))
        // continue;

        // reject rare cases of very soft particles that don't have enough e to
        // get reasonable virtuality
        // if (particle.pT() < 1.0 / sqrt(vir_factor))
        // continue;

        // if(particle.id()==22) cout<<"########this is a photon!######" <<endl;
        //  accept
      } else {  // FSR_on true: use Pythia vacuum shower instead of MATTER
        if (!particle.isFinal())
          continue;
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        if (fabs(particle.id()) > 5 &&
            (particle.id() != 21 && particle.id() != 22))
          continue;
      }

      VERBOSE(1) << MAGENTA << " particle from pythiagun id = " << particle.id()
                 << " pz = " << particle.pz() << " px = " << particle.px()
                 << " py = " << particle.py() << " E =  " << particle.e()
                 << " status = " << particle.status() << " idex "
                 << particle.index() << " Color " << particle.col() << " "
                 << particle.acol() << " Mothers " << particle.mother1() << " "
                 << particle.mother2() << " daughter " << particle.daughter1()
                 << " " << particle.daughter2();

      p62.push_back(particle);

    SkipParton:;
    }

    // if you want at least 2
    if (p62.size() < 2)
      continue;
    // if ( p62.size() < 1 ) continue;

    // Now have all candidates, sort them
    // sort by pt
    // std::sort(p62.begin(), p62.end(), greater_than_pt());
    // // check...
    // for (auto& p : p62 ) cout << p.pT() << endl;

    flag62 = true;

  } while (!flag62);

  // // Roll for a starting point
  // // See:
  // https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  // std::random_device device;
  // std::mt19937 engine(device()); // Seed the random number engine

  FourVector x_p;

  if (!ini) {
    JSINFO << BOLDYELLOW
           << "No initial state module, setting the starting location to "
              "0. Make sure to add e.g. trento before PythiaIsrGun.";
  } else {
    double t, x, y, z;
    bool pass = false;
    //  while (!pass)
    //  {
    ini->SampleABinaryCollisionPoint(t, x, y, z);
    // JSINFO << MAGENTA << " pass at time " << t << " with x = " << x << " y =
    // " << y << " z = " << z << "  ? " ;
    // cin >> pass;
    //}
    x_p.Set(x, y, z, t);
    ini->OutputHardCollisionPosition(t, x, y, z);
  }

  // Loop through particles
  // Accept them all

  double initial_state_label = -1;
  double final_state_label = 1;
  ini->pTHat.resize((p62.size()) / 4);
  int hCounter = 0;
  SetTotalMomentumPositive(0.0);
  SetTotalMomentumNegative(0.0);
  SetTotalMomentumFractionPositive(0.0);
  SetTotalMomentumFractionNegative(0.0);
  double TotalEnergyOfInitialStatePartons = 0.0;

  for (int np = 0; np < p62.size(); ++np)
  // for (int np = p62.size()-1; np >= 0 ; --np)
  {
    Pythia8::Particle &particle = p62.at(np);

    if (particle.status() == -21 || particle.status() == -31) {
      TotalEnergyOfInitialStatePartons += particle.e();
      if (particle.pz() >= 0.0) {
        SetTotalMomentumPositive(GetTotalMomentumPositive() + particle.e());
        SetTotalMomentumFractionPositive(GetTotalMomentumFractionPositive() +
                                         (particle.e() + particle.pz()) /
                                             (0.94 * eCM));
      } else {
        SetTotalMomentumNegative(GetTotalMomentumNegative() + particle.e());
        SetTotalMomentumFractionNegative(GetTotalMomentumFractionNegative() +
                                         (particle.e() - particle.pz()) /
                                             (0.94 * eCM));
      }
    }
  }

  VERBOSE(2) << "Negative Partons Momentum " << GetTotalMomentumNegative()
             << " Positive Partons Momentum " << GetTotalMomentumPositive()
             << " eCM " << eCM << " TotalEnergyOfInitialStatePartons = "
             << TotalEnergyOfInitialStatePartons;
  if (GetTotalMomentumFractionNegative() >= 1.0 ||
      GetTotalMomentumFractionPositive() >= 1.) {
    JSINFO
        << "Redoing Pythia Sampling Since MPI energy is larger than eCM/2.1 ";
    if (NSamplings < 1000) {
      goto ReDoSampling;
    }
    JSWARN << "Negative Partons Momentum Fraction "
           << GetTotalMomentumFractionNegative()
           << " Positive Partons Momentum Fraction "
           << GetTotalMomentumFractionPositive() << " eCM " << eCM
           << " TotalEnergyOfInitialStatePartons = "
           << TotalEnergyOfInitialStatePartons;
    throw std::runtime_error(
        "Pythia Isr Gun outputs more energy in the MPI partons than eCM");
  }

  for (int np = 0; np < p62.size(); ++np)
  // for (int np = p62.size()-1; np >= 0 ; --np)
  {
    Pythia8::Particle &particle = p62.at(np);

    // VERBOSE(7) << "Adding particle with pid = " << particle.id()
    //            << ", pT = " << particle.pT() << ", y = " << particle.y()
    //            << ", phi = " << particle.phi() << ", e = " << particle.e();

    // JSINFO<< MAGENTA << " at x=" << x_p.x() << ", y=" << x_p.y() << ", z=" <<
    // x_p.z() << ", t = " << x_p.t();

    int label = 0;
    int stat = 0;

    if (particle.status() == -21 || particle.status() == -31) {
      label = initial_state_label;
      initial_state_label--;
      stat = -1000;  // raw initial state status, must go to an initial state
                     // module
    }
    if (particle.status() == -23 || particle.status() == -33) {
      label = final_state_label;
      final_state_label++;
      stat = 1000;  // raw final state status, must go to a final state module
                    // with virtuality generation.
      if ((label - 1) % 2 == 0) {
        ini->pTHat[(label - 1) / 2] = particle.pT();
      }
    }

    FourVector p_p(particle.px(), particle.py(), particle.pz(), particle.e());

    auto ptn = make_shared<Parton>(label, particle.id(), stat, p_p, x_p);
    ptn->set_color(particle.col());
    ptn->set_anti_color(particle.acol());
    ptn->set_max_color(GetMax_ColorPerShower() * (np + 1));
    AddParton(ptn);
  }

  // Getting Number of hard partons
  int NPP = p62.size();
  SetMax_Color(GetMax_ColorPerShower() * NPP);
  FourVector Zeros(0, 0, 0, 0);

  ini->CollisionNegativeMomentum = std::vector<FourVector>(NPP / 2, Zeros);
  ini->CollisionPositiveMomentum = std::vector<FourVector>(NPP / 2, Zeros);
  ini->CollisionNegativeRotatedMomentum =
      std::vector<FourVector>(NPP / 2, Zeros);
  ini->CollisionPositiveRotatedMomentum =
      std::vector<FourVector>(NPP / 2, Zeros);
  ini->ClearHardPartonMomentum();

  VERBOSE(8) << GetNHardPartons();

  // REMARK: Check why this has to be called explictly, something wrong with
  // generic recursive execution!!????

  // JetScapeTask::ExecuteTasks();
}
