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

#include "PythiaGun.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#define MAGENTA "\033[35m"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<PythiaGun> PythiaGun::reg("PythiaGun");

PythiaGun::~PythiaGun() { VERBOSE(8); }

void PythiaGun::InitTask() {
  JSDEBUG << "Initialize PythiaGun";
  VERBOSE(8);

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

  ReadPtHatBins();
  ReadPartonYCut();
  ReadTrigger();
  ReadExtraHard(seed);
  const int nBins = GetNPtHatBins();

  // Window 0 keeps the seed as it is (a single window is unchanged).  The other
  // windows need their own streams: Pythia seeds from Random:seed only, so equal
  // seeds would repeat window 0's random numbers.  Pythia's seeds run from 1 to
  // 900,000,000; seed 0 (Pythia: from the clock) gives the others fresh ones.
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
  extraHardPythia_.reset();
  if (extraOn_)
    extraHardPythia_ =
        std::make_unique<Pythia8::Pythia>(settings, particleData, false);

  softQCDBin_.assign(nBins, false);
  yTried_.assign(CutsOn() ? nBins : 0, 0);
  yKept_.assign(CutsOn() ? nBins : 0, 0);
  for (int k = 0; k < nBins; ++k)
    ConfigurePythia(PythiaOf(k), k);

  activeBin_ = 0;
  pTHatMin = pTHatBins_[0].first;
  pTHatMax = pTHatBins_[0].second;
  softQCD = softQCDBin_[0];
  if (extraOn_)
    ConfigureExtraPythia();  // after window 0: FSR_on, eCM, the cutoffs

  std::ofstream sigma_printer;
  sigma_printer.open(printer, std::ios::trunc);
}

//! <pTHatBins>: pairs of numbers, "min max min max ...".  Absent or empty:
//! the one window pTHatMin .. pTHatMax.
void PythiaGun::ReadPtHatBins() {
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
      throw std::runtime_error("PythiaGun: <pTHatBins> takes numbers, got '" +
                               token + "'");
    }
  }
  if (edges.empty()) {
    pTHatBins_.emplace_back(pTHatMin, pTHatMax);
    return;
  }
  if (edges.size() % 2)
    throw std::runtime_error("PythiaGun: <pTHatBins> needs pairs 'min max', "
                             "got an odd number of values");
  for (size_t j = 0; j < edges.size(); j += 2) {
    if (!(edges[j] >= 0 && edges[j + 1] > edges[j]))
      throw std::runtime_error(
          "PythiaGun: every <pTHatBins> window needs 0 <= min < max");
    // ConfigurePythia hands the edges over with one decimal, as for pTHatMin.
    for (size_t e = j; e < j + 2; ++e)
      if (std::fabs(edges[e] * 10 - std::round(edges[e] * 10)) > 1e-6)
        throw std::runtime_error(
            "PythiaGun: <pTHatBins> edges take at most one decimal (Pythia "
            "gets them rounded to one), got " + std::to_string(edges[e]));
    pTHatBins_.emplace_back(edges[j], edges[j + 1]);
  }
  const int nBins = GetNPtHatBins();
  if (nBins == 1)
    return;

  // Settings that would apply to every window at once.
  std::string lines =
      GetXMLElementText({"Hard", "PythiaGun", "LinesToRead"}, false);
  if (lines.find("PhaseSpace:pTHat") != std::string::npos)
    throw std::runtime_error(
        "PythiaGun: <LinesToRead> sets PhaseSpace:pTHat..., which would "
        "override every <pTHatBins> window");

  // Hydro reuse: every hydro event should get the same number of jets per
  // window.  Event i uses window i mod K and a reuse group starts at a multiple
  // of nReuseHydro, so that needs nReuseHydro % K == 0.
  std::string reuse = GetXMLElementText({"setReuseHydro"}, false);
  int nReuse = GetXMLElementInt({"nReuseHydro"}, false);
  if (reuse.find("true") != std::string::npos && nReuse % nBins != 0)
    throw std::runtime_error(
        "PythiaGun: nReuseHydro = " + std::to_string(nReuse) +
        " is not a multiple of the " + std::to_string(nBins) +
        " <pTHatBins> windows, so the hydro events would not get the same "
        "number of jets per window");

  std::ostringstream msg;
  for (const auto &b : pTHatBins_)
    msg << " [" << b.first << ", " << b.second << "]";
  JSINFO << MAGENTA << "Pythia Gun with " << nBins
         << " pTHat windows (event i uses window i mod " << nBins
         << "):" << msg.str();
}

//! <partonYMax> (empty: no cut) and <partonYMode> (leading | both | any).
void PythiaGun::ReadPartonYCut() {
  partonYMax_ = 0.;
  std::stringstream text(
      GetXMLElementText({"Hard", "PythiaGun", "partonYMax"}, false));
  std::string token;
  if (text >> token) {
    size_t used = 0;
    try {
      partonYMax_ = std::stod(token, &used);
    } catch (const std::exception &) {
      used = 0;
    }
    if (used != token.size() || !(partonYMax_ > 0))
      throw std::runtime_error(
          "PythiaGun: <partonYMax> takes a rapidity > 0 (empty: no cut), got '" +
          token + "'");
  }
  std::stringstream mode(
      GetXMLElementText({"Hard", "PythiaGun", "partonYMode"}, false));
  partonYMode_ = "leading";
  mode >> partonYMode_;
  if (partonYMode_ != "leading" && partonYMode_ != "both" &&
      partonYMode_ != "any")
    throw std::runtime_error(
        "PythiaGun: <partonYMode> is leading, both or any, got '" +
        partonYMode_ + "'");
  if (partonYMax_ > 0)
    JSINFO << MAGENTA << "Pythia Gun: handed-over partons, two hardest, |y| < "
           << partonYMax_ << " (" << partonYMode_
           << "); sigma = sigmaGen x kept/tried";
}

namespace {

//! An optional cut value: empty -> 0 (no cut), else a number > 0.
double ReadOptionalCut(const std::string &text, const std::string &tag) {
  std::stringstream in(text);
  std::string token;
  if (!(in >> token))
    return 0.;
  size_t used = 0;
  double v = 0.;
  try {
    v = std::stod(token, &used);
  } catch (const std::exception &) {
    used = 0;
  }
  if (used != token.size() || !(v > 0))
    throw std::runtime_error("PythiaGun: <" + tag +
                             "> takes a number > 0 (empty: no cut), got '" +
                             token + "'");
  return v;
}

}  // namespace

//! <minPartons> (default 2), <trigger> (none | Z) and the Z acceptance
//! <triggerYMax>, <leptonPtMin>, <leptonEtaMax> (empty: no cut).
void PythiaGun::ReadTrigger() {
  minPartons_ = 2;
  std::stringstream mp(
      GetXMLElementText({"Hard", "PythiaGun", "minPartons"}, false));
  std::string token;
  if (mp >> token) {
    size_t used = 0;
    try {
      minPartons_ = std::stoi(token, &used);
    } catch (const std::exception &) {
      used = 0;
    }
    if (used != token.size() || minPartons_ < 1)
      throw std::runtime_error(
          "PythiaGun: <minPartons> takes an integer >= 1, got '" + token + "'");
  }
  std::stringstream trig(
      GetXMLElementText({"Hard", "PythiaGun", "trigger"}, false));
  trigger_ = "none";
  trig >> trigger_;
  if (trigger_ != "none" && trigger_ != "Z")
    throw std::runtime_error("PythiaGun: <trigger> is none or Z, got '" +
                             trigger_ + "'");
  triggerYMax_ = ReadOptionalCut(
      GetXMLElementText({"Hard", "PythiaGun", "triggerYMax"}, false),
      "triggerYMax");
  leptonPtMin_ = ReadOptionalCut(
      GetXMLElementText({"Hard", "PythiaGun", "leptonPtMin"}, false),
      "leptonPtMin");
  leptonEtaMax_ = ReadOptionalCut(
      GetXMLElementText({"Hard", "PythiaGun", "leptonEtaMax"}, false),
      "leptonEtaMax");
  if (trigger_ != "Z" &&
      (triggerYMax_ > 0 || leptonPtMin_ > 0 || leptonEtaMax_ > 0))
    throw std::runtime_error("PythiaGun: <triggerYMax>, <leptonPtMin> and "
                             "<leptonEtaMax> need <trigger>Z</trigger>");
  if (minPartons_ != 2)
    JSINFO << MAGENTA << "Pythia Gun: events need " << minPartons_
           << " parton(s) (leptons do not count)";
  if (trigger_ == "Z") {
    std::ostringstream cuts;
    if (triggerYMax_ > 0)
      cuts << ", |y_Z| < " << triggerYMax_;
    if (leptonPtMin_ > 0)
      cuts << ", lepton pT > " << leptonPtMin_;
    if (leptonEtaMax_ > 0)
      cuts << ", lepton |eta| < " << leptonEtaMax_;
    JSINFO << MAGENTA << "Pythia Gun: trigger Z (an e+e- or mu+mu- pair)"
           << cuts.str() << "; sigma = sigmaGen x kept/tried";
  }
}

namespace {

//! A seed from the job seed and a salt (splitmix64), as the windows' seeds.
uint64_t MixSeed(unsigned int seed, uint64_t salt) {
  uint64_t z = (static_cast<uint64_t>(seed) << 16) + salt;
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

//! An optional number: empty -> ``fallback``.
double ReadNumber(const std::string &text, const std::string &tag,
                  double fallback) {
  std::stringstream in(text);
  std::string token;
  if (!(in >> token))
    return fallback;
  size_t used = 0;
  double v = 0.;
  try {
    v = std::stod(token, &used);
  } catch (const std::exception &) {
    used = 0;
  }
  if (used != token.size())
    throw std::runtime_error("PythiaGun: <extraHard><" + tag +
                             "> takes a number, got '" + token + "'");
  return v;
}

//! 3dMCGlauber's sigma_NN^inel(sqrt s) [mb] (Glauber::
//! compute_NN_inelastic_cross_section), which its N_coll is made with.
double GlauberSigmaNN(double ecm) {
  const double ls = std::log(ecm * ecm);
  return (44.4 - 2.9 * ls + 0.33 * ls * ls) -
         (11.4 - 1.52 * ls + 0.13 * ls * ls);
}

}  // namespace

//! <extraHard>: <on>, <noTrigger>, <pTHatMin>, <sigmaNN>, <warmup>, <fixedN>,
//! <LinesToRead> (empty values: the defaults).
void PythiaGun::ReadExtraHard(unsigned int seed) {
  auto text = [this](const char *tag) {
    return GetXMLElementText({"Hard", "PythiaGun", "extraHard", tag}, false);
  };
  extraOn_ = ReadNumber(text("on"), "on", 0.) != 0.;
  extraNoTrigger_ = ReadNumber(text("noTrigger"), "noTrigger", 0.) != 0.;
  extraPtHatMin_ = ReadNumber(text("pTHatMin"), "pTHatMin", 15.);
  extraWarmup_ = static_cast<long>(ReadNumber(text("warmup"), "warmup", 20000.));
  extraFixedN_ = static_cast<int>(ReadNumber(text("fixedN"), "fixedN", -1.));
  const double ecm = GetXMLElementDouble({"Hard", "PythiaGun", "eCM"});
  extraSigmaNN_ = ReadNumber(text("sigmaNN"), "sigmaNN", 0.);
  if (!(extraSigmaNN_ > 0))
    extraSigmaNN_ = GlauberSigmaNN(ecm);
  nExtra_ = 0;
  extraMu_ = 0.;
  if (!extraOn_) {
    if (extraNoTrigger_)
      throw std::runtime_error(
          "PythiaGun: <extraHard><noTrigger> needs <extraHard><on>1");
    return;
  }
  if (!(extraPtHatMin_ > 0) || extraWarmup_ < 1)
    throw std::runtime_error("PythiaGun: <extraHard> needs <pTHatMin> > 0 and "
                             "<warmup> >= 1");
  if (extraNoTrigger_ &&
      (GetNPtHatBins() > 1 || partonYMax_ > 0 || trigger_ != "none"))
    throw std::runtime_error(
        "PythiaGun: <extraHard><noTrigger> makes no triggering scattering, so "
        "<pTHatBins>, <partonYMax> and <trigger> do not apply");
  // their own streams: the trigger's Pythia and the framework's engines are
  // untouched
  if (seed == 0) {
    std::random_device entropy;
    extraPythiaSeed_ = 1 + entropy() % 900000000u;
    extraRngSeed_ = entropy();
  } else {
    extraPythiaSeed_ =
        1 + static_cast<unsigned int>(MixSeed(seed, 0xE1A0) % 900000000ULL);
    extraRngSeed_ = static_cast<unsigned int>(MixSeed(seed, 0xE1A1));
  }
  extraRng_.seed(extraRngSeed_);
  JSINFO << MAGENTA << "Pythia Gun: extra hard scatterings above pTHat "
         << extraPtHatMin_ << " GeV, N ~ "
         << (extraFixedN_ >= 0 ? "fixed " + std::to_string(extraFixedN_)
                               : std::string("Poisson(N_coll sigma/sigma_NN)"))
         << ", sigma_NN " << extraSigmaNN_ << " mb"
         << (extraNoTrigger_ ? "; no triggering scattering" : "");
}

//! The extras' Pythia: HardQCD above <extraHard><pTHatMin>, ISR on, MPI off,
//! FSR as <FSR_on>; then <warmup> events for their cross section.
void PythiaGun::ConfigureExtraPythia() {
  Pythia8::Pythia &py = *extraHardPythia_;
  for (const char *s :
       {"Init:showProcesses = off", "Init:showChangedSettings = off",
        "Init:showMultipartonInteractions = off",
        "Init:showChangedParticleData = off", "Next:numberShowInfo = 0",
        "Next:numberShowProcess = 0", "Next:numberShowEvent = 0",
        "HadronLevel:Decay = off", "HadronLevel:all = off",
        "PartonLevel:ISR = on", "PartonLevel:MPI = off", "HardQCD:all = on",
        "Beams:idA = 2212", "Beams:idB = 2212", "Random:setSeed = on"})
    py.readString(s);
  py.readString(FSR_on ? "PartonLevel:FSR = on" : "PartonLevel:FSR = off");
  std::ostringstream v;
  v.precision(17);
  v << "PhaseSpace:pTHatMin = " << extraPtHatMin_;
  py.readString(v.str());
  py.readString("Random:seed = " + std::to_string(extraPythiaSeed_));
  v.str("");
  v << "Beams:eCM = " << eCM;
  py.readString(v.str());
  std::stringstream lines;
  lines << GetXMLElementText({"Hard", "PythiaGun", "extraHard", "LinesToRead"},
                             false);
  std::string s;
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue;
    py.readString(s);
  }
  if (!py.init())
    throw std::runtime_error("PythiaGun: the <extraHard> Pythia init() failed");
  for (long i = 0; i < extraWarmup_; ++i)
    py.next();
  extraSigma_ = py.info.sigmaGen();
  extraSigmaErr_ = py.info.sigmaErr();
  JSINFO << MAGENTA << "Pythia Gun: extra scatterings' sigma(pTHat > "
         << extraPtHatMin_ << " GeV) = " << extraSigma_ << " +- "
         << extraSigmaErr_ << " mb from " << extraWarmup_
         << " events (Pythia seed " << extraPythiaSeed_ << ")";
}

//! The two hardest partons to hand over (``sorted`` by pT; leptons skipped, at
//! least one parton) against |y| < partonYMax_: the hardest ("leading"), both
//! ("both") or either ("any").  With one parton every mode cuts on it.
bool PythiaGun::PassesPartonYCut(
    const std::vector<Pythia8::Particle> &sorted) const {
  const Pythia8::Particle *two[2] = {nullptr, nullptr};
  int n = 0;
  for (const auto &p : sorted) {
    if (Parton::isLepton(p.id()))
      continue;
    two[n++] = &p;
    if (n == 2)
      break;
  }
  if (n == 0)
    return false;
  const bool in0 = std::fabs(two[0]->y()) < partonYMax_;
  const bool in1 = n > 1 ? std::fabs(two[1]->y()) < partonYMax_ : in0;
  if (partonYMode_ == "leading")
    return in0;
  if (partonYMode_ == "both")
    return in0 && in1;
  return in0 || in1;
}

//! The Z: the opposite-charge, same-flavour lepton pair with the largest
//! invariant mass; triggerLeptons_ hardest first (``sorted`` is by pT).
bool PythiaGun::FindTriggerPair(const std::vector<Pythia8::Particle> &sorted) {
  triggerLeptons_.clear();
  double best = -1.;
  for (size_t i = 0; i < sorted.size(); ++i) {
    if (!Parton::isLepton(sorted[i].id()))
      continue;
    for (size_t j = i + 1; j < sorted.size(); ++j) {
      if (sorted[j].id() != -sorted[i].id())
        continue;
      const double m = (sorted[i].p() + sorted[j].p()).mCalc();
      if (m > best) {
        best = m;
        triggerLeptons_ = {sorted[i], sorted[j]};
      }
    }
  }
  return !triggerLeptons_.empty();
}

bool PythiaGun::PassesTriggerCut() const {
  const Pythia8::Vec4 z = triggerLeptons_[0].p() + triggerLeptons_[1].p();
  if (triggerYMax_ > 0 && !(std::fabs(z.rap()) < triggerYMax_))
    return false;
  for (const auto &l : triggerLeptons_) {
    if (leptonPtMin_ > 0 && !(l.pT() > leptonPtMin_))
      return false;
    if (leptonEtaMax_ > 0 && !(std::fabs(l.eta()) < leptonEtaMax_))
      return false;
  }
  return true;
}

double PythiaGun::GetSigmaErr(int bin) {
  const double e = InfoOf(bin).sigmaErr();
  if (!CutsOn())
    return e;
  const double s = InfoOf(bin).sigmaGen(), a = GetYAcceptance(bin);
  const long n = GetNYTried(bin);
  const double ea = n > 0 ? std::sqrt(a * (1. - a) / n) : 0.;
  return std::sqrt(e * a * e * a + s * ea * s * ea);
}

Pythia8::Pythia &PythiaGun::PythiaOf(int bin) {
  if (bin == 0)
    return *this;
  return *extraPythia_.at(bin - 1);
}

const Pythia8::Info &PythiaGun::InfoOf(int bin) {
  if (bin == 0)
    return info;
  return extraPythia_.at(bin - 1)->info;
}

//! One window's Pythia instance: the settings of pTHatMin/pTHatMax, with that
//! window's range and seed.  For window 0 this is the former InitTask body.
void PythiaGun::ConfigurePythia(Pythia8::Pythia &py, int bin) {
  const double binMin = pTHatBins_[bin].first;
  const double binMax = pTHatBins_[bin].second;

  // Show initialization at INFO level
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

  // No event record printout.
  py.readString("Next:numberShowInfo = 0");
  py.readString("Next:numberShowProcess = 0");
  py.readString("Next:numberShowEvent = 0");

  // For parsing text
  stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);
  numbf.setf(ios::showpoint);
  numbf.precision(1);
  stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

  std::string s = GetXMLElementText({"Hard", "PythiaGun", "name"});
  SetId(s);
  // cout << s << endl;

  // other Pythia settings
  py.readString("HadronLevel:Decay = off");
  py.readString("HadronLevel:all = off");
  py.readString("PartonLevel:ISR = on");
  py.readString("PartonLevel:MPI = on");
  // py.readString("PartonLevel:FSR = on");
  py.readString("PromptPhoton:all=on");
  py.readString("WeakSingleBoson:all=off");
  py.readString("WeakDoubleBoson:all=off");

  if (binMin < 0.01) {  // assuming low bin where softQCD should be used
    // running softQCD - inelastic nondiffrative (min-bias)
    py.readString("HardQCD:all = off");
    py.readString("SoftQCD:nonDiffractive = on");
    softQCDBin_[bin] = true;
  } else {                              // running normal hardQCD
    py.readString("HardQCD:all = on");  // will repeat this line in the xml for
                                        // demonstration
    //  py.readString("HardQCD:gg2ccbar = on"); // switch on heavy quark channel
    // py.readString("HardQCD:qqbar2ccbar = on");
    numbf.str("PhaseSpace:pTHatMin = ");
    numbf << binMin;
    py.readString(numbf.str());
    numbf.str("PhaseSpace:pTHatMax = ");
    numbf << binMax;
    py.readString(numbf.str());
    softQCDBin_[bin] = false;
  }

  // SC: read flag for FSR
  FSR_on = GetXMLElementInt({"Hard", "PythiaGun", "FSR_on"});
  if (FSR_on)
    py.readString("PartonLevel:FSR = on");
  else
    py.readString("PartonLevel:FSR = off");

  JSINFO << MAGENTA << "Pythia Gun with FSR_on: " << FSR_on;
  JSINFO << MAGENTA << "Pythia Gun with " << binMin << " < pTHat < " << binMax
         << (GetNPtHatBins() > 1 ? " (window " + std::to_string(bin) + ")"
                                 : std::string());

  // random seed (read in InitTask)
  py.readString("Random:setSeed = on");
  numbi.str("Random:seed = ");
  VERBOSE(7) << "Seeding pythia to " << seedBin_[bin];
  numbi << seedBin_[bin];
  py.readString(numbi.str());

  // Species
  py.readString("Beams:idA = 2212");
  py.readString("Beams:idB = 2212");

  // Energy
  eCM = GetXMLElementDouble({"Hard", "PythiaGun", "eCM"});
  numbf.str("Beams:eCM = ");
  numbf << eCM;
  py.readString(numbf.str());

  // Reading vir_factor from xml for MATTER
  vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});
  initial_virtuality_pT =
      GetXMLElementInt({"Eloss", "Matter", "initial_virtuality_pT"});
  if (vir_factor < rounding_error) {
    JSWARN << "vir_factor should not be zero or negative";
    exit(1);
  }

  softMomentumCutoff =
      GetXMLElementDouble({"Hard", "PythiaGun", "softMomentumCutoff"});

  std::stringstream lines;
  lines << GetXMLElementText({"Hard", "PythiaGun", "LinesToRead"}, false);
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

void PythiaGun::ExecuteTask() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();

  partonScatter_.clear();
  scatters_.clear();
  nExtra_ = 0;
  extraMu_ = 0.;
  if (extraNoTrigger_) {  // only the extras (ReadExtraHard: no windows, cuts)
    AddExtraScatterings();
    VERBOSE(8) << GetNHardPartons();
    return;
  }

  // This event's pTHat window: event i uses window i mod K.
  const int nBins = GetNPtHatBins();
  activeBin_ = nBins > 1 ? GetCurrentEvent() % nBins : 0;
  pTHatMin = pTHatBins_[activeBin_].first;
  pTHatMax = pTHatBins_[activeBin_].second;
  softQCD = softQCDBin_[activeBin_];
  Pythia8::Pythia &py = PythiaOf(activeBin_);
  if (nBins > 1)
    VERBOSE(1) << "pTHat window " << activeBin_ << ": " << pTHatMin
               << " < pTHat < " << pTHatMax;

  bool flag62 = false;
  vector<Pythia8::Particle> p62;
  triggerLeptons_.clear();
  constexpr long kMaxRejected = 100000;  //!< tries per event at the cuts
  long nRejected = 0;

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
    if (!printer.empty()) {
      std::ofstream sigma_printer;
      sigma_printer.open(printer, std::ios::out | std::ios::app);

      sigma_printer << "sigma = " << GetSigmaGen()
                    << " Err =  " << GetSigmaErr();
      if (nBins > 1)
        sigma_printer << " window " << activeBin_;
      sigma_printer << endl;
      // sigma_printer.close();

      //      JSINFO << BOLDYELLOW << " sigma = " << GetSigmaGen() << " sigma
      //      err = " << GetSigmaErr() << " printer = " << printer << " is " <<
      //      sigma_printer.is_open() ;
    };

    // pTarr[0]=0.0; pTarr[1]=0.0;
    // pindexarr[0]=0; pindexarr[1]=0;

    SelectPartons(py, p62);

    // at least <minPartons> (default 2); leptons do not count
    const long nPartons =
        std::count_if(p62.begin(), p62.end(), [](const Pythia8::Particle &p) {
          return !Parton::isLepton(p.id());
        });
    if (nPartons < minPartons_)
      continue;

    // Now have all candidates, sort them
    // sort by pt
    std::sort(p62.begin(), p62.end(), greater_than_pt());
    // // check...
    // for (auto& p : p62 ) cout << p.pT() << endl;

    // skipping event if softQCD is on & pThat exceeds max (where next bin is
    // HardQCD with this as pThatmin)
    if (softQCD && (py.info.pTHat() >= pTHatMax)) {
      continue;
    }

    // the cuts on what is handed over: rapidity (<partonYMax>), trigger
    if (CutsOn()) {
      ++yTried_[activeBin_];
      // a process that never passes (e.g. <trigger>Z without a Z in
      // <LinesToRead>) would regenerate forever
      if (++nRejected > kMaxRejected)
        throw std::runtime_error(
            "PythiaGun: no event passed <partonYMax> / <trigger> in " +
            std::to_string(kMaxRejected) +
            " tries; with <trigger>Z, does <LinesToRead> make Z -> e+e- or "
            "mu+mu-?");
      if (partonYMax_ > 0 && !PassesPartonYCut(p62))
        continue;
      if (trigger_ == "Z" && !(FindTriggerPair(p62) && PassesTriggerCut()))
        continue;
      ++yKept_[activeBin_];
    }

    flag62 = true;

  } while (!flag62);

  double xLoc[4] = {0.0, 0.0, 0.0, 0.0};

  if (!ini) {
    VERBOSE(1) << "No initial state module, setting the starting location to "
                  "0. Make sure to add e.g. trento before PythiaGun.";
  } else {
    double t, x, y, z;
    ini->SampleABinaryCollisionPoint(t, x, y, z);
    xLoc[1] = x;
    xLoc[2] = y;
  }

  // Accept them all
  AddScatterPartons(p62, xLoc[1], xLoc[2], 0, py.info.pTHat(), py.info.code());

  if (extraOn_)
    AddExtraScatterings();

  VERBOSE(8) << GetNHardPartons();
}

//! The particles of ``py``'s event to hand over, in event order.
void PythiaGun::SelectPartons(Pythia8::Pythia &py,
                              std::vector<Pythia8::Particle> &p62) {
  for (int parid = 0; parid < py.event.size(); parid++) {
    if (parid < 3)
      continue;  // 0, 1, 2: total event and beams
    Pythia8::Particle &particle = py.event[parid];

    // replacing diquarks with antiquarks (and anti-dq's with quarks)
    // the id is set to the heaviest quark in the diquark (except down quark)
    // this technically violates baryon number conservation over the entire
    // event also can violate electric charge conservation
    if ((std::abs(particle.id()) > 1100) &&
        (std::abs(particle.id()) < 6000) &&
        ((std::abs(particle.id()) / 10) % 10 == 0)) {
      if (particle.id() > 0) {
        particle.id(-1 * particle.id() / 1000);
      } else {
        particle.id(particle.id() / 1000);
      }
    }

    if (!FSR_on) {
      // leptons (a Z's decay: final, status 23 after the resonance decay,
      // not 62) go over as they are, like photons; not partons, so no
      // softMomentumCutoff
      if (Parton::isLepton(particle.id())) {
        if (particle.isFinal())
          p62.push_back(particle);
        continue;
      }
      // only accept particles after MPI
      if (particle.status() != 62)
        continue;
      // only accept gluons and quarks
      // Also accept Gammas to put into the hadron's list
      if (fabs(particle.id()) > 5 &&
          (particle.id() != 21 && particle.id() != 22))
        continue;

      // reject rare cases of very soft particles that don't have enough e to
      // get reasonable virtuality
      if (initial_virtuality_pT && (particle.pT() < softMomentumCutoff)) {
        // this cutoff was 1.0/sqrt(vir_factor) in versions < 3.6
        continue;
      } else if (!initial_virtuality_pT &&
                 (particle.pAbs() < softMomentumCutoff)) {
        continue;
      }

      // if(particle.id()==22) cout<<"########this is a photon!######" <<endl;
      //  accept
    } else {  // FSR_on true: use Pythia vacuum shower instead of MATTER
      if (!particle.isFinal())
        continue;
      // only accept gluons and quarks
      // Also accept Gammas and leptons to put into the hadron's list
      if (fabs(particle.id()) > 5 &&
          (particle.id() != 21 && particle.id() != 22) &&
          !Parton::isLepton(particle.id()))
        continue;
    }
    p62.push_back(particle);
  }
}

//! Hands ``sorted`` (by pT) over at (x, y), z = t = 0, as scattering
//! ``scatter``: label ``scatter``, colour tags offset by 100000 x scatter (none
//! for the trigger), and records it.
void PythiaGun::AddScatterPartons(const std::vector<Pythia8::Particle> &sorted,
                                  double x, double y, int scatter,
                                  double pthat, int code) {
  double xLoc[4] = {0.0, x, y, 0.0};
  const int offset = 100000 * scatter;
  ScatterRecord rec{scatter, pthat, code, x, y,
                    static_cast<int>(sorted.size()), 0., 0.};
  for (int np = 0; np < static_cast<int>(sorted.size()); ++np) {
    const Pythia8::Particle &particle = sorted.at(np);

    VERBOSE(7) << "Adding particle with pid = " << particle.id()
               << ", pT = " << particle.pT() << ", y = " << particle.y()
               << ", phi = " << particle.phi() << ", e = " << particle.e()
               << " at x=" << xLoc[1] << ", y=" << xLoc[2] << ", z=" << xLoc[3]
               << " (scattering " << scatter << ")";

    auto ptn = make_shared<Parton>(scatter, particle.id(), 0, particle.pT(),
                                   particle.eta(), particle.phi(), particle.e(),
                                   xLoc);
    ptn->set_color(particle.col() > 0 ? particle.col() + offset : particle.col());
    ptn->set_anti_color(particle.acol() > 0 ? particle.acol() + offset
                                            : particle.acol());
    ptn->set_max_color(offset + 1000 * (np + 1));
    AddParton(ptn);
    partonScatter_.push_back(scatter);
    rec.e += particle.e();
    rec.pt += particle.pT();
  }
  scatters_.push_back(rec);
}

//! The extras: N ~ Poisson(N_coll sigma/sigma_NN) (or <fixedN>), each at its
//! own binary collision, from their own Pythia and random stream.
void PythiaGun::AddExtraScatterings() {
  if (!ini)
    throw std::runtime_error(
        "PythiaGun: <extraHard> needs an initial state (3dMCGlauber) for "
        "N_coll and the binary collisions");
  const double ncoll = ini->GetNcoll();
  if (!(ncoll >= 0))
    throw std::runtime_error(
        "PythiaGun: <extraHard> needs an initial state that knows N_coll "
        "(3dMCGlauber), got " + std::to_string(ncoll));
  extraMu_ = ncoll * extraSigma_ / extraSigmaNN_;
  if (extraFixedN_ >= 0)
    nExtra_ = extraFixedN_;
  else if (extraMu_ > 0)
    nExtra_ = std::poisson_distribution<int>(extraMu_)(extraRng_);
  else
    nExtra_ = 0;

  struct greater_than_pt {
    inline bool operator()(const Pythia8::Particle &p1,
                           const Pythia8::Particle &p2) {
      return (p1.pT() > p2.pT());
    }
  };
  Pythia8::Pythia &py = *extraHardPythia_;
  std::vector<Pythia8::Particle> sel;
  for (int s = 1; s <= nExtra_; ++s) {
    int tries = 0;
    while (!py.next())
      if (++tries > 1000)
        throw std::runtime_error("PythiaGun: the <extraHard> Pythia failed "
                                 "1000 times in a row");
    sel.clear();
    SelectPartons(py, sel);
    std::sort(sel.begin(), sel.end(), greater_than_pt());
    double t, x, y, z;
    ini->SampleABinaryCollisionPointWith(extraRng_, t, x, y, z);
    AddScatterPartons(sel, x, y, s, py.info.pTHat(), py.info.code());
  }
  VERBOSE(3) << "Pythia Gun: " << nExtra_ << " extra scattering(s), mu = "
             << extraMu_ << " (N_coll " << ncoll << ")";
}
