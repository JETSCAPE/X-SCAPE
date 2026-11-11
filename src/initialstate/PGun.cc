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
// Parton Gun Test

#include "PGun.h"
#include "JetScapeParticles.h"
#include "Pythia8/Pythia.h"

#include <sstream>
#include <stdexcept>
#include <string>

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<PGun> PGun::reg("PGun");

Pythia8::Pythia PGun::InternalHelperPythia("IntentionallyEmpty", false);

PGun::PGun() : HardProcess() {
  fixed_pT = 0;
  parID = 21;
  SetId("PGun");
  VERBOSE(8);
}

PGun::~PGun() { VERBOSE(8); }

void PGun::InitTask() {
  JSDEBUG << "Initialize PGun Brick (Test) ...";
  VERBOSE(8);

  std::string s = GetXMLElementText({"Hard", "PGun", "name"});
  JSDEBUG << s << " to be initialized ...";

  fixed_pT = GetXMLElementDouble({"Hard", "PGun", "pT"});
  JSDEBUG << s << " with fixed pT = " << fixed_pT;
  JSINFO << "Parton Gun with fixed pT = " << fixed_pT;

  parID = GetXMLElementDouble({"Hard", "PGun", "parID"});
  JSINFO << "Parton Gun with parID = " << parID;

  // Optional: the parton's rapidity and phi (both 0 before these existed).
  if (!GetXMLElementText({"Hard", "PGun", "rapidity"}, false).empty())
    rapidity_ = GetXMLElementDouble({"Hard", "PGun", "rapidity"});
  if (!GetXMLElementText({"Hard", "PGun", "phi"}, false).empty())
    phi_ = GetXMLElementDouble({"Hard", "PGun", "phi"});

  // Optional: several (pT, rapidity) points, event i uses point i mod K.
  points_.clear();
  std::istringstream in(GetXMLElementText({"Hard", "PGun", "points"}, false));
  std::vector<double> v;
  std::string tok;
  while (in >> tok) {
    try {
      size_t used = 0;
      v.push_back(std::stod(tok, &used));
      if (used != tok.size())
        throw std::invalid_argument(tok);
    } catch (const std::exception &) {
      throw std::runtime_error("PGun: <points> takes numbers 'pT y pT y ...', got '" +
                               tok + "'");
    }
  }
  if (v.size() % 2)
    throw std::runtime_error("PGun: <points> needs pairs 'pT y', got " +
                             std::to_string(v.size()) + " numbers");
  for (size_t j = 0; j < v.size(); j += 2) {
    if (!(v[j] > 0))
      throw std::runtime_error("PGun: every <points> pT must be > 0");
    points_.emplace_back(v[j], v[j + 1]);
  }
  if (points_.empty()) {
    points_.emplace_back(fixed_pT, rapidity_);
  } else {
    const int K = GetNPoints();
    std::string reuse = GetXMLElementText({"setReuseHydro"}, false);
    int nReuse = GetXMLElementInt({"nReuseHydro"}, false);
    if (reuse.find("true") != std::string::npos && nReuse % K != 0)
      throw std::runtime_error(
          "PGun: nReuseHydro = " + std::to_string(nReuse) +
          " is not a multiple of the " + std::to_string(K) +
          " <points>, so the hydro events would not get the same number of "
          "events per point");
    std::ostringstream msg;
    for (const auto &p : points_)
      msg << " (" << p.first << ", " << p.second << ")";
    JSINFO << "Parton Gun with " << K << " (pT, y) points (event i uses point i mod "
           << K << "):" << msg.str();
  }
  JSINFO << "Parton Gun: rapidity = " << points_[0].second << ", phi = " << phi_
         << ", from the origin";
}

void PGun::ExecuteTask() {
  VERBOSE(2) << "Run Hard Process : " << GetId() << " ...";

  double p[4], xLoc[4];

  double pT, rapidity, phi, pseudorapidity;
  double eta_cut = 1.0;
  double tempRand;
  const double maxN = 1.0 * RAND_MAX;
  const double PI = 3.1415926;

  double ppx, ppy, ppz, pp0, mass;

  // for (int i=0;i<1;i++)
  //    {
  //      tempRand = rand()/maxN;
  //    if(tempRand < 0.25) parID = 21;
  //    else if(tempRand < 0.50) parID = 1;
  //    else if(tempRand < 0.75) parID = 2;
  //    else parID = 3;
  //    if (parID != 21) {
  //	 tempRand = rand()/maxN;
  //	 if(tempRand < 0.50) parID = -parID;
  //      }
  //     mass = 0.0;
  mass = InternalHelperPythia.particleData.m0(parID);
  // JSINFO << BOLDYELLOW << " Mass = " << mass ;
  // This event's (pT, rapidity): point i mod K (one point without <points>).
  const int K = GetNPoints();
  activePoint_ = K > 1 ? GetCurrentEvent() % K : 0;
  pT = points_[activePoint_].first;

  phi = 2.0 * PI * (rand() / maxN);  // kept: the rand() sequence stays the same
  rapidity = points_[activePoint_].second;
  phi = phi_;

  p[1] = pT * cos(phi);
  p[2] = pT * sin(phi);
  p[3] = sqrt(pT * pT + mass * mass) * sinh(rapidity);
  p[0] = sqrt(pT * pT + mass * mass) * cosh(rapidity);

  double p_abs = std::sqrt(pT * pT + p[3] * p[3]);
  if (std::abs(p_abs - p[3]) > rounding_error) {
    pseudorapidity = 0.5 * std::log((p_abs + p[3]) / (p_abs - p[3]));
  } else {
    JSWARN << "Particle in PGun has infinite pseudorapidity.";
  }

  // Roll for a starting point
  // See:
  // https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  for (int i = 0; i <= 3; i++) {
    xLoc[i] = 0.0;
  };

  if (!ini) {
    VERBOSE(1)
        << "No initial state module, setting the starting location to 0. ";
  } else {
    double t, x, y, z;
    ini->SampleABinaryCollisionPoint(t, x, y, z);
    xLoc[0] = t;
    xLoc[1] = x;
    xLoc[2] = y;
    xLoc[3] = z;
  }

  xLoc[0] = 0.0;
  xLoc[1] = 0.0;
  xLoc[2] = 0.0;
  xLoc[3] = 0.0;

  auto ptn =
      make_shared<Parton>(0, parID, 0, pT, pseudorapidity, phi, p[0], xLoc);
  ptn->set_color((parID > 0) ? 100 : 0);
  ptn->set_anti_color(((parID > 0) || (parID == 21)) ? 0 : 101);
  ptn->set_max_color(102);
  AddParton(ptn);

  VERBOSE(8) << GetNHardPartons();
}
