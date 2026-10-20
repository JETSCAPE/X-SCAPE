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

#include "Jewel.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>

#include "FluidDynamics.h"
#include "JetScapeLogger.h"
#include "JetScapeModuleMutex.h"
#include "JetScapeTaskSupport.h"
#include "jewel/jetscape_api.hpp"

using namespace Jetscape;
using namespace std;

RegisterJetScapeModule<Jewel> Jewel::reg("Jewel");

bool Jewel::jewel_initialised_ = false;

namespace {

// JEWEL works on global PYTHIA 6 common blocks: one shower at a time
std::mutex jewel_lock;
// module (clone) whose parton is being showered: its fluid signal is connected
Jewel *current_module = nullptr;

// JEWEL simulates the whole final-state evolution and cannot share partons
// with another energy loss module
class JewelMutex : public JetScapeModuleMutex {
 public:
  bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules) {
    for (auto module : modules) {
      if (module->GetId() != "Jewel") return false;
    }
    return true;
  }
};

// optional XML element (the defaults are in jetscape_main.xml)
double xml_double(std::initializer_list<const char *> path, double def) {
  auto *e = JetScapeXML::Instance()->GetElement(path, false);
  double v = def;
  if (e) e->QueryDoubleText(&v);
  return v;
}

// optional XML text element; empty elements give "" (GetXMLElementText would
// construct a std::string from the null pointer tinyxml returns for them)
std::string xml_text(std::initializer_list<const char *> path) {
  auto *e = JetScapeXML::Instance()->GetElement(path, false);
  const char *t = e ? e->GetText() : nullptr;
  std::string v = t ? t : "";
  v.erase(0, v.find_first_not_of(" \t\n\r"));
  v.erase(v.find_last_not_of(" \t\n\r") + 1);
  return v;
}

}  // namespace

Jewel::Jewel() {
  SetId("Jewel");
  VERBOSE(8);
  SetMutex(make_shared<JewelMutex>());
}

Jewel::~Jewel() { VERBOSE(8); }

void Jewel::InitTask() {
  JSINFO << "Initialize JEWEL ...";

  qmax_factor_ = GetXMLElementDouble({"Eloss", "Jewel", "qmax_factor"});

  if (jewel_initialised_) return;

  jewel::api::Settings s;
  s.param_file = xml_text({"Eloss", "Jewel", "param_file"});
  s.run_param_file = GetXMLElementText({"Eloss", "Jewel", "run_param_file"});
  s.seed = static_cast<long>(JetScapeTaskSupport::GetRandomSeed() % 900000000u);

  // parameters that JETSCAPE knows (they override the JEWEL parameter file)
  const double ptmax = xml_double({"Hard", "PythiaGun", "pTHatMax"}, -1.0);
  const double ecm = xml_double({"Hard", "PythiaGun", "eCM"}, 5020.0);
  std::ostringstream os;
  os << "SQRTS " << ecm;
  s.extra_params.push_back(os.str());
  os.str("");
  os << "PTMAX " << ptmax;
  s.extra_params.push_back(os.str());
  s.extra_params.push_back(
      "LOGFILE " + GetXMLElementText({"Eloss", "Jewel", "log_file"}));
  s.extra_params.push_back(
      std::string("KEEPRECOILS ") +
      (GetXMLElementInt({"Eloss", "Jewel", "keep_recoils"}) ? "T" : "F"));

  jewel::api::MediumSettings &m = s.medium;
  m.vacuum = GetXMLElementInt({"Eloss", "Jewel", "in_vac"}) != 0;
  m.tmax = GetXMLElementDouble({"Eloss", "Jewel", "Tmax"});
  m.tc = GetXMLElementDouble({"Eloss", "Jewel", "Tc"});
  m.nf = GetXMLElementInt({"Eloss", "Jewel", "nf"});
  m.mdfactor = GetXMLElementDouble({"Eloss", "Jewel", "mdfactor"});
  m.mdscalefac = GetXMLElementDouble({"Eloss", "Jewel", "mdscalefac"});
  m.pre_tau0 = GetXMLElementInt({"Eloss", "Jewel", "pre_tau0"});
  m.fluxmax = GetXMLElementDouble({"Eloss", "Jewel", "fluxmax"});
  m.etamax = GetXMLElementDouble({"Eloss", "Jewel", "etamax"});
  m.tau0 = GetXMLElementDouble({"Eloss", "tStart"});
  // medium interactions end at the end of the JETSCAPE energy-loss time
  // (e.g. the brick length) unless set explicitly
  m.lifetime = GetXMLElementDouble({"Eloss", "Jewel", "lifetime"});
  if (m.lifetime <= 0.0) m.lifetime = GetXMLElementDouble({"Eloss", "maxT"});

  // JEWEL reuses lookup tables found on disk, but they depend on the JEWEL
  // parameters and on the medium bounds (Tmax, Tc, Debye mass): name them after
  // a hash of everything that enters, so that changed settings get new tables
  {
    std::ostringstream key;
    if (!s.param_file.empty()) key << std::ifstream(s.param_file).rdbuf();
    for (const auto &l : s.extra_params) key << l << '\n';
    key << m.tmax << ' ' << m.tc << ' ' << m.nf << ' ' << m.mdfactor << ' '
        << m.mdscalefac << ' ' << m.etamax;
    std::uint64_t h = 1469598103934665603ull;  // FNV-1a
    for (unsigned char c : key.str()) h = (h ^ c) * 1099511628211ull;
    std::ostringstream tag;
    tag << std::hex << h;
    for (const char *t : {"splitint", "pdfs", "xsecs"}) {
      std::string label = std::string(t) == "splitint" ? "SPLITINTFILE"
                          : std::string(t) == "pdfs"   ? "PDFFILE"
                                                       : "XSECFILE";
      s.extra_params.push_back(label + " jewel_" + t + "_" + tag.str() +
                               ".dat");
    }
  }

  JSINFO << "JEWEL: parameter file '" << s.param_file << "', in_vac = "
         << m.vacuum << ", Tmax = " << m.tmax << " GeV, Tc = " << m.tc
         << " GeV, lifetime = " << m.lifetime << " fm, qmax = "
         << qmax_factor_ << " pT";

  // the medium is read from the JETSCAPE fluid in the same way as in LBT
  auto query = [](double t, double x, double y, double z) {
    jewel::api::MediumCell c;
    if (!current_module) return c;
    std::unique_ptr<FluidCellInfo> cell;
    // points outside the fluid's grid or lifetime: no medium (some fluid
    // modules throw there instead of returning an empty cell)
    try {
      current_module->GetHydroCellSignal(t, x, y, z, cell);
    } catch (const std::exception &) {
      return c;
    }
    if (cell) {
      c.temperature = cell->temperature;
      c.vx = cell->vx;
      c.vy = cell->vy;
      c.vz = cell->vz;
    }
    return c;
  };
  {
    std::lock_guard<std::mutex> lock(jewel_lock);
    jewel::api::initialise(s, query);
  }
  jewel_initialised_ = true;
  JSINFO << "JEWEL initialised (log in "
         << GetXMLElementText({"Eloss", "Jewel", "log_file"}) << ")";
}

void Jewel::DoEnergyLoss(double deltaT, double time, double Q2,
                         vector<Parton> &pIn, vector<Parton> &pOut) {
  for (int i = 0; i < pIn.size(); i++) {
    // photons pass through (as in LBT)
    if (pIn[i].pid() == photonid) {
      if (pIn[i].pstat() != 22) {
        pIn[i].set_stat(22);
        pOut.push_back(pIn[i]);
      }
      return;
    }
    // partons made by JEWEL (shower, recoils, holes) are final
    if (pIn[i].has_user_info<JewelUserInfo>()) return;
    // holes of another module or partons handed to the liquefier
    if (pIn[i].pstat() < 0) return;

    TakeResponsibilityFor(pIn[i]);

    // the cell query of the medium calls GetHydroCellSignal of this clone
    double hydro_tau0 = -1.0;
    GetHydroTau0Signal(hydro_tau0);

    const double p[4] = {pIn[i].e(), pIn[i].px(), pIn[i].py(), pIn[i].pz()};
    const double x[4] = {pIn[i].x_in().t(), pIn[i].x_in().x(),
                         pIn[i].x_in().y(), pIn[i].x_in().z()};
    const double qmax = qmax_factor_ * pIn[i].pt();

    vector<jewel::api::OutParton> out;
    bool ok;
    {
      std::lock_guard<std::mutex> lock(jewel_lock);
      current_module = this;
      jewel::api::next_event(hydro_tau0);
      ok = jewel::api::shower(pIn[i].pid(), p, x, qmax, out);
      current_module = nullptr;
    }
    if (!ok) {
      // JEWEL discarded the shower: keep the parton unmodified
      pOut.push_back(pIn[i]);
      pOut.back().set_user_info(new JewelUserInfo());
      return;
    }

    // colours: JEWEL numbers them from 601 (the incoming parton), shift them
    // into the colour range that the hard process reserved for this parton
    const unsigned int base = pIn[i].max_color();
    const int id_in = pIn[i].pid();
    auto colour = [&](int c) -> unsigned int {
      if (c <= 0) return 0;
      if (c == 601) return (id_in < 0 && id_in != 21) ? pIn[i].anti_color()
                                                      : pIn[i].color();
      if (c == 602 && id_in == 21) return pIn[i].anti_color();
      return base + static_cast<unsigned int>(c - 602);
    };

    unsigned int maxcol = base;
    for (const auto &o : out) {
      double po[4] = {o.p[0], o.p[1], o.p[2], o.p[3]};
      double xo[4] = {o.x[0], o.x[1], o.x[2], o.x[3]};
      pOut.push_back(Parton(0, o.id, o.status, FourVector(po), FourVector(xo)));
      Parton &np = pOut.back();
      np.set_user_info(new JewelUserInfo());
      if (o.status >= 0) {
        np.set_color(colour(o.col));
        np.set_anti_color(colour(o.acol));
        maxcol = std::max({maxcol, np.color(), np.anti_color()});
      }
      np.set_form_time(0.0);
    }
    for (auto &np : pOut) np.set_max_color(maxcol);
  }
}

void Jewel::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  auto f = w.lock();
  if (!f) return;
  const auto &st = jewel::api::stats();
  std::ostringstream os;
  os << "JEWEL showers " << st.n_showers << " discarded " << st.n_discarded
     << " scatterings " << st.n_scatterings << " cells above Tmax "
     << st.n_above_tmax << " rate queries above fluxmax " << st.n_above_fluxmax
     << " (max gamma(1-v.beta) " << st.max_flux << ") holes dropped "
     << st.n_holes_dropped;
  f->WriteComment(os.str());
}

void Jewel::FinishTask() {
  static bool done = false;
  if (done) return;
  done = true;
  const auto &st = jewel::api::stats();
  JSINFO << "JEWEL: " << st.n_showers << " showers (" << st.n_discarded
         << " discarded), " << st.n_scatterings << " stored scatterings, "
         << st.n_holes_dropped << " holes dropped (scattering undone), "
         << st.n_above_tmax << " cells above Tmax, " << st.n_above_fluxmax
         << " rate queries above fluxmax (max gamma(1-v.beta) = " << st.max_flux
         << ")";
}
