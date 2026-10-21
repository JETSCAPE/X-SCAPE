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

#include "JetScapeSignalManager.h"
#include "iColoredHadronization.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include "tinyxml2.h"
#include "MCGlauberGenStringWrapper.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>

using namespace Jetscape;
using namespace Pythia8;

// Register the module with the base class
RegisterJetScapeModule<iColoredHadronization> iColoredHadronization::reg(
    "iColoredHadronization");

Pythia8::Pythia iColoredHadronization::pythia("IntentionallyEmpty", false);

iColoredHadronization::iColoredHadronization() {
  SetId("ISRMyHadroTest");
  VERBOSE(8);
}

iColoredHadronization::~iColoredHadronization() { VERBOSE(8); }

void iColoredHadronization::InitTask() {
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
  // Optionally read in additional settings
  std::stringstream lines;
  lines << GetXMLElementText({"JetHadronization", "LinesToRead"}, false);
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue;  // skip empty lines
    JSINFO << "Also reading in: " << s;
    pythia.readString(s);
  }

  reseed_per_event_ =
      GetXMLElementInt({"JetHadronization", "reseed_per_event"}, false) == 1;
  if (reseed_per_event_) {
    JSINFO << "iColoredHadronization: reseeding Pythia every event";
  }
  store_input_ =
      GetXMLElementInt({"JetHadronization", "store_input"}, false) == 1;
  // <colour_duplicates>: fail (default; the event throws) or relabel
  std::string dup =
      GetXMLElementText({"JetHadronization", "colour_duplicates"}, false);
  dup.erase(0, dup.find_first_not_of(" \t\n"));
  dup.erase(dup.find_last_not_of(" \t\n") + 1);
  if (dup == "relabel")
    relabel_duplicates_ = true;
  else if (!dup.empty() && dup != "fail")
    throw std::runtime_error(
        "iColoredHadronization: <colour_duplicates> takes fail or relabel, "
        "not '" + dup + "'");
  if (relabel_duplicates_)
    JSINFO << "iColoredHadronization: duplicate colour tags are relabelled";
  // <colour_closure>: how the colour lines that end in the medium are closed
  std::string closure =
      GetXMLElementText({"JetHadronization", "colour_closure"}, false);
  closure.erase(0, closure.find_first_not_of(" \t\n"));
  closure.erase(closure.find_last_not_of(" \t\n") + 1);
  SetColourClosure(closure.empty() ? "reconnect" : closure);
  JSINFO << "iColoredHadronization: colour closure " << GetColourClosure();

  pythia.init();
}

void iColoredHadronization::SetColourClosure(const std::string &mode) {
  if (mode == "reconnect")
    closure_ = Closure::Reconnect;
  else if (mode == "reconnect_global")
    closure_ = Closure::ReconnectGlobal;
  else if (mode == "thermal")
    closure_ = Closure::Thermal;
  else if (mode == "fake")
    closure_ = Closure::Fake;
  else if (mode == "none")
    closure_ = Closure::None;
  else
    throw std::runtime_error(
        "iColoredHadronization: <colour_closure> takes reconnect, "
        "reconnect_global, thermal, fake or none, not '" + mode + "'");
}

std::string iColoredHadronization::GetColourClosure() const {
  switch (closure_) {
  case Closure::Reconnect:
    return "reconnect";
  case Closure::ReconnectGlobal:
    return "reconnect_global";
  case Closure::Thermal:
    return "thermal";
  case Closure::Fake:
    return "fake";
  case Closure::None:
    return "none";
  }
  return "?";
}

void iColoredHadronization::Reseed(unsigned int seed) {
  // Pythia's Rndm takes 0 < seed <= 900000000 (0 would mean "from the clock")
  pythia.rndm.init(static_cast<int>(seed % 900000000u) + 1);
  last_seed_ = seed;
}

void iColoredHadronization::ReseedForEvent() {
  if (has_next_seed_) {
    Reseed(next_seed_);
    has_next_seed_ = false;
  } else if (reseed_per_event_) {
    Reseed((*GetMt19937Generator())());
  }
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
  ReseedForEvent();
  // auto Particles = shower->GetFinalPartons();
  VERBOSE(2) << BOLDYELLOW << "iColored Hadronization partons ";
  std::vector<ShowerParton> partons;
  for (unsigned int ishower = 0; ishower < shower.size(); ++ishower) {
    JSDEBUG << "&&&&&&&&&&&&&&&&&&& there are " << shower.at(ishower).size()
            << " partons in the shower number " << ishower;
    for (unsigned int ipart = 0; ipart < shower.at(ishower).size(); ++ipart) {
      const auto &q = shower.at(ishower).at(ipart);
      double onshellE =
          pow(pow(q->px(), 2) + pow(q->py(), 2) + pow(q->pz(), 2), 0.5);

      if (q->pid() == 22) {
        VERBOSE(1) << BOLDYELLOW
                   << " photon found in colored hadronization with ";
        VERBOSE(1) << BOLDYELLOW << "px = " << q->px();
        // cin >> blurb;
      }

      partons.push_back({int(ishower), q->pid(), q->pstat(), int(q->color()),
                         int(q->anti_color()), q->px(), q->py(), q->pz(),
                         onshellE});

      VERBOSE(2) << BOLDYELLOW << " pid = " << q->pid()
                 << " plabel = " << q->plabel() << " pstat = " << q->pstat()
                 << " col = " << q->color() << " acol = " << q->anti_color();
    }
  }

  auto ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  auto Hard = JetScapeSignalManager::Instance()->GetHardProcessPointer().lock();

  auto MCGsecond = std::dynamic_pointer_cast<MCGlauberGenStringWrapper>(
      Hard->GetTaskList()[1]);
  auto Remnants = Hard->GetRemnants();
  if (2 * ini->pTHat.size() != Remnants.size()) {
    throw std::runtime_error(
        "Not enough remnants = " + std::to_string(Remnants.size()) +
        " Scattering = " + std::to_string(ini->pTHat.size()));
  }
  double NHardScatterings = double(ini->pTHat.size());

  for (unsigned int ipart = 0; ipart < Remnants.size(); ++ipart) {
    auto Rem = Remnants[ipart];
    double Pz, Px, Py, En;
    if (Rem.pz() >= 0) {
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

    // std::cout << "Px = " << Px << " Py = " << Py << " Pz = " << Pz << " En =
    // "<< En << " " << NHardScatterings <<  std::endl;

    double onshellE =
        pow(pow(Rem.px() + Px, 2) + pow(Rem.py() + Py, 2) + pow(Pz, 2), 0.5);
    partons.push_back({-1, Rem.pid(), Rem.pstat(), int(Rem.color()),
                       int(Rem.anti_color()), Rem.px() + Px, Rem.py() + Py, Pz,
                       onshellE});
    VERBOSE(2) << BOLDYELLOW << " pid = " << Rem.pid()
               << " plabel = " << Rem.plabel() << " pstat = " << Rem.pstat()
               << " col = " << Rem.color() << " acol = " << Rem.anti_color();
  }

  HadronizeList(partons, hOut);
  VERBOSE(2) << "There are " << hOut.size() << " Hadrons and " << pOut.size()
             << " partons after Hadronization";
  shower.clear();
}

void iColoredHadronization::HadronizePartons(
    const std::vector<ShowerParton> &partons,
    vector<shared_ptr<Hadron>> &hOut) {
  ReseedForEvent();
  HadronizeList(partons, hOut);
}

void iColoredHadronization::HadronizeList(
    const std::vector<ShowerParton> &partons,
    vector<shared_ptr<Hadron>> &hOut) {
  Event &event = pythia.event;
  event.reset();
  audit_ = ColourAudit();
  // the shower of every record entry (entry 0 is the system line; -1: a beam
  // remnant), for the closure
  std::vector<int> shower_of(1, -2);
  // Holes (pstat -1: the thermal partons a scattering took out of the medium)
  // and the liquefier's partons (-11, -13, -17) are not string partons. The
  // holes' colour tags tell the closure which lines end in the medium; their
  // energy is dropped and recorded.
  // The beam remnants (shower -1, whatever their status) always go to Pythia.
  std::vector<ShowerParton> holes;
  for (const auto &q : partons) {
    if (q.shower >= 0 && q.pstat == -1) {
      ++audit_.holes;
      audit_.hole_energy += q.e;
      holes.push_back(q);
      continue;
    }
    if (q.shower >= 0 && q.pstat < 0) {
      ++audit_.skipped;
      continue;
    }
    if (q.shower < 0)
      ++audit_.remnants;
    else
      ++audit_.partons;
    event.append(q.id, 23, q.col, q.acol, q.px, q.py, q.pz, q.e);
    shower_of.push_back(q.shower);
  }

  // Pair the colour tags by counting: every tag must occur exactly once as a
  // colour and once as an anticolour. cols/acols get one entry per missing
  // partner.
  std::vector<int> cols;
  std::vector<int> acols;
  audit_.duplicate_tags = CountColourTags(event, cols, acols);
  if (audit_.duplicate_tags > 0) {
    if (!relabel_duplicates_)
      throw std::runtime_error(
          "iColoredHadronization: " + std::to_string(audit_.duplicate_tags) +
          " duplicate colour tags (a tag used more than once as a colour or "
          "as an anticolour); <colour_duplicates>relabel</colour_duplicates> "
          "relabels them");
    audit_.relabelled_tags = RelabelDuplicateTags(event);
    cols.clear();
    acols.clear();
    CountColourTags(event, cols, acols);
  }
  audit_.unpaired_col = cols.size();
  audit_.unpaired_acol = acols.size();

  // Close the lines that end in the medium (or are broken otherwise),
  // <colour_closure>. A copy of the record before it, for the fallback.
  const bool has_open = !cols.empty() || !acols.empty();
  Event before;
  std::vector<int> shower_before;
  if (has_open) {
    before = event;
    shower_before = shower_of;
    CloseColourLines(event, shower_of, holes, closure_);
  }
  cols.clear();
  acols.clear();
  CountColourTags(event, cols, acols);
  if (cols.size() > 0 || acols.size() > 0) {
    for (auto col : cols)
      JSWARN << " col = " << col;
    for (auto col : acols)
      JSWARN << " acol = " << col;
    throw std::runtime_error("Unpaired colors sent to Pythia");
  }

  RecordInput();
  // A failed fragmentation leaves the input partons in the record as final
  // particles: never hand them on as hadrons. With a closure that added
  // partners or joined lines, try once more with soft end partons (B3).
  if (!pythia.next()) {
    const bool retry =
        has_open && closure_ != Closure::Fake && closure_ != Closure::None;
    if (!retry)
      throw std::runtime_error(
          "iColoredHadronization: Pythia could not hadronize the event (" +
          std::to_string(event.size() - 1) + " partons in the record)");
    ++audit_.retries;
    event = before;
    shower_of = shower_before;
    CloseColourLines(event, shower_of, holes, Closure::Fake);
    RecordInput();
    if (!pythia.next())
      throw std::runtime_error(
          "iColoredHadronization: Pythia could not hadronize the event, also "
          "with soft end partons (" +
          std::to_string(event.size() - 1) + " partons in the record)");
  }
  // event.list();

  CollectHadrons(hOut);
}

void iColoredHadronization::RecordInput() {
  if (!store_input_)
    return;
  // entry 0 is the system line event.reset() adds
  Event &event = pythia.event;
  last_input_.clear();
  for (int i = 1; i < event.size(); ++i) {
    last_input_.push_back({event[i].id(), event[i].col(), event[i].acol(),
                           event[i].px(), event[i].py(), event[i].pz(),
                           event[i].e(), event[i].m()});
  }
}

void iColoredHadronization::OpenEnds(const Event &event,
                                     std::map<int, int> &col_at,
                                     std::map<int, int> &acol_at) {
  std::map<int, int> ncol, nacol;
  for (int i = 1; i < event.size(); ++i) {
    if (event[i].col() != 0)
      ++ncol[event[i].col()];
    if (event[i].acol() != 0)
      ++nacol[event[i].acol()];
  }
  col_at.clear();
  acol_at.clear();
  for (int i = 1; i < event.size(); ++i) {
    const int c = event[i].col(), a = event[i].acol();
    if (c != 0 && nacol.count(c) == 0)
      col_at[c] = i;
    if (a != 0 && ncol.count(a) == 0)
      acol_at[a] = i;
  }
}

int iColoredHadronization::Reconnect(Event &event,
                                     const std::vector<int> &shower_of,
                                     bool within_shower) {
  // B2: join open colour ends with open anticolour ends, the pair with the
  // smallest invariant mass (the shortest string) first. The anticolour end
  // takes the colour end's tag. Never one gluon with itself.
  std::map<int, int> col_at, acol_at;
  OpenEnds(event, col_at, acol_at);
  struct Candidate {
    double m2;
    int ic, ia, tag;
  };
  std::vector<Candidate> cand;
  for (const auto &c : col_at)
    for (const auto &a : acol_at) {
      if (c.second == a.second)
        continue;
      if (within_shower && shower_of[c.second] != shower_of[a.second])
        continue;
      cand.push_back({(event[c.second].p() + event[a.second].p()).m2Calc(),
                      c.second, a.second, c.first});
    }
  std::sort(cand.begin(), cand.end(),
            [](const Candidate &x, const Candidate &y) {
              if (x.m2 != y.m2)
                return x.m2 < y.m2;
              return x.ic != y.ic ? x.ic < y.ic : x.ia < y.ia;
            });
  std::set<int> used_col, used_acol;  // record entries whose end is joined
  int joined = 0;
  for (const auto &k : cand) {
    if (used_col.count(k.ic) || used_acol.count(k.ia))
      continue;
    event[k.ia].acol(k.tag);
    used_col.insert(k.ic);
    used_acol.insert(k.ia);
    ++joined;
  }
  return joined;
}

void iColoredHadronization::InsertPartners(
    Event &event, std::vector<int> &shower_of,
    const std::vector<ShowerParton> &holes) {
  // B1: the medium partner of a hole, a parton of the conjugate colour charge
  // with the hole's tags swapped (the hole carries the incoming thermal
  // parton's tags) and the hole's momentum, for each hole whose lines are
  // still open. A gluon hole with one line open gets a light (anti)quark.
  std::map<int, int> col_at, acol_at;
  OpenEnds(event, col_at, acol_at);
  for (const auto &h : holes) {
    const bool c_open = h.col != 0 && col_at.count(h.col);
    const bool a_open = h.acol != 0 && acol_at.count(h.acol);
    if (!c_open && !a_open)
      continue;
    int id;
    if (c_open && a_open)
      id = 21;
    else if (h.id == 21)
      id = c_open ? -1 : 1;
    else
      id = -h.id;
    const int col = a_open ? h.acol : 0;
    const int acol = c_open ? h.col : 0;
    const double e = std::sqrt(h.px * h.px + h.py * h.py + h.pz * h.pz);
    event.append(id, 23, col, acol, h.px, h.py, h.pz, e);
    shower_of.push_back(h.shower);
    if (c_open)
      col_at.erase(h.col);
    if (a_open)
      acol_at.erase(h.acol);
    ++audit_.partners;
    audit_.partner_energy += e;
  }
}

void iColoredHadronization::InsertFakes(Event &event,
                                        std::vector<int> &shower_of) {
  // B3: a soft (1 GeV) light antiquark (quark) along the beam, on the side of
  // the parton with the open colour (anticolour), for every open end left.
  std::map<int, int> col_at, acol_at;
  OpenEnds(event, col_at, acol_at);
  for (const auto &c : col_at) {
    const double pz = event[c.second].pz() >= 0 ? 1. : -1.;
    event.append(-1, 23, 0, c.first, 0., 0., pz, 1.);
    shower_of.push_back(shower_of[c.second]);
    ++audit_.fakes;
    audit_.fake_energy += 1.;
  }
  for (const auto &a : acol_at) {
    const double pz = event[a.second].pz() >= 0 ? 1. : -1.;
    event.append(1, 23, a.first, 0, 0., 0., pz, 1.);
    shower_of.push_back(shower_of[a.second]);
    ++audit_.fakes;
    audit_.fake_energy += 1.;
  }
}

void iColoredHadronization::CloseColourLines(
    Event &event, std::vector<int> &shower_of,
    const std::vector<ShowerParton> &holes, Closure mode) {
  switch (mode) {
  case Closure::None:
    return;
  case Closure::Fake:
    InsertFakes(event, shower_of);
    return;
  case Closure::Reconnect:
    audit_.reconnected += Reconnect(event, shower_of, true);
    break;
  case Closure::ReconnectGlobal:
    audit_.reconnected += Reconnect(event, shower_of, false);
    break;
  case Closure::Thermal:
    break;
  }
  InsertPartners(event, shower_of, holes);
  // what is still open (lines broken without a hole of their own): across
  // showers, then soft end partons; both counted
  audit_.reconnected_cross += Reconnect(event, shower_of, false);
  InsertFakes(event, shower_of);
}

void iColoredHadronization::HadronizeInput(
    const std::vector<InputParton> &input, vector<shared_ptr<Hadron>> &hOut) {
  ReseedForEvent();
  Event &event = pythia.event;
  event.reset();
  for (const auto &q : input) {
    event.append(q.id, 23, q.col, q.acol, q.px, q.py, q.pz, q.e, q.m);
  }
  if (!pythia.next())
    throw std::runtime_error(
        "iColoredHadronization: Pythia could not hadronize the stored input (" +
        std::to_string(input.size()) + " partons)");
  CollectHadrons(hOut);
}

int iColoredHadronization::CountColourTags(const Event &event,
                                           std::vector<int> &cols,
                                           std::vector<int> &acols) {
  std::map<int, int> ncol, nacol;
  for (int i = 1; i < event.size(); ++i) {  // entry 0 is the system line
    if (event[i].col() != 0)
      ++ncol[event[i].col()];
    if (event[i].acol() != 0)
      ++nacol[event[i].acol()];
  }
  int duplicates = 0;
  for (const auto &t : ncol) {
    auto a = nacol.find(t.first);
    const int na = a == nacol.end() ? 0 : a->second;
    duplicates += std::max(0, t.second - 1);
    for (int k = na; k < t.second; ++k)
      cols.push_back(t.first);
  }
  for (const auto &t : nacol) {
    auto c = ncol.find(t.first);
    const int nc = c == ncol.end() ? 0 : c->second;
    duplicates += std::max(0, t.second - 1);
    for (int k = nc; k < t.second; ++k)
      acols.push_back(t.first);
  }
  return duplicates;
}

int iColoredHadronization::RelabelDuplicateTags(Event &event) {
  // The k-th colour and the k-th anticolour (in record order) with the same
  // tag become a pair with a new tag, for k >= 2. Which colour belongs to which
  // anticolour is not known here, so this is a repair, counted in the audit.
  int max_tag = 0;
  for (int i = 1; i < event.size(); ++i)
    max_tag = std::max(max_tag, std::max(event[i].col(), event[i].acol()));
  std::map<int, int> seen_col, seen_acol;  // occurrences so far, per tag
  std::map<std::pair<int, int>, int> new_tag;  // (tag, k) -> new tag
  int relabelled = 0;
  for (int i = 1; i < event.size(); ++i) {
    const int c = event[i].col();
    if (c != 0 && ++seen_col[c] > 1) {
      auto key = std::make_pair(c, seen_col[c]);
      if (!new_tag.count(key))
        new_tag[key] = ++max_tag;
      event[i].col(new_tag[key]);
      ++relabelled;
    }
    const int a = event[i].acol();
    if (a != 0 && ++seen_acol[a] > 1) {
      auto key = std::make_pair(a, seen_acol[a]);
      if (!new_tag.count(key))
        new_tag[key] = ++max_tag;
      event[i].acol(new_tag[key]);
      ++relabelled;
    }
  }
  return relabelled;
}

void iColoredHadronization::CollectHadrons(vector<shared_ptr<Hadron>> &hOut) {
  Event &event = pythia.event;
  unsigned int ip = hOut.size();
  for (unsigned int i = 0; i < event.size(); ++i) {
    if (!event[i].isFinal())
      continue;
    // if ( !event[i].isHadron() )  continue;
    if (fabs(event[i].eta()) > 20)
      continue;  // To prevent "nan" from propagating, very rare though

    double x[4] = {0, 0, 0, 0};
    hOut.push_back(make_shared<Hadron>(ip, event[i].id(), event[i].status(),
                                       event[i].pT(), event[i].eta(),
                                       event[i].phi(), event[i].e(), x));
    ++ip;
  }
}
