// Weighted spectra and jet substructure from JETSCAPE final-state hadron
// files (or named pipes), using fjcore anti-kT.
//
//   jetana out.spec file1_final_state_hadrons.dat [file2 ...]
//
// Every histogram stores sum(w*n) ("_h") and sum((w*n)^2) ("_h2") over events,
// together with sum(w), so files can be summed and yields are per unit weight.
//
// Spectra
//   fulljet : anti-kT R=0.4, all hadrons, |y_jet| < 2.8           (ATLAS)
//   chjet   : anti-kT R=0.4, charged hadrons, |eta_jet| < 0.5     (ALICE)
//   hadron  : charged hadrons, |eta| < 1                          (CMS)
//   sphenix : anti-kT R=0.4, all hadrons, |eta_jet| < 0.7         (sPHENIX-like)
// Radius dependence (full jets, holes subtracted within each jet's own R)
//   rs_specR02/R04/R06 : |eta_jet| < 1.1 - R, pT 15-60 GeV   (sPHENIX-like)
//   ls_specR02/R04/R06 : |y_jet| < 2.0, pT 30-316 GeV
// Substructure of full R=0.4 jets, in jet-pT bins of two selections
//   rs (RHIC): |eta_jet| < 0.7, pT 15-20, 20-30, 30-60 GeV
//   ls (LHC) : |eta_jet| < 1.6, pT 40-60, 60-80, 80-120, 120-200 GeV
//   <sel>_njet     : number of jets per jet-pT bin (per-jet normalisation)
//   <sel>_eec<k>   : sum_{i<j} pT_i pT_j / pT_jet^2 in log bins of dR_ij
//   <sel>_mass<k>  : jet mass
//   <sel>_ang05<k>, _ang1<k>, _ang2<k> : angularities
//                    lambda_a = sum_i pT_i (dR_i / R)^a / pT_jet
// Substructure uses charged constituents with pT > 1 GeV.
// Data-style observables (LHC), added for the lambda_t = 2.4 reruns; cuts follow
// the measurements approximately and can be refined from the kept event files:
//   xj_lead, xj<k>  : dijet x_J = pT2/pT1 (ATLAS-like): R = 0.4 full jets with
//                     |y| < 2.1, leading pT1 in 100-126-158-200-251-316 GeV,
//                     subleading = next jet with |y| < 2.1, pT2 > 25 GeV and
//                     dphi > 7pi/8 to the leading jet; xj_lead counts leading
//                     jets per pT1 bin, xj_pair the accepted pairs
//   ff_njet, ff_z<k>, ff_pt<k> : fragmentation (ATLAS-like): R = 0.4 full jets,
//                     |y| < 2.1, pT 126-158-200-251-316 GeV; charged particles
//                     pT > 1 GeV within dR < 0.4 of the axis (holes negative);
//                     z = pT cos(dR) / pT_jet
//   js_njet, js_rho : jet shape (CMS-like): R = 0.4 full jets pT > 120 GeV,
//                     0.3 < |eta| < 1.6; charged particles pT > 0.7 GeV within
//                     r < 1 of the axis (holes negative), pT/pT_jet summed in
//                     annuli of 0.05 (normalise to unit integral in the analysis)
//   chm_njet, chm_mass<k> : charged-jet mass (ALICE-like): anti-kT R = 0.4 on
//                     charged particles with pT > 0.15 GeV, |eta| < 0.5, charged
//                     jet pT 60-80-100-120 GeV, holes subtracted
// Hadrons with negative status (LBT holes) are not clustered; their
// 4-momenta within dR < R of a jet axis are subtracted from the jet
// (JETSCAPE hole subtraction), they enter substructure sums with negative
// weight, and they count negatively in the hadron spectrum.
#include "fjcore.hh"

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using fjcore::PseudoJet;

namespace {
const double R = 0.4;  // substructure and the named spectra above
const std::vector<std::pair<std::string, double>> kRadii = {
    {"R02", 0.2}, {"R04", 0.4}, {"R06", 0.6}};
const double kTrackPtMin = 1.0;
const std::set<int> kCharged = {211, 321, 2212, 3222, 3112, 3312, 3334, 11, 13};

struct Hist {
  std::vector<double> edges, h, h2, ev;  // ev: current-event accumulator
  Hist() {}
  explicit Hist(std::vector<double> e)
      : edges(e), h(e.size() - 1, 0.), h2(e.size() - 1, 0.), ev(e.size() - 1, 0.) {}
  void Fill(double x, double v = 1.) {
    if (x < edges.front() || x >= edges.back())
      return;
    for (size_t i = 0; i + 1 < edges.size(); i++)
      if (x < edges[i + 1]) {
        ev[i] += v;
        return;
      }
  }
  void EndEvent(double w) {
    for (size_t i = 0; i < ev.size(); i++) {
      double c = w * ev[i];
      h[i] += c;
      h2[i] += c * c;
      ev[i] = 0.;
    }
  }
};

std::vector<double> Linear(double lo, double hi, int n) {
  std::vector<double> e;
  for (int i = 0; i <= n; i++)
    e.push_back(lo + (hi - lo) * i / n);
  return e;
}

std::vector<double> Log(double lo, double hi, int n) {
  std::vector<double> e;
  for (int i = 0; i <= n; i++)
    e.push_back(lo * std::pow(hi / lo, double(i) / n));
  return e;
}

struct Particle {
  PseudoJet p;
  bool charged, hole;
};

struct Selection {
  std::string name;
  double etaMax;
  std::vector<double> ptBins;
};

const std::vector<Selection> kSelections = {
    {"rs", 0.7, {15, 20, 30, 60}},
    {"ls", 1.6, {40, 60, 80, 120, 200}}};
}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: jetana out.spec files...\n";
    return 1;
  }
  std::map<std::string, Hist> H;
  H["fulljet"] = Hist({30, 40, 50, 63, 80, 100, 112, 125, 141, 158, 177, 199,
                       223, 251, 281, 316});
  H["chjet"] = Hist({30, 40, 50, 60, 70, 85, 100, 120, 140});
  H["hadron"] = Hist({10, 14, 20, 28, 40, 56, 80, 115, 160});
  H["sphenix"] = Hist({15, 20, 25, 30, 35, 40, 50, 60});
  for (auto &r : kRadii) {
    H["rs_spec" + r.first] = Hist({15, 20, 25, 30, 35, 40, 50, 60});
    H["ls_spec" + r.first] = Hist({30, 40, 50, 63, 80, 100, 126, 158, 200, 251, 316});
  }
  for (auto &s : kSelections) {
    H[s.name + "_njet"] = Hist(s.ptBins);
    for (size_t k = 0; k + 1 < s.ptBins.size(); k++) {
      std::string t = std::to_string(k);
      H[s.name + "_eec" + t] = Hist(Log(0.005, 0.4, 20));
      H[s.name + "_mass" + t] = Hist(Linear(0., 30., 30));
      H[s.name + "_ang05" + t] = Hist(Linear(0., 0.8, 32));
      H[s.name + "_ang1" + t] = Hist(Linear(0., 0.5, 25));
      H[s.name + "_ang2" + t] = Hist(Linear(0., 0.3, 30));
    }
  }
  const std::vector<double> kXjBins = {100, 126, 158, 200, 251, 316};
  const std::vector<double> kFfBins = {126, 158, 200, 251, 316};
  const std::vector<double> kChmBins = {60, 80, 100, 120};
  H["xj_lead"] = Hist(kXjBins);
  H["xj_pair"] = Hist(kXjBins);
  H["ff_njet"] = Hist(kFfBins);
  H["js_njet"] = Hist({120, 1000});
  H["js_rho"] = Hist(Linear(0., 1., 20));
  H["chm_njet"] = Hist(kChmBins);
  for (size_t k = 0; k + 1 < kXjBins.size(); k++)
    H["xj" + std::to_string(k)] = Hist(Linear(0.2, 1.0, 16));
  for (size_t k = 0; k + 1 < kFfBins.size(); k++) {
    H["ff_z" + std::to_string(k)] = Hist(Log(0.01, 1.0, 20));
    H["ff_pt" + std::to_string(k)] = Hist(Log(1.0, 400., 20));
  }
  for (size_t k = 0; k + 1 < kChmBins.size(); k++)
    H["chm_mass" + std::to_string(k)] = Hist(Linear(0., 30., 30));
  auto bin_of = [](const std::vector<double> &b, double x) {
    for (size_t i = 0; i + 1 < b.size(); i++)
      if (x >= b[i] && x < b[i + 1])
        return int(i);
    return -1;
  };
  double sumw = 0.;
  long nev = 0;

  auto process = [&](const std::vector<Particle> &parts, double w) {
    // Full jets (with constituent bookkeeping) and charged jets.
    std::vector<PseudoJet> in, chin, holes, chholes, jaxis;
    for (size_t i = 0; i < parts.size(); i++) {
      PseudoJet p = parts[i].p;
      p.set_user_index(i);
      if (parts[i].hole) {
        holes.push_back(p);
        if (parts[i].charged)
          chholes.push_back(p);
      } else {
        in.push_back(p);
        if (parts[i].charged)
          chin.push_back(p);
      }
    }
    fjcore::JetDefinition def(fjcore::antikt_algorithm, R);
    fjcore::ClusterSequence cs(in, def);
    std::vector<PseudoJet> subjets;  // hole-subtracted jets (axis of the jet)
    for (auto &j : cs.inclusive_jets(5.0)) {
      PseudoJet sub = j;
      std::vector<PseudoJet> jholes;
      for (auto &h : holes)
        if (j.delta_R(h) < R) {
          sub -= h;
          jholes.push_back(h);
        }
      double pt = sub.pt();
      subjets.push_back(PseudoJet(sub.px(), sub.py(), sub.pz(), sub.E()));
      jaxis.push_back(j);  // clustered axis, for the track-in-jet observables
      if (std::abs(j.rap()) < 2.8)
        H["fulljet"].Fill(pt);
      if (std::abs(j.eta()) < 0.7)
        H["sphenix"].Fill(pt);

      for (auto &s : kSelections) {
        if (std::abs(j.eta()) >= s.etaMax)
          continue;
        int k = -1;
        for (size_t b = 0; b + 1 < s.ptBins.size(); b++)
          if (pt >= s.ptBins[b] && pt < s.ptBins[b + 1])
            k = b;
        if (k < 0)
          continue;
        std::string t = std::to_string(k);
        H[s.name + "_njet"].Fill(pt);
        // charged constituents above threshold, holes with negative sign
        std::vector<std::pair<PseudoJet, double>> trk;
        for (auto &c : j.constituents())
          if (parts[c.user_index()].charged && c.pt() > kTrackPtMin)
            trk.push_back({c, 1.});
        for (auto &h : jholes)
          if (parts[h.user_index()].charged && h.pt() > kTrackPtMin)
            trk.push_back({h, -1.});
        for (size_t a = 0; a < trk.size(); a++)
          for (size_t b = a + 1; b < trk.size(); b++)
            H[s.name + "_eec" + t].Fill(
                trk[a].first.delta_R(trk[b].first),
                trk[a].second * trk[b].second * trk[a].first.pt() *
                    trk[b].first.pt() / (pt * pt));
        double l05 = 0., l1 = 0., l2 = 0.;
        for (auto &q : trk) {
          double x = q.first.delta_R(j) / R, z = q.second * q.first.pt() / pt;
          l05 += z * std::sqrt(x);
          l1 += z * x;
          l2 += z * x * x;
        }
        H[s.name + "_ang05" + t].Fill(l05);
        H[s.name + "_ang1" + t].Fill(l1);
        H[s.name + "_ang2" + t].Fill(l2);
        double m2 = sub.m2();
        H[s.name + "_mass" + t].Fill(m2 > 0. ? std::sqrt(m2) : 0.);
      }
    }
    // ---- data-style observables on the R = 0.4 full jets
    {
      // dijet x_J
      int i1 = -1, i2 = -1;
      for (size_t i = 0; i < subjets.size(); i++) {
        if (std::abs(jaxis[i].rap()) >= 2.1)
          continue;
        if (i1 < 0 || subjets[i].pt() > subjets[i1].pt()) {
          i2 = i1;
          i1 = i;
        } else if (i2 < 0 || subjets[i].pt() > subjets[i2].pt()) {
          i2 = i;
        }
      }
      if (i1 >= 0) {
        double pt1 = subjets[i1].pt();
        int k = bin_of(kXjBins, pt1);
        if (k >= 0) {
          H["xj_lead"].Fill(pt1);
          if (i2 >= 0 && subjets[i2].pt() > 25. &&
              std::abs(jaxis[i1].delta_phi_to(jaxis[i2])) > 7. * M_PI / 8.) {
            H["xj_pair"].Fill(pt1);
            H["xj" + std::to_string(k)].Fill(subjets[i2].pt() / pt1);
          }
        }
      }
      for (size_t i = 0; i < subjets.size(); i++) {
        const PseudoJet &ax = jaxis[i];
        double pt = subjets[i].pt();
        // fragmentation functions
        int k = std::abs(ax.rap()) < 2.1 ? bin_of(kFfBins, pt) : -1;
        if (k >= 0) {
          H["ff_njet"].Fill(pt);
          std::string t = std::to_string(k);
          for (auto &q : parts) {
            if (!q.charged || q.p.pt() <= 1.0)
              continue;
            double dr = q.p.delta_R(ax);
            if (dr >= R)
              continue;
            double v = q.hole ? -1. : 1.;
            H["ff_z" + t].Fill(q.p.pt() * std::cos(dr) / pt, v);
            H["ff_pt" + t].Fill(q.p.pt(), v);
          }
        }
        // jet shape rho(r)
        double aeta = std::abs(ax.eta());
        if (pt > 120. && aeta > 0.3 && aeta < 1.6) {
          H["js_njet"].Fill(pt);
          for (auto &q : parts) {
            if (!q.charged || q.p.pt() <= 0.7)
              continue;
            double dr = q.p.delta_R(ax);
            if (dr < 1.0)
              H["js_rho"].Fill(dr, (q.hole ? -1. : 1.) * q.p.pt() / pt);
          }
        }
      }
      // ALICE-like charged-jet mass: charged particles with pT > 0.15 GeV
      std::vector<PseudoJet> ch015, chh015;
      for (auto &p : chin)
        if (p.pt() > 0.15)
          ch015.push_back(p);
      for (auto &p : chholes)
        if (p.pt() > 0.15)
          chh015.push_back(p);
      fjcore::ClusterSequence cm(ch015, def);
      for (auto &j : cm.inclusive_jets(5.0)) {
        if (std::abs(j.eta()) >= 0.5)
          continue;
        PseudoJet sub = j;
        for (auto &h : chh015)
          if (j.delta_R(h) < R)
            sub -= h;
        int k = bin_of(kChmBins, sub.pt());
        if (k < 0)
          continue;
        H["chm_njet"].Fill(sub.pt());
        double m2 = sub.m2();
        H["chm_mass" + std::to_string(k)].Fill(m2 > 0. ? std::sqrt(m2) : 0.);
      }
    }
    for (auto &r : kRadii) {
      fjcore::JetDefinition defR(fjcore::antikt_algorithm, r.second);
      fjcore::ClusterSequence csR(in, defR);
      for (auto &j : csR.inclusive_jets(5.0)) {
        PseudoJet sub = j;
        for (auto &h : holes)
          if (j.delta_R(h) < r.second)
            sub -= h;
        if (std::abs(j.eta()) < 1.1 - r.second)
          H["rs_spec" + r.first].Fill(sub.pt());
        if (std::abs(j.rap()) < 2.0)
          H["ls_spec" + r.first].Fill(sub.pt());
      }
    }
    fjcore::ClusterSequence chcs(chin, def);
    for (auto &j : chcs.inclusive_jets(5.0)) {
      PseudoJet sub = j;
      for (auto &h : chholes)
        if (j.delta_R(h) < R)
          sub -= h;
      if (std::abs(j.eta()) < 0.5)
        H["chjet"].Fill(sub.pt());
    }
    for (auto &q : parts)
      if (q.charged && std::abs(q.p.eta()) < 1.)
        H["hadron"].Fill(q.p.pt(), q.hole ? -1. : 1.);
    for (auto &kv : H)
      kv.second.EndEvent(w);
    sumw += w;
    nev++;
  };

  for (int a = 2; a < argc; a++) {
    std::ifstream f(argv[a]);
    std::string line;
    double w = -1.;
    std::vector<Particle> parts;
    while (std::getline(f, line)) {
      if (line.empty())
        continue;
      if (line[0] == '#') {
        std::istringstream ss(line);
        std::vector<std::string> tok;
        std::string t;
        while (ss >> t)
          tok.push_back(t);
        if (tok.size() > 4 && tok[1] == "Event") {
          if (w >= 0.)
            process(parts, w);
          w = std::atof(tok[4].c_str());
          parts.clear();
        }
        continue;
      }
      std::istringstream ss(line);
      int idx, pid, st;
      double E, px, py, pz;
      if (!(ss >> idx >> pid >> st >> E >> px >> py >> pz))
        continue;
      if (E <= 0.)
        continue;
      parts.push_back({PseudoJet(px, py, pz, E),
                       kCharged.count(std::abs(pid)) > 0, st < 0});
    }
    if (w >= 0.)
      process(parts, w);
  }

  std::ofstream out(argv[1]);
  out << "sumw " << sumw << " nev " << nev << "\n";
  for (auto &kv : H) {
    out << kv.first;
    for (double e : kv.second.edges)
      out << " " << e;
    out << "\n" << kv.first << "_h";
    for (double v : kv.second.h)
      out << " " << v;
    out << "\n" << kv.first << "_h2";
    for (double v : kv.second.h2)
      out << " " << v;
    out << "\n";
  }
  return 0;
}
