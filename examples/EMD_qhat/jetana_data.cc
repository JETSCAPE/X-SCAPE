// Data-matched jet observables from kept JETSCAPE final-state event files
// (plain or gzipped), binned exactly as the HEPData records they are compared
// with.  Output has the jetana format (sum(w*n), sum((w*n)^2), sum(w)).
//
//   jetana_data out.spec events1.dat.gz [events2.dat.gz ...]
//
// ATLAS dijet x_J (PLB 846 (2023) 138954, HEPData ins2075431)
//   R = 0.4 full jets (holes subtracted); the two highest-pT jets with
//   |y| < 2.8 must both have |y| < 2.1 and dphi > 7pi/8; pT1 in the HEPData
//   leading-jet bins; atl_xj<k> in the HEPData log x_J bins (0.316-1);
//   atl_xjlead counts leading jets per pT1 bin
// ATLAS D(z) (PRC 98 (2018) 024908, HEPData ins1673177)
//   R = 0.4 full jets, |y| < 2.1, pT 126-158-200-251-316-398 GeV; charged
//   particles pT > 1 GeV, dR < 0.4 to the jet axis (holes negative),
//   z = pT cos(dR) / pT_jet in the HEPData z bins of each pT bin; atl_ffnjet
//   counts jets; atl_ffzbkg<k>: underlying event from two cones at
//   phi_jet +- pi/2 (each weighted 1/2), to be subtracted per jet
// CMS jet shape (JHEP 05 (2018) 006, HEPData ins1658057)
//   R = 0.4 full jets pT > 120 GeV, |eta| < 1.6; charged particles
//   pT > 0.7 GeV, |eta| < 2.4, within dr < 1 (holes negative), pT/pT_jet in
//   the HEPData annuli; normalise to unit integral over 0 < dr < 1;
//   cms3_*: jets 0.3 < |eta| < 1.6 with an eta-reflected background cone
//   (cms3_rhobkg; a cross-check, not the CMS method);
//   cms_rhoall, cms_rhosb: the CMS method (JHEP 05 (2018) 006, in pp and
//   PbPb): particles at any eta (acceptance-corrected data), and the
//   background from the sideband 1.5 < |deta| < 2.5 as a function of dphi,
//   spread over each annulus by its deta length at that dphi; subtract
//   cms_rhosb from cms_rhoall
// ALICE charged-jet mass and angularities (arXiv:2411.03106, ins2845788)
//   anti-kT R = 0.2 on charged particles pT > 0.15 GeV, |eta_jet| < 0.7,
//   charged-jet pT 40-60-80-100-150 GeV (holes subtracted);
//   ali_m<k>, ali_l<a>_<k> (a = 1, 15, 2, 3 for alpha = 1, 1.5, 2, 3):
//   ungroomed, holes subtracted (4-momentum and negative angularity terms);
//   ali_mg<k>, ali_lg<a>_<k>: Soft Drop z_cut = 0.2, beta = 0 on the C/A
//   reclustered particles and holes of the jet, prongs and groomed jet as
//   signed sums (holes negative); ali_mgnh, ali_lgnh: the same grooming of
//   the particles alone; jets failing the grooming go to the [-0.1, 0)
//   "untagged" bin
//   lambda_a = sum_i pT_i / pT_jet (dR_i / R)^a, dR to the (groomed) axis.
// Events that contain quarks or gluons in the final state are skipped and
// counted ("skipped" in the output): ColoredHadronization does not check
// pythia.next(), and after its first junction failure in a run it writes the
// unhadronized partons of every later event.  No other sample has any.
// With JETANA_PIONMASS=1 in the environment, every charged particle gets the
// pion mass (E recomputed from its momentum), as in ALICE's track-based jets.
#include "fjcore.hh"

#include <cmath>
#include <cstdio>
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
const std::set<int> kCharged = {211, 321, 2212, 3222, 3112, 3312, 3334, 11, 13};

struct Hist {
  std::vector<double> edges, h, h2, ev;
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

std::vector<double> Log(double lo, double hi, int n) {
  std::vector<double> e;
  for (int i = 0; i <= n; i++)
    e.push_back(lo * std::pow(hi / lo, double(i) / n));
  return e;
}

int BinOf(const std::vector<double> &b, double x) {
  for (size_t i = 0; i + 1 < b.size(); i++)
    if (x >= b[i] && x < b[i + 1])
      return int(i);
  return -1;
}

struct Particle {
  PseudoJet p;
  bool charged, hole;
};

// HEPData binning
const std::vector<double> kXjPt = {100, 112, 126, 141, 158, 178, 200, 224, 251, 282, 316, 398};
const std::vector<double> kFfPt = {126, 158.49, 199.53, 251.19, 316.22, 398.1};
const std::vector<double> kCmsR = {0, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .6, .7, .8, 1.};
const std::vector<double> kAliPt = {40, 60, 80, 100, 150};
const std::vector<std::string> kAlpha = {"1", "15", "2", "3"};
const double kAlphaVal[] = {1., 1.5, 2., 3.};
// [pT bin][alpha] ungroomed and groomed angularity edges, ALICE Tables 1-32
const std::vector<std::vector<std::vector<double>>> kAng = {
    {{0, .05, .1, .15, .2, .25, .3, .35, .4, .55}, {0, .05, .08, .1, .14, .18, .2, .24, .3, .5},
     {0, .02, .04, .07, .1, .13, .16, .2, .3}, {0, .01, .02, .03, .04, .05, .06, .08, .1, .15, .25}},
    {{0, .1, .15, .2, .25, .3, .5}, {0, .05, .08, .12, .17, .25, .4},
     {0, .02, .04, .07, .1, .13, .2, .3}, {0, .01, .02, .04, .06, .08, .12, .25}},
    {{0, .1, .15, .2, .25, .3, .5}, {0, .04, .08, .12, .17, .35},
     {0, .02, .04, .07, .1, .25}, {0, .01, .02, .04, .07, .1}},
    {{0, .07, .1, .12, .14, .17, .2, .24, .3, .5}, {0, .03, .06, .09, .12, .16, .2, .26, .4},
     {0, .02, .03, .04, .05, .07, .1, .15, .3}, {0, .01, .02, .03, .04, .055, .075, .1, .2}}};
const std::vector<std::vector<std::vector<double>>> kAngG = {
    {{-.1, 0, .05, .1, .15, .2, .25, .3, .4, .55}, {-.1, 0, .02, .05, .08, .11, .15, .2, .3, .5},
     {-.1, 0, .02, .04, .07, .1, .13, .16, .25}, {-.1, 0, .005, .01, .02, .03, .05, .06, .08, .1, .15}},
    {{-.1, 0, .04, .08, .15, .2, .25, .3, .4, .55}, {-.1, 0, .015, .03, .06, .1, .15, .25, .4},
     {-.1, 0, .005, .01, .02, .03, .05, .07, .1, .2, .35}, {-.1, 0, .005, .01, .02, .03, .05, .07, .1, .15, .25}},
    {{-.1, 0, .05, .1, .2, .3, .45}, {-.1, 0, .02, .05, .1, .35},
     {-.1, 0, .005, .01, .03, .08, .25}, {-.1, 0, .001, .005, .02, .05, .1}},
    {{-.1, 0, .03, .05, .07, .09, .11, .15, .2, .3, .4, .55}, {-.1, 0, .007, .01, .015, .02, .04, .07, .1, .25, .45},
     {-.1, 0, .002, .004, .007, .01, .02, .04, .1, .35}, {-.1, 0, .001, .002, .005, .01, .02, .05, .1, .25}}};
// ALICE Tables 39-46 (Pb-Pb; the pp tables 33-38 share the first three)
const std::vector<std::vector<double>> kMass = {
    {0, 1.2, 1.6, 2, 2.4, 2.8, 3.2, 3.6, 4, 4.4, 5, 5.6, 6.2, 7},
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
    {0, 2, 4, 6, 8, 12},
    {0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 16}};
const std::vector<std::vector<double>> kMassG = {
    {-.1, 0, .3, .6, .9, 1.2, 1.6, 2, 2.4, 2.8, 3.2, 3.6, 4, 4.6, 5.2, 6, 7},
    {-.1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
    {-.1, 0, 1.8, 3.6, 5.8, 8, 12},
    {-.1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 16}};

// Signed 4-momentum of a C/A cluster: holes (user_index -1) count negatively.
PseudoJet Signed(const PseudoJet &j) {
  PseudoJet s(0., 0., 0., 0.);
  for (auto &c : j.constituents())
    s += c.user_index() < 0 ? PseudoJet(-c.px(), -c.py(), -c.pz(), -c.E()) : c;
  return s;
}

// Soft Drop (beta = 0) with holes: particles and holes inside the jet are
// reclustered together with C/A, and every prong and the groomed jet are the
// signed sums (holes negative).  gcons gets the groomed constituents with the
// sign in user_index.  False if untagged.
bool SoftDropSigned(const std::vector<PseudoJet> &cons, double zcut, PseudoJet &groomed,
                    std::vector<PseudoJet> &gcons) {
  if (cons.size() < 2)
    return false;
  fjcore::JetDefinition ca(fjcore::cambridge_algorithm, 10.);
  fjcore::ClusterSequence cs(cons, ca);
  std::vector<PseudoJet> js = cs.exclusive_jets(1);
  PseudoJet j = js[0], a, b;
  while (j.has_parents(a, b)) {
    double pa = Signed(a).pt(), pb = Signed(b).pt();
    if (Signed(a).E() > 0. && Signed(b).E() > 0. && std::min(pa, pb) / (pa + pb) > zcut) {
      groomed = Signed(j);
      gcons.clear();
      for (auto &c : j.constituents()) {
        gcons.push_back(PseudoJet(c.px(), c.py(), c.pz(), c.E()));
        gcons.back().set_user_index(c.user_index());
      }
      return true;
    }
    j = Signed(a).pt() > Signed(b).pt() ? a : b;
  }
  return false;
}

// Soft Drop (beta = 0) on C/A-reclustered constituents; false if untagged.
bool SoftDrop(const std::vector<PseudoJet> &cons, double zcut, PseudoJet &groomed,
              std::vector<PseudoJet> &gcons) {
  if (cons.size() < 2)
    return false;
  fjcore::JetDefinition ca(fjcore::cambridge_algorithm, 10.);
  fjcore::ClusterSequence cs(cons, ca);
  std::vector<PseudoJet> js = cs.exclusive_jets(1);
  PseudoJet j = js[0], a, b;
  while (j.has_parents(a, b)) {
    double pa = a.pt(), pb = b.pt();
    if (std::min(pa, pb) / (pa + pb) > zcut) {
      groomed = PseudoJet(j.px(), j.py(), j.pz(), j.E());
      gcons.clear();
      for (auto &c : j.constituents())
        gcons.push_back(PseudoJet(c.px(), c.py(), c.pz(), c.E()));
      return true;
    }
    j = pa > pb ? a : b;
  }
  return false;
}
}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: jetana_data out.spec files...\n";
    return 1;
  }
  std::map<std::string, Hist> H;
  H["atl_xjlead"] = Hist(kXjPt);
  H["atl_xjpair"] = Hist(kXjPt);
  for (size_t k = 0; k + 1 < kXjPt.size(); k++)
    H["atl_xj" + std::to_string(k)] = Hist(Log(std::pow(10., -0.5), 1., 10));
  H["atl_ffnjet"] = Hist(kFfPt);
  // z bins differ between the HEPData pT bins
  std::vector<double> zFine = Log(0.01, 1., 20), zCoarse = {
      0.0068, 0.01, 0.0158, 0.0251, 0.0398, 0.0631, 0.1, 0.1585, 0.2512, 0.3981, 0.631, 1.};
  std::vector<std::vector<double>> zEdges = {{0.00794328}, {0.00794328}, {0.00630957, 0.00794328},
                                             {0.0046}, {0.0032, 0.0046}};
  for (size_t k = 0; k < zEdges.size(); k++) {
    for (double e : (k < 3 ? zFine : zCoarse))
      zEdges[k].push_back(e);
    H["atl_ffz" + std::to_string(k)] = Hist(zEdges[k]);
    // underlying event: two cones at phi_jet +- pi/2, same rapidity (x 1/2)
    H["atl_ffzbkg" + std::to_string(k)] = Hist(zEdges[k]);
  }
  H["cms_njet"] = Hist({120, 10000});
  H["cms_rho"] = Hist(kCmsR);
  H["cms_rhopos"] = Hist(kCmsR);  // diagnostics: particles only, holes only
  H["cms_rhohole"] = Hist(kCmsR);
  // CMS eta-reflection: jets 0.3 < |eta| < 1.6, background from the cone at
  // (-eta_jet, phi_jet), as subtracted in the measurement
  H["cms_rhoall"] = Hist(kCmsR);
  H["cms_rhosb"] = Hist(kCmsR);
  H["cms3_njet"] = Hist({120, 10000});
  H["cms3_rho"] = Hist(kCmsR);
  H["cms3_rhobkg"] = Hist(kCmsR);
  H["ali_njet"] = Hist(kAliPt);
  for (size_t k = 0; k + 1 < kAliPt.size(); k++) {
    std::string t = std::to_string(k);
    H["ali_m" + t] = Hist(kMass[k]);
    H["ali_mg" + t] = Hist(kMassG[k]);
    H["ali_mgnh" + t] = Hist(kMassG[k]);  // grooming without holes (v1)
    for (size_t a = 0; a < kAlpha.size(); a++) {
      H["ali_l" + kAlpha[a] + "_" + t] = Hist(kAng[k][a]);
      H["ali_lg" + kAlpha[a] + "_" + t] = Hist(kAngG[k][a]);
      H["ali_lgnh" + kAlpha[a] + "_" + t] = Hist(kAngG[k][a]);
    }
  }
  double sumw = 0.;
  long nev = 0, nskip = 0;
  const char *pm = std::getenv("JETANA_PIONMASS");
  const bool pionMass = pm && std::string(pm) == "1";

  auto process = [&](const std::vector<Particle> &parts, double w) {
    std::vector<PseudoJet> in, holes, ch, chholes;
    for (auto &q : parts) {
      if (q.hole) {
        holes.push_back(q.p);
        if (q.charged && q.p.pt() > 0.15)
          chholes.push_back(q.p);
      } else {
        in.push_back(q.p);
        if (q.charged && q.p.pt() > 0.15)
          ch.push_back(q.p);
      }
    }
    // ---- R = 0.4 full jets, holes subtracted
    fjcore::JetDefinition def4(fjcore::antikt_algorithm, 0.4);
    fjcore::ClusterSequence cs(in, def4);
    std::vector<PseudoJet> ax, sub;
    for (auto &j : fjcore::sorted_by_pt(cs.inclusive_jets(10.))) {
      PseudoJet s = j;
      for (auto &h : holes)
        if (j.delta_R(h) < 0.4)
          s -= h;
      ax.push_back(j);
      sub.push_back(s);
    }
    // ATLAS x_J
    {
      int i1 = -1, i2 = -1;
      for (size_t i = 0; i < sub.size(); i++) {
        if (std::abs(ax[i].rap()) >= 2.8)
          continue;
        if (i1 < 0 || sub[i].pt() > sub[i1].pt()) {
          i2 = i1;
          i1 = i;
        } else if (i2 < 0 || sub[i].pt() > sub[i2].pt()) {
          i2 = i;
        }
      }
      if (i1 >= 0 && std::abs(ax[i1].rap()) < 2.1) {
        double pt1 = sub[i1].pt();
        int k = BinOf(kXjPt, pt1);
        if (k >= 0) {
          H["atl_xjlead"].Fill(pt1);
          if (i2 >= 0 && std::abs(ax[i2].rap()) < 2.1 &&
              std::abs(ax[i1].delta_phi_to(ax[i2])) > 7. * M_PI / 8.) {
            double xj = sub[i2].pt() / pt1;
            if (xj >= std::pow(10., -0.5)) {
              H["atl_xjpair"].Fill(pt1);
              H["atl_xj" + std::to_string(k)].Fill(xj);
            }
          }
        }
      }
    }
    for (size_t i = 0; i < sub.size(); i++) {
      double pt = sub[i].pt();
      // ATLAS D(z)
      int k = std::abs(ax[i].rap()) < 2.1 ? BinOf(kFfPt, pt) : -1;
      if (k >= 0) {
        H["atl_ffnjet"].Fill(pt);
        std::string t = "atl_ffz" + std::to_string(k);
        for (auto &q : parts) {
          if (!q.charged || q.p.pt() <= 1.0)
            continue;
          double dr = q.p.delta_R(ax[i]);
          if (dr < 0.4)
            H[t].Fill(q.p.pt() * std::cos(dr) / pt, q.hole ? -1. : 1.);
          for (double sgn : {-1., 1.}) {
            PseudoJet perp = fjcore::PtYPhiM(1., ax[i].rap(), ax[i].phi() + sgn * M_PI / 2.);
            double db = q.p.delta_R(perp);
            if (db < 0.4)
              H["atl_ffzbkg" + std::to_string(k)].Fill(q.p.pt() * std::cos(db) / pt,
                                                        0.5 * (q.hole ? -1. : 1.));
          }
        }
      }
      // CMS rho(dr), eta-reflection selection
      double aeta = std::abs(ax[i].eta());
      if (pt > 120. && aeta > 0.3 && aeta < 1.6) {
        H["cms3_njet"].Fill(pt);
        PseudoJet refl = fjcore::PtYPhiM(1., -ax[i].eta(), ax[i].phi());
        for (auto &q : parts) {
          if (!q.charged || q.p.pt() <= 0.7 || std::abs(q.p.eta()) >= 2.4)
            continue;
          double v = (q.hole ? -1. : 1.) * q.p.pt() / pt;
          double dr = q.p.delta_R(ax[i]), db = q.p.delta_R(refl);
          if (dr < 1.0)
            H["cms3_rho"].Fill(dr, v);
          if (db < 1.0)
            H["cms3_rhobkg"].Fill(db, v);
        }
      }
      // CMS rho(dr)
      if (pt > 120. && std::abs(ax[i].eta()) < 1.6) {
        H["cms_njet"].Fill(pt);
        Hist &all = H["cms_rhoall"], &sb = H["cms_rhosb"];
        for (auto &q : parts) {
          if (!q.charged || q.p.pt() <= 0.7)
            continue;
          double v = (q.hole ? -1. : 1.) * q.p.pt() / pt;
          double dr = q.p.delta_R(ax[i]);
          if (dr < 1.0)
            all.Fill(dr, v);
          double deta = std::abs(q.p.eta() - ax[i].eta());
          double dphi = std::abs(q.p.delta_phi_to(ax[i]));
          if (deta > 1.5 && deta < 2.5 && dphi < 1.0) {
            // background per unit deta (sideband width 2) times the deta
            // length of each annulus at this dphi
            for (size_t b = 0; b + 1 < kCmsR.size(); b++) {
              double r1 = kCmsR[b], r2 = kCmsR[b + 1];
              if (dphi >= r2)
                continue;
              double len = 2. * (std::sqrt(r2 * r2 - dphi * dphi) -
                                 (dphi < r1 ? std::sqrt(r1 * r1 - dphi * dphi) : 0.));
              sb.Fill(0.5 * (r1 + r2), v * len / 2.);
            }
          }
        }
        for (auto &q : parts) {
          if (!q.charged || q.p.pt() <= 0.7 || std::abs(q.p.eta()) >= 2.4)
            continue;
          double dr = q.p.delta_R(ax[i]);
          if (dr < 1.0) {
            H["cms_rho"].Fill(dr, (q.hole ? -1. : 1.) * q.p.pt() / pt);
            H[q.hole ? "cms_rhohole" : "cms_rhopos"].Fill(dr, q.p.pt() / pt);
          }
        }
      }
    }
    // ---- ALICE R = 0.2 charged jets
    const double R2 = 0.2;
    fjcore::JetDefinition def2(fjcore::antikt_algorithm, R2);
    fjcore::ClusterSequence cs2(ch, def2);
    for (auto &j : cs2.inclusive_jets(10.)) {
      if (std::abs(j.eta()) >= 0.7)
        continue;
      PseudoJet s = j;
      std::vector<PseudoJet> jh;
      for (auto &h : chholes)
        if (j.delta_R(h) < R2) {
          s -= h;
          jh.push_back(h);
        }
      double pt = s.pt();
      int k = BinOf(kAliPt, pt);
      if (k < 0)
        continue;
      std::string t = std::to_string(k);
      H["ali_njet"].Fill(pt);
      double m2 = s.m2();
      H["ali_m" + t].Fill(m2 > 0. ? std::sqrt(m2) : 0.);
      std::vector<PseudoJet> cons;
      for (auto &c : j.constituents())
        cons.push_back(PseudoJet(c.px(), c.py(), c.pz(), c.E()));
      for (size_t a = 0; a < kAlpha.size(); a++) {
        double l = 0.;
        for (auto &c : cons)
          l += c.pt() * std::pow(c.delta_R(j) / R2, kAlphaVal[a]);
        for (auto &h : jh)
          l -= h.pt() * std::pow(h.delta_R(j) / R2, kAlphaVal[a]);
        H["ali_l" + kAlpha[a] + "_" + t].Fill(l / pt);
      }
      PseudoJet g;
      std::vector<PseudoJet> gc;
      if (SoftDrop(cons, 0.2, g, gc)) {
        double mg2 = g.m2();
        H["ali_mgnh" + t].Fill(mg2 > 0. ? std::sqrt(mg2) : 0.);
        for (size_t a = 0; a < kAlpha.size(); a++) {
          double l = 0.;
          for (auto &c : gc)
            l += c.pt() * std::pow(c.delta_R(g) / R2, kAlphaVal[a]);
          H["ali_lgnh" + kAlpha[a] + "_" + t].Fill(l / g.pt());
        }
      } else {
        H["ali_mgnh" + t].Fill(-0.05);
        for (size_t a = 0; a < kAlpha.size(); a++)
          H["ali_lgnh" + kAlpha[a] + "_" + t].Fill(-0.05);
      }
      // with holes (default)
      std::vector<PseudoJet> scons;
      for (auto &c : cons) {
        scons.push_back(c);
        scons.back().set_user_index(1);
      }
      for (auto &h : jh) {
        scons.push_back(PseudoJet(h.px(), h.py(), h.pz(), h.E()));
        scons.back().set_user_index(-1);
      }
      if (SoftDropSigned(scons, 0.2, g, gc)) {
        double mg2 = g.m2();
        H["ali_mg" + t].Fill(mg2 > 0. ? std::sqrt(mg2) : 0.);
        for (size_t a = 0; a < kAlpha.size(); a++) {
          double l = 0.;
          for (auto &c : gc)
            l += (c.user_index() < 0 ? -1. : 1.) * c.pt() *
                 std::pow(c.delta_R(g) / R2, kAlphaVal[a]);
          H["ali_lg" + kAlpha[a] + "_" + t].Fill(l / g.pt());
        }
      } else {
        H["ali_mg" + t].Fill(-0.05);
        for (size_t a = 0; a < kAlpha.size(); a++)
          H["ali_lg" + kAlpha[a] + "_" + t].Fill(-0.05);
      }
    }
    for (auto &kv : H)
      kv.second.EndEvent(w);
    sumw += w;
    nev++;
  };

  for (int a = 2; a < argc; a++) {
    std::string fn = argv[a];
    bool gz = fn.size() > 3 && fn.substr(fn.size() - 3) == ".gz";
    FILE *fp = gz ? popen(("gzip -dc '" + fn + "'").c_str(), "r") : fopen(fn.c_str(), "r");
    if (!fp) {
      std::cerr << "cannot open " << fn << "\n";
      return 1;
    }
    char buf[4096];
    double w = -1.;
    bool partons = false;
    std::vector<Particle> parts;
    while (fgets(buf, sizeof buf, fp)) {
      if (buf[0] == '#') {
        std::istringstream ss(buf);
        std::vector<std::string> tok;
        std::string t;
        while (ss >> t)
          tok.push_back(t);
        if (tok.size() > 4 && tok[1] == "Event") {
          if (w >= 0.)
            partons ? (void)nskip++ : process(parts, w);
          w = std::atof(tok[4].c_str());
          parts.clear();
          partons = false;
        }
        continue;
      }
      int idx, pid, st;
      double E, px, py, pz;
      if (std::sscanf(buf, "%d %d %d %lf %lf %lf %lf", &idx, &pid, &st, &E, &px, &py, &pz) != 7)
        continue;
      if (pid == 21 || (pid != 0 && std::abs(pid) <= 6))
        partons = true;
      if (E <= 0.)
        continue;
      bool charged = kCharged.count(std::abs(pid)) > 0;
      if (pionMass && charged)
        E = std::sqrt(px * px + py * py + pz * pz + 0.13957 * 0.13957);
      parts.push_back({PseudoJet(px, py, pz, E), charged, st < 0});
    }
    if (w >= 0.)
      partons ? (void)nskip++ : process(parts, w);
    gz ? pclose(fp) : fclose(fp);
  }

  std::ofstream out(argv[1]);
  out << "sumw " << sumw << " nev " << nev << " skipped " << nskip << "\n";
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
