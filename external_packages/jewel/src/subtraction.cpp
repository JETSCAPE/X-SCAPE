// JEWEL 2.6.0 C++ port -- event-wide constituent subtraction of thermal
// momenta (DOSUBTRACTION). Modernised from subtraction.C of the JEWEL 2.6.0
// distribution; the algorithm, including its quirks, is unchanged:
//   * `lstart` evaluates to 2 for every process (see NOTE below),
//   * particle/thermal pairs are ordered with qsort (as in the original, so
//     that ties are resolved identically),
//   * the species of left-over thermal momenta is drawn with C rand().
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <numbers>
#include <vector>

namespace {

constexpr int kLines = 23000;
using KArray = int (*)[kLines];
using PArray = double (*)[kLines];

struct Mom4 {
  double pt, mdelta, phi, y;
};
struct Particle {
  double pt, mdelta, phi, y;
  int id;
};
struct Pair {
  double dR;
  std::size_t ipart, itherm;
};

constexpr double sqr(double x) { return x * x; }

int comparePairs(const void* a, const void* b) {
  const auto* d1 = static_cast<const Pair*>(a);
  const auto* d2 = static_cast<const Pair*>(b);
  if (d1->dR < d2->dR) return -1;
  if (d1->dR > d2->dR) return 1;
  return 0;
}

void printLine(int i, double a, double b, const double* p0, const double* p1, const double* p2,
               const double* p3, const double* p4) {
  std::cout << i << ": " << a << " - " << b << " -> (" << *p0 << " , " << *p1 << " , " << *p2
            << " ; " << *p3 << ") m = " << *p4 << '\n';
}

std::vector<Particle> readEvent(KArray k, PArray p, int n, int lstart) {
  std::vector<Particle> parts;
  for (int i = lstart; i < n; ++i) {
    if (k[0][i] < 6 || k[0][i] == 17) {
      if (p[3][i] - std::fabs(p[2][i]) > 0.) {
        const double pT = std::sqrt(sqr(p[0][i]) + sqr(p[1][i]));
        const double mass2 = sqr(p[3][i]) - sqr(pT) - sqr(p[2][i]);
        parts.push_back({.pt = pT,
                         .mdelta = std::sqrt(mass2 + sqr(pT)) - pT,
                         .phi = std::atan2(p[1][i], p[0][i]),
                         .y = 0.5 * std::log((p[3][i] + p[2][i]) / (p[3][i] - p[2][i])),
                         .id = k[1][i]});
      } else {
        printLine(i, k[0][i], k[1][i], &p[0][i], &p[1][i], &p[2][i], &p[3][i], &p[4][i]);
      }
    }
  }
  return parts;
}

std::vector<Mom4> readThermals(PArray sc, int nscatcen) {
  std::vector<Mom4> thermals;
  for (int i = 0; i < nscatcen; ++i) {
    if (sc[3][i] - std::fabs(sc[2][i]) > 0.) {
      const double pT = std::sqrt(sqr(sc[0][i]) + sqr(sc[1][i]));
      const double mass2 = sqr(sc[3][i]) - sqr(pT) - sqr(sc[2][i]);
      thermals.push_back({.pt = pT,
                          .mdelta = std::sqrt(mass2 + sqr(pT)) - pT,
                          .phi = std::atan2(sc[1][i], sc[0][i]),
                          .y = 0.5 * std::log((sc[3][i] + sc[2][i]) / (sc[3][i] - sc[2][i]))});
    } else {
      printLine(i, sc[0][i], sc[1][i], &sc[0][i], &sc[1][i], &sc[2][i], &sc[3][i], &sc[4][i]);
    }
  }
  return thermals;
}

// All particle-thermal distances, sorted in dR.
std::vector<Pair> buildPairs(const std::vector<Particle>& parts, const std::vector<Mom4>& thermals) {
  std::vector<Pair> dists;
  dists.reserve(parts.size() * thermals.size());
  for (std::size_t i = 0; i < parts.size(); ++i) {
    for (std::size_t j = 0; j < thermals.size(); ++j) {
      double dphi = std::abs(parts[i].phi - thermals[j].phi);
      if (dphi > std::numbers::pi) dphi = 2. * std::numbers::pi - dphi;
      dists.push_back({std::sqrt(sqr(dphi) + sqr(parts[i].y - thermals[j].y)), i, j});
    }
  }
  if (!dists.empty()) std::qsort(dists.data(), dists.size(), sizeof(Pair), comparePairs);
  return dists;
}

// Go through all particle-thermal pairs and re-distribute pT and mass.
void doSubtraction(const std::vector<Pair>& dists, std::vector<Particle>& parts,
                   std::vector<Mom4>& thermals, double dRmax) {
  for (const Pair& d : dists) {
    if (d.dR > dRmax) break;
    Particle& part = parts[d.ipart];
    Mom4& th = thermals[d.itherm];
    if (part.pt > th.pt) {
      part.pt -= th.pt;
      th.pt = 0.;
    } else {
      th.pt -= part.pt;
      part.pt = 0.;
    }
    if (part.mdelta > th.mdelta) {
      part.mdelta -= th.mdelta;
      th.mdelta = 0.;
    } else {
      th.mdelta -= part.mdelta;
      part.mdelta = 0.;
    }
  }
}

template <typename T>
void storeLine(KArray k, PArray p, int line, int id, const T& m) {
  k[0][line] = 1;
  k[1][line] = id;
  k[2][line] = 0;
  k[3][line] = 0;
  k[4][line] = 0;
  p[0][line] = m.pt * std::cos(m.phi);
  p[1][line] = m.pt * std::sin(m.phi);
  p[2][line] = (m.pt + m.mdelta) * std::sinh(m.y);
  p[3][line] = (m.pt + m.mdelta) * std::cosh(m.y);
  p[4][line] = std::sqrt(sqr(m.pt + m.mdelta) - sqr(m.pt));
}

// Collect resulting 4-momenta to get the subtracted event.
void collectResult(const std::vector<Particle>& parts, const std::vector<Mom4>& thermals, KArray k,
                   PArray p, int* pyjetsN, int lstart) {
  int n = lstart;
  // NOTE(port): as in subtraction.C, the loop starts at index lstart of the
  // particle list (not of the event record).
  for (std::size_t i = static_cast<std::size_t>(lstart); i < parts.size(); ++i) {
    if (parts[i].pt > 0.) {
      ++n;
      storeLine(k, p, n - 1, parts[i].id, parts[i]);
    }
  }
  for (const Mom4& th : thermals) {
    if (th.pt > 0.) {
      const double r = double(std::rand()) / double(RAND_MAX);
      int thermid = 111;
      if (r < 1. / 3.)
        thermid = 211;
      else if (r < 2. / 3.)
        thermid = -211;
      ++n;
      if (n == kLines) {
        *pyjetsN = 0;
        return;
      }
      storeLine(k, p, n - 1, thermid, th);
    }
  }
  *pyjetsN = n;
}

}  // namespace

extern "C" void subtract_thmom_(int pyjetsK[5][kLines], double pyjetsP[5][kLines], int* pyjetsN,
                                double scatcens[5][kLines], int* nscatcen, double* Rmax,
                                char* proc) {
  // NOTE(port): the original conditions are `if (strcmp(proc,"PPYJ")) lstart=1;
  // if (strcmp(proc,"PPZJ") || strcmp(proc,"PPWJ")) lstart=2;`. strcmp returns
  // non-zero for *different* strings, and proc can never equal both "PPZJ" and
  // "PPWJ", so lstart is 2 for every process. Kept for identical results.
  int lstart = 0;
  if (std::strcmp(proc, "PPYJ")) lstart = 1;
  if (std::strcmp(proc, "PPZJ") || std::strcmp(proc, "PPWJ")) lstart = 2;
  std::vector<Particle> parts = readEvent(pyjetsK, pyjetsP, *pyjetsN, lstart);
  std::vector<Mom4> thermals = readThermals(scatcens, *nscatcen);
  const std::vector<Pair> dists = buildPairs(parts, thermals);
  doSubtraction(dists, parts, thermals, *Rmax);
  collectResult(parts, thermals, pyjetsK, pyjetsP, pyjetsN, lstart);
}
