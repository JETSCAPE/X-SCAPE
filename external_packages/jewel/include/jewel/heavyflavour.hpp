// JEWEL C++ -- heavy-flavour extension (vacuum part).
//
// Switched on with HEAVYFLAVOUR T in the parameter file (default F: original JEWEL, bit-identical).
// With it on, charm and bottom quarks keep their mass through the final-state shower:
//
//  * Conventions.  P(l,5) of a heavy quark is its invariant mass M, M^2 = m^2 + t with t >= 0 the
//    off-shellness; a heavy quark is on shell when M = m.  Everything in the JEWEL core that is a
//    virtuality (Sudakov evolution variable, formation time 2E/t, angular-ordering estimates)
//    uses t, everything that is kinematics uses M.  getmass() keeps returning the virtuality
//    sqrt(t) (0 = no further splitting); callers convert with invmass().  For light partons
//    m = 0 and all helpers return the original JEWEL expressions.
//
//  * Q -> Q g uses the quasi-collinear massive splitting function (Catani, Dittmaier, Trocsanyi,
//    hep-ph/0011222)  P = C_F [ (1+z^2)/(1-z) - 2 m^2 / t ],  z = energy fraction of the quark;
//    it is <= the massless kernel, so it enters as an extra veto in getmass, and produces the
//    dead cone (radiation suppressed for theta < m/E).
//
//  * g -> Q Qbar (Q = c, b) is added to the gluon Sudakov with  P = T_R [ z^2 + (1-z)^2 + 2 m^2/t ]
//    above threshold.  Light g -> q qbar keeps JEWEL's original weight.
//
//  * z limits.  JEWEL's  eps = 1/2 - 1/2 sqrt(1 - Q0^2/t) sqrt(1 - t/E^2)  is exactly the
//    two-body phase space of a decay into two daughters of mass Q0/2.  The heavy version uses the
//    same formula with daughter masses sqrt(m^2 + Q0^2/4) (heavy quark) and Q0/2 (gluon).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "jewel/commons.hpp"
#include "jewel/pythia6.hpp"

namespace jewel {

namespace cb::hf {
inline bool on = false;           // HEAVYFLAVOUR
inline bool stable = false;       // HFSTABLE: weakly decaying c and b hadrons are not decayed
inline double hqds = 0.0;         // HQDS: 2 pi T D_s of the heavy-quark Langevin transport (0 = off)
inline double mc = 0.0, mb = 0.0; // charm and bottom masses = PYTHIA PMAS(4,1), PMAS(5,1)
// flavour of a g -> q qbar splitting decided in getmass (1, 2, 3: d, u, s; 4: c cbar; 5: b bbar;
// 0: not decided by getmass, e.g. HEAVYFLAVOUR F or g -> g g, then g -> q qbar gives d dbar as in JEWEL)
inline F1<int, kLines> qqbarfl;
// flavour belonging to the most recent (zdec, qqbardec) decided by a medium routine (getnewmass, medium
// getmass calls); stored into qqbarfl when the line is created (dokinematics, makeinsplit)
inline int pendingfl = 0;
// restricted hard process (PROCESS PPQQ / PPCC / PPBB): flavour-creation 2 -> 2 only, e.g. to compare with
// Craft, Lee, Mecaj, Moult, arXiv:2210.09311 (Pythia gg, qqbar -> Q Qbar).  Internally the run is a PPJJ run.
//   "QQ": q qbar -> q' qbar' and g g -> q qbar (ISUB 12, 53), light flavours d, u, s only
//   "CC", "BB": g g, q qbar -> c cbar / b bbar (ISUB 81, 82, MSTP(7) = 4 / 5)
inline std::string hardproc;
}  // namespace cb::hf

// mass of a parton in the JEWEL shower: m_c, m_b for heavy quarks if HEAVYFLAVOUR, otherwise 0
inline double hfmass(int kf) noexcept {
  if (!cb::hf::on) return 0.0;
  const int a = std::abs(kf);
  if (a == 4) return cb::hf::mc;
  if (a == 5) return cb::hf::mb;
  return 0.0;
}
inline bool isheavy(int kf) noexcept { return hfmass(kf) > 0.0; }

// invariant mass from virtuality q = sqrt(t)
// Signed conventions (as JEWEL's P(l,5) < 0 for spacelike partons): the off-shellness is
// t = sign(P5) P5^2 - m^2 and the signed virtuality v = sign(t) sqrt(|t|); for m = 0 both are the identity.
// invariant mass (signed P5) from the signed virtuality v
inline double invmass(double q, double m) noexcept {
  if (m == 0.0) return q;
  const double M2 = (q < 0.0 ? -q * q : q * q) + m * m;
  return M2 < 0.0 ? -std::sqrt(-M2) : std::sqrt(M2);
}
// virtuality sqrt(t) from invariant mass M (0 on shell)
inline double virtuality(double M, double m) noexcept {
  if (m == 0.0) return M;
  const double t = (M < 0.0 ? -M * M : M * M) - m * m;
  return t < 0.0 ? -std::sqrt(-t) : std::sqrt(t);
}
// off-shellness of line l (JEWEL: P(l,5) > 0)
inline bool offshell(int l) noexcept {
  const double m = hfmass(py::pyjets::K(l, 2));
  return m == 0.0 ? py::pyjets::P(l, 5) > 0.0 : py::pyjets::P(l, 5) > m * (1.0 + 1e-9);
}
// virtuality squared t of line l (JEWEL: P(l,5)**2)
inline double virt2(int l) noexcept {
  const double m = hfmass(py::pyjets::K(l, 2));
  const double p5 = py::pyjets::P(l, 5);
  return m == 0.0 ? p5 * p5 : std::max(p5 * p5 - m * m, 0.0);
}

// Allowed range [zlo, zhi] of the energy fraction z of daughter 1 for a parton of energy e and
// virtuality squared t (invariant mass^2 = t + m^2) splitting into daughters of mass m1, m2.
// Returns false if the splitting is kinematically closed.  For m = 0, m1 = m2 = Q0/2 this is
// JEWEL's [eps, 1-eps].
inline bool hfzlimits(double t, double m, double e, double m1, double m2, double& zlo, double& zhi) {
  const double M2 = t + m * m;
  if (M2 >= e * e) return false;
  const double a = m1 * m1 / M2, b = m2 * m2 / M2;
  const double lam = (1.0 - (m1 + m2) * (m1 + m2) / M2) * (1.0 - (m1 - m2) * (m1 - m2) / M2);
  if (lam < 0.0) return false;
  const double beta = std::sqrt(1.0 - M2 / (e * e));
  const double c = 0.5 * (1.0 + a - b), h = 0.5 * beta * std::sqrt(lam);
  zlo = c - h;
  zhi = c + h;
  return zlo < zhi;
}
// diagnostic switches (environment JEWEL_HFTEST, comma separated; testing only):
//   nothin  - no initial-state dead-cone thinning (hfisweight = 1)
//   nodress - Q -> Q g daughters undressed (quark m, gluon 0) instead of sqrt(m^2 + Q0^2/4), Q0/2
//   norad   - no medium-induced radiation of heavy quarks (elastic energy loss only)
//   stats   - print counters of heavy-quark momentum transfers and b-quark scatterings at exit (stderr)
inline bool hftest(const char* flag) {
  const char* e = std::getenv("JEWEL_HFTEST");
  return e && std::strstr(e, flag);
}
// z limits for Q -> Q(z) g(1-z) with JEWEL's resolution cut Q0
inline bool hfzlimitsqq(double t, double m, double e, double& zlo, double& zhi) {
  static const bool nodress = hftest("nodress");
  const double q02 = cb::param::q0 * cb::param::q0 / 4.0;
  if (nodress) return hfzlimits(t, m, e, m, 0.0, zlo, zhi);
  return hfzlimits(t, m, e, std::sqrt(m * m + q02), std::sqrt(q02), zlo, zhi);
}
// z limits for g -> Q(z) Qbar(1-z)
inline bool hfzlimitsgq(double t, double m, double e, double& zlo, double& zhi) {
  const double mq = std::sqrt(m * m + cb::param::q0 * cb::param::q0 / 4.0);
  return hfzlimits(t, 0.0, e, mq, mq, zlo, zhi);
}

// Is z outside the allowed range for a parton with virtuality q, quark mass m (0 for gluons and
// light quarks), energy e?  fl = g -> q qbar flavour of a gluon (0 light / none).  For m = 0 and
// fl = 0 this is JEWEL's (z < eps) || (z > 1 - eps) with the eps passed in.
inline bool hfzout(double q, double m, double e, double z, double eps, int fl = 0) {
  if (m == 0.0 && hfmass(fl) == 0.0) return (z < eps) || (z > (1.0 - eps));
  double zlo, zhi;
  const bool open = (m > 0.0) ? hfzlimitsqq(q * q, m, e, zlo, zhi) : hfzlimitsgq(q * q, hfmass(fl), e, zlo, zhi);
  return !open || (z < zlo) || (z > zhi);
}

// debugging aid (JEWEL_HFDEBUG): report heavy-flavour lines whose stored (P5, zd) cannot be
// realised with on-shell daughters at their actual energy
void hfcheck(const char* stage);
// Langevin transport of final-state c and b quarks (HQDS > 0), hfdiffusion.cpp
void hfdiffuse();

// transverse momentum squared of a splitting with virtuality t, z = energy fraction of daughter 1:
//   Q -> Q g:      kt^2 = z(1-z) t - (1-z)^2 m^2
//   g -> Q Qbar:   kt^2 = z(1-z) t - m^2
inline double hfkt2qq(double t, double m, double z) { return z * (1.0 - z) * t - (1.0 - z) * (1.0 - z) * m * m; }
inline double hfkt2gq(double t, double m, double z) { return z * (1.0 - z) * t - m * m; }
// opening-angle estimate theta = kt / (z (1-z) E) (JEWEL: sqrt(t / (z(1-z))) / E for m = 0)
inline double hftheta(double t, double m, double z, double e) {
  if (m == 0.0) return std::sqrt(t) / (std::sqrt(z * (1.0 - z)) * e);
  return std::sqrt(std::max(hfkt2qq(t, m, z), 0.0)) / (z * (1.0 - z) * e);
}

// e+e- -> gamma*/Z -> Q Qbar g matrix-element correction for the first emission of a heavy quark (JEWEL's
// EEJJ correction is massless and kept for light quarks).  Q Qbar g energy fractions x1 (emitting quark),
// x2 (other quark), r = m/sqrt(s).  ME: PYTHIA's massive three-jet matrix element PYMAEL, class
// V/A -> q qbar with the Z-pole vector fraction alpha = v^2/(v^2 + a^2), normalised like JEWEL's massless
// (x1^2 + x2^2)/((1-x1)(1-x2)).  Shower: the massive kernel of both legs in the same variables,
//   PS_1 = [(1+z^2)/x3 - 2 r^2/((1-x2)(2-x2))] / (1-x2),  z = x1/(2-x2).
// ME/(PS_1 + PS_2) <= 1 over the whole phase space for c and b at the Z (checked numerically).
inline double hfmeweight(double x1, double x2, double r, int kf) {
  const double sw2 = py::pydat1::PARU(102);
  const double ef = (std::abs(kf) % 2 == 0) ? 2.0 / 3.0 : -1.0 / 3.0;
  const double af = (ef > 0.0) ? 1.0 : -1.0, vf = af - 4.0 * ef * sw2;
  const double alpha = vf * vf / (vf * vf + af * af);
  const int ni = 13;  // class 2 (V/A -> q qbar), vector/axial mixture alpha
  const double me = f77::pymael_(&ni, &x1, &x2, &r, &r, &alpha);
  auto ps = [r](double a, double b) {
    const double z = a / (2.0 - b), x3 = 2.0 - a - b;
    return std::max((1.0 + z * z) / x3 - 2.0 * r * r / ((1.0 - b) * (2.0 - b)), 0.0) / (1.0 - b);
  };
  const double p = ps(x1, x2) + ps(x2, x1);
  return p > 0.0 ? me / p : 1.0;
}

// Initial-state (spacelike) emission off a heavy quark in the medium, Q(p) -> Q*(x p) + g((1-x) p): the massive
// quasi-collinear kernel is C_F[(1+x^2)/(1-x) - 2 x m^2/|t|] with |t| = (kt^2 + (1-x)^2 m^2)/(1-x) the
// off-shellness of Q* (Catani, Dittmaier, Trocsanyi).  Ratio to the massless kernel, used to thin JEWEL's
// initial-state medium emissions of heavy quarks (dead cone); 0 if kt^2 < 0.
inline double hfisweight(double x, double t, double m) {
  static const bool nothin = hftest("nothin");
  if (nothin) return 1.0;
  if (t <= (1.0 - x) * m * m) return 0.0;
  return std::clamp(1.0 - 2.0 * x * (1.0 - x) * m * m / ((1.0 + x * x) * t), 0.0, 1.0);
}
// signed mass^2 of a line (P5 < 0 denotes a spacelike parton)
inline double signedm2(double p5) { return p5 < 0.0 ? -p5 * p5 : p5 * p5; }

// massive splitting functions (without alpha_s / 2pi)
inline double pqqmassive(double z, double t, double m) {
  return 4.0 / 3.0 * ((1.0 + z * z) / (1.0 - z) - 2.0 * m * m / t);
}
inline double pgqmassive(double z, double t, double m) {
  return 0.5 * (z * z + (1.0 - z) * (1.0 - z) + 2.0 * m * m / t);
}

}  // namespace jewel
