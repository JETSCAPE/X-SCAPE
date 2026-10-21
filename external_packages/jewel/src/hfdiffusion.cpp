// Heavy-flavour extension: Langevin transport of charm and bottom quarks in the medium (HQDS > 0).
//
// JEWEL's explicit 2->2 scatterings contain the leading-order perturbative elastic energy loss; for slow heavy
// quarks this is too weak to describe B and D suppression (see VALIDATION.md).  This adds a non-perturbative
// diffusion component, as in LIDO / Langevin-based heavy-quark models, characterised by the spatial diffusion
// coefficient 2 pi T D_s = HQDS (temperature independent):
//   momentum diffusion  kappa = 2 T^2 / D_s = 4 pi T^3 / HQDS            (isotropic, kappa_L = kappa_T)
//                       optionally falling off at high momentum (HQP0 = p0 > 0, rest-frame momentum p):
//                       kappa(p) = kappa / (1 + (p / p0)^2)
//   drag                eta_D = kappa / (2 T E) - kappa'(p) / (2 p)      (Einstein relation, pre-point Ito;
//                       the second term = kappa(p) / (p0^2 (1 + (p/p0)^2)) keeps the thermal distribution exact)
// so that heavy quarks relax towards the thermal distribution.  Each step is done in the local fluid rest frame:
//   p -> p - eta_D p dt + xi,   <xi_i xi_j> = kappa dt delta_ij.
// It is applied after the parton evolution (before hadronisation) to every final-state c and b quark, over its
// whole path in the medium: the path is reconstructed from the quark's ancestry (scatterings and Q -> Q g
// splittings of the same quark), and the medium temperature and flow are taken along it; the momentum that
// diffuses is the final one.  The energy exchanged with the medium is not tracked (like soft recoils).
//
// Recombination (HQRECO = N > 0; hadronisation only): instantaneous coalescence of a heavy quark Q with a thermal
// light antiquark at the point where Q leaves the QGP (last step with T > T_c along its path), as in the Duke
// Langevin and LIDO models.  Gaussian meson Wigner function in the pair relative momentum k (pair rest frame),
//   W(k) = exp(-k^2 sigma^2),  sigma^2 = 1/(mu omega),  omega from <r^2> = (3/2)(m_Q^2 + m_q^2)/((m_Q + m_q)^2 mu omega)
// (HQRECOR2 = <r^2> [fm^2], default 0.18), light quarks Boltzmann at the local temperature with constituent masses
// 0.3 (u, d) and 0.475 GeV (s).  The coalescence probability at fluid-frame momentum p*,
//   P(p*) = N sum_f I_f(p*) / sum_f I_f(0),   I_f(p*) = int d^3q exp(-E_q/T) W(k),
// is 1 at rest for N = 1 (all slow heavy quarks recombine).  If Q coalesces, the antiquark flavour is chosen with
// weight I_f, its momentum sampled from exp(-E_q/T) W(k); the meson (pseudoscalar : vector = 1 : 3, D/D*, B/B*;
// no baryons) gets p = p_Q + p_q, E from its mass, and is added after the string fragmentation (PYTHIA decays it
// as usual).  The partner light quark of the thermal q qbar pair, with a thermal momentum, takes Q's place in the
// jet's colour system, so the strings are unchanged in topology.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "jewel/commons.hpp"
#include "jewel/heavyflavour.hpp"
#include "jewel/medium.hpp"
#include "jewel/pythia6.hpp"

namespace jewel {

namespace {

constexpr double kHbarc = 0.1973269804;  // GeV fm
constexpr double kPi = 3.141592653589793;
constexpr double kStep = 0.05;           // lab-frame time step [fm]

struct Meson {
  int kf;
  double p[3];
};
std::vector<Meson> pending;  // coalesced mesons of the current event (added after the string fragmentation)

constexpr double kMq[3] = {0.3, 0.3, 0.475};  // constituent masses of d, u, s [GeV]

// relative momentum of two particles in their rest frame
double krel(double e1, const double p1[3], double m1, double e2, const double p2[3], double m2) {
  const double s = (e1 + e2) * (e1 + e2) - (p1[0] + p2[0]) * (p1[0] + p2[0]) - (p1[1] + p2[1]) * (p1[1] + p2[1]) -
                   (p1[2] + p2[2]) * (p1[2] + p2[2]);
  const double a = s - (m1 + m2) * (m1 + m2), b = s - (m1 - m2) * (m1 - m2);
  return (a > 0.0 && b > 0.0) ? std::sqrt(a * b) / (2.0 * std::sqrt(s)) : 0.0;
}

double sigma2(double mQ, double mq) {  // Wigner-function width^2 [GeV^-2] from <r^2>
  const double r2 = cb::hf::hqrecor2 / (kHbarc * kHbarc);
  const double mu = mQ * mq / (mQ + mq);
  const double omega = 1.5 * (mQ * mQ + mq * mq) / ((mQ + mq) * (mQ + mq) * mu * r2);
  return 1.0 / (mu * omega);
}

// I_f(p*) = int d^3q exp(-E_q/T) W(k) for a heavy quark of mass mQ and momentum ps along z (2D quadrature)
double coalint(double ps, double mQ, double mq, double T) {
  const double s2 = sigma2(mQ, mq), eQ = std::sqrt(ps * ps + mQ * mQ);
  const double pQ[3] = {0.0, 0.0, ps};
  const int nq = 80, nc = 40;
  const double qmax = 15.0 * T + 3.0;
  double sum = 0.0;
  for (int i = 0; i < nq; ++i) {
    const double q = (i + 0.5) * qmax / nq, eq = std::sqrt(q * q + mq * mq);
    double ang = 0.0;
    for (int j = 0; j < nc; ++j) {
      const double c = -1.0 + (j + 0.5) * 2.0 / nc, sn = std::sqrt(1.0 - c * c);
      const double pq[3] = {q * sn, 0.0, q * c};
      const double k = krel(eQ, pQ, mQ, eq, pq, mq);
      ang += std::exp(-k * k * s2) * 2.0 / nc;
    }
    sum += q * q * std::exp(-eq / T) * ang * qmax / nq;
  }
  return sum;
}

double boltzmann(double m, double T) {  // |q| from q^2 exp(-sqrt(q^2+m^2)/T)
  for (;;) {
    const double q = -T * std::log(std::max(py::pyr(0), 1e-300) * std::max(py::pyr(0), 1e-300) *
                                   std::max(py::pyr(0), 1e-300));
    if (py::pyr(0) < std::exp(-(std::sqrt(q * q + m * m) - q) / T)) return q;
  }
}

void isotropic(double q, double v[3]) {
  const double c = 2.0 * py::pyr(0) - 1.0, sn = std::sqrt(1.0 - c * c), ph = 2.0 * kPi * py::pyr(0);
  v[0] = q * sn * std::cos(ph);
  v[1] = q * sn * std::sin(ph);
  v[2] = q * c;
}

double gauss() {  // standard normal (Box-Muller) from PYTHIA's generator
  double u1 = 0.0;
  do {
    u1 = py::pyr(0);
  } while (u1 <= 0.0);
  const double u2 = py::pyr(0);
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
}

// Lorentz boost of the 4-vector (e, p) by velocity b (to the frame moving with -b)
void boost(double b[3], double& e, double p[3]) {
  const double b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
  if (b2 <= 0.0) return;
  const double g = 1.0 / std::sqrt(1.0 - b2);
  const double bp = b[0] * p[0] + b[1] * p[1] + b[2] * p[2];
  const double g2 = (g - 1.0) / b2;
  for (int i = 0; i < 3; ++i) p[i] += g2 * bp * b[i] + g * b[i] * e;
  e = g * (e + bp);
}

}  // namespace

void hfcoalesce(int l, double T, const double b[3]);

// diagnostic (JEWEL_HFTEST=stats): per heavy quark, lab time spent at T > 0, mean T, energy change by diffusion
static void hfstat_diff(int aq, double tin, double tsum, double de) {
  static const bool on = hftest("stats");
  struct S { long n[2] = {0, 0}; double t[2] = {0, 0}, ts[2] = {0, 0}, e[2] = {0, 0};
    ~S() { for (int f = 0; f < 2; ++f) if (n[f])
             std::fprintf(stderr, "HFSTATS diffusion %s: quarks %ld  <time at T>0> = %.3f fm  <T> = %.3f GeV  <dE> = %.3f GeV\n",
                          f ? "b" : "c", n[f], t[f] / n[f], t[f] > 0 ? ts[f] / t[f] : 0.0, e[f] / n[f]); } };
  static S s;
  if (!on) return;
  const int f = (aq == 5) ? 1 : 0;
  ++s.n[f]; s.t[f] += tin; s.ts[f] += tsum; s.e[f] += de;
}

// diagnostic (JEWEL_HFTEST=stats): coalesced fraction vs fluid-frame momentum, for c and b
static void hfstat_reco(int aq, double ps, bool coal) {
  static const bool on = hftest("stats");
  struct S { long n[2][12] = {}, c[2][12] = {};
    ~S() { for (int f = 0; f < 2; ++f) { std::fprintf(stderr, "HFSTATS coalescence %s: p*[GeV] frac (n):", f ? "b" : "c");
             for (int i = 0; i < 12; ++i) if (n[f][i]) std::fprintf(stderr, " %d-%d:%.3f(%ld)", i, i + 1, double(c[f][i]) / n[f][i], n[f][i]);
             std::fprintf(stderr, "\n"); } } };
  static S s;
  if (!on) return;
  const int f = (aq == 5) ? 1 : 0, i = std::min(11, int(ps));
  ++s.n[f][i];
  if (coal) ++s.c[f][i];
}

void hfdiffuse() {
  using namespace py::pyjets;
  using namespace cb::time;
  const double x = cb::hf::hqds, p0 = cb::hf::hqp0;
  // recombination only in the stand-alone program (JEWEL's strings); in library mode (JETSCAPE) the
  // framework's hadronisation (e.g. HybridHadronization) recombines the heavy quarks
  const bool reco = (cb::hf::hqreco > 0.0) && cb::evrecord::hadro && !cb::lib::active;
  pending.clear();
  if ((x <= 0.0) && !reco) return;
  for (int l = 1; l <= N; ++l) {
    // final-state partons: status 1 (still active) or 4 (evolution finished, cascade.cpp)
    if (((K(l, 1) != 1) && (K(l, 1) != 4)) || !isheavy(K(l, 2))) continue;
    // ancestry of the same heavy quark (earliest first); stop at a different flavour or a time-ordering break
    int chain[1000];
    int n = 0;
    chain[n++] = l;
    for (int a = K(l, 3); a > 0 && n < 1000; a = K(a, 3)) {
      if (K(a, 2) != K(l, 2) || mv(a, 4) > mv(chain[n - 1], 4) + 1e-9) break;
      chain[n++] = a;
    }
    double p[3] = {P(l, 1), P(l, 2), P(l, 3)};
    const double m = P(l, 5);
    double pos[3] = {0.0, 0.0, 0.0};
    double tex = 0.0, bex[3] = {0.0, 0.0, 0.0};  // temperature and fluid velocity where Q leaves the QGP
    double tin = 0.0, tsum = 0.0;                  // diagnostics: time at T > 0 and its T-weighted sum
    double xex[4] = {-1.0, 0.0, 0.0, 0.0};          // (t, x, y, z) where Q leaves the QGP
    const double e0 = P(l, 4);
    // segments k = n-1 (earliest) ... 0 (final line); the final segment runs to the end of the medium
    for (int k = n - 1; k >= 0; --k) {
      const int c = chain[k];
      const double t0 = mv(c, 4);
      const double t1 = (k > 0) ? mv(chain[k - 1], 4) : cb::param::ltime;
      if (t1 <= t0) continue;
      double v[3];
      for (int i = 0; i < 3; ++i) {
        pos[i] = mv(c, i + 1);
        v[i] = (k > 0) ? P(c, i + 1) / P(c, 4) : 0.0;
      }
      for (double t = t0; t < t1; t += kStep) {
        const double dt = std::min(kStep, t1 - t);
        const double tm = t + 0.5 * dt;
        if (k == 0) {
          const double e = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + m * m);
          for (int i = 0; i < 3; ++i) v[i] = p[i] / e;
        }
        double xm[3];
        for (int i = 0; i < 3; ++i) xm[i] = pos[i] + 0.5 * dt * v[i];
        const double temp = medium().gettemp(xm[0], xm[1], xm[2], tm);
        if (temp > 0.0) {
          double u[5];
          medium().avscatcen(xm[0], xm[1], xm[2], tm, u[0], u[1], u[2], u[3], u[4]);
          double b[3] = {0.0, 0.0, 0.0};
          if (u[3] > 0.0)
            for (int i = 0; i < 3; ++i) b[i] = u[i] / u[3];
          tex = temp;
          xex[0] = t + dt;
          for (int i = 0; i < 3; ++i) xex[i + 1] = pos[i] + dt * v[i];
          tin += dt;
          tsum += dt * temp;
          for (int i = 0; i < 3; ++i) bex[i] = b[i];
          if (x <= 0.0) {
            for (int i = 0; i < 3; ++i) pos[i] += dt * v[i];
            continue;
          }
          // to the fluid rest frame
          double e = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + m * m);
          const double elab = e;
          double mb[3] = {-b[0], -b[1], -b[2]};
          boost(mb, e, p);
          const double dtr = dt * e / elab;                                   // rest-frame time step [fm]
          double kappa = 4.0 * kPi * temp * temp * temp / x / kHbarc;  // GeV^2 / fm
          double eta = 0.0;                                             // 1 / fm
          if (p0 > 0.0) {
            const double f = 1.0 + (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) / (p0 * p0);
            kappa /= f;
            eta = kappa / (p0 * p0 * f);
          }
          eta += kappa / (2.0 * temp * e);
          const double sig = std::sqrt(kappa * dtr);
          for (int i = 0; i < 3; ++i) p[i] += -eta * p[i] * dtr + sig * gauss();
          e = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + m * m);
          boost(b, e, p);
        }
        for (int i = 0; i < 3; ++i) pos[i] += dt * v[i];
      }
    }
    P(l, 1) = p[0];
    P(l, 2) = p[1];
    P(l, 3) = p[2];
    P(l, 4) = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + m * m);
    hfstat_diff(std::abs(K(l, 2)), tin, tsum, P(l, 4) - e0);
    // library mode (JETSCAPE): the heavy quark is handed to the framework's hadronisation at the point where it
    // leaves the QGP (e.g. HybridHadronization recombines it with thermal partons near that point), not at
    // the vertex of its last interaction inside the medium
    if (cb::lib::active && (xex[0] > 0.0)) {
      mv(l, 4) = xex[0];
      for (int i = 0; i < 3; ++i) mv(l, i + 1) = xex[i + 1];
    }
    if (reco && (tex > 0.0)) hfcoalesce(l, tex, bex);
  }
}

// coalescence of the heavy quark in line l with a thermal antiquark (temperature T, fluid velocity b)
void hfcoalesce(int l, double T, const double b[3]) {
  using namespace py::pyjets;
  const int kq = K(l, 2), sgn = (kq > 0) ? 1 : -1, aq = std::abs(kq);
  const double mQ = (hfmass(kq) > 0.0) ? P(l, 5) : 0.0;
  if (mQ <= 0.0) return;
  // heavy quark in the fluid rest frame
  double pQ[3] = {P(l, 1), P(l, 2), P(l, 3)}, eQ = P(l, 4);
  double mb[3] = {-b[0], -b[1], -b[2]};
  boost(mb, eQ, pQ);
  const double ps = std::sqrt(pQ[0] * pQ[0] + pQ[1] * pQ[1] + pQ[2] * pQ[2]);
  // normalisation sum_f I_f(0), cached per (flavour, temperature)
  static double tcache[2] = {-1.0, -1.0}, norm[2] = {0.0, 0.0};
  const int ic = (aq == 5) ? 1 : 0;
  if (std::fabs(tcache[ic] - T) > 1e-6) {
    norm[ic] = 0.0;
    for (int f = 0; f < 3; ++f) norm[ic] += coalint(0.0, mQ, kMq[f], T);
    tcache[ic] = T;
  }
  double wf[3], wsum = 0.0;
  for (int f = 0; f < 3; ++f) wsum += (wf[f] = coalint(ps, mQ, kMq[f], T));
  const bool coal = py::pyr(0) < cb::hf::hqreco * wsum / norm[ic];
  hfstat_reco(aq, ps, coal);
  if (!coal) return;  // fragments through the strings
  // antiquark flavour (1 d, 2 u, 3 s) and momentum
  double r = py::pyr(0) * wsum;
  int f = 0;
  while (f < 2 && r > wf[f]) r -= wf[f++];
  const double mq = kMq[f], s2 = sigma2(mQ, mq);
  double q[3];
  int ntry = 0;
  for (;; ++ntry) {
    if (ntry > 100000) return;
    isotropic(boltzmann(mq, T), q);
    const double eq = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + mq * mq);
    const double k = krel(std::sqrt(ps * ps + mQ * mQ), pQ, mQ, eq, q, mq);
    if (py::pyr(0) < std::exp(-k * k * s2)) break;
  }
  // meson: Q qbar_f; PDG code (pseudoscalar 1/4, vector 3/4)
  const int fl = f + 1;  // 1 d, 2 u, 3 s
  int kf = 0;
  if (aq == 4) kf = (fl == 2) ? 421 : (fl == 1 ? 411 : 431);    // c ubar = D0, c dbar = D+, c sbar = Ds+
  else kf = -((fl == 2) ? 521 : (fl == 1 ? 511 : 531));        // b ubar = B-, b dbar = B0bar, b sbar = Bs0bar
  if (py::pyr(0) < 0.75) kf += (kf > 0) ? 2 : -2;
  kf *= sgn;
  Meson mes{kf, {0.0, 0.0, 0.0}};
  double eq = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + mq * mq);
  double bb[3] = {b[0], b[1], b[2]};
  boost(bb, eq, q);  // antiquark to the lab
  for (int i = 0; i < 3; ++i) mes.p[i] = P(l, i + 1) + q[i];
  pending.push_back(mes);
  // partner light quark of the thermal pair replaces Q in the colour system
  double qp[3];
  isotropic(boltzmann(mq, T), qp);
  double ep = std::sqrt(qp[0] * qp[0] + qp[1] * qp[1] + qp[2] * qp[2]);  // massless, like JEWEL's light partons
  boost(bb, ep, qp);
  K(l, 2) = sgn * fl;
  for (int i = 0; i < 3; ++i) P(l, i + 1) = qp[i];
  P(l, 4) = ep;
  P(l, 5) = 0.0;
}

// after the string fragmentation: add the coalesced mesons and let PYTHIA decay them
void hfrecoadd() {
  using namespace py::pyjets;
  if (pending.empty()) return;
  for (const Meson& m : pending) {
    ++N;
    K(N, 1) = 1;
    K(N, 2) = m.kf;
    K(N, 3) = 0;
    K(N, 4) = 0;
    K(N, 5) = 0;
    const double mm = py::pymass(m.kf);
    for (int i = 0; i < 3; ++i) P(N, i + 1) = m.p[i];
    P(N, 4) = std::sqrt(m.p[0] * m.p[0] + m.p[1] * m.p[1] + m.p[2] * m.p[2] + mm * mm);
    P(N, 5) = mm;
    for (int i = 1; i <= 5; ++i) cb::time::mv(N, i) = 0.0;
  }
  pending.clear();
  py::pyexec();
}

}  // namespace jewel
