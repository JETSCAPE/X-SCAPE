// Heavy-flavour extension: Langevin transport of charm and bottom quarks in the medium (HQDS > 0).
//
// JEWEL's explicit 2->2 scatterings contain the leading-order perturbative elastic energy loss; for slow heavy
// quarks this is too weak to describe B and D suppression (see VALIDATION.md).  This adds a non-perturbative
// diffusion component, as in LIDO / Langevin-based heavy-quark models, characterised by the spatial diffusion
// coefficient 2 pi T D_s = HQDS (temperature independent):
//   momentum diffusion  kappa = 2 T^2 / D_s = 4 pi T^3 / HQDS            (isotropic, kappa_L = kappa_T)
//   drag                eta_D = kappa / (2 T E)                          (Einstein relation, pre-point Ito)
// so that heavy quarks relax towards the thermal distribution.  Each step is done in the local fluid rest frame:
//   p -> p - eta_D p dt + xi,   <xi_i xi_j> = kappa dt delta_ij.
// It is applied after the parton evolution (before hadronisation) to every final-state c and b quark, over its
// whole path in the medium: the path is reconstructed from the quark's ancestry (scatterings and Q -> Q g
// splittings of the same quark), and the medium temperature and flow are taken along it; the momentum that
// diffuses is the final one.  The energy exchanged with the medium is not tracked (like soft recoils).

#include <cmath>

#include "jewel/commons.hpp"
#include "jewel/heavyflavour.hpp"
#include "jewel/medium.hpp"
#include "jewel/pythia6.hpp"

namespace jewel {

namespace {

constexpr double kHbarc = 0.1973269804;  // GeV fm
constexpr double kPi = 3.141592653589793;
constexpr double kStep = 0.05;           // lab-frame time step [fm]

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

void hfdiffuse() {
  using namespace py::pyjets;
  using namespace cb::time;
  const double x = cb::hf::hqds;
  if (x <= 0.0) return;
  for (int l = 1; l <= N; ++l) {
    if (K(l, 1) != 1 || !isheavy(K(l, 2))) continue;
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
          // to the fluid rest frame
          double e = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + m * m);
          const double elab = e;
          double mb[3] = {-b[0], -b[1], -b[2]};
          boost(mb, e, p);
          const double dtr = dt * e / elab;                                   // rest-frame time step [fm]
          const double kappa = 4.0 * kPi * temp * temp * temp / x / kHbarc;  // GeV^2 / fm
          const double eta = kappa / (2.0 * temp * e);                        // 1 / fm
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
  }
}

}  // namespace jewel
