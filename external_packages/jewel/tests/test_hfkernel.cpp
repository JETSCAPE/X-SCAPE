// test_hfkernel: emission density of JEWEL's (massive) Q -> Q g Sudakov veto algorithm in isolation.
//
// A quark of fixed energy E is evolved with getmass() as an eikonal source: every accepted virtuality t is
// one emission, the next one is generated below it (no recoil).  The density of emissions in ln(1/theta),
// theta = kt / (z (1-z) E) the q-g opening angle, is written for a light quark and for the heavy one.
// Their ratio tests the massive kernel against the dead-cone factor (theta^2/(theta^2+theta0^2))^2 times
// the change of the phase-space limits (see analysis/hfplot_ratio.py).
//
// usage: test_hfkernel <E> <m_heavy> <nsources> > out.txt
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "jewel/commons.hpp"
#include "jewel/heavyflavour.hpp"
#include "jewel/routines.hpp"

int main(int argc, char** argv) {
  using namespace jewel;
  const double E = argc > 1 ? std::atof(argv[1]) : 20.0;
  const double m = argc > 2 ? std::atof(argv[2]) : 4.8;
  const long n = argc > 3 ? std::atol(argv[3]) : 200000;
  cb::param::q0 = 1.5;
  cb::param::lps = 0.4;
  cb::param::nf = 3;
  cb::param::angord = 1;
  cb::alphasfac::ptfac = 1.0;
  cb::hf::on = true;
  cb::hf::mc = 1.5;
  cb::hf::mb = 4.8;
  constexpr int nb = 40;
  const double lmin = 0.0, lmax = 5.0;
  std::vector<double> hl(nb, 0.0), hh(nb, 0.0);
  for (int pass = 0; pass < 2; ++pass) {
    const double mq = pass == 0 ? 0.0 : m;
    auto& h = pass == 0 ? hl : hh;
    for (long i = 0; i < n; ++i) {
      double qmax = E, zdec = 0.0;
      bool qq = false;
      for (;;) {
        const double q = getmass(0.0, qmax, -1.0, E, "QQ", qmax, false, zdec, qq, mq);
        if (q <= 0.0) break;
        const double t = q * q, z = zdec;
        const double kt2 = (mq > 0.0) ? hfkt2qq(t, mq, z) : z * (1.0 - z) * t;
        if (kt2 > 0.0) {
          const double theta = std::sqrt(kt2) / (z * (1.0 - z) * E);
          const double l = std::log(1.0 / theta);
          const int ib = static_cast<int>(std::floor((l - lmin) / (lmax - lmin) * nb));
          if (ib >= 0 && ib < nb) h[ib] += 1.0;
        }
        qmax = q;
      }
    }
  }
  std::printf("# E=%g m=%g sources=%ld  columns: ln(1/theta) dN_light/dl dN_heavy/dl ratio ratio_err\n", E, m, n);
  const double dl = (lmax - lmin) / nb;
  for (int b = 0; b < nb; ++b) {
    const double r = hl[b] > 0 ? hh[b] / hl[b] : 0.0;
    const double er = (hl[b] > 0 && hh[b] > 0) ? r * std::sqrt(1.0 / hl[b] + 1.0 / hh[b]) : 0.0;
    std::printf("%.4f %.6g %.6g %.4f %.4f\n", lmin + (b + 0.5) * dl, hl[b] / (n * dl), hh[b] / (n * dl), r, er);
  }
  return 0;
}
