// Check of the JETSCAPE interface of JEWEL (jetscape_api) without JETSCAPE:
// showers single light quarks of energy E0 at the origin, moving along x, in a
// static brick (T for 0 < t <= L, no flow) and reports
//   * the energy given to the medium, E0 - sum E(status 0), to compare with the
//     stand-alone JEWEL brick (medium-simple BRICK T, e+e- -> q qbar at
//     sqrt(s) = 2 E0, recoils not kept: 13.16 +- 0.08 GeV for E0 = 100 GeV,
//     T = 0.3 GeV, L = 4 fm, TC 0.15, HEAVYFLAVOUR T; HEAVYFLAVOUR.md),
//   * energy conservation with recoils and holes:
//     sum E(status 0) + sum E(status 1) - sum E(status -1) = E(in)
//     (the three-momentum is not conserved exactly: the initial virtuality is
//     given to the parton at fixed energy).
//
//   test_jetscape_brick <jewel param file> <nshowers> <E0> <T> <L> <seed> [vac]
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "jewel/jetscape_api.hpp"

int main(int argc, char** argv) {
  if (argc < 7) {
    std::fprintf(stderr, "usage: %s params nshowers E0 T L seed [vac]\n", argv[0]);
    return 1;
  }
  const std::string params = argv[1];
  const int n = std::atoi(argv[2]);
  const double e0 = std::atof(argv[3]), temp = std::atof(argv[4]), len = std::atof(argv[5]);
  const long seed = std::atol(argv[6]);
  const bool vac = argc > 7 && std::string(argv[7]) == "vac";

  jewel::api::Settings s;
  s.param_file = params;
  s.run_param_file = "jewel_params_run.dat";
  s.seed = seed;
  s.medium.tmax = temp;
  s.medium.tc = 0.15;
  s.medium.tau0 = 0.0;
  s.medium.lifetime = len;
  s.medium.etamax = 10.0;
  s.medium.vacuum = vac;
  jewel::api::initialise(s, [temp](double t, double, double, double) {
    jewel::api::MediumCell c;
    if (t > 0.0) c.temperature = temp;
    return c;
  });

  double sum = 0, sum2 = 0, maxviol = 0, sumviol = 0, sumrec = 0, sumhole = 0;
  long nrec = 0, nhole = 0, nok = 0;
  std::vector<jewel::api::OutParton> out;
  for (int i = 0; i < n; ++i) {
    const int id = 1 + (i % 3);  // d, u, s
    const double p[4] = {e0, e0, 0.0, 0.0};
    const double x[4] = {0.0, 0.0, 0.0, 0.0};
    jewel::api::next_event(0.0);
    if (!jewel::api::shower(id, p, x, 2.0 * e0, out)) continue;
    ++nok;
    double esh = 0, bal[4] = {-p[0], -p[1], -p[2], -p[3]};
    for (const auto& o : out) {
      const double w = (o.status == -1) ? -1.0 : 1.0;
      for (int k = 0; k < 4; ++k) bal[k] += w * o.p[k];
      if (o.status == 0) esh += o.p[0];
      if (o.status == 1) {
        ++nrec;
        sumrec += o.p[0];
      }
      if (o.status == -1) {
        ++nhole;
        sumhole += o.p[0];
      }
    }
    const double de = e0 - esh;
    sum += de;
    sum2 += de * de;
    maxviol = std::max(maxviol, std::fabs(bal[0]));
    sumviol += bal[0];
  }
  const double mean = sum / nok, err = std::sqrt((sum2 / nok - mean * mean) / nok);
  const auto& st = jewel::api::stats();
  std::printf("showers %ld (discarded %ld), E0 = %g GeV, T = %g GeV, L = %g fm%s\n", nok, st.n_discarded,
              e0, temp, len, vac ? " (vacuum)" : "");
  std::printf("energy given to the medium: %.3f +- %.3f GeV\n", mean, err);
  std::printf("recoils per shower %.3f (<E> %.3f GeV), holes per shower %.3f (<E> %.3f GeV)\n",
              double(nrec) / nok, nrec ? sumrec / nrec : 0.0, double(nhole) / nok,
              nhole ? sumhole / nhole : 0.0);
  std::printf("E(0) + E(1) - E(-1) - E(in): mean %.3g GeV, max |.| %.3g GeV\n", sumviol / nok, maxviol);
  return 0;
}
