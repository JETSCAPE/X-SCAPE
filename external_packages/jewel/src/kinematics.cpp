// JEWEL 2.6.0 C++ port -- final-state scattering, new masses, momentum transfers and 2->2 kinematics
// (DOFISTATESCAT, GETNEWMASS, GETPNORAD1, GETQVEC, DOKINEMATICS).
// Translated from jewel-2.6.0.f:5080-6458.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "jewel/heavyflavour.hpp"
#include "jewel/commons.hpp"
#include "jewel/farray.hpp"
#include "jewel/fortran_compat.hpp"
#include "jewel/io.hpp"
#include "jewel/medium.hpp"
#include "jewel/pythia6.hpp"
#include "jewel/routines.hpp"

namespace jewel {


namespace {
const char* fbool(bool b) { return b ? "T" : "F"; }
}  // namespace

// ***********************************************************************
// ***	  subroutine dofistatescat
// ***********************************************************************
void dofistatescat(const int& l, const double& tnow, const double& dtleft, double& deltat,
                   double& newmass, bool& overq0, double& z, bool& qqbar) {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::ftimefac;
  using namespace cb::param;
  using namespace cb::disc;
  using namespace cb::coherent;
  using namespace cb::angor;
  // C--local variables
  constexpr int countmax = 100;  // DATA COUNTMAX/100/
  int count2 = 0;
  double tleft = 0, deltal = 0, q2old = 0, tsum = 0, rtau = 0, lambda = 0, phi = 0, tauest = 0,
         low = 0;
  F1<double, 4> qsumvecold;
  F1<double, 3> pos;
  bool noscat = false;
  [[maybe_unused]] char typ = ' ';

  deltal = 0.0;
  // current virtuality of the line (= P(l,5) for light partons; sqrt(M^2 - m^2) for heavy quarks).  newmass
  // is a virtuality throughout (the target handed to dokinematics).
  const double cv = virtuality(P(l, 5), hfmass(K(l, 2)));

  if (-qsum2 > (1.0 + 1.e-8) * ipow<2>(P(l, 4)))
    io::log(" DOFISTATESCAT has a problem: {} {}", -qsum2, ipow<2>(P(l, 4)));

  if (K(l, 2) == 21) {
    typ = 'G';
  } else {
    typ = 'Q';
  }
  low = ipow<2>(q0) / ipow<2>(scalefacm);

  tsum = 0.0;
  phi = 0.0;
  deltat = 0.0;

  // C--check for radiation with first (given) momentum transfer
  q2old = 0.0;
  if (overq0 || (-qsum2 > low)) {
    const double a2 = ipow<2>(scalefacm) * qsum2;
    const double a3 = ipow<2>(scalefacm) * q2old;
    const double res = getnewmass(l, a2, a3, newmass, false, 1.0, z, qqbar);
    newmass = res;
    overq0 = true;
  } else {
    newmass = cv;
  }

  rtau = py::pyr(0);

  tleft = dtleft;
l222:
  if ((newmass == 0.0) || (newmass == cv)) {
    tauest = tleft;
  } else {
    tauest = std::max(ftfac * (1. - phi) * r4(0.2f) * P(l, 4) / ipow<2>(newmass), 0.0);
  }
  if (tauest > 0.0) {
    lambda = 1.0 / tauest;
    tauest = -std::log(1.0 - rtau) / lambda;
  }
  {
    const double tpos = tnow + tsum;
    getposition(l, tpos, pos(1), pos(2), pos(3));
  }
  {
    const double tst = tnow + tsum;
    const double dtm = std::min(tauest, tleft);
    noscat = !getdeltat(l, tst, dtm, deltal, false);
  }
  if (!noscat) {
    // C--do scattering
    nend = nend + 1;
    if (nend > countmax) {
      nend = nend - 1;
      goto l218;
    }
    if (nstart == 0) nstart = 1;
    tsum = tsum + deltal;
    tleft = tleft - deltal;
    if ((newmass != 0.0) && (newmass != cv)) {
      phi = phi + 5. * deltal * ipow<2>(newmass) / (1. * P(l, 4));
    }
    allqs(nend, 6) = tnow + tsum;
    qsumvecold(1) = qsumvec(1);
    qsumvecold(2) = qsumvec(2);
    qsumvecold(3) = qsumvec(3);
    qsumvecold(4) = qsumvec(4);
    q2old = qsum2;
    // C--get new momentum transfer
    count2 = 0;
  l219: {
    const double dtq = tnow + tsum - mv(l, 4);
    getqvec(l, nend, dtq, 1.0, false);
  }
    if (-qsum2 > ipow<2>(P(l, 4))) {
      qsumvec(1) = qsumvecold(1);
      qsumvec(2) = qsumvecold(2);
      qsumvec(3) = qsumvecold(3);
      qsumvec(4) = qsumvecold(4);
      qsum2 = q2old;
      if (count2 < 100) {
        count2 = count2 + 1;
        goto l219;
      } else {
        allqs(nend, 1) = 0.0;
        allqs(nend, 2) = 0.0;
        allqs(nend, 3) = 0.0;
        allqs(nend, 4) = 0.0;
        allqs(nend, 5) = 0.0;
        nend = nend - 1;
      }
    }
    // C--figure out new virtuality
    if (overq0 || (-qsum2 > low)) {
      const double a2 = ipow<2>(scalefacm) * qsum2;
      const double a3 = ipow<2>(scalefacm) * q2old;
      const double res = getnewmass(l, a2, a3, newmass, false, 1.0, z, qqbar);
      newmass = res;
      overq0 = true;
    }
    goto l222;
  }
  // C--no more scattering
l218:
  if ((ipow<2>(newmass) > low) && (newmass != cv)) {
    if ((tleft < tauest) || (py::pyr(0) > 1.0 / (nend * 1.0))) {
      if (nend == countmax) {
        deltat = tsum;
      } else if (tleft < tauest) {
        deltat = tsum + tleft;
      } else {
        deltat = tsum + tauest;
      }
      newmass = cv;
    } else {
      deltat = tsum + tauest;
    }
  } else {
    deltat = 0.0;
    nstart = 1;
    nend = 1;
    qsum2 = allqs(nend, 1);
    qsumvec(1) = allqs(nend, 2);
    qsumvec(2) = allqs(nend, 3);
    qsumvec(3) = allqs(nend, 4);
    qsumvec(4) = allqs(nend, 5);
    if (-allqs(nend, 1) > low) {
      overq0 = true;
    } else {
      overq0 = false;
    }
    newmass = cv;
  }
  return;
}

// ***********************************************************************
// ***	  function getnewmass
// ***********************************************************************
// Heavy quark (HEAVYFLAVOUR T): same logic as below, but in terms of the virtuality sqrt(t), t = M^2 - m^2
// ('mass' and the result are virtualities), with the kinematic bound sqrt(x^2 E^2 - m^2).
//  * final-state route (in = false): massive Q -> Q g kernel and z limits (getmass with mq)
//  * initial-state route (in = true): JEWEL's integrated in-medium splitting tables, exactly as for light quarks
//    (the mass suppression of these emissions is applied by the thinning in doinstatescat), so that the heavy
//    path reduces to the light one for m -> 0
// The no-splitting probabilities used for pkeep are the massless ones (approximation).
static double getnewmass_heavy(const int& l, const double& q2, const double& qold2, const double& mass,
                               const bool& in, const double& x, double& zdec, bool& qqbardec) {
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::param;
  const double m = hfmass(K(l, 2));
  const double ep = x * P(l, 4);
  const double vmax = std::sqrt(std::max(ep * ep - m * m, 0.0));  // kinematic bound on the virtuality
  const double cv = virtuality(P(l, 5), m);
  double res = 0.0, qa = 1.0, z = 1.0;
  static const bool norad = hftest("norad");  // diagnostic: no medium-induced radiation of heavy quarks
  if (norad || ep < q0 || vmax <= 0.0) {
    zdec = 0.0;
    qqbardec = false;
    return 0.0;
  }
  if (-q2 < ipow<2>(q0)) return 0.0;
  const std::string typ = "QQ";
  if (std::sqrt(-qold2) <= q0) {
    res = getmass(0.0, std::sqrt(-q2), -1.0, ep, typ, in ? vmax : std::sqrt(-q2), in, zdec, qqbardec, m);
    return std::min(res, vmax);
  }
  if (std::max(cv, mass) > 0.0) {
    if (-q2 > -qold2) {
      const double zold = zdec;
      const bool qqbarold = qqbardec;
      const double qtmp = getmass(0.0, std::sqrt(-q2), -1.0, ep, typ, std::sqrt(-q2), in, zdec, qqbardec, m);
      if (qtmp < std::sqrt(-qold2)) {
        res = mass;
        zdec = zold;
        qqbardec = qqbarold;
      } else {
        res = qtmp;
      }
    } else {
      const double pnosplit1 = getsudakov(std::sqrt(-qold2), qa, q0, z, ep, typ, mv(l, 4), in);
      const double pnosplit2 = getsudakov(std::sqrt(-q2), qa, q0, z, ep, typ, mv(l, 4), in);
      const double pkeep = (1. - pnosplit2) / (1. - pnosplit1);
      if (py::pyr(0) < pkeep) {
        if (cv < std::sqrt(-q2)) {
          res = mass;
        } else {
          int ntry = 0;
          do {
            res = getmass(q0, std::sqrt(-q2), -1.0, ep, typ, std::sqrt(-q2), in, zdec, qqbardec, m);
          } while ((res == 0.0) && (++ntry < 1000));
        }
      } else {
        res = 0.0;
        zdec = 0.0;
        qqbardec = false;
      }
    }
  } else {
    if (-q2 > (-qold2) * r4(1.0000001f)) {
      res = getmass(0.0, std::sqrt(-q2), -1.0, ep, typ, vmax, in, zdec, qqbardec, m);
      if (res < std::sqrt(-qold2)) {
        res = 0.0;
        zdec = 0.0;
        qqbardec = false;
      }
    } else {
      res = 0.0;
      zdec = 0.0;
      qqbardec = false;
    }
  }
  return std::min(res, vmax);
}

double getnewmass(const int& l, const double& q2, const double& qold2, const double& mass,
                  const bool& in, const double& x, double& zdec, bool& qqbardec) {
  if (isheavy(py::pyjets::K(l, 2))) return getnewmass_heavy(l, q2, qold2, mass, in, x, zdec, qqbardec);
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::angor;
  using namespace cb::param;
  // C--local variables
  double getnewmass = 0, pnosplit1 = 0, pnosplit2 = 0, z = 0, qa = 0, pkeep = 0, qtmp = 0,
         zold = 0;
  bool qqbarold = false;
  std::string typ;

  if (x * P(l, 4) < q0) {
    getnewmass = 0.0;
    zdec = 0.0;
    qqbardec = false;
    return getnewmass;
  }
  if (-q2 < ipow<2>(q0)) {
    getnewmass = 0.0;
    return getnewmass;
  }
  if (K(l, 2) == 21) {
    typ = "GC";
  } else {
    typ = "QQ";
  }
  if (std::sqrt(-qold2) <= q0) {
    if (in) {
      getnewmass = getmass(0.0, std::sqrt(-q2), -1.0, x * P(l, 4), typ, x * P(l, 4), in, zdec,
                           qqbardec, 0.0, &cb::hf::pendingfl);
    } else {
      getnewmass = getmass(0.0, std::sqrt(-q2), -1.0, P(l, 4), typ, std::sqrt(-q2), in, zdec,
                           qqbardec, 0.0, &cb::hf::pendingfl);
    }
    getnewmass = std::min(getnewmass, x * P(l, 4));
    return getnewmass;
  }
  z = 1.0;
  qa = 1.0;
  if (std::max(P(l, 5), mass) > 0.0) {
    if (-q2 > -qold2) {
      zold = zdec;
      qqbarold = qqbardec;
      qtmp = getmass(0.0, std::sqrt(-q2), -1.0, x * P(l, 4), typ, std::sqrt(-q2), in, zdec,
                     qqbardec, 0.0, &cb::hf::pendingfl);
      if (qtmp < std::sqrt(-qold2)) {
        getnewmass = mass;
        zdec = zold;
        qqbardec = qqbarold;
      } else {
        getnewmass = qtmp;
      }
    } else {
      pnosplit1 = getsudakov(std::sqrt(-qold2), qa, q0, z, x * P(l, 4), typ, mv(l, 4), in);
      pnosplit2 = getsudakov(std::sqrt(-q2), qa, q0, z, x * P(l, 4), typ, mv(l, 4), in);
      pkeep = (1. - pnosplit2) / (1. - pnosplit1);
      if (py::pyr(0) < pkeep) {
        if (P(l, 5) < std::sqrt(-q2)) {
          getnewmass = mass;
        } else {
          do {  // label 55
            getnewmass = getmass(q0, std::sqrt(-q2), -1.0, x * P(l, 4), typ, std::sqrt(-q2), in,
                                 zdec, qqbardec, 0.0, &cb::hf::pendingfl);
          } while ((getnewmass == 0.0) && (x * P(l, 4) > q0));
        }
      } else {
        getnewmass = 0.0;
        zdec = 0.0;
        qqbardec = false;
      }
    }
  } else {
    if (-q2 > (-qold2) * r4(1.0000001f)) {
      getnewmass = getmass(0.0, std::sqrt(-q2), -1.0, x * P(l, 4), typ, x * P(l, 4), in, zdec,
                           qqbardec, 0.0, &cb::hf::pendingfl);
      if (getnewmass < std::sqrt(-qold2)) {
        getnewmass = 0.0;
        zdec = 0.0;
        qqbardec = false;
      }
    } else {
      getnewmass = 0.0;
      zdec = 0.0;
      qqbardec = false;
    }
  }
  getnewmass = std::min(getnewmass, x * P(l, 4));
  return getnewmass;
}

// ***********************************************************************
// ***	  function getpnorad1
// ***********************************************************************
double getpnorad1(const int& line, const double& x, const double& y, const double& z,
                  const double& t) {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace cb::param;
  // C--local variables
  double getpnorad1 = 0, up = 0, low = 0, ccol = 0, sigmatot = 0, ms1 = 0, md1 = 0, shat = 0,
         pcms2 = 0;
  F1<double, 5> avmom;

  md1 = medium().getmd(x, y, z, t);
  medium().avscatcen(x, y, z, t, avmom(1), avmom(2), avmom(3), avmom(4), avmom(5));
  ms1 = avmom(5);
  shat = ipow<2>(avmom(5)) + ipow<2>(P(line, 5)) +
         2. * (avmom(4) * P(line, 4) - avmom(1) * P(line, 1) - avmom(2) * P(line, 2) -
               avmom(3) * P(line, 3));
  pcms2 = ipow<2>(shat + ipow<2>(P(line, 5)) - ipow<2>(ms1)) / (4. * shat) - ipow<2>(P(line, 5));
  up = 4. * pcms2;
  low = ipow<2>(q0) / ipow<2>(scalefacm);
  if ((up <= low) || (P(line, 4) < q0 / scalefacm)) {
    getpnorad1 = 1.0;
    return getpnorad1;
  }
  if (K(line, 2) == 21) {
    ccol = 3. / 2.;
    // C--probability for no initial state radiation
    sigmatot = getsscat(P(line, 4), P(line, 1), P(line, 2), P(line, 3), P(line, 5), 0.0, 'G', 'C',
                        x, y, z, t, 0);
    if (sigmatot == 0.0) {
      getpnorad1 = -1.0;
      return getpnorad1;
    }
    getpnorad1 = (ccol * (scatprimfunc(low, md1) - scatprimfunc(0.0, md1)) +
                  getxsecint(up, md1, "GB")) /
                 sigmatot;
  } else {
    ccol = r4(2.f / 3.f);
    // C--probability for no initial state radiation
    sigmatot = getsscat(P(line, 4), P(line, 1), P(line, 2), P(line, 3), P(line, 5), 0.0, 'Q', 'C',
                        x, y, z, t, 0);
    if (sigmatot == 0.0) {
      getpnorad1 = 1.0;
      return getpnorad1;
    }
    getpnorad1 = (ccol * (scatprimfunc(low, md1) - scatprimfunc(0.0, md1)) +
                  getxsecint(up, md1, "QB")) /
                 sigmatot;
  }
  if ((getpnorad1 < -1.e-4) || (getpnorad1 > 1.0 + 1.e-4)) {
    io::log(" error: P_norad= {} {} {} {} {} {} {}", getpnorad1, P(line, 4), P(line, 5), low, up,
            K(line, 2), md1);
  }
  return getpnorad1;
}

// ***********************************************************************
// ***	  subroutine getqvec
// ***********************************************************************
// diagnostic counters for heavy-quark momentum transfers (JEWEL_HFTEST=stats): calls, |t| > Q0^2, t reset to 0
static void hfstat_qvec(double mt, double maxt) {
  static const bool on = hftest("stats");
  struct S { long n = 0, over = 0, zero = 0; double sumt = 0, summax = 0;
    ~S() { if (n) std::fprintf(stderr, "HFSTATS getqvec heavy: n=%ld |t|>Q0^2: %.4f t->0: %.4f <|t|>=%.3f <maxt>=%.3f\n",
                                 n, double(over) / n, double(zero) / n, sumt / n, summax / n); } };
  static S s;
  if (!on) return;
  if (mt < 0.0) { ++s.zero; return; }
  ++s.n; s.sumt += mt; s.summax += maxt;
  if (mt > ipow<2>(cb::param::q0)) ++s.over;
}

// diagnostic: scatterings applied to b-quark lines in dokinematics (JEWEL_HFTEST=stats)
static void hfstat_kin(int kf, int what, double de) {
  static const bool on = hftest("stats");
  struct S { long acc[2] = {0, 0}, rej[2] = {0, 0}; double de[2] = {0, 0};
    ~S() { const char* nm[2] = {"m=0(light)", "b"};
      for (int i = 0; i < 2; ++i) if (acc[i] + rej[i])
        std::fprintf(stderr, "HFSTATS dokinematics %s: accepted %ld rejected %ld <dE_lab>/scat = %.3f\n", nm[i], acc[i],
                     rej[i], acc[i] ? de[i] / acc[i] : 0.0); } };
  static S s;
  if (!on) return;
  const int i = (std::abs(kf) == 5) ? 1 : (std::abs(kf) <= 3 ? 0 : -1);
  if (i < 0) return;
  if (what) { ++s.acc[i]; s.de[i] += de; } else ++s.rej[i];
}

void getqvec(const int& l, const int& j, const double& dt, const double& x, const bool& recalc) {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::coherent;
  using namespace cb::disc;
  using namespace cb::param;
  // C--local variables
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/
  constexpr int countmax = 1000;            // DATA COUNTMAX/1000/
  int counter = 0;
  [[maybe_unused]] int count2 = 0;
  double xsc = 0, ysc = 0, zsc = 0, tsc = 0, t = 0, pt = 0, maxt = 0, phi2 = 0, phi = 0,
         theta = 0, pt2 = 0, xi = 0, md = 0, shat = 0, pcms2 = 0;
  F1<double, 4> newmom;
  F1<double, 3> beta;
  F1<double, 5> savemom, avmom;
  [[maybe_unused]] char typs = ' ';

  if (j > 10000) {
    discard = true;
    return;
  }

  counter = 0;
  count2 = 0;

  if (recalc) {
    if (j != 1) io::log(" ERROR in GETQVEC: Nend != 1 for recalc = true");
    qsumvec(1) = 0.0;
    qsumvec(2) = 0.0;
    qsumvec(3) = 0.0;
    qsumvec(4) = 0.0;
    qsum2 = 0.0;
  }

  if (recalc) {
    xsc = scatcentres(j, 7);
    ysc = scatcentres(j, 8);
    zsc = scatcentres(j, 9);
    tsc = scatcentres(j, 10);
  } else {
    xsc = mv(l, 1) + dt * P(l, 1) / P(l, 4);
    ysc = mv(l, 2) + dt * P(l, 2) / P(l, 4);
    zsc = mv(l, 3) + dt * P(l, 3) / P(l, 4);
    tsc = mv(l, 4) + dt;
  }
  md = medium().getmd(xsc, ysc, zsc, tsc);

  medium().avscatcen(xsc, ysc, zsc, tsc, avmom(1), avmom(2), avmom(3), avmom(4), avmom(5));

  for (int i = 1; i <= 5; ++i) {
    savemom(i) = P(l, i);
  }

  xi = std::sqrt(std::max(ipow<2>(x) * ipow<2>(P(l, 4)), ipow<2>(P(l, 5))) - ipow<2>(P(l, 5))) /
       py::pyp(l, 8);
  P(l, 1) = xi * P(l, 1);
  P(l, 2) = xi * P(l, 2);
  P(l, 3) = xi * P(l, 3);
  P(l, 4) = std::max(x * P(l, 4), P(l, 5));

  shat = ipow<2>(avmom(5)) + ipow<2>(savemom(5)) +
         2. * (avmom(4) * savemom(4) - avmom(1) * savemom(1) - avmom(2) * savemom(2) -
               avmom(3) * savemom(3));
  pcms2 = ipow<2>(shat + ipow<2>(savemom(5)) - ipow<2>(avmom(5))) / (4. * shat) -
          ipow<2>(savemom(5));
  maxt = 4. * pcms2;
  if (recalc) maxt = std::min(4. * pcms2, ipow<2>(P(l, 4)));

l444:
  if (recalc) {
    K(1, 2) = fint(scatcentres(j, 1));
    P(1, 1) = scatcentres(j, 2);
    P(1, 2) = scatcentres(j, 3);
    P(1, 3) = scatcentres(j, 4);
    P(1, 4) = scatcentres(j, 5);
    P(1, 5) = scatcentres(j, 6);
    mv(1, 1) = scatcentres(j, 7);
    mv(1, 2) = scatcentres(j, 8);
    mv(1, 3) = scatcentres(j, 9);
    mv(1, 4) = scatcentres(j, 10);
  } else {
    medium().getscatterer(xsc, ysc, zsc, tsc, K(1, 2), P(1, 1), P(1, 2), P(1, 3), P(1, 4),
                          P(1, 5));
    mv(1, 1) = xsc;
    mv(1, 2) = ysc;
    mv(1, 3) = zsc;
    mv(1, 4) = tsc;
    typs = 'Q';
    if (K(1, 2) == 21) typs = 'G';

    scatcentres(j, 1) = K(1, 2);
    scatcentres(j, 2) = P(1, 1);
    scatcentres(j, 3) = P(1, 2);
    scatcentres(j, 4) = P(1, 3);
    scatcentres(j, 5) = P(1, 4);
    scatcentres(j, 6) = P(1, 5);
    scatcentres(j, 7) = mv(1, 1);
    scatcentres(j, 8) = mv(1, 2);
    scatcentres(j, 9) = mv(1, 3);
    scatcentres(j, 10) = mv(1, 4);
  }
  K(1, 1) = 1;
  // C--transform to scattering centre's rest frame and rotate such that parton momentum is in z-direction
  beta(1) = P(1, 1) / P(1, 4);
  beta(2) = P(1, 2) / P(1, 4);
  beta(3) = P(1, 3) / P(1, 4);
  py::pyrobo(l, l, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
  py::pyrobo(1, 1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
  theta = py::pyp(l, 13);
  phi = py::pyp(l, 15);
  py::pyrobo(l, l, 0.0, -phi, 0.0, 0.0, 0.0);
  py::pyrobo(1, 1, 0.0, -phi, 0.0, 0.0, 0.0);
  py::pyrobo(l, l, -theta, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(1, 1, -theta, 0.0, 0.0, 0.0, 0.0);
  // C--pick a t from differential scattering cross section
  const bool hq = isheavy(K(l, 2));
  if (hq && !recalc) {
    // heavy quark: kinematic limit |t| <= 4 p_cms^2 of the actual (thermal) scattering centre, which is at rest
    // here, instead of JEWEL's average centre at rest with mass m_s.  For slow heavy projectiles the average-centre
    // bound is far too small (b at E = 10 GeV, T = 0.3 GeV: 2.9 vs ~7 GeV^2); light partons are unchanged.
    const double ms = P(1, 5), mp = P(l, 5);
    const double sh = ms * ms + mp * mp + 2.0 * ms * P(l, 4);
    maxt = std::max(ipow<2>(sh + mp * mp - ms * ms) / (4.0 * sh) - mp * mp, 0.0) * 4.0;
  }
  if (recalc) {
    t = std::max(allqs(j, 1), -maxt);
  } else {
    t = -gett(0.0, maxt, md);
  }
  if (hq && !recalc) hfstat_qvec(-t, maxt);
l202:
  newmom(4) = P(l, 4) + t / (2. * P(1, 5));
  newmom(3) = (t - 2. * ipow<2>(P(l, 5)) + 2. * P(l, 4) * newmom(4)) / (2. * P(l, 3));
  pt2 = ipow<2>(newmom(4)) - ipow<2>(newmom(3)) - ipow<2>(P(l, 5));
  if (std::fabs(pt2) < 1.e-10) pt2 = 0.0;
  if (t == 0.0) pt2 = 0.0;
  if (pt2 < 0.0) {
    if (hq && !recalc) hfstat_qvec(-1.0, 0.0);
    t = 0.0;
    goto l202;
  }
  pt = std::sqrt(pt2);
  phi2 = py::pyr(0) * 2 * pi;
  newmom(1) = pt * gf_cos(phi2);
  newmom(2) = pt * gf_sin(phi2);
  P(1, 1) = newmom(1) - P(l, 1);
  P(1, 2) = newmom(2) - P(l, 2);
  P(1, 3) = newmom(3) - P(l, 3);
  P(1, 4) = newmom(4) - P(l, 4);
  P(1, 5) = 0.0;
  // C--transformation to lab
  py::pyrobo(l, l, theta, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(1, 1, theta, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(l, l, 0.0, phi, 0.0, 0.0, 0.0);
  py::pyrobo(1, 1, 0.0, phi, 0.0, 0.0, 0.0);
  py::pyrobo(l, l, 0.0, 0.0, beta(1), beta(2), beta(3));
  py::pyrobo(1, 1, 0.0, 0.0, beta(1), beta(2), beta(3));
  allqs(j, 1) = t;
  allqs(j, 2) = P(1, 1);
  allqs(j, 3) = P(1, 2);
  allqs(j, 4) = P(1, 3);
  allqs(j, 5) = P(1, 4);
  // NOTE(port): the Fortran adds ALLQS(NEND,*) here (not ALLQS(J,*)); kept as is.
  qsumvec(1) = qsumvec(1) + allqs(nend, 2);
  qsumvec(2) = qsumvec(2) + allqs(nend, 3);
  qsumvec(3) = qsumvec(3) + allqs(nend, 4);
  qsumvec(4) = qsumvec(4) + allqs(nend, 5);
  qsum2 = ipow<2>(qsumvec(4)) - ipow<2>(qsumvec(1)) - ipow<2>(qsumvec(2)) - ipow<2>(qsumvec(3));
  if (qsum2 > 0.0) {
    qsumvec(1) = qsumvec(1) - allqs(nend, 2);
    qsumvec(2) = qsumvec(2) - allqs(nend, 3);
    qsumvec(3) = qsumvec(3) - allqs(nend, 4);
    qsumvec(4) = qsumvec(4) - allqs(nend, 5);
    qsum2 = ipow<2>(qsumvec(4)) - ipow<2>(qsumvec(1)) - ipow<2>(qsumvec(2)) - ipow<2>(qsumvec(3));
    if (counter > countmax) {
      io::log(" GETQVEC unable to find q vector");
      allqs(j, 1) = 0.0;
      allqs(j, 2) = 0.0;
      allqs(j, 3) = 0.0;
      allqs(j, 4) = 0.0;
      allqs(j, 5) = 0.0;
    } else {
      counter = counter + 1;
      goto l444;
    }
  }
  for (int i = 1; i <= 5; ++i) {
    P(l, i) = savemom(i);
  }
  K(1, 1) = 13;
}

// ***********************************************************************
// ***	  subroutine dokinematics
// ***********************************************************************
void dokinematics(const int& l, [[maybe_unused]] const int& lold, const int& n1, int& n2,
                  const double& newm, bool& retrysplit, const double& time, const double& x,
                  double& z, const bool& qqbar, const bool& forcohscat) {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::ftimefac;
  using namespace cb::colour;
  using namespace cb::param;
  using namespace cb::disc;
  using namespace cb::angor;
  using namespace cb::coherent;
  using namespace cb::weight;
  using namespace cb::storescatcen;
  using namespace cb::storedummies;
  // C--local variables
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/
  constexpr double pdummy = 1.e-6;          // data pdummy/1.d-6/
  int line = 0, j = 0, j_end = 0, dir = 0, nold = 0, colmaxold = 0, statold = 0, nscatcenold = 0,
      lpart = 0;
  double theta = 0, phi = 0, phi2 = 0, maxt = 0, t = 0, newmass = 0, deltam = 0, dm = 0, ttot = 0,
         dmleft = 0, lambda = 0, endtime = 0, tmp = 0, m32 = 0, newm2 = 0, shat = 0, theta2 = 0,
         e3new = 0, e4new = 0, p32 = 0, p42 = 0, p3old = 0, mass2 = 0, enew = 0, pt2 = 0, pt = 0,
         pl = 0, m12 = 0, firsttime = 0, ys = 0, p3boost = 0, pboost = 0, m42 = 0, localt = 0,
         oldt = 0, qmass2 = 0;
  double m4 = 0, z4 = 0;
  F1<double, 3> beta;
  std::string typ;
  bool qqbardec = false, rejectt = false, reshuffle = false, softrec = false, splitrec = false,
       isrecoil = false;
  [[maybe_unused]] bool redokin = false;

  lpart = l;

  if ((N + 2 * (n2 - n1 + 1)) > 22990) {
    io::log(" event too long for event record");
    discard = true;
    return;
  }

  if (K(l, 1) == 2) {
    isrecoil = true;
  } else {
    isrecoil = false;
  }

  firsttime = mv(l, 5);

  redokin = false;

  // heavy quark (HEAVYFLAVOUR T): newm is the target virtuality; internally the target is the invariant
  // mass.  vir() / tgt() convert mass -> virtuality / virtuality -> mass; both are the identity for m = 0.
  const double mqk = hfmass(K(l, 2));
  auto vir = [mqk](double p5) { return (mqk > 0.0) ? virtuality(p5, mqk) : p5; };
  auto tgt = [mqk](double v) { return (mqk > 0.0) ? invmass(v, mqk) : v; };

  newm2 = tgt(newm);
  static const bool hfdbg = std::getenv("JEWEL_HFDEBUG") != nullptr;
  int nretry = 0;
  nold = N;
  nscatcenold = nscatcen;
  colmaxold = colmax;
  statold = K(l, 1);
l204:
  deltam = newm2 - P(l, 5);
  dmleft = deltam;

  ttot = 0.0;
  for (int jj = n1, jj_end = n2; jj <= jj_end; ++jj) {
    ttot = ttot + allqs(jj, 1);
  }

  line = l;

  for (j = n1, j_end = n2; j <= j_end; ++j) {
    splitrec = false;
    // C--projectile type
    if (K(line, 2) == 21) {
      typ = "GC";
      if (py::pyr(0) < 0.5) {
        dir = 1;
      } else {
        dir = -1;
      }
    } else {
      typ = "QQ";
      dir = 0;
    }
    K(1, 1) = 6;
    K(1, 2) = fint(scatcentres(j, 1));
    P(1, 1) = scatcentres(j, 2);
    P(1, 2) = scatcentres(j, 3);
    P(1, 3) = scatcentres(j, 4);
    P(1, 4) = scatcentres(j, 5);
    P(1, 5) = scatcentres(j, 6);
    mv(1, 1) = scatcentres(j, 7);
    mv(1, 2) = scatcentres(j, 8);
    mv(1, 3) = scatcentres(j, 9);
    mv(1, 4) = scatcentres(j, 10);
    t = allqs(j, 1);
    if (t == 0.0) {
      rejectt = true;
      continue;  // goto 222
    } else {
      rejectt = false;
    }

    if (ttot == 0.0) {
      dm = 0.0;
    } else {
      if (dmleft < 0.0) {
        dm = std::max(dmleft * t / ttot * 1.5, dmleft);
      } else {
        dm = std::min(dmleft * t / ttot * 1.5, dmleft);
      }
    }
    ttot = ttot - allqs(j, 1);

    // C--transform to c.m.s. and rotate such that parton momentum is in z-direction
    beta(1) = (P(1, 1) + P(line, 1)) / (P(1, 4) + P(line, 4));
    beta(2) = (P(1, 2) + P(line, 2)) / (P(1, 4) + P(line, 4));
    beta(3) = (P(1, 3) + P(line, 3)) / (P(1, 4) + P(line, 4));
    if ((beta(1) > 1.0) || (beta(2) > 1.0) || (beta(3) > 1.0) ||
        (std::sqrt(ipow<2>(beta(1)) + ipow<2>(beta(2)) + ipow<2>(beta(3))) > 1.0)) {
      reshuffle = false;
    } else {
      reshuffle = true;
    }
  l205:
    if (!reshuffle) {
      beta(1) = P(1, 1) / P(1, 4);
      beta(2) = P(1, 2) / P(1, 4);
      beta(3) = P(1, 3) / P(1, 4);
      py::pyrobo(line, line, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
      py::pyrobo(1, 1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
      theta = py::pyp(line, 13);
      phi = py::pyp(line, 15);
      py::pyrobo(line, line, 0.0, -phi, 0.0, 0.0, 0.0);
      py::pyrobo(1, 1, 0.0, -phi, 0.0, 0.0, 0.0);
      py::pyrobo(line, line, -theta, 0.0, 0.0, 0.0, 0.0);
      py::pyrobo(1, 1, -theta, 0.0, 0.0, 0.0, 0.0);

      if (kinmode == 0) {
        m42 = 0.0;
      } else if (kinmode == 1) {
        m42 = ipow<2>(P(1, 5));
      } else {
        if (scalefacm * std::sqrt(-t) > q0) {
          m4 = getmass(0.0, scalefacm * std::sqrt(-t), -1.0, P(1, 4), typ, P(1, 4), false, z4,
                       qqbardec);
          if (m4 > 0.0) splitrec = true;
          m42 = ipow<2>(m4);
        } else {
          m42 = ipow<2>(P(1, 5));
        }
      }
      if (t == 0.0) m42 = ipow<2>(P(1, 5));
      maxt = -2. * P(1, 5) * P(line, 4) - ipow<2>(P(1, 5)) + m42;
      if (t < maxt) {
        t = 0.0;
        rejectt = true;
        dm = 0.0;
        m42 = ipow<2>(P(1, 5));
        // C--transformation to lab
        py::pyrobo(line, line, theta, 0.0, 0.0, 0.0, 0.0);
        py::pyrobo(line, line, 0.0, phi, 0.0, 0.0, 0.0);
        py::pyrobo(line, line, 0.0, 0.0, beta(1), beta(2), beta(3));
        py::pyrobo(1, 1, theta, 0.0, 0.0, 0.0, 0.0);
        py::pyrobo(1, 1, 0.0, phi, 0.0, 0.0, 0.0);
        py::pyrobo(1, 1, 0.0, 0.0, beta(1), beta(2), beta(3));
        continue;  // goto 222
      }
      // (JEWEL: -P5^2, correct for spacelike projectiles only; heavy quarks use the signed mass^2)
      m12 = (mqk > 0.0) ? signedm2(P(line, 5)) : -ipow<2>(P(line, 5));
    l203:
      newmass = P(line, 5) + dm;
      if (newmass < 0.0) {
        m32 = -ipow<2>(newmass);
      } else {
        m32 = ipow<2>(newmass);
      }
      if ((deltam == 0.0) && isrecoil && (kinmode == 3)) {
        localt = medium().gettemp(mv(1, 1), mv(1, 2), mv(1, 3), mv(1, 4));
        oldt = medium().gettemp(mv(l, 1), mv(l, 2), mv(l, 3), mv(l, 4));
        if (localt > 0.0) {
          m32 = ipow<2>(P(l, 5) * localt / oldt);
          newm2 = std::sqrt(m32);
        }
      }
      if (t == 0.0) {
        enew = P(line, 4);
      } else {
        enew = P(line, 4) + (t + ipow<2>(P(1, 5)) - m42) / (2. * P(1, 5));
      }
      pl = (t + 2. * P(line, 4) * enew - m12 - m32) / (2. * P(line, 3));
      pt2 = ipow<2>(enew) - ipow<2>(pl) - m32;
      if (t == 0.0) pt2 = 0.0;
      if (std::fabs(pt2) < 1.e-8) pt2 = 0.0;
      if (pt2 < 0.0) {
        if (splitrec) {
          m4 = getmass(0.0, std::sqrt(m42), -1.0, P(1, 4), typ, P(1, 4), false, z4, qqbardec);
          if (m4 == 0.0) splitrec = false;
          m42 = ipow<2>(m4);
          goto l203;
        }
        if (dm != 0.0) {
          dm = 0.0;
          goto l203;
        } else {
          io::log(" This should not have happened: pt^2<0!");
          io::log(" {} {} {} {}", t, enew, pl, pt2);
          t = 0.0;
          m42 = ipow<2>(P(1, 5));
          rejectt = true;
          // C--transformation to lab
          py::pyrobo(line, line, theta, 0.0, 0.0, 0.0, 0.0);
          py::pyrobo(line, line, 0.0, phi, 0.0, 0.0, 0.0);
          py::pyrobo(line, line, 0.0, 0.0, beta(1), beta(2), beta(3));
          py::pyrobo(1, 1, theta, 0.0, 0.0, 0.0, 0.0);
          py::pyrobo(1, 1, 0.0, phi, 0.0, 0.0, 0.0);
          py::pyrobo(1, 1, 0.0, 0.0, beta(1), beta(2), beta(3));
          continue;  // goto 222
        }
      }
      pt = std::sqrt(pt2);
      phi2 = py::pyr(0) * 2. * pi;
      N = N + 2;
      lpart = N;
      P(N, 1) = pt * gf_cos(phi2);
      P(N, 2) = pt * gf_sin(phi2);
      P(N, 3) = pl;
      P(N, 4) = enew;
      P(N, 5) = fsign(std::sqrt(std::fabs(m32)), newmass);
      // !---------------------------------
      P(N - 1, 1) = P(1, 1) + P(line, 1) - P(N, 1);
      P(N - 1, 2) = P(1, 2) + P(line, 2) - P(N, 2);
      P(N - 1, 3) = P(1, 3) + P(line, 3) - P(N, 3);
      P(N - 1, 4) = P(1, 4) + P(line, 4) - P(N, 4);
      mass2 = ipow<2>(P(N - 1, 4)) - ipow<2>(P(N - 1, 1)) - ipow<2>(P(N - 1, 2)) -
              ipow<2>(P(N - 1, 3));
      // NOTE(port): Fortran reads `mass2.gt.-1.-6`, i.e. -1.-6 = -7.0 (probably meant -1.d-6); kept.
      if ((mass2 < 0.0) && (mass2 > -1.0 - 6)) mass2 = 0.0;
      if (mass2 < 0.0)
        io::log(" messed up scattering centres mass^2:  {} {}", mass2, ipow<2>(P(1, 5)));
      P(N - 1, 5) = std::sqrt(mass2);
      if (std::fabs(P(N - 1, 5) - std::sqrt(m42)) > 1.e-6)
        io::log(" messed up scattering centres mass (no rs):  {} {} {} {} {}", P(N - 1, 5),
                P(1, 5), P(l, 5), std::sqrt(m42), fbool(rejectt));
      io::log_flush();
      // !---------------------------------
      // !        P(N-1,1)=P(1,1)
      // !        P(N-1,2)=P(1,2)
      // !        P(N-1,3)=P(1,3)
      // !        P(N-1,4)=P(1,4)
      // !        P(N-1,5)=P(1,5)
      // !---------------------------------
    } else {
      py::pyrobo(line, line, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
      py::pyrobo(1, 1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
      if ((P(1, 4) < 0.0) || (P(line, 4) < 0.0)) {
        py::pyrobo(1, 1, 0.0, 0.0, beta(1), beta(2), beta(3));
        py::pyrobo(line, line, 0.0, 0.0, beta(1), beta(2), beta(3));
        reshuffle = false;
        goto l205;
      }
      theta = py::pyp(line, 13);
      phi = py::pyp(line, 15);
      py::pyrobo(line, line, 0.0, -phi, 0.0, 0.0, 0.0);
      py::pyrobo(1, 1, 0.0, -phi, 0.0, 0.0, 0.0);
      py::pyrobo(line, line, -theta, 0.0, 0.0, 0.0, 0.0);
      py::pyrobo(1, 1, -theta, 0.0, 0.0, 0.0, 0.0);
      shat = ipow<2>(P(1, 4) + P(line, 4));
      p3old = P(line, 3);

      maxt = -4. * ipow<2>(P(line, 3));
      if (t < maxt) {
        t = 0.0;
        rejectt = true;
        // C--transformation to lab
        py::pyrobo(line, line, theta, 0.0, 0.0, 0.0, 0.0);
        py::pyrobo(line, line, 0.0, phi, 0.0, 0.0, 0.0);
        py::pyrobo(line, line, 0.0, 0.0, beta(1), beta(2), beta(3));
        py::pyrobo(1, 1, theta, 0.0, 0.0, 0.0, 0.0);
        py::pyrobo(1, 1, 0.0, phi, 0.0, 0.0, 0.0);
        py::pyrobo(1, 1, 0.0, 0.0, beta(1), beta(2), beta(3));
        continue;  // goto 222
      }
      theta2 = std::acos(1.0 + t / (2. * ipow<2>(P(line, 3))));
      phi2 = py::pyr(0) * 2. * pi;
      N = N + 2;
      lpart = N;
      P(N, 1) = P(line, 3) * gf_sin(theta2) * gf_cos(phi2);
      P(N, 2) = P(line, 3) * gf_sin(theta2) * gf_sin(phi2);
      P(N, 3) = P(line, 3) * gf_cos(theta2);
      P(N, 4) = P(line, 4);
      P(N, 5) = P(line, 5);
      // !---------------------------------
      P(N - 1, 1) = P(1, 1) + P(line, 1) - P(N, 1);
      P(N - 1, 2) = P(1, 2) + P(line, 2) - P(N, 2);
      P(N - 1, 3) = P(1, 3) + P(line, 3) - P(N, 3);
      P(N - 1, 4) = P(1, 4) + P(line, 4) - P(N, 4);
      mass2 = ipow<2>(P(N - 1, 4)) - ipow<2>(P(N - 1, 1)) - ipow<2>(P(N - 1, 2)) -
              ipow<2>(P(N - 1, 3));
      // NOTE(port): Fortran reads `mass2.gt.-1.-6`, i.e. -1.-6 = -7.0 (probably meant -1.d-6); kept.
      if ((mass2 < 0.0) && (mass2 > -1.0 - 6)) mass2 = 0.0;
      if (mass2 < 0.0)
        io::log(" messed up scattering centres mass^2:  {} {}", mass2, ipow<2>(P(1, 5)));
      P(N - 1, 5) = std::sqrt(mass2);
      if (std::fabs(P(N - 1, 5) - P(1, 5)) > 1.e-6)
        io::log(" messed up scattering centres mass:  {} {} {}", P(N - 1, 5), P(1, 5), P(l, 5));
      io::log_flush();
      // !---------------------------------
      // !        P(N-1,1)=P(1,1)
      // !        P(N-1,2)=P(1,2)
      // !        P(N-1,3)=P(1,3)
      // !        P(N-1,4)=P(1,4)
      // !        P(N-1,5)=P(1,5)
      // !---------------------------------
    }
    // C--outgoing projectile
    K(N, 1) = K(line, 1);
    if (isrecoil && (newm > vir(P(l, 5))) && (vir(P(N, 5)) > q0) && (j == n2)) {
      K(N, 1) = 1;
    }
    K(N, 2) = K(line, 2);
    // !	 K(N,3)=L
    K(N, 3) = line;
    K(N, 4) = 0;
    K(N, 5) = 0;
    // C--take care of incoming projectile
    if (K(line, 1) == 1) {
      K(line, 1) = 12;
    } else {
      K(line, 1) = 14;
    }
    K(line, 4) = N - 1;
    K(line, 5) = N;
    // C--temporary status code, will be overwritten later
    K(N - 1, 1) = 3;
    K(N - 1, 2) = 21;
    K(N - 1, 3) = 0;
    K(N - 1, 4) = 0;
    K(N - 1, 5) = 0;

    if (reshuffle) {
      // C--adjust mass and re-shuffle momenta

      if (kinmode == 0) {
        m42 = 0.0;
      } else if (kinmode == 1) {
        m42 = ipow<2>(P(1, 5));
      } else {
        if (scalefacm * std::sqrt(-t) > q0) {
          m4 = getmass(0.0, scalefacm * std::sqrt(-t), -1.0, P(1, 4), typ, P(1, 4), false, z4,
                       qqbardec);
          if (m4 > 0.0) splitrec = true;
          m42 = ipow<2>(m4);
        } else {
          m42 = ipow<2>(P(1, 5));
        }
      }
    l206:
      newmass = P(N, 5) + dm;
      if (newmass < 0.0) {
        m32 = -ipow<2>(newmass);
      } else {
        m32 = ipow<2>(newmass);
      }
      if ((deltam == 0.0) && isrecoil && (kinmode == 3)) {
        localt = medium().gettemp(mv(1, 1), mv(1, 2), mv(1, 3), mv(1, 4));
        oldt = medium().gettemp(mv(l, 1), mv(l, 2), mv(l, 3), mv(l, 4));
        if (localt > 0.0) {
          m32 = ipow<2>(P(l, 5) * localt / oldt);
          newm2 = std::sqrt(m32);
        }
      }
      if (t == 0.0) m42 = ipow<2>(P(1, 5));
      e3new = (shat + m32 - m42) / (2.0 * std::sqrt(shat));
      e4new = (shat - m32 + m42) / (2.0 * std::sqrt(shat));
      p32 = ipow<2>(e3new) - m32;
      p42 = ipow<2>(e4new) - m42;
      if ((p32 < 0.0) || (p42 < 0.0) || (e3new < 0.0) || (e4new < 0.0)) {
        if (m42 == 0.0) {
          p42 = 1.e-4;
        } else {
          p42 = 0.0;
        }
        e4new = std::sqrt(p42 + m42);
        e3new = std::sqrt(shat) - e4new;
        p32 = ipow<2>(e4new) - m42;
        m32 = ipow<2>(e3new) - ipow<2>(e4new) + m42;
        if ((e3new < 0.0) || (e4new < 0.0)) {
          if (splitrec) {
            m4 = getmass(0.0, std::sqrt(m42), -1.0, P(1, 4), typ, P(1, 4), false, z4, qqbardec);
            if (m4 == 0.0) splitrec = false;
            m42 = ipow<2>(m4);
            goto l206;
          }
          if (dm != 0.0) {
            dm = 0.0;
            goto l206;
          }
          m42 = ipow<2>(P(1, 5));
          e3new = P(N, 4);
          e4new = P(N - 1, 4);
          p32 = ipow<2>(p3old);
          p42 = ipow<2>(p3old);
          if (P(N, 5) < 0.0) {
            m32 = -ipow<2>(P(N, 5));
          } else {
            m32 = ipow<2>(P(N, 5));
          }
        }
      }
      P(N, 1) = std::sqrt(p32) * P(N, 1) / p3old;
      P(N, 2) = std::sqrt(p32) * P(N, 2) / p3old;
      P(N, 3) = std::sqrt(p32) * P(N, 3) / p3old;
      P(N, 4) = e3new;
      P(N, 5) = fsign(std::sqrt(std::fabs(m32)), newmass);
      tmp = ipow<2>(P(N, 4)) - ipow<2>(P(N, 1)) - ipow<2>(P(N, 2)) - ipow<2>(P(N, 3)) - m32;
      if (std::fabs(tmp) > 1.e-6)
        io::log(" Oups, messed up projectiles mass (rs): {} {} {} {} {} {}", tmp, m32, P(N, 5), dm,
                m42, p32);
      // !---------------------------------
      P(N - 1, 1) = std::sqrt(p42) * P(N - 1, 1) / p3old;
      P(N - 1, 2) = std::sqrt(p42) * P(N - 1, 2) / p3old;
      P(N - 1, 3) = std::sqrt(p42) * P(N - 1, 3) / p3old;
      P(N - 1, 4) = e4new;
      P(N - 1, 5) = std::sqrt(m42);
      tmp = ipow<2>(P(N - 1, 4)) - ipow<2>(P(N - 1, 1)) - ipow<2>(P(N - 1, 2)) -
            ipow<2>(P(N - 1, 3)) - ipow<2>(P(N - 1, 5));
      if (std::fabs(tmp) > 1.e-6)
        io::log(" Oups, messed up scattering centres mass (rs): {} {} {} {} {} {} {}", tmp, p3old,
                P(N - 1, 1), P(N - 1, 2), P(N - 1, 3), P(N - 1, 4), P(N - 1, 5));
      if ((std::fabs(P(N, 1) + P(N - 1, 1)) > 1.e-6) || (std::fabs(P(N, 2) + P(N - 1, 2)) > 1.e-6) ||
          (std::fabs(P(N, 3) + P(N - 1, 3)) > 1.e-6)) {
        io::log(" Oups, momentum not conserved (rs) {} {} {}", P(N, 1) + P(N - 1, 1),
                P(N, 2) + P(N - 1, 2), P(N, 3) + P(N - 1, 3));
        io::log(" {} {} {} {}", m42, dm, e3new, e4new);
      }
      // !---------------------------------
      // !        P(N-1,1)=P(1,1)
      // !        P(N-1,2)=P(1,2)
      // !        P(N-1,3)=P(1,3)
      // !        P(N-1,4)=P(1,4)
      // !        P(N-1,5)=P(1,5)
      // !---------------------------------
    }

    // C--transformation to lab
    py::pyrobo(N - 1, N, theta, 0.0, 0.0, 0.0, 0.0);
    py::pyrobo(line, line, theta, 0.0, 0.0, 0.0, 0.0);
    py::pyrobo(N - 1, N, 0.0, phi, 0.0, 0.0, 0.0);
    py::pyrobo(line, line, 0.0, phi, 0.0, 0.0, 0.0);
    py::pyrobo(N - 1, N, 0.0, 0.0, beta(1), beta(2), beta(3));
    py::pyrobo(line, line, 0.0, 0.0, beta(1), beta(2), beta(3));
    py::pyrobo(1, 1, theta, 0.0, 0.0, 0.0, 0.0);
    py::pyrobo(1, 1, 0.0, phi, 0.0, 0.0, 0.0);
    py::pyrobo(1, 1, 0.0, 0.0, beta(1), beta(2), beta(3));
    if (!allhad) {
      K(N - 1, 1) = 13;
      softrec = false;
    } else {
      // C--boost to fluid rest frame
      ys = 0.5 * std::log((mv(1, 4) + mv(1, 3)) / (mv(1, 4) - mv(1, 3)));
      p3boost = std::sinh(-ys) * P(N - 1, 4) + std::cosh(-ys) * P(N - 1, 3);
      pboost = std::sqrt(ipow<2>(p3boost) + ipow<2>(P(N - 1, 1)) + ipow<2>(P(N - 1, 2)));
      localt = medium().gettemp(mv(1, 1), mv(1, 2), mv(1, 3), mv(1, 4));
      if (pboost < (recsoftcut * 3. * localt)) {
        softrec = true;
        K(N - 1, 1) = 13;
      } else {
        softrec = false;
        if (scatrecoil && (pboost > (rechardcut * 3. * localt))) {
          K(N - 1, 1) = 2;
        } else {
          K(N - 1, 1) = 3;
        }
      }
    }
    if (rejectt) K(N - 1, 1) = 11;
    // C--outgoing projectile
    if (allhad && (!rejectt) && (!softrec)) {
      if (K(N, 2) == 21) {
        if (dir == 1) {
          trip(N) = colmax + 1;
          anti(N) = anti(line);
        } else {
          trip(N) = trip(line);
          anti(N) = colmax + 1;
        }
      } else if (K(N, 2) > 0) {
        trip(N) = colmax + 1;
        anti(N) = 0;
      } else {
        trip(N) = 0;
        anti(N) = colmax + 1;
      }
      colmax = colmax + 1;
    } else {
      trip(N) = trip(line);
      anti(N) = anti(line);
    }
    // C--outgoing scattering centre
    if (allhad && (!rejectt) && (!softrec)) {
      if ((K(N, 2) > 0) && (dir >= 0)) {
        trip(N - 1) = trip(line);
        anti(N - 1) = trip(N);
      } else {
        trip(N - 1) = anti(N);
        anti(N - 1) = anti(line);
      }
    } else {
      trip(N - 1) = 0;
      anti(N - 1) = 0;
    }
    // C--splitting variables for outgoing projectile
    if ((K(N, 1) == 2) && (z == 0.0)) {
      zd(N) = -1.0;
      thetad(N) = -1.0;
    } else {
      zd(N) = z;
      if (z == 0.0) {
        thetad(N) = 0.0;
      } else {
        thetad(N) = (mqk > 0.0) ? hftheta(virt2(N), mqk, z, P(N, 4))
                                : P(N, 5) / (std::sqrt(z * (1.0 - z)) * P(N, 4));
      }
    }
    qqbard(N) = qqbar;
    if (cb::hf::on) cb::hf::qqbarfl(N) = qqbar ? cb::hf::pendingfl : 0;
    aoreq(N) = false;
    // C--splitting variables for outgoing scattering centre
    if (splitrec) {
      if (K(N - 1, 1) == 2) K(N - 1, 1) = 1;
      thetad(N - 1) = P(N - 1, 5) / (std::sqrt(z4 * (1. - z4)) * P(N - 1, 4));
      zd(N - 1) = z4;
      qqbard(N - 1) = qqbardec;
      if (cb::hf::on) cb::hf::qqbarfl(N - 1) = 0;
    } else {
      thetad(N - 1) = -1.0;
      zd(N - 1) = -1.0;
      qqbard(N - 1) = false;
    }
    aoreq(N - 1) = false;
    mv(N, 4) = mv(1, 4);
    mv(N - 1, 4) = mv(1, 4);
    // C--set the production vertices: x_mother + (tprod - tprod_mother) * beta_mother
    mv(N - 1, 1) = mv(line, 1) + (mv(N - 1, 4) - mv(line, 4)) * P(line, 1) /
                                     std::max(py::pyp(line, 8), P(line, 4));
    mv(N - 1, 2) = mv(line, 2) + (mv(N - 1, 4) - mv(line, 4)) * P(line, 2) /
                                     std::max(py::pyp(line, 8), P(line, 4));
    mv(N - 1, 3) = mv(line, 3) + (mv(N - 1, 4) - mv(line, 4)) * P(line, 3) /
                                     std::max(py::pyp(line, 8), P(line, 4));
    mv(N, 1) = mv(line, 1) +
               (mv(N, 4) - mv(line, 4)) * P(line, 1) / std::max(py::pyp(line, 8), P(line, 4));
    mv(N, 2) = mv(line, 2) +
               (mv(N, 4) - mv(line, 4)) * P(line, 2) / std::max(py::pyp(line, 8), P(line, 4));
    mv(N, 3) = mv(line, 3) +
               (mv(N, 4) - mv(line, 4)) * P(line, 3) / std::max(py::pyp(line, 8), P(line, 4));
    if (P(N - 1, 5) > P(1, 5)) {
      lambda = 1.0 / (ftfac * r4(0.2f) * P(N - 1, 4) / ipow<2>(P(N - 1, 5)));
      mv(N - 1, 5) = mv(N - 1, 4) - std::log(1.0 - py::pyr(0)) / lambda;
    } else {
      mv(N - 1, 5) = 0.0;
    }
    if (j < n2) {
      mv(N, 5) = scatcentres(j + 1, 10);
    } else {
      if (offshell(N)) {
        if (deltam == 0.0) {
          if (forcohscat) {
            endtime = time;
          } else {
            endtime = firsttime;
          }
        } else {
          lambda = 1.0 / (ftfac * P(N, 4) * r4(0.2f) / virt2(N));
          endtime = scatcentres(j, 10) - std::log(1.0 - py::pyr(0)) / lambda;
        }
        mv(N, 5) = endtime;
      } else {
        mv(N, 5) = 0.0;
      }
    }
    mv(line, 5) = allqs(j, 6);

    // C--store scattering centre before interaction in separate common block
    if ((writescatcen || dosubtraction) && (!rejectt) && (nscatcen < maxnscatcen)) {
      nscatcen = nscatcen + 1;
      if (nscatcen > maxnscatcen) {
        io::log(" WARNING: no room left to store further scattering centres");
        goto l230;
      }
      // position of the scattering (line 1 holds the scattering centre)
      for (int c = 1; c <= 4; ++c) scatpos(nscatcen, c) = cb::time::mv(1, c);
      // library mode: soft recoils (RECSOFTCUT) leave this slot unfilled; mark it empty
      // instead of keeping the momentum from an earlier shower
      if (cb::lib::active) for (int c = 1; c <= 5; ++c) scatcen(nscatcen, c) = 0.0;
      // library mode: the recoil of this scattering is line N-1 (holes are only kept
      // while their recoil exists, see jetscape_api.cpp)
      if (cb::lib::active) cb::lib::scatrec(nscatcen) = N - 1;
      if (recmode == 0) {
        if (!softrec) {
          scatflav(nscatcen) = K(1, 2);
          scatcen(nscatcen, 1) = P(1, 1);
          scatcen(nscatcen, 2) = P(1, 2);
          scatcen(nscatcen, 3) = P(1, 3);
          scatcen(nscatcen, 4) = P(1, 4);
          scatcen(nscatcen, 5) = P(1, 5);
          dummies(nscatcen, 1) = scatcen(nscatcen, 1) * pdummy / scatcen(nscatcen, 4);
          dummies(nscatcen, 2) = scatcen(nscatcen, 2) * pdummy / scatcen(nscatcen, 4);
          dummies(nscatcen, 3) = scatcen(nscatcen, 3) * pdummy / scatcen(nscatcen, 4);
          dummies(nscatcen, 4) = pdummy;
          dummies(nscatcen, 5) =
              std::sqrt(ipow<2>(dummies(nscatcen, 4)) - ipow<2>(dummies(nscatcen, 1)) -
                        ipow<2>(dummies(nscatcen, 2)) - ipow<2>(dummies(nscatcen, 3)));
          // C--------------------
        }
      } else if (recmode == 1) {
        if (!softrec) {
          scatflav(nscatcen) = K(1, 2);
          scatcen(nscatcen, 1) = P(1, 1);
          scatcen(nscatcen, 2) = P(1, 2);
          scatcen(nscatcen, 3) = P(1, 3);
          scatcen(nscatcen, 4) = P(1, 4);
          scatcen(nscatcen, 5) = P(1, 5);
          dummies(nscatcen, 1) = P(N - 1, 1) * pdummy / P(N - 1, 4);
          dummies(nscatcen, 2) = P(N - 1, 2) * pdummy / P(N - 1, 4);
          dummies(nscatcen, 3) = P(N - 1, 3) * pdummy / P(N - 1, 4);
          dummies(nscatcen, 4) = pdummy;
          dummies(nscatcen, 5) =
              std::sqrt(ipow<2>(dummies(nscatcen, 4)) - ipow<2>(dummies(nscatcen, 1)) -
                        ipow<2>(dummies(nscatcen, 2)) - ipow<2>(dummies(nscatcen, 3)));
        }
      } else if (recmode == 2) {
        scatflav(nscatcen) = K(1, 2);
        scatcen(nscatcen, 1) = P(N - 1, 1) - P(1, 1);
        scatcen(nscatcen, 2) = P(N - 1, 2) - P(1, 2);
        scatcen(nscatcen, 3) = P(N - 1, 3) - P(1, 3);
        scatcen(nscatcen, 4) = P(N - 1, 4) - P(1, 4);
        qmass2 = ipow<2>(scatcen(nscatcen, 4)) - ipow<2>(scatcen(nscatcen, 1)) -
                 ipow<2>(scatcen(nscatcen, 2)) - ipow<2>(scatcen(nscatcen, 3));
        scatcen(nscatcen, 5) = fsign(std::sqrt(std::fabs(qmass2)), qmass2);
        dummies(nscatcen, 1) = scatcen(nscatcen, 1) * pdummy / scatcen(nscatcen, 4);
        dummies(nscatcen, 2) = scatcen(nscatcen, 2) * pdummy / scatcen(nscatcen, 4);
        dummies(nscatcen, 3) = scatcen(nscatcen, 3) * pdummy / scatcen(nscatcen, 4);
        dummies(nscatcen, 4) = pdummy;
        dummies(nscatcen, 5) =
            -std::sqrt(ipow<2>(dummies(nscatcen, 1)) + ipow<2>(dummies(nscatcen, 2)) +
                       ipow<2>(dummies(nscatcen, 3)) - ipow<2>(dummies(nscatcen, 4)));
        if (scatcen(nscatcen, 4) < 0.0) {
          dummies(nscatcen, 1) = -1. * dummies(nscatcen, 1);
          dummies(nscatcen, 2) = -1. * dummies(nscatcen, 2);
          dummies(nscatcen, 3) = -1. * dummies(nscatcen, 3);
        }
      } else if (recmode == 3) {
        if (softrec) {
          scatflav(nscatcen) = K(1, 2);
          scatcen(nscatcen, 1) = P(N - 1, 1) - P(1, 1);
          scatcen(nscatcen, 2) = P(N - 1, 2) - P(1, 2);
          scatcen(nscatcen, 3) = P(N - 1, 3) - P(1, 3);
          scatcen(nscatcen, 4) = P(N - 1, 4) - P(1, 4);
          qmass2 = ipow<2>(scatcen(nscatcen, 4)) - ipow<2>(scatcen(nscatcen, 1)) -
                   ipow<2>(scatcen(nscatcen, 2)) - ipow<2>(scatcen(nscatcen, 3));
          scatcen(nscatcen, 5) = fsign(std::sqrt(std::fabs(qmass2)), qmass2);
          dummies(nscatcen, 1) = scatcen(nscatcen, 1) * pdummy / scatcen(nscatcen, 4);
          dummies(nscatcen, 2) = scatcen(nscatcen, 2) * pdummy / scatcen(nscatcen, 4);
          dummies(nscatcen, 3) = scatcen(nscatcen, 3) * pdummy / scatcen(nscatcen, 4);
          dummies(nscatcen, 4) = pdummy;
          dummies(nscatcen, 5) =
              -std::sqrt(ipow<2>(dummies(nscatcen, 1)) + ipow<2>(dummies(nscatcen, 2)) +
                         ipow<2>(dummies(nscatcen, 3)) - ipow<2>(dummies(nscatcen, 4)));
          if (scatcen(nscatcen, 4) < 0.0) {
            dummies(nscatcen, 1) = -1. * dummies(nscatcen, 1);
            dummies(nscatcen, 2) = -1. * dummies(nscatcen, 2);
            dummies(nscatcen, 3) = -1. * dummies(nscatcen, 3);
          }
        } else {
          scatflav(nscatcen) = K(1, 2);
          scatcen(nscatcen, 1) = P(1, 1);
          scatcen(nscatcen, 2) = P(1, 2);
          scatcen(nscatcen, 3) = P(1, 3);
          scatcen(nscatcen, 4) = P(1, 4);
          scatcen(nscatcen, 5) = P(1, 5);
          dummies(nscatcen, 1) = scatcen(nscatcen, 1) * pdummy / scatcen(nscatcen, 4);
          dummies(nscatcen, 2) = scatcen(nscatcen, 2) * pdummy / scatcen(nscatcen, 4);
          dummies(nscatcen, 3) = scatcen(nscatcen, 3) * pdummy / scatcen(nscatcen, 4);
          dummies(nscatcen, 4) = pdummy;
          dummies(nscatcen, 5) =
              std::sqrt(ipow<2>(dummies(nscatcen, 4)) - ipow<2>(dummies(nscatcen, 1)) -
                        ipow<2>(dummies(nscatcen, 2)) - ipow<2>(dummies(nscatcen, 3)));
          // C--------------------
        }
      }
    }
  l230:
    if (!isrecoil) hfstat_kin(K(line, 2), rejectt ? 0 : 1, P(line, 4) - P(N, 4));

    dmleft = dmleft - (P(N, 5) - P(line, 5));
    line = N;
    tmp = std::fabs(ipow<2>(P(N, 4)) - ipow<2>(P(N, 1)) - ipow<2>(P(N, 2)) - ipow<2>(P(N, 3))) -
          ipow<2>(P(N, 5));
    if (std::fabs(tmp) >= 1.e-6)
      io::log(" dokinematics 4-momentum test failed: {} {} {} {} {} {}", tmp, j, P(l, 5),
              P(line, 5), P(N, 5), fbool(reshuffle));
  }  // 222 CONTINUE
  // C--set decay time for last outgoing parton
  if (N > nold) {
    if (offshell(N)) {
      if (deltam == 0.0) {
        if (forcohscat) {
          endtime = time;
        } else {
          endtime = firsttime;
        }
      } else {
        lambda = 1.0 / (ftfac * P(N, 4) * r4(0.2f) / virt2(N));
        endtime = mv(N, 4) - std::log(1.0 - py::pyr(0)) / lambda;
      }
      mv(N, 5) = endtime;
    } else {
      mv(N, 5) = 0.0;
    }
    // NOTE(port): SCATCENTRES(J,10) is read with J = N2+1 after the loop (as in Fortran; only logged).
    if ((!rejectt) && (mv(N, 5) < mv(N, 4) - 1e-6) && (mv(N, 5) > 1e-6))
      io::log(" unhealthy:  {} {} {} {} {} {} {} {} {} {} {}", N, mv(N, 4), scatcentres(j, 10),
              endtime, time, n2, l, P(N, 5), deltam, firsttime, fbool(forcohscat));
  }

  if (P(lpart, 5) < 0.0) {
    retrysplit = true;
    return;
  }
  if (std::fabs(P(lpart, 5) / newm2 - 1.0) > 1.e-6) {
    // !        write(logfid,*)'dokinematics did not reach target mass',
    // !     &	p(lpart,5),newm2,newm,rejectt
    // C--if reached mass is reasonable but lower than target, then accept it nevertheless
    if ((x == 1.0) && (vir(P(lpart, 5)) > q0)) goto l302;
    if (x != 1.0) retrysplit = true;
    redokin = true;
    ++nretry;
    if (hfdbg && (nretry % 1000 == 0))
      io::log(" HFDEBUG dokinematics retry {}: kf={} K1={} P5(l)={} reached={} newm={} newm2={} x={} n1={} n2={} E={}", nretry,
              K(l, 2), K(l, 1), P(l, 5), P(lpart, 5), newm, newm2, x, n1, n2, P(l, 4));
    N = nold;
    colmax = colmaxold;
    nscatcen = nscatcenold;
    K(l, 1) = statold;
    // heavy quark: a target mass that is kinematically out of reach (e.g. bringing a spacelike heavy quark
    // back to M = m with too little momentum transfer) would be retried forever; give up after 10 attempts
    // and let the caller redo the interaction without the splitting (retrysplit), as JEWEL does elsewhere
    if ((mqk > 0.0) && (nretry >= 10)) {
      retrysplit = true;
      return;
    }
    if (P(l, 5) < 0.0) {
      newm2 = tgt(0.0);
    } else {
      if ((vir(P(l, 5)) < q0) && (K(l, 1) != 14)) {
        if (newm2 == tgt(0.0)) {
          retrysplit = true;
          return;
        } else if ((newm2 == tgt(newm)) && (newm != q0 + 1.e-6)) {
          newm2 = tgt(q0 + 1.e-6);
        } else {
          newm2 = tgt(0.0);
        }
      } else {
        newm2 = P(l, 5);
        if (K(l, 1) == 14) z = 0.0;
      }
      for (int i = n1, i_end = n2; i <= i_end; ++i) {
        if (allqs(i, 1) < 0.0) {
          n2 = i;
          break;  // goto 301
        }
      }
      // 301 continue
    }
    goto l204;
  }
l302:;
}

}  // namespace jewel
