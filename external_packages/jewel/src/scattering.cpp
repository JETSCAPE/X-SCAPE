// JEWEL 2.6.0 C++ port -- initial-state scattering (DOINSTATESCAT) and coherent
// scattering of colour dipoles (docohscattering). Translated from jewel-2.6.0.f:4132-5079.
#include <algorithm>
#include <cmath>
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

// ***********************************************************************
// ***	  subroutine doinstatescat
// ***********************************************************************
void doinstatescat(int& l, double& x, int& typi, double& q, const double& tstart, double& deltat,
                   bool& overq0, const bool& retrysplit) {
  // C--Common block of Pythia
  using namespace py::pyjets;
  // C--time common block
  using namespace cb::time;
  // C--factor in front of formation times
  using namespace cb::ftimefac;
  // C--Parameter common block
  using namespace cb::param;
  // C--coherence parameters
  using namespace cb::cohparams;
  // C--discard event flag
  using namespace cb::disc;
  // C--variables for coherent scattering
  using namespace cb::coherent;
  // C--information needed to go back to before coherent scattering
  using namespace cb::revertcohscat;
  // C--variables for angular ordering
  using namespace cb::angor;
  // C--identifier of file for hepmc output and logfile
  using namespace cb::hepmcid;

  // C--local variables
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/
  constexpr int countmax = 10000;           // DATA COUNTMAX/10000/
  int counter, count2;
  double deltal, pnorad, weight, low, fmax, sigmatot, pfchange, tnow, tleft, xmax, tsum,
      q2old, tmax, dt, xsc, ysc, zsc, tsc, ms1, md1, lambda, rtau, phi, tauest, zdum, dtfirst;
  [[maybe_unused]] double tmaxnew, qold;
  [[maybe_unused]] double nrejl, naccl;
  F1<double, 4> qsumvecold;
  bool fchange, norad, noscat, qqbardum;
  char typ;
  [[maybe_unused]] std::string typ2;
  double splitangle, splittime, dperp, qperp, qdotp, q2;
  double mql = 0.0, cvl = 0.0;  // heavy flavour: quark mass and current virtuality of line l

  dtfirst = -1.0;
  tsum = 0.0;
  naccl = 0.0;
  nrejl = 0.0;
  didcohscat = false;
  // C--go back to last splitting of parton in line L
  // NOTE(port): if K(l,3)==0 the Fortran reads K(0,1) (out of bounds, i.e. the word before K in
  // PYJETS = NPAD); the FArray access reproduces the same address. Kept as in the original.
  if (K(K(l, 3), 1) == 11) {
    splitangle = thetad(K(l, 3));
    splittime = mv(K(l, 3), 5);
  } else {
    splitangle = -1.0;
    splittime = 0.0;
  }

  counter = 0;

  xsc = mv(l, 1) + (tstart - mv(l, 4)) * P(l, 1) / P(l, 4);
  ysc = mv(l, 2) + (tstart - mv(l, 4)) * P(l, 2) / P(l, 4);
  zsc = mv(l, 3) + (tstart - mv(l, 4)) * P(l, 3) / P(l, 4);
  tsc = tstart;
  md1 = medium().getmd(xsc, ysc, zsc, tsc);
  ms1 = medium().getms(xsc, ysc, zsc, tsc);
  tleft = deltat;
  tnow = tstart;

  if (md1 <= 1.e-4 || ms1 <= 1.e-4) {
    const double tempv = medium().gettemp(xsc, ysc, zsc, tsc);
    const double neffv = medium().getneff(xsc, ysc, zsc, tsc);
    io::log(" problem! {} {}", tempv, neffv);
  }

  // C--check for scattering
l120:
  noscat = !getdeltat(l, tnow, tleft, dt, false);
  if (dtfirst < 0.0) dtfirst = dt;
  if (noscat && (!retrysplit)) {
    // C--allow for coherent scattering of colour dipole before first resolved scattering
    nstart = 0;
    nend = 0;
    docohscattering(l, tstart, deltat, splitangle, splittime);
    goto l116;
  }

  // C--decide whether there will be radiation
  pnorad = getpnorad1(l, xsc, ysc, zsc, tsc);
  if ((py::pyr(0) < pnorad) || (P(l, 4) < r4(1.001f) * q0)) {
    norad = true;
  } else {
    norad = false;
  }

  // C--decide whether q or g is to be scattered
  xsc = mv(l, 1) + (tstart + tsum + dt - mv(l, 4)) * P(l, 1) / P(l, 4);
  ysc = mv(l, 2) + (tstart + tsum + dt - mv(l, 4)) * P(l, 2) / P(l, 4);
  zsc = mv(l, 3) + (tstart + tsum + dt - mv(l, 4)) * P(l, 3) / P(l, 4);
  tsc = tstart + tsum + dt;
  if (K(l, 2) == 21) {
    typ = 'G';
    typ2 = "GC";
    sigmatot = getsscat(P(l, 4), P(l, 1), P(l, 2), P(l, 3), P(l, 5), q0, 'G', 'C', xsc, ysc, zsc,
                        tsc, 0);
    if ((sigmatot == 0.0) || (pnorad == 1.0)) {
      pfchange = 0.0;
    } else {
      pfchange = getsscat(P(l, 4), P(l, 1), P(l, 2), P(l, 3), P(l, 5), q0, 'G', 'Q', xsc, ysc,
                          zsc, tsc, 0) /
                 sigmatot;
    }
    sigmatot = getsscat(P(l, 4), P(l, 1), P(l, 2), P(l, 3), P(l, 5), 0.0, 'G', 'C', xsc, ysc, zsc,
                        tsc, 0);
  } else {
    typ = 'Q';
    typ2 = "QQ";
    sigmatot = getsscat(P(l, 4), P(l, 1), P(l, 2), P(l, 3), P(l, 5), q0, 'Q', 'C', xsc, ysc, zsc,
                        tsc, 0);
    if ((sigmatot == 0.0) || (pnorad == 1.0)) {
      pfchange = 0.0;
    } else {
      pfchange = getsscat(P(l, 4), P(l, 1), P(l, 2), P(l, 3), P(l, 5), q0, 'Q', 'G', xsc, ysc,
                          zsc, tsc, 0) /
                 sigmatot;
    }
    sigmatot = getsscat(P(l, 4), P(l, 1), P(l, 2), P(l, 3), P(l, 5), 0.0, 'Q', 'C', xsc, ysc, zsc,
                        tsc, 0);
  }
  if ((pfchange < -1.e-4) || (pfchange > 1.0 + 1.e-4)) {
    io::log(" error: flavour change probability= {} for {}", pfchange, typ);
  }
  if (py::pyr(0) < pfchange) {
    fchange = true;
  } else {
    fchange = false;
  }
  if (norad) fchange = false;
  // heavy quarks keep their identity (no Q -> g flavour change; the random number above is still drawn)
  mql = hfmass(K(l, 2));
  if (mql > 0.0) fchange = false;
  cvl = virtuality(P(l, 5), mql);  // current virtuality (= P(l,5) for light partons)
  // C--set TYPI
  if (typ == 'G') {
    if (fchange) {
      typi = fint(fsign(2.0, py::pyr(0) - 0.5));
    } else {
      typi = K(l, 2);
    }
  } else {
    if (fchange) {
      typi = 21;
    } else {
      typi = K(l, 2);
    }
  }
  low = ipow<2>(q0) / ipow<2>(scalefacm);
  tmax = 4.0 * (ipow<2>(P(l, 4)) - ipow<2>(P(l, 5)));
  xmax = 1.0 - ipow<2>(q0) / (ipow<2>(scalefacm) * 4.0 * tmax);

  if (sigmatot == 0.0) goto l116;

  rtau = py::pyr(0);

  // C--generate a trial emission
  // C--pick a x value from splitting function
l112:
  counter = counter + 1;
  if (typ == 'G') {
    if (fchange) {
      x = generatez(0.0, 0.0, 1.0 - xmax, "QG");
    } else {
      x = generatez(0.0, 0.0, 1.0 - xmax, "GG");
    }
  } else {
    if (fchange) {
      x = 1.0 - generatez(0.0, 0.0, 1.0 - xmax, "QQ");
    } else {
      x = generatez(0.0, 0.0, 1.0 - xmax, "QQ");
    }
  }
  if (norad) x = 1.0;
  // C--initialisation
  tmaxnew = ipow<2>(x * P(l, 4));
  phi = 0.0;
  qsumvec(1) = 0.0;
  qsumvec(2) = 0.0;
  qsumvec(3) = 0.0;
  qsumvec(4) = 0.0;
  qsum2 = -1.e-10;
  overq0 = false;
  q = cvl;
  qold = cvl;
  tauest = deltat;
  // C--generate first momentum transfer
  deltal = dt;
  nstart = 1;
  nend = 1;
  tnow = tnow + deltal;
  tsum = tsum + deltal;
  tleft = tleft - deltal;
  allqs(nend, 6) = tnow;
  q2old = qsum2;
  // C--get new momentum transfer
  count2 = 0;
l118:
  getqvec(l, nend, tnow - mv(l, 4), x, false);
  if (-qsum2 > ipow<2>(P(l, 4))) {
    qsumvec(1) = 0.0;
    qsumvec(2) = 0.0;
    qsumvec(3) = 0.0;
    qsumvec(4) = 0.0;
    qsum2 = q2old;
    if (count2 < 100) {
      count2 = count2 + 1;
      goto l118;
    } else {
      allqs(nend, 1) = 0.0;
      allqs(nend, 2) = 0.0;
      allqs(nend, 3) = 0.0;
      allqs(nend, 4) = 0.0;
      allqs(nend, 5) = 0.0;
      nend = 0;
    }
  }
  // C--if we have a scattering, check whether momentum transfer is above resolution scale
  if ((nend == 1) && (splitangle > 0.0)) {
    qdotp = allqs(1, 1) * P(l, 1) + allqs(1, 2) * P(l, 2) + allqs(1, 3) * P(l, 3);
    q2 = ipow<2>(allqs(1, 1)) + ipow<2>(allqs(1, 2)) + ipow<2>(allqs(1, 3));
    qperp = std::sqrt(q2 - ipow<2>(qdotp) / py::pyp(l, 7));
    dperp = cohlengthfac * splitangle * (tnow - splittime);
    // C--factor 0.2 to convert 1/fm to GeV
    if ((qperp < (r4(0.2f) / dperp)) && (py::pyr(0) < pcohrej)) {
      qsumvec(1) = 0.0;
      qsumvec(2) = 0.0;
      qsumvec(3) = 0.0;
      qsumvec(4) = 0.0;
      qsum2 = 0.0;
      allqs(1, 1) = 0.0;
      allqs(1, 2) = 0.0;
      allqs(1, 3) = 0.0;
      allqs(1, 4) = 0.0;
      allqs(1, 5) = 0.0;
      nstart = 0;
      nend = 0;
      goto l120;
    }
  }
  // C--update OVERQ0
  if (nend > 0) {
    if (-allqs(nend, 1) > low) overq0 = true;
  }
  // C--get new virtuality
  if (overq0 && (!norad)) {
    q = getnewmass(l, ipow<2>(scalefacm) * qsum2, ipow<2>(scalefacm) * q2old, 0.0, true, x, zdum,
                   qqbardum);
  } else {
    q = 0.0;
  }

  // C--allow for coherent scattering of colour dipole before first resolved scattering
  docohscattering(l, tstart, tsum, splitangle, splittime);
  if (didcohscat) getqvec(l, nend, tnow - mv(l, 4), x, true);

  // C--estimate formation time
l111:
  if ((q == 0.0) || (q == cvl)) {
    tauest = deltat;
  } else {
    tauest = std::max(ftfac * (1.0 - phi) * r4(0.2f) * x * P(l, 4) / ipow<2>(q), 0.0);
  }
  if (tauest > 0.0) {
    lambda = 1.0 / tauest;
    tauest = -std::log(1.0 - rtau) / lambda;
  }

  // C--find number, position and momentum transfers of further scatterings
  {
    const double dtmaxarg = std::min(tleft, tauest);
    noscat = !getdeltat(l, tnow, dtmaxarg, deltal, false);
  }
  noscat = true;
  if ((!noscat) && (!norad)) {
    // C--add a momentum transfer
    nend = nend + 1;
    if (nend >= 100) {
      nend = nend - 1;
      goto l114;
    }
    tnow = tnow + deltal;
    tsum = tsum + deltal;
    tleft = tleft - deltal;
    // C--update phase
    if ((q != 0.0) && (q != cvl)) {
      phi = phi + 5.0 * deltal * ipow<2>(q) / (1.0 * x * P(l, 4));
    }
    // C--get new momentum transfer
    allqs(nend, 6) = tnow;
    q2old = qsum2;
    qsumvecold(1) = qsumvec(1);
    qsumvecold(2) = qsumvec(2);
    qsumvecold(3) = qsumvec(3);
    qsumvecold(4) = qsumvec(4);
    count2 = 0;
  l119:
    getqvec(l, nend, tnow - mv(l, 4), x, false);
    if (-qsum2 > ipow<2>(P(l, 4))) {
      qsumvec(1) = qsumvecold(1);
      qsumvec(2) = qsumvecold(2);
      qsumvec(3) = qsumvecold(3);
      qsumvec(4) = qsumvecold(4);
      qsum2 = q2old;
      if (count2 < 100) {
        count2 = count2 + 1;
        goto l119;
      } else {
        allqs(nend, 1) = 0.0;
        allqs(nend, 2) = 0.0;
        allqs(nend, 3) = 0.0;
        allqs(nend, 4) = 0.0;
        allqs(nend, 5) = 0.0;
        nend = nend - 1;
      }
    }
    // C--update OVERQ0
    if ((nend > 0) && ((-qsum2 > low) || (-allqs(nend, 1) > low))) overq0 = true;
    // C--get new virtuality
    qold = q;
    if (overq0 && (!norad)) {
      q = getnewmass(l, ipow<2>(scalefacm) * qsum2, ipow<2>(scalefacm) * q2old, 0.0, true, x,
                     zdum, qqbardum);
    } else {
      q = 0.0;
    }
    goto l111;
  }

  // C--do reweighting
l114:
  tmaxnew = ipow<2>(x) * ipow<2>(P(l, 4));
  if (norad) {
    weight = 1.0;
    q = 0.0;
    x = 1.0;
  } else if ((-qsum2 < low) || (q == 0.0)) {
    weight = 0.0;
  } else if (-qsum2 > ipow<2>(P(l, 4))) {
    weight = 0.0;
  } else {
    if (typ == 'G') {
      fmax = 2.0 * std::log(-ipow<2>(scalefacm) * qsum2 / ipow<2>(q0)) *
             alphas(ipow<2>(q0) / 4.0, lps) / (2.0 * pi);
      if (qsum2 == 0.0) {
        weight = 0.0;
        norad = true;
      } else {
        if (fchange) {
          weight = 2.0 * getpdf(x, scalefacm * std::sqrt(-qsum2), "QG") / (pqg(x) * fmax);
          if ((weight > 1.0 + 1.e-4) || (weight < -1.e-4)) {
            const double sq = std::sqrt(-qsum2);
            const double pdfv = getpdf(x, scalefacm * std::sqrt(-qsum2), "QG");
            io::log(" x,sqrt(qsum^2),getpdf,fmax: {} {} {} qg {}", x, sq, pdfv, fmax);
          }
        } else {
          weight = getpdf(x, scalefacm * std::sqrt(-qsum2), "GG") / (pgg(x) * fmax);
          if ((weight > 1.0 + 1.e-4) || (weight < -1.e-4)) {
            const double sq = std::sqrt(-qsum2);
            const double pdfv = getpdf(x, scalefacm * std::sqrt(-qsum2), "GG");
            io::log(" x,sqrt(qsum^2),getpdf,fmax: {} {} {} gg {}", x, sq, pdfv, fmax);
          }
        }
      }
    } else {
      fmax = std::log(-ipow<2>(scalefacm) * qsum2 / ipow<2>(q0)) *
             alphas(ipow<2>(q0) / 4.0, lps) / (2.0 * pi);
      if (qsum2 == 0.0) {
        weight = 0.0;
        norad = true;
      } else {
        if (fchange) {
          weight = getpdf(x, scalefacm * std::sqrt(-qsum2), "GQ") / (pgq(x) * fmax);
          if ((weight > 1.0 + 1.e-4) || (weight < -1.e-4)) {
            const double sq = std::sqrt(-qsum2);
            const double pdfv = getpdf(x, scalefacm * std::sqrt(-qsum2), "GQ");
            io::log(" x,sqrt(qsum^2),getpdf:,fmax {} {} {} gq {}", x, sq, pdfv, fmax);
          }
        } else {
          weight = getpdf(x, scalefacm * std::sqrt(-qsum2), "QQ") / (pqq(x) * fmax);
          if ((weight > 1.0 + 1.e-4) || (weight < -1.e-4)) {
            const double sq = std::sqrt(-qsum2);
            const double pdfv = getpdf(x, scalefacm * std::sqrt(-qsum2), "QQ");
            io::log(" x,sqrt(qsum^2),getpdf,fmax: {} {} {} qq {}", x, sq, pdfv, fmax);
          }
        }
      }
    }
  }
  if ((weight > 1.0 + 1.e-4) || (weight < -1.e-4)) io::log(" error: weight= {}", weight);
  // 115
  if (py::pyr(0) > weight) {
    if (counter < countmax) {
      // C--undo coherent scattering, go back and retry
      if (didcohscat) {
        N = nbef;
        // NOTE(port): JEWEL 2.6.0 keeps the scattering centres stored by the undone coherent
        // scattering (with DOSUBTRACTION they are subtracted although their recoils are gone).
        // Kept in the stand-alone program; as a library they would be holes without recoils.
        if (cb::lib::active) cb::storescatcen::nscatcen = cb::lib::nscatbef;
        l = lbef;
        K(l, 1) = klbef1;
        K(lsib, 1) = klsib1;
        K(l, 4) = 0;
        K(l, 5) = 0;
        K(lsib, 4) = 0;
        K(lsib, 5) = 0;
        mv(l, 5) = decaytl;
        mv(lsib, 5) = decaytsib;
        didcohscat = false;
      }
      tleft = deltat;
      tnow = tstart;
      tsum = 0.0;
      dt = dtfirst;
      goto l112;
    } else {
      q = 0.0;
      x = 1.0;
      nend = nstart;
      qsum2 = allqs(nend, 1);
      qsumvec(1) = allqs(nend, 2);
      qsumvec(2) = allqs(nend, 3);
      qsumvec(3) = allqs(nend, 4);
      qsumvec(4) = allqs(nend, 5);
      typi = K(l, 2);
      if (-allqs(nend, 1) > low) {
        overq0 = true;
      } else {
        overq0 = false;
      }
      deltat = allqs(nend, 6) - tstart;
      tnow = allqs(1, 6);
      return;
    }
  }
  // heavy quark: thin the accepted initial-state emission with the massive/massless kernel ratio (dead cone);
  // a thinned emission becomes an elastic scattering, as in the fallback above
  if ((mql > 0.0) && (!norad) && (x < 1.0) &&
      (py::pyr(0) > hfisweight(x, -ipow<2>(scalefacm) * qsum2, mql))) {
    q = 0.0;
    x = 1.0;
    nend = nstart;
    qsum2 = allqs(nend, 1);
    qsumvec(1) = allqs(nend, 2);
    qsumvec(2) = allqs(nend, 3);
    qsumvec(3) = allqs(nend, 4);
    qsumvec(4) = allqs(nend, 5);
    typi = K(l, 2);
    overq0 = (-allqs(nend, 1) > low);
    deltat = allqs(nend, 6) - tstart;
    tnow = allqs(1, 6);
    return;
  }
  // C--found meaningful configuration, now do final checks
  // C--check if phase is unity and weight with 1/Nscat
  if (nend == 0) {
    q = 0.0;
    x = 1.0;
    typi = K(l, 2);
    nstart = 0;
    nend = 0;
    qsum2 = 0.0;
    qsumvec(1) = 0.0;
    qsumvec(2) = 0.0;
    qsumvec(3) = 0.0;
    qsumvec(4) = 0.0;
    overq0 = false;
    return;
  }
  if (((tleft < tauest) || (py::pyr(0) > 1.0 / (nend * 1.0))) && (!norad)) {
    q = 0.0;
    x = 1.0;
    nend = nstart;
    qsum2 = allqs(nend, 1);
    qsumvec(1) = allqs(nend, 2);
    qsumvec(2) = allqs(nend, 3);
    qsumvec(3) = allqs(nend, 4);
    qsumvec(4) = allqs(nend, 5);
    typi = K(l, 2);
    if (-allqs(nend, 1) > low) {
      overq0 = true;
    } else {
      overq0 = false;
    }
    deltat = allqs(nend, 6) - tstart;
    tnow = allqs(1, 6);
  } else {
    if (!norad) {
      tleft = tleft - tauest;
      tnow = tnow + tauest;
      tsum = tsum + tauest;
    }
    deltat = tsum;
  }
  return;
  // C--exit in case of failure
l116:
  q = 0.0;
  x = 1.0;
  nstart = 0;
  nend = 0;
  qsumvec(1) = 0.0;
  qsumvec(2) = 0.0;
  qsumvec(3) = 0.0;
  qsumvec(4) = 0.0;
  qsum2 = 0.0;
  overq0 = false;
  typi = K(l, 2);
  return;
}

// ***********************************************************************
// ***	  subroutine docohscattering
// ***********************************************************************
void docohscattering(int& l, const double& tstart, const double& tsum, const double& splitangle,
                     const double& splittime) {
  // C--identifier of file for hepmc output and logfile
  using namespace cb::hepmcid;
  // C--Common block of Pythia
  using namespace py::pyjets;
  // C--time common block
  using namespace cb::time;
  // C--keep track of times until which partons have been evolved
  using namespace cb::evoltimes;
  // C--colour index common block
  using namespace cb::colour;
  // C--Parameter common block
  using namespace cb::param;
  // C--variables for angular ordering
  using namespace cb::angor;
  // C--coherence parameters
  using namespace cb::cohparams;
  // C--variables for coherent scattering
  using namespace cb::coherent;
  // C--information needed to go back to before coherent scattering
  using namespace cb::revertcohscat;
  // C--event weight
  using namespace cb::weight;
  // C--histograms for internal analysis
  using namespace cb::histograms;

  // C--local variables
  int mother, sib, count2, sib1, sib2, dir, nold;
  [[maybe_unused]] int nscatl;
  [[maybe_unused]] int lcohorig;
  int bint, binthat;
  double deltat, dt, z, phisplit, phim, thetam, tnow, tleft, qdotp, q2, qperp, dperp, pz, pts, mb,
      mc, p1, p2, decayt1, decayt2;
  F1<double, 4> pos1, pos2, mom1, mom2;
  bool scat, quark, qqbar;

  count2 = 0;
  nscatl = 0;
  deltat = tsum;

  didcohscat = false;

  if (!cohscat) return;

  tnow = tstart;
  tleft = deltat;
  // C--first figure out whether sibling has scattered or had a chance to scatter
  mother = K(l, 3);
  // NOTE(port): mother==0 would read K(0,3) (= K(23000,2) in column-major storage); kept as original.
  if (K(mother, 3) == 0) return;
  // NOTE(port): Fortran `call flush()` (flush all units); only the JEWEL log is a Fortran unit here.
  io::log_flush();
  if (K(mother, 4) == l) {
    sib = K(mother, 5);
  } else {
    sib = K(mother, 4);
  }
  if (sib == 0) return;
  sib1 = std::min(sib, l);
  sib2 = std::max(sib, l);
  if (K(mother, 2) == 21) {
    quark = false;
    if (K(sib1, 2) == 21) {
      qqbar = false;
    } else {
      qqbar = true;
    }
  } else {
    quark = true;
    qqbar = false;
  }
  if (quark || qqbar) {
    dir = 0;
  } else {
    if (trip(sib1) == trip(mother)) {
      dir = 1;
    } else {
      dir = -1;
    }
  }
  // C--if mother did not decay by splitting or sibling is not active any more coherent scattering was already generated
  if ((K(mother, 1) != 11) || (K(sib, 1) != 1) || (mv(sib, 5) < tstart) || (splitangle < 0.0)) {
    return;
  }
  // C--if sibling will decay soon reduce available time for coherent scattering accordingly
  if ((mv(sib, 5) > 0.0) && (mv(sib, 5) < tstart + deltat)) deltat = mv(sib, 5) - tstart;
  if (deltat < 0.0) {
    io::log(" ERROR in docohscattering: deltat < 0 {} {} {} {} {}", deltat, l, mv(sib, 5), tstart,
            tstart + deltat);
    return;
  }
  // C--store information needed to undo coherent scattering
  didcohscat = true;
  nbef = N;
  cb::lib::nscatbef = cb::storescatcen::nscatcen;
  lbef = l;
  lsib = sib;
  klbef1 = K(l, 1);
  klsib1 = K(sib, 1);
  decaytl = mv(l, 5);
  decaytsib = mv(sib, 5);
  decayt1 = mv(sib1, 5);
  decayt2 = mv(sib2, 5);
  // C--there can be coherent scattering -> undo the splitting and save z and phi
  for (int i = 1; i <= 4; ++i) {
    mom1(i) = P(sib1, i);
    mom2(i) = P(sib2, i);
    pos1(i) = mv(sib1, i);
    pos2(i) = mv(sib2, i);
  }
  z = P(sib2, 4) / P(mother, 4);
  // C--rotate such that mother points in z-direction
  thetam = py::pyp(mother, 13);
  phim = py::pyp(mother, 15);
  py::pyrobo(sib2, sib2, 0.0, -phim, 0.0, 0.0, 0.0);
  py::pyrobo(sib2, sib2, -thetam, 0.0, 0.0, 0.0, 0.0);
  phisplit = py::pyp(sib2, 15);
  py::pyrobo(sib2, sib2, thetam, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(sib2, sib2, 0.0, phim, 0.0, 0.0, 0.0);
  mv(sib1, 5) = mv(sib1, 4);
  mv(sib2, 5) = mv(sib2, 4);
  N = N + 1;
  P(N, 1) = P(mother, 1);
  P(N, 2) = P(mother, 2);
  P(N, 3) = P(mother, 3);
  P(N, 4) = P(mother, 4);
  P(N, 5) = P(mother, 5);
  mv(N, 1) = mv(l, 1);
  mv(N, 2) = mv(l, 2);
  mv(N, 3) = mv(l, 3);
  mv(N, 4) = mv(l, 4);
  mv(N, 5) = 0.0;
  trip(N) = trip(mother);
  anti(N) = anti(mother);
  K(N, 1) = 1;
  K(N, 2) = K(mother, 2);
  K(N, 3) = sib1;
  K(l, 1) = 20;
  K(sib, 1) = 20;
  K(l, 4) = N;
  K(l, 5) = N;
  K(sib, 4) = N;
  K(sib, 5) = N;
  thetad(N) = thetad(mother);
  zd(N) = zd(mother);
  qqbard(N) = qqbard(mother);
  if (cb::hf::on) cb::hf::qqbarfl(N) = cb::hf::qqbarfl(mother);
  trip(N) = trip(mother);
  anti(N) = anti(mother);
  lcohorig = N;
  // C--check whether coherent state (=mother) wants to scatter
l120:
  scat = getdeltat(N, tnow, tleft, dt, true);
  if (scat) {
    tnow = tnow + dt;
    tleft = tleft - dt;
    nold = nend;
    nend = nend + 1;
  l118:
    getqvec(N, nend, tnow - mv(N, 4), 1.0, false);
    if (-allqs(nend, 1) > ipow<2>(P(l, 4))) {
      if (count2 < 100) {
        count2 = count2 + 1;
        goto l118;
      } else {
        allqs(nend, 1) = 0.0;
        allqs(nend, 2) = 0.0;
        allqs(nend, 3) = 0.0;
        allqs(nend, 4) = 0.0;
        allqs(nend, 5) = 0.0;
        nend = nend - 1;
      }
    }
    // C--if we have a scattering, check whether momentum transfer is above resolution scale
    if (nend > nold) {
      qdotp = allqs(nend, 1) * P(l, 1) + allqs(nend, 2) * P(l, 2) + allqs(nend, 3) * P(l, 3);
      q2 = ipow<2>(allqs(nend, 1)) + ipow<2>(allqs(nend, 2)) + ipow<2>(allqs(nend, 3));
      qperp = std::sqrt(q2 - ipow<2>(qdotp) / py::pyp(l, 7));
      dperp = cohlengthfac * splitangle * (tnow - splittime);
      // C--factor 0.2 to convert 1/fm to GeV
      if (qperp > (r4(0.2f) / dperp)) {
        allqs(nend, 1) = 0.0;
        allqs(nend, 2) = 0.0;
        allqs(nend, 3) = 0.0;
        allqs(nend, 4) = 0.0;
        allqs(nend, 5) = 0.0;
        nend = nend - 1;
        goto l120;
      } else {
        if ((scatcentres(nend, 10) > tmin) && (scatcentres(nend, 10) < timax)) {
          bint = fint((scatcentres(nend, 10) - tmin) * ntbinsaccrej / (timax - tmin)) + 1;
          histotcoh(bint) = histotcoh(bint) + evweight;
        }
        if ((-allqs(nend, 1) > thatmin) && (-allqs(nend, 1) < thatmax)) {
          binthat = fint((-allqs(nend, 1) - thatmin) * nthatbinsaccrej / (thatmax - thatmin)) + 1;
          histothatcoh(binthat) = histothatcoh(binthat) + evweight;
        }
      }
    }
    // C--now do the kinemtatics for the coherent scatterings
    if (nend > nold) {
      nscatl = nscatl + 1;
      {
        bool retrysplit_tmp = false;
        // target: the current virtuality (heavy quarks: sqrt(M^2 - m^2)) and the line's own g -> q qbar flavour
        cb::hf::pendingfl = cb::hf::qqbarfl(N);
        dokinematics(N, N, nend, nend, virtuality(P(N, 5), hfmass(K(N, 2))), retrysplit_tmp, tnow, 1.0, zd(N),
                     qqbard(N), true);
      }
      // C--update positions and momenta of incoherent state
      // C--make temporary lines in event record for calculation of momenta
      // C--rotate such that coherent state points in z-direction
      thetam = py::pyp(N, 13);
      phim = py::pyp(N, 15);
      py::pyrobo(N, N, 0.0, -phim, 0.0, 0.0, 0.0);
      py::pyrobo(N, N, -thetam, 0.0, 0.0, 0.0, 0.0);
      mb = P(sib2, 5);
      mc = P(sib1, 5);
      pz = (2.0 * z * ipow<2>(P(N, 4)) - ipow<2>(P(N, 5)) - ipow<2>(mb) + ipow<2>(mc)) /
           (2.0 * P(N, 3));
      pts = ipow<2>(z) * (ipow<2>(P(N, 4))) - ipow<2>(pz) - ipow<2>(mb);
      if (pts < 0.0) {
        // C--kinematics doesn't work out -> have to undo last scattering
        N = N - 2;
        K(N, 1) = 1;
        K(N, 4) = 0;
        K(N, 5) = 0;
        mv(N, 5) = 0.0;
        goto l502;
      }
      N = N + 2;
      P(N - 1, 1) = -std::sqrt(pts) * gf_cos(phisplit);
      P(N - 1, 2) = -std::sqrt(pts) * gf_sin(phisplit);
      P(N - 1, 3) = P(N - 2, 3) - pz;
      P(N - 1, 4) = (1.0 - z) * P(N - 2, 4);
      P(N - 1, 5) = mc;
      P(N, 1) = std::sqrt(pts) * gf_cos(phisplit);
      P(N, 2) = std::sqrt(pts) * gf_sin(phisplit);
      P(N, 3) = pz;
      P(N, 4) = z * P(N - 2, 4);
      P(N, 5) = mb;
      K(N - 1, 1) = 1;
      K(N, 1) = 1;
      py::pyrobo(N - 2, N, thetam, 0.0, 0.0, 0.0, 0.0);
      py::pyrobo(N - 2, N, 0.0, phim, 0.0, 0.0, 0.0);
      p1 = std::sqrt(ipow<2>(mom1(1)) + ipow<2>(mom1(2)) + ipow<2>(mom1(3)));
      p2 = std::sqrt(ipow<2>(mom2(1)) + ipow<2>(mom2(2)) + ipow<2>(mom2(3)));
      pos1(1) = pos1(1) + dt * mom1(1) / std::max(p1, mom1(4));
      pos1(2) = pos1(2) + dt * mom1(2) / std::max(p1, mom1(4));
      pos1(3) = pos1(3) + dt * mom1(3) / std::max(p1, mom1(4));
      pos1(4) = pos1(4) + dt;
      pos2(1) = pos2(1) + dt * mom2(1) / std::max(p2, mom2(4));
      pos2(2) = pos2(2) + dt * mom2(2) / std::max(p2, mom2(4));
      pos2(3) = pos2(3) + dt * mom2(3) / std::max(p2, mom2(4));
      pos2(4) = pos2(4) + dt;
      for (int i = 1; i <= 4; ++i) {
        mom1(i) = P(N - 1, i);
        mom2(i) = P(N, i);
      }
      N = N - 2;
      // C--remove scattering from list
    l502:
      nend = nend - 1;
      if (nend == 1) {
        qsum2 = allqs(nend, 1);
        qsumvec(1) = allqs(nend, 2);
        qsumvec(2) = allqs(nend, 3);
        qsumvec(3) = allqs(nend, 4);
        qsumvec(4) = allqs(nend, 5);
      } else {
        qsum2 = 0.0;
        qsumvec(1) = 0.0;
        qsumvec(2) = 0.0;
        qsumvec(3) = 0.0;
        qsumvec(4) = 0.0;
      }
    }
    if (tleft > 0.0) goto l120;
    // C--no (more) scattering
  }
  // C--split coherent state up again
  tnow = tstart + deltat;
  // C--rotate such that coherent state points in z-direction
  thetam = py::pyp(N, 13);
  phim = py::pyp(N, 15);
  py::pyrobo(N, N, 0.0, -phim, 0.0, 0.0, 0.0);
  py::pyrobo(N, N, -thetam, 0.0, 0.0, 0.0, 0.0);
  mb = P(sib2, 5);
  mc = P(sib1, 5);
  pz = (2.0 * z * ipow<2>(P(N, 4)) - ipow<2>(P(N, 5)) - ipow<2>(mb) + ipow<2>(mc)) /
       (2.0 * P(N, 3));
  pts = ipow<2>(z) * (ipow<2>(P(N, 4))) - ipow<2>(pz) - ipow<2>(mb);
  if (pts < 0.0) io::log(" ERROR: docohscattering has pt^2 =  {}", pts);
  mv(N, 5) = tnow;
  K(N, 1) = 21;
  K(N, 4) = N + 1;
  K(N, 5) = N + 2;
  N = N + 2;
  P(N - 1, 1) = -std::sqrt(pts) * std::cos(phisplit);
  P(N - 1, 2) = -std::sqrt(pts) * std::sin(phisplit);
  P(N - 1, 3) = P(N - 2, 3) - pz;
  P(N - 1, 4) = (1.0 - z) * P(N - 2, 4);
  P(N - 1, 5) = mc;
  P(N, 1) = std::sqrt(pts) * std::cos(phisplit);
  P(N, 2) = std::sqrt(pts) * std::sin(phisplit);
  P(N, 3) = pz;
  P(N, 4) = z * P(N - 2, 4);
  P(N, 5) = mb;
  K(N - 1, 1) = 1;
  K(N - 1, 2) = K(sib1, 2);
  K(N - 1, 3) = N - 2;
  K(N, 1) = 1;
  K(N, 2) = K(sib2, 2);
  K(N, 3) = N - 2;
  K(N - 1, 4) = 0;
  K(N - 1, 5) = 0;
  K(N, 4) = 0;
  K(N, 5) = 0;
  mv(N - 1, 4) = mv(N - 2, 5);
  mv(N, 4) = mv(N - 2, 5);
  mv(N - 1, 5) = decayt1;
  mv(N, 5) = decayt2;
  p1 = std::sqrt(ipow<2>(mom1(1)) + ipow<2>(mom1(2)) + ipow<2>(mom1(3)));
  p2 = std::sqrt(ipow<2>(mom2(1)) + ipow<2>(mom2(2)) + ipow<2>(mom2(3)));
  mv(N - 1, 1) = pos1(1) + (mv(N - 1, 4) - mv(N - 2, 4)) * mom1(1) / std::max(p1, mom1(4));
  mv(N - 1, 2) = pos1(2) + (mv(N - 1, 4) - mv(N - 2, 4)) * mom1(2) / std::max(p1, mom1(4));
  mv(N - 1, 3) = pos1(3) + (mv(N - 1, 4) - mv(N - 2, 4)) * mom1(3) / std::max(p1, mom1(4));
  mv(N, 1) = pos2(1) + (mv(N, 4) - mv(N - 2, 4)) * mom2(1) / std::max(p2, mom2(4));
  mv(N, 2) = pos2(2) + (mv(N, 4) - mv(N - 2, 4)) * mom2(2) / std::max(p2, mom2(4));
  mv(N, 3) = pos2(3) + (mv(N, 4) - mv(N - 2, 4)) * mom2(3) / std::max(p2, mom2(4));
  thetad(N - 1) = thetad(sib1);
  zd(N - 1) = zd(sib1);
  qqbard(N - 1) = qqbard(sib1);
  if (cb::hf::on) cb::hf::qqbarfl(N - 1) = cb::hf::qqbarfl(sib1);
  thetad(N) = thetad(sib2);
  zd(N) = zd(sib2);
  qqbard(N) = qqbard(sib2);
  if (cb::hf::on) cb::hf::qqbarfl(N) = cb::hf::qqbarfl(sib2);
  aoreq(N - 1) = true;
  aoreq(N) = true;
  if (qqbar) {
    trip(N - 1) = 0;
    anti(N - 1) = anti(N - 2);
  } else {
    if ((K(N - 2, 2) > 0) && (dir >= 0)) {
      trip(N - 1) = trip(N - 2);
      anti(N - 1) = colmax + 1;
    } else {
      trip(N - 1) = colmax + 1;
      anti(N - 1) = anti(N - 2);
    }
    colmax = colmax + 1;
  }
  if (quark) {
    if (K(N, 2) > 0) {
      trip(N) = anti(N - 1);
      anti(N) = 0;
    } else {
      trip(N) = 0;
      anti(N) = trip(N - 1);
    }
  } else if (qqbar) {
    trip(N) = trip(N - 2);
    anti(N) = 0;
  } else {
    if (dir == 1) {
      trip(N) = anti(N - 1);
      anti(N) = anti(N - 2);
    } else {
      trip(N) = trip(N - 2);
      anti(N) = trip(N - 1);
    }
  }
  py::pyrobo(N - 2, N, thetam, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(N - 2, N, 0.0, phim, 0.0, 0.0, 0.0);
  if (l == sib1) {
    starttime(N - 1) = starttime(l);
    starttime(N) = tnow;
    l = N - 1;
  } else {
    starttime(N) = starttime(l);
    starttime(N - 1) = tnow;
    l = N;
  }

  return;
}

}  // namespace jewel
