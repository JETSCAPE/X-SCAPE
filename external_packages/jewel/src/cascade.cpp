// JEWEL 2.6.0 C++ port -- parton cascade: MAKECASCADE, MAKEBRANCH, MAKESPLITTING,
// REJECTSPLITTING, MAKEINSPLIT. Translated from jewel-2.6.0.f:2902-4131.
#include <algorithm>
#include <cmath>
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

// ***********************************************************************
// ***	  subroutine makecascade
// ***********************************************************************
void makecascade() {
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::param;
  using namespace cb::disc;

  // C--local variables
  int nold;
  bool cont;

l10:
  nold = N;
  cont = false;
  for (int i = 2, i_end = nold; i <= i_end; ++i) {
    if (i > N) goto l10;
    // C--check if parton may evolve, i.e. do splitting or scattering
    if ((K(i, 1) == 1) || (K(i, 1) == 2)) {
      cont = true;
      makebranch(i);
      if (discard) return;  // GOTO 12
    }
  }
  if (cont) goto l10;
  // 12 END
}

// ***********************************************************************
// ***	  subroutine makebranch
// ***********************************************************************
void makebranch(const int& l) {
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::param;
  using namespace cb::disc;
  using namespace cb::angor;
  using namespace cb::evoltimes;
  using namespace cb::coherent;
  using namespace cb::weight;
  using namespace cb::hepmcid;
  using namespace cb::storescatcen;

  // C--local variables
  int line, nold = 0, typi, lineold = 0, oldstcode = 0, nscatcenold = 0, nold2;
  // NOTE(port): possibly uninitialised in Fortran (LKINE, NEWMASS are only
  // unset on the path X>1 after DOINSTATESCAT, which cannot occur); set to 0.
  int lkine = 0;
  double newmass = 0.0;
  double formtime, tleft, tsum, deltat, q, zdec = 0.0, x, olddecayt = 0.0, eps;
  bool overq0, qqbardec = false, radiation, retrysplit, medindtmp;
  [[maybe_unused]] bool roomleft;
  [[maybe_unused]] std::string typ;

  line = l;
  nstart = 0;
  nend = 0;
  // C--sanity check of starting time
  if (starttime(line) > -1.e5) {
    if (starttime(line) < mv(line, 4) - 1.e-6) {
      io::log(" Error in MAKEBRANCH: starting time before production time {} {} {}", line,
              starttime(line), mv(line, 4));
      pevrec(3, false);
    }
    if ((starttime(line) > mv(line, 5) + 1.e-6) && (mv(line, 5) > 0.0)) {
      io::log(" Error in MAKEBRANCH: starting time after decay time {} {} {}", line,
              starttime(line), mv(line, 5));
      pevrec(3, false);
    }
  }
  if (starttime(line) < -1.e5) {
    if ((mv(line, 4) < 0.0) && (mv(line, 5) > 0.0)) {
      starttime(line) = 0.0;
    } else if ((mv(line, 4) < 0.0) && (!offshell(line))) {
      starttime(line) = 0.0;
    } else {
      starttime(line) = mv(line, 4);
    }
  }
  tsum = 0.0;
  qsum2 = 0.0;
  qsumvec(1) = 0.0;
  qsumvec(2) = 0.0;
  qsumvec(3) = 0.0;
  qsumvec(4) = 0.0;
  retrysplit = false;
  medindtmp = false;
  x = 0.0;
  q = 0.0;
  typi = 0;

  // label 20 (not a jump target)
  if (discard) return;
  if ((N > 20000) && compress) roomleft = compressevent(line);
  if (((K(line, 1) == 1) && offshell(line)) || ((K(line, 1) == 2) && (zd(line) > 0.0))) {
    if (medind(line)) {
      formtime = starttime(line);
    } else {
      formtime = std::min(mv(line, 5), ltime);
    }
    radiation = true;
  } else {
    formtime = ltime;
    radiation = false;
  }
  tleft = formtime - starttime(line);
  if (K(line, 2) == 21) {
    typ = "G";
  } else {
    typ = "Q";
  }
  medind(line) = false;

  if ((tleft <= 1.e-6) || (starttime(line) < 0.0)) {
    // C--no scattering
    if (radiation) {
      // C--if there is radiation associated with the parton then form it now
      nold = N;
      nscatcenold = nscatcen;
      makesplitting(line);
      if (discard) return;
      // C--no new daughters if splitting was rejected due to angular ordering
      if (nold < N) {
        // C--set the production vertices: x_mother + (tprod - tprod_mother) * beta_mother
        mv(N - 1, 1) = mv(line, 1) +
                       (mv(N - 1, 4) - mv(line, 4)) * P(line, 1) / std::max(py::pyp(line, 8), P(line, 4));
        mv(N - 1, 2) = mv(line, 2) +
                       (mv(N - 1, 4) - mv(line, 4)) * P(line, 2) / std::max(py::pyp(line, 8), P(line, 4));
        mv(N - 1, 3) = mv(line, 3) +
                       (mv(N - 1, 4) - mv(line, 4)) * P(line, 3) / std::max(py::pyp(line, 8), P(line, 4));
        mv(N, 1) = mv(line, 1) +
                   (mv(N, 4) - mv(line, 4)) * P(line, 1) / std::max(py::pyp(line, 8), P(line, 4));
        mv(N, 2) = mv(line, 2) +
                   (mv(N, 4) - mv(line, 4)) * P(line, 2) / std::max(py::pyp(line, 8), P(line, 4));
        mv(N, 3) = mv(line, 3) +
                   (mv(N, 4) - mv(line, 4)) * P(line, 3) / std::max(py::pyp(line, 8), P(line, 4));
        line = N;
        if ((mv(N - 1, 4) < 0.0) && (mv(N - 1, 5) > 0.0)) {
          starttime(N - 1) = 0.0;
        } else {
          starttime(N - 1) = mv(N - 1, 4);
        }
        if ((mv(N, 4) < 0.0) && (mv(N, 5) > 0.0)) {
          starttime(line) = 0.0;
        } else {
          starttime(line) = mv(N, 4);
        }
      }
      // GOTO 21
    } else {
      starttime(line) = formtime;
      // GOTO 21
    }
  } else {
    // C--do scattering
    // C--find delta t for the scattering
    deltat = tleft;
    overq0 = false;
    doinstatescat(line, x, typi, q, starttime(line) + tsum, deltat, overq0, false);
    tsum = tsum + deltat;
    tleft = tleft - deltat;
    // C--do initial state splitting if there is one
    nold = N;
    lineold = line;
    oldstcode = K(line, 1);
    zdec = zd(line);
    qqbardec = qqbard(line);
    nscatcenold = nscatcen;
    olddecayt = mv(line, 5);
  l25:
    if (x < 1.0) {
      makeinsplit(line, x, qsum2, q, typi, starttime(line) + tsum, deltat);
      if (discard) return;
      if (x < 1.0) {
        starttime(N - 1) = starttime(line);
        starttime(N) = starttime(line);
        line = N;
        lkine = N;
        if (K(line, 2) == 21) {
          newmass = getmass(0.0, scalefacm * std::sqrt(-qsum2), -1.0, P(line, 4), "GC",
                            std::sqrt(-qsum2), false, zdec, qqbardec, 0.0, &cb::hf::pendingfl);
          if (zdec > 0.0) {
            thetad(line) = newmass / (std::sqrt(zdec * (1. - zdec)) * P(line, 4));
          } else {
            thetad(line) = 0.0;
          }
        } else {
          // heavy quark: massive kernel; newmass is the virtuality (target for dokinematics)
          const double mql = hfmass(K(line, 2));
          newmass = getmass(0.0, scalefacm * std::sqrt(-qsum2), -1.0, P(line, 4), "QQ",
                            std::sqrt(-qsum2), false, zdec, qqbardec, mql, &cb::hf::pendingfl);
          if (zdec > 0.0) {
            thetad(line) = (mql > 0.0) ? hftheta(newmass * newmass, mql, zdec, P(line, 4))
                                       : newmass / (std::sqrt(zdec * (1. - zdec)) * P(line, 4));
          } else {
            thetad(line) = 0.0;
          }
        }
        zd(line) = zdec;
        qqbard(line) = qqbardec;
        if (cb::hf::on) cb::hf::qqbarfl(line) = qqbardec ? cb::hf::pendingfl : 0;
      } else {
        lkine = line;
        nend = nstart;
        qsum2 = allqs(nend, 1);
        qsumvec(1) = allqs(nend, 2);
        qsumvec(2) = allqs(nend, 3);
        qsumvec(3) = allqs(nend, 4);
        qsumvec(4) = allqs(nend, 5);
        if (-allqs(nend, 1) > ipow<2>(q0) / ipow<2>(scalefacm)) {
          overq0 = true;
        } else {
          overq0 = false;
        }
        tleft = starttime(line) + tsum + tleft - allqs(1, 6);
        tsum = allqs(1, 6) - starttime(line);
      }
    }
    if (x == 1.0) {
      newmass = 0.0;
      if (nend > 0) {
        dofistatescat(line, starttime(line) + tsum, tleft, deltat, newmass, overq0, zdec, qqbardec);
        // (newmass is a virtuality; compare with the current one)
        if (newmass > (virtuality(P(line, 5), hfmass(K(line, 2))) * (1.0 + 1.e-6))) {
          medindtmp = true;
        } else {
          medindtmp = false;
          zdec = zd(line);
          qqbardec = qqbard(line);
          cb::hf::pendingfl = cb::hf::qqbarfl(line);
        }
        tsum = tsum + deltat;
        tleft = tleft - deltat;
        lkine = line;
      }
    }
    // C--do kinematics
    retrysplit = false;
    if (nend > 0) {
      [[maybe_unused]] int nendold = nend;
      nold2 = N;
      dokinematics(lkine, lineold, nstart, nend, newmass, retrysplit, starttime(lkine) + tsum, x, zdec,
                   qqbardec, false);
      if (medindtmp && (zd(N) > 0.0) && isheavy(K(N, 2))) {
        // heavy quark: massive z limits and opening angle, with the virtuality
        const double mqn = hfmass(K(N, 2)), t = virt2(N);
        double zlo, zhi;
        if ((t > 0.0) && hfzlimitsqq(t, mqn, P(N, 4), zlo, zhi)) {
          while ((zd(N) < zlo) || (zd(N) > zhi))
            zd(N) = generatez(0.0, 0.0, std::max(std::min(zlo, 1.0 - zhi), 1e-6), "QQ");
          thetad(N) = hftheta(t, mqn, zd(N), P(N, 4));
        } else {
          zd(N) = 0.0;
          thetad(N) = 0.0;
        }
      } else if (medindtmp && (zd(N) > 0.0)) {
        eps = 0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(P(N, 5))) *
                        std::sqrt(1.0 - ipow<2>(P(N, 5)) / ipow<2>(P(N, 4)));
        if ((zd(N) < eps) || (zd(N) > (1.0 - eps))) {
          if (K(N, 2) < 21) {
            zd(N) = generatez(0.0, 0.0, eps, "QQ");
          } else {
            if (qqbard(N)) {
              zd(N) = generatez(0.0, 0.0, eps, "GQ");
            } else {
              zd(N) = generatez(0.0, 0.0, eps, "GG");
            }
          }
        }
        thetad(N) = P(N, 5) / (std::sqrt(zd(N) * (1.0 - zd(N))) * P(N, 4));
      }
      if (retrysplit) {
        tleft = starttime(line) + tsum + tleft - allqs(1, 6);
        tsum = allqs(1, 6) - starttime(line);
        if (x < 1.0) {
          nend = nstart;
          qsum2 = allqs(nend, 1);
          qsumvec(1) = allqs(nend, 2);
          qsumvec(2) = allqs(nend, 3);
          qsumvec(3) = allqs(nend, 4);
          qsumvec(4) = allqs(nend, 5);
          typi = K(l, 2);
          if (-allqs(nend, 1) > ipow<2>(q0) / ipow<2>(scalefacm)) {
            overq0 = true;
          } else {
            overq0 = false;
          }
          N = nold;
          line = lineold;
          x = 1.0;
          K(line, 1) = oldstcode;
          nscatcen = nscatcenold;
          mv(line, 5) = olddecayt;
          goto l25;
        } else {
          starttime(N) = starttime(line) + tsum;
          if (N > nold2) line = N;
          medind(line) = medindtmp;
          tsum = 0.0;
        }
      } else {
        if (N > nold2) {
          if (newmass == virtuality(P(lkine, 5), hfmass(K(lkine, 2)))) {
            starttime(N) = starttime(line) + tsum;
          } else {
            // !	       starttime(n) = mv(n,4)
            if (mv(N, 5) > 0.0) {
              starttime(N) = std::min(mv(N, 5), starttime(line) + tsum);
            } else {
              starttime(N) = starttime(line) + tsum;
            }
          }
          line = N;
        } else {
          if (newmass == virtuality(P(lkine, 5), hfmass(K(lkine, 2)))) {
            starttime(line) = starttime(line) + tsum;
          } else {
            // !	       starttime(line) = mv(line,4)
            if (mv(line, 5) > 0.0) {
              starttime(line) = std::min(mv(line, 5), starttime(line) + tsum);
            } else {
              starttime(line) = starttime(line) + tsum;
            }
          }
        }
        medind(line) = medindtmp;
        tsum = 0.0;
      }
    } else {
      starttime(line) = starttime(line) + tsum;
      tsum = 0.0;
    }
  }
  // label 21
  if (((K(line, 1) == 1) && offshell(line)) || ((K(line, 1) == 2) && (zd(line) > 0.0)) ||
      (starttime(line) < ltime)) {
    return;
  }
  if ((K(line, 1) == 1) && (!offshell(line))) K(line, 1) = 4;
  if ((K(line, 1) == 2) && (zd(line) < 0.0)) K(line, 1) = 5;
}

// ***********************************************************************
// ***	  subroutine makesplitting
// ***********************************************************************
void makesplitting(const int& l) {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::ftimefac;
  using namespace cb::colour;
  using namespace cb::param;
  using namespace cb::disc;
  using namespace cb::angor;
  using namespace cb::evoltimes;
  using namespace cb::weight;

  // C--local variables
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/
  int dir, mother, sib1, sib2, daughter;
  [[maybe_unused]] int premother;
  double phiq, bmax1, cmax1, pts, mb = 0.0, mc = 0.0, pz, z, lambda, zdecb = 0.0, zdecc = 0.0,
                                   theta, thetain, phiin, thetaa, thetab, thetac;
  F1<double, 3> xdec;
  bool quark, qqbar, qqbardecb = false, qqbardecc = false, aoreject, newmb, newmc;
  bool hasscattered, firstdaughter, hasdecayed;
  [[maybe_unused]] bool firstsplit, aoreql;
  // heavy flavour: quark masses of daughters b (quark / quark from g->qqbar) and c (gluon /
  // antiquark), their invariant masses, g->qqbar flavour of this line and of the daughters.
  // mb, mc are virtualities as returned by getmass; for light partons Mb = mb, Mc = mc.
  double mqb = 0.0, mqc = 0.0, Mb = 0.0, Mc = 0.0;
  int hfl = 0, flb = 0, flc = 0;

  // C--on-shell partons cannot split
  if ((!offshell(l)) || (K(l, 1) == 11) || (K(l, 1) == 12) || (K(l, 1) == 13) || (K(l, 1) == 14) ||
      (K(l, 1) == 3) || (zd(l) < 0.0))
    return;

  if ((N + 2) > 22990) {
    io::log(" event too long for event record");
    discard = true;
    return;
  }

  mother = K(l, 3);
  premother = l;
  while (K(mother, 1) == 12) {
    mother = K(mother, 3);
    premother = K(premother, 3);
  }
  sib1 = K(mother, 4);
  sib2 = K(mother, 5);
  if ((sib1 == 0) && (sib2 == 0)) {
    firstsplit = true;
  } else {
    firstsplit = false;
  }
  aoreql = false;

  // C--rotate such that momentum points in z-direction
  thetain = py::pyp(l, 13);
  phiin = py::pyp(l, 15);
  py::pyrobo(l, l, 0.0, -phiin, 0.0, 0.0, 0.0);
  py::pyrobo(l, l, -thetain, 0.0, 0.0, 0.0, 0.0);

  if (angord == 3) {
    xdec(1) = mv(l, 1) + (mv(l, 5) - mv(l, 4)) * P(l, 1) / P(l, 4);
    xdec(2) = mv(l, 2) + (mv(l, 5) - mv(l, 4)) * P(l, 2) / P(l, 4);
    xdec(3) = mv(l, 3) + (mv(l, 5) - mv(l, 4)) * P(l, 3) / P(l, 4);
    if (medium().gettemp(xdec(1), xdec(2), xdec(3), mv(l, 5)) > 0.0) {
      theta = -1.0;
    } else {
      theta = thetad(l);
    }
  } else if (angord == 4) {
    theta = thetad(l);
  } else {
    theta = 0.0;
  }

  // C--quark or gluon?
  if (K(l, 2) == 21) {
    quark = false;
  } else {
    quark = true;
    qqbar = false;
  }
  // C--if gluon decide on kind of splitting
  qqbar = qqbard(l);
  if (cb::hf::on) {
    if (qqbar) hfl = cb::hf::qqbarfl(l);
    mqb = quark ? hfmass(K(l, 2)) : hfmass(hfl);
    mqc = (qqbar && !quark) ? hfmass(hfl) : 0.0;
  }
  // C--if g->gg splitting decide on colour order
  if (quark || qqbar) {
    dir = 0;
  } else {
    if (py::pyr(0) < 0.5) {
      dir = 1;
    } else {
      dir = -1;
    }
  }
  z = zd(l);
  if (z == 0.0) {
    io::log(" makesplitting: z=0 {} {}", l, P(l, 5));
    discard = true;
    return;
  }
  // C--maximum virtualities for daughters
  bmax1 = virtuality(std::min(P(l, 5), z * P(l, 4)), mqb);
  cmax1 = virtuality(std::min(P(l, 5), (1. - z) * P(l, 4)), mqc);
  newmb = true;
  newmc = true;
l30:
  if (newmb) {
    // C--generate mass of quark or gluon (particle b) from Sudakov FF
    if (quark || qqbar) {
      mb = getmass(0.0, bmax1, theta, z * P(l, 4), "QQ", bmax1, false, zdecb, qqbardecb, mqb, &flb);
    } else {
      mb = getmass(0.0, bmax1, theta, z * P(l, 4), "GC", bmax1, false, zdecb, qqbardecb, mqb, &flb);
    }
    Mb = invmass(mb, mqb);
  }
  if (newmc) {
    // C--generate mass gluon (particle c) from Sudakov FF
    if (quark || (!qqbar)) {
      mc = getmass(0.0, cmax1, theta, (1. - z) * P(l, 4), "GC", cmax1, false, zdecc, qqbardecc, mqc, &flc);
    } else {
      mc = getmass(0.0, cmax1, theta, (1. - z) * P(l, 4), "QQ", cmax1, false, zdecc, qqbardecc, mqc, &flc);
    }
    Mc = invmass(mc, mqc);
  }
  // C--quark (parton b) momentum
  // label 182 (not a jump target)
  pz = (2. * z * ipow<2>(P(l, 4)) - ipow<2>(P(l, 5)) - ipow<2>(Mb) + ipow<2>(Mc)) / (2. * P(l, 3));
  pts = ipow<2>(z) * (ipow<2>(P(l, 4))) - ipow<2>(pz) - ipow<2>(Mb);
  // C--if kinematics doesn't work out, generate new virtualities
  // C     for daughters
  if ((pts < 0.0) || ((Mb + Mc) > P(l, 5))) {
    if ((mb == 0.0) && (mc == 0.0)) {
      io::log(" makesplitting: no kinematically allowed splitting found, will reject event.");
      if (std::getenv("JEWEL_HFDEBUG"))
        io::log(" HFDEBUG kf={} qqbar={} hfl={} M={} E={} z={} mqb={} mqc={} Mb={} Mc={} pts={} K1={} mother={}",
                K(l, 2), qqbar, hfl, P(l, 5), P(l, 4), z, mqb, mqc, Mb, Mc, pts, K(l, 1), K(l, 3));
      discard = true;
      goto l31;
    } else if (mb == 0.0) {
      newmb = false;
      newmc = true;
      cmax1 = mc;
      goto l30;
    } else if (mc == 0.0) {
      newmb = true;
      newmc = false;
      bmax1 = mb;
      goto l30;
    } else {
      if (py::pyr(0) > 0.5) {
        newmb = false;
        newmc = true;
        cmax1 = mc;
      } else {
        newmb = true;
        newmc = false;
        bmax1 = mb;
      }
      goto l30;
    }
  }

  // C---START of ANGORD=5 mode
  if (angord == 5) {
    // C--update decay angle of parton in line L with actual angle (including masses of daughters)
    thetaa = std::asin(std::sqrt(pts / (ipow<2>(z) * ipow<2>(P(l, 4)) - ipow<2>(Mb)))) +
             std::asin(std::sqrt(pts / (ipow<2>(1. - z) * ipow<2>(P(l, 4)) - ipow<2>(Mc))));
    thetad(l) = thetaa;
  }
  // C--check if splitting should respect angular ordering and reject if necessary
  if ((angord == 1) || (angord == 5)) {
    // C--find out whether parton has scattered
    hasscattered = false;
    mother = K(l, 3);
    firstdaughter = (K(mother, 4) == l);
    while (K(mother, 1) == 12) {
      if ((K(mother, 1) == 12) && (K(K(mother, 4), 1) != 11)) {
        hasscattered = true;
      }
      firstdaughter = (K(K(mother, 3), 4) == mother);
      mother = K(mother, 3);
    }
    if (hasscattered) {
      // C--current parton has scattered -> no angular ordering
      aoreject = false;
    } else {
      sib1 = K(mother, 4);
      sib2 = K(mother, 5);
      // C--now we check whether this is the first splitting
      if (!aoreq(l)) {
        aoreject = false;
      } else if (thetad(mother) < 0.0) {
        aoreject = false;
      } else {
        // C--check whether sibling has scattered or decayed
        if (firstdaughter) {
          daughter = sib2;
        } else {
          daughter = sib1;
        }
        hasscattered = false;
        while (K(daughter, 1) == 12) {
          if ((K(daughter, 1) == 12) && (K(K(K(daughter, 3), 4), 1) != 11)) {
            hasscattered = true;
          }
          daughter = K(daughter, 5);
        }
        if (K(daughter, 1) == 11) {
          hasdecayed = true;
        } else {
          hasdecayed = false;
        }
        if (hasdecayed && (!hasscattered)) aoreql = true;
        if (hasscattered || hasdecayed) {
          // C--the sibling has scattered or already decayed (in which case both siblings have been checked)
          aoreject = false;
        } else {
          // C--none of the siblings has scattered or decayed
          aoreql = true;
          if ((thetad(sib1) > thetad(mother)) || (thetad(sib2) > thetad(mother))) {
            aoreject = true;
          } else {
            aoreject = false;
          }
          // C--end if (sibling has scattered or decayed)
        }
        // C--end if (first splitting)
      }
      // C--end if parton has scattered
    }
    if (aoreject) {
      py::pyrobo(l, l, thetain, 0.0, 0.0, 0.0, 0.0);
      py::pyrobo(l, l, 0.0, phiin, 0.0, 0.0, 0.0);
      rejectsplitting(l, mother);
      return;
    }
  }
  // C---END OF ANGOD=1 or 5 mode

  if (zdecb > 0.0) {
    thetab = (mqb > 0.0) ? hftheta(mb * mb, mqb, zdecb, z * P(l, 4))
                         : mb / (std::sqrt(zdecb * (1.0 - zdecb)) * z * P(l, 4));
  } else {
    thetab = 0.0;
  }
  if (zdecc > 0.0) {
    thetac = (mqc > 0.0) ? hftheta(mc * mc, mqc, zdecc, (1.0 - z) * P(l, 4))
                         : mc / (std::sqrt(zdecc * (1.0 - zdecc)) * (1.0 - z) * P(l, 4));
  } else {
    thetac = 0.0;
  }
  // C--force angular ordering also in medium with new way of doing AO
  if (angord == 2) {
    if ((thetab > thetad(l)) && (thetad(l) > 0.0)) {
      newmb = true;
      newmc = false;
      bmax1 = mb;
      goto l30;
    }
    if ((thetac > thetad(l)) && (thetad(l) > 0.0)) {
      newmc = true;
      newmb = false;
      cmax1 = mc;
      goto l30;
    }
  }
  N = N + 2;
  // C--take care of first daughter (radiated gluon or antiquark)
  K(N - 1, 1) = 1;
  if (qqbar) {
    K(N - 1, 2) = (hfl > 0) ? -hfl : -1;
    trip(N - 1) = 0;
    anti(N - 1) = anti(l);
  } else {
    K(N - 1, 2) = 21;
    if ((K(l, 2) > 0) && (dir >= 0)) {
      trip(N - 1) = trip(l);
      anti(N - 1) = colmax + 1;
    } else {
      trip(N - 1) = colmax + 1;
      anti(N - 1) = anti(l);
    }
    colmax = colmax + 1;
  }
  K(N - 1, 3) = l;
  K(N - 1, 4) = 0;
  K(N - 1, 5) = 0;
  P(N - 1, 4) = (1 - z) * P(l, 4);
  P(N - 1, 5) = Mc;
  thetad(N - 1) = thetac;
  zd(N - 1) = zdecc;
  qqbard(N - 1) = qqbardecc;
  if (cb::hf::on) cb::hf::qqbarfl(N - 1) = flc;
  aoreq(N - 1) = true;
  // C--take care of second daughter (final quark or gluon or quark from
  // C	 gluon splitting)
  // !	K(N,1)=K(L,1)
  K(N, 1) = 1;
  if (quark) {
    K(N, 2) = K(l, 2);
    if (K(N, 2) > 0) {
      trip(N) = anti(N - 1);
      anti(N) = 0;
    } else {
      trip(N) = 0;
      anti(N) = trip(N - 1);
    }
  } else if (qqbar) {
    K(N, 2) = (hfl > 0) ? hfl : 1;
    trip(N) = trip(l);
    anti(N) = 0;
  } else {
    K(N, 2) = 21;
    if (dir == 1) {
      trip(N) = anti(N - 1);
      anti(N) = anti(l);
    } else {
      trip(N) = trip(l);
      anti(N) = trip(N - 1);
    }
  }
  K(N, 3) = l;
  K(N, 4) = 0;
  K(N, 5) = 0;
  P(N, 3) = pz;
  P(N, 4) = z * P(l, 4);
  P(N, 5) = Mb;
  thetad(N) = thetab;
  zd(N) = zdecb;
  qqbard(N) = qqbardecb;
  if (cb::hf::on) cb::hf::qqbarfl(N) = flb;
  aoreq(N) = true;
  // C--azimuthal angle
  phiq = 2 * pi * py::pyr(0);
  P(N, 1) = std::sqrt(pts) * gf_cos(phiq);
  P(N, 2) = std::sqrt(pts) * gf_sin(phiq);
  // C--gluon momentum
  P(N - 1, 1) = P(l, 1) - P(N, 1);
  P(N - 1, 2) = P(l, 2) - P(N, 2);
  P(N - 1, 3) = P(l, 3) - P(N, 3);
  mv(N - 1, 4) = mv(l, 5);
  // formation time from the virtuality (= P(.,5) for light partons)
  if (mc > 0.0) {
    lambda = 1.0 / (ftfac * P(N - 1, 4) * r4(0.2f) / ipow<2>(mc));
    mv(N - 1, 5) = mv(l, 5) - std::log(1.0 - py::pyr(0)) / lambda;
  } else {
    mv(N - 1, 5) = 0.0;
  }
  mv(N, 4) = mv(l, 5);
  if (mb > 0.0) {
    lambda = 1.0 / (ftfac * P(N, 4) * r4(0.2f) / ipow<2>(mb));
    mv(N, 5) = mv(l, 5) - std::log(1.0 - py::pyr(0)) / lambda;
  } else {
    mv(N, 5) = 0.0;
  }
  // C--take care of initial quark (or gluon)
  if (K(l, 1) == 2) {
    K(l, 1) = 13;
  } else {
    K(l, 1) = 11;
  }
  K(l, 4) = N - 1;
  K(l, 5) = N;
l31:
  py::pyrobo(l, l, thetain, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(l, l, 0.0, phiin, 0.0, 0.0, 0.0);
  py::pyrobo(N - 1, N, thetain, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(N - 1, N, 0.0, phiin, 0.0, 0.0, 0.0);
}

// ***********************************************************************
// ***	  subroutine rejectsplitting
// ***********************************************************************
void rejectsplitting(const int& line, const int& mother) {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::param;
  using namespace cb::disc;
  using namespace cb::angor;
  using namespace cb::evoltimes;

  // C--local variables
  int sib1, sib2, ltemp;
  double mbnew, mcnew, zdecb = 0.0, zdecc = 0.0, bmax1, cmax1, theta, phi, pts, pznew, ptsnew, deltat1,
                       deltat2, mb, mc;
  bool qqbardecb = false, qqbardecc = false, newmb, newmc;
  // heavy flavour: quark masses of the siblings; mb, mc, mbnew, mcnew are virtualities,
  // Mbnew, Mcnew invariant masses (identical for light partons)
  double mqb, mqc, Mbnew, Mcnew;
  int flb = 0, flc = 0;
  [[maybe_unused]] bool bviolao, cviolao;

  // C--find sibling
  sib1 = K(mother, 4);
  sib2 = K(mother, 5);
  // C--if particle had rejected scatterings set status codes of these lines to -1
  ltemp = sib1;
  while (K(ltemp, 1) == 12) {
    K(K(ltemp, 4), 1) = -1;
    K(K(ltemp, 5), 1) = -1;
    ltemp = K(ltemp, 3);
  }
  ltemp = sib2;
  while (K(ltemp, 1) == 12) {
    K(K(ltemp, 4), 1) = -1;
    K(K(ltemp, 5), 1) = -1;
    ltemp = K(ltemp, 3);
  }
  if (K(sib1, 1) == 12) K(sib1, 1) = 1;
  if (K(sib2, 1) == 12) K(sib2, 1) = 1;
  // C--rotate such that mother points in z-direction (for kinematics)
  theta = py::pyp(mother, 13);
  phi = py::pyp(mother, 15);
  py::pyrobo(mother, mother, 0.0, -phi, 0.0, 0.0, 0.0);
  py::pyrobo(mother, mother, -theta, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(sib1, sib2, 0.0, -phi, 0.0, 0.0, 0.0);
  py::pyrobo(sib1, sib2, -theta, 0.0, 0.0, 0.0, 0.0);
  pts = ipow<2>(P(sib1, 1)) + ipow<2>(P(sib1, 2));
  mqb = hfmass(K(sib2, 2));
  mqc = hfmass(K(sib1, 2));
  mb = virtuality(P(sib2, 5), mqb);
  mc = virtuality(P(sib1, 5), mqc);
  mbnew = mb;
  mcnew = mc;
  Mbnew = P(sib2, 5);
  Mcnew = P(sib1, 5);
  newmb = (thetad(sib2) > thetad(mother));
  newmc = (thetad(sib1) > thetad(mother));
  bviolao = newmb;
  cviolao = newmc;
  bmax1 = virtuality(std::min(P(sib2, 5), P(sib2, 4)), mqb);
  cmax1 = virtuality(std::min(P(sib1, 5), P(sib1, 4)), mqc);
  // C--find new virtual masses for the legs that violate angular ordering or kinematics
l40:
  if (newmc) {
    // C--are we dealing with a quark or gluon?
    if (K(sib1, 2) == 21) {
      mcnew = getmass(0.0, mc, thetad(mother), P(sib1, 4), "GC", cmax1, false, zdecc, qqbardecc, mqc, &flc);
    } else {
      mcnew = getmass(0.0, mc, thetad(mother), P(sib1, 4), "QQ", cmax1, false, zdecc, qqbardecc, mqc, &flc);
    }
    Mcnew = invmass(mcnew, mqc);
    P(sib1, 5) = Mcnew;
    zd(sib1) = zdecc;
    qqbard(sib1) = qqbardecc;
    if (cb::hf::on) cb::hf::qqbarfl(sib1) = flc;
    if (zdecc > 0.0) {
      thetad(sib1) = (mqc > 0.0) ? hftheta(mcnew * mcnew, mqc, zdecc, P(sib1, 4))
                                 : mcnew / (std::sqrt(zdecc * (1. - zdecc)) * P(sib1, 4));
    } else {
      thetad(sib1) = 0.0;
    }
  }
  if (newmb) {
    // C--are we dealing with a quark or gluon?
    if (K(sib2, 2) == 21) {
      mbnew = getmass(0.0, mb, thetad(mother), P(sib2, 4), "GC", bmax1, false, zdecb, qqbardecb, mqb, &flb);
    } else {
      mbnew = getmass(0.0, mb, thetad(mother), P(sib2, 4), "QQ", bmax1, false, zdecb, qqbardecb, mqb, &flb);
    }
    Mbnew = invmass(mbnew, mqb);
    P(sib2, 5) = Mbnew;
    zd(sib2) = zdecb;
    qqbard(sib2) = qqbardecb;
    if (cb::hf::on) cb::hf::qqbarfl(sib2) = flb;
    if (zdecb > 0.0) {
      thetad(sib2) = (mqb > 0.0) ? hftheta(mbnew * mbnew, mqb, zdecb, P(sib2, 4))
                                 : mbnew / (std::sqrt(zdecb * (1. - zdecb)) * P(sib2, 4));
    } else {
      thetad(sib2) = 0.0;
    }
  }
  pznew = (2. * zd(mother) * ipow<2>(P(mother, 4)) - ipow<2>(P(mother, 5)) - ipow<2>(Mbnew) +
           ipow<2>(Mcnew)) /
          (2. * P(mother, 3));
  ptsnew = ipow<2>(zd(mother)) * ipow<2>(P(mother, 4)) - ipow<2>(pznew) - ipow<2>(Mbnew);
  if (ptsnew < 0.0) {
    if ((mbnew == 0.0) && (mcnew == 0.0)) {
      io::log(" WARNING: cannot find valid kinematics after AO rejection of splitting -> will discard event. {} {}",
              line, ptsnew);
      discard = true;
      return;
    } else if (mbnew == 0.0) {
      cmax1 = mcnew;
      newmb = false;
      newmc = true;
      goto l40;
    } else if (mcnew == 0.0) {
      bmax1 = mbnew;
      newmc = false;
      newmb = true;
      goto l40;
    } else {
      if (py::pyr(0) > 0.5) {
        newmb = false;
        newmc = true;
        cmax1 = mcnew;
      } else {
        newmb = true;
        newmc = false;
        bmax1 = mbnew;
      }
      goto l40;
    }
  }
  P(sib2, 1) = P(sib2, 1) * std::sqrt(ptsnew / pts);
  P(sib2, 2) = P(sib2, 2) * std::sqrt(ptsnew / pts);
  P(sib2, 3) = pznew;
  P(sib1, 1) = P(mother, 1) - P(sib2, 1);
  P(sib1, 2) = P(mother, 2) - P(sib2, 2);
  P(sib1, 3) = P(mother, 3) - P(sib2, 3);
  // C--update lifetimes
  if (newmc) {
    starttime(sib1) = mv(sib1, 5);
  }
  if (newmb) {
    starttime(sib2) = mv(sib2, 5);
  }
  if (mcnew == 0.0) {
    mv(sib1, 5) = 0.0;
  } else {
    deltat1 = mv(sib1, 5) - mv(sib1, 4);
    mv(sib1, 5) = mv(sib1, 4) + deltat1 * ipow<2>(mc) / ipow<2>(mcnew);
  }
  if (mbnew == 0.0) {
    mv(sib2, 5) = 0.0;
  } else {
    deltat2 = mv(sib2, 5) - mv(sib2, 4);
    mv(sib2, 5) = mv(sib2, 4) + deltat2 * ipow<2>(mb) / ipow<2>(mbnew);
  }
  // C--rotate back
  py::pyrobo(mother, mother, theta, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(mother, mother, 0.0, phi, 0.0, 0.0, 0.0);
  py::pyrobo(sib1, sib2, theta, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(sib1, sib2, 0.0, phi, 0.0, 0.0, 0.0);
}

// ***********************************************************************
// ***	  subroutine makeinsplit
// ***********************************************************************
void makeinsplit(const int& l, double& x, const double& tsum, const double& virt, const int& typi,
                 const double& time, const double& taurad) {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::ftimefac;
  using namespace cb::colour;
  using namespace cb::angor;
  using namespace cb::disc;
  using namespace cb::param;
  using namespace cb::weight;

  // C--local variables
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/
  int dir;
  [[maybe_unused]] int nold;
  double mb2, mc2, pz, kt2, theta, phi, phiq, lambda, zdec = 0.0;
  bool qqbardec = false;
  std::string typ2, typc;

  if ((N + 2) > 22990) {
    io::log(" event too long for event record");
    discard = true;
    return;
  }

  if (K(l, 2) == 21) {
    if (typi == 21) {
      typ2 = "GG";
      typc = "GC";
    } else {
      typ2 = "QG";
      typc = "QQ";
    }
  } else {
    if (typi == 21) {
      typ2 = "GQ";
      typc = "QQ";
    } else {
      typ2 = "QQ";
      typc = "GC";
    }
  }

  // C--if g->gg decide on colour configuration
  if (typ2 == "GG") {
    if (py::pyr(0) < 0.5) {
      dir = 1;
    } else {
      dir = -1;
    }
  } else {
    dir = 0;
  }

  mb2 = ipow<2>(virt);
  mb2 = ipow<2>(P(l, 5)) - mb2;
  // !	MB2=-VIRT**2
  mc2 = ipow<2>(getmass(0.0, scalefacm * std::sqrt(-tsum), -1.0, (1. - x) * P(l, 4), typc,
                        (1. - x) * P(l, 4), false, zdec, qqbardec, 0.0, &cb::hf::pendingfl));

  // C--rotate such that momentum points in z-direction
  nold = N;
  theta = py::pyp(l, 13);
  phi = py::pyp(l, 15);
  py::pyrobo(l, l, 0.0, -phi, 0.0, 0.0, 0.0);
  py::pyrobo(l, l, -theta, 0.0, 0.0, 0.0, 0.0);
  pz = (2 * x * ipow<2>(P(l, 4)) - ipow<2>(P(l, 5)) - mb2 + mc2) / (2 * P(l, 3));
  kt2 = ipow<2>(x) * (ipow<2>(P(l, 4))) - ipow<2>(pz) - mb2;
  if (kt2 < 0.0) {
    mc2 = 0.0;
    if (K(l, 1) == 2) zdec = -1.0;
    pz = (2 * x * ipow<2>(P(l, 4)) - ipow<2>(P(l, 5)) - mb2 + mc2) / (2 * P(l, 3));
    kt2 = ipow<2>(x) * (ipow<2>(P(l, 4))) - ipow<2>(pz) - mb2;
    if (kt2 < 0.0) {
      py::pyrobo(l, l, theta, 0.0, 0.0, 0.0, 0.0);
      py::pyrobo(l, l, 0.0, phi, 0.0, 0.0, 0.0);
      x = 1.0;
      return;
    }
  }
  N = N + 2;
  // C--take care of first daughter (radiated gluon or antiquark)
  // !	K(N-1,1)=K(L,1)
  K(N - 1, 1) = 1;
  if (typ2 == "QG") {
    K(N - 1, 2) = -typi;
    if (K(N - 1, 2) > 0) {
      trip(N - 1) = trip(l);
      anti(N - 1) = 0;
    } else {
      trip(N - 1) = 0;
      anti(N - 1) = anti(l);
    }
  } else if (typ2 == "GQ") {
    K(N - 1, 2) = K(l, 2);
    if (K(N - 1, 2) > 0) {
      trip(N - 1) = colmax + 1;
      anti(N - 1) = 0;
    } else {
      trip(N - 1) = 0;
      anti(N - 1) = colmax + 1;
    }
    colmax = colmax + 1;
  } else {
    K(N - 1, 2) = 21;
    if ((K(l, 2) > 0) && (dir >= 0)) {
      trip(N - 1) = trip(l);
      anti(N - 1) = colmax + 1;
    } else {
      trip(N - 1) = colmax + 1;
      anti(N - 1) = anti(l);
    }
    colmax = colmax + 1;
  }
  K(N - 1, 3) = l;
  K(N - 1, 4) = 0;
  K(N - 1, 5) = 0;
  P(N - 1, 4) = (1. - x) * P(l, 4);
  P(N - 1, 5) = std::sqrt(mc2);
  // C--take care of second daughter (final quark or gluon or quark from
  // C	 gluon splitting)
  // !	K(N,1)=K(L,1)
  K(N, 1) = 1;
  if (typ2 == "QG") {
    K(N, 2) = typi;
    if (K(N, 2) > 0) {
      trip(N) = trip(l);
      anti(N) = 0;
    } else {
      trip(N) = 0;
      anti(N) = anti(l);
    }
  } else if (typi != 21) {
    K(N, 2) = K(l, 2);
    if (K(N, 2) > 0) {
      trip(N) = anti(N - 1);
      anti(N) = 0;
    } else {
      trip(N) = 0;
      anti(N) = trip(N - 1);
    }
  } else {
    K(N, 2) = 21;
    if (K(N - 1, 2) == 21) {
      if (dir == 1) {
        trip(N) = anti(N - 1);
        anti(N) = anti(l);
      } else {
        trip(N) = trip(l);
        anti(N) = trip(N - 1);
      }
    } else if (K(N - 1, 2) > 0) {
      trip(N) = trip(l);
      anti(N) = trip(N - 1);
    } else {
      trip(N) = anti(N - 1);
      anti(N) = anti(l);
    }
  }
  K(N, 3) = l;
  K(N, 4) = 0;
  K(N, 5) = 0;
  P(N, 3) = pz;
  P(N, 4) = x * P(l, 4);
  if (mb2 < 0.0) {
    P(N, 5) = -std::sqrt(-mb2);
  } else {
    P(N, 5) = std::sqrt(mb2);
  }
  // C--azimuthal angle
  phiq = 2 * pi * py::pyr(0);
  P(N, 1) = std::sqrt(kt2) * gf_cos(phiq);
  P(N, 2) = std::sqrt(kt2) * gf_sin(phiq);
  // C--gluon momentum
  P(N - 1, 1) = P(l, 1) - P(N, 1);
  P(N - 1, 2) = P(l, 2) - P(N, 2);
  P(N - 1, 3) = P(l, 3) - P(N, 3);
  mv(l, 5) = time - taurad;
  mv(N - 1, 4) = mv(l, 5);
  if (P(N - 1, 5) > 0.0) {
    lambda = 1.0 / (ftfac * P(N - 1, 4) * r4(0.2f) / ipow<2>(P(N - 1, 5)));
    mv(N - 1, 5) = mv(l, 5) - std::log(1.0 - py::pyr(0)) / lambda;
  } else {
    mv(N - 1, 5) = 0.0;
  }
  mv(N, 4) = mv(l, 5);
  if (offshell(N)) {  // (heavy quark: timelike if M > m)
    mv(N, 5) = time;
  } else {
    mv(N, 5) = 0.0;
  }
  zd(N - 1) = zdec;
  qqbard(N - 1) = qqbardec;
  if (cb::hf::on) {
    cb::hf::qqbarfl(N - 1) = qqbardec ? cb::hf::pendingfl : 0;
    cb::hf::qqbarfl(N) = 0;
  }
  aoreq(N - 1) = false;
  if (zdec == 0.0) {
    thetad(N - 1) = -1.0;
  } else {
    thetad(N - 1) = P(N - 1, 5) / (std::sqrt(zdec * (1.0 - zdec)) * P(N - 1, 4));
  }
  thetad(N) = -1.0;
  zd(N) = 0.0;
  qqbard(N) = false;
  aoreq(N) = false;
  // C--take care of initial quark (or gluon)
  if (K(l, 1) == 2) {
    K(l, 1) = 13;
  } else {
    K(l, 1) = 11;
  }
  K(l, 4) = N - 1;
  K(l, 5) = N;
  py::pyrobo(l, l, theta, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(N - 1, N, theta, 0.0, 0.0, 0.0, 0.0);
  py::pyrobo(l, l, 0.0, phi, 0.0, 0.0, 0.0);
  py::pyrobo(N - 1, N, 0.0, phi, 0.0, 0.0, 0.0);

  // C--set the production vertices: x_mother + (tprod - tprod_mother) * beta_mother
  mv(N - 1, 1) = mv(l, 1) + (mv(N - 1, 4) - mv(l, 4)) * P(l, 1) / std::max(py::pyp(l, 8), P(l, 4));
  mv(N - 1, 2) = mv(l, 2) + (mv(N - 1, 4) - mv(l, 4)) * P(l, 2) / std::max(py::pyp(l, 8), P(l, 4));
  mv(N - 1, 3) = mv(l, 3) + (mv(N - 1, 4) - mv(l, 4)) * P(l, 3) / std::max(py::pyp(l, 8), P(l, 4));
  mv(N, 1) = mv(l, 1) + (mv(N, 4) - mv(l, 4)) * P(l, 1) / std::max(py::pyp(l, 8), P(l, 4));
  mv(N, 2) = mv(l, 2) + (mv(N, 4) - mv(l, 4)) * P(l, 2) / std::max(py::pyp(l, 8), P(l, 4));
  mv(N, 3) = mv(l, 3) + (mv(N, 4) - mv(l, 4)) * P(l, 3) / std::max(py::pyp(l, 8), P(l, 4));
}

}  // namespace jewel
