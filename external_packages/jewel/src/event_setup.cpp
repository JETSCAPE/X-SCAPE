// JEWEL 2.6.0 C++ port -- event set-up: interpretation of the PYTHIA event, special lines,
// colour indices, initial virtualities, string formation and clean-up.
// Translated from jewel-2.6.0.f:1405-2901.
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
// ***	  subroutine interpretpythiaevent
// ***********************************************************************
void interpretpythiaevent() {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace cb::param;
  using namespace cb::jetpoint;
  using namespace cb::evrecord;
  using namespace cb::ftimefac;
  using namespace cb::time;
  using namespace cb::angor;

  // C--local parameters
  int nnew = 0, lstart = 0, lprev = 0, l1 = 0, ntmp = 0, mother = 0, slen = 0;
  // NOTE(port): possibly uninitialised in Fortran (l2: set in loop 302 only if a
  // matching line is found; l3: set in loop 410 only if a line is found).
  int l2 = 0, l3 = 0;
  double lambda = 0, pt2 = 0, eps = 0;
  bool onlyzeros = false, recomb = false;

  // C--special treatment for Jeweling intial state radiation (currently only available for di-jets)
  if (collider != "EEJJ") {
    nnew = offset;
    // C--find vector bosons and their decay products and move them up first
    if ((collider == "PPZJ") || (collider == "PPZQ") || (collider == "PPZG") ||
        (collider == "PPWJ") || (collider == "PPWQ") || (collider == "PPWG")) {
      for (int i = 101, i_end = 106; i <= i_end; ++i) {
        if ((K(i, 2) == 23) || (std::abs(K(i, 2)) == 24)) {
          nnew = nnew + 1;
          copyline(i, nnew, 0);
          K(nnew, 3) = K(i, 3);
          K(nnew, 4) = nnew + 1;
          K(nnew, 5) = nnew + 2;
          K(i, 1) = 11;
          // C--now come the decay leptons
          nnew = nnew + 1;
          copyline(K(i, 4), nnew, 0);
          K(nnew, 3) = K(K(i, 4), 3);
          K(K(i, 4), 1) = 11;
          nnew = nnew + 1;
          copyline(K(i, 5), nnew, 0);
          K(nnew, 3) = K(K(i, 5), 3);
          K(K(i, 5), 1) = 11;
        }
      }
    }
    if ((collider == "PPYJ") || (collider == "PPYQ") || (collider == "PPYG")) {
      for (int i = 101, i_end = 106; i <= i_end; ++i) {
        if (K(i, 2) == 22) {
          nnew = nnew + 1;
          copyline(i, nnew, 0);
          K(nnew, 3) = K(i, 3);
          K(i, 1) = 11;
        }
      }
    }
    // C--find other non-strongly interacting particles and move them up
    for (int i = 107, i_end = N; i <= i_end; ++i) {
      if ((K(i, 1) <= 2) && (!isparton(K(i, 2)))) {
        nnew = nnew + 1;
        copyline(i, nnew, 0);
        K(i, 1) = 11;
      }
    }
    // C--walk backwards from ME to beam and set production vertices and times
    for (int i = 0, i_end = 1; i <= i_end; ++i) {
      lstart = 105 + i;
      lambda = 1.0 / (ftfac * P(lstart, 4) * r4(0.2f) / ipow<2>(P(lstart, 5)));
      if (isrscat) {
        mv(lstart, 5) = 0.0;
      } else {
        mv(lstart, 5) = -1.e6;
      }
      mv(lstart, 4) = mv(lstart, 5) + std::log(1.0 - py::pyr(0)) / lambda;
      mv(lstart, 1) = x0 + mv(lstart, 4) * P(lstart, 1) / std::max(py::pyp(lstart, 8), P(lstart, 4));
      mv(lstart, 2) = y0 + mv(lstart, 4) * P(lstart, 2) / std::max(py::pyp(lstart, 8), P(lstart, 4));
      mv(lstart, 3) = mv(lstart, 4) * P(lstart, 3) / std::max(py::pyp(lstart, 8), P(lstart, 4));
      lprev = lstart;
      l1 = K(lprev, 3);
      while (K(lprev, 3) > 4) {
        for (int jj = 105, jj_end = N; jj <= jj_end; ++jj) {
          if (K(jj, 3) == l1) l2 = jj;
        }
        if ((K(l1, 3) == 3) || (K(l1, 3) == 4)) {
          mv(l1, 1) = 1.e6;
          mv(l1, 2) = 1.e6;
          if (K(l1, 3) == 3) {
            mv(l1, 3) = 1.e6;
          } else {
            mv(l1, 3) = -1.e6;
          }
          mv(l1, 4) = -1.e6;
          mv(l1, 5) = mv(lprev, 4);
        } else {
          lambda = 1.0 / (ftfac * P(l1, 4) * r4(0.2f) / ipow<2>(P(l1, 5)));
          mv(l1, 4) = mv(lprev, 4) + std::log(1.0 - py::pyr(0)) / lambda;
          mv(l1, 5) = mv(lprev, 4);
          mv(l1, 1) = mv(lprev, 1) + (mv(l1, 4) - mv(l1, 5)) * P(l1, 1) / std::max(py::pyp(l1, 8), P(l1, 4));
          mv(l1, 2) = mv(lprev, 2) + (mv(l1, 4) - mv(l1, 5)) * P(l1, 2) / std::max(py::pyp(l1, 8), P(l1, 4));
          mv(l1, 3) = mv(lprev, 3) + (mv(l1, 4) - mv(l1, 5)) * P(l1, 3) / std::max(py::pyp(l1, 8), P(l1, 4));
        }
        mv(l2, 4) = mv(lprev, 4);
        if (P(l2, 5) == py::pymass(K(l2, 2))) {
          mv(l2, 5) = ltime;
        } else {
          lambda = 1.0 / (ftfac * P(l2, 4) * r4(0.2f) / virt2(l2));
          mv(l2, 5) = mv(lprev, 4) - std::log(1.0 - py::pyr(0)) / lambda;
        }
        mv(l2, 1) = mv(lprev, 1);
        mv(l2, 2) = mv(lprev, 2);
        mv(l2, 3) = mv(lprev, 3);
        lprev = l1;
        l1 = K(lprev, 3);
      }
    }
    // C--order all singlets from triplet to antitriplet
    for (int i = lstart, i_end = N; i <= i_end; ++i) {
      if ((K(i, 1) == 2) && ((K(i, 2) < 0) || isdiquark(K(i, 2)))) {
        ntmp = N;
        for (int jj = i, jj_end = N; jj <= jj_end; ++jj) {
          ntmp = ntmp + 1;
          copyline(jj, ntmp, 1);
          K(ntmp, 3) = K(jj, 3);
          if (K(jj, 1) == 1) break;  // goto 403
        }
        // 403
        for (int jj = 0, jj_end = ntmp - N - 1; jj <= jj_end; ++jj) {
          copyline(ntmp - jj, i + jj, 1);
          K(i + jj, 3) = K(ntmp - jj, 3);
          if (jj == ntmp - N - 1) {
            K(i + jj, 1) = 1;
          } else {
            K(i + jj, 1) = 2;
          }
        }
      }
    }
    // C--undo final state splittings
    lstart = 0;
    for (int i = 107, i_end = N; i <= i_end; ++i) {
      if ((lstart == 0) && (K(i, 1) <= 2)) lstart = i;
      if (K(i, 1) <= 2) {
        for (int jj = 1; jj <= 5; ++jj) {
          mv(i, jj) = mv(K(i, 3), jj);
        }
        K(i, 3) = K(K(i, 3), 3);
      }
    }
    recomb = true;
    while (recomb) {
      recomb = false;
      for (int i = lstart, i_end = N; i <= i_end; ++i) {
        if (K(i, 3) > 106) {
          onlyzeros = true;
          for (int jj = i + 1, jj_end = N; jj <= jj_end; ++jj) {
            if (K(i, 3) == K(jj, 3)) {
              l1 = i;
              l2 = jj;
              recomb = true;
              goto l308;
            } else {
              if ((jj > i + 1) && (K(jj, 1) > 0)) onlyzeros = false;
            }
          }
        }
      }
    l308:
      if (recomb) {
        if (onlyzeros) {
          if (cb::hf::on) cb::hf::qqbarfl(l1) = 0;
          if (K(K(l1, 3), 2) == 21) {
            zd(l1) = P(l1, 4) / P(K(l1, 3), 4);
            if (K(l1, 2) != 21) {
              qqbard(l1) = true;
              // heavy flavour: keep the flavour of PYTHIA's g -> Q Qbar
              // (all flavours: PYTHIA already chose them, light ones included)
              if (std::abs(K(l1, 2)) <= 5) cb::hf::qqbarfl(l1) = std::abs(K(l1, 2));
            } else {
              qqbard(l1) = false;
            }
          } else {
            if (K(l1, 2) == 21) {
              zd(l1) = P(l2, 4) / P(K(l2, 3), 4);
            } else {
              zd(l1) = P(l1, 4) / P(K(l1, 3), 4);
            }
            qqbard(l1) = false;
          }
          for (int i = 1; i <= 5; ++i) {
            P(l1, i) = P(K(l1, 3), i);
          }
          K(l1, 2) = K(K(l1, 3), 2);
          K(l1, 3) = K(K(l1, 3), 3);
          if (K(l1, 2) == 21) K(l1, 1) = 2;
          if (K(l2, 1) == 1) K(l1, 1) = 1;
          K(l2, 1) = 0;
        } else {
          // C--have to change order in which singlets appear
          if (K(l1, 1) == 2) {
            ntmp = N;
            for (int i = l1, i_end = N; i <= i_end; ++i) {
              ntmp = ntmp + 1;
              copyline(i, ntmp, 1);
              K(ntmp, 3) = K(i, 3);
              if (K(i, 1) == 1) break;  // goto 406
            }
            // 406
            slen = ntmp - N;
            for (int i = l1 + slen, i_end = l2; i <= i_end; ++i) {
              copyline(i, i - slen, 1);
              K(i - slen, 3) = K(i, 3);
            }
            for (int i = 1, i_end = slen; i <= i_end; ++i) {
              copyline(N + i, l2 - slen + i, 1);
              K(l2 - slen + i, 3) = K(N + i, 3);
            }
          } else {
            for (int i = 1, i_end = l1 - lstart; i <= i_end; ++i) {
              if ((K(l1 - i, 2) != 21) && (K(l1 - i, 1) != 0)) {
                l3 = l1 - i;
                break;  // goto 411
              }
            }
            // 411
            ntmp = N;
            for (int i = l3, i_end = l1; i <= i_end; ++i) {
              ntmp = ntmp + 1;
              copyline(i, ntmp, 1);
              K(ntmp, 3) = K(i, 3);
            }
            slen = ntmp - N;
            for (int i = l3 + slen, i_end = l2 - 1; i <= i_end; ++i) {
              copyline(i, i - slen, 1);
              K(i - slen, 3) = K(i, 3);
            }
            for (int i = 1, i_end = slen; i <= i_end; ++i) {
              copyline(N + i, l2 - slen - 1 + i, 1);
              K(l2 - slen - 1 + i, 3) = K(N + i, 3);
            }
          }
        }
      }
    }
    // C--copy remaining lines to top of event record
    for (int jj = lstart, jj_end = N; jj <= jj_end; ++jj) {
      if (K(jj, 1) > 0) {
        if ((K(jj, 3) > 106) && (P(jj, 5) != py::pymass(K(jj, 2)))) {
          for (int i = 1; i <= 5; ++i) {
            mv(jj, i) = mv(K(jj, 3), i);
          }
        }
        if ((K(jj, 3) == 1) || (K(jj, 3) == 2)) {
          mv(jj, 1) = x0;
          mv(jj, 2) = y0;
          mv(jj, 3) = 0.0;
          mv(jj, 4) = 0.0;
          mv(jj, 5) = ltime;
        }
        nnew = nnew + 1;
        copyline(jj, nnew, 1);
        if (K(jj, 3) < 100) {
          K(nnew, 3) = K(jj, 3);
        } else {
          mother = K(jj, 3);
          while (mother > 100) {
            mother = K(mother, 3);
          }
          K(nnew, 3) = mother;
        }
        zd(nnew) = zd(jj);
        qqbard(nnew) = qqbard(jj);
        aoreq(nnew) = false;
        if (cb::hf::on) cb::hf::qqbarfl(nnew) = cb::hf::qqbarfl(jj);
        if (zd(nnew) > 0.0) {
          thetad(nnew) = isheavy(K(nnew, 2))
                             ? hftheta(virt2(nnew), hfmass(K(nnew, 2)), zd(nnew), P(nnew, 4))
                             : P(nnew, 5) / (std::sqrt(zd(nnew) * (1.0 - zd(nnew))) * P(nnew, 4));
        } else {
          thetad(nnew) = 0.0;
        }
      }
    }
    N = nnew;
    npart = N - offset;
  }

  // C--extra check to see whether all splittings are kinematically allowed
  for (int jj = offset, jj_end = N; jj <= jj_end; ++jj) {
    // heavy g -> Q Qbar flavour (light flavours need no special treatment here)
    const int hfl = (cb::hf::on && (K(jj, 2) == 21) && qqbard(jj) && isheavy(cb::hf::qqbarfl(jj))) ? cb::hf::qqbarfl(jj) : 0;
    if (isheavy(K(jj, 2)) || (hfl > 0)) {
      // heavy flavour: same check with the massive z limits.  These splittings come from
      // PYTHIA's shower history and may lie below JEWEL's Q0 resolution; then only the
      // kinematic limits for on-shell daughters apply.  z is redrawn inside the range.
      const double t = virt2(jj), e = P(jj, 4);
      const double m = (hfl > 0) ? hfmass(hfl) : hfmass(K(jj, 2));
      double zlo = 0.0, zhi = 0.0;
      bool open;
      if (hfl > 0) {
        open = hfzlimitsgq(t, m, e, zlo, zhi) || hfzlimits(t, 0.0, e, m, m, zlo, zhi);
      } else {
        open = hfzlimitsqq(t, m, e, zlo, zhi) || hfzlimits(t, m, e, m, 0.0, zlo, zhi);
      }
      if ((t > 0.0) && open && ((zd(jj) < zlo) || (zd(jj) > zhi))) {
        do {
          zd(jj) = generatez(0.0, 0.0, std::max(std::min(zlo, 1.0 - zhi), 1e-6), (hfl > 0) ? "GQ" : "QQ");
        } while ((zd(jj) < zlo) || (zd(jj) > zhi));
      }
    } else if (P(jj, 5) > 0.0) {
      eps = 0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(P(jj, 5))) *
                      std::sqrt(1.0 - ipow<2>(P(jj, 5)) / ipow<2>(P(jj, 4)));
      if ((zd(jj) < eps) || (zd(jj) > 1.0 - eps)) {
        if (K(jj, 2) == 21) {
          if (qqbard(jj)) {
            zd(jj) = generatez(0.0, 0.0, eps, "GQ");
          } else {
            zd(jj) = generatez(0.0, 0.0, eps, "GG");
          }
        } else {
          zd(jj) = generatez(0.0, 0.0, eps, "QQ");
        }
      }
      pt2 = ipow<2>(zd(jj)) * ipow<2>(P(jj, 4)) -
            ipow<2>(2.0 * zd(jj) * ipow<2>(P(jj, 4)) - ipow<2>(P(jj, 5))) /
                (4.0 * (ipow<2>(P(jj, 4)) - ipow<2>(P(jj, 5))));
      if (pt2 < -1.e-4) {
      }
    }
  }
  // C--end special treatment
}

void hfcheck(const char* stage) {
  using namespace py::pyjets;
  using namespace cb::angor;
  static const bool dbg = std::getenv("JEWEL_HFDEBUG") != nullptr;
  if (!dbg || !cb::hf::on) return;
  for (int j = 1; j <= N; ++j) {
    if (K(j, 1) != 1 && K(j, 1) != 2) continue;
    if (zd(j) <= 0.0) continue;
    const int hfl = (K(j, 2) == 21 && qqbard(j) && isheavy(cb::hf::qqbarfl(j))) ? cb::hf::qqbarfl(j) : 0;
    if (!isheavy(K(j, 2)) && hfl == 0) continue;
    const double m = hfl ? hfmass(hfl) : hfmass(K(j, 2));
    const double t = virt2(j);
    double zlo, zhi;
    const bool open = hfl ? hfzlimits(t, 0.0, P(j, 4), m, m, zlo, zhi) : hfzlimits(t, m, P(j, 4), m, 0.0, zlo, zhi);
    if (!open || zd(j) < zlo || zd(j) > zhi)
      io::log(" HFCHECK {} line {} kf={} K1={} K3={} hfl={} M={} E={} z={} open={} zlo={} zhi={}", stage, j,
              K(j, 2), K(j, 1), K(j, 3), hfl, P(j, 5), P(j, 4), zd(j), open, open ? zlo : 0.0, open ? zhi : 0.0);
  }
}

// ***********************************************************************
// ***	  subroutine findspeclines
// ***********************************************************************
void findspeclines() {
  using namespace py::pyjets;
  using namespace cb::evrecord;
  using namespace cb::speclines;

  // c--preparatory steps
  if ((collider == "PPZJ") || (collider == "PPZQ") || (collider == "PPZG") ||
      (collider == "PPWJ") || (collider == "PPWQ") || (collider == "PPWG")) {
    lme1orig = 7;
    lme2orig = 8;
    if (std::abs(K(7, 2)) > 21) {
      lv = 7;
    } else {
      lv = 8;
    }
  } else {
    lme1orig = offset - 1;
    lme2orig = offset;
  }
  for (int ipart = offset + 1, ipart_end = offset + npart; ipart <= ipart_end; ++ipart) {
    // C--find decay leptons in V+jet events
    if ((collider == "PPZJ") || (collider == "PPZQ") || (collider == "PPZG") ||
        (collider == "PPWJ") || (collider == "PPWQ") || (collider == "PPWG")) {
      if (K(ipart, 3) == offset - 1) llep1 = ipart;
      if (K(ipart, 3) == offset) llep2 = ipart;
    }
    if (K(ipart, 3) == lme1orig) {
      lme1 = ipart;
      if (K(ipart, 2) == 21) {
        type1 = "GC";
      } else {
        type1 = "QQ";
      }
    } else if (K(ipart, 3) == lme2orig) {
      lme2 = ipart;
      if (K(ipart, 2) == 21) {
        type2 = "GC";
      } else {
        type2 = "QQ";
      }
    }
  }
}

// ***********************************************************************
// ***	  subroutine setcolour
// ***********************************************************************
void setcolour() {
  using namespace py::pyjets;
  using namespace cb::evrecord;
  using namespace cb::colour;
  using namespace cb::speclines;

  // C--local variables
  // NOTE(port): possibly uninitialised in Fortran (firsttrip is read in the
  // middle/end-of-singlet branches before any beginning-of-singlet assignment
  // if the record does not start with a singlet beginning).
  bool firsttrip = false;

  for (int ipart = offset + 1, ipart_end = offset + npart; ipart <= ipart_end; ++ipart) {
    if (K(ipart, 3) != (lme1orig) && (K(ipart, 3) != lme2orig)) {
      trip(ipart) = 0;
      anti(ipart) = 0;
      // C--assign colour indices
    }
    if (K(ipart, 1) == 2) {
      if (K(ipart - 1, 1) == 2) {
        // C--in middle of colour singlet
        if (firsttrip) {
          trip(ipart) = colmax + 1;
          anti(ipart) = trip(ipart - 1);
        } else {
          trip(ipart) = anti(ipart - 1);
          anti(ipart) = colmax + 1;
        }
        colmax = colmax + 1;
      } else {
        // C--beginning of colour singlet
        if (((std::abs(K(ipart, 2)) < 10) && (K(ipart, 2) > 0)) ||
            (isdiquark(K(ipart, 2)) && (K(ipart, 2) < 0))) {
          trip(ipart) = colmax + 1;
          anti(ipart) = 0;
          firsttrip = true;
        } else {
          trip(ipart) = 0;
          anti(ipart) = colmax + 1;
          firsttrip = false;
        }
        colmax = colmax + 1;
      }
    }
    if ((!ishadron(K(ipart, 2))) && (K(ipart, 1) == 1)) {
      // C--end of colour singlet
      if (firsttrip) {
        trip(ipart) = 0;
        anti(ipart) = trip(ipart - 1);
      } else {
        trip(ipart) = anti(ipart - 1);
        anti(ipart) = 0;
      }
    }
  }
}

// ***********************************************************************
// ***	  subroutine findvirtualities
// ***********************************************************************
void findvirtualities(const int& line1, const int& line2, const double& x0, const double& y0,
                      double& z1, double& z2, const bool& mpiscat) {
  using namespace py::pyjets;
  using namespace py::pypars;
  using namespace cb::time;
  using namespace cb::angor;
  using namespace cb::evrecord;
  using namespace cb::param;
  using namespace cb::ftimefac;
  using namespace cb::speclines;

  // C--local variables
  int l1 = 0, l2 = 0;
  double qmax1 = 0, qmax2 = 0, emax = 0, theta1 = 0, theta2 = 0, etot = 0, m1 = 0, m2 = 0,
         p21 = 0, p22 = 0, weight = 0, q1 = 0, q2 = 0, enew1 = 0, enew2 = 0, eps1 = 0, eps2 = 0,
         pold = 0, lambda = 0, thetam1 = 0, thetam2 = 0, theta = 0;
  [[maybe_unused]] double en1 = 0, en2 = 0;
  double x1 = 0, x2 = 0, x3 = 0, meweight = 0, psweight = 0;
  F1<double, 3> beta;
  bool which1 = false, which2 = false;
  // heavy flavour: quark masses of the two partons (0 for light partons / HEAVYFLAVOUR F),
  // g -> q qbar flavours; q1, q2 stay virtualities, m1, m2 are invariant masses
  double mq1 = 0, mq2 = 0;
  int fl1 = 0, fl2 = 0;

  // C--find virtualities and adapt four-vectors
  l1 = lme1;
  l2 = lme2;
  if (mpiscat) {
    thetam1 = -1.0;
    thetam2 = -1.0;
  } else {
    thetam1 = getthetam(lme1);
    thetam2 = getthetam(lme2);
  }
  if (((collider == "PPZJ") || (collider == "PPZQ") || (collider == "PPZG") ||
       (collider == "PPWJ") || (collider == "PPWQ") || (collider == "PPWG")) &&
      (!mpiscat)) {
    if (std::abs(K(l1, 2)) > 21) {
      qmax1 = 0.0;
      qmax2 = std::sqrt(PARI(18) + ipow<2>(P(l1, 5)));
    } else {
      qmax1 = std::sqrt(PARI(18) + ipow<2>(P(l2, 5)));
      qmax2 = 0.0;
    }
    emax = P(l1, 4) + P(l2, 4);
    theta1 = -1.0;
    theta2 = -1.0;
  } else if ((collider == "PPJJ") || (collider == "PPYJ") || (collider == "PPYQ") ||
             (collider == "PPYG") || mpiscat) {
    if (mpiscat) {
      l1 = line1;
      l2 = line2;
    }
    if ((K(l1, 1) == 4) || (mpiscat && (!mpifsr))) {
      qmax1 = 0.0;
    } else {
      if (mpiscat && isparton(K(l1, 2))) {
        qmax1 = py::pyp(l1, 10) * std::exp(r4(0.3f) * std::abs(py::pyp(l1, 17) - py::pyp(l2, 17)) / 2.0) / 2.0;
      } else {
        qmax1 = PARI(17);
      }
    }
    if ((K(l2, 1) == 4) || (mpiscat && (!mpifsr))) {
      qmax2 = 0.0;
    } else {
      if (mpiscat && isparton(K(l2, 2))) {
        qmax2 = py::pyp(l2, 10) * std::exp(r4(0.3f) * std::abs(py::pyp(l1, 17) - py::pyp(l2, 17)) / 2.0) / 2.0;
      } else {
        qmax2 = PARI(17);
      }
    }
    emax = P(l1, 4) + P(l2, 4);
    theta1 = -1.0;
    theta2 = -1.0;
  }
  en1 = P(l1, 4);
  en2 = P(l2, 4);
  beta(1) = (P(l1, 1) + P(l2, 1)) / (P(l1, 4) + P(l2, 4));
  beta(2) = (P(l1, 2) + P(l2, 2)) / (P(l1, 4) + P(l2, 4));
  beta(3) = (P(l1, 3) + P(l2, 3)) / (P(l1, 4) + P(l2, 4));
  py::pyrobo(l1, l1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
  py::pyrobo(l2, l2, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
  etot = P(l1, 4) + P(l2, 4);
  py::pyrobo(l1, l1, 0.0, 0.0, beta(1), beta(2), beta(3));
  py::pyrobo(l2, l2, 0.0, 0.0, beta(1), beta(2), beta(3));
  if (collider == "EEJJ") {
    qmax1 = etot;
    qmax2 = etot;
    emax = P(l1, 4) + P(l2, 4);
    theta1 = -1.0;
    theta2 = -1.0;
  }
  // C--find virtuality
  mq1 = hfmass(K(l1, 2));
  mq2 = hfmass(K(l2, 2));
  q1 = getmass(0.0, qmax1, theta1, emax, type1, emax, false, z1, which1, mq1, &fl1);
  q2 = getmass(0.0, qmax2, theta2, emax, type2, emax, false, z2, which2, mq2, &fl2);
l182:
  if (std::abs(K(l1, 2)) > 21) {
    m1 = P(l1, 5);
  } else {
    m1 = invmass(q1, mq1);
  }
  if (std::abs(K(l2, 2)) > 21) {
    m2 = P(l2, 5);
  } else {
    m2 = invmass(q2, mq2);
  }
  enew1 = etot / 2.0 + (ipow<2>(m1) - ipow<2>(m2)) / (2.0 * etot);
  enew2 = etot / 2.0 - (ipow<2>(m1) - ipow<2>(m2)) / (2.0 * etot);
  p21 = ipow<2>(etot / 2.0 + (ipow<2>(m1) - ipow<2>(m2)) / (2.0 * etot)) - ipow<2>(m1);
  p22 = ipow<2>(etot / 2.0 - (ipow<2>(m1) - ipow<2>(m2)) / (2.0 * etot)) - ipow<2>(m2);
  weight = 1.0;
  if ((py::pyr(0) > weight) || (p21 < 0.0) || (p22 < 0.0) || (enew1 < 0.0) || (enew2 < 0.0)) {
    if (q1 > q2) {
      q1 = getmass(0.0, q1, theta1, emax, type1, emax, false, z1, which1, mq1, &fl1);
    } else {
      q2 = getmass(0.0, q2, theta2, emax, type2, emax, false, z2, which2, mq2, &fl2);
    }
    goto l182;
  }
  py::pyrobo(l1, l1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
  py::pyrobo(l2, l2, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
  pold = py::pyp(l1, 8);
  P(l1, 1) = P(l1, 1) * std::sqrt(p21) / pold;
  P(l1, 2) = P(l1, 2) * std::sqrt(p21) / pold;
  P(l1, 3) = P(l1, 3) * std::sqrt(p21) / pold;
  P(l1, 4) = enew1;
  P(l1, 5) = m1;
  pold = py::pyp(l2, 8);
  P(l2, 1) = P(l2, 1) * std::sqrt(p22) / pold;
  P(l2, 2) = P(l2, 2) * std::sqrt(p22) / pold;
  P(l2, 3) = P(l2, 3) * std::sqrt(p22) / pold;
  P(l2, 4) = enew2;
  P(l2, 5) = m2;
  py::pyrobo(l1, l1, 0.0, 0.0, beta(1), beta(2), beta(3));
  py::pyrobo(l2, l2, 0.0, 0.0, beta(1), beta(2), beta(3));
  // C--correct for overestimated energy and impose angular ordering
  if (q1 > 0.0) {
    eps1 = 0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(q1)) *
                     std::sqrt(1.0 - ipow<2>(q1) / ipow<2>(P(l1, 4)));
    if (hfzout(q1, mq1, P(l1, 4), z1, eps1, which1 ? fl1 : 0)) {
      q1 = getmass(0.0, q1, theta1, emax, type1, emax, false, z1, which1, mq1, &fl1);
      goto l182;
    }
    if (aofirstsplit && (thetam1 > 0.0)) {
      theta = (mq1 > 0.0) ? hftheta(q1 * q1, mq1, z1, P(l1, 4))
                           : q1 / (std::sqrt(z1 * (1.0 - z1)) * P(l1, 4));
      if (theta > thetam1) {
        q1 = getmass(0.0, q1, theta1, emax, type1, emax, false, z1, which1, mq1, &fl1);
        goto l182;
      }
    }
  }
  if (q2 > 0.0) {
    eps2 = 0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(q2)) *
                     std::sqrt(1.0 - ipow<2>(q2) / ipow<2>(P(l2, 4)));
    if (hfzout(q2, mq2, P(l2, 4), z2, eps2, which2 ? fl2 : 0)) {
      q2 = getmass(0.0, q2, theta2, emax, type2, emax, false, z2, which2, mq2, &fl2);
      goto l182;
    }
    if (aofirstsplit && (thetam2 > 0.0)) {
      theta = (mq2 > 0.0) ? hftheta(q2 * q2, mq2, z2, P(l2, 4))
                           : q2 / (std::sqrt(z2 * (1.0 - z2)) * P(l2, 4));
      if (theta > thetam2) {
        q2 = getmass(0.0, q2, theta2, emax, type2, emax, false, z2, which2, mq2, &fl2);
        goto l182;
      }
    }
  }

  // C--correct to ME for first parton
  if (collider == "EEJJ") {
    beta(1) = (P(l1, 1) + P(l2, 1)) / (P(l1, 4) + P(l2, 4));
    beta(2) = (P(l1, 2) + P(l2, 2)) / (P(l1, 4) + P(l2, 4));
    beta(3) = (P(l1, 3) + P(l2, 3)) / (P(l1, 4) + P(l2, 4));
    py::pyrobo(l1, l1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
    py::pyrobo(l2, l2, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
    if (q1 > 0.0) {
      // C--generate z value
      x1 = z1 * (ipow<2>(etot) + ipow<2>(q1)) / ipow<2>(etot);
      x2 = (ipow<2>(etot) - ipow<2>(q1)) / ipow<2>(etot);
      x3 = (1.0 - z1) * (ipow<2>(etot) + ipow<2>(q1)) / ipow<2>(etot);
      if (mq1 > 0.0) {
        // heavy quark: massive ME correction (x1: emitting quark, x2: the other one)
        weight = hfmeweight(x1, x2, mq1 / etot, K(l1, 2));
      } else {
        psweight = (1.0 - x1) * (1.0 + ipow<2>(x1 / (2.0 - x2))) / x3 +
                   (1.0 - x2) * (1.0 + ipow<2>(x2 / (2.0 - x1))) / x3;
        meweight = ipow<2>(x1) + ipow<2>(x2);
        weight = meweight / psweight;
      }
      if (py::pyr(0) > weight) {
        // 184
        q1 = getmass(0.0, q1, theta1, emax, type1, emax, false, z1, which1, mq1, &fl1);
      }
    }
    // C--correct to ME for second parton
    if (q2 > 0.0) {
      // C--generate z value
      x1 = (ipow<2>(etot) - ipow<2>(q2)) / ipow<2>(etot);
      x2 = z2 * (ipow<2>(etot) + ipow<2>(q2)) / ipow<2>(etot);
      x3 = (1.0 - z2) * (ipow<2>(etot) + ipow<2>(q2)) / ipow<2>(etot);
      if (mq2 > 0.0) {
        weight = hfmeweight(x2, x1, mq2 / etot, K(l2, 2));
      } else {
        psweight = (1.0 - x1) * (1.0 + ipow<2>(x1 / (2.0 - x2))) / x3 +
                   (1.0 - x2) * (1.0 + ipow<2>(x2 / (2.0 - x1))) / x3;
        meweight = ipow<2>(x1) + ipow<2>(x2);
        weight = meweight / psweight;
      }
      if (py::pyr(0) > weight) {
        // 185
        q2 = getmass(0.0, q2, theta2, emax, type2, emax, false, z2, which2, mq2, &fl2);
      }
    }
  l186:
    m1 = invmass(q1, mq1);
    m2 = invmass(q2, mq2);
    enew1 = etot / 2.0 + (ipow<2>(m1) - ipow<2>(m2)) / (2.0 * etot);
    enew2 = etot / 2.0 - (ipow<2>(m1) - ipow<2>(m2)) / (2.0 * etot);
    p21 = ipow<2>(etot / 2.0 + (ipow<2>(m1) - ipow<2>(m2)) / (2.0 * etot)) - ipow<2>(m1);
    p22 = ipow<2>(etot / 2.0 - (ipow<2>(m1) - ipow<2>(m2)) / (2.0 * etot)) - ipow<2>(m2);
    pold = py::pyp(l1, 8);
    P(l1, 1) = P(l1, 1) * std::sqrt(p21) / pold;
    P(l1, 2) = P(l1, 2) * std::sqrt(p21) / pold;
    P(l1, 3) = P(l1, 3) * std::sqrt(p21) / pold;
    P(l1, 4) = enew1;
    P(l1, 5) = m1;
    pold = py::pyp(l2, 8);
    P(l2, 1) = P(l2, 1) * std::sqrt(p22) / pold;
    P(l2, 2) = P(l2, 2) * std::sqrt(p22) / pold;
    P(l2, 3) = P(l2, 3) * std::sqrt(p22) / pold;
    P(l2, 4) = enew2;
    P(l2, 5) = m2;
    py::pyrobo(l1, l1, 0.0, 0.0, beta(1), beta(2), beta(3));
    py::pyrobo(l2, l2, 0.0, 0.0, beta(1), beta(2), beta(3));
    // C--correct for overestimated energy
    if (q1 > 0.0) {
      eps1 = 0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(q1)) *
                       std::sqrt(1.0 - ipow<2>(q1) / ipow<2>(P(l1, 4)));
      if (hfzout(q1, mq1, P(l1, 4), z1, eps1, which1 ? fl1 : 0)) {
        q1 = getmass(0.0, q1, theta1, emax, type1, emax, false, z1, which1, mq1, &fl1);
        py::pyrobo(l1, l1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
        py::pyrobo(l2, l2, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
        goto l186;
      }
    }
    if (q2 > 0.0) {
      eps2 = 0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(q2)) *
                       std::sqrt(1.0 - ipow<2>(q2) / ipow<2>(P(l2, 4)));
      if (hfzout(q2, mq2, P(l2, 4), z2, eps2, which2 ? fl2 : 0)) {
        q2 = getmass(0.0, q2, theta2, emax, type2, emax, false, z2, which2, mq2, &fl2);
        py::pyrobo(l1, l1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
        py::pyrobo(l2, l2, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
        goto l186;
      }
    }
  }

  // C--transfer recoil to decay leptons in V+jet
  if (((collider == "PPZJ") || (collider == "PPZQ") || (collider == "PPZG") ||
       (collider == "PPWJ") || (collider == "PPWQ") || (collider == "PPWG")) &&
      (!mpiscat)) {
    beta(1) = P(lv, 1) / P(lv, 4);
    beta(2) = P(lv, 2) / P(lv, 4);
    beta(3) = P(lv, 3) / P(lv, 4);
    py::pyrobo(llep1, llep1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
    py::pyrobo(llep2, llep2, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
    if (std::abs(K(l1, 2)) > 21) {
      beta(1) = P(l1, 1) / P(l1, 4);
      beta(2) = P(l1, 2) / P(l1, 4);
      beta(3) = P(l1, 3) / P(l1, 4);
    } else {
      beta(1) = P(l2, 1) / P(l2, 4);
      beta(2) = P(l2, 2) / P(l2, 4);
      beta(3) = P(l2, 3) / P(l2, 4);
    }
    py::pyrobo(llep1, llep1, 0.0, 0.0, beta(1), beta(2), beta(3));
    py::pyrobo(llep2, llep2, 0.0, 0.0, beta(1), beta(2), beta(3));
  }

  if (z1 > 0.0) {
    thetad(l1) = (mq1 > 0.0) ? hftheta(q1 * q1, mq1, z1, P(l1, 4))
                             : P(l1, 5) / (std::sqrt(z1 * (1.0 - z1)) * P(l1, 4));
  } else {
    thetad(l1) = 0.0;
  }
  if (z2 > 0.0) {
    thetad(l2) = (mq2 > 0.0) ? hftheta(q2 * q2, mq2, z2, P(l2, 4))
                             : P(l2, 5) / (std::sqrt(z2 * (1.0 - z2)) * P(l2, 4));
  } else {
    thetad(l2) = 0.0;
  }
  zd(l1) = z1;
  zd(l2) = z2;
  qqbard(l1) = which1;
  qqbard(l2) = which2;
  if (cb::hf::on) {
    cb::hf::qqbarfl(l1) = fl1;
    cb::hf::qqbarfl(l2) = fl2;
  }
  aoreq(l1) = false;
  aoreq(l2) = false;

  mv(l1, 1) = x0;
  mv(l1, 2) = y0;
  mv(l1, 3) = 0.0;
  mv(l1, 4) = 0.0;
  if (offshell(l1)) {
    lambda = 1.0 / (ftfac * P(l1, 4) * r4(0.2f) / ipow<2>(q1));
    mv(l1, 5) = -std::log(1.0 - py::pyr(0)) / lambda;
  } else {
    mv(l1, 5) = ltime;
  }

  mv(l2, 1) = x0;
  mv(l2, 2) = y0;
  mv(l2, 3) = 0.0;
  mv(l2, 4) = 0.0;
  if (offshell(l2)) {
    lambda = 1.0 / (ftfac * P(l2, 4) * r4(0.2f) / ipow<2>(q2));
    mv(l2, 5) = -std::log(1.0 - py::pyr(0)) / lambda;
  } else {
    mv(l2, 5) = ltime;
  }
}

// ***********************************************************************
// ***	  subroutine findmpivirtualities
// ***********************************************************************
void findmpivirtualities(const double& x0, const double& y0) {
  using namespace py::pyjets;
  using namespace cb::evrecord;
  using namespace cb::param;
  using namespace cb::pythiaparams;
  using namespace cb::speclines;

  // C--local variables
  double z1 = 0, z2 = 0;

  if (!mpi) return;
  // !TODO baryons can appear in beam remnants and have k(line,3) = 0 -> have to take care of this when finding mpi scatterings!
  for (int i = offset + 1, i_end = offset + npart; i <= i_end; ++i) {
    // C--check whether this is an mpi parton
    if (K(i, 3) == 0) {
      // C--find partner parton from same scattering
      for (int j = i + 1, j_end = offset + npart; j <= j_end; ++j) {
        if (std::abs(py::pyp(i, 10) - py::pyp(j, 10)) < 1.e-4) {
          findvirtualities(i, j, x0, y0, z1, z2, true);
          break;  // goto 200
        }
      }
    }
  }
}

// ***********************************************************************
// ***	  function getthetam
// ***********************************************************************
double getthetam(const int& line) {
  using namespace py::pyjets;
  using namespace cb::colour;
  using namespace cb::hepmcid;

  // C--local variables
  int flavcode = 0;
  // NOTE(port): possibly uninitialised in Fortran (l1/l2 are only set if a
  // colour partner is found in loops 100-102).
  int l1 = 0, l2 = 0;
  double theta1 = 0, theta2 = 0, p1dotp2 = 0;

  flavcode = K(line, 2);
  double getthetam = -1.0;
  if ((std::abs(flavcode) > 6) && (flavcode != 21)) return getthetam;
  if (flavcode == 21) {
    for (int i = 1, i_end = N; i <= i_end; ++i) {
      if (K(i, 1) < 10) {
        if (anti(i) == trip(line)) l1 = i;
        if (trip(i) == anti(line)) l2 = i;
      }
    }
    if (K(l1, 3) > 2) {
      p1dotp2 = P(line, 1) * P(l1, 1) + P(line, 2) * P(l1, 2) + P(line, 3) * P(l1, 3);
      theta1 = std::acos(p1dotp2 / (py::pyp(line, 8) * py::pyp(l1, 8)));
    } else {
      theta1 = 1.0e6;
    }
    if (K(l2, 3) > 2) {
      p1dotp2 = P(line, 1) * P(l2, 1) + P(line, 2) * P(l2, 2) + P(line, 3) * P(l2, 3);
      theta2 = std::acos(p1dotp2 / (py::pyp(line, 8) * py::pyp(l2, 8)));
    } else {
      theta2 = 1.0e6;
    }
    getthetam = std::min(theta1, theta2);
    if (getthetam > 1.0e5) getthetam = -1.0;
  } else if (flavcode < 0) {
    for (int i = 1, i_end = N; i <= i_end; ++i) {
      if (K(i, 1) < 10) {
        if (trip(i) == anti(line)) l1 = i;
      }
    }
    if (K(l1, 3) > 2) {
      p1dotp2 = P(line, 1) * P(l1, 1) + P(line, 2) * P(l1, 2) + P(line, 3) * P(l1, 3);
      getthetam = std::acos(p1dotp2 / (py::pyp(line, 8) * py::pyp(l1, 8)));
    }
  } else {
    for (int i = 1, i_end = N; i <= i_end; ++i) {
      if (K(i, 1) < 10) {
        if (anti(i) == trip(line)) l1 = i;
      }
    }
    if (K(l1, 3) > 2) {
      p1dotp2 = P(line, 1) * P(l1, 1) + P(line, 2) * P(l1, 2) + P(line, 3) * P(l1, 3);
      getthetam = std::acos(p1dotp2 / (py::pyp(line, 8) * py::pyp(l1, 8)));
    }
  }
  return getthetam;
}

// ***********************************************************************
// ***	  subroutine makestrings
// ***********************************************************************
void makestrings(const int& which) {
  using namespace cb::hepmcid;
  if (which == 0) {
    makestrings_vac();
  } else if (which == 1) {
    makestrings_minl();
  } else {
    io::log(" error: unknown hadronisation type in MAKESTRINGS");
  }
}

// ***********************************************************************
// ***	  subroutine makestrings_vac
// ***********************************************************************
void makestrings_vac() {
  using namespace cb::hepmcid;
  using namespace py::pyjets;
  using namespace cb::param;
  using namespace cb::colour;
  using namespace cb::disc;

  // C--local variables
  int nold = 0, i = 0, j = 0, lquark = 0, lloose = 0, nold1 = 0;
  [[maybe_unused]] int lmatch = 0;
  [[maybe_unused]] bool roomleft = false;
  // DATA EADDEND/10.d0/ (unused)

  i = 0;
  if (compress) roomleft = compressevent(i);
  nold1 = N;
  // C--remove all active lines that are leptons, gammas, hadrons etc.
  for (int i_end = (i = 1, nold1); i <= i_end; ++i) {
    if ((K(i, 1) == 4) && (trip(i) == 0) && (anti(i) == 0)) {
      // C--copy line to end of event record
      N = N + 1;
      if (N > 22990) {
        io::log(" event too long for event record");
        discard = true;
        return;
      }
      K(N, 1) = 11;
      K(N, 2) = K(i, 2);
      K(N, 3) = i;
      K(N, 4) = 0;
      K(N, 5) = 0;
      P(N, 1) = P(i, 1);
      P(N, 2) = P(i, 2);
      P(N, 3) = P(i, 3);
      P(N, 4) = P(i, 4);
      P(N, 5) = P(i, 5);
      K(i, 1) = 17;
      K(i, 4) = N;
      K(i, 5) = N;
      trip(N) = trip(i);
      anti(N) = anti(i);
    }
  }
  nold = N;
  // C--first do strings with existing (anti)triplets
  // C--find string end (=quark or antiquark)
l43:
  lquark = 0;
  for (int i_end = (i = 1, nold); i <= i_end; ++i) {
    if ((K(i, 1) == 11) || (K(i, 1) == 12) || (K(i, 1) == 13) || (K(i, 1) == 14)) K(i, 1) = 17;
    if (((K(i, 1) == 1) || (K(i, 1) == 3) || (K(i, 1) == 4) || (K(i, 1) == 5)) &&
        ((K(i, 2) < 6) || isdiquark(K(i, 2)))) {
      lquark = i;
      goto l41;
    }
  }
  goto l50;
l41:
  // C--copy string end to end of event record
  N = N + 1;
  if (N > 22990) {
    io::log(" event too long for event record");
    discard = true;
    return;
  }
  K(N, 1) = 2;
  K(N, 2) = K(lquark, 2);
  K(N, 3) = lquark;
  K(N, 4) = 0;
  K(N, 5) = 0;
  P(N, 1) = P(lquark, 1);
  P(N, 2) = P(lquark, 2);
  P(N, 3) = P(lquark, 3);
  P(N, 4) = P(lquark, 4);
  P(N, 5) = P(lquark, 5);
  K(lquark, 1) = 16;
  K(lquark, 4) = N;
  K(lquark, 5) = N;
  trip(N) = trip(lquark);
  anti(N) = anti(lquark);
  // C--append matching colour partner
  lmatch = 0;
  for (j = 1; j <= 10000000; ++j) {
    for (int i_end = (i = 1, nold); i <= i_end; ++i) {
      if (((K(i, 1) == 1) || (K(i, 1) == 3) || (K(i, 1) == 4) || (K(i, 1) == 5)) &&
          (((trip(i) == anti(N)) && (trip(i) != 0)) || ((anti(i) == trip(N)) && (anti(i) != 0)))) {
        N = N + 1;
        if (N > 22990) {
          io::log(" event too long for event record");
          discard = true;
          return;
        }
        K(N, 2) = K(i, 2);
        K(N, 3) = i;
        K(N, 4) = 0;
        K(N, 5) = 0;
        P(N, 1) = P(i, 1);
        P(N, 2) = P(i, 2);
        P(N, 3) = P(i, 3);
        P(N, 4) = P(i, 4);
        P(N, 5) = P(i, 5);
        trip(N) = trip(i);
        anti(N) = anti(i);
        K(i, 1) = 16;
        K(i, 4) = N;
        K(i, 5) = N;
        if (K(i, 2) == 21) {
          K(N, 1) = 2;
          goto l44;
        } else {
          K(N, 1) = 1;
          goto l43;
        }
      }
    }
    // C--no matching colour partner found
    io::log(
        " Error in MAKESTRINGS_VAC: failed to reconstruct colour singlet system, will discard event {} {} {} {}",
        N, trip(N), anti(N), i);
    pevrec(2, true);
    discard = true;
    return;
  l44:;
  }
  // C--now take care of purely gluonic remainder system
  // C-----------------------------------------
  // C--find gluon where anti-triplet is not matched
l50:
  lloose = 0;
  for (int i_end = (i = 1, nold); i <= i_end; ++i) {
    if (((K(i, 1) == 1) || (K(i, 1) == 3) || (K(i, 1) == 4) || (K(i, 1) == 5))) {
      for (int j_end = (j = 1, nold); j <= j_end; ++j) {
        // NOTE(port): the Fortran tests K(I,1) (not K(J,1)) here; kept as is.
        if (((K(i, 1) == 1) || (K(i, 1) == 3) || (K(i, 1) == 4) || (K(i, 1) == 5))) {
          if (anti(i) == trip(j)) goto l45;
        }
      }
      lloose = i;
      goto l47;
    }
  l45:;
  }
  goto l51;
l47:
  // C--generate artificial triplet end
  io::log(" Error in MAKESTRINGS_VAC: failed to reconstruct colour singlet system, will discard event");
  discard = true;
  return;
  // NOTE(port): the code below (copy loose gluon, loop 48/49) is unreachable in
  // the Fortran original (it follows an unconditional RETURN); translated for completeness.
  // C--copy loose gluon to end of event record
  N = N + 1;
  if (N > 22990) {
    io::log(" event too long for event record");
    discard = true;
    return;
  }
  K(N, 1) = 2;
  K(N, 2) = K(lloose, 2);
  K(N, 3) = lloose;
  K(N, 4) = 0;
  K(N, 5) = 0;
  P(N, 1) = P(lloose, 1);
  P(N, 2) = P(lloose, 2);
  P(N, 3) = P(lloose, 3);
  P(N, 4) = P(lloose, 4);
  P(N, 5) = P(lloose, 5);
  K(lloose, 1) = 16;
  K(lloose, 4) = N;
  K(lloose, 5) = N;
  trip(N) = trip(lloose);
  anti(N) = anti(lloose);
  // C--append matching colour partner
  lmatch = 0;
  for (j = 1; j <= 10000000; ++j) {
    for (int i_end = (i = 1, nold); i <= i_end; ++i) {
      if (((K(i, 1) == 1) || (K(i, 1) == 3) || (K(i, 1) == 4) || (K(i, 1) == 5)) &&
          (anti(i) == trip(N))) {
        N = N + 1;
        if (N > 22990) {
          io::log(" event too long for event record");
          discard = true;
          return;
        }
        K(N, 2) = K(i, 2);
        K(N, 3) = i;
        K(N, 4) = 0;
        K(N, 5) = 0;
        P(N, 1) = P(i, 1);
        P(N, 2) = P(i, 2);
        P(N, 3) = P(i, 3);
        P(N, 4) = P(i, 4);
        P(N, 5) = P(i, 5);
        trip(N) = trip(i);
        anti(N) = anti(i);
        K(i, 1) = 16;
        K(i, 4) = N;
        K(i, 5) = N;
        K(N, 1) = 2;
        goto l48;
      }
    }
    // C--no matching colour partner found, add artificial end point
    io::log(" Error in MAKESTRINGS_VAC: failed to reconstruct colour singlet system, will discard event");
    discard = true;
    return;
  l48:;
  }
l51:
  cleanup(nold1);
}

// ***********************************************************************
// ***	  subroutine makestrings_minl
// ***********************************************************************
void makestrings_minl() {
  using namespace py::pyjets;
  using namespace cb::param;
  using namespace cb::colour;
  using namespace cb::hepmcid;
  using namespace cb::disc;
  using namespace cb::rapmax;

  // C--local variables
  int nold = 0, i = 0, lmax = 0, lmin = 0, lend = 0;
  double emax = 0, minv = 0, mmin = 0, z = 0, dir = 0;
  [[maybe_unused]] bool roomleft = false;
  constexpr double mcut = 1.e8;       // DATA MCUT/1.d8/
  constexpr double eaddend = 10.0;    // DATA EADDEND/10.d0/
  constexpr double deltaymax = 4.0;   // data deltaymax/4.d0/

  i = 0;
  if (compress) roomleft = compressevent(i);
  // C--label all active lines that are leptons, gammas, hadrons etc.
  for (int i_end = (i = 1, N); i <= i_end; ++i) {
    if ((K(i, 1) <= 4) && (trip(i) == 0) && (anti(i) == 0)) {
      K(i, 1) = 19;
    }
  }
  nold = N;
  // C--find most energetic unfragmented parton in event
l43:
  emax = 0;
  lmax = 0;
  for (int i_end = (i = 1, nold); i <= i_end; ++i) {
    if ((K(i, 1) == 11) || (K(i, 1) == 12) || (K(i, 1) == 13) || (K(i, 1) == 14)) K(i, 1) = 17;
    if (std::abs(py::pyp(i, 17)) > etamax + 1.0) K(i, 1) = 17;
    if (((K(i, 1) == 1) || (K(i, 1) == 3) || (K(i, 1) == 4) || (K(i, 1) == 5)) && (P(i, 4) > emax)) {
      emax = P(i, 4);
      lmax = i;
    }
  }
  // C--if there is none, we are done
  if (lmax == 0) goto l50;
  // C--check if highest energy parton is (anti)quark or gluon
  if (K(lmax, 2) == 21) {
    // C--split gluon in qqbar pair and store one temporarily in line 1
    // C--make new line in event record for string end
    N = N + 2;
    if (N > 22990) {
      io::log(" event too long for event record");
      discard = true;
      return;
    }
    if ((N - 2) > nold) {
      // NOTE(port): this shift moves K and P only, not trip/anti (original behaviour kept).
      for (int j = nold, j_end = N - 3; j <= j_end; ++j) {
        K(N + nold - j, 1) = K(N + nold - j - 2, 1);
        K(N + nold - j, 2) = K(N + nold - j - 2, 2);
        if (K(N + nold - j - 2, 3) > nold) {
          K(N + nold - j, 3) = K(N + nold - j - 2, 3) + 2;
        } else {
          K(N + nold - j, 3) = K(N + nold - j - 2, 3);
        }
        K(N + nold - j, 4) = 0;
        K(N + nold - j, 5) = 0;
        P(N + nold - j, 1) = P(N + nold - j - 2, 1);
        P(N + nold - j, 2) = P(N + nold - j - 2, 2);
        P(N + nold - j, 3) = P(N + nold - j - 2, 3);
        P(N + nold - j, 4) = P(N + nold - j - 2, 4);
        P(N + nold - j, 5) = P(N + nold - j - 2, 5);
        if (K(N + nold - j - 2, 3) > 0) {
          K(K(N + nold - j - 2, 3), 4) = K(K(N + nold - j - 2, 3), 4) + 2;
          K(K(N + nold - j - 2, 3), 5) = K(K(N + nold - j - 2, 3), 5) + 2;
        }
      }
    }
    nold = nold + 2;
    K(lmax, 1) = 18;
    z = generatez(0.0, 0.0, 1.e-3, "QG");
    if (z > 0.5) {
      K(nold - 1, 2) = 1;
      K(nold, 2) = -1;
    } else {
      z = 1.0 - z;
      K(nold - 1, 2) = -1;
      K(nold, 2) = 1;
    }
    K(nold - 1, 1) = 1;
    K(nold - 1, 3) = lmax;
    K(nold - 1, 4) = 0;
    K(nold - 1, 5) = 0;
    P(nold - 1, 1) = (1.0 - z) * P(lmax, 1);
    P(nold - 1, 2) = (1.0 - z) * P(lmax, 2);
    P(nold - 1, 3) = (1.0 - z) * P(lmax, 3);
    P(nold - 1, 4) = (1.0 - z) * P(lmax, 4);
    P(nold - 1, 5) = P(lmax, 5);
    K(nold, 1) = 1;
    K(nold, 3) = lmax;
    K(nold, 4) = 0;
    K(nold, 5) = 0;
    P(nold, 1) = z * P(lmax, 1);
    P(nold, 2) = z * P(lmax, 2);
    P(nold, 3) = z * P(lmax, 3);
    P(nold, 4) = z * P(lmax, 4);
    P(nold, 5) = P(lmax, 5);
    K(lmax, 1) = 18;
    K(lmax, 4) = nold - 1;
    K(lmax, 5) = nold;
    lmax = nold;
  }
  N = N + 1;
  if (N > 22990) {
    io::log(" event too long for event record");
    discard = true;
    return;
  }
  K(N, 1) = 2;
  K(N, 2) = K(lmax, 2);
  K(N, 3) = lmax;
  K(N, 4) = 0;
  K(N, 5) = 0;
  P(N, 1) = P(lmax, 1);
  P(N, 2) = P(lmax, 2);
  P(N, 3) = P(lmax, 3);
  P(N, 4) = P(lmax, 4);
  P(N, 5) = P(lmax, 5);
  K(lmax, 1) = 16;
  K(lmax, 4) = N;
  K(lmax, 5) = N;
  lend = lmax;
  // C--find closest partner
l42:
  mmin = 1.e10;
  lmin = 0;
  for (int i_end = (i = 1, nold); i <= i_end; ++i) {
    if (((K(i, 1) == 1) || (K(i, 1) == 3) || (K(i, 1) == 4) || (K(i, 1) == 5)) &&
        ((K(i, 2) == 21) ||
         ((static_cast<double>(K(i, 2) * K(lend, 2)) < 0.0) && (K(i, 3) != K(lend, 3)))) &&
        (std::abs(py::pyp(i, 17) - py::pyp(lend, 17)) < deltaymax)) {
      minv = P(i, 4) * P(lmax, 4) - P(i, 1) * P(lmax, 1) - P(i, 2) * P(lmax, 2) - P(i, 3) * P(lmax, 3);
      if ((minv < mmin) && (minv > 0.0) && (minv < mcut)) {
        mmin = minv;
        lmin = i;
      }
    }
  }
  // C--if no closest partner can be found, generate artificial end point for string
  if (lmin == 0) {
    N = N + 1;
    if (N > 22990) {
      io::log(" event too long for event record");
      discard = true;
      return;
    }
    K(N, 1) = 1;
    K(N, 2) = -K(lend, 2);
    K(N, 3) = 0;
    K(N, 4) = 0;
    K(N, 5) = 0;
    P(N, 1) = 0.0;
    P(N, 2) = 0.0;
    if (py::pyr(0) < 0.5) {
      dir = 1.0;
    } else {
      dir = -1.0;
    }
    P(N, 3) = dir * eaddend;
    P(N, 4) = eaddend;
    P(N, 5) = 0.0;
    goto l43;
  } else {
    // C--else build closest partner in string
    N = N + 1;
    if (N > 22990) {
      io::log(" event too long for event record");
      discard = true;
      return;
    }
    K(N, 2) = K(lmin, 2);
    K(N, 3) = lmin;
    K(N, 4) = 0;
    K(N, 5) = 0;
    P(N, 1) = P(lmin, 1);
    P(N, 2) = P(lmin, 2);
    P(N, 3) = P(lmin, 3);
    P(N, 4) = P(lmin, 4);
    P(N, 5) = P(lmin, 5);
    K(lmin, 1) = 16;
    K(lmin, 4) = N;
    K(lmin, 5) = N;
    if (K(lmin, 2) == 21) {
      K(N, 1) = 2;
      lmax = lmin;
      goto l42;
    } else {
      K(N, 1) = 1;
      goto l43;
    }
  }
l50:
  for (int i_end = (i = 1, nold); i <= i_end; ++i) {
    if ((K(i, 1) == 19) && (std::abs(py::pyp(i, 17)) < etamax + 1.0)) {
      // C--copy line to end of event record
      N = N + 1;
      if (N > 22990) {
        io::log(" event too long for event record");
        discard = true;
        return;
      }
      K(N, 1) = 4;
      K(N, 2) = K(i, 2);
      K(N, 3) = i;
      K(N, 4) = 0;
      K(N, 5) = 0;
      P(N, 1) = P(i, 1);
      P(N, 2) = P(i, 2);
      P(N, 3) = P(i, 3);
      P(N, 4) = P(i, 4);
      P(N, 5) = P(i, 5);
      K(i, 1) = 17;
      K(i, 4) = N;
      K(i, 5) = N;
      trip(N) = trip(i);
      anti(N) = anti(i);
    }
  }
  cleanup(nold);
}

// ***********************************************************************
// ***	  subroutine cleanup
// ***********************************************************************
void cleanup(const int& nfirst) {
  using namespace py::pyjets;

  // C--local variables
  int nlast = N;
  for (int i = 1, i_end = nlast - nfirst; i <= i_end; ++i) {
    for (int j = 1; j <= 5; ++j) {
      K(i, j) = K(nfirst + i, j);
      P(i, j) = P(nfirst + i, j);
      V(i, j) = V(nfirst + i, j);
    }
    K(i, 3) = 0;
  }
  N = nlast - nfirst;
}

}  // namespace jewel
