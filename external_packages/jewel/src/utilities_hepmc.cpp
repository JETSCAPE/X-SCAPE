// JEWEL 2.6.0 C++ port -- GETDELTAT, utilities, event-record helpers, HepMC output, logo and time.
// Translated from jewel-2.6.0.f:8423-9697.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "jewel/jfmt.hpp"
#include <string>
#include <string_view>

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

// Fortran Ew.d / E0.d edit descriptor as produced by gfortran:
// [-]0.<d digits>E(+|-)dd, three-digit exponents as (+|-)ddd without the 'E'.
// Rounding is round-to-nearest on the exact binary value (gfortran uses the C
// library printf for the digits, as std::format does). Negative zero keeps its
// sign (gfortran default -fsign-zero). w == 0 gives the minimal width.
std::string fedit(double x, int w, int d) {
  std::string s;
  if (std::isnan(x)) {
    s = "NaN";
  } else if (std::isinf(x)) {
    s = x < 0 ? "-Inf" : "Inf";
  } else {
    const bool neg = std::signbit(x);
    const double ax = std::fabs(x);
    std::string digits;
    int e = 0;
    if (ax == 0.0) {
      digits.assign(static_cast<std::size_t>(d), '0');
      e = 0;
    } else {
      // d significant digits: "D.DDDDDe+XX"
      const std::string t = jewel::jfmt::format("{:.{}e}", ax, d - 1);
      const std::size_t epos = t.find('e');
      digits = t.substr(0, 1);
      if (d > 1) digits += t.substr(2, epos - 2);
      e = std::stoi(t.substr(epos + 1)) + 1;
    }
    const int ae = std::abs(e);
    const char esign = e < 0 ? '-' : '+';
    // gfortran: E0.d always prints a three-digit exponent (E+000); Ew.d uses
    // two digits and drops the letter for |exponent| > 99.
    std::string ex = w == 0     ? jewel::jfmt::format("E{}{:03}", esign, ae)
                     : ae <= 99 ? jewel::jfmt::format("E{}{:02}", esign, ae)
                                : jewel::jfmt::format("{}{:03}", esign, ae);
    s = std::string(neg ? "-" : "") + "0." + digits + ex;
  }
  if (w > 0) {
    if (static_cast<int>(s.size()) > w) return std::string(static_cast<std::size_t>(w), '*');
    return std::string(static_cast<std::size_t>(w) - s.size(), ' ') + s;
  }
  return s;
}
std::string e06(double x) { return fedit(x, 0, 6); }
std::string e08(double x) { return fedit(x, 0, 8); }

void hepmc_line(const std::string& line) {
  std::fputs(line.c_str(), io::hepmc());
  std::fputc('\n', io::hepmc());
}

// 5000 FORMAT(A2,2(I0,' '),3(E0.6,' '),7(I0,' '),E0.6)
void w5000(std::string_view a, int i1, int i2, double r1, double r2, double r3, int j1, int j2, int j3,
           int j4, int j5, int j6, int j7, double r4v) {
  hepmc_line(jewel::jfmt::format("{}{} {} {} {} {} {} {} {} {} {} {} {} {}", a, i1, i2, e06(r1), e06(r2), e06(r3),
                         j1, j2, j3, j4, j5, j6, j7, e06(r4v)));
}
// 5100 FORMAT(A2,E0.6,' ',E0.6)
void w5100(std::string_view a, double r1, double r2) {
  hepmc_line(jewel::jfmt::format("{}{} {}", a, e06(r1), e06(r2)));
}
// 5200 FORMAT(A2,9(I0,' '),3(E0.6,' '),E0.6)
void w5200(std::string_view a, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9,
           double r1, double r2, double r3, double r4v) {
  hepmc_line(jewel::jfmt::format("{}{} {} {} {} {} {} {} {} {} {} {} {} {}", a, i1, i2, i3, i4, i5, i6, i7, i8, i9,
                         e06(r1), e06(r2), e06(r3), e06(r4v)));
}
// 5300 FORMAT(A2,2(I0,' '),5(E0.6,' '),I0,' ',I0)
void w5300(std::string_view a, int i1, int i2, double r1, double r2, double r3, double r4v, double r5,
           int j1, int j2) {
  hepmc_line(jewel::jfmt::format("{}{} {} {} {} {} {} {} {} {}", a, i1, i2, e06(r1), e06(r2), e06(r3), e06(r4v),
                         e06(r5), j1, j2));
}
// 5400 FORMAT(A2,8(I0,' '),I0)
void w5400(std::string_view a, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9) {
  hepmc_line(jewel::jfmt::format("{}{} {} {} {} {} {} {} {} {}", a, i1, i2, i3, i4, i5, i6, i7, i8, i9));
}
// 5500 FORMAT(A2,2(I0,' '),5(E0.8,' '),4(I0,' '),I0)
void w5500(std::string_view a, int i1, int i2, double r1, double r2, double r3, double r4v, double r5,
           int j1, int j2, int j3, int j4, int j5) {
  hepmc_line(jewel::jfmt::format("{}{} {} {} {} {} {} {} {} {} {} {} {}", a, i1, i2, e08(r1), e08(r2), e08(r3),
                         e08(r4v), e08(r5), j1, j2, j3, j4, j5));
}
// WRITE(J,'(A2,I2,A5)')'N ',1,'"0"'
void wnline() { hepmc_line(jewel::jfmt::format("{}{:2}{:>5}", "N ", 1, "\"0\"")); }
// WRITE(J,'(A)')'U GEV MM'
void wuline() { hepmc_line("U GEV MM"); }

}  // namespace

// ***********************************************************************
// ***	  function getdeltat
// ***********************************************************************
bool getdeltat(const int& line, const double& tstart, const double& dtmax1, double& deltat, const bool& elas) {
  using namespace py::pyjets;
  using namespace cb::param;
  using namespace cb::time;
  using namespace cb::errline;
  using namespace cb::scalfac;
  bool getdeltat_ = false;
  double dtmax, sigmamax, neffmax, linvmax, r, toff, xs, ys, zs, ts, xstart, ystart, zstart, weight, neff,
      sigma, tau, sigmamin, neffmin;
  int nnull;
  char ptype;
  bool stopnow;

  // C--initialization
  getdeltat_ = false;
  deltat = 0.0;
  dtmax = dtmax1;
  if (K(line, 2) == 21) {
    ptype = 'G';
  } else {
    ptype = 'Q';
  }

  nnull = 0;
  stopnow = false;

  // C--check for upper bound from plasma lifetime
  if ((tstart + dtmax) >= ltime) dtmax = ltime - tstart;
  if (dtmax < 0.0) return getdeltat_;

  // C--calculate time relative to production of the considered parton
  toff = tstart - mv(line, 4);
  xstart = mv(line, 1) + toff * P(line, 1) / P(line, 4);
  ystart = mv(line, 2) + toff * P(line, 2) / P(line, 4);
  zstart = mv(line, 3) + toff * P(line, 3) / P(line, 4);

  if (elas) {
    sigmamax = getsscat(P(line, 4), P(line, 1), P(line, 2), P(line, 3), P(line, 5), 0.0, ptype, ptype, xstart,
                        ystart, zstart, tstart, 1);
    sigmamin = getsscat(P(line, 4), P(line, 1), P(line, 2), P(line, 3), P(line, 5), 0.0, ptype, ptype, xstart,
                        ystart, zstart, tstart, 2);
  } else {
    // C--calculate upper limit for density*cross section
    sigmamax = getsscat(P(line, 4), P(line, 1), P(line, 2), P(line, 3), P(line, 5), 0.0, ptype, 'C', xstart,
                        ystart, zstart, tstart, 1);
    sigmamin = getsscat(P(line, 4), P(line, 1), P(line, 2), P(line, 3), P(line, 5), 0.0, ptype, 'C', xstart,
                        ystart, zstart, tstart, 2);
  }
  neffmax = medium().getneffmax();
  neffmin = medium().getnatmdmin();
  linvmax = xsecscalefac * 5.0 * std::max(neffmin * sigmamax, neffmax * sigmamin);
  if (linvmax == 0.0) return getdeltat_;

  for (int i = 1; i <= 1000000; ++i) {
    deltat = deltat - std::log(py::pyr(0)) / linvmax;
    xs = xstart + deltat * P(line, 1) / P(line, 4);
    ys = ystart + deltat * P(line, 2) / P(line, 4);
    zs = zstart + deltat * P(line, 3) / P(line, 4);
    ts = tstart + deltat;
    if (ts < zs) {
      tau = -1.0;
    } else {
      tau = std::sqrt(ipow<2>(ts) - ipow<2>(zs));
    }
    neff = medium().getneffprobe(xs, ys, zs, ts, P(line, 1) / P(line, 4), P(line, 2) / P(line, 4),
                                 P(line, 3) / P(line, 4));
    if ((tau > 1.0) && (neff == 0.0)) {
      if (nnull > 4) {
        stopnow = true;
      } else {
        nnull = nnull + 1;
      }
    } else {
      nnull = 0;
    }
    if ((deltat > dtmax) || stopnow) {
      deltat = dtmax;
      return getdeltat_;
    }
    if (neff > 0.0) {
      if (elas) {
        sigma = getsscat(P(line, 4), P(line, 1), P(line, 2), P(line, 3), P(line, 5), 0.0, ptype, ptype, xs, ys,
                         zs, ts, 0) *
                xsecscalefac;
      } else {
        sigma = getsscat(P(line, 4), P(line, 1), P(line, 2), P(line, 3), P(line, 5), 0.0, ptype, 'C', xs, ys,
                         zs, ts, 0) *
                xsecscalefac;
      }
    } else {
      sigma = 0.0;
    }
    weight = 5.0 * neff * sigma / linvmax;
    if (weight > 1.0 + 1e-6) {
      if (line != errl) {
        io::log(" error in GETDELTAT: weight > 1 {} {} {} {}", io::flist(weight),
                io::flist(neff * sigma / (neffmax * sigmamin)), io::flist(neff * sigma / (neffmin * sigmamax)),
                io::flist(P(line, 4)));
        errl = line;
      }
    }
    r = py::pyr(0);
    if (r < weight) {
      getdeltat_ = true;
      return getdeltat_;
    }
  }  // 333
  return getdeltat_;
}

// ***********************************************************************
// ***	  function poissonian
// ***********************************************************************
int poissonian(const double& lambda) {
  int n;
  double disc, p, u, v;
  constexpr double pi = 3.141592653589793;  // data pi/3.141592653589793d0/
  int poissonian_;

  if (lambda > 745.0) {
    u = py::pyr(0);
    v = py::pyr(0);
    poissonian_ = fint(std::sqrt(lambda) * std::sqrt(-2.0 * std::log(u)) * std::cos(2.0 * pi * v) + lambda);
  } else {
    disc = std::exp(-lambda);
    p = 1.0;
    n = 0;
    for (;;) {  // 800
      p = p * py::pyr(0);
      if (p > disc) {
        n = n + 1;
        continue;  // goto 800
      }
      break;
    }
    poissonian_ = n;
  }
  return poissonian_;
}

// ***********************************************************************
// ***	  subroutine makemassless
// ***********************************************************************
void makemassless(const int& l3, const int& l1, const int& l2) {
  using namespace py::pyjets;
  F1<double, 3> beta;
  double etot, e1new, e2new, e3new, p1new, p1new2, m1, m2;

  beta(1) = (P(l1, 1) + P(l2, 1)) / (P(l1, 4) + P(l2, 4));
  beta(2) = (P(l1, 2) + P(l2, 2)) / (P(l1, 4) + P(l2, 4));
  beta(3) = (P(l1, 3) + P(l2, 3)) / (P(l1, 4) + P(l2, 4));
  py::pyrobo(l1, l1, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
  py::pyrobo(l2, l2, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
  py::pyrobo(l3, l3, 0.0, 0.0, -beta(1), -beta(2), -beta(3));
  etot = P(l1, 4) + P(l2, 4) + P(l3, 4);
  e3new = py::pyp(l3, 8);
  m1 = P(l1, 5);
  m2 = P(l2, 5);
  p1new2 = (ipow<2>(ipow<2>(etot - e3new) - ipow<2>(m1) - ipow<2>(m2)) - 4.0 * ipow<2>(m1) * ipow<2>(m2)) /
           (4.0 * ipow<2>(etot - e3new));
  if (p1new2 < 0.0) io::log(" negative mass squared: {}", io::flist(p1new2));
  p1new = std::sqrt(std::fabs(p1new2));
  e1new = std::sqrt(p1new2 + ipow<2>(m1));
  e2new = std::sqrt(p1new2 + ipow<2>(m2));
  P(l1, 1) = P(l1, 1) * p1new / py::pyp(l1, 8);
  P(l1, 2) = P(l1, 2) * p1new / py::pyp(l1, 8);
  P(l1, 3) = P(l1, 3) * p1new / py::pyp(l1, 8);
  P(l1, 4) = e1new;
  P(l2, 1) = P(l2, 1) * p1new / py::pyp(l2, 8);
  P(l2, 2) = P(l2, 2) * p1new / py::pyp(l2, 8);
  P(l2, 3) = P(l2, 3) * p1new / py::pyp(l2, 8);
  P(l2, 4) = e2new;
  P(l3, 4) = e3new;
  P(l3, 5) = 0.0;
  py::pyrobo(l1, l1, 0.0, 0.0, beta(1), beta(2), beta(3));
  py::pyrobo(l2, l2, 0.0, 0.0, beta(1), beta(2), beta(3));
  py::pyrobo(l3, l3, 0.0, 0.0, beta(1), beta(2), beta(3));
}

// ***********************************************************************
// ***	  function ishadron
// ***********************************************************************
bool ishadron(const int& id) {
  bool ishadron_;
  if (std::abs(id) < 100) {
    ishadron_ = false;
  } else {
    // INT(ABS(ID)/10.): REAL*4 division
    if (fint(static_cast<double>(static_cast<float>(std::abs(id)) / 10.0f)) % 10 == 0) {
      ishadron_ = false;
    } else {
      ishadron_ = true;
    }
  }
  return ishadron_;
}

// ***********************************************************************
// ***	  function isdiquark
// ***********************************************************************
bool isdiquark(const int& id) {
  bool isdiquark_;
  if (std::abs(id) < 1000) {
    isdiquark_ = false;
  } else {
    if ((id / 10) % 10 == 0) {
      isdiquark_ = true;
    } else {
      isdiquark_ = false;
    }
  }
  return isdiquark_;
}

// ***********************************************************************
// ***	  function islepton
// ***********************************************************************
bool islepton(const int& id) {
  bool islepton_;
  if ((std::abs(id) == 11) || (std::abs(id) == 13) || (std::abs(id) == 15)) {
    islepton_ = true;
  } else {
    islepton_ = false;
  }
  return islepton_;
}

// ***********************************************************************
// ***	  function isparton
// ***********************************************************************
bool isparton(const int& id) {
  bool isparton_;
  if ((std::abs(id) < 6) || (id == 21) || isdiquark(id)) {
    isparton_ = true;
  } else {
    isparton_ = false;
  }
  return isparton_;
}

// ***********************************************************************
// ***	  function isqorg
// ***********************************************************************
bool isqorg(const int& id) {
  bool isqorg_;
  if ((std::abs(id) < 6) || (id == 21)) {
    isqorg_ = true;
  } else {
    isqorg_ = false;
  }
  return isqorg_;
}

// ***********************************************************************
// ***	  function isprimstring
// ***********************************************************************
bool isprimstring(const int& l) {
  using namespace py::pyjets;
  bool isprimstring_;

  if ((K(l, 2) != 91) && (K(l, 2) != 92)) {
    isprimstring_ = false;
    return isprimstring_;
  }
  isprimstring_ = true;
  if (K(l, 3) == 0) return isprimstring_;
  if (K(K(l, 3), 3) == 0) return isprimstring_;
  if ((K(K(l, 3), 3) == 0) || (isparton(K(K(K(l, 3), 3), 2)))) {
    isprimstring_ = true;
  } else {
    isprimstring_ = false;
  }
  return isprimstring_;
}

// ***********************************************************************
// ***	  function issecstring
// ***********************************************************************
bool issecstring(const int& l) {
  using namespace py::pyjets;
  bool issecstring_;

  if ((K(l, 2) != 91) && (K(l, 2) != 92)) {
    issecstring_ = false;
    return issecstring_;
  }
  if (isprimstring(l)) {
    issecstring_ = false;
    return issecstring_;
  }
  issecstring_ = false;
  if (K(l, 3) == 0) return issecstring_;
  if (K(K(l, 3), 3) == 0) return issecstring_;
  if (isparton(K(K(K(l, 3), 3), 2))) {
    issecstring_ = false;
  } else {
    issecstring_ = true;
  }
  return issecstring_;
}

// ***********************************************************************
// ***	  function isprimhadron
// ***********************************************************************
bool isprimhadron(const int& l) {
  using namespace py::pyjets;
  bool isprimhadron_;

  isprimhadron_ = false;
  if (K(l, 3) == 0) return isprimhadron_;
  if (((K(K(l, 3), 2) == 91) || (K(K(l, 3), 2) == 92)) && isprimstring(K(l, 3)) && (!isparton(K(l, 2)))) {
    isprimhadron_ = true;
  } else {
    isprimhadron_ = false;
  }
  if (K(l, 1) == 17) isprimhadron_ = true;
  return isprimhadron_;
}

// ***********************************************************************
// ***	  subroutine getposition
// ***********************************************************************
void getposition(const int& line, const double& time, double& posx, double& posy, double& posz) {
  using namespace py::pyjets;
  using namespace cb::time;

  posx = mv(line, 1) + (time - mv(line, 4)) * P(line, 1) / std::max(py::pyp(line, 8), P(line, 4));
  posy = mv(line, 2) + (time - mv(line, 4)) * P(line, 2) / std::max(py::pyp(line, 8), P(line, 4));
  posz = mv(line, 3) + (time - mv(line, 4)) * P(line, 3) / std::max(py::pyp(line, 8), P(line, 4));
}

// ***********************************************************************
// ***	  function compressevent
// library mode: origin of line l (0 shower, 1 medium recoil, 2 deleted), from its
// mother chain; a chain ends at a line without mother: the documentation line
// of the hard parton (line 2), a recoil, or a line whose origin COMPRESSEVENT kept
int lib_origin(int l) {
  using namespace py::pyjets;
  bool deleted = false;
  for (int guard = 0; guard < cb::kLines; ++guard) {
    if (K(l, 1) == -1) deleted = true;
    if (K(l, 3) <= 0) break;
    l = K(l, 3);
  }
  int o = cb::lib::origin(l);
  if (o < 0) o = (l == 2) ? 0 : 1;
  return deleted ? 2 : o;
}

// ***********************************************************************
bool compressevent(int& l1) {
  using namespace py::pyjets;
  using namespace cb::angor;
  using namespace cb::time;
  using namespace cb::colour;
  int nold, nnew = 0, nstart;
  bool compressevent_;

  nold = N;

  {
    bool found = false;
    for (int i = 2, i_end = nold; i <= i_end; ++i) {
      if (((K(i, 1) == 11) || (K(i, 1) == 12) || (K(i, 1) == 13) || (K(i, 1) == 14)) && (i != l1)) {
        nnew = i;
        found = true;
        break;  // goto 778
      }
    }  // 777
    if (!found) {
      compressevent_ = false;
      return compressevent_;
    }
  }
  // 778
  // library mode: the moved lines lose their mother pointers; keep their origin
  // (shower / recoil / deleted) and the recoil lines of the stored scattering centres
  static F1<int, cb::kLines> newline;
  if (cb::lib::active) {
    for (int i = 2; i <= nold; ++i) cb::lib::origin(i) = lib_origin(i);
    for (int i = 1; i <= nold; ++i) newline(i) = (i < nnew) ? i : 0;
  }
  nstart = nnew;
  for (int i = nstart, i_end = nold; i <= i_end; ++i) {
    if (((K(i, 1) != 11) && (K(i, 1) != 12) && (K(i, 1) != 13) && (K(i, 1) != 14)) || (i == l1)) {
      for (int j = 1; j <= 5; ++j) {
        P(nnew, j) = P(i, j);
        V(nnew, j) = V(i, j);
        mv(nnew, j) = mv(i, j);
      }  // 780
      trip(nnew) = trip(i);
      anti(nnew) = anti(i);
      zd(nnew) = zd(i);
      thetad(nnew) = thetad(i);
      qqbard(nnew) = qqbard(i);
      cb::hf::qqbarfl(nnew) = cb::hf::qqbarfl(i);
      aoreq(nnew) = aoreq(i);
      K(nnew, 1) = K(i, 1);
      K(nnew, 2) = K(i, 2);
      K(nnew, 3) = 0;
      K(nnew, 4) = 0;
      K(nnew, 5) = 0;
      if (l1 == i) l1 = nnew;
      if (cb::lib::active) {
        cb::lib::origin(nnew) = cb::lib::origin(i);
        newline(i) = nnew;
      }
      nnew = nnew + 1;
    }
  }  // 779
  N = nnew - 1;
  if (cb::lib::active) {
    using cb::storescatcen::nscatcen;
    for (int j = 1; j <= nscatcen; ++j) {
      const int r = cb::lib::scatrec(j);
      if (r > 0) cb::lib::scatrec(j) = (newline(r) > 0) ? newline(r) : -2;
    }
    for (int i = nnew; i <= nold; ++i) cb::lib::origin(i) = -1;
  }
  if ((nold - N) <= 10) {
    compressevent_ = false;
  } else {
    compressevent_ = true;
  }
  for (int i = nnew, i_end = nold; i <= i_end; ++i) {
    for (int j = 1; j <= 5; ++j) {
      K(i, j) = 0;
      P(i, j) = 0.0;
      V(i, j) = 0.0;
      mv(i, j) = 0.0;
    }  // 782
    trip(i) = 0;
    anti(i) = 0;
    zd(i) = 0.0;
    thetad(i) = 0.0;
    qqbard(i) = false;
    cb::hf::qqbarfl(i) = 0;
    // NOTE(port): aoreq(i) is not reset here in the Fortran either; kept as is.
  }  // 781
  if (N > 23000) io::log(" Error in compressevent: n = {}", N);
  if (l1 > N) io::log(" Error in compressevent: l1 = {}", l1);
  io::log_flush();
  return compressevent_;
}

// ***********************************************************************
// ***	  subroutine pevrec
// ***********************************************************************
void pevrec(const int& num, const bool& col) {
  using namespace py::pyjets;
  using namespace cb::angor;
  using namespace cb::time;
  using namespace cb::colour;
  // 5000 FORMAT(I6,' (',I6,',',I6')    ['3I6,' ]  {',I5,I3,' } ',2E14.6,l2)

  for (int i = 1, i_end = N; i <= i_end; ++i) {
    V(i, 1) = mv(i, 1);
    V(i, 2) = mv(i, 2);
    V(i, 3) = mv(i, 3);
    V(i, 4) = mv(i, 4);
    V(i, 5) = mv(i, 5);
    if (col)
      io::log("{:6} ({:6},{:6})    [{:6}{:6}{:6} ]  {{{:5}{:3} }} {}{}{:>2}", i, trip(i), anti(i), K(i, 3),
              K(i, 4), K(i, 5), K(i, 2), K(i, 1), fedit(zd(i), 14, 6), fedit(thetad(i), 14, 6),
              aoreq(i) ? "T" : "F");
  }  // 202
  py::pylist(num);
}

// ***********************************************************************
// ***	  subroutine putparticlesonshell
// ***********************************************************************
void putparticlesonshell() {
  using namespace py::pyjets;
  using namespace py::pypars;
  using namespace cb::storescatcen;
  using namespace cb::evrecord;
  using namespace cb::hepmcid;
  int lstart;
  F1<double, 4> p4totlab, p4totcms;
  // NOTE(port): NEWENS(23000) is a >64 kB local, i.e. static in gfortran. It is
  // indexed from LSTART, which is 0 for COLLIDER='PPJJ'; NEWENS(0) is out of
  // bounds in the Fortran (overwrites whatever static precedes it). Here the
  // array is given the index range 0:23000 so that element 0 is well defined.
  static FArray<double, Dim<0, 23000>> newens;
  double xmt, et, x, x2, f0, g0, acc;
  constexpr double relacc = 1.e-6;  // data relacc/1.d-6/

  if (putonshell == 0) {
    return;
  } else if (putonshell == 1) {
    for (int i = 1, i_end = N; i <= i_end; ++i) {
      P(i, 4) = std::sqrt(py::pyp(i, 7) + ipow<2>(py::pymass(K(i, 2))));
      P(i, 5) = py::pymass(K(i, 2));
    }  // 700
  } else if (putonshell == 2) {
    if ((collider == "EEJJ") || (collider == "PPDY")) {
      return;
    } else if (collider == "PPJJ") {
      // NOTE(port): LSTART=0 makes the loops below touch line 0 of PYJETS
      // (P(0,j), K(0,2)), which is outside the Fortran array bounds and aliases
      // the preceding storage (K(22999..23000,5) / P(23000,j-1) / K(23000,1)).
      // The FArray index arithmetic yields the same addresses as gfortran, so
      // the behaviour is reproduced, but it is formally out of bounds
      // (will trip JEWEL_BOUNDS_CHECK). Kept as in the original.
      lstart = 0;
    } else if ((collider == "PPYJ") || (collider == "PPYQ") || (collider == "PPYG")) {
      lstart = 1;
    } else {
      lstart = 2;
    }
    // C--boost to c.m.s.
    p4totlab(1) = 0.0;
    p4totlab(2) = 0.0;
    p4totlab(3) = 0.0;
    p4totlab(4) = 0.0;
    for (int i = lstart, i_end = N; i <= i_end; ++i) {
      p4totlab(1) = p4totlab(1) + P(i, 1);
      p4totlab(2) = p4totlab(2) + P(i, 2);
      p4totlab(3) = p4totlab(3) + P(i, 3);
      p4totlab(4) = p4totlab(4) + P(i, 4);
    }  // 701
    py::pyrobo(lstart, N, 0.0, 0.0, -p4totlab(1) / p4totlab(4), -p4totlab(2) / p4totlab(4),
               -p4totlab(3) / p4totlab(4));
    // C--first make all momenta massless
    p4totcms(1) = 0.0;
    p4totcms(2) = 0.0;
    p4totcms(3) = 0.0;
    p4totcms(4) = 0.0;
    xmt = 0.0;
    for (int i = lstart, i_end = N; i <= i_end; ++i) {
      xmt = xmt + std::sqrt(std::max(ipow<2>(P(i, 4)) - py::pyp(i, 7), 0.0));
      p4totcms(1) = p4totcms(1) + P(i, 1);
      p4totcms(2) = p4totcms(2) + P(i, 2);
      p4totcms(3) = p4totcms(3) + P(i, 3);
      p4totcms(4) = p4totcms(4) + P(i, 4);
    }  // 702
    et = p4totcms(4);
    x = 1.0 / std::sqrt(1.0 - ipow<2>(xmt / et));
    acc = relacc * et;
    for (int j = 1; j <= 10; ++j) {
      f0 = -et;
      g0 = 0.0;
      x2 = ipow<2>(x);
      for (int i = lstart, i_end = N; i <= i_end; ++i) {
        newens(i) = std::sqrt(x2 * py::pyp(i, 7));
        f0 = f0 + newens(i);
        g0 = g0 + py::pyp(i, 7) / newens(i);
      }  // 704
      if (std::fabs(f0) < acc) break;  // goto 705
      x = x - f0 / (x * g0);
    }  // 703
    // 705
    for (int i = lstart, i_end = N; i <= i_end; ++i) {
      P(i, 1) = x * P(i, 1);
      P(i, 2) = x * P(i, 2);
      P(i, 3) = x * P(i, 3);
      P(i, 4) = newens(i);
      P(i, 5) = std::sqrt(std::fabs(ipow<2>(P(i, 4)) - py::pyp(i, 7)));
    }  // 706
    // C--now put all particles on their mass shell
    p4totcms(1) = 0.0;
    p4totcms(2) = 0.0;
    p4totcms(3) = 0.0;
    p4totcms(4) = 0.0;
    xmt = 0.0;
    for (int i = lstart, i_end = N; i <= i_end; ++i) {
      xmt = xmt + py::pymass(K(i, 2));
      p4totcms(1) = p4totcms(1) + P(i, 1);
      p4totcms(2) = p4totcms(2) + P(i, 2);
      p4totcms(3) = p4totcms(3) + P(i, 3);
      p4totcms(4) = p4totcms(4) + P(i, 4);
    }  // 707
    if (p4totcms(4) > xmt) {
      et = p4totcms(4);
      x = 1.0 / std::sqrt(1.0 - ipow<2>(xmt / et));
      acc = relacc * et;
      for (int j = 1; j <= 10; ++j) {
        f0 = -et;
        g0 = 0.0;
        x2 = ipow<2>(x);
        for (int i = lstart, i_end = N; i <= i_end; ++i) {
          newens(i) = std::sqrt(ipow<2>(py::pymass(K(i, 2))) + x2 * ipow<2>(P(i, 4)));
          f0 = f0 + newens(i);
          g0 = g0 + ipow<2>(P(i, 4)) / newens(i);
        }  // 709
        if (std::fabs(f0) < acc) break;  // goto 710
        x = x - f0 / (x * g0);
      }  // 708
      // 710
      for (int i = lstart, i_end = N; i <= i_end; ++i) {
        P(i, 1) = x * P(i, 1);
        P(i, 2) = x * P(i, 2);
        P(i, 3) = x * P(i, 3);
        P(i, 4) = newens(i);
        P(i, 5) = std::sqrt(std::fabs(ipow<2>(P(i, 4)) - py::pyp(i, 7)));
      }  // 711
    } else {
      io::log(" WARNING: not enough energy available to put particles on shell.");
    }
    py::pyrobo(lstart, N, 0.0, 0.0, p4totlab(1) / p4totlab(4), p4totlab(2) / p4totlab(4),
               p4totlab(3) / p4totlab(4));
    for (int i = lstart, i_end = N; i <= i_end; ++i) {
      P(i, 4) = std::sqrt(py::pyp(i, 7) + ipow<2>(py::pymass(K(i, 2))));
      P(i, 5) = py::pymass(K(i, 2));
    }  // 712
  }
}

// ***********************************************************************
// ***	  subroutine converttohepmc
// ***********************************************************************
// The Fortran writes to unit J (= HPMCFID); the C++ HepMC file is io::hepmc().
void converttohepmc([[maybe_unused]] const int& j, const int& evnum, const int& pid, const std::string& beam1,
                    const std::string& beam2) {
  using namespace py::pyjets;
  using namespace py::pypars;
  using namespace cb::param;
  using namespace cb::evrecord;
  using namespace cb::storescatcen;
  using namespace cb::storedummies;
  // NOTE(port): the Fortran declares Rsub LOGICAL in this routine's copy of
  // /STORESCATCEN/ (DOUBLE PRECISION elsewhere). Rsub is not referenced in
  // CONVERTTOHEPMC, so the mismatch has no effect; the double cb::storescatcen::rsub is used.
  //
  // NOTE(port): CODELIST(25000) is a 100 kB local array, which gfortran places
  // in static storage (> 64 kB limit), so its contents survive between calls
  // (and start as zero). Entries that are read before being set in the current
  // event therefore see values from earlier events; reproduced with `static`.
  static F1<int, 25000> codelist;
  int pbarcode, vbarcode, i, nfirst, nvertex, codefirst;
  // NOTE(port): possibly uninitialised in Fortran (NSTART is only set if some
  // line has K(I,3).NE.0 in the hadronised branch); initialised to 0.
  int nstart = 0;
  double centr;
  const double mproton = r4(0.9383f);   // data mproton/0.9383/
  const double mneutron = r4(0.9396f);  // data mneutron/0.9396/

  pbarcode = 0;
  vbarcode = 0;

  centr = medium().getcentrality();

  if (shorthepmc) {
    // C--short output
    if (collider == "EEJJ") {
      nvertex = 3;
      pbarcode = 5;
    } else {
      nvertex = 1;
      pbarcode = 2;
    }
    nfirst = 0;
    for (i = 1; i <= N; ++i) {
      if (((K(i, 1) < 6) || (K(i, 1) == 17))) nfirst = nfirst + 1;
    }  // 131
    if (writescatcen) nfirst = nfirst + nscatcen;
    if (writedummies) nfirst = nfirst + nscatcen;

    w5000("E ", evnum, -1, 0.0, 0.0, 0.0, 0, 0, nvertex, 1, 2, 0, 1, PARI(10));
    wnline();
    wuline();
    w5100("C ", PARI(1) * 1.e9, 0.0);
    w5200("H ", 0, 0, 0, 0, 0, 0, 0, 0, 0, centr, 0.0, 0.0, 0.0);
    w5300("F ", 0, 0, -1.0, -1.0, -1.0, -1.0, -1.0, 0, 0);
    // C--write out vertex line
    if (collider == "EEJJ") {
      w5400("V ", -1, 0, 0, 0, 0, 0, 2, 1, 0);
      w5500("P ", 1, -11, 0.0, 0.0, sqrts / 2.0, sqrts / 2.0, r4(0.00051f), 2, 0, 0, -1, 0);
      w5500("P ", 2, 11, 0.0, 0.0, -sqrts / 2.0, sqrts / 2.0, r4(0.00051f), 2, 0, 0, -1, 0);
      w5500("P ", 3, 23, 0.0, 0.0, 0.0, sqrts, r4(91.2f), 2, 0, 0, -2, 0);
      w5400("V ", -2, 0, 0, 0, 0, 0, 0, 2, 0);
      w5500("P ", 4, pid, sqrts / 2.0, 0.0, 0.0, sqrts / 2.0, 0.0, 2, 0, 0, -3, 0);
      w5500("P ", 5, -pid, -sqrts / 2.0, 0.0, 0.0, sqrts / 2.0, 0.0, 2, 0, 0, -3, 0);
      w5400("V ", -3, 0, 0, 0, 0, 0, 0, nfirst, 0);
    } else {
      w5400("V ", -1, 0, 0, 0, 0, 0, 2, nfirst, 0);
      if (beam1 == "p+") {
        w5500("P ", 1, 2212, 0.0, 0.0, std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mproton)), sqrts / 2.0, mproton,
              4, 0, 0, -1, 0);
      } else {
        w5500("P ", 1, 2212, 0.0, 0.0, std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mneutron)), sqrts / 2.0,
              mneutron, 4, 0, 0, -1, 0);
      }
      if (beam2 == "p+") {
        w5500("P ", 2, 2212, 0.0, 0.0, -std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mproton)), sqrts / 2.0,
              mproton, 4, 0, 0, -1, 0);
      } else {
        w5500("P ", 2, 2212, 0.0, 0.0, -std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mneutron)), sqrts / 2.0,
              mneutron, 4, 0, 0, -1, 0);
      }
    }
    // C--write out scattering centres
    if (writescatcen) {
      for (i = 1; i <= nscatcen; ++i) {
        pbarcode = pbarcode + 1;
        w5500("P ", pbarcode, scatflav(i), scatcen(i, 1), scatcen(i, 2), scatcen(i, 3), scatcen(i, 4),
              scatcen(i, 5), 3, 0, 0, 0, 0);
      }  // 133
    }
    // C--write out dummy particles
    if (writedummies) {
      for (i = 1; i <= nscatcen; ++i) {
        pbarcode = pbarcode + 1;
        w5500("P ", pbarcode, 111, dummies(i, 1), dummies(i, 2), dummies(i, 3), dummies(i, 4), 0.0, 1, 0, 0, 0,
              0);
      }  // 137
    }
    // C--write out particle lines
    for (i = 1; i <= N; ++i) {
      if (((K(i, 1) < 6) || (K(i, 1) == 17))) {
        pbarcode = pbarcode + 1;
        if ((K(i, 1) == 3) || (K(i, 1) == 5)) {
          w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 4, 0, 0, 0, 0);
        } else {
          w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 1, 0, 0, 0, 0);
        }
      }
    }  // 132

  } else {
    // C--long output
    if (hadro) {
      // C--hadronised events
      nfirst = 0;
      if (collider == "EEJJ") {
        nvertex = 3;
      } else {
        nvertex = 1;
      }
      for (i = 1; i <= N; ++i) {
        if (K(i, 3) != 0) {
          nstart = i;
          break;  // GOTO 124
        }
      }  // 123
      // 124
      // !	    nstart=0

      for (i = nstart + 1; i <= N; ++i) {
        if (isprimhadron(i)) nfirst = nfirst + 1;
        if ((ishadron(K(i, 2)) || (std::abs(K(i, 2)) == 15)) && (K(i, 4) != 0)) nvertex = nvertex + 1;
      }  // 126
      // 127

      if (writescatcen) nfirst = nfirst + nscatcen;
      if (writedummies) nfirst = nfirst + nscatcen;

      w5000("E ", evnum, -1, 0.0, 0.0, 0.0, 0, 0, nvertex, 1, 2, 0, 1, PARI(10));
      wnline();
      wuline();
      w5100("C ", PARI(1) * 1.e9, 0.0);
      w5200("H ", 0, 0, 0, 0, 0, 0, 0, 0, 0, centr, 0.0, 0.0, 0.0);
      w5300("F ", 0, 0, -1.0, -1.0, -1.0, -1.0, -1.0, 0, 0);

      // C--write out vertex line
      if (collider == "EEJJ") {
        vbarcode = -3;
        pbarcode = 5;
      } else {
        vbarcode = -1;
        pbarcode = 2;
      }
      if (collider == "EEJJ") {
        w5400("V ", -1, 0, 0, 0, 0, 0, 2, 1, 0);
        w5500("P ", 1, -11, 0.0, 0.0, sqrts / 2.0, sqrts / 2.0, r4(0.00051f), 2, 0, 0, -1, 0);
        w5500("P ", 2, 11, 0.0, 0.0, -sqrts / 2.0, sqrts / 2.0, r4(0.00051f), 2, 0, 0, -1, 0);
        w5500("P ", 3, 23, 0.0, 0.0, 0.0, sqrts, r4(91.2f), 2, 0, 0, -2, 0);
        w5400("V ", -2, 0, 0, 0, 0, 0, 0, 2, 0);
        w5500("P ", 4, pid, sqrts / 2.0, 0.0, 0.0, sqrts / 2.0, 0.0, 2, 0, 0, -3, 0);
        w5500("P ", 5, -pid, -sqrts / 2.0, 0.0, 0.0, sqrts / 2.0, 0.0, 2, 0, 0, -3, 0);
        w5400("V ", vbarcode, 0, 0, 0, 0, 0, 0, nfirst, 0);
      } else {
        w5400("V ", -1, 0, 0, 0, 0, 0, 2, nfirst, 0);
        if (beam1 == "p+") {
          w5500("P ", 1, 2212, 0.0, 0.0, std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mproton)), sqrts / 2.0,
                mproton, 4, 0, 0, -1, 0);
        } else {
          w5500("P ", 1, 2212, 0.0, 0.0, std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mneutron)), sqrts / 2.0,
                mneutron, 4, 0, 0, -1, 0);
        }
        if (beam2 == "p+") {
          w5500("P ", 2, 2212, 0.0, 0.0, -std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mproton)), sqrts / 2.0,
                mproton, 4, 0, 0, -1, 0);
        } else {
          w5500("P ", 2, 2212, 0.0, 0.0, -std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mneutron)), sqrts / 2.0,
                mneutron, 4, 0, 0, -1, 0);
        }
      }

      codefirst = nfirst + pbarcode;

      // C--write out scattering centres
      if (writescatcen) {
        for (i = 1; i <= nscatcen; ++i) {
          pbarcode = pbarcode + 1;
          w5500("P ", pbarcode, scatflav(i), scatcen(i, 1), scatcen(i, 2), scatcen(i, 3), scatcen(i, 4),
                scatcen(i, 5), 3, 0, 0, 0, 0);
        }  // 134
      }
      // C--write out dummy particles
      if (writedummies) {
        for (i = 1; i <= nscatcen; ++i) {
          pbarcode = pbarcode + 1;
          w5500("P ", pbarcode, 111, dummies(i, 1), dummies(i, 2), dummies(i, 3), dummies(i, 4), 0.0, 1, 0, 0,
                0, 0);
        }  // 138
      }

      // C--first write out all particles coming directly from string or cluster decays
      for (i = nstart + 1; i <= N; ++i) {
        if (!isprimhadron(i)) {
          continue;  // GOTO 125
        } else {
          if (pbarcode == codefirst) break;  // GOTO 130
          pbarcode = pbarcode + 1;
          // C--write out particle line
          if (K(i, 4) > 0) {
            vbarcode = vbarcode - 1;
            codelist(i) = vbarcode;
            w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 2, 0, 0, vbarcode, 0);
          } else {
            w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 1, 0, 0, 0, 0);
          }
        }
      }  // 125
      // 130
      // C--now write out all other particles and vertices
      for (i = nstart + 1; i <= N; ++i) {
        if (isprimhadron(i) || isprimstring(i)) continue;  // goto 129
        if (isparton(K(i, 2))) {
          if (ishadron(K(K(i, 3), 2))) codelist(i) = codelist(K(i, 3));
          continue;  // goto 129
        }
        if (issecstring(i)) {
          codelist(i) = codelist(K(i, 3));
          continue;  // goto 129
        }
        pbarcode = pbarcode + 1;
        if ((K(i, 3) != K(i - 1, 3))) {
          // C--write out vertex line
          w5400("V ", codelist(K(i, 3)), 0, 0, 0, 0, 0, 0, K(K(i, 3), 5) - K(K(i, 3), 4) + 1, 0);
        }
        // C--write out particle line
        if (K(i, 4) > 0) {
          vbarcode = vbarcode - 1;
          codelist(i) = vbarcode;
          w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 2, 0, 0, vbarcode, 0);
        } else {
          w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 1, 0, 0, 0, 0);
        }
      }  // 129

    } else {
      // C---------------------------------------------------------------------------------------
      // C--partonic events
      // !	    call pevrec(2,.false.)
      // C--hadronised events
      nfirst = 0;
      if (collider == "EEJJ") {
        nvertex = 3;
      } else {
        nvertex = 1;
      }

      for (i = 9; i <= N; ++i) {
        if ((K(i, 3) == 1) || (K(i, 3) == 2) || (K(i, 3) == 7) || (K(i, 3) == 8)) nfirst = nfirst + 1;
        if (K(i, 4) != 0) nvertex = nvertex + 1;
      }  // 150
      nstart = 9 + nfirst;

      if (writescatcen) nfirst = nfirst + nscatcen;
      if (writedummies) nfirst = nfirst + nscatcen;

      w5000("E ", evnum, -1, 0.0, 0.0, 0.0, 0, 0, nvertex, 1, 2, 0, 1, PARI(10));
      wnline();
      wuline();
      w5100("C ", PARI(1) * 1.e9, 0.0);
      w5200("H ", 0, 0, 0, 0, 0, 0, 0, 0, 0, centr, 0.0, 0.0, 0.0);
      w5300("F ", 0, 0, -1.0, -1.0, -1.0, -1.0, -1.0, 0, 0);

      // C--write out vertex line
      if (collider == "EEJJ") {
        vbarcode = -3;
        pbarcode = 5;
      } else {
        vbarcode = -1;
        pbarcode = 2;
      }
      if (collider == "EEJJ") {
        w5400("V ", -1, 0, 0, 0, 0, 0, 2, 1, 0);
        w5500("P ", 1, -11, 0.0, 0.0, sqrts / 2.0, sqrts / 2.0, r4(0.00051f), 2, 0, 0, -1, 0);
        w5500("P ", 2, 11, 0.0, 0.0, -sqrts / 2.0, sqrts / 2.0, r4(0.00051f), 2, 0, 0, -1, 0);
        w5500("P ", 3, 23, 0.0, 0.0, 0.0, sqrts, r4(91.2f), 2, 0, 0, -2, 0);
        w5400("V ", -2, 0, 0, 0, 0, 0, 0, 2, 0);
        w5500("P ", 4, pid, sqrts / 2.0, 0.0, 0.0, sqrts / 2.0, 0.0, 2, 0, 0, -3, 0);
        w5500("P ", 5, -pid, -sqrts / 2.0, 0.0, 0.0, sqrts / 2.0, 0.0, 2, 0, 0, -3, 0);
        w5400("V ", vbarcode, 0, 0, 0, 0, 0, 0, nfirst, 0);
      } else {
        w5400("V ", -1, 0, 0, 0, 0, 0, 2, nfirst, 0);
        if (beam1 == "p+") {
          w5500("P ", 1, 2212, 0.0, 0.0, std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mproton)), sqrts / 2.0,
                mproton, 4, 0, 0, -1, 0);
        } else {
          w5500("P ", 1, 2212, 0.0, 0.0, std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mneutron)), sqrts / 2.0,
                mneutron, 4, 0, 0, -1, 0);
        }
        if (beam2 == "p+") {
          w5500("P ", 2, 2212, 0.0, 0.0, -std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mproton)), sqrts / 2.0,
                mproton, 4, 0, 0, -1, 0);
        } else {
          w5500("P ", 2, 2212, 0.0, 0.0, -std::sqrt(ipow<2>(sqrts) / 4.0 - ipow<2>(mneutron)), sqrts / 2.0,
                mneutron, 4, 0, 0, -1, 0);
        }
      }

      // C--write out scattering centres
      if (writescatcen) {
        for (i = 1; i <= nscatcen; ++i) {
          pbarcode = pbarcode + 1;
          w5500("P ", pbarcode, scatflav(i), scatcen(i, 1), scatcen(i, 2), scatcen(i, 3), scatcen(i, 4),
                scatcen(i, 5), 3, 0, 0, 0, 0);
        }  // 151
      }
      // C--write out dummy particles
      if (writedummies) {
        for (i = 1; i <= nscatcen; ++i) {
          pbarcode = pbarcode + 1;
          w5500("P ", pbarcode, 111, dummies(i, 1), dummies(i, 2), dummies(i, 3), dummies(i, 4), 0.0, 1, 0, 0,
                0, 0);
        }  // 152
      }

      // C--write out outgoing particles of first vertex
      for (i = 9; i <= nstart - 1; ++i) {
        pbarcode = pbarcode + 1;
        // C--write out particle line
        if (K(i, 4) > 0) {
          vbarcode = vbarcode - 1;
          codelist(i) = vbarcode;
          w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 2, 0, 0, vbarcode, 0);
        } else {
          w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 1, 0, 0, 0, 0);
        }
      }  // 154

      // C--now write out all other particles and vertices
      for (i = nstart; i <= N; ++i) {
        pbarcode = pbarcode + 1;
        if (K(i, 3) == 0) {
          // C--write out vertex line - scattering
          w5400("V ", codelist(K(i + 1, 3)), K(K(i + 1, 3), 1), 0, 0, 0, 0, 0,
                K(K(i + 1, 3), 5) - K(K(i + 1, 3), 4) + 1, 0);
        } else if ((K(i, 3) != K(i - 1, 3)) && (K(i - 1, 3) != 0)) {
          // C--write out vertex line - splitting
          w5400("V ", codelist(K(i, 3)), K(K(i, 3), 1), 0, 0, 0, 0, 0, K(K(i, 3), 5) - K(K(i, 3), 4) + 1, 0);
        }
        // C--write out particle line
        if (K(i, 4) > 0) {
          vbarcode = vbarcode - 1;
          codelist(i) = vbarcode;
          w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 2, 0, 0, vbarcode, 0);
        } else {
          if ((K(i, 1) == 3) || (K(i, 1) == 5)) {
            w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 4, 0, 0, 0, 0);
          } else if ((K(i, 1) == 11) && (K(i, 3) == 0)) {
            w5500("P ", pbarcode, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0);
          } else {
            w5500("P ", pbarcode, K(i, 2), P(i, 1), P(i, 2), P(i, 3), P(i, 4), P(i, 5), 1, 0, 0, 0, 0);
          }
        }
      }  // 153
    }
  }
  std::fflush(io::hepmc());  // call flush(j)
}

// ***********************************************************************
// ***	  subroutine copyline
// ***********************************************************************
void copyline(const int& nfr, const int& nto, const int& mode) {
  using namespace py::pyjets;
  using namespace cb::time;
  using namespace cb::colour;
  using namespace cb::angor;

  if (mode >= 0) {
    K(nto, 1) = K(nfr, 1);
    K(nto, 2) = K(nfr, 2);
    for (int j = 1; j <= 5; ++j) {
      P(nto, j) = P(nfr, j);
    }  // 100
    zd(nto) = zd(nfr);
    thetad(nto) = thetad(nfr);
    qqbard(nto) = qqbard(nfr);
    cb::hf::qqbarfl(nto) = cb::hf::qqbarfl(nfr);
    aoreq(nto) = aoreq(nfr);
  }
  if (mode >= 1) {
    for (int j = 1; j <= 5; ++j) {
      mv(nto, j) = mv(nfr, j);
    }  // 101
  }
  if (mode == 2) {
    trip(nto) = trip(nfr);
    anti(nto) = anti(nfr);
  }
}

// ***********************************************************************
// ***	  subroutine printlogo
// ***********************************************************************
// The logo and copyright/program messages are reproduced verbatim (MCnet
// guidelines). Each Fortran `write(fid,*)'...'` is a list-directed record,
// which gfortran starts with one blank; an empty `write(fid,*)` gives an
// empty record.
void printlogo(const int& fid) {
  constexpr std::string_view logo[] = {
      "                   _________________________________________                  ",
      "                  |                                         |                 ",
      "                  |  JJJJJ  EEEEE  W       W  EEEEE  L      |                  ",
      "                  |      J  E      W       W  E      L      |                  ",
      " _________________|      J  EEE     W  W  W   EEE    L      |_________________ ",
      "|                 |  J   J  E       W W W W   E      L      |                 |",
      "|                 |   JJJ   EEEEE    W   W    EEEEE  LLLLL  |                 |",
      "|                 |_________________________________________|                 |",
      "|                                                                             |",
      "|                            this is JEWEL 2.6.0                              |",
      "|                                                                             |",
      "| Copyright Korinna C. Zapp (2026)  [Korinna.Zapp@fysik.lu.se]                |",
      "|                                                                             |",
      "| The JEWEL homepage is jewel.hepforge.org                                    |",
      "|                                                                             |",
      "| The medium model was partly implemented by Jochen Klein.                    |",
      "| Rithya (Raghav) Kunnawalkam Elayavalli helped with the implementation of    |",
      "| the V+jet processes.                                                        |",
      "|                                                                             |",
      "| Please cite JHEP 1303 (2013) 080 [arXiv:1212.1599] for physics and          |",
      "| EPJC 74 (2014) no.2, 2762 [arXiv:1311.0048]  for the code.                  |",
      "| The reference for V+jet processes is EPJC 76 (2016) no.12 695               |",
      "| [arXiv:1608.03099] and for recoil effects it is JHEP 07 (2017) 141          |",
      "| [arXiv:1707.01539] and EPJ C 82 (2022) no.11,  1010 [arXiv:2207.14814].     |",
      "| Quenching of initial state radiation is discussed in PLB 835 (2022) 137567  |",
      "| [arXiv:2208.00813] and multi-parton interactions in arXiv:2407.03066.       |",
      "| Colour coherence is introduced in arXiv:2604.13932.                         |",
      "|                                                                             |",
      "| JEWEL relies heavily on PYTHIA 6 for the event generation. The modified     |",
      "| version of PYTHIA 6.4.25 that is shipped with JEWEL is, however, not an     |",
      "| official PYTHIA release and must not be used for anything else. Please      |",
      "| refer to results as \"JEWEL+PYTHIA\".                                         |",
      "|                                                                             |",
      "| JEWEL also uses code provided by S. Zhang and J. M. Jing                    |",
      "| (Computation of Special Functions, John Wiley & Sons, New York, 1996 and    |",
      "| http://jin.ece.illinois.edu) for computing the exponential integral Ei(x).  |",
      "|                                                                             |",
      "|_____________________________________________________________________________|",
  };

  io::fortran_write(fid, "");
  for (const std::string_view l : logo) io::fortran_write(fid, std::string(" ") + std::string(l));
  // Addition of the C++ port (not part of the original logo).
  io::fortran_write(fid, " C++ port of JEWEL 2.6.0 (translated from the original Fortran; PYTHIA 6.4.25 unchanged)");
  io::fortran_write(fid, "");
  io::fortran_write(fid, "");
}

// ***********************************************************************
// ***	  subroutine printtime
// ***********************************************************************
// 1000 format (i2.2, '.', i2.2, '.', i4.4, ', ', i2.2, ':', i2.2, ':', i2.2 )
// IDATE(date): date(1)=day, date(2)=month, date(3)=year; ITIME(time):
// hour, minute, second -- both in local time.
void printtime() {
  const std::time_t now = std::time(nullptr);
  std::tm lt{};
  localtime_r(&now, &lt);
  io::log_line(jewel::jfmt::format("{:02}.{:02}.{:04}, {:02}:{:02}:{:02}", lt.tm_mday, lt.tm_mon + 1, lt.tm_year + 1900,
                           lt.tm_hour, lt.tm_min, lt.tm_sec));
}

}  // namespace jewel
