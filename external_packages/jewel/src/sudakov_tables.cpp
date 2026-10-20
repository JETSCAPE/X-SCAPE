// JEWEL 2.6.0 C++ port -- Sudakov factors, splitting/pdf/cross-section integrals, lookup tables,
// exponential integral, splitting functions, alpha_s and the adaptive Runge-Kutta integrator.
// Translated from jewel-2.6.0.f:6459-8422.
//
// Note on single precision: in Fortran `2.*NF` (REAL*4 literal times INTEGER) is REAL*4, so the
// ubiquitous factor `(33.-2.*NF)` and `(11.-2.*NF/3.)` are evaluated in single precision and only
// then promoted to DOUBLE PRECISION. This is reproduced explicitly with float arithmetic below.

#include <algorithm>
#include <cmath>
#include <string>

#include "jewel/commons.hpp"
#include "jewel/farray.hpp"
#include "jewel/fortran_compat.hpp"
#include "jewel/heavyflavour.hpp"
#include "jewel/io.hpp"
#include "jewel/medium.hpp"
#include "jewel/pythia6.hpp"
#include "jewel/routines.hpp"

namespace jewel {


namespace {
// CHARACTER*2 literals passed as arguments (hoisted to avoid per-call allocation; same values)
const std::string kQQ = "QQ";
const std::string kQG = "QG";
const std::string kGQ = "GQ";
const std::string kGG = "GG";
const std::string kGC = "GC";
const std::string kQA = "QA";
const std::string kGA = "GA";
const std::string kGB = "GB";
}  // namespace

// ***********************************************************************
// ***	  function getproba
// ***********************************************************************
double getproba(const double& qi, const double& qf, const double& qaa, const double& zaa,
                const double& ebb, const std::string& type, const double& t1, const bool& ins2) {
  // C--variables for Sudakov integration
  using namespace cb::sudaint;

  qa = qaa;
  za2 = zaa;
  eb = ebb;
  typ = type;
  t = t1;
  instate = ins2;
  const double s = getsudakov(qi, qaa, qf, zaa, ebb, type, t1, ins2);
  return s * deriv(qf, 1);
}

// ***********************************************************************
// ***	  function getsudakov
// ***********************************************************************
double getsudakov(const double& qmax1, const double& qa1, const double& qb1, const double& za1,
                  const double& eb1, const std::string& type3, const double& t2, const bool& ins) {
  // C--Parameter common block
  using namespace cb::param;
  // C--variables for Sudakov integration
  using namespace cb::sudaint;
  // C--local variables
  constexpr double epsi = 1.e-4;  // DATA EPSI/1.d-4/
  double getsudakov_ = 0.0;

  double qb2 = qb1;
  if (ins) {
    if (qb2 < q0) io::log(" error: Q < Q0 {} {}", qb2, qmax1);
    if (qb2 < (q0 + 1.e-10)) qb2 = qb2 + 1.e-10;
  } else {
    if (qb2 < q0) io::log(" error: Q < min {} {}", qb2, qmax1);
    if (qb2 < (q0 + 1.e-10)) qb2 = qb2 + 1.e-10;
  }
  if (qb2 >= (qmax1 - 1.e-10)) {
    getsudakov_ = 1.0;
  } else {
    if (ins) {
      getsudakov_ = getinsudafast(qb1, qmax1, type3);
    } else {
      qa = qa1;
      za2 = za1;
      eb = eb1;
      typ = type3;
      t = t2;
      instate = false;
      const double hfirst = r4(0.01f) * (qmax1 - qb1);
      double ystart = 0.0;
      odeint(ystart, qb2, qmax1, epsi, hfirst, 0.0, 1);
      getsudakov_ = std::exp(-ystart);
    }
  }
  return getsudakov_;
}

// ***********************************************************************
// ***	  function getinsudakov
// ***********************************************************************
double getinsudakov(const double& qb, const double& qmax1, const std::string& type3) {
  // C--Parameter common block
  using namespace cb::param;
  // C--variables for Sudakov integration
  using namespace cb::sudaint;
  // C--local variables
  constexpr double epsi = 1.e-4;  // DATA EPSI/1.d-4/
  double getinsudakov_ = 0.0;

  double qb1 = qb;
  if (qb1 < q0) io::log(" error: Q < Q0 {} {}", qb1, qmax1);
  if (qb1 < (q0 + 1.e-12)) qb1 = qb1 + 1.e-12;
  if (qb1 >= (qmax1 - 1.e-12)) {
    getinsudakov_ = 1.0;
  } else {
    typ = type3;
    const double hfirst = r4(0.01f) * (qmax1 - qb1);
    double ystart = 0.0;
    odeint(ystart, qb1, qmax1, epsi, hfirst, 0.0, 6);
    getinsudakov_ = std::exp(-ystart);
  }
  return getinsudakov_;
}

// ***********************************************************************
// ***	  function deriv
// ***********************************************************************
double deriv(const double& xval, const int& w4) {
  // C--Parameter common block
  using namespace cb::param;
  // C--variables for splitting function integration
  using namespace cb::intsplitf;
  // C--variables for Sudakov integration
  using namespace cb::sudaint;
  // C--variables for pdf integration
  using namespace cb::pdfintv;
  // C--variables for cross section integration
  using namespace cb::xsecv;
  // C--local variables
  constexpr double pi = 3.141592653589793;
  double deriv_ = 0.0;

  if (w4 == 1) {
    // C--Sudakov integration
    if (instate) {
      deriv_ = 2.0 * getinspliti(xval, typ) / xval;
    } else {
      deriv_ = 2.0 * getspliti(qa, xval, za2, eb, typ) / xval;
    }
  } else if (w4 == 2) {
    // C--P(q->qg) integration
    deriv_ = (1.0 + fm) * alphas(xval * (1.0 - xval) * qquad / 1.0, lps) * pqq(xval) / (2.0 * pi);
  } else if (w4 == 3) {
    // C--P(g->gg) integration
    deriv_ = (1.0 + fm) * alphas(xval * (1.0 - xval) * qquad / 1.0, lps) * pgg(xval) / (2.0 * pi);
  } else if (w4 == 4) {
    // C--P(g->qq) integration
    deriv_ = (1.0 + fm) * alphas(xval * (1 - xval) * qquad / 1.0, lps) * pqg(xval) / (2.0 * pi);
  } else if (w4 == 5) {
    deriv_ = std::exp(-xval) / xval;
  } else if (w4 == 6) {
    deriv_ = 2.0 * getinspliti(xval, typ) / xval;
  } else if (w4 == 7) {
    const double s = getinsudafast(xval, xmax, kQQ);
    deriv_ = 2.0 * s * alphas((1.0 - z) * ipow<2>(xval) / 1.0, lps) * pqq(z) / (2.0 * pi * xval);
  } else if (w4 == 8) {
    const double s = getinsudafast(xval, xmax, kGC);
    deriv_ = 2.0 * s * alphas((1.0 - z) * ipow<2>(xval) / 1.0, lps) * pgq(z) / (2.0 * pi * xval);
  } else if (w4 == 9) {
    const double s = getinsudafast(xval, xmax, kQQ);
    deriv_ = 2.0 * s * alphas((1.0 - z) * ipow<2>(xval) / 1.0, lps) * pqg(z) / (2.0 * pi * xval);
  } else if (w4 == 10) {
    // NOTE(port): the Fortran continuation lines read `...*ALPHAS(...)*` / `*2.*PGG(Z)...`, i.e.
    // `ALPHAS(...)**2.*PGG(Z)` -- alpha_s squared instead of a factor 2 (as in GETPDFXINTEXACT,
    // W4=24). Kept as in the original; `**2.` (real exponent 2.0) is folded to x*x by gfortran.
    const double s = getinsudafast(xval, xmax, kGC);
    const double as = alphas((1.0 - z) * ipow<2>(xval) / 1.0, lps);
    deriv_ = 2.0 * s * (as * as) * pgg(z) / (2.0 * pi * xval);
  } else if (w4 == 11) {
    const double s = getinspliti(scalefacm * std::sqrt(xval), kGQ);
    const double sp = scatprimfunc(xval, mdx);
    deriv_ = 3.0 * s * sp / (2.0 * xval);
  } else if (w4 == 12) {
    const double s = getinspliti(scalefacm * std::sqrt(xval), kQG);
    const double sp = scatprimfunc(xval, mdx);
    deriv_ = 2.0 * s * sp / (3.0 * xval);
  } else if (w4 == 13) {
    const double s = getinsudafast(qlow, scalefacm * std::sqrt(xval), kGC);
    const double as = alphas(xval + ipow<2>(mdx), lps);
    deriv_ = s * 3.0 * 2.0 * pi * ipow<2>(as) / (2.0 * ipow<2>(xval + ipow<2>(mdx)));
  } else if (w4 == 14) {
    const double s = getinsudafast(qlow, scalefacm * std::sqrt(xval), kQQ);
    const double as = alphas(xval + ipow<2>(mdx), lps);
    deriv_ = s * 2.0 * 2.0 * pi * ipow<2>(as) / (3.0 * ipow<2>(xval + ipow<2>(mdx)));
  } else if (w4 == 21) {
    const double s = getinsudafast(xval, xmax, kQQ);
    const double sp = getinspliti(xval, kQQ);
    deriv_ = 2.0 * s * sp / xval;
  } else if (w4 == 22) {
    const double s = getinsudafast(xval, xmax, kGC);
    const double sp = getinspliti(xval, kGQ);
    deriv_ = 2.0 * s * sp / xval;
  } else if (w4 == 23) {
    const double s = getinsudafast(xval, xmax, kQQ);
    const double sp = getinspliti(xval, kQG);
    deriv_ = 2.0 * s * sp / xval;
  } else if (w4 == 24) {
    const double s = getinsudafast(xval, xmax, kGC);
    const double sp = getinspliti(xval, kGG);
    deriv_ = 2.0 * s * 2.0 * sp / xval;
  } else {
    const int w = w4 - 100;
    deriv_ = medium().medderiv(xval, w);
  }
  return deriv_;
}

// ***********************************************************************
// ***	  function getspliti
// ***********************************************************************
double getspliti(const double& qa, const double& qb, const double& zeta, const double& eb,
                 const std::string& type1) {
  // C--Parameter common block
  using namespace cb::param;
  // C--splitting integral
  using namespace cb::splitint;
  // C--variables for splitting function integration
  using namespace cb::intsplitf;
  // C--number of extrapolations in tables
  using namespace cb::extrapolations;
  // C--local variables
  F1<double, 2> x1a, x2a, yb;
  F2<double, 2, 2> ya;
  double low, y, a, b;
  // NOTE(port): possibly uninitialised in Fortran (result for an unknown TYPE1)
  double getspliti_ = 0.0;
  double splitintgg = 0.0, splitintqg = 0.0;

  ntotspliti = ntotspliti + 1;
  if (qb > qmax) {
    noverspliti = noverspliti + 1;
    if (noverspliti <= 25) io::log(" WARNING in getspliti: need to extrapolate:  {} {}", qb, qmax);
  }

  // C--find boundaries for z integration
  if ((angord > 0) && (zeta != 1.0)) {
    low = std::max(0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(qb)) *
                             std::sqrt(1.0 - ipow<2>(qb) / ipow<2>(eb)),
                   0.5 - 0.5 * std::sqrt(1.0 - 4.0 * ipow<2>(qb) * (1.0 - zeta) / (zeta * ipow<2>(qa))));
  } else {
    low = 0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(qb)) *
                      std::sqrt(1.0 - ipow<2>(qb) / ipow<2>(eb));
  }
  // C--find values in array
  const int qlmax = fint((qb - qval(1)) * npoint / (qval(1000) - qval(1)) + 1);
  int qline = std::max(qlmax, 1);
  qline = std::min(qline, npoint);
  const int zlmax = fint((std::log(low) - std::log(zmval(1))) * npoint /
                             (std::log(zmval(1000)) - std::log(zmval(1))) +
                         1);
  int zline = std::max(zlmax, 1);
  zline = std::min(zline, npoint);
  if ((qline > 999) || (zline > 999) || (qline < 1) || (zline < 1)) {
    io::log(" ERROR in GETSPLITI: line number out of bound {} {}", qline, zline);
  }
  if ((type1 == "GG") || (type1 == "GC")) {
    for (int i = 1; i <= 2; ++i) {
      x1a(i) = qval(qline - 1 + i);
      x2a(i) = zmval(zline - 1 + i);
      for (int j = 1; j <= 2; ++j) {
        ya(i, j) = splitiggv(qline - 1 + i, zline - 1 + j);
      }
    }
    for (int i = 1; i <= 2; ++i) {
      a = (ya(i, 2) - ya(i, 1)) / (x2a(2) - x2a(1));
      b = ya(i, 1) - a * x2a(1);
      yb(i) = a * low + b;
    }
    if (x1a(1) == x1a(2)) {
      y = (yb(1) + yb(2)) / 2.0;
    } else {
      a = (yb(2) - yb(1)) / (x1a(2) - x1a(1));
      b = yb(1) - a * x1a(1);
      y = a * qb + b;
    }
    if (type1 == "GG") {
      getspliti_ = std::min(y, 10.0);
    } else {
      splitintgg = std::min(y, 10.0);
    }
  }
  if ((type1 == "QG") || (type1 == "GC")) {
    for (int i = 1; i <= 2; ++i) {
      x1a(i) = qval(qline - 1 + i);
      x2a(i) = zmval(zline - 1 + i);
      for (int j = 1; j <= 2; ++j) {
        ya(i, j) = splitiqgv(qline - 1 + i, zline - 1 + j);
      }
    }
    for (int i = 1; i <= 2; ++i) {
      a = (ya(i, 2) - ya(i, 1)) / (x2a(2) - x2a(1));
      b = ya(i, 1) - a * x2a(1);
      yb(i) = a * low + b;
    }
    if (x1a(1) == x1a(2)) {
      y = (yb(1) + yb(2)) / 2.0;
    } else {
      a = (yb(2) - yb(1)) / (x1a(2) - x1a(1));
      b = yb(1) - a * x1a(1);
      y = a * qb + b;
    }
    if (type1 == "QG") {
      getspliti_ = nf * std::min(y, 10.0);
    } else {
      splitintqg = nf * std::min(y, 10.0);
    }
  }
  if (type1 == "QQ") {
    for (int i = 1; i <= 2; ++i) {
      x1a(i) = qval(qline - 1 + i);
      x2a(i) = zmval(zline - 1 + i);
      for (int j = 1; j <= 2; ++j) {
        ya(i, j) = splitiqqv(qline - 1 + i, zline - 1 + j);
      }
    }
    for (int i = 1; i <= 2; ++i) {
      a = (ya(i, 2) - ya(i, 1)) / (x2a(2) - x2a(1));
      b = ya(i, 1) - a * x2a(1);
      yb(i) = a * low + b;
    }
    if (x1a(1) == x1a(2)) {
      y = (yb(1) + yb(2)) / 2.0;
    } else {
      a = (yb(2) - yb(1)) / (x1a(2) - x1a(1));
      b = yb(1) - a * x1a(1);
      y = a * qb + b;
    }
    getspliti_ = std::min(y, 10.0);
  }
  if (type1 == "GC") getspliti_ = splitintgg + splitintqg;
  return getspliti_;
}

// ***********************************************************************
// ***	  function getinspliti
// ***********************************************************************
double getinspliti(const double& qb, const std::string& type1) {
  // C--Parameter common block
  using namespace cb::param;
  // C--local variables
  constexpr double pi = 3.141592653589793;
  double low, y;
  // NOTE(port): possibly uninitialised in Fortran (result for an unknown TYPE1)
  double getinspliti_ = 0.0;
  double splitintgg = 0.0, splitintqg = 0.0;
  // (33.-2.*NF) is REAL*4 in Fortran
  const double b0 = static_cast<double>(33.f - 2.f * static_cast<float>(nf));

  // C--find boundaries for z integration
  const double up = 1.0 - ipow<2>(q0) / (4.0 * ipow<2>(qb));
  if ((type1 == "GG") || (type1 == "GC")) {
    low = 1.0 - up;
    if (up <= low) {
      getinspliti_ = 0.0;
      return getinspliti_;
    }
    const double lnlow = std::log((1.0 - low) * ipow<2>(qb) / ipow<2>(lps));
    const double lnup = std::log((1.0 - up) * ipow<2>(qb) / ipow<2>(lps));
    const double e1 = ei(lnlow);
    const double e2 = ei(2.0 * lnlow);
    const double e3 = ei(3.0 * lnlow);
    const double e4 = ei(lnup);
    const double e5 = ei(2.0 * lnup);
    const double e6 = ei(3.0 * lnup);
    y = 2.0 *
        (std::log(lnlow) - ipow<2>(lps) * e1 / ipow<2>(qb) + ipow<4>(lps) * e2 / ipow<4>(qb) -
         ipow<6>(lps) * e3 / ipow<6>(qb) - std::log(lnup) + ipow<2>(lps) * e4 / ipow<2>(qb) -
         ipow<4>(lps) * e5 / ipow<4>(qb) + ipow<6>(lps) * e6 / ipow<6>(qb) + low - std::log(low) -
         up + std::log(up)) *
        3.0 * 12.0 * pi / (2.0 * pi * b0);
    if (type1 == "GG") {
      getinspliti_ = y;
    } else {
      splitintgg = y;
    }
  }
  if ((type1 == "QG") || (type1 == "GC")) {
    low = 0.0;
    if (up <= low) {
      getinspliti_ = 0.0;
      return getinspliti_;
    }
    const double lnlow = std::log((1.0 - low) * ipow<2>(qb) / ipow<2>(lps));
    const double lnup = std::log((1.0 - up) * ipow<2>(qb) / ipow<2>(lps));
    const double e1 = ei(3.0 * lnlow);
    const double e2 = ei(2.0 * lnlow);
    const double e3 = ei(lnlow);
    const double e4 = ei(3.0 * lnup);
    const double e5 = ei(2.0 * lnup);
    const double e6 = ei(lnup);
    y = (2.0 * ipow<6>(lps) * e1 / ipow<6>(qb) - 2.0 * ipow<4>(lps) * e2 / ipow<4>(qb) +
         2.0 * ipow<2>(lps) * e3 / ipow<2>(qb) - 2.0 * ipow<6>(lps) * e4 / ipow<6>(qb) +
         2.0 * ipow<4>(lps) * e5 / ipow<4>(qb) - 2.0 * ipow<2>(lps) * e6 / ipow<2>(qb)) *
        12.0 * pi / (2.0 * 2.0 * pi * b0);
    if (type1 == "QG") {
      getinspliti_ = nf * y;
    } else {
      splitintqg = nf * y;
    }
  }
  if (type1 == "QQ") {
    low = 0.0;
    if (up <= low) {
      getinspliti_ = 0.0;
      return getinspliti_;
    }
    const double lnlow = std::log((1.0 - low) * ipow<2>(qb) / ipow<2>(lps));
    const double lnup = std::log((1.0 - up) * ipow<2>(qb) / ipow<2>(lps));
    const double e1 = ei(lnlow);
    const double e2 = ei(2.0 * lnlow);
    const double e3 = ei(lnup);
    const double e4 = ei(2.0 * lnup);
    y = (2.0 * std::log(lnlow) - 2.0 * ipow<2>(lps) * e1 / ipow<2>(qb) +
         ipow<4>(lps) * e2 / ipow<4>(qb) - 2.0 * std::log(lnup) +
         2.0 * ipow<2>(lps) * e3 / ipow<2>(qb) - ipow<4>(lps) * e4 / ipow<4>(qb)) *
        4.0 * 12.0 * pi / (3.0 * 2.0 * pi * b0);
    getinspliti_ = y;
  }
  if (type1 == "GQ") {
    low = 1.0 - up;
    if (up <= low) {
      getinspliti_ = 0.0;
      return getinspliti_;
    }
    y = (ipow<2>(up) / 2.0 - 2.0 * up + 2.0 * std::log(up) - ipow<2>(low) / 2.0 + 2.0 * low -
         2.0 * std::log(low)) *
        4.0 * 12.0 * pi / (3.0 * 2.0 * pi * b0 * std::log(ipow<2>(qb) / ipow<2>(lps)));
    getinspliti_ = y;
  }
  if (type1 == "GC") getinspliti_ = splitintgg + splitintqg;
  return getinspliti_;
}

// ***********************************************************************
// ***	  function getpdf
// ***********************************************************************
double getpdf(const double& x, const double& q, const std::string& typ) {
  // C--pdf common block
  using namespace cb::pdfs;
  // C--Parameter common block
  using namespace cb::param;
  // C--variables for pdf integration
  using namespace cb::pdfintv;
  // C--local variables
  constexpr double epsi = 1.e-4;  // DATA EPSI/1.d-4/
  double qlow, qhigh, ystart, hfirst;
  double getpdf_ = 0.0;

  if ((x < 0.0) || (x > 1.0) || (q < q0)) {
    io::log(" error in GETPDF: parameter out of bound {} {}", x, q);
    getpdf_ = 0.0;
    return getpdf_;
  }

  if (typ == "QQ") {
    z = x;
    xmax = q;
    // C--f_q^q
    qlow = std::max(q0, q0 / (2.0 * std::sqrt(1.0 - x)));
    qhigh = q;
    if ((qlow >= qhigh * (1.0 - 1.e-10)) || (x > 1.0 - 1.e-10)) {
      ystart = 0.0;
    } else {
      hfirst = r4(0.01f) * (qhigh - qlow);
      ystart = 0.0;
      odeint(ystart, qlow, qhigh, epsi, hfirst, 0.0, 7);
    }
    getpdf_ = ystart;
  } else if (typ == "GQ") {
    z = x;
    xmax = q;
    // C--f_q^g
    qlow = std::max(q0, std::max(q0 / (2.0 * std::sqrt(x)), q0 / (2.0 * std::sqrt(1.0 - x))));
    qhigh = q;
    if ((qlow >= qhigh * (1.0 - 1.e-10)) || (x < 0.0 + 1.e-10) || (x > 1.0 - 1.e-10)) {
      ystart = 0.0;
    } else {
      hfirst = r4(0.01f) * (qhigh - qlow);
      ystart = 0.0;
      odeint(ystart, qlow, qhigh, epsi, hfirst, 0.0, 8);
    }
    getpdf_ = ystart;
  } else if (typ == "QG") {
    z = x;
    xmax = q;
    // C--f_q^g
    qlow = std::max(q0, q0 / (2.0 * std::sqrt(1.0 - x)));
    qhigh = q;
    if ((qlow >= qhigh * (1.0 - 1.e-10)) || (x > 1.0 - 1.e-10)) {
      ystart = 0.0;
    } else {
      hfirst = r4(0.01f) * (qhigh - qlow);
      ystart = 0.0;
      odeint(ystart, qlow, qhigh, epsi, hfirst, 0.0, 9);
    }
    getpdf_ = ystart;
  } else if (typ == "GG") {
    z = x;
    xmax = q;
    // C--f_q^q
    qlow = std::max(q0, std::max(q0 / (2.0 * std::sqrt(x)), q0 / (2.0 * std::sqrt(1.0 - x))));
    qhigh = q;
    if ((qlow >= qhigh * (1.0 - 1.e-10)) || (x < 0.0 + 1.e-10) || (x > 1.0 - 1e-10)) {
      ystart = 0.0;
    } else {
      hfirst = r4(0.01f) * (qhigh - qlow);
      ystart = 0.0;
      odeint(ystart, qlow, qhigh, epsi, hfirst, 0.0, 10);
    }
    getpdf_ = ystart;
  } else {
    io::log(" error: pdf-type {} does not exist", typ);
    getpdf_ = 0.0;
  }
  return getpdf_;
}

// ***********************************************************************
// ***	  function getpdfxint
// ***********************************************************************
double getpdfxint(const double& q, const std::string& typ) {
  // C--pdf common block
  using namespace cb::pdfs;
  // C--number of extrapolations in tables
  using namespace cb::extrapolations;
  // C--local variables
  // NOTE(port): possibly uninitialised in Fortran (XA, YA for an unknown TYP)
  F1<double, 2> xa, ya;

  ntotpdf = ntotpdf + 1;
  if (ipow<2>(q) > qinqx(1, 1000)) {
    noverpdf = noverpdf + 1;
    if (noverpdf <= 25)
      io::log(" WARNING in getpdfxint: need to extrapolate:  {} {}", ipow<2>(q), qinqx(1, 1000));
  }

  const int q2close = fint((std::log(ipow<2>(q)) - std::log(qinqx(1, 1))) * 999.0 /
                               (std::log(qinqx(1, 1000)) - std::log(qinqx(1, 1))) +
                           1);
  int q2line = std::max(q2close, 1);
  q2line = std::min(q2line, 999);
  if ((q2line > 999) || (q2line < 1)) {
    io::log(" ERROR in GETPDFXINT: line number out of bound {}", q2line);
  }

  if (typ == "QQ") {
    for (int j = 1; j <= 2; ++j) {
      xa(j) = qinqx(1, q2line - 1 + j);
      ya(j) = qinqx(2, q2line - 1 + j);
    }
  } else if (typ == "GQ") {
    for (int j = 1; j <= 2; ++j) {
      xa(j) = ginqx(1, q2line - 1 + j);
      ya(j) = ginqx(2, q2line - 1 + j);
    }
  } else if (typ == "QG") {
    for (int j = 1; j <= 2; ++j) {
      xa(j) = qingx(1, q2line - 1 + j);
      ya(j) = qingx(2, q2line - 1 + j);
    }
  } else if (typ == "GG") {
    for (int j = 1; j <= 2; ++j) {
      xa(j) = gingx(1, q2line - 1 + j);
      ya(j) = gingx(2, q2line - 1 + j);
    }
  } else {
    io::log(" error in GETPDFXINT: unknown integral type {}", typ);
  }
  const double a = (ya(2) - ya(1)) / (xa(2) - xa(1));
  const double b = ya(1) - a * xa(1);
  const double y = a * ipow<2>(q) + b;
  return y;
}

// ***********************************************************************
// ***	  subroutine getpdfxintexact
// ***********************************************************************
double getpdfxintexact(const double& q, const std::string& typ) {
  // C--Parameter common block
  using namespace cb::param;
  // C--variables for pdf integration
  using namespace cb::pdfintv;
  // C--local variables
  constexpr double epsi = 1.e-4;  // DATA EPSI/1.d-4/

  const double hfirst = 0.01;
  double ystart = 0.0;
  xmax = q;
  z = 0.0;
  if (typ == "QQ") {
    odeint(ystart, q0, q, epsi, hfirst, 0.0, 21);
  } else if (typ == "QG") {
    odeint(ystart, q0, q, epsi, hfirst, 0.0, 23);
  } else if (typ == "GQ") {
    odeint(ystart, q0, q, epsi, hfirst, 0.0, 22);
  } else if (typ == "GG") {
    odeint(ystart, q0, q, epsi, hfirst, 0.0, 24);
  }
  return ystart;
}

// ***********************************************************************
// ***	  function getxsecint
// ***********************************************************************
double getxsecint(const double& tm, const double& md, const std::string& typ2) {
  // C--Parameter common block
  using namespace cb::param;
  // C--cross secttion common block
  using namespace cb::xsecs;
  // C--variables for cross section integration
  using namespace cb::xsecv;
  // C--number of extrapolations in tables
  using namespace cb::extrapolations;
  // C--local variables
  // NOTE(port): possibly uninitialised in Fortran (X1A, X2A, YA for an unknown TYP2)
  F1<double, 2> x1a, x2a, yb;
  F2<double, 2, 2> ya;
  double y, a, b;

  ntotxsec = ntotxsec + 1;
  if (tm > intq1(1000, 101)) {
    noverxsec = noverxsec + 1;
    // NOTE(port): the original tests noverpdf (not noverxsec) here; kept.
    if (noverpdf <= 25)
      io::log(" WARNING in getxsecint: need to extrapolate:  {} {}", tm, intq1(1000, 101));
  }

  const int tclose = fint((std::log(tm) - std::log(intq1(1, 101))) * 999.0 /
                              (std::log(intq1(1000, 101)) - std::log(intq1(1, 101))) +
                          1);
  int tline = std::max(tclose, 1);
  tline = std::min(tline, 999);
  const int mdclose = fint((md - intq1(1001, 1)) * 99.0 / (intq1(1001, 100) - intq1(1001, 1)) + 1);
  int mdline = std::max(mdclose, 1);
  mdline = std::min(mdline, 99);
  if ((tline > 999) || (mdline > 99) || (tline < 1) || (mdline < 1)) {
    io::log(" ERROR in GETXSECINT: line number out of bound {} {}", tline, mdline);
  }

  if (typ2 == "QA") {
    // C--first quark integral
    for (int i = 1; i <= 2; ++i) {
      x1a(i) = intq1(1001, mdline - 1 + i);
      x2a(i) = intq1(tline - 1 + i, 101);
      for (int j = 1; j <= 2; ++j) {
        ya(i, j) = intq1(tline - 1 + j, mdline - 1 + i);
      }
    }
  } else if (typ2 == "QB") {
    // C--second quark integral
    for (int i = 1; i <= 2; ++i) {
      x1a(i) = intq2(1001, mdline - 1 + i);
      x2a(i) = intq2(tline - 1 + i, 101);
      for (int j = 1; j <= 2; ++j) {
        ya(i, j) = intq2(tline - 1 + j, mdline - 1 + i);
      }
    }
  } else if (typ2 == "GA") {
    // C--first gluon integral
    for (int i = 1; i <= 2; ++i) {
      x1a(i) = intg1(1001, mdline - 1 + i);
      x2a(i) = intg1(tline - 1 + i, 101);
      for (int j = 1; j <= 2; ++j) {
        ya(i, j) = intg1(tline - 1 + j, mdline - 1 + i);
      }
    }
  } else if (typ2 == "GB") {
    // C--second gluon integral
    for (int i = 1; i <= 2; ++i) {
      x1a(i) = intg2(1001, mdline - 1 + i);
      x2a(i) = intg2(tline - 1 + i, 101);
      for (int j = 1; j <= 2; ++j) {
        ya(i, j) = intg2(tline - 1 + j, mdline - 1 + i);
      }
    }
  } else {
    io::log(" error in GETXSECINT: unknown integral type {}", typ2);
  }
  for (int i = 1; i <= 2; ++i) {
    a = (ya(i, 2) - ya(i, 1)) / (x2a(2) - x2a(1));
    b = ya(i, 1) - a * x2a(1);
    yb(i) = a * tm + b;
  }
  if (x1a(1) == x1a(2)) {
    y = yb(1);
  } else {
    a = (yb(2) - yb(1)) / (x1a(2) - x1a(1));
    b = yb(1) - a * x1a(1);
    y = a * md + b;
  }
  return y;
}

// ***********************************************************************
// ***	  function getinsudafast
// ***********************************************************************
double getinsudafast(const double& q1, const double& q2, const std::string& typ) {
  // C--Parameter common block
  using namespace cb::param;
  // C--local variables
  double getinsudafast_;

  if (q2 <= q1) {
    getinsudafast_ = 1.0;
  } else if (q1 <= q0) {
    getinsudafast_ = getinsudared(q2, typ);
  } else {
    const double s2 = getinsudared(q2, typ);
    getinsudafast_ = s2 / getinsudared(q1, typ);
  }
  if (getinsudafast_ > 1.0) getinsudafast_ = 1.0;
  if (getinsudafast_ < (-1.e-10)) {
    io::log(" ERROR: GETINSUDAFAST < 0: {} for {}   {}   {}", getinsudafast_, q1, q2, typ);
  }
  if (getinsudafast_ < 0.0) getinsudafast_ = 0.0;
  return getinsudafast_;
}

// ***********************************************************************
// ***	  function getinsudared
// ***********************************************************************
double getinsudared(const double& q, const std::string& typ2) {
  // C--Parameter common block
  using namespace cb::param;
  // C--Sudakov common block
  using namespace cb::insuda;
  // C--number of extrapolations in tables
  using namespace cb::extrapolations;
  // C--local variables
  // NOTE(port): possibly uninitialised in Fortran (XA, YA for an unknown TYP2)
  F1<double, 2> xa, ya;

  ntotsuda = ntotsuda + 1;
  if (q > sudaqq(1000, 1)) {
    noversuda = noversuda + 1;
    if (noversuda <= 25)
      io::log(" WARNING in getinsudared: need to extrapolate:  {} {}", q, sudaqq(1000, 1));
  }

  const int qclose = fint((std::log(q) - std::log(sudaqq(1, 1))) * 999.0 /
                              (std::log(sudaqq(1000, 1)) - std::log(sudaqq(1, 1))) +
                          1);
  int qbin = std::max(qclose, 1);
  qbin = std::min(qbin, 999);
  if ((qbin > 999) || (qbin < 1)) {
    io::log(" ERROR in GETINSUDARED: line number out of bound {}", qbin);
  }
  if (typ2 == "QQ") {
    for (int i = 1; i <= 2; ++i) {
      xa(i) = sudaqq(qbin - 1 + i, 1);
      ya(i) = sudaqq(qbin - 1 + i, 2);
    }
  } else if (typ2 == "QG") {
    for (int i = 1; i <= 2; ++i) {
      xa(i) = sudaqg(qbin - 1 + i, 1);
      ya(i) = sudaqg(qbin - 1 + i, 2);
    }
  } else if (typ2 == "GG") {
    for (int i = 1; i <= 2; ++i) {
      xa(i) = sudagg(qbin - 1 + i, 1);
      ya(i) = sudagg(qbin - 1 + i, 2);
    }
  } else if (typ2 == "GC") {
    for (int i = 1; i <= 2; ++i) {
      xa(i) = sudagc(qbin - 1 + i, 1);
      ya(i) = sudagc(qbin - 1 + i, 2);
    }
  } else {
    io::log(" error in GETINSUDARED: unknown type {}", typ2);
  }
  const double a = (ya(2) - ya(1)) / (xa(2) - xa(1));
  const double b = ya(1) - a * xa(1);
  const double y = a * q + b;
  double getinsudared_ = y;
  if (getinsudared_ < (-1.e-10)) {
    io::log(" ERROR: GETINSUDARED < 0: {} {} {}", getinsudared_, q, typ2);
  }
  if (getinsudared_ < 0.0) getinsudared_ = 0.0;
  return getinsudared_;
}

// ***********************************************************************
// ***	  function getsscat
// ***********************************************************************
double getsscat(const double& en, const double& px, const double& py, const double& pz,
                const double& mp, const double& lw, const char& type1, const char& type2,
                const double& x, const double& y, const double& z, const double& t,
                const int& mode) {
  // C--Parameter common block
  using namespace cb::param;
  // C--variables for cross section integration
  using namespace cb::xsecv;
  // C--local variables
  double up, ccol, low, mdeb, pcms2, shat, pproj, psct;
  F1<double, 5> avmom;
  double getsscat_;

  if (type1 == 'Q') {
    ccol = r4(2.f / 3.f);
  } else {
    ccol = 3.0 / 2.0;
  }
  if (mode == 0) {
    mdeb = medium().getmd(x, y, z, t);
    medium().avscatcen(x, y, z, t, avmom(1), avmom(2), avmom(3), avmom(4), avmom(5));
    shat = ipow<2>(avmom(5)) + ipow<2>(mp) +
           2.0 * (avmom(4) * en - avmom(1) * px - avmom(2) * py - avmom(3) * pz);
    if (shat == 0.0) {
      getsscat_ = 0.0;
      return getsscat_;
    }
    pcms2 = ipow<2>(shat + ipow<2>(mp) - ipow<2>(avmom(5))) / (4.0 * shat) - ipow<2>(mp);
    up = 4.0 * pcms2;
  } else {
    if (mode == 1) {
      mdeb = medium().getmdmin();
    } else {
      mdeb = medium().getmdmax();
    }
    medium().maxscatcen(avmom(1), avmom(2), avmom(3), avmom(4), avmom(5));
    psct = std::sqrt(ipow<2>(avmom(1)) + ipow<2>(avmom(2)) + ipow<2>(avmom(3)));
    pproj = std::sqrt(ipow<2>(px) + ipow<2>(py) + ipow<2>(pz));
    shat = ipow<2>(avmom(5)) + ipow<2>(mp) + 2.0 * (en * avmom(4) + pproj * psct);
    pcms2 = ipow<2>(shat + ipow<2>(mp) - ipow<2>(avmom(5))) / (4.0 * shat) - ipow<2>(mp);
    up = 4.0 * pcms2;
  }
  low = ipow<2>(lw);
  if (low > up) {
    getsscat_ = 0.0;
    return getsscat_;
  }
  if ((type2 == 'C') || ((type1 == 'Q') && (type2 == 'Q')) || ((type1 == 'G') && (type2 == 'G'))) {
    const double s1 = scatprimfunc(up, mdeb);
    const double s2 = scatprimfunc(low, mdeb);
    getsscat_ = ccol * (s1 - s2);
  } else {
    getsscat_ = 0.0;
  }
  low = ipow<2>(q0) / ipow<2>(scalefacm);
  if (up > low) {
    if (type1 == 'Q') {
      if ((type2 == 'C') || (type2 == 'G')) {
        const double pdf = getpdfxint(scalefacm * std::sqrt(up), kGQ);
        const double sp = scatprimfunc(up, mdeb);
        getsscat_ = getsscat_ + pdf * 3.0 * sp / 2.0;
        getsscat_ = getsscat_ - getxsecint(up, mdeb, kQA);
      }
    } else {
      if ((type2 == 'C') || (type2 == 'G')) {
        const double s1 = scatprimfunc(up, mdeb);
        const double s2 = scatprimfunc(low, mdeb);
        const double xs = getxsecint(up, mdeb, kGB);
        getsscat_ = getsscat_ + ccol * (s1 - s2) - xs;
      }
      if ((type2 == 'C') || (type2 == 'Q')) {
        const double pdf = getpdfxint(scalefacm * std::sqrt(up), kQG);
        const double sp = scatprimfunc(up, mdeb);
        getsscat_ = getsscat_ + 2.0 * pdf * 2.0 * sp / 3.0;
        getsscat_ = getsscat_ - 2.0 * getxsecint(up, mdeb, kGA);
      }
    }
  }
  if (getsscat_ < -1.e-4) {
    io::log(" error: cross section < 0 {} for {} {} {} {} {} {} {} {}", getsscat_, en, mp, lw,
            type1, type2, ipow<2>(lw), up, mode);
  }
  getsscat_ = std::max(getsscat_, 0.0);
  return getsscat_;
}

// ***********************************************************************
// ***	  function getmass
// ***********************************************************************
double getmass(const double& qbmin, const double& qbmax, const double& theta, const double& ep,
               const std::string& type, const double& max2, const bool& ins, double& zdec,
               bool& qqbardec, const double mq, int* qqfl) {
  // C--Common block of Pythia
  using namespace py::pyjets;
  using namespace py::pydat1;
  using namespace py::pydat3;
  // C--Parameter common block
  using namespace cb::param;
  // C--time common block
  using namespace cb::time;
  // C--factor in front of alphas argument
  using namespace cb::alphasfac;
  // C--local variables
  constexpr double pi = 3.141592653589793;
  double q2min, alphmax, pref, q2max, gmax, arg, cand, eps, trueval, oest, weight, r, rz,
      thetanew, r2, rmin;
  // NOTE(port): possibly uninitialised in Fortran: Z and TRUEEPS are only set in the .NOT.INS
  // branch but are read in the 'weight> 1' message and (Z) in the ANGORD=3/4 theta check.
  double z = 0.0, trueeps = 0.0;
  double getmass_ = 0.0;
  // heavy flavour (HEAVYFLAVOUR T): weights of the gluon channels at the accepted point
  double hfal = 0.0, hfwc = 0.0, hfwb = 0.0;
  if (qqfl) *qqfl = 0;

  q2min = ipow<2>(q0);

  alphmax = alphas(3.0 * ptfac * q2min / 16.0, lps);

  if (type == "QQ") {
    pref = 4.0 * alphmax / (3.0 * 2.0 * pi);
  } else {
    pref = 29.0 * alphmax / (8.0 * 2.0 * pi);
  }

  // C--check if phase space available, return 0.d0 otherwise
  if ((qbmax <= qbmin) || (ep < qbmin)) {
    getmass_ = 0.0;
    zdec = 0.0;
    qqbardec = false;
    return getmass_;
  }

  q2max = ipow<2>(qbmax);
l21:
  if (q2max - ipow<2>(qbmin) < r4(1e-4f)) {
    getmass_ = qbmin;
    zdec = 0.5;
    if (type == "QQ") {
      qqbardec = false;
    } else if (cb::hf::on) {
      const int nl = std::clamp(nf, 1, 3);
      const double r3 = py::pyr(0) * (nl * pqg(0.5) + pgg(0.5));
      qqbardec = (r3 >= pgg(0.5));
      if (qqfl && qqbardec) *qqfl = 1 + std::min(static_cast<int>(nl * (r3 - pgg(0.5)) / (nl * pqg(0.5))), nl - 1);
    } else {
      if (py::pyr(0) < pqg(0.5) / (pqg(0.5) + pgg(0.5))) {
        qqbardec = true;
      } else {
        qqbardec = false;
      }
    }
    return getmass_;
  }
  gmax = pref * ipow<2>(std::log(q2min / (4.0 * q2max)));
  if (qbmin > 0.0) {
    rmin = std::exp(pref * ipow<2>(std::log(q2min / (4.0 * ipow<2>(qbmin)))) - gmax);
  } else {
    rmin = 0.0;
  }

  r = py::pyr(0) * (1.0 - rmin) + rmin;
  arg = gmax + std::log(r);
  if (arg < 0.0) {
    getmass_ = 0.0;
    zdec = 0.0;
    qqbardec = false;
    return getmass_;
  }
  cand = q2min * std::exp(std::sqrt(arg / pref)) / 4.0;
  eps = q2min / (4.0 * cand);

  if ((cand < q2min) || (cand < ipow<2>(qbmin))) {
    getmass_ = 0.0;
    zdec = 0.0;
    qqbardec = false;
    return getmass_;
  }

  if ((cand > ipow<2>(max2)) || (cand + mq * mq > ipow<2>(ep))) {
    q2max = cand;
    goto l21;
  }

  if (ins) {
    trueval = getinspliti(std::sqrt(cand), type);
    oest = -2.0 * pref * std::log(eps);
    weight = trueval / oest;
  } else if (cb::hf::on) {
    // heavy-flavour shower: massive Q -> Q g, and g -> c cbar, b bbar added to the gluon
    trueeps = 0.5 - 0.5 * std::sqrt(1.0 - q2min / cand) * std::sqrt(1.0 - cand / ipow<2>(ep));
    rz = py::pyr(0);
    z = 1.0 - std::pow(eps, rz);
    const double asmin = 3.0 * ptfac * q2min / 16.0;  // scale of alphmax
    if (type == "QQ") {
      double zlo = trueeps, zhi = 1.0 - trueeps;
      if (mq > 0.0 && !hfzlimitsqq(cand, mq, ep, zlo, zhi)) zlo = zhi = 2.0;
      if ((z < zlo) || (z > zhi)) {
        weight = 0.0;
      } else if (mq > 0.0) {
        trueval = alphas(std::max(ptfac * hfkt2qq(cand, mq, z), asmin), lps) *
                  std::max(pqqmassive(z, cand, mq), 0.0) / (2.0 * pi);
        weight = trueval / (2.0 * pref / (1.0 - z));
      } else {
        trueval = alphas(ptfac * z * (1.0 - z) * cand, lps) * pqq(z) / (2.0 * pi);
        weight = trueval / (2.0 * pref / (1.0 - z));
      }
    } else {
      if (py::pyr(0) < r4(17.f / 29.f)) z = 1.0 - z;
      if ((z < trueeps) || (z > (1.0 - trueeps))) {
        weight = 0.0;
      } else {
        hfal = alphas(ptfac * z * (1.0 - z) * cand, lps);
        double zlo, zhi;
        hfwc = hfwb = 0.0;
        if (hfzlimitsgq(cand, cb::hf::mc, ep, zlo, zhi) && (z > zlo) && (z < zhi))
          hfwc = alphas(std::max(ptfac * hfkt2gq(cand, cb::hf::mc, z), asmin), lps) *
                 pgqmassive(z, cand, cb::hf::mc);
        if (hfzlimitsgq(cand, cb::hf::mb, ep, zlo, zhi) && (z > zlo) && (z < zhi))
          hfwb = alphas(std::max(ptfac * hfkt2gq(cand, cb::hf::mb, z), asmin), lps) *
                 pgqmassive(z, cand, cb::hf::mb);
        // light g -> q qbar with the physical weight N_f T_R (JEWEL 2.6.0: one flavour's worth, PQG)
        trueval = (hfal * (pgg(z) + std::clamp(nf, 1, 3) * pqg(z)) + hfwc + hfwb) / (2.0 * pi);
        oest = alphmax * (17.0 / (4.0 * z) + 3.0 / (1.0 - z)) / (2.0 * pi);
        weight = trueval / oest;
      }
    }
  } else {
    // C--find true z interval
    trueeps = 0.5 - 0.5 * std::sqrt(1.0 - q2min / cand) * std::sqrt(1.0 - cand / ipow<2>(ep));
    if (trueeps < eps) io::log(" error in getmass: true eps < eps {} {}", trueeps, eps);
    rz = py::pyr(0);
    z = 1.0 - std::pow(eps, rz);
    if ((z < trueeps) || (z > (1.0 - trueeps))) {
      weight = 0.0;
    } else {
      if (type == "QQ") {
        trueval = alphas(ptfac * z * (1.0 - z) * cand, lps) * pqq(z) / (2.0 * pi);
        oest = 2.0 * pref / (1.0 - z);
        weight = trueval / oest;
      } else {
        if (py::pyr(0) < r4(17.f / 29.f)) z = 1.0 - z;
        trueval = alphas(ptfac * z * (1.0 - z) * cand, lps) * (pgg(z) + pqg(z)) / (2.0 * pi);
        oest = alphmax * (17.0 / (4.0 * z) + 3.0 / (1.0 - z)) / (2.0 * pi);
        weight = trueval / oest;
      }
    }
  }
  if (weight > 1.0)
    io::log(" problem in getmass: weight> 1 {} {} {} {} {} {} {}", weight, type, eps, trueeps, z,
            cand, ep);
  r2 = py::pyr(0);
  if (r2 > weight) {
    q2max = cand;
    goto l21;
  } else {
    if ((angord == 3) || (angord == 4)) {
      thetanew = (mq > 0.0) ? hftheta(cand, mq, z, ep) : std::sqrt(cand / (z * (1.0 - z))) / ep;
      if ((theta > 0.0) && (thetanew > theta)) {
        q2max = cand;
        goto l21;
      }
    }
    getmass_ = std::sqrt(cand);
    if (!ins) {
      zdec = z;
      if (type == "QQ") {
        qqbardec = false;
      } else if (cb::hf::on) {
        // choose g -> g g, q qbar (light), c cbar, b bbar according to their weights
        const double wgg = hfal * pgg(z), wqq = hfal * std::clamp(nf, 1, 3) * pqg(z);
        const double r3 = py::pyr(0) * (wgg + wqq + hfwc + hfwb);
        qqbardec = (r3 >= wgg);
        if (qqfl) {
          if (r3 < wgg + wqq) {
            // light g -> q qbar: d, u, s with equal weight (massless), from the same random number
            const int nl = std::clamp(nf, 1, 3);
            *qqfl = (r3 >= wgg) ? 1 + std::min(static_cast<int>(nl * (r3 - wgg) / wqq), nl - 1) : 0;
          } else {
            *qqfl = (r3 < wgg + wqq + hfwc) ? 4 : 5;
          }
        }
      } else {
        if (py::pyr(0) < pqg(z) / (pqg(z) + pgg(z))) {
          qqbardec = true;
        } else {
          qqbardec = false;
        }
      }
    }
  }
  return getmass_;
}

// ***********************************************************************
// ***	  function generatez
// ***********************************************************************
double generatez(const double& ti, const double& ea, const double& epsi, const std::string& type) {
  // C--Parameter common block
  using namespace cb::param;
  // C--local variables
  double eps, x, r, help, r1;
  double generatez_ = 0.0;

  if (ti == 0.0) {
    eps = epsi;
  } else {
    eps = std::max(0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ti) * std::sqrt(1.0 - ti / ipow<2>(ea)),
                   epsi);
  }
  if (eps > 0.5) {
    generatez_ = 0.5;
    return generatez_;  // GOTO 61
  }
  for (;;) {  // 60
    r = py::pyr(0);
    if (type == "QQ") {
      x = 1.0 - (1.0 - eps) * std::pow(eps / (1.0 - eps), r);
      r = py::pyr(0);
      if (r < ((1.0 + ipow<2>(x)) / 2.0)) {
        generatez_ = x;
        break;
      }
    } else if (type == "GG") {
      x = 1.0 / (1.0 + std::pow((1.0 - eps) / eps, 1.0 - 2.0 * r));
      r = py::pyr(0);
      help = ((1.0 - x) / x + x / (1.0 - x) + x * (1.0 - x)) / (1.0 / (1.0 - x) + 1.0 / x);
      if (r < help) {
        generatez_ = x;
        break;
      }
    } else {
      r = py::pyr(0) * (1.0 - 2.0 * eps) + eps;
      r1 = py::pyr(0) / 2.0;
      help = 0.5 * (ipow<2>(r) + ipow<2>(1.0 - r));
      if (r1 < help) {
        generatez_ = r;
        break;
      }
    }
  }
  return generatez_;  // 61
}

// ***********************************************************************
// ***	  function scatprimfunc
// ***********************************************************************
double scatprimfunc(const double& t, const double& mdeb) {
  // C--Parameter common block
  using namespace cb::param;
  // C--local variables
  constexpr double pi = 3.141592653589793;
  // (33.-2.*NF)**2 is REAL*4 in Fortran
  const float b0f = 33.f - 2.f * static_cast<float>(nf);
  const double b0sq = static_cast<double>(b0f * b0f);

  const double e = ei(-std::log((t + ipow<2>(mdeb)) / ipow<2>(lqcd)));
  return 2.0 * pi * ipow<2>(12.0 * pi) *
         (-(e / ipow<2>(lqcd)) -
          1.0 / ((t + ipow<2>(mdeb)) * std::log((t + ipow<2>(mdeb)) / ipow<2>(lqcd)))) /
         b0sq;
}

// ***********************************************************************
// ***	  function intpqq
// ***********************************************************************
double intpqq(const double& z, const double& q) {
  // C--Parameter common block
  using namespace cb::param;
  // ((33.-2.*NF)*3.) is REAL*4 in Fortran
  const double d = static_cast<double>((33.f - 2.f * static_cast<float>(nf)) * 3.f);

  return 6.0 * 4.0 * (-2.0 * std::log(std::log(ipow<2>(q) / ipow<2>(lps)) + std::log(1.0 - z))) / d;
}

// ***********************************************************************
// ***	  function intpgglow
// ***********************************************************************
double intpgglow(const double& z, const double& q) {
  // C--Parameter common block
  using namespace cb::param;
  // (33.-2.*NF) is REAL*4 in Fortran
  const double d = static_cast<double>(33.f - 2.f * static_cast<float>(nf));

  return 6.0 * 3.0 * (std::log(std::log(ipow<2>(q) / ipow<2>(lps)) + std::log(z))) / d;
}

// ***********************************************************************
// ***	  function intpgghigh
// ***********************************************************************
double intpgghigh(const double& z, const double& q) {
  // C--Parameter common block
  using namespace cb::param;
  // (33.-2.*NF) is REAL*4 in Fortran
  const double d = static_cast<double>(33.f - 2.f * static_cast<float>(nf));

  return -(6.0 * 3.0 * (std::log(std::log(ipow<2>(q) / ipow<2>(lps)) + std::log(1.0 - z))) / d);
}

// ***********************************************************************
// ***	  function intpqglow
// ***********************************************************************
double intpqglow(const double& z, const double& q) {
  // C--Parameter common block
  using namespace cb::param;
  // ((33.-2.*NF)*2.) is REAL*4 in Fortran
  const double d = static_cast<double>((33.f - 2.f * static_cast<float>(nf)) * 2.f);

  const double l = std::log(ipow<2>(q) / ipow<2>(lps)) + std::log(z);
  const double e1 = ei(l);
  const double e2 = ei(2.0 * l);
  const double e3 = ei(3.0 * l);
  return 6.0 *
         (ipow<2>(lps) * e1 / ipow<2>(q) - 2.0 * ipow<4>(lps) * e2 / ipow<4>(q) +
          2.0 * ipow<6>(lps) * e3 / ipow<6>(q)) /
         d;
}

// ***********************************************************************
// ***	  function intpqghigh
// ***********************************************************************
double intpqghigh(const double& z, const double& q) {
  // C--Parameter common block
  using namespace cb::param;
  // ((33.-2.*NF)*2.) is REAL*4 in Fortran
  const double d = static_cast<double>((33.f - 2.f * static_cast<float>(nf)) * 2.f);

  const double l = std::log(ipow<2>(q) / ipow<2>(lps)) + std::log(1.0 - z);
  const double e1 = ei(l);
  const double e2 = ei(2.0 * l);
  const double e3 = ei(3.0 * l);
  return -(6.0 *
           (ipow<2>(lps) * e1 / ipow<2>(q) - 2.0 * ipow<4>(lps) * e2 / ipow<4>(q) +
            2.0 * ipow<6>(lps) * e3 / ipow<6>(q)) /
           d);
}

// ***********************************************************************
// ***	  function gett
// ***********************************************************************
double gett(const double& mint, const double& maxt, const double& mdeb) {
  // C--Parameter common block
  using namespace cb::param;
  // C--local variables
  double tmin, tmax, r1, r2, t;

  tmax = maxt + ipow<2>(mdeb);
  tmin = mint + ipow<2>(mdeb);
  if (tmin > tmax) {
    return 0.0;
  }
  for (;;) {  // 20
    r1 = py::pyr(0);
    t = tmax * tmin / (tmax + r1 * (tmin - tmax));
    r2 = py::pyr(0);
    if (r2 < ipow<2>(alphas(t, lqcd)) / ipow<2>(alphas(tmin, lqcd))) {
      return t - ipow<2>(mdeb);
    }
  }
}

// ***********************************************************************
// ***	  function ei
// ***********************************************************************
double ei(const double& x) {
  // C--exponential integral for negative arguments
  using namespace cb::expint;
  // C--local variables
  int lmax, line;
  F1<double, 2> xa, ya;
  double y, a, b;
  // DATA EPSI/1.e-5/ (unused)

  if (std::fabs(x) > valmax) io::log(" warning: value out of array in Ei(x) {} {}", x, valmax);

  if (x >= 0.0) {
    lmax = fint(x * nval / valmax);
    line = std::max(lmax, 1);
    line = std::min(line, 999);
    if ((line > 999) || (line < 1)) {
      io::log(" ERROR in EI: line number out of bound {}", line);
    }
    for (int k = 1; k <= 2; ++k) {
      xa(k) = eixs(1, line - 1 + k);
      ya(k) = eixs(3, line - 1 + k);
    }
    a = (ya(2) - ya(1)) / (xa(2) - xa(1));
    b = ya(1) - a * xa(1);
    y = a * x + b;
  } else {
    lmax = fint(-x * nval / valmax);
    line = std::max(lmax, 1);
    line = std::min(line, 999);
    if ((line > 999) || (line < 1)) {
      io::log(" ERROR in EI: line number out of bound {}", line);
    }
    for (int k = 1; k <= 2; ++k) {
      xa(k) = eixs(1, line - 1 + k);
      ya(k) = eixs(2, line - 1 + k);
    }
    a = (ya(2) - ya(1)) / (xa(2) - xa(1));
    b = ya(1) - a * xa(1);
    y = -(a * x) + b;
  }
  return y;
}

// ***********************************************************************
// ***	  function pqq
// ***********************************************************************
double pqq(const double& z) { return 4.0 * (1.0 + ipow<2>(z)) / (3.0 * (1.0 - z)); }

// ***********************************************************************
// ***	  function pgq
// ***********************************************************************
double pgq(const double& z) { return 4.0 * (1.0 + ipow<2>(1.0 - z)) / (3.0 * z); }

// ***********************************************************************
// ***	  function pgg
// ***********************************************************************
double pgg(const double& z) { return 3.0 * ((1.0 - z) / z + z / (1.0 - z) + z * (1.0 - z)); }

// ***********************************************************************
// ***	  function pqg
// ***********************************************************************
double pqg(const double& z) { return 0.5 * (ipow<2>(z) + ipow<2>(1.0 - z)); }

// ***********************************************************************
// ***	  function alphas
// ***********************************************************************
double alphas(const double& t, const double& lambda) {
  // C--Parameter common block
  using namespace cb::param;
  // C--local variables
  constexpr double pi = 3.141592653589793;
  // (11.-2.*NF/3.) is REAL*4 in Fortran
  const double b0 = static_cast<double>(11.f - 2.f * static_cast<float>(nf) / 3.f);

  return 4.0 * pi / (b0 * std::log(t / ipow<2>(lambda)));
}

// ***********************************************************************
// ***	  subroutine splitfncint
// ***********************************************************************
void splitfncint(const double& emax) {
  // C--Parameter common block
  using namespace cb::param;
  // C--splitting integral
  using namespace cb::splitint;
  // C--variables for splitting function integration
  using namespace cb::intsplitf;
  // C--max rapidity
  using namespace cb::rapmax;
  // C--local variables
  constexpr double zmmax = 0.5;   // DATA ZMMAX/0.5/
  constexpr int nstep = 999;      // DATA NSTEP/999/
  constexpr double epsi = 1.e-5;  // DATA EPSI/1.d-5/
  double hfirst, ystart, lnzmmin, lnzmmax, zm, zm2, q, shat, pcms2;
  F1<double, 5> avmom;

  medium().maxscatcen(avmom(1), avmom(2), avmom(3), avmom(4), avmom(5));
  shat = ipow<2>(avmom(5)) +
         2.0 * emax * (avmom(4) + std::sqrt(ipow<2>(avmom(1)) + ipow<2>(avmom(2)) + ipow<2>(avmom(3))));
  pcms2 = ipow<2>(shat - ipow<2>(avmom(5))) / (4.0 * shat);
  qmax = std::sqrt(scalefacm * 4.0 * pcms2);

  zmmin = q0 / emax;

  lnzmmin = std::log(zmmin);
  lnzmmax = std::log(zmmax);

  npoint = nstep;

  for (int i = 1, i_end = nstep + 1; i <= i_end; ++i) {
    q = (i - 1) * (qmax - q0) / nstep + q0;
    qval(i) = q;
    qquad = ipow<2>(q);
    for (int j = 1, j_end = nstep + 1; j <= j_end; ++j) {
      zm = std::exp((j - 1) * (lnzmmax - lnzmmin) / nstep + lnzmmin);
      zmval(j) = zm;
      if (ipow<2>(q) < ipow<2>(q0)) {
        zm2 = 0.5;
      } else {
        zm2 = 0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(q));
      }
      zm = std::max(zm, zm2);
      if (zm == 0.5) {
        splitiqqv(i, j) = 0.0;
        splitiggv(i, j) = 0.0;
        splitiqgv(i, j) = 0.0;
      } else {
        ystart = 0.0;
        hfirst = r4(0.01f);
        fm = 0.0;
        odeint(ystart, zm, 1.0 - zm, epsi, hfirst, 0.0, 2);
        splitiqqv(i, j) = ystart;
        ystart = 0.0;
        hfirst = r4(0.01f);
        fm = 0.0;
        odeint(ystart, zm, 1.0 - zm, epsi, hfirst, 0.0, 3);
        splitiggv(i, j) = ystart;
        ystart = 0.0;
        hfirst = r4(0.01f);
        fm = 0.0;
        odeint(ystart, zm, 1.0 - zm, epsi, hfirst, 0.0, 4);
        splitiqgv(i, j) = ystart;
      }
    }
  }
}

// ***********************************************************************
// ***	  subroutine pdfint
// ***********************************************************************
void pdfint(const double& emax) {
  // C--Parameter common block
  using namespace cb::param;
  // C--pdf common block
  using namespace cb::pdfs;
  // C--variables for pdf integration
  using namespace cb::pdfintv;
  // C--max rapidity
  using namespace cb::rapmax;
  // C--local variables
  double q2, q2max, deltaq2, shat, pcms2;
  F1<double, 5> avmom;

  medium().maxscatcen(avmom(1), avmom(2), avmom(3), avmom(4), avmom(5));
  shat = ipow<2>(avmom(5)) +
         2.0 * emax * (avmom(4) + std::sqrt(ipow<2>(avmom(1)) + ipow<2>(avmom(2)) + ipow<2>(avmom(3))));
  pcms2 = ipow<2>(shat - ipow<2>(avmom(5))) / (4.0 * shat);
  q2max = scalefacm * 4.0 * pcms2;

  deltaq2 = std::log(q2max) - std::log(ipow<2>(q0));
  qinqx(1, 1) = ipow<2>(q0);
  ginqx(1, 1) = ipow<2>(q0);
  qingx(1, 1) = ipow<2>(q0);
  gingx(1, 1) = ipow<2>(q0);
  qinqx(2, 1) = 0.0;
  ginqx(2, 1) = 0.0;
  qingx(2, 1) = 0.0;
  gingx(2, 1) = 0.0;
  for (int j = 2; j <= 1000; ++j) {
    q2 = std::exp((j - 1) * deltaq2 / 999.0 + std::log(ipow<2>(q0)));
    qinqx(1, j) = q2;
    ginqx(1, j) = q2;
    qingx(1, j) = q2;
    gingx(1, j) = q2;
    qinqx(2, j) = getpdfxintexact(std::sqrt(q2), kQQ);
    ginqx(2, j) = getpdfxintexact(std::sqrt(q2), kGQ);
    qingx(2, j) = getpdfxintexact(std::sqrt(q2), kQG);
    gingx(2, j) = getpdfxintexact(std::sqrt(q2), kGG);
  }
}

// ***********************************************************************
// ***	  subroutine xsecint
// ***********************************************************************
void xsecint(const double& emax) {
  // C--Parameter common block
  using namespace cb::param;
  // C--cross secttion common block
  using namespace cb::xsecs;
  // C--variables for cross section integration
  using namespace cb::xsecv;
  // C--max rapidity
  using namespace cb::rapmax;
  // C--local variables
  constexpr double epsi = 1.e-4;  // DATA EPSI/1.d-4/
  double tmax, tmaxmax, deltatmax, ystart, mdmin, mdmax, deltamd, shat, pcms2;
  double hfirst;
  F1<double, 5> avmom;

  medium().maxscatcen(avmom(1), avmom(2), avmom(3), avmom(4), avmom(5));
  shat = ipow<2>(avmom(5)) +
         2.0 * emax * (avmom(4) + std::sqrt(ipow<2>(avmom(1)) + ipow<2>(avmom(2)) + ipow<2>(avmom(3))));
  pcms2 = ipow<2>(shat - ipow<2>(avmom(5))) / (4.0 * shat);
  tmaxmax = scalefacm * 4.0 * pcms2;
  deltatmax =
      (std::log(tmaxmax) - std::log(ipow<2>(q0) * (1.0 + 1.e-6) / ipow<2>(scalefacm))) / 999.0;
  mdmin = medium().getmdmin();
  mdmax = std::max(mdmin, medium().getmdmax());
  deltamd = (mdmax - mdmin) / 99.0;

  for (int j = 1; j <= 1000; ++j) {
    tmax = std::exp((j - 1) * deltatmax + std::log(ipow<2>(q0) * (1.0 + 1.e-6) / ipow<2>(scalefacm)));
    intq1(j, 101) = tmax;
    intq2(j, 101) = tmax;
    intg1(j, 101) = tmax;
    intg2(j, 101) = tmax;
    for (int k = 1; k <= 100; ++k) {
      mdx = mdmin + (k - 1) * deltamd;
      intq1(1001, k) = mdx;
      intq2(1001, k) = mdx;
      intg1(1001, k) = mdx;
      intg2(1001, k) = mdx;
      if (tmax < ipow<2>(q0) / ipow<2>(scalefacm)) {
        intq1(j, k) = 0.0;
        intq2(j, k) = 0.0;
        intg1(j, k) = 0.0;
        intg2(j, k) = 0.0;
      } else {
        // C--first quark integral
        qlow = q0;
        hfirst = r4(0.01f) * (tmax - ipow<2>(q0) / ipow<2>(scalefacm));
        ystart = 0.0;
        odeint(ystart, ipow<2>(q0) / ipow<2>(scalefacm), tmax, epsi, hfirst, 0.0, 11);
        intq1(j, k) = ystart;
        // C--second quark integral
        qlow = q0;
        hfirst = r4(0.01f) * (tmax - ipow<2>(q0) / ipow<2>(scalefacm));
        ystart = 0.0;
        odeint(ystart, ipow<2>(q0) / ipow<2>(scalefacm), tmax, epsi, hfirst, 0.0, 14);
        intq2(j, k) = ystart;
        // C--first gluon integral
        qlow = q0;
        ystart = 0.0;
        odeint(ystart, ipow<2>(q0) / ipow<2>(scalefacm), tmax, epsi, hfirst, 0.0, 12);
        intg1(j, k) = ystart;
        // C--second gluon integral
        qlow = q0;
        ystart = 0.0;
        odeint(ystart, ipow<2>(q0) / ipow<2>(scalefacm), tmax, epsi, hfirst, 0.0, 13);
        intg2(j, k) = ystart;
      }
    }
  }
}

// ***********************************************************************
// ***	  function insudaint
// ***********************************************************************
void insudaint(const double& emax) {
  // C--Parameter common block
  using namespace cb::param;
  // C--Sudakov common block
  using namespace cb::insuda;
  // C--max rapidity
  using namespace cb::rapmax;
  // C--local variables
  double qmax, q, delta, shat, pcms2;
  F1<double, 5> avmom;

  medium().maxscatcen(avmom(1), avmom(2), avmom(3), avmom(4), avmom(5));
  shat = ipow<2>(avmom(5)) +
         2.0 * emax * (avmom(4) + std::sqrt(ipow<2>(avmom(1)) + ipow<2>(avmom(2)) + ipow<2>(avmom(3))));
  pcms2 = ipow<2>(shat - ipow<2>(avmom(5))) / (4.0 * shat);
  qmax = std::sqrt(scalefacm * 4.0 * pcms2);
  delta = (std::log(3.0 * qmax) - std::log(ipow<2>(q0) * (1.0 + 1.e-6))) / 999.0;
  for (int i = 1; i <= 1000; ++i) {
    q = std::exp((i - 1) * delta + std::log(ipow<2>(q0) * (1.0 + 1.e-6)));
    sudaqq(i, 1) = q;
    sudaqg(i, 1) = q;
    sudagg(i, 1) = q;
    sudagc(i, 1) = q;
    sudaqq(i, 2) = getinsudakov(q0, q, kQQ);
    sudaqg(i, 2) = getinsudakov(q0, q, kQG);
    sudagg(i, 2) = getinsudakov(q0, q, kGG);
    sudagc(i, 2) = getinsudakov(q0, q, kGC);
  }
}

// ***********************************************************************
// ***	  function eixint
// ***********************************************************************
void eixint() {
  // C--exponential integral for negative arguments
  using namespace cb::expint;
  // C-local variables
  constexpr double epsi = 1.e-6;  // DATA EPSI/1.d-6/
  double x, hfirst, ystart, ei_val;

  nval = 1000;
  valmax = 55.0;

  for (int i = 1, i_end = nval; i <= i_end; ++i) {
    x = i * valmax / (nval * 1.0);
    eixs(1, i) = x;
    // C--do negative arguments first
    ystart = 0.0;
    hfirst = r4(0.01f);
    odeint(ystart, x, 1000.0, epsi, hfirst, 0.0, 5);
    eixs(2, i) = -ystart;
    // C--now do the positive arguments
    eix(x, ei_val);
    eixs(3, i) = ei_val;
  }
}

// ***********************************************************************
// ***	  function odeint
// ***********************************************************************
void odeint(double& ystart, const double& a, const double& b, const double& eps, const double& h1,
            const double& hmin, const int& w1) {
  // C--local variables
  constexpr int nmax = 100000;  // data nmax/100000/
  double x, h, y, dydx, yscale;
  // NOTE(port): possibly uninitialised in Fortran: HNEW/HDID stay unset if RKSTEPPER returns
  // early ('step size not significant') in the first step.
  double hdid = 0.0, hnew = 0.0;

  x = a;
  y = ystart;
  h = fsign(h1, b - a);
  for (int nstep = 1; nstep <= nmax; ++nstep) {
    dydx = deriv(x, w1);
    yscale = std::fabs(y) + std::fabs(h * dydx) + r4(1.e-25f);
    if (((x + h - b) * h) > 0.0) h = b - x;
    rkstepper(x, y, dydx, h, hdid, hnew, yscale, eps, w1);
    if ((x - b) * h >= 0) {
      ystart = y;
      return;
    }
    h = hnew;
    if (std::fabs(h) < std::fabs(hmin)) {
      io::log(" Error in odeint: stepsize too small {} {} {} {} {}", w1, ystart, a, b, h1);
      return;
    }
  }
  io::log(" Error in odeint: too many steps {} {} {} {} {}", w1, ystart, a, b, h1);
}

// ***********************************************************************
// ***	  function rkstepper
// ***********************************************************************
void rkstepper(double& x, double& y, const double& dydx, const double& htest, double& hdid,
               double& hnew, const double& yscale, const double& eps, const int& w1) {
  // C--local variables
  constexpr double powerdown = 0.25;       // data powerdown/0.25/
  constexpr double powerup = r4(0.2f);     // data powerup/0.2/
  constexpr double safety = r4(0.9f);      // data safety/0.9/
  constexpr double maxdown = 10.0;         // data maxdown/10./
  constexpr double maxup = 5.0;            // data maxup/5./
  double yhalf, y1, y2, dydxhalf, xnew, delta, err, h, fac;
  bool reject;

  reject = false;
  h = htest;
l10:
  xnew = x + h;
  if (x == xnew) {
    io::log(" Error in rkstepper: step size not significant");
    return;
  }
  yhalf = rk4step(x, y, dydx, h / 2.0, w1);
  dydxhalf = deriv(x + h / 2.0, w1);
  y2 = rk4step(x + h / 2.0, yhalf, dydxhalf, h / 2.0, w1);
  y1 = rk4step(x, y, dydx, h, w1);
  delta = y2 - y1;
  err = std::fabs(delta) / (yscale * eps);
  if (err > 1.0) {
    reject = true;
    fac = std::max(1.0 / maxdown, safety / std::pow(err, powerdown));
    h = h * fac;
    goto l10;
  } else {
    if (reject) {
      hnew = h;
    } else {
      fac = std::min(maxup, safety / std::pow(err, powerup));
      hnew = fac * h;
    }
    x = xnew;
    y = y2 + delta / 15.0;
    hdid = h;
  }
}

// ***********************************************************************
// ***	  function rk4step
// ***********************************************************************
double rk4step(const double& x, const double& y, const double& dydx, const double& h,
               const int& w1) {
  const double k1 = h * dydx;
  const double k2 = h * deriv(x + h / 2.0, w1);
  const double k4 = h * deriv(x + h, w1);
  const double yout = y + k1 / 6.0 + 2.0 * k2 / 3.0 + k4 / 6.0;
  return yout;
}

}  // namespace jewel
