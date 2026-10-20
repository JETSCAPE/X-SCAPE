// JEWEL 2.6.0 C++ port -- simple medium model (Bjorken expanding, Woods-Saxon
// nuclear geometry). Translated from medium-simple.f:61-829.
//
// The medium's own COMMON blocks (MEDPARAM, MEDPARAMINT, RAPMAX2, BOOSTMED,
// MDFAC, THICKFNC, CROSSSEC, LOGFILE, INTEG) are not shared with the core and
// are therefore private data members of MediumSimple (Fortran names in lower
// case). NPART, NTHICK, CALCTA and CALCXSECTION are private member functions.
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include "jewel/jfmt.hpp"
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>

#include "jewel/farray.hpp"
#include "jewel/fortran_compat.hpp"
#include "jewel/io.hpp"
#include "jewel/medium.hpp"
#include "jewel/pythia6.hpp"
#include "jewel/routines.hpp"

namespace jewel {

namespace {

// ---------------------------------------------------------------------------
// Helpers emulating the gfortran list-directed READ(BUFFER,*,IOSTAT=IOS) X of
// a single item from an internal record. On any error / end of record / null
// value the variable is left unchanged (as gfortran does); the IOSTAT value is
// irrelevant in MEDINIT (it is overwritten by the next READ '(A)').

// First list-directed item of the record, with a repeat count "r*" removed.
// Returns nullopt for end of record, a null value or a slash.
std::optional<std::string> ld_first_item(std::string_view s) {
  auto is_blank = [](char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; };
  std::size_t i = 0;
  while (i < s.size() && is_blank(s[i])) ++i;
  if (i == s.size()) return std::nullopt;             // end of record
  if (s[i] == ',' || s[i] == '/') return std::nullopt;  // null value / slash
  std::size_t j = i;
  while (j < s.size() && !is_blank(s[j]) && s[j] != ',' && s[j] != '/') ++j;
  std::string tok(s.substr(i, j - i));
  const auto star = tok.find('*');
  if (star != std::string::npos) {
    // repeat count r*c: r must be a positive integer constant
    if (star == 0) return std::nullopt;
    for (std::size_t k = 0; k < star; ++k)
      if (tok[k] < '0' || tok[k] > '9') return std::nullopt;
    tok = tok.substr(star + 1);
    if (tok.empty()) return std::nullopt;  // r* = r null values
  }
  return tok;
}

bool is_digit(char c) { return c >= '0' && c <= '9'; }

// REAL / DOUBLE PRECISION item (gfortran validates, then converts with strtod).
void ld_read(std::string_view rec, double& x) {
  const auto item = ld_first_item(rec);
  if (!item) return;
  const std::string& t = *item;
  std::string n;
  std::size_t k = 0;
  if (k < t.size() && (t[k] == '+' || t[k] == '-')) n += t[k++];
  if (k < t.size() && !is_digit(t[k]) && t[k] != '.') {
    // Inf / Infinity / NaN
    std::string low;
    for (std::size_t m = k; m < t.size(); ++m)
      low += static_cast<char>((t[m] >= 'A' && t[m] <= 'Z') ? t[m] - 'A' + 'a' : t[m]);
    if (low == "inf" || low == "infinity" || low.starts_with("nan")) {
      x = std::strtod((n + (low.starts_with("nan") ? std::string("nan") : low)).c_str(), nullptr);
    }
    return;
  }
  int ndig = 0;
  while (k < t.size() && is_digit(t[k])) { n += t[k++]; ++ndig; }
  if (k < t.size() && t[k] == '.') {
    n += t[k++];
    while (k < t.size() && is_digit(t[k])) { n += t[k++]; ++ndig; }
  }
  if (ndig == 0) return;  // bad real
  if (k < t.size()) {
    const char c = t[k];
    if (c == 'e' || c == 'E' || c == 'd' || c == 'D' || c == 'q' || c == 'Q') {
      ++k;
    } else if (c != '+' && c != '-') {
      return;  // bad real
    }
    n += 'e';
    if (k < t.size() && (t[k] == '+' || t[k] == '-')) n += t[k++];
    int nexp = 0;
    while (k < t.size() && is_digit(t[k])) { n += t[k++]; ++nexp; }
    if (nexp == 0 || k != t.size()) return;  // bad real
  }
  x = std::strtod(n.c_str(), nullptr);
}

// INTEGER item
void ld_read(std::string_view rec, int& iv) {
  const auto item = ld_first_item(rec);
  if (!item) return;
  std::string_view t = *item;
  if (!t.empty() && t.front() == '+') t.remove_prefix(1);
  int v = 0;
  const auto [p, ec] = std::from_chars(t.data(), t.data() + t.size(), v);
  if (ec != std::errc() || p != t.data() + t.size() || t.empty()) return;  // bad integer
  iv = v;
}

// LOGICAL item: optional '.', then T or F (rest of the item is ignored)
void ld_read(std::string_view rec, bool& lv) {
  const auto item = ld_first_item(rec);
  if (!item) return;
  std::string_view t = *item;
  if (!t.empty() && t.front() == '.') t.remove_prefix(1);
  if (t.empty()) return;
  if (t.front() == 't' || t.front() == 'T') lv = true;
  else if (t.front() == 'f' || t.front() == 'F') lv = false;
}

// Fortran CHARACTER comparison ignores trailing blanks.
std::string_view rtrim(std::string_view s) {
  while (!s.empty() && s.back() == ' ') s.remove_suffix(1);
  return s;
}

const char* flog(bool b) { return b ? "T" : "F"; }

// ---------------------------------------------------------------------------
class MediumSimple final : public Medium {
 public:
  void medinit(const std::string& file, const int& id, const double& etam, const int& mass) override;
  void mednextevt() override;
  double getcentrality() override;
  void pickvtx(double& x, double& y) override;
  void setb(const double& bval) override;
  void getscatterer(const double& x, const double& y, const double& z, const double& t, int& type,
                    double& px, double& py, double& pz, double& e, double& ms) override;
  void avscatcen(const double& x, const double& y, const double& z, const double& t, double& px,
                 double& py, double& pz, double& e, double& m) override;
  void maxscatcen(double& px, double& py, double& pz, double& e, double& m) override;
  double getmd(const double& x1, const double& y1, const double& z1, const double& t1) override;
  double getms(const double& x2, const double& y2, const double& z2, const double& t2) override;
  double getneff(const double& x3, const double& y3, const double& z3, const double& t3) override;
  double gettemp(const double& x4, const double& y4, const double& z4, const double& t4) override;
  double gettempmax() override;
  double getmdmax() override;
  double getmdmin() override;
  double getmsmax() override;
  double getnatmdmin() override;
  double getltimemax() override;
  double getneffmax() override;
  double medderiv(const double& xval, const int& w) override;

 private:
  double npart(const double& xx1, const double& yy1, const double& xx2, const double& yy2);
  double nthick(const double& x1, const double& y1);
  void calcta();
  void calcxsection();

  // C--medium parameters  COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,NF
  double centrmin = 0, centrmax = 0, breal = 0, centr = 0;
  double rau = 10.0;  // DATA RAU/10./
  int nf = 0;
  // C--internal medium parameters
  // COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,N0,SIGMANN,A,WOODSSAXON
  double taui = 0, ti = 0, tc = 0;
  double d3 = 0.9;     // DATA D3/0.9d0/
  double zeta3 = 1.2;  // DATA ZETA3/1.2d0/
  double d = 0, n0 = 0, sigmann = 0;
  int a = 0;
  bool woodssaxon = false;
  // C--max rapidity  common/rapmax2/etamax2
  double etamax2 = 0;
  // C--longitudinal boost of momentum distribution  common/boostmed/boost
  bool boost = false;
  // C--factor to vary Debye mass  COMMON/MDFAC/MDFACTOR,MDSCALEFAC
  double mdfactor = 0, mdscalefac = 0;
  // Debye mass m_D = mdcoef * T (floored at MDFACTOR):
  //   MDMODE 0 (default, JEWEL 2.6.0): mdcoef = MDSCALEFAC * 3
  //   MDMODE 1: m_D^2 = (4 pi ALPHASMED / 3) (Nc + NF/2) T^2, NF = thermal light flavours of this file
  //             (same form as e.g. JETSCAPE, arXiv:2208.00983, below Eq. 6).  ALPHASMED defaults to the value
  //             that reproduces MDSCALEFAC * 3 for NF = 3, so the existing tune is kept.
  int mdmode = 0;
  // static brick for energy-loss tests (BRICK T): T = TI for 0 < t <= LBRICK [fm], everywhere in space,
  // production point at the origin, scattering centres at rest (no longitudinal flow)
  bool brick = false;
  double lbrick = 5.0;
  double alphasmed = -1.0, mdcoef = 0;
  // C--nuclear thickness function  COMMON /THICKFNC/ RMAX,TA(100,2)
  double rmax = 0;
  F2<double, 100, 2> ta;
  // C--geometrical cross section  COMMON /CROSSSEC/ IMPMAX,CROSS(200,3)
  double impmax = 0;
  F2<double, 200, 3> cross;
  // C--identifier of log file  common/logfile/logfid
  int logfid = 0;
  // C--variables for integration  COMMON/INTEG/B,R
  double b = 0, r = 0;
};

// ---------------------------------------------------------------------------
void MediumSimple::medinit(const std::string& file_arg, const int& id, const double& etam,
                           const int& mass) {
  // CHARACTER*80 FILE (gfortran trims trailing blanks of file names)
  const std::string file(rtrim(std::string_view(file_arg).substr(0, 80)));
  // C--local variables: LUN=77 / IOS are replaced by a std::ifstream and
  // stream state (IOS of the internal READs is never used).

  etamax2 = etam;
  logfid = id;

  // C--default settings
  taui = 0.4;
  ti = 0.59;
  tc = 0.17;
  woodssaxon = true;
  centrmin = 0.0;
  centrmax = 10.0;
  nf = 3;
  a = mass;
  n0 = 0.17;
  d = 0.54;
  sigmann = r4(6.2f);
  mdfactor = 0.45;
  mdscalefac = 1.0;
  boost = true;

  // C--read settings from file
  io::fortran_write(logfid, "");
  std::error_code ec;
  const bool fileexist = !file.empty() && std::filesystem::exists(file, ec);
  if (fileexist) {
    io::fortran_write(logfid, jewel::jfmt::format(" Reading medium parameters from {}", file));
    // OPEN(unit=LUN,file=FILE,status='old',err=10)
    std::ifstream in;
    if (!std::filesystem::is_directory(file, ec)) in.open(file);
    if (in.is_open()) {
      for (int i = 1; i <= 1000; ++i) {
        // READ(LUN, '(A)', iostat=ios) BUFFER   (CHARACTER*100)
        std::string line;
        if (!std::getline(in, line)) break;  // if (ios.ne.0) goto 30
        std::string buffer = line.substr(0, 100);
        buffer.resize(100, ' ');
        const char firstchar = buffer[0];
        if (firstchar == '#') continue;
        // POS=SCAN(BUFFER,' ')
        const auto blank = buffer.find(' ');
        const std::size_t pos = (blank == std::string::npos) ? 0 : blank + 1;
        const std::string label = buffer.substr(0, pos);
        buffer = buffer.substr(pos);
        const std::string_view lab = rtrim(label);
        if (lab == "TAUI") {
          ld_read(buffer, taui);
        } else if (lab == "TI") {
          ld_read(buffer, ti);
        } else if (lab == "TC") {
          ld_read(buffer, tc);
        } else if (lab == "WOODSSAXON") {
          ld_read(buffer, woodssaxon);
        } else if (lab == "CENTRMIN") {
          ld_read(buffer, centrmin);
        } else if (lab == "CENTRMAX") {
          ld_read(buffer, centrmax);
        } else if (lab == "NF") {
          ld_read(buffer, nf);
        } else if (lab == "N0") {
          ld_read(buffer, n0);
        } else if (lab == "D") {
          ld_read(buffer, d);
        } else if (lab == "SIGMANN") {
          ld_read(buffer, sigmann);
        } else if (lab == "MDFACTOR") {
          ld_read(buffer, mdfactor);
        } else if (lab == "MDSCALEFAC") {
          ld_read(buffer, mdscalefac);
        } else if (lab == "BRICK") {
          ld_read(buffer, brick);
        } else if (lab == "LBRICK") {
          ld_read(buffer, lbrick);
        } else if (lab == "MDMODE") {
          ld_read(buffer, mdmode);
        } else if (lab == "ALPHASMED") {
          ld_read(buffer, alphasmed);
        } else {
          io::fortran_write(logfid, jewel::jfmt::format(" unknown label {}", lab));
        }
      }
      // 30 close(LUN,status='keep')
      in.close();
      io::fortran_write(logfid, " ...done");
    } else {
      // 10
      io::fortran_write(logfid,
                        " Could not open medium parameter file, will run with default settings.");
    }
  } else {
    io::fortran_write(logfid, " No medium parameter file found, will run with default settings.");
  }

  // 40
  io::fortran_write(logfid, " using parameters:");
  io::fortran_write(logfid, jewel::jfmt::format(" TAUI       = {}", io::flist(taui)));
  io::fortran_write(logfid, jewel::jfmt::format(" TI         = {}", io::flist(ti)));
  io::fortran_write(logfid, jewel::jfmt::format(" TC         = {}", io::flist(tc)));
  io::fortran_write(logfid, jewel::jfmt::format(" WOODSSAXON = {}", flog(woodssaxon)));
  io::fortran_write(logfid, jewel::jfmt::format(" CENTRMIN   = {}", io::flist(centrmin)));
  io::fortran_write(logfid, jewel::jfmt::format(" CENTRMAX   = {}", io::flist(centrmax)));
  io::fortran_write(logfid, jewel::jfmt::format(" NF         = {}", nf));
  io::fortran_write(logfid, jewel::jfmt::format(" A          = {}", a));
  io::fortran_write(logfid, jewel::jfmt::format(" N0         = {}", io::flist(n0)));
  io::fortran_write(logfid, jewel::jfmt::format(" D          = {}", io::flist(d)));
  io::fortran_write(logfid, jewel::jfmt::format(" SIGMANN    = {}", io::flist(sigmann)));
  io::fortran_write(logfid, jewel::jfmt::format(" MDFACTOR   = {}", io::flist(mdfactor)));
  io::fortran_write(logfid, jewel::jfmt::format(" MDSCALEFAC = {}", io::flist(mdscalefac)));
  if (brick) {
    boost = false;
    io::fortran_write(logfid, jewel::jfmt::format(" BRICK      = T: T = TI = {} GeV for 0 < t <= LBRICK = {} fm", ti, lbrick));
  }
  if (mdmode == 1) {
    constexpr double pi = 3.141592653589793;
    if (alphasmed < 0.0) alphasmed = ipow<2>(3.0 * mdscalefac) * 3.0 / (4.0 * pi * (3.0 + 3.0 / 2.0));
    mdcoef = std::sqrt(4.0 * pi * alphasmed / 3.0 * (3.0 + nf / 2.0));
    io::fortran_write(logfid, jewel::jfmt::format(" MDMODE     = 1: m_D^2 = 4 pi alpha_s/3 (Nc + NF/2) T^2"));
    io::fortran_write(logfid, jewel::jfmt::format(" ALPHASMED  = {}  ->  m_D = {} T", alphasmed, mdcoef));
  } else {
    mdcoef = mdscalefac * 3.0;
  }
  io::fortran_write(logfid, "");
  io::fortran_write(logfid, "");
  io::fortran_write(logfid, "");

  // C--calculate T_A(x,y)
  calcta();
  // C--calculate geometrical cross section
  calcxsection();
}

// ---------------------------------------------------------------------------
void MediumSimple::mednextevt() {
  // C--local variables
  int i;
  double rr, b1, b2;

  // C--pick an impact parameter
  rr = (py::pyr(0) * (centrmax - centrmin) + centrmin) / 100.0;
  i = 0;
  for (int j = 1; j <= 200; ++j) {
    if ((rr - cross(j, 3) / cross(200, 3)) >= 0.0) {
      i = i + 1;
    } else {
      break;  // goto 132
    }
  }
  // 132
  b1 = static_cast<double>(i - 1) * 0.1;
  b2 = static_cast<double>(i) * 0.1;
  // NOTE(port): for i=200 (r>=1) the Fortran reads CROSS(201,3), i.e. past the
  // end of the CROSSSEC common block; same out-of-range access here (i=0 maps to
  // CROSS(200,2) in both). Not reachable for CENTRMAX<=100.
  breal = (b2 * (cross(i, 3) / cross(200, 3) - rr) + b1 * (rr - cross(i + 1, 3) / cross(200, 3))) /
          (cross(i, 3) / cross(200, 3) - cross(i + 1, 3) / cross(200, 3));
  centr = rr;
}

double MediumSimple::getcentrality() { return centr; }

// ---------------------------------------------------------------------------
void MediumSimple::pickvtx(double& x, double& y) {
  // C--local variables
  double x1, x2, y1, y2, z, xval, yval, zval;

  if (brick) {
    x = 0.0;
    y = 0.0;
    return;
  }

  x1 = breal / 2.0 - rau;
  x2 = rau - breal / 2.0;
  y1 = -std::sqrt(4.0 * (rau * rau) - breal * breal) / 2.0;
  y2 = std::sqrt(4.0 * (rau * rau) - breal * breal) / 2.0;
  for (;;) {
    // 131
    xval = py::pyr(0) * (x2 - x1) + x1;
    yval = py::pyr(0) * (y2 - y1) + y1;
    if ((nthick(xval - breal / 2.0, yval) == 0.0) || nthick(xval + breal / 2.0, yval) == 0.0)
      continue;
    zval = py::pyr(0) * nthick(-breal / 2.0, 0.0) * nthick(breal / 2.0, 0.0);
    z = nthick(xval - breal / 2.0, yval) * nthick(xval + breal / 2.0, yval);
    if (zval > z) continue;
    break;
  }
  x = xval;
  y = yval;
}

void MediumSimple::setb(const double& bval) { breal = bval; }

// ---------------------------------------------------------------------------
void MediumSimple::getscatterer(const double& x, const double& y, const double& z, const double& t,
                                int& type, double& px, double& py, double& pz, double& e,
                                double& ms) {
  // C--local variables
  double md, temp;
  double rr, pmax, wt, theta, phi, p, ys, pz2, e2;
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/

  rr = py::pyr(0);
  if (rr < (static_cast<double>(24.f * static_cast<float>(nf)) * d3 / 3.0) /
               (static_cast<double>(24.f * static_cast<float>(nf)) * d3 / 3.0 +
                48.0 * zeta3 / 2.0)) {
    type = 2;
    ms = getms(x, y, z, t);
  } else {
    type = 21;
    ms = getmd(x, y, z, t);
  }
  temp = gettemp(x, y, z, t);
  [[maybe_unused]] double tau = std::sqrt(t * t - z * z);
  if (boost) {
    ys = 0.5 * std::log((t + z) / (t - z));
  } else {
    ys = 0.0;
  }
  pmax = 10.0 * temp;

  if (temp < 1.e-2) {
    io::fortran_write(logfid, " asking for a scattering centre without medium:");
    io::fortran_write(logfid, jewel::jfmt::format(" at (x,y,z,t)= {} {} {} {}", io::flist(x), io::flist(y),
                                          io::flist(z), io::flist(t)));
    io::fortran_write(logfid, " making one up to continue but something is wrong!");
    type = 21;
    px = 0.0;
    py = 0.0;
    pz = 0.0;
    const double zero = 0.0;
    ms = getms(zero, zero, zero, zero);
    md = getmd(zero, zero, zero, zero);
    static_cast<void>(md);
    e = std::sqrt(px * px + py * py + pz * pz + ms * ms);
    return;
  }

  do {
    // 10
    p = std::pow(py::pyr(0), r4(0.3333333f)) * pmax;
    e2 = std::sqrt(p * p + ms * ms);
    if (type == 2) {
      wt = (std::exp(ms / temp) - 1.0) / (std::exp(e2 / temp) - 1.0);
    } else {
      wt = (std::exp(ms / temp) + 1.0) / (std::exp(e2 / temp) + 1.0);
    }
    if (wt > 1.0) io::fortran_write(logfid, jewel::jfmt::format(" Error in getscatterer: weight = {}", io::flist(wt)));
    if (wt < 0.0) io::fortran_write(logfid, jewel::jfmt::format(" Error in getscatterer: weight = {}", io::flist(wt)));
  } while (py::pyr(0) > wt);  // goto 10
  phi = py::pyr(0) * 2.0 * pi;
  theta = -std::acos(2.0 * py::pyr(0) - 1.0) + pi;
  px = p * gf_sin(theta) * gf_cos(phi);
  py = p * gf_sin(theta) * gf_sin(phi);
  pz2 = p * gf_cos(theta);
  e = std::cosh(ys) * e2 + std::sinh(ys) * pz2;
  pz = std::sinh(ys) * e2 + std::cosh(ys) * pz2;
}

// ---------------------------------------------------------------------------
void MediumSimple::avscatcen(const double& x, const double& y, const double& z, const double& t,
                             double& px, double& py, double& pz, double& e, double& m) {
  // C--local variables
  double ys;

  if (std::fabs(z) >= t) {
    px = 0.0;
    py = 0.0;
    pz = 0.0;
    e = 0.0;
    m = 0.0;
    return;
  }
  if (boost) {
    ys = 0.5 * std::log((t + z) / (t - z));
    if ((z == 0.0) && (t == 0.0)) ys = 0.0;
    if (ys > etamax2) ys = etamax2;
    if (ys < -etamax2) ys = -etamax2;
  } else {
    ys = 0.0;
  }
  m = getms(x, y, z, t);
  e = m * std::cosh(ys);
  px = 0.0;
  py = 0.0;
  pz = m * std::sinh(ys);
}

// ---------------------------------------------------------------------------
void MediumSimple::maxscatcen(double& px, double& py, double& pz, double& e, double& m) {
  // C--local variables
  double ys;
  if (boost) {
    ys = etamax2;
  } else {
    ys = 0.0;
  }
  m = getmsmax();
  e = m * std::cosh(ys);
  px = 0.0;
  py = 0.0;
  pz = m * std::sinh(ys);
}

// ---------------------------------------------------------------------------
double MediumSimple::getmd(const double& x1, const double& y1, const double& z1, const double& t1) {
  double getmd_v = mdcoef * gettemp(x1, y1, z1, t1);
  getmd_v = std::max(getmd_v, mdfactor);
  return getmd_v;
}

double MediumSimple::getms(const double& x2, const double& y2, const double& z2, const double& t2) {
  // SQRT(2.) is a REAL*4 square root
  return getmd(x2, y2, z2, t2) / static_cast<double>(std::sqrt(2.0f));
}

double MediumSimple::getneff(const double& x3, const double& y3, const double& z3,
                             const double& t3) {
  // C--   local variables
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/
  if (brick)
    return (static_cast<double>(12.f * static_cast<float>(nf)) * d3 * 2.0 / 3.0 + 16.0 * zeta3 * 3.0 / 2.0) *
           ipow<3>(gettemp(x3, y3, z3, t3)) / ipow<2>(pi);
  const double tau = std::sqrt(t3 * t3 - z3 * z3);
  const double cosheta = t3 / tau;
  double getneff_v = (static_cast<double>(12.f * static_cast<float>(nf)) * d3 * 2.0 / 3.0 +
                      16.0 * zeta3 * 3.0 / 2.0) *
                     ipow<3>(gettemp(x3, y3, z3, t3)) / ipow<2>(pi);
  getneff_v = getneff_v / cosheta;
  return getneff_v;
}

// ---------------------------------------------------------------------------
double MediumSimple::gettemp(const double& x4, const double& y4, const double& z4,
                             const double& t4) {
  // C--local variables
  double tau, eps0, epsin, tempin, ys;
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/
  // REAL*4 sub-expression (16.*8.+7.*2.*6.*NF)
  const double dof = static_cast<double>(16.f * 8.f + 7.f * 2.f * 6.f * static_cast<float>(nf));

  double gettemp_v = 0.0;

  if (brick) return ((t4 > 0.0) && (t4 <= lbrick)) ? ti : 0.0;
  if (std::fabs(z4) > t4) return gettemp_v;

  tau = std::sqrt(t4 * t4 - z4 * z4);
  // C--check for overlap region
  if ((nthick(x4 - breal / 2.0, y4) == 0.0) || nthick(x4 + breal / 2.0, y4) == 0.0)
    return gettemp_v;

  ys = 0.5 * std::log((t4 + z4) / (t4 - z4));
  if (std::fabs(ys) > etamax2) return gettemp_v;
  // C--determine initial temperature at transverse position
  if (woodssaxon) {
    eps0 = dof * ipow<2>(pi) * ipow<4>(ti) / 240.0;
    //         EPSIN=EPS0*NPART(X4-BREAL/2.,Y4,X4+BREAL/2.,Y4)
    //     &        *PI*RAU**2/(2.*A)
    const double zero = 0.0;
    epsin = eps0 * npart(x4 - breal / 2.0, y4, x4 + breal / 2.0, y4) / npart(zero, zero, zero, zero);
    tempin = std::pow(epsin * 240.0 / (ipow<2>(pi) * dof), 0.25);
  } else {
    tempin = ti;
  }
  // C--calculate temperature if before initial time
  if (tau <= taui) {
    gettemp_v = tempin * tau / taui;
  } else {
    // C--evolve temperature
    gettemp_v = tempin * std::pow(taui / tau, r4(0.3333f));
  }
  if (gettemp_v < tc) gettemp_v = 0.0;
  return gettemp_v;
}

// ---------------------------------------------------------------------------
double MediumSimple::gettempmax() {
  if (brick) return ti;
  const double zero = 0.0;
  return gettemp(zero, zero, zero, taui);
}

double MediumSimple::getmdmax() {
  double getmdmax_v = mdcoef * gettempmax();
  getmdmax_v = std::max(getmdmax_v, mdfactor);
  return getmdmax_v;
}

double MediumSimple::getmdmin() {
  double getmdmin_v = mdcoef * tc;
  getmdmin_v = std::max(getmdmin_v, mdfactor);
  return getmdmin_v;
}

double MediumSimple::getmsmax() { return getmdmax() / std::sqrt(2.0); }

double MediumSimple::getnatmdmin() {
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/
  // C--local variables
  const double tt = getmdmin() / mdcoef;
  return (static_cast<double>(12.f * static_cast<float>(nf)) * d3 * 2.0 / 3.0 +
          16.0 * zeta3 * 3.0 / 2.0) *
         ipow<3>(tt) / ipow<2>(pi);
}

double MediumSimple::getltimemax() {
  if (brick) return lbrick;
  return taui * ipow<3>(gettempmax() / tc) * std::cosh(etamax2);
}

double MediumSimple::getneffmax() {
  // C--   local variables
  constexpr double pi = 3.141592653589793;  // DATA PI/3.141592653589793d0/
  return (static_cast<double>(12.f * static_cast<float>(nf)) * d3 * 2.0 / 3.0 +
          16.0 * zeta3 * 3.0 / 2.0) *
         ipow<3>(gettempmax()) / ipow<2>(pi);
}

// ---------------------------------------------------------------------------
double MediumSimple::npart(const double& xx1, const double& yy1, const double& xx2,
                           const double& yy2) {
  return nthick(xx1, yy1) * (1.0 - std::exp(-sigmann * nthick(xx2, yy2))) +
         nthick(xx2, yy2) * (1.0 - std::exp(-sigmann * nthick(xx1, yy1)));
}

// ---------------------------------------------------------------------------
double MediumSimple::nthick(const double& x1, const double& y1) {
  // Locals R, B, C (R and B hide the INTEG members, as in the Fortran where
  // NTHICK does not declare COMMON/INTEG/).
  int line, lmin;
  F1<double, 4> xa, ya;
  double r, c, b, nthick_v;

  r = std::sqrt(x1 * x1 + y1 * y1);
  if (r > ta(100, 1)) {
    nthick_v = 0.0;
  } else {
    line = fint(r * 99.0 / ta(100, 1) + 1.0);
    lmin = std::max(line, 1);
    lmin = std::min(lmin, 99);
    if ((r < ta(lmin, 1)) || (r > ta(lmin + 1, 1)))
      io::fortran_write(logfid, jewel::jfmt::format(" {} {} {} {} {}", line, lmin, io::flist(r),
                                            io::flist(ta(lmin, 1)), io::flist(ta(lmin + 1, 1))));
    xa(1) = ta(lmin, 1);
    xa(2) = ta(lmin + 1, 1);
    ya(1) = ta(lmin, 2);
    ya(2) = ta(lmin + 1, 2);
    c = (ya(2) - ya(1)) / (xa(2) - xa(1));
    b = ya(1) - c * xa(1);
    nthick_v = c * r + b;
  }
  return nthick_v;
}

// ---------------------------------------------------------------------------
void MediumSimple::calcta() {
  // C--local variables
  int nsteps;
  double eps, hfirst, y;

  nsteps = 100;
  eps = r4(1.e-4f);
  hfirst = 0.1;

  // INTEGER**REAL*4 is evaluated in single precision (A converted to REAL*4)
  {
    const float af = static_cast<float>(a);
    r = static_cast<double>(1.12f * std::pow(af, 0.33333f) - 0.86f * std::pow(af, -0.33333f));
  }
  rmax = 2.0 * r;

  for (int i = 1; i <= nsteps; ++i) {
    // C--set transverse position
    b = static_cast<double>(i - 1) * 2.0 * r / static_cast<double>(nsteps);
    y = 0.0;
    // C--integrate along longitudinal line
    const double lo = -(2.0 * r), hi = 2.0 * r, hmin = 0.0;
    const int w = 101;
    odeint(y, lo, hi, eps, hfirst, hmin, w);
    ta(i, 1) = b;
    ta(i, 2) = y;
  }
}

// ---------------------------------------------------------------------------
void MediumSimple::calcxsection() {
  // C--local variables (B, X, Y are locals; B hides the INTEG member as in Fortran)
  double b, p, prod, x, y, pprev;

  pprev = 0.0;
  for (int ib = 1; ib <= 200; ++ib) {
    b = 0.1 * static_cast<double>(ib);
    prod = 1.0;
    for (int ix = 1; ix <= 100; ++ix) {
      for (int iy = 1; iy <= 100; ++iy) {
        x = -20.0 + static_cast<double>(ix) * 0.4;
        y = -20.0 + static_cast<double>(iy) * 0.4;
        prod = prod * std::pow(std::exp(-nthick(x + b / 2.0, y) * sigmann),
                               0.16 * nthick(x - b / 2.0, y));
      }
    }
    p = (1.0 - prod) * 8.8 / 14.0 * b;
    cross(ib, 1) = b;
    cross(ib, 2) = p;
    if (ib == 1) {
      cross(ib, 3) = 0.0;
    } else {
      cross(ib, 3) = cross(ib - 1, 3) + (p + pprev) / 2.0 * r4(0.1f);
    }
    pprev = p;
  }
  impmax = r4(19.95f);
}

// ---------------------------------------------------------------------------
double MediumSimple::medderiv(const double& xval, const int& w) {
  if (w == 1) {
    // C--XVAL corresponds to z-coordinate
    return n0 / (1.0 + std::exp((std::sqrt(b * b + xval * xval) - r) / d));
  } else {
    return 0.0;
  }
}

}  // namespace

std::unique_ptr<Medium> make_medium_simple() { return std::make_unique<MediumSimple>(); }

}  // namespace jewel
