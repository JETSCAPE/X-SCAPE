// JEWEL 2.6.0 C++ port -- interface for use as a final-state shower in JETSCAPE
// (see jetscape_api.hpp).
#include "jewel/jetscape_api.hpp"

#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <utility>

#include "jewel/commons.hpp"
#include "jewel/fortran_compat.hpp"
#include "jewel/heavyflavour.hpp"
#include "jewel/io.hpp"
#include "jewel/medium.hpp"
#include "jewel/pythia6.hpp"
#include "jewel/routines.hpp"

namespace jewel::api {

namespace {

constexpr double kPi = 3.141592653589793;
constexpr double kD3 = 0.9, kZeta3 = 1.2;  // as in medium-simple.f

Stats g_stats;

// ---------------------------------------------------------------------------
// Medium from an external fluid (temperature and flow velocity per cell).
// Same Debye mass, thermal masses and density as medium-simple.f; the Bjorken
// boost of the scattering centres is replaced by a boost with the local fluid
// velocity, the 1/cosh(eta) flux factor of GETNEFF by gamma (1 - v.beta_probe).
class MediumJetscape final : public Medium {
 public:
  MediumJetscape(MediumSettings s, CellQuery q) : s_(s), query_(std::move(q)) {
    mdcoef_ = s_.mdscalefac * 3.0;
  }

  void set_tau0(double tau0) {
    s_.tau0 = tau0;
    cached_ = false;
  }
  void set_vertex(double x, double y) {
    vx_ = x;
    vy_ = y;
  }

  void medinit(const std::string&, const int& id, const double&, const int&) override {
    const int logfid = id;
    io::fortran_write(logfid, "");
    io::fortran_write(logfid, " JETSCAPE medium (temperature and flow from the JETSCAPE fluid):");
    io::fortran_write(logfid, jfmt::format(" TMAX       = {}", s_.tmax));
    io::fortran_write(logfid, jfmt::format(" TC         = {}", s_.tc));
    io::fortran_write(logfid, jfmt::format(" NF         = {}", s_.nf));
    io::fortran_write(logfid, jfmt::format(" MDFACTOR   = {}", s_.mdfactor));
    io::fortran_write(logfid, jfmt::format(" MDSCALEFAC = {}", s_.mdscalefac));
    io::fortran_write(logfid, jfmt::format(" LIFETIME   = {} fm", s_.lifetime));
    io::fortran_write(logfid, jfmt::format(" TAU0       = {} fm (PRE_TAU0 = {})", s_.tau0, s_.pre_tau0));
    io::fortran_write(logfid, jfmt::format(" FLUXMAX    = {}", s_.fluxmax));
    io::fortran_write(logfid, jfmt::format(" FLUX_MODE  = {}", s_.flux_mode));
    io::fortran_write(logfid, jfmt::format(" ETAMAX     = {}", s_.etamax));
    io::fortran_write(logfid, jfmt::format(" VACUUM     = {}", s_.vacuum ? "T" : "F"));
    io::fortran_write(logfid, "");
  }
  void mednextevt() override {}
  double getcentrality() override { return -1.0; }
  void pickvtx(double& x, double& y) override {
    x = vx_;
    y = vy_;
  }
  void setb(const double&) override {}

  void getscatterer(const double& x, const double& y, const double& z, const double& t, int& type,
                    double& px, double& py, double& pz, double& e, double& ms) override {
    // flavour of the scattering centre (quark or gluon), as in medium-simple.f
    const double rr = py::pyr(0);
    const double wq = 24.0 * s_.nf * kD3 / 3.0;
    if (rr < wq / (wq + 48.0 * kZeta3 / 2.0)) {
      type = 2;
      ms = getms(x, y, z, t);
    } else {
      type = 21;
      ms = getmd(x, y, z, t);
    }
    const MediumCell& c = cell(x, y, z, t);
    const double temp = c.temperature;
    if (temp < 1.e-2) {
      io::log(" asking for a scattering centre without medium:");
      io::log(" at (x,y,z,t)= {} {} {} {}", x, y, z, t);
      io::log(" making one up to continue but something is wrong!");
      type = 21;
      px = py = pz = 0.0;
      ms = s_.mdfactor / std::sqrt(2.0);
      e = ms;
      return;
    }
    // thermal momentum in the fluid rest frame (Bose / Fermi weight)
    const double pmax = 10.0 * temp;
    double p, e2, wt;
    do {
      p = std::cbrt(py::pyr(0)) * pmax;
      e2 = std::sqrt(p * p + ms * ms);
      if (type == 2) {
        wt = (std::exp(ms / temp) - 1.0) / (std::exp(e2 / temp) - 1.0);
      } else {
        wt = (std::exp(ms / temp) + 1.0) / (std::exp(e2 / temp) + 1.0);
      }
    } while (py::pyr(0) > wt);
    const double phi = py::pyr(0) * 2.0 * kPi;
    const double theta = -std::acos(2.0 * py::pyr(0) - 1.0) + kPi;
    double q[4] = {e2, p * std::sin(theta) * std::cos(phi), p * std::sin(theta) * std::sin(phi),
                   p * std::cos(theta)};
    boost(c, q);  // to the lab frame, as LBT does for its thermal partons
    e = q[0];
    px = q[1];
    py = q[2];
    pz = q[3];
  }

  // average scattering centre: at rest in the fluid frame, mass m_s
  void avscatcen(const double& x, const double& y, const double& z, const double& t, double& px,
                 double& py, double& pz, double& e, double& m) override {
    if (std::fabs(z) >= t) {
      px = py = pz = e = m = 0.0;
      return;
    }
    const MediumCell& c = cell(x, y, z, t);
    m = getms(x, y, z, t);
    const double g = gamma(c);
    e = m * g;
    px = m * g * c.vx;
    py = m * g * c.vy;
    pz = m * g * c.vz;
  }
  // bound used for the cross-section tables (as medium-simple.f: boosted with ETAMAX)
  void maxscatcen(double& px, double& py, double& pz, double& e, double& m) override {
    m = getmsmax();
    e = m * std::cosh(s_.etamax);
    px = 0.0;
    py = 0.0;
    pz = m * std::sinh(s_.etamax);
  }

  double getmd(const double& x, const double& y, const double& z, const double& t) override {
    return std::max(mdcoef_ * gettemp(x, y, z, t), s_.mdfactor);
  }
  double getms(const double& x, const double& y, const double& z, const double& t) override {
    return getmd(x, y, z, t) / static_cast<double>(std::sqrt(2.0f));
  }
  // density in the fluid rest frame
  double getneff(const double& x, const double& y, const double& z, const double& t) override {
    return density(gettemp(x, y, z, t));
  }
  // scattering rate per unit lab time / cross section: n u.p / E = n gamma (1 - v.beta)
  double getneffprobe(const double& x, const double& y, const double& z, const double& t,
                      const double& bx, const double& by, const double& bz) override {
    const MediumCell& c = cell(x, y, z, t);
    if (c.temperature <= 0.0) return 0.0;
    double flux;
    if (s_.flux_mode == 1 && t > std::fabs(z)) {
      // medium-2D.f: transverse velocity at eta_s = 0 is v_lab,T cosh eta_s for boost-invariant flow
      const double ch = t / std::sqrt(t * t - z * z);
      const double ux = c.vx * ch, uy = c.vy * ch, u = std::sqrt(ux * ux + uy * uy);
      const double pt = std::sqrt(bx * bx + by * by);
      const double yt = u < 1.0 ? std::atanh(u) : 0.0;
      const double cth = (u > 0.0 && pt > 0.0) ? (ux * bx + uy * by) / (u * pt) : 0.0;
      flux = (std::cosh(yt) - std::sinh(yt) * cth) / ch;
    } else {
      flux = gamma(c) * (1.0 - (c.vx * bx + c.vy * by + c.vz * bz));
    }
    if (flux > g_stats.max_flux) g_stats.max_flux = flux;
    if (flux > s_.fluxmax) ++g_stats.n_above_fluxmax;
    return density(c.temperature) * flux;
  }
  double gettemp(const double& x, const double& y, const double& z, const double& t) override {
    return cell(x, y, z, t).temperature;
  }
  double gettempmax() override { return s_.tmax; }
  double getmdmax() override { return std::max(mdcoef_ * s_.tmax, s_.mdfactor); }
  double getmdmin() override { return std::max(mdcoef_ * s_.tc, s_.mdfactor); }
  double getmsmax() override { return getmdmax() / std::sqrt(2.0); }
  // bounds of the rate in GETDELTAT: like getneffprobe they include the flow
  // factor gamma (1 - v.beta) <= fluxmax (without it, cells just above T_c that
  // move towards the parton exceed the bound: "weight > 1")
  double getnatmdmin() override { return density(getmdmin() / mdcoef_) * s_.fluxmax; }
  double getltimemax() override { return s_.lifetime; }
  double getneffmax() override { return density(s_.tmax) * s_.fluxmax; }
  double medderiv(const double&, const int&) override { return 0.0; }

 private:
  double density(double temp) const {
    return (12.0 * s_.nf * kD3 * 2.0 / 3.0 + 16.0 * kZeta3 * 3.0 / 2.0) * temp * temp * temp /
           (kPi * kPi);
  }
  static double gamma(const MediumCell& c) {
    const double v2 = c.vx * c.vx + c.vy * c.vy + c.vz * c.vz;
    return v2 < 1.0 ? 1.0 / std::sqrt(1.0 - v2) : 1.0;
  }
  // q = (E, px, py, pz) from the fluid rest frame to the lab frame
  static void boost(const MediumCell& c, double q[4]) {
    const double v2 = c.vx * c.vx + c.vy * c.vy + c.vz * c.vz;
    if (v2 <= 0.0 || v2 >= 1.0) return;
    const double g = 1.0 / std::sqrt(1.0 - v2);
    const double vq = c.vx * q[1] + c.vy * q[2] + c.vz * q[3];
    const double f = (g - 1.0) * vq / v2 + g * q[0];
    q[0] = g * (q[0] + vq);
    q[1] += f * c.vx;
    q[2] += f * c.vy;
    q[3] += f * c.vz;
  }

  // cell at (t, x, y, z); the JEWEL core asks for the same point several times in a row
  const MediumCell& cell(double x, double y, double z, double t) {
    if (cached_ && t == ct_ && x == cx_ && y == cy_ && z == cz_) return cell_;
    ct_ = t;
    cx_ = x;
    cy_ = y;
    cz_ = z;
    cached_ = true;
    cell_ = MediumCell{};
    if (s_.vacuum || t <= std::fabs(z) || t > s_.lifetime) return cell_;
    const double etas = 0.5 * std::log((t + z) / (t - z));
    if (std::fabs(etas) > s_.etamax) return cell_;
    const double tau = std::sqrt(t * t - z * z);
    if (tau < s_.tau0) {
      if (s_.pre_tau0 != 1) return cell_;
      // linear rise before tau0 at the same (x, y, eta_s)
      cell_ = query_(s_.tau0 * std::cosh(etas), x, y, s_.tau0 * std::sinh(etas));
      cell_.temperature *= tau / s_.tau0;
    } else {
      cell_ = query_(t, x, y, z);
    }
    if (cell_.temperature > s_.tmax) ++g_stats.n_above_tmax;
    if (cell_.temperature < s_.tc) cell_ = MediumCell{};
    return cell_;
  }

  MediumSettings s_;
  CellQuery query_;
  double mdcoef_ = 3.0;
  double vx_ = 0.0, vy_ = 0.0;
  bool cached_ = false;
  double ct_ = 0, cx_ = 0, cy_ = 0, cz_ = 0;
  MediumCell cell_;
};

MediumJetscape* g_medium = nullptr;

void write_param_lines(const Settings& s, const std::string& file) {
  std::ofstream out(file);
  if (!out) throw std::runtime_error("JEWEL: cannot write " + file);
  if (!s.param_file.empty()) {
    std::ifstream in(s.param_file);
    if (!in) throw std::runtime_error("JEWEL: cannot read parameter file " + s.param_file);
    std::string line;
    while (std::getline(in, line)) out << line << '\n';
  }
  out << "# set by the JETSCAPE module (override the lines above)\n";
  for (const auto& l : s.extra_params) out << l << '\n';
}

}  // namespace

void initialise(const Settings& settings, CellQuery cell_query) {
  auto m = std::make_unique<MediumJetscape>(settings.medium, std::move(cell_query));
  g_medium = m.get();
  set_medium(std::move(m));
  // JEWEL reads a private copy of the parameter file, which is then kept as
  // run_param_file: jobs started at the same time in one directory never read
  // each other's (or a half-written) file
  const std::string mine = settings.run_param_file + ".tmp" + std::to_string(::getpid());
  write_param_lines(settings, mine);
  init_library(mine);
  if (std::rename(mine.c_str(), settings.run_param_file.c_str()) != 0) std::remove(mine.c_str());

  using namespace py::pydat1;
  using namespace py::pydatr;
  using namespace cb::param;
  using namespace cb::storescatcen;
  // PYTHIA 6 settings of INITPYTHIA that the shower uses (no hard process is generated here)
  PARJ(81) = lps;
  PARJ(82) = q0;
  MSTU(4) = 23000;
  MSTU(5) = 10000;
  // scattering centres are needed for the holes; no event-wide subtraction
  dosubtraction = false;
  writescatcen = allhad;
  // random numbers: PYR restarts from MRPY(1) when MRPY(2) = 0
  MRPY(1) = static_cast<int>(settings.seed % 900000000L);
  MRPY(2) = 0;
  io::log(" JETSCAPE: PYR seed {}", MRPY(1));
  io::log_flush();
}

void next_event(double tau0) {
  if (g_medium && tau0 > 0.0) g_medium->set_tau0(tau0);
}

bool shower(int id, const double p[4], const double x[4], double qmax, std::vector<OutParton>& out) {
  using namespace py::pyjets;
  using namespace cb::param;
  using namespace cb::disc;
  using namespace cb::angor;
  using namespace cb::evoltimes;
  using namespace cb::ftimefac;
  using namespace cb::time;
  using namespace cb::colour;
  using namespace cb::jetpoint;
  using namespace cb::evrecord;
  using namespace cb::storescatcen;
  using namespace cb::speclines;

  out.clear();
  ++g_stats.n_showers;

  // C--reset the event record (as GENEVENT)
  N = 0;
  colmax = 600;
  discard = false;
  for (int i = 1; i <= cb::kLines; ++i) {
    for (int j = 1; j <= 5; ++j) mv(i, j) = 0.0;
    zd(i) = 0.0;
    thetad(i) = 0.0;
    qqbard(i) = false;
    cb::hf::qqbarfl(i) = 0;
    starttime(i) = -1.e6;
    medind(i) = false;
    cb::lib::origin(i) = -1;
  }
  nscatcen = 0;
  ltime = medium().getltimemax();
  x0 = x[1];
  y0 = x[2];
  if (g_medium) g_medium->set_vertex(x0, y0);

  // line 1 is JEWEL's work space for the scattering centre; line 2 documents the
  // incoming parton (status 21, like the hard-process lines of a PYTHIA event), so
  // that the mother-chain walks in the shower stop there; the parton is line 3
  N = 3;
  offset = 2;
  npart = 1;
  for (int l = 1; l <= 3; ++l)
    for (int j = 1; j <= 5; ++j) {
      K(l, j) = 0;
      P(l, j) = 0.0;
      V(l, j) = 0.0;
    }
  const int ldoc = 2, l = 3;
  K(ldoc, 1) = 21;
  K(ldoc, 2) = id;
  P(ldoc, 1) = p[1];
  P(ldoc, 2) = p[2];
  P(ldoc, 3) = p[3];
  P(ldoc, 4) = p[0];
  K(l, 1) = 1;
  K(l, 2) = id;
  K(l, 3) = ldoc;
  P(l, 1) = p[1];
  P(l, 2) = p[2];
  P(l, 3) = p[3];
  P(l, 4) = p[0];
  lme1 = lme2 = l;
  lme1orig = lme2orig = ldoc;
  type1 = type2 = (id == 21) ? "GC" : "QQ";
  trip(1) = anti(1) = trip(ldoc) = anti(ldoc) = 0;
  if (id == 21) {
    trip(l) = ++colmax;
    anti(l) = ++colmax;
  } else if (id > 0) {
    trip(l) = ++colmax;
    anti(l) = 0;
  } else {
    trip(l) = 0;
    anti(l) = ++colmax;
  }

  // C--virtuality of the parton (FINDVIRTUALITIES for a single parton: the energy is kept,
  // the three-momentum is reduced; no colour partner, so no angular constraint)
  const double e = p[0];
  const double mq = hfmass(id);
  double z = 0.0;
  bool which = false;
  int fl = 0;
  const std::string typ = (id == 21) ? "GC" : "QQ";
  double q = getmass(0.0, qmax, -1.0, e, typ, e, false, z, which, mq, &fl);
  for (int tries = 0;; ++tries) {
    const double m = invmass(q, mq);
    if (tries > 1000) {
      q = 0.0;
      z = 0.0;
      which = false;
      break;
    }
    if (m >= e) {
      q = getmass(0.0, q, -1.0, e, typ, e, false, z, which, mq, &fl);
      continue;
    }
    if (q > 0.0) {
      const double eps = 0.5 - 0.5 * std::sqrt(1.0 - ipow<2>(q0) / ipow<2>(q)) *
                                   std::sqrt(1.0 - ipow<2>(q) / ipow<2>(e));
      if (hfzout(q, mq, e, z, eps, which ? fl : 0)) {
        q = getmass(0.0, q, -1.0, e, typ, e, false, z, which, mq, &fl);
        continue;
      }
    }
    break;
  }
  {
    const double m = invmass(q, mq);
    const double pold = std::sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]);
    const double pnew = std::sqrt(std::max(e * e - m * m, 0.0));
    if (pold > 0.0) {
      P(l, 1) *= pnew / pold;
      P(l, 2) *= pnew / pold;
      P(l, 3) *= pnew / pold;
    }
    P(l, 5) = m;
  }
  if (z > 0.0) {
    thetad(l) = (mq > 0.0) ? hftheta(q * q, mq, z, e) : P(l, 5) / (std::sqrt(z * (1.0 - z)) * e);
  } else {
    thetad(l) = 0.0;
  }
  zd(l) = z;
  qqbard(l) = which;
  if (cb::hf::on) cb::hf::qqbarfl(l) = fl;
  aoreq(l) = false;
  // production point and formation time
  mv(l, 1) = x[1];
  mv(l, 2) = x[2];
  mv(l, 3) = x[3];
  mv(l, 4) = x[0];
  for (int j = 1; j <= 4; ++j) mv(ldoc, j) = mv(l, j);
  mv(ldoc, 5) = x[0];
  if (offshell(l)) {
    const double lambda = 1.0 / (ftfac * P(l, 4) * r4(0.2f) / ipow<2>(q));
    mv(l, 5) = x[0] - std::log(1.0 - py::pyr(0)) / lambda;
  } else {
    mv(l, 5) = ltime;
  }

  // C--develop parton shower
  makecascade();
  if (discard) {
    ++g_stats.n_discarded;
    io::log(" JETSCAPE: JEWEL discarded the shower of a parton with id {} and E = {}", id, e);
    return false;
  }
  // heavy flavour: Langevin transport of final-state heavy quarks
  if (cb::hf::on && (cb::hf::hqds > 0.0)) hfdiffuse();
  if (!allhad) {
    for (int i = 1; i <= N; ++i)
      if (K(i, 1) == 3) K(i, 1) = 22;
  }

  // final partons: shower (descends from line 2) or recoil (descends from a scattering
  // centre); lines descending from a scattering that REJECTSPLITTING undid (status -1)
  // are dropped
  for (int i = l; i <= N; ++i) {
    const int st = K(i, 1);
    if (st < 1 || st > 5) continue;
    const int origin = lib_origin(i);
    if (origin == 2) continue;
    OutParton o;
    o.id = K(i, 2);
    o.status = origin;
    o.p[0] = P(i, 4);
    o.p[1] = P(i, 1);
    o.p[2] = P(i, 2);
    o.p[3] = P(i, 3);
    o.mass = P(i, 5);
    o.x[0] = mv(i, 4);
    o.x[1] = mv(i, 1);
    o.x[2] = mv(i, 2);
    o.x[3] = mv(i, 3);
    o.col = trip(i);
    o.acol = anti(i);
    out.push_back(o);
  }
  // thermal momenta of the scattering centres before the scattering
  if (allhad) {
    for (int j = 1; j <= nscatcen; ++j) {
      if (scatcen(j, 4) <= 0.0) continue;
      // no hole without its recoil: REJECTSPLITTING can undo a scattering (its recoil
      // line gets status -1) but leaves the scattering centre stored
      const int r = cb::lib::scatrec(j);
      if (r == 0 || (r > 0 && K(r, 1) == -1)) {
        ++g_stats.n_holes_dropped;
        continue;
      }
      // JEWEL_LIBDEBUG: recoils that do not reach the output although their hole does
      if (std::getenv("JEWEL_LIBDEBUG") && r > 0 &&
          (K(r, 1) < 1 || K(r, 1) > 5 || lib_origin(r) == 2))
        io::log(" LIBDEBUG hole {} recoil line {} K={} {} {} {} {} origin {} E_hole {} E_rec {} N {}", j, r,
                K(r, 1), K(r, 2), K(r, 3), K(r, 4), K(r, 5), lib_origin(r), scatcen(j, 4), P(r, 4), N);
      OutParton o;
      o.id = scatflav(j);
      o.status = -1;
      o.p[0] = scatcen(j, 4);
      o.p[1] = scatcen(j, 1);
      o.p[2] = scatcen(j, 2);
      o.p[3] = scatcen(j, 3);
      o.mass = scatcen(j, 5);
      o.x[0] = scatpos(j, 4);
      o.x[1] = scatpos(j, 1);
      o.x[2] = scatpos(j, 2);
      o.x[3] = scatpos(j, 3);
      out.push_back(o);
    }
  }
  g_stats.n_scatterings += nscatcen;
  return true;
}

const Stats& stats() { return g_stats; }

}  // namespace jewel::api
