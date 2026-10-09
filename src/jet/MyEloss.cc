/*************************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 *
 *************************************************************************************/

#include "MyEloss.h"
#include "EMDQhat.h"
#include "HTLQhat.h"
#include "FluidDynamics.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

#include <array>
#include <vector>
#include <cmath>
#include <string>

#define MAGENTA "\033[35m"

// Register the module with the base class
RegisterJetScapeModule<MyEloss> MyEloss::reg("MyEloss");

namespace {
const double CA = 3.;
const double CF = 4. / 3.;
// Keep the emission probability per sub-step below this value.
const double maxProbPerStep = 0.05;
}  // namespace

MyEloss::MyEloss() {
  SetId("MyEloss");
  VERBOSE(8);
  SetMutex(make_shared<MyElossMutex>());
}

MyEloss::~MyEloss() {
  VERBOSE(8);
#ifdef MYELOSS_DIAG
  std::cerr << "MYELOSS_DIAG acted=" << n_acted << " deferred_to_other_module="
            << n_deferred << std::endl;
#endif
}

void MyEloss::InitTask() {
  JSINFO << "Initialize MyEloss (BDMPS-Z + EMD qhat) ...";

  tStart = GetXMLElementDouble({"Eloss", "tStart"});
  in_vac = GetXMLElementInt({"Eloss", "MyEloss", "in_vac"});
  Q0 = GetXMLElementDouble({"Eloss", "MyEloss", "Q0"});
  T0 = GetXMLElementDouble({"Eloss", "MyEloss", "T0"});
  lambda_t = GetXMLElementDouble({"Eloss", "MyEloss", "lambda_t"});
  qhat_K = GetXMLElementDouble({"Eloss", "MyEloss", "qhat_K"});
  conformal_Tref = GetXMLElementDouble({"Eloss", "MyEloss", "conformal_Tref"});
  qhat_model = GetXMLElementInt({"Eloss", "MyEloss", "qhat_model"});
  htl_alphas = GetXMLElementDouble({"Eloss", "MyEloss", "htl_alphas"});
  if (qhat_model != 0 && qhat_model != 1)
    throw std::runtime_error("MyEloss: qhat_model must be 0 (EMD) or 1 (HTL)");
  alphas = GetXMLElementDouble({"Eloss", "MyEloss", "alphas"});
  omega_min = GetXMLElementDouble({"Eloss", "MyEloss", "omega_min"});
  z_max = GetXMLElementDouble({"Eloss", "MyEloss", "z_max"});
  theta_max = GetXMLElementDouble({"Eloss", "MyEloss", "theta_max"});
  E_cut = GetXMLElementDouble({"Eloss", "MyEloss", "E_cut"});
  radiation_on = GetXMLElementInt({"Eloss", "MyEloss", "radiation_on"});
  broadening_on = GetXMLElementInt({"Eloss", "MyEloss", "broadening_on"});
  formation_time_on =
      GetXMLElementInt({"Eloss", "MyEloss", "formation_time_on"});
  recoil_on = GetXMLElementInt({"Eloss", "MyEloss", "recoil_on"});
  recoil_Q0 = GetXMLElementDouble({"Eloss", "MyEloss", "recoil_Q0"});
  el_alphas = GetXMLElementDouble({"Eloss", "MyEloss", "el_alphas"});
  if (recoil_on)
    BuildKickTable();

  JSINFO << MAGENTA << "MyEloss: qhat_model="
         << (qhat_model == 0 ? "EMD" : "HTL") << " htl_alphas=" << htl_alphas;
  JSINFO << MAGENTA << "MyEloss: in_vac=" << in_vac << " Q0=" << Q0
         << (Q0 < 0 ? " (dynamical)" : "") << " T0=" << T0
         << " lambda_t=" << lambda_t << " qhat_K=" << qhat_K
         << " conformal_Tref=" << conformal_Tref
         << " alphas=" << alphas << " omega_min=" << omega_min
         << " z_max=" << z_max << " theta_max=" << theta_max
         << " E_cut=" << E_cut << " radiation=" << radiation_on
         << " broadening=" << broadening_on
         << " formation_time=" << formation_time_on << " recoil=" << recoil_on
         << " recoil_Q0=" << recoil_Q0 << " el_alphas=" << el_alphas;
  JSINFO << MAGENTA
         << "MyEloss: gluon qhat [GeV^2/fm] at T=0.2/0.3/0.4 GeV, E=100 GeV = "
         << QhatGluon(0.2, 100.) / 0.1973 << " / "
         << QhatGluon(0.3, 100.) / 0.1973 << " / "
         << QhatGluon(0.4, 100.) / 0.1973;
}

void MyEloss::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  auto f = w.lock();
  if (!f)
    return;
  f->WriteComment("ElossModule Parton List: " + GetId());
}

void MyEloss::ClearTask() {}

double MyEloss::QhatGluon(double T, double E) const {
  if (qhat_model == 1)
    return qhat_K * HTLQhat::GluonQhat(T, E, htl_alphas);
  return qhat_K * EMDQhat::GluonQhat(T, lambda_t, conformal_Tref);
}

double MyEloss::OmegaMax(double E, double M2, double qhat_g, double L) const {
  double wmax = std::min(z_max * E, E - std::sqrt(M2));
  if (formation_time_on)
    wmax = std::min(wmax, 0.5 * qhat_g * L * L);  // t_f(omega) <= L
  return wmax;
}

bool MyEloss::TryEmission(double P[3], double M2, double CR, double qhat_g,
                          double L, double dt, double k[4]) {
  auto &rng = *GetMt19937Generator();
  double pmod = std::sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
  double E = std::sqrt(pmod * pmod + M2);
  double wmax = OmegaMax(E, M2, qhat_g, L);
  if (wmax <= omega_min || qhat_g <= 0.)
    return false;

  // Majorant: theta-integrated kernel with theta_max -> infinity,
  // (alpha_s C_R / pi) sqrt(qhat) omega^{-3/2}. Thinning below makes the
  // accepted emissions follow the truncated kernel exactly.
  double umin = 1. / std::sqrt(wmax), umax = 1. / std::sqrt(omega_min);
  double rate = alphas * CR / M_PI * std::sqrt(qhat_g) * 2. * (umax - umin);
  if (ZeroOne(rng) >= rate * dt)
    return false;

  // omega^{-3/2} <=> uniform in omega^{-1/2}
  double u = umin + ZeroOne(rng) * (umax - umin);
  double omega = 1. / (u * u);

  // Angular part theta^3 exp(-a theta^2), a = omega^2 / sqrt(2 omega qhat):
  // x = a theta^2 follows x exp(-x), i.e. Gamma(2,1), truncated at x_max.
  double a = omega * omega / std::sqrt(2. * omega * qhat_g);
  double xmax = a * theta_max * theta_max;
  double accept = 1. - std::exp(-xmax) * (1. + xmax);
  if (ZeroOne(rng) >= accept)
    return false;
  double x;
  do {
    x = -std::log((1. - ZeroOne(rng)) * (1. - ZeroOne(rng)));
  } while (x > xmax);
  double theta = std::sqrt(x / a);
  double phi = 2. * M_PI * ZeroOne(rng);

  double n[3] = {P[0] / pmod, P[1] / pmod, P[2] / pmod};
  double e1[3], e2[3];
  OrthonormalBasis(n, e1, e2);
  double st = std::sin(theta), ct = std::cos(theta);
  for (int a3 = 0; a3 < 3; a3++)
    k[a3 + 1] = omega * (st * std::cos(phi) * e1[a3] +
                         st * std::sin(phi) * e2[a3] + ct * n[a3]);
  k[0] = omega;

  // Parent: energy conserved locally (E -> E - omega), recoils against the
  // gluon's transverse momentum; the medium absorbs the momentum mismatch.
  double Enew = E - omega;
  if (Enew * Enew <= M2)
    return false;
  double Pn[3] = {P[0] - k[1], P[1] - k[2], P[2] - k[3]};
  double pn = std::sqrt(Pn[0] * Pn[0] + Pn[1] * Pn[1] + Pn[2] * Pn[2]);
  if (pn <= 0.)
    return false;
  double scale = std::sqrt(Enew * Enew - M2) / pn;
  for (int a3 = 0; a3 < 3; a3++)
    P[a3] = Pn[a3] * scale;
  return true;
}

void MyEloss::Broaden(double P[3], double qhat_R, double dt) {
  // <k_T^2> = qhat_R dt, Gaussian, transverse to the parton; |P| is kept so
  // the kick is a pure rotation (energy conserving in the fluid frame).
  auto &rng = *GetMt19937Generator();
  double pmod = std::sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
  if (pmod <= 0. || qhat_R <= 0.)
    return;
  double n[3] = {P[0] / pmod, P[1] / pmod, P[2] / pmod};
  double e1[3], e2[3];
  OrthonormalBasis(n, e1, e2);
  double sigma = std::sqrt(0.5 * qhat_R * dt);
  double kx = sigma * Gauss(rng), ky = sigma * Gauss(rng);
  double Pn[3];
  for (int a3 = 0; a3 < 3; a3++)
    Pn[a3] = P[a3] + kx * e1[a3] + ky * e2[a3];
  double pn = std::sqrt(Pn[0] * Pn[0] + Pn[1] * Pn[1] + Pn[2] * Pn[2]);
  for (int a3 = 0; a3 < 3; a3++)
    P[a3] = Pn[a3] * pmod / pn;
}

// In units of T the elastic kinematics depend only on x = E/T (m_D is
// proportional to T), so the correction is tabulated once in ln x.
void MyEloss::BuildKickTable() {
  const int nf = 3, nx = 30, nsample = 4000;
  const double mu2 = 4. * M_PI * el_alphas * (1. + nf / 6.);  // T = 1
  std::mt19937 gen(12345);
  std::uniform_real_distribution<double> U(0., 1.);
  kick_lnx.clear();
  kick_corr.clear();
  for (int ix = 0; ix < nx; ix++) {
    double lnx = std::log(3.) + ix * (std::log(1.e4) - std::log(3.)) / (nx - 1);
    double E = std::exp(lnx);
    double sum_kt2 = 0.;
    int n = 0;
    for (int is = 0; is < nsample; is++) {
      double k = -std::log((1. - U(gen)) * (1. - U(gen)) * (1. - U(gen)));
      double ct = 2. * U(gen) - 1., st = std::sqrt(1. - ct * ct);
      double ph = 2. * M_PI * U(gen);
      double kt[4] = {k, k * st * std::cos(ph), k * st * std::sin(ph), k * ct};
      double tot[4] = {E + kt[0], kt[1], kt[2], E + kt[3]};  // jet along z
      double s = tot[0] * tot[0] - tot[1] * tot[1] - tot[2] * tot[2] - tot[3] * tot[3];
      if (s <= mu2)
        continue;
      double u = U(gen);
      double inv = 0.5 / mu2 - u * (0.5 / mu2 - 1. / (s + mu2));
      double absT = 1. / inv - mu2;
      double beta[3] = {tot[1] / tot[0], tot[2] / tot[0], tot[3] / tot[0]};
      double mbeta[3] = {-beta[0], -beta[1], -beta[2]};
      double pc[4] = {E, 0., 0., E};
      Boost(beta, pc);
      double pcm = std::sqrt(pc[1] * pc[1] + pc[2] * pc[2] + pc[3] * pc[3]);
      double cth = 1. - absT / (2. * pcm * pcm);
      if (cth < -1.)
        continue;
      double sth = std::sqrt(std::max(0., 1. - cth * cth));
      double nn[3] = {pc[1] / pcm, pc[2] / pcm, pc[3] / pcm}, e1[3], e2[3];
      OrthonormalBasis(nn, e1, e2);
      double phi = 2. * M_PI * U(gen);
      double pn[4] = {pc[0], 0., 0., 0.};
      for (int a3 = 0; a3 < 3; a3++)
        pn[a3 + 1] = pcm * (cth * nn[a3] + sth * (std::cos(phi) * e1[a3] +
                                                  std::sin(phi) * e2[a3]));
      Boost(mbeta, pn);
      sum_kt2 += pn[1] * pn[1] + pn[2] * pn[2];
      n++;
    }
    // <|t|> implied by the rate normalisation in ElasticSplit (T = 1)
    double X = std::max(3.70 * E, 2. * mu2);
    double num = std::log((X + mu2) / (2. * mu2)) + mu2 / (X + mu2) - 0.5;
    double den = 0.5 / mu2 - 1. / (X + mu2);
    double corr = (n > 0 && sum_kt2 > 0.) ? (num / den) / (sum_kt2 / n) : 1.;
    kick_lnx.push_back(lnx);
    kick_corr.push_back(corr);
  }
  JSINFO << MAGENTA << "MyEloss: hard-rate correction <|t|>/<q_T^2> at E/T = 3, "
         << "100, 1e4: " << kick_corr.front() << ", "
         << KickCorrection(100., 1.) << ", " << kick_corr.back();
}

double MyEloss::KickCorrection(double E, double T) const {
  if (kick_lnx.empty() || T <= 0.)
    return 1.;
  double lx = std::log(E / T);
  if (lx <= kick_lnx.front())
    return kick_corr.front();
  if (lx >= kick_lnx.back())
    return kick_corr.back();
  size_t i = 0;
  while (i + 2 < kick_lnx.size() && kick_lnx[i + 1] < lx)
    i++;
  double f = (lx - kick_lnx[i]) / (kick_lnx[i + 1] - kick_lnx[i]);
  return kick_corr[i] + f * (kick_corr[i + 1] - kick_corr[i]);
}

void MyEloss::ElasticSplit(double E, double T, double qhat_R,
                           double &qhat_soft, double &gamma_hard, double &mu2,
                           double &qmax2) const {
  const int nf = 3;
  mu2 = 4. * M_PI * el_alphas * T * T * (1. + nf / 6.);
  // Each scattering can transfer at most |t| = s = 2 E k (1 - cos) to its
  // thermal partner. The normalisation uses the typical reach: for a
  // Boltzmann partner (<ln k> = ln T + psi(3)) and isotropic angles
  // (<ln(1-cos)> = ln 2 - 1), <ln s> = ln(2 E T) + 0.6159, i.e. s ~ 3.70 E T.
  qmax2 = std::max(3.70 * E * T, 2. * mu2);
  // int_0^X x A/(x+a)^2 dx = A [ln((X+a)/a) + a/(X+a) - 1]
  double I = std::log((qmax2 + mu2) / mu2) + mu2 / (qmax2 + mu2) - 1.;
  double A = qhat_R / I;
  qhat_soft = A * (std::log(2.) - 0.5);
  // rate such that the jet's <q_T^2> per unit length is the hard part of qhat
  gamma_hard = A * (0.5 / mu2 - 1. / (qmax2 + mu2)) * KickCorrection(E, T);
}

bool MyEloss::HardElastic(double P[3], double M2, double T, double mu2,
                          double qmax2, double kRec[4], double kHole[4],
                          int &pid, double &absT) {
  auto &rng = *GetMt19937Generator();
  double E = std::sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2] + M2);
  for (int tries = 0; tries < 20; tries++) {
    // thermal partner: Boltzmann k^2 exp(-k/T), isotropic
    double k = -T * std::log((1. - ZeroOne(rng)) * (1. - ZeroOne(rng)) *
                             (1. - ZeroOne(rng)));
    double ct = 2. * ZeroOne(rng) - 1., st = std::sqrt(1. - ct * ct);
    double ph = 2. * M_PI * ZeroOne(rng);
    double kt[4] = {k, k * st * std::cos(ph), k * st * std::sin(ph), k * ct};
    double tot[4] = {E + k, P[0] + kt[1], P[1] + kt[2], P[2] + kt[3]};
    double s = tot[0] * tot[0] - tot[1] * tot[1] - tot[2] * tot[2] -
               tot[3] * tot[3];
    if (s <= M2)
      continue;
    // |t| from A/(x+mu2)^2 on [mu2, tmax], tmax = this partner's reach
    double tmax = s - M2;
    if (tmax <= mu2)
      continue;
    double u = ZeroOne(rng);
    double inv = 0.5 / mu2 - u * (0.5 / mu2 - 1. / (tmax + mu2));
    absT = 1. / inv - mu2;
    // CM frame of jet + thermal parton
    double beta[3] = {tot[1] / tot[0], tot[2] / tot[0], tot[3] / tot[0]};
    double mbeta[3] = {-beta[0], -beta[1], -beta[2]};
    double pc[4] = {E, P[0], P[1], P[2]};
    Boost(beta, pc);
    double pcm = std::sqrt(pc[1] * pc[1] + pc[2] * pc[2] + pc[3] * pc[3]);
    double cth = 1. - absT / (2. * pcm * pcm);
    if (cth < -1.)
      continue;  // |t| not reachable with this partner
    double sth = std::sqrt(std::max(0., 1. - cth * cth));
    double n[3] = {pc[1] / pcm, pc[2] / pcm, pc[3] / pcm}, e1[3], e2[3];
    OrthonormalBasis(n, e1, e2);
    double phi = 2. * M_PI * ZeroOne(rng);
    double pn[4] = {pc[0], 0., 0., 0.}, kn[4];
    for (int a3 = 0; a3 < 3; a3++)
      pn[a3 + 1] = pcm * (cth * n[a3] + sth * (std::cos(phi) * e1[a3] +
                                               std::sin(phi) * e2[a3]));
    double Ecm = std::sqrt(s);
    kn[0] = Ecm - pc[0];
    for (int a3 = 0; a3 < 3; a3++)
      kn[a3 + 1] = -pn[a3 + 1];
    Boost(mbeta, pn);
    Boost(mbeta, kn);
    for (int a3 = 0; a3 < 3; a3++)
      P[a3] = pn[a3 + 1];
    for (int a3 = 0; a3 < 4; a3++) {
      kRec[a3] = kn[a3];
      kHole[a3] = kt[a3];
    }
    // thermal species: gluon vs light (anti)quark by Boltzmann degeneracy,
    // g_g = 16, g_q = 4 N_c N_f * 3/4 for nf = 3
    double fg = 16. / (16. + 27.);
    if (ZeroOne(rng) < fg)
      pid = 21;
    else {
      int fl = 1 + (int)(3. * ZeroOne(rng));
      pid = ZeroOne(rng) < 0.5 ? fl : -fl;
    }
    return true;
  }
  return false;
}

void MyEloss::Boost(const double v[3], double p[4]) {
  double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  if (v2 < 1e-12)
    return;
  if (v2 >= 1.)
    v2 = 1. - 1e-10;
  double gamma = 1. / std::sqrt(1. - v2);
  double vp = v[0] * p[1] + v[1] * p[2] + v[2] * p[3];
  double c = (gamma - 1.) * vp / v2 - gamma * p[0];
  p[0] = gamma * (p[0] - vp);
  for (int a3 = 0; a3 < 3; a3++)
    p[a3 + 1] += c * v[a3];
}

void MyEloss::OrthonormalBasis(const double n[3], double e1[3], double e2[3]) {
  // e1 perpendicular to n, built from the coordinate axis least aligned with n.
  double ref[3] = {0., 0., 0.};
  int imin = 0;
  for (int a3 = 1; a3 < 3; a3++)
    if (std::abs(n[a3]) < std::abs(n[imin]))
      imin = a3;
  ref[imin] = 1.;
  double dot = ref[0] * n[0] + ref[1] * n[1] + ref[2] * n[2];
  for (int a3 = 0; a3 < 3; a3++)
    e1[a3] = ref[a3] - dot * n[a3];
  double m = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  for (int a3 = 0; a3 < 3; a3++)
    e1[a3] /= m;
  e2[0] = n[1] * e1[2] - n[2] * e1[1];
  e2[1] = n[2] * e1[0] - n[0] * e1[2];
  e2[2] = n[0] * e1[1] - n[1] * e1[0];
}

void MyEloss::DoEnergyLoss(double deltaT, double time, double Q2,
                           vector<Parton> &pIn, vector<Parton> &pOut) {
  if (in_vac)
    return;

  GetHydroTau0Signal(tStart);

  for (int i = 0; i < pIn.size(); i++) {
    if (pIn[i].pid() == photonid || pIn[i].pstat() < 0)
      continue;

    // Propagate from the last update point to the current time.
    double E = pIn[i].e();
    double p[4] = {E, pIn[i].px(), pIn[i].py(), pIn[i].pz()};
    double pmod = std::sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]);
    if (pmod < 1e-6)
      continue;
    double w[3] = {p[1] / E, p[2] / E, p[3] / E};
    double t0 = pIn[i].x_in().t();
    double x[4] = {time, pIn[i].x_in().x() + (time - t0) * w[0],
                   pIn[i].x_in().y() + (time - t0) * w[1],
                   pIn[i].x_in().z() + (time - t0) * w[2]};
    if (x[0] * x[0] <= x[3] * x[3])
      continue;
    double tau = std::sqrt(x[0] * x[0] - x[3] * x[3]);
    if (tau < tStart)
      continue;

    std::unique_ptr<FluidCellInfo> cell;
    GetHydroCellSignal(x[0], x[1], x[2], x[3], cell);
    if (!GetJetSignalConnected()) {
      JSWARN << "Couldn't find a hydro module attached!";
      throw std::runtime_error(
          "Please attach a hydro module or set in_vac to 1 in the XML file");
    }
    double T = cell->temperature;
    if (T < T0)
      continue;

    bool isGluon = (pIn[i].pid() == 21);
    double CR = isGluon ? CA : CF;
    double qhat_g = QhatGluon(T, E);
    double qhat_R = qhat_g * CR / CA;

    // Hand-off from MATTER: act only below the switching virtuality.
    double Q0loc = Q0;
    if (Q0 < 0.)
      Q0loc = std::max(1., std::sqrt(std::sqrt(2. * E * qhat_R * std::sqrt(2.))));
    if (pIn[i].t() > Q0loc * Q0loc + rounding_error)
      continue;

    // Fluid rest frame
    double v[3] = {cell->vx, cell->vy, cell->vz};
    double mv[3] = {-v[0], -v[1], -v[2]};
    double pl[4] = {p[0], p[1], p[2], p[3]};
    Boost(v, pl);
    if (pl[0] < E_cut)
      continue;

    // MATTER runs first and may still own this parton: it samples the medium
    // at the parton's last vertex, we at its current position.
    if (pIn[i].GetControlled()) {
#ifdef MYELOSS_DIAG
      n_deferred++;
#endif
      continue;
    }
#ifdef MYELOSS_DIAG
    n_acted++;
#endif

    TakeResponsibilityFor(pIn[i]);

    double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    double gamma = 1. / std::sqrt(std::max(1e-10, 1. - v2));
    double dt_lrf =
        deltaT * gamma * (1. - (v[0] * w[0] + v[1] * w[1] + v[2] * w[2]));
    double dt = dt_lrf * fmToGeVinv;

    double l_med = 0.;
    if (pIn[i].has_user_info<MyElossUserInfo>())
      l_med = pIn[i].user_info<MyElossUserInfo>().l_med();
    double L0 = l_med * fmToGeVinv;

    double P[3] = {pl[1], pl[2], pl[3]};
    double M2 = pl[0] * pl[0] - P[0] * P[0] - P[1] * P[1] - P[2] * P[2];
    if (M2 < 0.) {
      // Space-like input: put it on shell at fixed energy.
      double s = pl[0] / std::sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
      for (int a3 = 0; a3 < 3; a3++)
        P[a3] *= s;
      M2 = 0.;
    }

    // Elastic split of qhat_R into soft broadening and hard scatterings.
    double qhat_soft = qhat_R, gamma_hard = 0., mu2 = 0., qmax2 = 0.;
    if (broadening_on && recoil_on)
      ElasticSplit(pl[0], T, qhat_R, qhat_soft, gamma_hard, mu2, qmax2);

    // Split the step so that each emission or scattering trial has a small
    // probability.
    int nsub = 1;
    double rate = gamma_hard;
    if (radiation_on) {
      double wmax = OmegaMax(pl[0], M2, qhat_g, L0 + dt);
      if (wmax > omega_min)
        rate += alphas * CR / M_PI * std::sqrt(qhat_g) * 2. *
                (1. / std::sqrt(omega_min) - 1. / std::sqrt(wmax));
    }
    nsub = std::max(1, (int)std::ceil(rate * dt / maxProbPerStep));
    double dts = dt / nsub;

    vector<std::array<double, 4>> gluons;
    // hard recoils kept as jet partons: {E, px, py, pz, |t|, pid}; holes
    vector<std::array<double, 6>> recoils;
    vector<std::array<double, 5>> holes;
    for (int s = 0; s < nsub; s++) {
      if (radiation_on) {
        double k[4];
        if (TryEmission(P, M2, CR, qhat_g, L0 + (s + 1) * dts, dts, k))
          gluons.push_back({k[0], k[1], k[2], k[3]});
      }
      if (broadening_on)
        Broaden(P, qhat_soft, dts);
      if (gamma_hard > 0. &&
          ZeroOne(*GetMt19937Generator()) < gamma_hard * dts) {
        double kr[4], kh[4], absT;
        int pid;
        if (HardElastic(P, M2, T, mu2, qmax2, kr, kh, pid, absT) &&
            absT > recoil_Q0 * recoil_Q0) {
          recoils.push_back({kr[0], kr[1], kr[2], kr[3], absT, (double)pid});
          holes.push_back({kh[0], kh[1], kh[2], kh[3], (double)pid});
        }
        // otherwise the recoil and the hole become part of the medium
      }
    }

    // Back to the lab frame and hand the partons to the framework.
    double fx[4] = {x[0], x[1], x[2], x[3]};
    double velocity_jet[4] = {1.0, pIn[i].jet_v().x(), pIn[i].jet_v().y(),
                              pIn[i].jet_v().z()};

    double pf[4] = {std::sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2] + M2),
                    P[0], P[1], P[2]};
    Boost(mv, pf);
    FourVector xVec;
    pOut.push_back(Parton(pIn[i].plabel(), pIn[i].pid(), pIn[i].pstat(),
                          FourVector(pf[1], pf[2], pf[3], pf[0]), xVec));
    pOut.back().set_x(fx);
    pOut.back().set_color(pIn[i].color());
    pOut.back().set_anti_color(pIn[i].anti_color());
    pOut.back().set_user_info(new MyElossUserInfo(l_med + dt_lrf));
    pOut.back().set_jet_v(velocity_jet);
    pOut.back().set_mean_form_time();
    pOut.back().set_form_time(pOut.back().mean_form_time());

    for (auto &g : gluons) {
      double kf[4] = {g[0], g[1], g[2], g[3]};
      Boost(mv, kf);
      pOut.push_back(
          Parton(0, 21, 0, FourVector(kf[1], kf[2], kf[3], kf[0]), xVec));
      pOut.back().set_x(fx);
      pOut.back().set_user_info(new MyElossUserInfo(0.));
      pOut.back().set_jet_v(velocity_jet);
      pOut.back().set_mean_form_time();
      pOut.back().set_form_time(pOut.back().mean_form_time());
    }

    // The shower graph only records partons when a module returns two or
    // more of them, so a step without emission or kept recoil also returns
    // the momentum given to the medium as a status -13 ("missing") parton,
    // as AdSCFT does. Otherwise the broadening kick (and the elastic energy
    // loss) after a parton's last emission would never reach the output.
    // Hadronisation ignores status -13.
    if (gluons.empty() && recoils.empty()) {
      FourVector pVecM(pIn[i].px() - pf[1], pIn[i].py() - pf[2],
                       pIn[i].pz() - pf[3], pIn[i].e() - pf[0]);
      pOut.push_back(Parton(0, 21, -13, pVecM, xVec));
      pOut.back().set_x(fx);
    }

    // Hard recoils: jet partons (status 0) with virtuality |t| at fixed
    // energy, so MATTER showers them while |t| > Q0^2 and MyEloss takes over
    // below. Holes: status -1, as for LBT (subtracted in jet observables).
    for (size_t r = 0; r < recoils.size(); r++) {
      auto &k = recoils[r];
      double kr[4] = {k[0], k[1], k[2], k[3]};
      double pm = std::sqrt(kr[1] * kr[1] + kr[2] * kr[2] + kr[3] * kr[3]);
      if (kr[0] * kr[0] > k[4] && pm > 0.) {
        double sc = std::sqrt(kr[0] * kr[0] - k[4]) / pm;
        for (int a3 = 1; a3 < 4; a3++)
          kr[a3] *= sc;
      }
      Boost(mv, kr);
      pOut.push_back(Parton(0, (int)k[5], 0,
                            FourVector(kr[1], kr[2], kr[3], kr[0]), xVec));
      pOut.back().set_x(fx);
      pOut.back().set_user_info(new MyElossUserInfo(0.));
      pOut.back().set_jet_v(velocity_jet);
      pOut.back().set_mean_form_time();
      pOut.back().set_form_time(pOut.back().mean_form_time());

      auto &hh = holes[r];
      double kh[4] = {hh[0], hh[1], hh[2], hh[3]};
      Boost(mv, kh);
      pOut.push_back(Parton(0, (int)hh[4], -1,
                            FourVector(kh[1], kh[2], kh[3], kh[0]), xVec));
      pOut.back().set_x(fx);
      pOut.back().set_jet_v(velocity_jet);
    }
  }
}
