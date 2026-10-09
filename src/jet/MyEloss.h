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

// MyEloss: low-virtuality BDMPS-Z ("BDMPS-lite") medium-induced radiation with
// a temperature dependent holographic (EMD) qhat, plus transverse momentum
// broadening of every parton it handles. Meant to run after MATTER, taking
// over partons once their virtuality drops below Q0 (fixed, or the dynamical
// Q0^2 = sqrt(2 sqrt(2) E qhat_R) used by MATTER when Q0 < 0).
//
// Emission kernel per unit length in the fluid rest frame:
//   dN / (domega dtheta dL) = alpha omega theta^3 t_f exp(-theta^2 omega^2 /
//   (qhat t_f)),   t_f = sqrt(2 omega / qhat)
// Integrating over theta gives the BDMPS spectrum
//   omega dN / (domega dL) = sqrt(2) alpha sqrt(qhat / omega),
// so alpha = alpha_s C_R / (sqrt(2) pi) reproduces
//   omega dN / (domega dL) = (alpha_s C_R / pi) sqrt(qhat_A / omega).

#ifndef MYELOSS_H
#define MYELOSS_H

#include "JetEnergyLossModule.h"
#include "JetScapeModuleMutex.h"
#include <random>
#include <vector>

using namespace Jetscape;

// Path length [fm, fluid rest frame] travelled while handled by MyEloss.
// Enters the formation-time constraint t_f(omega) <= L, i.e.
// omega <= omega_c = qhat L^2 / 2.
class MyElossUserInfo : public fjcore::PseudoJet::UserInfoBase {
 public:
  MyElossUserInfo(double l_med) : _l_med(l_med){};
  double l_med() const { return _l_med; }
  double _l_med;
  ~MyElossUserInfo(){};
};

// MyEloss handles the same low-virtuality phase as LBT, MARTINI and AdS/CFT.
class MyElossMutex : public JetScapeModuleMutex {
 public:
  bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules) {
    for (auto module : modules) {
      string name = module->GetId();
      if (name == "LBT" || name == "Martini" || name == "AdSCFT")
        return false;
    }
    return true;
  }
};

class MyEloss : public JetEnergyLossModule<MyEloss> {
 public:
  MyEloss();
  virtual ~MyEloss();

  void InitTask();
  void ClearTask();
  void WriteTask(weak_ptr<JetScapeWriter> w);

  void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton> &pIn,
                    vector<Parton> &pOut);

 private:
  // Gluon qhat [GeV^3] at temperature T [GeV] for a parton of energy E [GeV]
  // (E only enters the HTL model).
  double QhatGluon(double T, double E) const;
  // Upper edge of the emission-energy window [GeV] for a parton of energy E
  // and mass^2 M2 with in-medium path length L [GeV^-1].
  double OmegaMax(double E, double M2, double qhat_g, double L) const;
  // One emission trial over dt [GeV^-1]; on success fills k (lrf 4-momentum
  // of the gluon) and updates the parent's 3-momentum P.
  bool TryEmission(double P[3], double M2, double CR, double qhat_g, double L,
                   double dt, double k[4]);
  void Broaden(double P[3], double qhat_R, double dt);

  // Elastic recoil option. The transverse-kick distribution
  //   dGamma/dq^2 = A / (q^2 + mu^2)^2,  0 < q^2 < qmax^2 ~ 3.7 E T,
  // mu^2 = m_D^2, is normalised to the total qhat_R. Kicks below mu stay
  // Gaussian broadening (qhat_soft); kicks above mu are discrete elastic
  // scatterings off thermal partons with rate gamma_hard.
  void ElasticSplit(double E, double T, double qhat_R, double &qhat_soft,
                    double &gamma_hard, double &mu2, double &qmax2) const;
  // One hard elastic scattering in the fluid frame: a thermal partner is
  // drawn first, then |t| from the tail A/(q^2+mu^2)^2 between mu^2 and that
  // partner's reach s. Updates P, returns recoil and hole 4-momenta, the
  // thermal parton's id and |t|. False if no partner with s > mu^2 was found.
  // (qmax2 sets only the rate normalisation; see ElasticSplit.)
  bool HardElastic(double P[3], double M2, double T, double mu2, double qmax2,
                   double kRec[4], double kHole[4], int &pid, double &absT);

  static void Boost(const double v[3], double p[4]);  // lab -> frame moving with v
  static void OrthonormalBasis(const double n[3], double e1[3], double e2[3]);

  bool in_vac;
  double tStart;         // start of quenching [fm/c]
  double Q0;             // switching virtuality [GeV]; < 0 -> dynamical
  double T0;             // no quenching below T0 [GeV]
  double lambda_t;       // 't Hooft coupling in qhat = Q(T) sqrt(lambda) T^3
  double qhat_K;         // overall qhat multiplier
  double conformal_Tref; // > 0: conformal control, Q(T) frozen at this T [GeV]
  int qhat_model;        // 0: EMD holographic, 1: JETSCAPE HTL (as in MATTER)
  double htl_alphas;     // fixed alpha_s of the HTL qhat
  double alphas;         // fixed alpha_s of the radiation kernel
  double omega_min;      // soft cut on emitted gluon energy [GeV]
  double z_max;          // omega <= z_max E
  double theta_max;      // emission angle cut [rad]
  double E_cut;          // partons below this lrf energy stream freely [GeV]
  bool radiation_on;
  bool broadening_on;
  bool formation_time_on;
  bool recoil_on;        // hard elastic scatterings with recoils and holes
  double recoil_Q0;      // |t| > recoil_Q0^2: recoil kept as jet parton with
                         // virtuality |t| (showered by MATTER); otherwise the
                         // recoil and hole become part of the medium
  double el_alphas;      // alpha_s in the Debye mass m_D^2 = 4 pi a T^2 (1+nf/6)
  // Hard-rate correction vs ln(E/T): (normalisation <|t|>) / <q_T^2 of the
  // jet per scattering>, from a small Monte Carlo at initialisation, so that
  // the jet's transverse momentum squared per unit length equals qhat.
  std::vector<double> kick_lnx, kick_corr;
  void BuildKickTable();
  double KickCorrection(double E, double T) const;

#ifdef MYELOSS_DIAG
  long n_acted = 0, n_deferred = 0;  // parton-steps handled / left to MATTER
#endif

  std::uniform_real_distribution<double> ZeroOne{0., 1.};
  std::normal_distribution<double> Gauss{0., 1.};

  // Allows the registration of the module so that it is available to be used
  // by the Jetscape framework.
  static RegisterJetScapeModule<MyEloss> reg;
};

#endif  // MYELOSS_H
