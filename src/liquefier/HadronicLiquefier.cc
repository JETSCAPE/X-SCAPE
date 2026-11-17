/*******************************************************************************
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
 ******************************************************************************/
// ------------------------------------------------------------
// This is a hadronic liquefier for the JETSCAPE framework
// It implements a Gaussian smearing kernel and a covariant one
// ------------------------------------------------------------

#include "HadronicLiquefier.h"
#include <algorithm>
#include <cmath>
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

namespace Jetscape {
HadronicLiquefier::HadronicLiquefier() { InitializeParameters(); }

HadronicLiquefier::HadronicLiquefier(bool covariant, double sigma_transverse,
                                     double sigma_longitudinal, double xMax,
                                     double yMax, double zMax, int Nx, int Ny,
                                     int Nz, bool hydro_Cartesian) {
  covariant_smearing_ = covariant;
  sigma_transverse_ = sigma_transverse;
  sigma_longitudinal_ = sigma_longitudinal;
  xMax_ = xMax;
  yMax_ = yMax;
  zMax_ = zMax;
  Nx_ = Nx;
  Ny_ = Ny;
  Nz_ = Nz;
  dx_ = 2. * xMax_ / (Nx_ - 1.);
  dy_ = 2. * yMax_ / (Ny_ - 1.);
  dz_ = 2. * zMax_ / (Nz_ - 1.);
  hydro_Cartesian_ = hydro_Cartesian;
  gamma_factor_max_ = 5.0;
}

void HadronicLiquefier::InitializeParameters() {
  int cov_kernel = JetScapeXML::Instance()->GetElementInt(
      {"Liquefier", "HadronicLiquefier", "covariant_kernel"});
  if (cov_kernel != 0 && cov_kernel != 1) {
    JSWARN << "The covariant_kernel parameter is not 0 or 1.";
    exit(1);
  }
  covariant_smearing_ = cov_kernel;

  sigma_transverse_ = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "HadronicLiquefier", "sigma_transverse"});
  sigma_longitudinal_ = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "HadronicLiquefier", "sigma_longitudinal"});

  gamma_factor_max_ = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "HadronicLiquefier", "gamma_factor_max"});

  // get the grid specifications from the initial state module, which are
  // also used in the hydro evolution
  xMax_ = JetScapeXML::Instance()->GetElementDouble({"IS", "grid_max_x"});
  yMax_ = JetScapeXML::Instance()->GetElementDouble({"IS", "grid_max_y"});
  zMax_ = JetScapeXML::Instance()->GetElementDouble({"IS", "grid_max_z"});
  dx_ = JetScapeXML::Instance()->GetElementDouble({"IS", "grid_step_x"});
  dy_ = JetScapeXML::Instance()->GetElementDouble({"IS", "grid_step_y"});
  dz_ = JetScapeXML::Instance()->GetElementDouble({"IS", "grid_step_z"});

  Nx_ = int(std::ceil(2 * xMax_ / dx_));
  Ny_ = int(std::ceil(2 * yMax_ / dy_));
  Nz_ = int(std::ceil(2 * zMax_ / dz_));

  if (Nz_ == 1 || dz_ == 0.) {
    JSWARN << "HadronicLiquefier: The longitudinal grid range is not set up "
              "properly.";
  }

  hydro_Cartesian_ = false;
  std::string strCartesianHydro =
      JetScapeXML::Instance()->GetElementText({"Hydro", "CartesianHydro"});
  if ((int)strCartesianHydro.find("true") >= 0) {
    hydro_Cartesian_ = true;
    JSINFO << "BulkDynamicsManager set up for run with Cartesian hydro ...";
  }
}

double HadronicLiquefier::smearing_kernel_covariant_Milne(
    const double x_diff, const double y_diff, const double eta_diff,
    const double ux, const double uy, const double ueta, const double tau,
    const double gamma) const {
  const double N = (tau * gamma) / (pow(M_PI, 1.5) * sigma_transverse_ *
                                    sigma_transverse_ * sigma_transverse_);
  // compute the squared distance and the scalar product r*u
  const double dr_squared =
      (x_diff * x_diff + y_diff * y_diff + eta_diff * eta_diff * tau * tau);
  const double dr_dot_u = (x_diff * ux + y_diff * uy + eta_diff * ueta * tau);

  return N * exp(-(dr_squared + dr_dot_u * dr_dot_u) /
                 (sigma_transverse_ * sigma_transverse_));
}

double HadronicLiquefier::smearing_kernel_covariant_Cartesian(
    const double x_diff, const double y_diff, const double z_diff,
    const double ux, const double uy, const double uz,
    const double gamma) const {
  const double N = gamma / (pow(M_PI, 1.5) * sigma_transverse_ *
                            sigma_transverse_ * sigma_transverse_);
  // compute the squared distance and the scalar product r*u
  const double dr_squared =
      (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
  const double dr_dot_u = (x_diff * ux + y_diff * uy + z_diff * uz);

  return N * exp(-(dr_squared + dr_dot_u * dr_dot_u) /
                 (sigma_transverse_ * sigma_transverse_));
}

double HadronicLiquefier::smearing_kernel_gaussian(
    const double x_diff, const double y_diff, const double eta_diff) const {
  const double N_eta = 1. / (sqrt(M_PI) * sigma_longitudinal_);
  const double N_trans = 1. / (M_PI * sigma_transverse_ * sigma_transverse_);

  const double r_trans_squared = x_diff * x_diff + y_diff * y_diff;
  const double deta_squared = eta_diff * eta_diff;

  const double exp_trans =
      N_trans * exp(-r_trans_squared / (sigma_transverse_ * sigma_transverse_));
  const double exp_eta =
      N_eta * exp(-deta_squared / (sigma_longitudinal_ * sigma_longitudinal_));
  return exp_trans * exp_eta;
}

double HadronicLiquefier::compute_drop_kernel_normalization(
    const double tau, const HadronDroplet &drop_i) const {
  const double skip_dis_x = skip_n_sigma_transverse_ * sigma_transverse_;
  double skip_dis_eta;
  if (covariant_smearing_) {
    // The covariant kernel has only one extension parameter
    skip_dis_eta = skip_n_sigma_transverse_ * sigma_transverse_;
  } else {
    skip_dis_eta = skip_n_sigma_longitudinal_ * sigma_longitudinal_;
  }

  if (covariant_smearing_ && !hydro_Cartesian_) {
    skip_dis_eta = skip_dis_x / tau;
  }
  double normalization = 0.;

  auto xmu_i = drop_i.get_xmu();
  auto pmu_i = drop_i.get_pmu();

  // Only the cells within the skip distances can contribute (the tests
  // below are unchanged): index ranges with one cell of margin, the same
  // cells in the same order as the loop over the whole grid.
  auto range = [](double c, double half, double lo, double d, int n,
                  int &i0, int &i1) {
    i0 = std::max(0, static_cast<int>(std::floor((c - half - lo) / d)) - 1);
    i1 = std::min(n - 1, static_cast<int>(std::ceil((c + half - lo) / d)) + 1);
  };
  int ix0, ix1, iy0, iy1, ie0, ie1;
  range(xmu_i[1], skip_dis_x, -xMax_, dx_, Nx_, ix0, ix1);
  range(xmu_i[2], skip_dis_x, -yMax_, dy_, Ny_, iy0, iy1);
  range(xmu_i[3], skip_dis_eta, -zMax_, dz_, Nz_, ie0, ie1);
  for (int ieta = ie0; ieta <= ie1; ieta++) {
    for (int ix = ix0; ix <= ix1; ix++) {
      for (int iy = iy0; iy <= iy1; iy++) {
        double eta = -zMax_ + ieta * dz_;
        double x = -xMax_ + ix * dx_;
        double y = -yMax_ + iy * dy_;

        double x_diff = x - xmu_i[1];
        if (abs(x_diff) > skip_dis_x) {
          continue;
        }
        double y_diff = y - xmu_i[2];
        if (abs(y_diff) > skip_dis_x) {
          continue;
        }
        const double eta_s = xmu_i[3];  // Cartesian z, Milne eta
        double eta_diff = eta - xmu_i[3];
        if (abs(eta_diff) > skip_dis_eta) {
          continue;
        }

        const double mass = sqrt(pmu_i[0] * pmu_i[0] - pmu_i[1] * pmu_i[1] -
                                 pmu_i[2] * pmu_i[2] - pmu_i[3] * pmu_i[3]);
        if (mass <= 1e-16) {
          continue;  // skip massless particles
        }
        const double mT =
            sqrt(mass * mass + pmu_i[1] * pmu_i[1] + pmu_i[2] * pmu_i[2]);
        double rapidity;
        if (mT >= 1e-16) {
          rapidity = std::asinh(pmu_i[3] / mT);
        } else {
          rapidity = 0.5 * log((pmu_i[0] + pmu_i[3]) / (pmu_i[0] - pmu_i[3]));
        }

        if (covariant_smearing_ && !hydro_Cartesian_) {
          const double ux = pmu_i[1] / mass;
          const double uy = pmu_i[2] / mass;
          const double peta = mT * sinh(rapidity - eta_s);
          const double ueta = peta / mass;
          const double gamma = mT * cosh(rapidity - eta_s) / mass;

          normalization += smearing_kernel_covariant_Milne(
              x_diff, y_diff, eta_diff, ux, uy, ueta, tau, gamma);
        } else if (!covariant_smearing_ && !hydro_Cartesian_) {
          normalization += smearing_kernel_gaussian(x_diff, y_diff, eta_diff);
        } else if (covariant_smearing_ && hydro_Cartesian_) {
          const double ux = pmu_i[1] / mass;
          const double uy = pmu_i[2] / mass;
          const double uz = pmu_i[3] / mass;
          const double gamma = pmu_i[0] / mass;  // Cartesian: the lab frame

          normalization += smearing_kernel_covariant_Cartesian(
              x_diff, y_diff, eta_diff, ux, uy, uz, gamma);
        } else if (!covariant_smearing_ && hydro_Cartesian_) {
          normalization += smearing_kernel_gaussian(x_diff, y_diff, eta_diff);
        }
      }
    }
  }
  return normalization * dx_ * dy_ * dz_;
}

void HadronicLiquefier::get_source_energy(const double tau, const double x,
                                          const double y, const double eta,
                                          std::array<double, 4> &jmu) const {
  jmu = {0., 0., 0., 0.};
  const double skip_dis_x = skip_n_sigma_transverse_ * sigma_transverse_;
  double skip_dis_eta;
  if (covariant_smearing_) {
    // The covariant kernel has only one extension parameter
    skip_dis_eta = skip_n_sigma_transverse_ * sigma_transverse_;
  } else {
    skip_dis_eta = skip_n_sigma_longitudinal_ * sigma_longitudinal_;
  }

  if (covariant_smearing_ && !hydro_Cartesian_) {
    skip_dis_eta = skip_dis_x / tau;
  }

  double value_kernel = 0.;
  // the step's droplets (prepare_active_droplets) or all, in list order
  const std::size_t n_loop =
      active_valid_ ? active_idx_.size() : hadron_droplets_list.size();
  for (std::size_t k_loop = 0; k_loop < n_loop; k_loop++) {
    const auto &drop_i =
        hadron_droplets_list[active_valid_ ? active_idx_[k_loop] : k_loop];
    auto xmu_i = drop_i.get_xmu();
    if (dtau_ > 0 &&
        ((xmu_i[0] < tau - 0.5 * dtau_) || (xmu_i[0] >= tau + 0.5 * dtau_))) {
      continue;
    }

    double x_diff = x - xmu_i[1];
    if (abs(x_diff) > skip_dis_x) {
      continue;
    }
    double y_diff = y - xmu_i[2];
    if (abs(y_diff) > skip_dis_x) {
      continue;
    }
    const double eta_s = xmu_i[3];  // Cartesian z, Milne eta
    double eta_diff = eta - eta_s;
    if (abs(eta_diff) > skip_dis_eta) {
      continue;
    }

    auto pmu_i = drop_i.get_pmu();
    const double mass = sqrt(pmu_i[0] * pmu_i[0] - pmu_i[1] * pmu_i[1] -
                             pmu_i[2] * pmu_i[2] - pmu_i[3] * pmu_i[3]);
    if (mass <= 1e-16) {
      continue;  // skip massless particles
    }
    const double mT =
        sqrt(mass * mass + pmu_i[1] * pmu_i[1] + pmu_i[2] * pmu_i[2]);
    double rapidity;
    if (mT >= 1e-16) {
      rapidity = std::asinh(pmu_i[3] / mT);
    } else {
      rapidity = 0.5 * log((pmu_i[0] + pmu_i[3]) / (pmu_i[0] - pmu_i[3]));
    }

    if (covariant_smearing_ && !hydro_Cartesian_) {
      const double ux = pmu_i[1] / mass;
      const double uy = pmu_i[2] / mass;
      const double peta = mT * sinh(rapidity - eta_s);
      const double ueta = peta / mass;
      const double gamma = mT * cosh(rapidity - eta_s) / mass;

      value_kernel = smearing_kernel_covariant_Milne(x_diff, y_diff, eta_diff,
                                                     ux, uy, ueta, tau, gamma) /
                     tau;
    } else if (!covariant_smearing_ && !hydro_Cartesian_) {
      value_kernel = smearing_kernel_gaussian(x_diff, y_diff, eta_diff) / tau;
    } else if (covariant_smearing_ && hydro_Cartesian_) {
      const double ux = pmu_i[1] / mass;
      const double uy = pmu_i[2] / mass;
      const double uz = pmu_i[3] / mass;
      const double gamma = pmu_i[0] / mass;  // Cartesian: the lab frame

      value_kernel = smearing_kernel_covariant_Cartesian(
          x_diff, y_diff, eta_diff, ux, uy, uz, gamma);
    } else if (!covariant_smearing_ && hydro_Cartesian_) {
      value_kernel = smearing_kernel_gaussian(x_diff, y_diff, eta_diff);
    }
    double norm = drop_i.get_normalization();
    if (norm <= 0.0 || std::isnan(norm)) {
      JSWARN << "Normalization is zero or NaN for a hadron droplet. Skipping.";
      continue;
    }
    value_kernel /= norm;

    if (hydro_Cartesian_) {
      // Cartesian: (t, x, y, z) components, p^mu itself
      jmu[0] += value_kernel * pmu_i[0];
      jmu[1] += value_kernel * pmu_i[1];
      jmu[2] += value_kernel * pmu_i[2];
      jmu[3] += value_kernel * pmu_i[3];
    } else {
      // Milne: p^tau = mT cosh(y - eta), tau p^eta = mT sinh(y - eta)
      jmu[0] += value_kernel * mT * cosh(rapidity - eta);
      jmu[1] += value_kernel * pmu_i[1];
      jmu[2] += value_kernel * pmu_i[2];
      jmu[3] += value_kernel * mT * sinh(rapidity - eta);
    }
  }
}

double HadronicLiquefier::get_source_quantity(const double tau, const double x,
                                              const double y, const double eta,
                                              const QuantityType qtype) const {
  double result = 0.;
  const double skip_dis_x = skip_n_sigma_transverse_ * sigma_transverse_;
  double skip_dis_eta;
  if (covariant_smearing_) {
    // The covariant kernel has only one extension parameter
    skip_dis_eta = skip_n_sigma_transverse_ * sigma_transverse_;
  } else {
    skip_dis_eta = skip_n_sigma_longitudinal_ * sigma_longitudinal_;
  }

  if (covariant_smearing_ && !hydro_Cartesian_) {
    skip_dis_eta = skip_dis_x / tau;
  }

  double value_kernel = 0.;
  // the step's droplets (prepare_active_droplets) or all, in list order
  const std::size_t n_loop =
      active_valid_ ? active_idx_.size() : hadron_droplets_list.size();
  for (std::size_t k_loop = 0; k_loop < n_loop; k_loop++) {
    const auto &drop_i =
        hadron_droplets_list[active_valid_ ? active_idx_[k_loop] : k_loop];
    int quantity_smear = 0;
    if (qtype == BARYON_NUMBER) {
      quantity_smear = drop_i.get_baryon_number();
    } else if (qtype == ELECTRIC_CHARGE) {
      quantity_smear = drop_i.get_electric_charge();
    } else if (qtype == STRANGENESS) {
      quantity_smear = drop_i.get_strangeness();
    } else {
      JSWARN << "The quantity to smear is not implemented.";
      exit(1);
    }
    if (quantity_smear == 0) {
      continue;
    }

    auto xmu_i = drop_i.get_xmu();
    if (dtau_ > 0 &&
        ((xmu_i[0] < tau - 0.5 * dtau_) || (xmu_i[0] >= tau + 0.5 * dtau_))) {
      continue;
    }

    double x_diff = x - xmu_i[1];
    if (abs(x_diff) > skip_dis_x) {
      continue;
    }
    double y_diff = y - xmu_i[2];
    if (abs(y_diff) > skip_dis_x) {
      continue;
    }
    const double eta_s = xmu_i[3];  // Cartesian z, Milne eta
    double eta_diff = eta - xmu_i[3];
    if (abs(eta_diff) > skip_dis_eta) {
      continue;
    }

    auto pmu_i = drop_i.get_pmu();
    const double mass = sqrt(pmu_i[0] * pmu_i[0] - pmu_i[1] * pmu_i[1] -
                             pmu_i[2] * pmu_i[2] - pmu_i[3] * pmu_i[3]);
    if (mass <= 1e-16) {
      continue;  // skip massless particles
    }
    const double mT =
        sqrt(mass * mass + pmu_i[1] * pmu_i[1] + pmu_i[2] * pmu_i[2]);
    double rapidity;
    if (mT >= 1e-16) {
      rapidity = std::asinh(pmu_i[3] / mT);
    } else {
      rapidity = 0.5 * log((pmu_i[0] + pmu_i[3]) / (pmu_i[0] - pmu_i[3]));
    }

    if (covariant_smearing_ && !hydro_Cartesian_) {
      const double ux = pmu_i[1] / mass;
      const double uy = pmu_i[2] / mass;
      const double peta = mT * sinh(rapidity - eta_s);
      const double ueta = peta / mass;
      const double gamma = mT * cosh(rapidity - eta_s) / mass;

      value_kernel = smearing_kernel_covariant_Milne(x_diff, y_diff, eta_diff,
                                                     ux, uy, ueta, tau, gamma) /
                     tau;
    } else if (!covariant_smearing_ && !hydro_Cartesian_) {
      value_kernel = smearing_kernel_gaussian(x_diff, y_diff, eta_diff) / tau;
    } else if (covariant_smearing_ && hydro_Cartesian_) {
      const double ux = pmu_i[1] / mass;
      const double uy = pmu_i[2] / mass;
      const double uz = pmu_i[3] / mass;
      const double gamma = pmu_i[0] / mass;  // Cartesian: the lab frame

      value_kernel = smearing_kernel_covariant_Cartesian(
          x_diff, y_diff, eta_diff, ux, uy, uz, gamma);
    } else if (!covariant_smearing_ && hydro_Cartesian_) {
      value_kernel = smearing_kernel_gaussian(x_diff, y_diff, eta_diff);
    }
    double norm = drop_i.get_normalization();
    if (norm <= 0.0 || std::isnan(norm)) {
      JSWARN << "Normalization is zero or NaN for a hadron droplet. Skipping.";
      continue;
    }
    value_kernel /= norm;
    result += value_kernel * quantity_smear;
  }
  return result;
}

double HadronicLiquefier::get_source_rhob(const double tau, const double x,
                                          const double y,
                                          const double eta) const {
  QuantityType qtype = BARYON_NUMBER;
  return get_source_quantity(tau, x, y, eta, qtype);
}

double HadronicLiquefier::get_source_rhoq(const double tau, const double x,
                                          const double y,
                                          const double eta) const {
  QuantityType qtype = ELECTRIC_CHARGE;
  return get_source_quantity(tau, x, y, eta, qtype);
}

double HadronicLiquefier::get_source_rhos(const double tau, const double x,
                                          const double y,
                                          const double eta) const {
  QuantityType qtype = STRANGENESS;
  return get_source_quantity(tau, x, y, eta, qtype);
}

void HadronicLiquefier::get_sources_on_grid(const double tau,
                                            const std::vector<double> &x,
                                            const std::vector<double> &y,
                                            const std::vector<double> &eta,
                                            std::vector<int> &cells,
                                            std::vector<double> &vals) const {
  cells.clear();
  vals.clear();
  const int nx = static_cast<int>(x.size());
  const int ny = static_cast<int>(y.size());
  const int neta = static_cast<int>(eta.size());
  if (nx == 0 || ny == 0 || neta == 0) return;

  // The same box sizes as get_source_energy() and get_source_quantity().
  const double skip_dis_x = skip_n_sigma_transverse_ * sigma_transverse_;
  double skip_dis_eta;
  if (covariant_smearing_) {
    skip_dis_eta = skip_n_sigma_transverse_ * sigma_transverse_;
  } else {
    skip_dis_eta = skip_n_sigma_longitudinal_ * sigma_longitudinal_;
  }
  if (covariant_smearing_ && !hydro_Cartesian_) {
    skip_dis_eta = skip_dis_x / tau;
  }

  // The droplets the cell-by-cell functions use at this time, in their
  // order, with what does not depend on the cell computed once (the same
  // expressions, so the same values).
  struct Drop {
    double x, y, eta_s, mass, mT, rapidity, norm, ux, uy, uz, ueta, gamma;
    double px, py, E, pz;
    int B, Q, S;
  };
  std::vector<Drop> drops;
  // the step's droplets (prepare_active_droplets) or all, in list order
  const std::size_t n_loop =
      active_valid_ ? active_idx_.size() : hadron_droplets_list.size();
  for (std::size_t k_loop = 0; k_loop < n_loop; k_loop++) {
    const auto &drop_i =
        hadron_droplets_list[active_valid_ ? active_idx_[k_loop] : k_loop];
    auto xmu_i = drop_i.get_xmu();
    if (dtau_ > 0 &&
        ((xmu_i[0] < tau - 0.5 * dtau_) || (xmu_i[0] >= tau + 0.5 * dtau_))) {
      continue;
    }
    auto pmu_i = drop_i.get_pmu();
    Drop d;
    d.mass = sqrt(pmu_i[0] * pmu_i[0] - pmu_i[1] * pmu_i[1] -
                  pmu_i[2] * pmu_i[2] - pmu_i[3] * pmu_i[3]);
    if (d.mass <= 1e-16) continue;  // skip massless particles
    d.mT = sqrt(d.mass * d.mass + pmu_i[1] * pmu_i[1] + pmu_i[2] * pmu_i[2]);
    if (d.mT >= 1e-16) {
      d.rapidity = std::asinh(pmu_i[3] / d.mT);
    } else {
      d.rapidity = 0.5 * log((pmu_i[0] + pmu_i[3]) / (pmu_i[0] - pmu_i[3]));
    }
    d.norm = drop_i.get_normalization();
    if (d.norm <= 0.0 || std::isnan(d.norm)) {
      JSWARN << "Normalization is zero or NaN for a hadron droplet. Skipping.";
      continue;
    }
    d.x = xmu_i[1];
    d.y = xmu_i[2];
    d.eta_s = xmu_i[3];  // Cartesian z, Milne eta
    d.px = pmu_i[1];
    d.py = pmu_i[2];
    d.ux = pmu_i[1] / d.mass;
    d.uy = pmu_i[2] / d.mass;
    d.uz = pmu_i[3] / d.mass;
    const double peta = d.mT * sinh(d.rapidity - d.eta_s);
    d.ueta = peta / d.mass;
    d.gamma = hydro_Cartesian_ ? pmu_i[0] / d.mass   // Cartesian: lab frame
                               : d.mT * cosh(d.rapidity - d.eta_s) / d.mass;
    d.E = pmu_i[0];
    d.pz = pmu_i[3];
    d.B = drop_i.get_baryon_number();
    d.Q = drop_i.get_electric_charge();
    d.S = drop_i.get_strangeness();
    drops.push_back(d);
  }
  if (drops.empty()) return;

  // [first, last] index of the coordinates within [lo, hi], padded; the
  // exact tests of the cell-by-cell functions follow per cell
  const double pad =
      1e-9 * (std::max({std::abs(x.front()), std::abs(x.back()),
                        std::abs(y.front()), std::abs(y.back()), 1.0}));
  auto index_range = [](const std::vector<double> &v, double lo, double hi,
                        int &i0, int &i1) {
    i0 = static_cast<int>(std::lower_bound(v.begin(), v.end(), lo) - v.begin());
    i1 = static_cast<int>(std::upper_bound(v.begin(), v.end(), hi) -
                          v.begin()) - 1;
  };

  std::vector<std::vector<int>> slice_cells(neta);
  std::vector<std::vector<double>> slice_vals(neta);
  constexpr int NV = 7;

#pragma omp parallel
  {
    // this thread's (x, y) plane: the sums, and which cells have one
    std::vector<double> acc;
    std::vector<char> hit;
    std::vector<int> touched;
#pragma omp for schedule(dynamic)
    for (int ie = 0; ie < neta; ie++) {
      const double eta_c = eta[ie];
      touched.clear();
      for (const Drop &d : drops) {
        const double eta_diff = eta_c - d.eta_s;
        if (abs(eta_diff) > skip_dis_eta) continue;
        // per (droplet, slice): the longitudinal factors of jmu (Milne)
        const double ch = hydro_Cartesian_ ? 0. : cosh(d.rapidity - eta_c);
        const double sh = hydro_Cartesian_ ? 0. : sinh(d.rapidity - eta_c);
        int ix0, ix1, iy0, iy1;
        index_range(x, d.x - skip_dis_x - pad, d.x + skip_dis_x + pad, ix0,
                    ix1);
        index_range(y, d.y - skip_dis_x - pad, d.y + skip_dis_x + pad, iy0,
                    iy1);
        for (int ix = ix0; ix <= ix1; ix++) {
          const double x_diff = x[ix] - d.x;
          if (abs(x_diff) > skip_dis_x) continue;
          for (int iy = iy0; iy <= iy1; iy++) {
            const double y_diff = y[iy] - d.y;
            if (abs(y_diff) > skip_dis_x) continue;
            double value_kernel = 0.;
            if (covariant_smearing_ && !hydro_Cartesian_) {
              value_kernel =
                  smearing_kernel_covariant_Milne(x_diff, y_diff, eta_diff,
                                                  d.ux, d.uy, d.ueta, tau,
                                                  d.gamma) /
                  tau;
            } else if (!covariant_smearing_ && !hydro_Cartesian_) {
              value_kernel =
                  smearing_kernel_gaussian(x_diff, y_diff, eta_diff) / tau;
            } else if (covariant_smearing_ && hydro_Cartesian_) {
              value_kernel = smearing_kernel_covariant_Cartesian(
                  x_diff, y_diff, eta_diff, d.ux, d.uy, d.uz, d.gamma);
            } else {
              value_kernel = smearing_kernel_gaussian(x_diff, y_diff, eta_diff);
            }
            value_kernel /= d.norm;

            if (acc.empty()) {
              acc.assign(static_cast<std::size_t>(NV) * nx * ny, 0.0);
              hit.assign(static_cast<std::size_t>(nx) * ny, 0);
            }
            const int c = ix * ny + iy;
            if (!hit[c]) {
              hit[c] = 1;
              touched.push_back(c);
            }
            double *a = &acc[static_cast<std::size_t>(NV) * c];
            if (hydro_Cartesian_) {   // as get_source_energy
              a[0] += value_kernel * d.E;
              a[1] += value_kernel * d.px;
              a[2] += value_kernel * d.py;
              a[3] += value_kernel * d.pz;
            } else {
              a[0] += value_kernel * d.mT * ch;
              a[1] += value_kernel * d.px;
              a[2] += value_kernel * d.py;
              a[3] += value_kernel * d.mT * sh;
            }
            // get_source_quantity() skips the droplets without the charge
            if (d.B != 0) a[4] += value_kernel * d.B;
            if (d.Q != 0) a[5] += value_kernel * d.Q;
            if (d.S != 0) a[6] += value_kernel * d.S;
          }
        }
      }
      // hand this slice's cells over and reset the plane
      std::sort(touched.begin(), touched.end());
      auto &sc = slice_cells[ie];
      auto &sv = slice_vals[ie];
      sc.reserve(3 * touched.size());
      sv.reserve(NV * touched.size());
      for (const int c : touched) {
        sc.push_back(c / ny);
        sc.push_back(c % ny);
        sc.push_back(ie);
        double *a = &acc[static_cast<std::size_t>(NV) * c];
        for (int i = 0; i < NV; i++) {
          sv.push_back(a[i]);
          a[i] = 0.0;
        }
        hit[c] = 0;
      }
    }
  }
  for (int ie = 0; ie < neta; ie++) {
    cells.insert(cells.end(), slice_cells[ie].begin(), slice_cells[ie].end());
    vals.insert(vals.end(), slice_vals[ie].begin(), slice_vals[ie].end());
  }
}

bool HadronicLiquefier::gamma_factor_too_large_check(
    const Hadron &hadron) const {
  auto p_init = hadron.p_in();
  // compute the hadron mass and rapidity
  std::array<double, 4> pmu_i = {
      static_cast<double>(p_init.t()), static_cast<double>(p_init.x()),
      static_cast<double>(p_init.y()), static_cast<double>(p_init.z())};

  const double mass = sqrt(pmu_i[0] * pmu_i[0] - pmu_i[1] * pmu_i[1] -
                           pmu_i[2] * pmu_i[2] - pmu_i[3] * pmu_i[3]);
  if (mass <= 1e-16) {
    return true;  // skip massless particles
  }
  auto x_init = hadron.x_in();
  const double eta_s =
      0.5 * log((x_init.t() + x_init.z()) / (x_init.t() - x_init.z()));
  const double mT =
      sqrt(mass * mass + pmu_i[1] * pmu_i[1] + pmu_i[2] * pmu_i[2]);
  double rapidity;
  if (mT >= 1e-16) {
    rapidity = std::asinh(pmu_i[3] / mT);
  } else {
    rapidity = 0.5 * log((pmu_i[0] + pmu_i[3]) / (pmu_i[0] - pmu_i[3]));
  }

  double gamma = 1.0;
  if (!hydro_Cartesian_) {
    const double ux = pmu_i[1] / mass;
    const double uy = pmu_i[2] / mass;
    const double peta = mT * sinh(rapidity - eta_s);
    const double ueta = peta / mass;
    gamma = mT * cosh(rapidity - eta_s) / mass;
  } else {
    const double ux = pmu_i[1] / mass;
    const double uy = pmu_i[2] / mass;
    const double uz = pmu_i[3] / mass;
    gamma = pmu_i[0] / mass;  // Cartesian: the lab frame
  }

  // Check if the gamma factor is larger than the maximum allowed value
  if (gamma > gamma_factor_max_) {
    VERBOSE(3) << "HadronicLiquefier: Gamma factor too large for hadron with "
               << "mass " << mass << ", rapidity " << rapidity << ", eta_s "
               << eta_s << ", and gamma " << gamma
               << ". Not added to hydro sources. Maximum allowed gamma = "
               << gamma_factor_max_;
    return true;  // Too large gamma factor, do not add to hydro sources
  }
  return false;  // Gamma factor is acceptable
}

std::vector<shared_ptr<Hadron>> HadronicLiquefier::add_hydro_sources_hadrons(
    std::vector<Hadron> &hIn) {
  // Create vector to store hadrons not added to the hydro sources
  std::vector<shared_ptr<Hadron>> hadrons_not_added;

  // Create droplets from the hadrons
  for (const auto &hadron : hIn) {
    // Check if the gamma factor is too large for covariant smearing
    if (covariant_smearing_ && gamma_factor_too_large_check(hadron)) {
      hadrons_not_added.push_back(make_shared<Hadron>(hadron));
      continue;
    }
    auto x_init = hadron.x_in();
    auto p_init = hadron.p_in();
    std::array<double, 4> x_hadron = {0.0, 0.0, 0.0, 0.0};
    if (hydro_Cartesian_) {
      x_hadron = {
          static_cast<double>(x_init.t()), static_cast<double>(x_init.x()),
          static_cast<double>(x_init.y()), static_cast<double>(x_init.z())};
    } else {
      const double tau =
          sqrt(x_init.t() * x_init.t() - x_init.z() * x_init.z());
      const double eta_s =
          0.5 * log((x_init.t() + x_init.z()) / (x_init.t() - x_init.z()));
      x_hadron = {static_cast<double>(tau), static_cast<double>(x_init.x()),
                  static_cast<double>(x_init.y()), static_cast<double>(eta_s)};
    }
    std::array<double, 4> p_hadron = {
        static_cast<double>(p_init.t()), static_cast<double>(p_init.x()),
        static_cast<double>(p_init.y()), static_cast<double>(p_init.z())};

    HadronDroplet hadron_droplet =
        HadronDroplet(x_hadron, p_hadron, hadron.baryon_number(),
                      hadron.charge(), hadron.strangeness());

    double norm =
        compute_drop_kernel_normalization(x_hadron[0], hadron_droplet);
    hadron_droplet.set_normalization(norm);
    active_valid_ = false;
    hadron_droplets_list.push_back(hadron_droplet);
  }
  return hadrons_not_added;
}

Jetscape::real HadronicLiquefier::get_dropletlist_total_energy() const {
  Jetscape::real total_E = 0.0;
  for (const auto &drop_i : hadron_droplets_list) {
    total_E += drop_i.get_pmu()[0];
  }
  return (total_E);
}

Jetscape::real HadronicLiquefier::get_dropletlist_net_baryon_number() const {
  Jetscape::real net_baryon_number = 0.0;
  for (const auto &drop_i : hadron_droplets_list) {
    net_baryon_number += drop_i.get_baryon_number();
  }
  return (net_baryon_number);
}

Jetscape::real HadronicLiquefier::get_dropletlist_net_electric_charge() const {
  Jetscape::real net_electric_charge = 0.0;
  for (const auto &drop_i : hadron_droplets_list) {
    net_electric_charge += drop_i.get_electric_charge();
  }
  return (net_electric_charge);
}

Jetscape::real HadronicLiquefier::get_dropletlist_net_strangeness() const {
  Jetscape::real net_strangeness = 0.0;
  for (const auto &drop_i : hadron_droplets_list) {
    net_strangeness += drop_i.get_strangeness();
  }
  return (net_strangeness);
}

void HadronicLiquefier::ClearTask() {
  hadron_droplets_list.clear();
  active_valid_ = false;
}

};  // namespace Jetscape