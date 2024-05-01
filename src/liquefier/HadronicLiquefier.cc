/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

namespace Jetscape {
HadronicLiquefier::HadronicLiquefier() {
  InitializeParameters();
}

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
  dx_ = 2. * xMax_ / (Nx_-1.);
  dy_ = 2. * yMax_ / (Ny_-1.);
  dz_ = 2. * zMax_ / (Nz_-1.);
  hydro_Cartesian_ = hydro_Cartesian;
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

  // get the grid specifications from the initial state module, which are
  // also used in the hydro evolution
  Nx_ = ini->GetXSize();
  Ny_ = ini->GetYSize();
  Nz_ = ini->GetZSize();
  dx_ = ini->GetXStep();
  dy_ = ini->GetYStep();
  dz_ = ini->GetZStep();
  xMax_ = ini->GetXMax();
  yMax_ = ini->GetYMax();
  zMax_ = ini->GetZMax();
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

  for (int ieta = 0; ieta < Nz_; ieta++) {
    for (int ix = 0; ix < Nx_; ix++) {
      for (int iy = 0; iy < Ny_; iy++) {
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
        const double eta_s =
            0.5 * log((xmu_i[0] + xmu_i[3]) / (xmu_i[0] - xmu_i[3]));
        double eta_diff;
        if (!hydro_Cartesian_) {
          eta_diff = eta - eta_s;
        } else {
          eta_diff = eta - xmu_i[3];
        }
        if (abs(eta_diff) > skip_dis_eta) {
          continue;
        }

        const double mass = sqrt(pmu_i[0] * pmu_i[0] - pmu_i[1] * pmu_i[1] -
                                 pmu_i[2] * pmu_i[2] - pmu_i[3] * pmu_i[3]);
        const double rapidity =
            0.5 * log((pmu_i[0] + pmu_i[3]) / (pmu_i[0] - pmu_i[3]));
        const double mT =
            sqrt(mass * mass + pmu_i[1] * pmu_i[1] + pmu_i[2] * pmu_i[2]);

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
          const double gamma = mT * cosh(rapidity - eta_s) / mass;

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

void HadronicLiquefier::get_source_energy(
    const double tau, const double x, const double y, const double eta,
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
  for (const auto &drop_i : hadron_droplets_list) {
    auto xmu_i = drop_i.get_xmu();
    auto pmu_i = drop_i.get_pmu();

    double x_diff = x - xmu_i[1];
    if (abs(x_diff) > skip_dis_x) {
      continue;
    }
    double y_diff = y - xmu_i[2];
    if (abs(y_diff) > skip_dis_x) {
      continue;
    }
    double eta_diff;
    const double eta_s =
        0.5 * log((xmu_i[0] + xmu_i[3]) / (xmu_i[0] - xmu_i[3]));
    if (!hydro_Cartesian_) {
      eta_diff = eta - eta_s;
    } else {
      eta_diff = eta - xmu_i[3];
    }
    if (abs(eta_diff) > skip_dis_eta) {
      continue;
    }

    const double mass = sqrt(pmu_i[0] * pmu_i[0] - pmu_i[1] * pmu_i[1] -
                             pmu_i[2] * pmu_i[2] - pmu_i[3] * pmu_i[3]);
    const double rapidity =
        0.5 * log((pmu_i[0] + pmu_i[3]) / (pmu_i[0] - pmu_i[3]));
    const double mT =
        sqrt(mass * mass + pmu_i[1] * pmu_i[1] + pmu_i[2] * pmu_i[2]);

    if (covariant_smearing_ && !hydro_Cartesian_) {
      const double ux = pmu_i[1] / mass;
      const double uy = pmu_i[2] / mass;
      const double peta = mT * sinh(rapidity - eta_s);
      const double ueta = peta / mass;
      const double gamma = mT * cosh(rapidity - eta_s) / mass;

      value_kernel = smearing_kernel_covariant_Milne(x_diff, y_diff, eta_diff,
                                                     ux, uy, ueta, tau, gamma);
    } else if (!covariant_smearing_ && !hydro_Cartesian_) {
      value_kernel = smearing_kernel_gaussian(x_diff, y_diff, eta_diff);
    } else if (covariant_smearing_ && hydro_Cartesian_) {
      const double ux = pmu_i[1] / mass;
      const double uy = pmu_i[2] / mass;
      const double uz = pmu_i[3] / mass;
      const double gamma = mT * cosh(rapidity - eta_s) / mass;

      value_kernel = smearing_kernel_covariant_Cartesian(
          x_diff, y_diff, eta_diff, ux, uy, uz, gamma);
    } else if (!covariant_smearing_ && hydro_Cartesian_) {
      value_kernel = smearing_kernel_gaussian(x_diff, y_diff, eta_diff);
    }
    value_kernel /= drop_i.get_normalization();

    jmu[0] += value_kernel * mT * cosh(rapidity - eta);
    jmu[1] += value_kernel * pmu_i[1];
    jmu[2] += value_kernel * pmu_i[2];
    jmu[3] += value_kernel * mT * sinh(rapidity - eta);
  }
}

double HadronicLiquefier::get_source_quantity(const double tau,
                                              const double x,
                                              const double y,
                                              const double eta,
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
  for (const auto &drop_i : hadron_droplets_list) {
    int quantity_smear = 0;
    if (quantity_smear == BARYON_NUMBER) {
      quantity_smear = drop_i.get_baryon_number();
    } else if (quantity_smear == ELECTRIC_CHARGE) {
      quantity_smear = drop_i.get_electric_charge();
    } else if (quantity_smear == STRANGENESS) {
      quantity_smear = drop_i.get_strangeness();
    } else {
      JSWARN << "The quantity to smear is not implemented.";
      exit(1);
    }
    if (quantity_smear == 0) {
      continue;
    }

    auto xmu_i = drop_i.get_xmu();
    auto pmu_i = drop_i.get_pmu();

    double x_diff = x - xmu_i[1];
    if (abs(x_diff) > skip_dis_x) {
      continue;
    }
    double y_diff = y - xmu_i[2];
    if (abs(y_diff) > skip_dis_x) {
      continue;
    }
    double eta_diff;
    const double eta_s =
        0.5 * log((xmu_i[0] + xmu_i[3]) / (xmu_i[0] - xmu_i[3]));
    if (hydro_Cartesian_) {
      eta_diff = eta - xmu_i[3];
    } else {
      eta_diff = eta - eta_s;
    }
    if (abs(eta_diff) > skip_dis_eta) {
      continue;
    }

    const double mass = sqrt(pmu_i[0] * pmu_i[0] - pmu_i[1] * pmu_i[1] -
                             pmu_i[2] * pmu_i[2] - pmu_i[3] * pmu_i[3]);
    const double rapidity =
        0.5 * log((pmu_i[0] + pmu_i[3]) / (pmu_i[0] - pmu_i[3]));
    const double mT =
        sqrt(mass * mass + pmu_i[1] * pmu_i[1] + pmu_i[2] * pmu_i[2]);

    if (covariant_smearing_ && !hydro_Cartesian_) {
      const double ux = pmu_i[1] / mass;
      const double uy = pmu_i[2] / mass;
      const double peta = mT * sinh(rapidity - eta_s);
      const double ueta = peta / mass;
      const double gamma = mT * cosh(rapidity - eta_s) / mass;

      value_kernel = smearing_kernel_covariant_Milne(x_diff, y_diff, eta_diff,
                                                     ux, uy, ueta, tau, gamma);
    } else if (!covariant_smearing_ && !hydro_Cartesian_) {
      value_kernel = smearing_kernel_gaussian(x_diff, y_diff, eta_diff);
    } else if (covariant_smearing_ && hydro_Cartesian_) {
      const double ux = pmu_i[1] / mass;
      const double uy = pmu_i[2] / mass;
      const double uz = pmu_i[3] / mass;
      const double gamma = mT * cosh(rapidity - eta_s) / mass;

      value_kernel = smearing_kernel_covariant_Cartesian(
          x_diff, y_diff, eta_diff, ux, uy, uz, gamma);
    } else if (!covariant_smearing_ && hydro_Cartesian_) {
      value_kernel = smearing_kernel_gaussian(x_diff, y_diff, eta_diff);
    }
    value_kernel /= drop_i.get_normalization();
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

void HadronicLiquefier::add_hydro_sources_hadrons(const double tau,
                                                  std::vector<Hadron> &hIn) {
  // Create droplets from the hadrons
  for (const auto &hadron : hIn) {
    auto p_init = hadron.p_in();
    auto x_init = hadron.x_in();

    std::array<double, 4> x_hadron = {
        static_cast<double>(x_init.t()),
        static_cast<double>(x_init.x()),
        static_cast<double>(x_init.y()),
        static_cast<double>(x_init.z())};
    std::array<double, 4> p_hadron = {
        static_cast<double>(p_init.t()),
        static_cast<double>(p_init.x()),
        static_cast<double>(p_init.y()),
        static_cast<double>(p_init.z())};

    int baryon_number = hadron.baryon_number();
    int electric_charge = hadron.charge();
    int strangeness = hadron.strangeness();
    HadronDroplet hadron_droplet = HadronDroplet(x_hadron, p_hadron, 
                          baryon_number, electric_charge, strangeness);

    double norm = compute_drop_kernel_normalization(tau, hadron_droplet);
    hadron_droplet.set_normalization(norm);
    hadron_droplets_list.push_back(hadron_droplet);
  }
}

void HadronicLiquefier::ClearTask() { hadron_droplets_list.clear(); }

}; // namespace Jetscape