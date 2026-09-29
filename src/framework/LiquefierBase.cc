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
#include "LiquefierBase.h"
#include "JetScapeXML.h"
#include <math.h>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace Jetscape {

LiquefierBase::LiquefierBase()
    : hydro_source_abs_err(1e-10),
      e_mom_rel_tol(5e-2),
      drop_stat(-11),
      miss_stat(-13),
      neg_stat(-17) {
  GetHydroCellSignalConnected = false;
}

// Whether (tau, x, y, eta) is in the future light cone of the droplet; the
// same check get_source() makes before calling smearing_kernel().
bool LiquefierBase::in_light_cone(Jetscape::real tau, Jetscape::real x,
                                  Jetscape::real y, Jetscape::real eta,
                                  const Droplet &drop_i) const {
  const auto x_drop = drop_i.get_xmu();
  double ds2 = tau * tau + x_drop[0] * x_drop[0] -
               2.0 * tau * x_drop[0] * cosh(eta - x_drop[3]) -
               (x - x_drop[1]) * (x - x_drop[1]) -
               (y - x_drop[2]) * (y - x_drop[2]);
  return tau >= x_drop[0] && ds2 >= 0.0;
}

void LiquefierBase::get_source(Jetscape::real tau, Jetscape::real x,
                               Jetscape::real y, Jetscape::real eta,
                               std::array<Jetscape::real, 4> &jmu) const {
  jmu = {0.0, 0.0, 0.0, 0.0};
  // Inside a prepared window only the droplets that can contribute there are
  // visited, in their original order; the others would add exactly zero.
  const bool pruned = (tau >= active_tau_lo_ && tau <= active_tau_hi_);
  const std::size_t n_loop =
      pruned ? active_droplets_.size() : dropletlist.size();
  for (std::size_t k = 0; k < n_loop; k++) {
    const std::size_t idx = pruned ? active_droplets_[k] : k;
    const auto &drop_i = dropletlist[idx];
    const double norm = droplet_norm_[idx];
    if (norm == 0.0)
      continue;  // deposited as a point (below) or lost off the grid
    if (in_light_cone(tau, x, y, eta, drop_i)) {
      std::array<Jetscape::real, 4> jmu_i = {0.0, 0.0, 0.0, 0.0};
      smearing_kernel(tau, x, y, eta, drop_i, jmu_i);
      if (norm == 1.0) {  // not normalized: bit-identical to before
        for (int i = 0; i < 4; i++)
          jmu[i] += jmu_i[i];
      } else {
        for (int i = 0; i < 4; i++)
          jmu[i] += static_cast<Jetscape::real>(norm * jmu_i[i]);
      }
    }
  }
  if (point_deposits_.empty())
    return;
  const auto &g = hydro_grid_;
  for (const auto &pd : point_deposits_) {
    if (std::abs(tau - pd.tau_q) < 0.5 * pd.dtau &&
        std::abs(x - pd.x) < 0.5 * g.dx && std::abs(y - pd.y) < 0.5 * g.dy &&
        std::abs(eta - pd.eta) < 0.5 * g.deta) {
      // p in tau-eta components, per unit tau dtau dx dy deta: the hydro adds
      // tau dtau dV J^mu at this query time, i.e. exactly the droplet's p^mu
      const auto p = dropletlist[pd.idx].get_pmu();
      const double ch = cosh(pd.eta), sh = sinh(pd.eta);
      const double w = 1.0 / (pd.tau_q * pd.dtau * g.dx * g.dy * g.deta);
      jmu[0] += static_cast<Jetscape::real>(w * (p[0] * ch - p[3] * sh));
      jmu[1] += static_cast<Jetscape::real>(w * p[1]);
      jmu[2] += static_cast<Jetscape::real>(w * p[2]);
      jmu[3] += static_cast<Jetscape::real>(w * (p[3] * ch - p[0] * sh));
    }
  }
}

//! If vertex does not conserve energy and momentum,
//! a p_missing will be added to the pOut list
void LiquefierBase::check_energy_momentum_conservation(
    const std::vector<Parton> &pIn, std::vector<Parton> &pOut) {
  FourVector p_init(0., 0., 0., 0.);
  for (const auto &iparton : pIn) {
    auto temp = iparton.p_in();
    if (iparton.pstat() == -1) {
      p_init -= temp;
    } else {
      p_init += temp;
    }
  }

  FourVector p_final(0., 0., 0., 0.);
  FourVector x_final(0., 0., 0., 0.);
  for (const auto &iparton : pOut) {
    auto temp = iparton.p_in();
    if (iparton.pstat() == -1) {
      p_final -= temp;
    } else {
      p_final += temp;
    }
    x_final = iparton.x_in();
  }

  FourVector p_missing = p_init;
  p_missing -= p_final;
  if (std::abs(p_missing.t()) > hydro_source_abs_err ||
      std::abs(p_missing.x()) > hydro_source_abs_err ||
      std::abs(p_missing.y()) > hydro_source_abs_err ||
      std::abs(p_missing.z()) > hydro_source_abs_err) {
    // Energy-loss kinematics (on-shell E recomputed from p, soft partons cut)
    // miss by far more than 1e-10 GeV, so warn only above e_mom_rel_tol.
    const double dmax = std::max(
        {std::abs(p_missing.t()), std::abs(p_missing.x()),
         std::abs(p_missing.y()), std::abs(p_missing.z())});
    if (dmax > e_mom_rel_tol * std::abs(p_init.t())) {
      JSWARN << "A vertex does not conserve energy momentum!";
      JSWARN << "E = " << p_missing.t() << " GeV, px = " << p_missing.x()
             << " GeV, py = " << p_missing.y() << " GeV, pz = "
             << p_missing.z() << " GeV (E_in = " << p_init.t() << " GeV).";
    }
    Parton parton_miss(0, 21, miss_stat, p_missing, x_final);
    pOut.push_back(parton_miss);
  }
}

void LiquefierBase::filter_partons(std::vector<Parton> &pOut) {
  // threshold_energy_switch = 1, use e_threshold
  // threshold_energy_switch = 0, use |e_threshold|*T
  threshold_energy_switch = JetScapeXML::Instance()->GetElementInt(
      {"Liquefier", "threshold_energy_switch"});
  if (threshold_energy_switch != 0 && threshold_energy_switch != 1) {
    JSWARN << "threshold_energy_switch should be 0 or 1, but it is "
           << threshold_energy_switch;
    exit(1);
  }
  e_threshold = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "e_threshold"});  // GeV

  for (auto &iparton : pOut) {
    if (iparton.pstat() == miss_stat)
      continue;

    // ignore photons
    if (iparton.isPhoton(iparton.pid()))
      continue;

    // ignore heavy quarks
    if (std::abs(iparton.pid()) == 4 || std::abs(iparton.pid()) == 5)
      continue;

    if (iparton.pstat() == -1) {
      // remove negative particles from parton list
      // iparton.set_stat(drop_stat);
      iparton.set_stat(neg_stat);
      continue;
    }

    // for positive particles, including jet partons and recoil partons
    auto tLoc = iparton.x_in().t();
    auto xLoc = iparton.x_in().x();
    auto yLoc = iparton.x_in().y();
    auto zLoc = iparton.x_in().z();
    std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
    GetHydroCellSignal(tLoc, xLoc, yLoc, zLoc, check_fluid_info_ptr);
    auto vxLoc = check_fluid_info_ptr->vx;
    auto vyLoc = check_fluid_info_ptr->vy;
    auto vzLoc = check_fluid_info_ptr->vz;
    auto beta2 = vxLoc * vxLoc + vyLoc * vyLoc + vzLoc * vzLoc;
    auto gamma = 1.0 / sqrt(1.0 - beta2);
    auto E_boosted = gamma * (iparton.e() - iparton.p(1) * vxLoc -
                              iparton.p(2) * vyLoc - iparton.p(3) * vzLoc);
    if (threshold_energy_switch) {
      // drop partons with energy smaller than e_threshold
      // (in the local rest frame) from parton list
      if (E_boosted < e_threshold) {
        iparton.set_stat(drop_stat);
        continue;
      }
    } else {
      // drop partons with energy smaller than |e_threshold|*T
      // (in the local rest frame) from parton list
      auto tempLoc = check_fluid_info_ptr->temperature;
      if (E_boosted < std::abs(e_threshold) * tempLoc) {
        iparton.set_stat(drop_stat);
        continue;
      }
    }
  }
}

void LiquefierBase::add_hydro_sources(std::vector<Parton> &pIn,
                                      std::vector<Parton> &pOut) {
  if (pOut.size() == 0) {
    // the process is freestreaming, ignore
    filter_partons(pIn);
    return;
  }
  check_energy_momentum_conservation(pIn, pOut);
  filter_partons(pOut);

  FourVector p_final;
  FourVector p_init;
  FourVector x_final;
  FourVector x_init;
  // cout << "debug, mid ......." << pOut.size() << endl;
  //  use energy conservation to deterime the source term
  const auto weight_init = 1.0;
  for (const auto &iparton : pIn) {
    auto temp = iparton.p_in();
    p_init += temp;
    x_init = iparton.x_in();
  }

  auto weight_final = 0.0;
  for (const auto &iparton : pOut) {
    if (iparton.pstat() == drop_stat)
      continue;
    if (iparton.pstat() == miss_stat)
      continue;
    if (iparton.pstat() == neg_stat)
      continue;
    auto temp = iparton.p_in();
    p_final += temp;
    x_final = iparton.x_in();
    weight_final = 1.0;
  }

  if (std::abs(p_init.t() - p_final.t()) / p_init.t() > hydro_source_abs_err ||
      std::abs(p_init.x() - p_final.x()) / p_init.t() > hydro_source_abs_err ||
      std::abs(p_init.y() - p_final.y()) / p_init.t() > hydro_source_abs_err ||
      std::abs(p_init.z() - p_final.z()) / p_init.t() > hydro_source_abs_err) {
    auto droplet_t = ((x_final.t() * weight_final + x_init.t() * weight_init) /
                      (weight_final + weight_init));
    auto droplet_x = ((x_final.x() * weight_final + x_init.x() * weight_init) /
                      (weight_final + weight_init));
    auto droplet_y = ((x_final.y() * weight_final + x_init.y() * weight_init) /
                      (weight_final + weight_init));
    auto droplet_z = ((x_final.z() * weight_final + x_init.z() * weight_init) /
                      (weight_final + weight_init));

    auto droplet_tau = sqrt(droplet_t * droplet_t - droplet_z * droplet_z);
    auto droplet_eta =
        (0.5 * log((droplet_t + droplet_z) / (droplet_t - droplet_z)));
    // A droplet at the hard vertex (t = z = 0) has eta = 0/0: NaN, which makes
    // every kernel value NaN, so the droplet never deposits.  At tau = 0 any
    // finite eta is the same point; use 0.
    if (!(droplet_tau > 0.0)) {
      droplet_tau = 0.0;
      droplet_eta = 0.0;
    }
    auto droplet_E = p_init.t() - p_final.t();
    auto droplet_px = p_init.x() - p_final.x();
    auto droplet_py = p_init.y() - p_final.y();
    auto droplet_pz = p_init.z() - p_final.z();

    std::array<Jetscape::real, 4> droplet_xmu = {
        static_cast<Jetscape::real>(droplet_tau),
        static_cast<Jetscape::real>(droplet_x),
        static_cast<Jetscape::real>(droplet_y),
        static_cast<Jetscape::real>(droplet_eta)};
    std::array<Jetscape::real, 4> droplet_pmu = {
        static_cast<Jetscape::real>(droplet_E),
        static_cast<Jetscape::real>(droplet_px),
        static_cast<Jetscape::real>(droplet_py),
        static_cast<Jetscape::real>(droplet_pz)};
    Droplet drop_i(droplet_xmu, droplet_pmu);
    add_a_droplet(drop_i);
  }
}

void LiquefierBase::ClearTask() {
  dropletlist.clear();
  droplet_norm_.clear();
  droplet_flux_.clear();
  reset_normalization();
  invalidate_active_droplets();
}

// ── normalization on the hydro grid ───────────────────────────────────────
void LiquefierBase::reset_normalization() {
  droplet_norm_.assign(dropletlist.size(), 1.0);
  droplet_flux_.assign(dropletlist.size(), -1.0);
  point_deposits_.clear();
  n_normalized_ = n_point_ = n_lost_ = 0;
  E_droplets_ = E_sampled_ = 0.;
  flux_min_ = 1e300;
  flux_max_ = 0.;
}

void LiquefierBase::set_hydro_grid(const HydroGrid &grid) {
  if (grid == hydro_grid_)
    return;
  hydro_grid_ = grid;
  reset_normalization();
}

double LiquefierBase::sampled_flux(const Droplet &drop_i, double tau,
                                   double dtau, double *tau_q_out) const {
  const auto &g = hydro_grid_;
  const double t_dep = deposit_time(drop_i);
  const auto xd = drop_i.get_xmu();
  // unit energy at rest in the lab: smearing_kernel() then returns
  // (n.j) (cosh eta, 0, 0, -sinh eta), whose lab energy is n.j
  const Droplet unit(xd, {1.0, 0.0, 0.0, 0.0});
  const double dV = g.dx * g.dy * g.deta;
  // The hydro's query times near here are tau + k dtau (both Runge-Kutta
  // stages; each gets weight dtau in total).  The kernel deposits at those
  // within deposit_half_width() of the deposit time: all of them are summed
  // (smearing_kernel() makes the exact cut).  tau_q_out is the nearest one.
  const double hw = deposit_half_width() + 1e-6;
  const double reach = std::max(hw, 0.5 * dtau + 1e-6);
  double flux = 0.0, best = -1.0, best_d = 1e300;
  for (int k = -3; k <= 6; k++) {
    const double tq = tau + k * dtau;
    if (tq <= 0.0)
      continue;
    const double d = std::abs(tq - t_dep);
    if (d <= reach && d < best_d) {
      best_d = d;
      best = tq;
    }
    if (d > hw)
      continue;
    const Jetscape::real tq_r = static_cast<Jetscape::real>(tq);
    double s = 0.0;
#pragma omp parallel for reduction(+ : s) schedule(dynamic)
    for (int ie = 0; ie < g.neta; ie++) {
      const Jetscape::real eta_r =
          static_cast<Jetscape::real>(g.eta_min + ie * g.deta);
      // transverse radius of the droplet's future light cone in this slice
      const double R2 = tq_r * tq_r + xd[0] * xd[0] -
                        2.0 * tq_r * xd[0] * cosh(eta_r - xd[3]);
      if (!(R2 >= 0.0))
        continue;
      const double R = std::sqrt(R2);
      const int ix0 = std::max(0, (int)std::floor((xd[1] - R - g.x_min) / g.dx) - 1);
      const int ix1 = std::min(g.nx - 1, (int)std::ceil((xd[1] + R - g.x_min) / g.dx) + 1);
      const int iy0 = std::max(0, (int)std::floor((xd[2] - R - g.y_min) / g.dy) - 1);
      const int iy1 = std::min(g.ny - 1, (int)std::ceil((xd[2] + R - g.y_min) / g.dy) + 1);
      const double ch = cosh(eta_r), sh = sinh(eta_r);
      for (int ix = ix0; ix <= ix1; ix++) {
        const Jetscape::real x = static_cast<Jetscape::real>(g.x_min + ix * g.dx);
        for (int iy = iy0; iy <= iy1; iy++) {
          const Jetscape::real y = static_cast<Jetscape::real>(g.y_min + iy * g.dy);
          if (!in_light_cone(tq_r, x, y, eta_r, drop_i))
            continue;
          std::array<Jetscape::real, 4> j = {0.0, 0.0, 0.0, 0.0};
          smearing_kernel(tq_r, x, y, eta_r, unit, j);
          s += ch * j[0] + sh * j[3];
        }
      }
    }
    flux += s * tq * dtau * dV;
  }
  if (tau_q_out)
    *tau_q_out = best;
  return flux;
}

void LiquefierBase::normalize_active_droplets(double tau, double dtau) {
  if (!normalize_on_grid_ || !hydro_grid_.valid() || !normalizable_on_grid() ||
      dtau <= 0.0)
    return;
  const auto &g = hydro_grid_;
  auto visit = [&](std::size_t idx) {
    if (droplet_flux_[idx] >= 0.0)
      return;  // already normalized
    const Droplet &d = dropletlist[idx];
    double tq = -1.0;
    const double flux = sampled_flux(d, tau, dtau, &tq);
    if (tq < 0.0)
      return;  // its query time is not within reach of this step yet
    droplet_flux_[idx] = flux;
    const auto p = d.get_pmu();
    if (flux > 1e-6) {
      droplet_norm_[idx] = 1.0 / flux;
      n_normalized_++;
      E_droplets_ += p[0];
      E_sampled_ += flux * p[0];
      flux_min_ = std::min(flux_min_, flux);
      flux_max_ = std::max(flux_max_, flux);
      return;
    }
    // the kernel misses every cell centre: deposit it whole into the nearest
    // cell at its rapidity (a droplet at tau_d = 0 has eta = 0/0: use 0)
    droplet_norm_[idx] = 0.0;
    const auto xd = d.get_xmu();
    const double eta_d = std::isfinite(xd[3]) ? xd[3] : 0.0;
    const long ix = std::lround((xd[1] - g.x_min) / g.dx);
    const long iy = std::lround((xd[2] - g.y_min) / g.dy);
    const long ie = std::lround((eta_d - g.eta_min) / g.deta);
    if (ix < 0 || ix >= g.nx || iy < 0 || iy >= g.ny || ie < 0 || ie >= g.neta) {
      n_lost_++;
      return;
    }
    point_deposits_.push_back({static_cast<int>(idx), tq, dtau,
                               g.x_min + ix * g.dx, g.y_min + iy * g.dy,
                               g.eta_min + ie * g.deta});
    n_point_++;
  };
  if (active_droplets_prepared()) {
    for (int idx : active_droplets_)
      visit(static_cast<std::size_t>(idx));
  } else {
    for (std::size_t idx = 0; idx < dropletlist.size(); idx++)
      visit(idx);
  }
}

std::string LiquefierBase::normalization_summary() const {
  if (!normalize_on_grid_)
    return "normalization on the hydro grid: off";
  if (!hydro_grid_.valid())
    return "normalization on the hydro grid: no grid set (not normalized)";
  char buf[400];
  snprintf(buf, sizeof(buf),
           "normalized %d droplets on the hydro grid (sampled flux %.3g .. "
           "%.3g; their sampled energy %.4g GeV -> %.4g GeV), %d deposited "
           "into one cell (flux 0), %d lost off the grid",
           n_normalized_, n_normalized_ ? flux_min_ : 0.0,
           n_normalized_ ? flux_max_ : 0.0, E_sampled_, E_droplets_, n_point_,
           n_lost_);
  return buf;
}

void LiquefierBase::prepare_active_droplets(double tau_lo, double tau_hi) {
  active_droplets_.clear();
  for (std::size_t i = 0; i < dropletlist.size(); i++) {
    if (droplet_may_contribute(dropletlist[i], tau_lo, tau_hi)) {
      active_droplets_.push_back(static_cast<int>(i));
    }
  }
  active_tau_lo_ = tau_lo;
  active_tau_hi_ = tau_hi;
}

Jetscape::real LiquefierBase::get_dropletlist_total_energy() const {
  Jetscape::real total_E = 0.0;
  for (const auto &drop_i : dropletlist) {
    total_E += drop_i.get_pmu()[0];
  }
  return (total_E);
}

};  // namespace Jetscape
