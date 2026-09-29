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

#include "CausalLiquefier.h"
#include "LiquefierBase.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

using namespace Jetscape;

// check coordinate transformation
TEST(CausalLiquifierTest, TEST_COORDINATES) {
  CausalLiquefier lqf(0.3, 0.3, 0.3, 0.3);

  // for the transformation from tau-eta to t-z (configuration space)
  EXPECT_DOUBLE_EQ(0.0, lqf.get_t(0.0, 0.5));
  EXPECT_DOUBLE_EQ(5.0, lqf.get_t(5.0, 0.0));
  EXPECT_DOUBLE_EQ(0.0, lqf.get_z(0.0, 0.5));
  EXPECT_DOUBLE_EQ(0.0, lqf.get_z(5.0, 0.0));

  for (double tau = 0.1; tau < 0.3; tau += 0.1) {
    for (double eta = 0.0; eta < 0.2; eta += 0.1) {
      double t = lqf.get_t(tau, eta);
      double z = lqf.get_z(tau, eta);
      EXPECT_DOUBLE_EQ(tau * tau, t * t - z * z);
      EXPECT_DOUBLE_EQ(exp(2.0 * eta), (t + z) / (t - z));
    }
  }

  // for the transformation from t-z to tau-eta (momentum)
  EXPECT_DOUBLE_EQ(0.0, lqf.get_ptau(0.0, 0.0, 0.5));
  EXPECT_DOUBLE_EQ(5.0, lqf.get_ptau(5.0, 0.0, 0.0));
  EXPECT_DOUBLE_EQ(0.0, lqf.get_peta(0.0, 0.0, 0.5));
  EXPECT_DOUBLE_EQ(5.0, lqf.get_peta(0.0, 5.0, 0.0));
}

// check causality
TEST(CausalLiquifierTest, TEST_CAUSALITY) {
  CausalLiquefier lqf(0.3, 0.3, 0.3, 0.3);

  // check zero source outside of the cousal area
  double time = lqf.tau_delay;
  double r_bound = lqf.c_diff * time;
  for (int scale = 1.0; scale < 5.0; scale++) {
    double r_out = r_bound * scale;
    EXPECT_DOUBLE_EQ(0.0, lqf.rho_smooth(time, r_out));
    EXPECT_DOUBLE_EQ(0.0, lqf.rho_delta(time, r_out));
    EXPECT_DOUBLE_EQ(0.0, lqf.kernel_rho(time, r_out));
  }

  // check zero source before the deposition time
  time = 0.99 * (lqf.tau_delay - 0.5 * lqf.dtau);
  std::array<Jetscape::real, 4> jmu = {0.0, 0.0, 0.0, 0.0};
  std::array<Jetscape::real, 4> droplet_xmu = {0.0, 0.0, 0.0, 0.0};
  std::array<Jetscape::real, 4> droplet_pmu = {1.0, 1.0, 0.0, 0.0};
  Droplet drop_i(droplet_xmu, droplet_pmu);
  lqf.smearing_kernel(time, 0.0, 0.0, 0.0, drop_i, jmu);
  EXPECT_DOUBLE_EQ(0.0, jmu[0]);

  // check zero source after the deposition time
  time = 1.01 * (lqf.tau_delay + 0.5 * lqf.dtau);
  lqf.smearing_kernel(time, 0.0, 0.0, 0.0, drop_i, jmu);
  EXPECT_DOUBLE_EQ(0.0, jmu[0]);
}

// check conservation law
TEST(CausalLiquifierTest, TEST_CONSERVATION) {
  double dr = 0.005;  // in [fm]
  double dt = 0.005;  // in [fm]

  CausalLiquefier lqf(dt, dr, dr, dr);

  for (double t = 0.3; t < 3.0; t += dt) {
    double integrated_value = 0.0;
    for (double r = 0.5 * dr; r < 1.5 * t; r += dr) {
      integrated_value += 4.0 * M_PI * r * r * dr * lqf.kernel_rho(t, r);
    }
    // require conservation within 5%
    EXPECT_NEAR(1.0, integrated_value, 0.05);
  }
}

// check conservation law
TEST(CausalLiquifierTest, TEST_GRID_CARTESIAN_CONSERVATION) {
  double ll[3] = {0.05, 0.1, 0.3};

  for (int i = 0; i < 0; i++) {
    double dt = 0.1;
    double l = ll[i];
    double dx = l;
    double dy = l;
    double dz = l;
    int n_cells = 5.0 / l;
    double dV = dx * dx * dz;

    CausalLiquefier lqf(dt, dx, dy, dz);

    // std::o/Users/yasukitachibana/Dropbox
    // (Personal)/Codes/Release2.2CandidateLiquefierUpdate/examples/unittests/causal_liquifier.ccfstream
    // ofs; string filename = "liq_cons_test_dx" + std::to_string(int(l*1000)) +
    // "_tr3000_d2400_delta300.txt"; ofs.open(filename.c_str(),
    // std::ios_base::out);

    for (double t = 0.5; t < 3.0; t += dt) {
      double integrated_value = 0.0;
      for (int ix = 0; ix < n_cells; ix++) {
        double x = (ix - 0.5 * double(n_cells)) * dx;
        for (int iy = 0; iy < n_cells; iy++) {
          double y = (iy - 0.5 * double(n_cells)) * dy;
          for (int iz = 0; iz < n_cells; iz++) {
            double z = (iz - 0.5 * double(n_cells)) * dz;

            double r = sqrt(x * x + y * y + z * z);
            integrated_value += lqf.kernel_rho(t, r) * dV;
          }
        }
      }
      // require conservation within 5%
      // ofs << t << " " << integrated_value <<"\n";
      EXPECT_NEAR(1.0, integrated_value, 0.05);
      // JSINFO <<"t="<<t <<" " <<integrated_value;
    }
    // ofs.close();
  }
}

// check conservation law
TEST(CausalLiquifierTest, TEST_GRID_TAU_ETA_CONSERVATION) {
  std::array<Jetscape::real, 4> x_in = {1.0, 0.0, 0.0, 1.0};  // in tau-eta
  std::array<Jetscape::real, 4> p_in = {1.0, 1.0, 1.0, 1.0};  // in Cartesian
  Droplet a_drop(x_in, p_in);

  double dtau = 0.1;
  double dl = 0.05;
  double dx = dl;
  double dy = dl;
  double deta = 0.05;
  int n_xy = 8.0 / dl;
  int n_eta = 5.0 / deta;

  CausalLiquefier lqf(dtau, dx, dy, deta);

  //    std::ofstream ofs;
  //    string filename =
  //    "liq_cons_tau_eta_test_tdep"+std::to_string(int(x_in[0]))+"_etadep"+std::to_string(int(x_in[3]))+"_dxy"
  //    + std::to_string(int(dl*1000)) +
  //    + "_deta" + std::to_string(int(deta*1000)) + "_tr100_d80_delta100.txt";
  //
  //    ofs.open(filename.c_str(), std::ios_base::out);

  double dt = 0.5;
  for (double t = 0.2; t < 4.0; t += dt) {
    double integrated_value = 0.0;
    double tau_delay = t;
    lqf.set_t_delay(tau_delay);
    // in tau-eta
    std::array<Jetscape::real, 4> total_pmu = {0.0, 0.0, 0.0, 0.0};
    std::array<Jetscape::real, 4> x_hydro = {0.0, 0.0, 0.0, 0.0};
    x_hydro[0] = x_in[0] + tau_delay;
    double dvolume = x_hydro[0] * dx * dy * deta;

    for (int ix = 0; ix < n_xy; ix++) {
      x_hydro[1] = (ix - 0.5 * double(n_xy)) * dx;
      for (int iy = 0; iy < n_xy; iy++) {
        x_hydro[2] = (iy - 0.5 * double(n_xy)) * dy;
        for (int ieta = 0; ieta < n_eta; ieta++) {
          x_hydro[3] = (ieta - 0.5 * double(n_eta)) * deta;

          std::array<Jetscape::real, 4> jmu = {0.0, 0.0, 0.0, 0.0};

          lqf.smearing_kernel(x_hydro[0], x_hydro[1], x_hydro[2], x_hydro[3],
                              a_drop, jmu);

          total_pmu[0] +=
              dtau * (jmu[0] * cosh(x_hydro[3]) + jmu[3] * sinh(x_hydro[3])) *
              dvolume;
          total_pmu[1] += dtau * jmu[1] * dvolume;
          total_pmu[2] += dtau * jmu[2] * dvolume;
          total_pmu[3] +=
              dtau * (jmu[0] * sinh(x_hydro[3]) + jmu[3] * cosh(x_hydro[3])) *
              dvolume;

          integrated_value += dtau * jmu[1] * dvolume;
        }
      }
    }
    //        ofs << t << " " << integrated_value <<"\n";
    //        JSINFO
    //        << total_pmu[0] << " "
    //        << total_pmu[1] << " "
    //        << total_pmu[2] << " "
    //        << total_pmu[3];
    //        ofs.close();
    if (t > 3.) {
      EXPECT_NEAR(1.0, integrated_value, 0.05);
    }
  }
}

// ── normalization on the hydro grid ─────────────────────────────────────────
// MUSIC samples get_source() once per cell centre, at the query times tau_n and
// tau_n + dtau of each step (Runge-Kutta, weight dtau/2 each), and adds
// tau dtau dx dy deta J^mu per cell.  With a production grid (0.3 fm, deta 0.2)
// the point-sampled kernel does not sum to the droplet's four-momentum (6.19x,
// 0.98x, 0.0125x for the droplets below, taken from a 0-10% Au+Au run); after
// LiquefierBase::normalize_active_droplets() it does.

namespace {

// the production liquefier: <dtau> 0.02, tau_delay 1, time_relax 0.1,
// d_diff 0.08, width_delta 0.1
void set_production_parameters(CausalLiquefier &lqf) {
  lqf.tau_delay = 1.0;
  lqf.time_relax = 0.1;
  lqf.d_diff = 0.08;
  lqf.width_delta = 0.1;
  lqf.c_diff = sqrt(lqf.d_diff / lqf.time_relax);
  lqf.gamma_relax = 0.5 / lqf.time_relax;
}

// MUSIC's grid for the IS grid grid_max_x 15, grid_step_x 0.3, grid_max_z 6,
// grid_step_z 0.2 (cell centres -size/2 + i d)
HydroGrid music_grid(double dx = 0.3, int nx = 100) {
  HydroGrid g;
  g.nx = g.ny = nx;
  g.neta = 60;
  g.dx = g.dy = dx;
  g.deta = 0.2;
  g.x_min = g.y_min = -0.5 * nx * dx;
  g.eta_min = -6.0;
  return g;
}

// Run the liquefier through the steps around the droplet's deposit as MUSIC
// does, and return the lab four-momentum MUSIC receives.
std::array<double, 4> deposited(CausalLiquefier &lqf, const Droplet &drop,
                                const HydroGrid &g, double tau0 = 0.4) {
  const double dt = lqf.dtau;
  const long k_dep = std::lround((lqf.deposit_time(drop) - tau0) / dt);
  std::array<double, 4> P = {0., 0., 0., 0.};
  const double dV = g.dx * g.dy * g.deta;
  for (long n = k_dep - 4; n <= k_dep + 3; n++) {
    const double tau_n = tau0 + n * dt;
    lqf.prepare_active_droplets(tau_n - dt, tau_n + 2. * dt);
    lqf.set_hydro_grid(g);
    lqf.normalize_active_droplets(tau_n, dt);
    for (int stage = 0; stage < 2; stage++) {
      const double tau = tau_n + stage * dt;
      for (int ie = 0; ie < g.neta; ie++) {
        const double eta = g.eta_min + ie * g.deta;
        const double ch = cosh(eta), sh = sinh(eta);
        for (int ix = 0; ix < g.nx; ix++) {
          for (int iy = 0; iy < g.ny; iy++) {
            std::array<Jetscape::real, 4> j = {0., 0., 0., 0.};
            lqf.get_source(tau, g.x_min + ix * g.dx, g.y_min + iy * g.dy, eta,
                           j);
            const double w = 0.5 * dt * tau * dV;
            P[0] += w * (ch * j[0] + sh * j[3]);
            P[1] += w * j[1];
            P[2] += w * j[2];
            P[3] += w * (sh * j[0] + ch * j[3]);
          }
        }
      }
    }
  }
  return P;
}

struct Case {
  std::array<Jetscape::real, 4> x, p;
  double raw_flux;
};

// tau_d, x, y, eta_d and E, px, py, pz of three production droplets
const Case kCases[] = {
    {{1.05349f, 1.06941f, -3.74047f, -2.23919f},
     {12.7132f, 2.63322f, -0.447677f, -12.3924f}, 6.1897},  // large |eta_s|
    {{0.133141f, 0.0040877f, -3.68418f, -0.498067f},
     {1.68015f, 0.582329f, -1.41034f, -0.473678f}, 0.9775},  // resolved
    {{8.4334f, 4.29754f, 4.26602f, 0.530885f},
     {3.41514f, 2.47391f, 1.43562f, 1.81253f}, 0.012520},  // late tau
};

}  // namespace

// Without normalization MUSIC receives raw_flux x p (the bug); with it, p.
TEST(CausalLiquifierTest, TEST_NORMALIZED_ON_HYDRO_GRID) {
  const HydroGrid g = music_grid();
  for (const auto &c : kCases) {
    const Droplet drop(c.x, c.p);
    for (bool on : {false, true}) {
      CausalLiquefier lqf(0.02, 0.3, 0.3, 0.2);
      set_production_parameters(lqf);
      lqf.set_normalize_on_hydro_grid(on);
      lqf.add_a_droplet(drop);
      const auto P = deposited(lqf, drop, g);
      const double f = on ? 1.0 : c.raw_flux;
      for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(f * c.p[i], P[i], 2e-3 * std::abs(c.p[0]) * std::max(f, 1.0))
            << "component " << i << ", normalize " << on << ", eta_d " << c.x[3];
      }
      if (on) {
        EXPECT_NEAR(c.raw_flux, lqf.get_droplet_flux(0), 1e-3 * c.raw_flux);
      }
    }
  }
}

// All droplets of one step together: each is normalized on its own.
TEST(CausalLiquifierTest, TEST_NORMALIZED_SEVERAL_DROPLETS) {
  const HydroGrid g = music_grid();
  CausalLiquefier lqf(0.02, 0.3, 0.3, 0.2);
  set_production_parameters(lqf);
  // the large-|eta| droplet and a copy moved to eta 0, same deposit time
  Case a = kCases[0], b = kCases[0];
  b.x[3] = 0.0f;
  lqf.add_a_droplet(Droplet(a.x, a.p));
  lqf.add_a_droplet(Droplet(b.x, b.p));
  const auto P = deposited(lqf, Droplet(a.x, a.p), g);
  for (int i = 0; i < 4; i++)
    EXPECT_NEAR(a.p[i] + b.p[i], P[i], 4e-3 * a.p[0]) << "component " << i;
}

// A kernel that misses every cell centre (3 fm cells) is deposited whole into
// the nearest cell; one off the grid is lost and deposits nothing.
TEST(CausalLiquifierTest, TEST_NORMALIZED_POINT_AND_LOST) {
  const HydroGrid g = music_grid(3.0, 10);  // centres at -15 + 3 i
  const Case c = {{1.0f, -1.5f, -1.5f, 0.0f}, {5.0f, 1.0f, -2.0f, 0.5f}, 0.};
  {
    CausalLiquefier lqf(0.02, 3.0, 3.0, 0.2);
    set_production_parameters(lqf);
    const Droplet drop(c.x, c.p);
    lqf.add_a_droplet(drop);
    const auto P = deposited(lqf, drop, g);
    EXPECT_EQ(0.0, lqf.get_droplet_flux(0));
    for (int i = 0; i < 4; i++)
      EXPECT_NEAR(c.p[i], P[i], 1e-4 * c.p[0]) << "component " << i;
  }
  {
    CausalLiquefier lqf(0.02, 0.3, 0.3, 0.2);
    set_production_parameters(lqf);
    Case far = c;
    far.x[1] = 100.0f;
    const Droplet drop(far.x, far.p);
    lqf.add_a_droplet(drop);
    const auto P = deposited(lqf, drop, music_grid());
    for (int i = 0; i < 4; i++)
      EXPECT_EQ(0.0, P[i]);
    EXPECT_NE(std::string::npos,
              lqf.normalization_summary().find("1 lost off the grid"));
  }
}

// A droplet at the hard vertex (tau_d = 0; add_hydro_sources() stores eta = 0
// there, 0/0 before) deposits normally and, normalized, exactly.
TEST(CausalLiquifierTest, TEST_NORMALIZED_AT_THE_VERTEX) {
  const HydroGrid g = music_grid();
  const Case c = {{0.0f, 0.0268f, 1.2035f, 0.0f}, {3.2287f, -2.3416f, 0.4344f, 1.7589f}, 0.};
  CausalLiquefier lqf(0.02, 0.3, 0.3, 0.2);
  set_production_parameters(lqf);
  const Droplet drop(c.x, c.p);
  lqf.add_a_droplet(drop);
  const auto P = deposited(lqf, drop, g);
  EXPECT_GT(lqf.get_droplet_flux(0), 0.5);   // the kernel itself deposits it
  for (int i = 0; i < 4; i++)
    EXPECT_NEAR(c.p[i], P[i], 2e-3 * c.p[0]) << "component " << i;
}
