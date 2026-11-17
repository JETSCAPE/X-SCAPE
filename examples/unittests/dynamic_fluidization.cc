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

#include <cmath>

#include "DynamicFluidization.h"
#include "gtest/gtest.h"

using namespace Jetscape;

// T^{mu nu} = (e + P) u^mu u^nu - P g^{mu nu}, boosted: the Landau match must
// give e and u back
TEST(DynamicFluidizationTest, LANDAU_MATCH_BOOSTED_IDEAL_FLUID) {
  const double e = 2.3, P = 0.6;
  const double ux = 0.8, uy = -0.3, uz = 1.7;
  const double u0 = std::sqrt(1. + ux * ux + uy * uy + uz * uz);
  const double u[4] = {u0, ux, uy, uz};
  const double g[4] = {1., -1., -1., -1.};
  std::array<std::array<double, 4>, 4> T{};
  for (int a = 0; a < 4; a++)
    for (int b = 0; b < 4; b++)
      T[a][b] = (e + P) * u[a] * u[b] - (a == b ? P * g[a] : 0.);
  double e_out = 0.;
  std::array<double, 4> u_out{};
  ASSERT_TRUE(DynamicFluidization::LandauMatch(T, e_out, u_out));
  EXPECT_NEAR(e_out, e, 1e-9);
  for (int a = 0; a < 4; a++) EXPECT_NEAR(u_out[a], u[a], 1e-8);
}

// A single hadron: its own contribution is subtracted, it never fluidizes
TEST(DynamicFluidizationTest, SINGLE_HADRON_DOES_NOT_FLUIDIZE) {
  DynamicFluidization df(0.5, 0.01);
  std::vector<DynamicFluidization::Candidate> h = {
      {1., 0., 0., 0., 0.938, 0., 0., 0., 0.938, true}};
  auto f = df.Evaluate(h);
  ASSERT_EQ(f.size(), 1u);
  EXPECT_FALSE(f[0]);
  // its rest-frame energy density is its own: m / (2 pi sigma^2)^{3/2}
  EXPECT_NEAR(df.LastEnergyDensity()[0],
              0.938 / std::pow(2. * M_PI * 0.25, 1.5), 1e-12);
}

// A dense cluster fluidizes; a far-away hadron and an ineligible one in the
// cluster do not; every hadron (also the ineligible) contributes to T
TEST(DynamicFluidizationTest, DENSE_CLUSTER_FLUIDIZES) {
  DynamicFluidization df(0.5, 0.5);
  std::vector<DynamicFluidization::Candidate> h;
  for (int i = 0; i < 20; i++) {
    const double x = 0.05 * (i % 5), y = 0.05 * (i / 5), z = 0.;
    const double px = 0.2 * std::cos(i), pz = 0.2 * std::sin(i);
    const double m = 0.14;
    h.push_back({1., x, y, z, std::sqrt(m * m + px * px + pz * pz), px, 0.,
                 pz, m, i != 3});
  }
  h.push_back({1., 20., 0., 0., 0.938, 0., 0., 0., 0.938, true});
  auto f = df.Evaluate(h);
  for (int i = 0; i < 20; i++) {
    if (i == 3) {
      EXPECT_FALSE(f[i]);
    } else {
      EXPECT_TRUE(f[i]) << "hadron " << i << " e = "
                        << df.LastEnergyDensity()[i];
    }
  }
  EXPECT_FALSE(f[20]);
}

// Lorentz invariance of the criterion: boosting the whole cluster along z
// (positions contracted at equal lab time do not matter for a static
// cluster at one point) leaves the rest-frame energy density unchanged
TEST(DynamicFluidizationTest, REST_FRAME_ENERGY_DENSITY_BOOST_INVARIANT) {
  DynamicFluidization df(0.5, 0.5);
  const double m = 1.0;
  std::vector<DynamicFluidization::Candidate> rest = {
      {0., 0., 0., 0., m, 0., 0., 0., m, true},
      {0., 0., 0., 0., m, 0., 0., 0., m, true}};
  df.Evaluate(rest);
  const double e_rest = df.LastEnergyDensity()[0];
  const double gam = 3., bg = std::sqrt(gam * gam - 1.);
  std::vector<DynamicFluidization::Candidate> moving = {
      {0., 0., 0., 0., gam * m, 0., 0., bg * m, m, true},
      {0., 0., 0., 0., gam * m, 0., 0., bg * m, m, true}};
  df.Evaluate(moving);
  EXPECT_NEAR(df.LastEnergyDensity()[0], e_rest, 1e-9 * e_rest);
}
