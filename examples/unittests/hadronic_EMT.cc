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

#include "HadronicEMT.h"
#include "gtest/gtest.h"

using namespace Jetscape;

// check coordinate transformation functions (same as in causal_liquifier.cc test)
TEST(HadronicEMTTest, TEST_COORDINATES){

  HadronicEMT hEMT(0.5,0.5,1);

  // for the transformation from tau-eta to t-z (configuration space)
  EXPECT_DOUBLE_EQ(0.0, hEMT.get_t(0.0,0.5));
  EXPECT_DOUBLE_EQ(5.0, hEMT.get_t(5.0,0.0));
  EXPECT_DOUBLE_EQ(0.0, hEMT.get_z(0.0,0.5));
  EXPECT_DOUBLE_EQ(0.0, hEMT.get_z(5.0,0.0));

  for(double tau=0.1;tau<0.3;tau+=0.1){
    for(double eta=0.0;eta<0.2;eta+=0.1){
      double t=hEMT.get_t(tau,eta);
      double z=hEMT.get_z(tau,eta);
      EXPECT_DOUBLE_EQ(tau*tau, t*t-z*z);
      EXPECT_DOUBLE_EQ(exp(2.0*eta), (t+z)/(t-z));
    }
  }

  // for the transformation from t-z to tau-eta (momentum)
  EXPECT_DOUBLE_EQ(0.0, hEMT.get_ptau(0.0,0.0,0.5));
  EXPECT_DOUBLE_EQ(5.0, hEMT.get_ptau(5.0,0.0,0.0));
  EXPECT_DOUBLE_EQ(0.0, hEMT.get_peta(0.0,0.0,0.5));
  EXPECT_DOUBLE_EQ(5.0, hEMT.get_peta(0.0,5.0,0.0));
}

TEST(HadronicEMTTest, TEST_LANDAU_MATCHING){
  // check the Landau matching for a simple case
  const std::array<std::array<double, 4>, 4> Tmn = {{{1.0, 0.0, 0.0, 0.0},
                                                    {0.0, 1./3., 0.0, 0.0},
                                                    {0.0, 0.0, 1./3., 0.0},
                                                    {0.0, 0.0, 0.0, 1./3.}}};
  double e = 0.0;
  std::array<double, 4> umu = {0.0, 0.0, 0.0, 0.0};

  HadronicEMT hEMT(0.5,0.5,1);
  hEMT.ComputeEnergyDensityAndFlowVelocity(Tmn, e, umu);

  EXPECT_NEAR(1.0,e,1e-4);
  EXPECT_NEAR(1.0,umu[0],1e-4);
  EXPECT_NEAR(0.0,umu[1],1e-4);
  EXPECT_NEAR(0.0,umu[2],1e-4);
  EXPECT_NEAR(0.0,umu[3],1e-4);

  // check the Landau matching for a more complicated case
  // This tensor is created with T^{\mu\nu} = (e + P) u^{\mu} u^{\nu} - P g^{\mu\nu}
  // e = 1.0, P = e/3, u^{\mu} = (u^0,0.1,0.1,0.1), u^0 = sqrt(1.0 + 0.1^2 + 0.1^2 + 0.1^2) = 1.01488916
  const std::array<std::array<double, 4>, 4> Tmn1 = {{{1.04, 0.13531855, 0.13531855, 0.13531855},
                                                    {0.13531855, 0.34666667, 0.01333333, 0.01333333},
                                                    {0.13531855, 0.01333333, 0.34666667, 0.01333333},
                                                    {0.13531855, 0.01333333, 0.01333333, 0.34666667}}};
  double e1 = 0.0;
  std::array<double, 4> umu1 = {0.0, 0.0, 0.0, 0.0};

  hEMT.ComputeEnergyDensityAndFlowVelocity(Tmn1, e1, umu1);

  EXPECT_NEAR(1.0,e1,1e-4);
  EXPECT_NEAR(1.01488916,umu1[0],1e-4);
  EXPECT_NEAR(0.1,umu1[1],1e-4);
  EXPECT_NEAR(0.1,umu1[2],1e-4);
  EXPECT_NEAR(0.1,umu1[3],1e-4);
}

TEST(HadronicEMTTest, TEST_TMUNU){
  // create fake hadrons
  std::vector<Hadron> hadron_list;
  unsigned int nparticles = 10;
  for (unsigned int ipart = 0; ipart < nparticles; ipart++) {
    const int hadron_label = 0;
    const int hadron_status = 11;
    const int hadron_id = 111;
    const double hadron_mass = 1.0;
    const double energy = hadron_mass;
    FourVector hadron_p(0.0, 0.0, 0.0, energy);
    FourVector hadron_x(0.0, 0.0, 0.0, 0.0);

    // create a JETSCAPE Hadron
    hadron_list.push_back(Hadron(hadron_label, hadron_id, hadron_status, 
                                hadron_p, hadron_x, hadron_mass));
  }

  std::unique_ptr<BulkMediaInfo> bulk_info_ptr;

  HadronicEMT hEMT_cov(0.5,0.5,1);
  hEMT_cov.GetBulkInfo(0.0,0.0,0.0,0.0,bulk_info_ptr,hadron_list);

  // value of the kernel is 1.4367
  EXPECT_NEAR(nparticles*1.4367,bulk_info_ptr->energy_density,1e-4);
  EXPECT_NEAR(0.0,bulk_info_ptr->vx,1e-4);
  EXPECT_NEAR(0.0,bulk_info_ptr->vy,1e-4);
  EXPECT_NEAR(0.0,bulk_info_ptr->vz,1e-4);

  
  // place half of the particles outside of the 5 sigma range
  // create fake hadrons
  std::vector<Hadron> hadron_list1;
  unsigned int nparticles1 = 5;
  for (unsigned int ipart = 0; ipart < nparticles1; ipart++) {
    const int hadron_label = 0;
    const int hadron_status = 11;
    const int hadron_id = 111;
    const double hadron_mass = 1.0;
    const double energy = hadron_mass;
    FourVector hadron_p(0.0, 0.0, 0.0, energy);
    FourVector hadron_x(0.0, 0.0, 0.0, 0.0);

    // create a JETSCAPE Hadron
    hadron_list1.push_back(Hadron(hadron_label, hadron_id, hadron_status, 
                                hadron_p, hadron_x, hadron_mass));
  }
  for (unsigned int ipart = 0; ipart < nparticles1; ipart++) {
    const int hadron_label = 0;
    const int hadron_status = 11;
    const int hadron_id = 111;
    const double hadron_mass = 1.0;
    const double energy = hadron_mass;
    FourVector hadron_p(0.0, 0.0, 0.0, energy);
    FourVector hadron_x(0.0, 50.0, 0.0, 0.0);

    // create a JETSCAPE Hadron
    hadron_list1.push_back(Hadron(hadron_label, hadron_id, hadron_status, 
                                hadron_p, hadron_x, hadron_mass));
  }

  hEMT_cov.GetBulkInfo(0.0,0.0,0.0,0.0,bulk_info_ptr,hadron_list1);

  // value of the kernel is 1.4367
  EXPECT_NEAR(nparticles1*1.4367,bulk_info_ptr->energy_density,1e-4);
  EXPECT_NEAR(0.0,bulk_info_ptr->vx,1e-4);
  EXPECT_NEAR(0.0,bulk_info_ptr->vy,1e-4);
  EXPECT_NEAR(0.0,bulk_info_ptr->vz,1e-4);


  // test the first setup with the Gaussian kernel
  HadronicEMT hEMT_gauss(0.5,0.5,0);
  hEMT_gauss.GetBulkInfo(0.0,0.0,0.0,0.0,bulk_info_ptr,hadron_list);

  // value of the kernel is 1.4367
  EXPECT_NEAR(nparticles*1.4367,bulk_info_ptr->energy_density,1e-4);
  EXPECT_NEAR(0.0,bulk_info_ptr->vx,1e-4);
  EXPECT_NEAR(0.0,bulk_info_ptr->vy,1e-4);
  EXPECT_NEAR(0.0,bulk_info_ptr->vz,1e-4);
}