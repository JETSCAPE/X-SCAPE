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

  HadronicEMT hEMT(0.3,0.3,0.3);

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

TEST(HadronicEMTTest, TEST_TMUNU){

  // create fake hadrons
  std::vector<Hadron> hadron_list;
  unsigned int nparticles = 10;
  for (unsigned int ipart = 0; ipart < nparticles; ipart++) {
    const int hadron_label = 0;
    const int hadron_status = 11;
    const int hadron_id = 111;
    const double hadron_mass = 1.0;
    const double pz = 1.0;
    const double energy = std::sqrt(hadron_mass*hadron_mass + pz*pz);
    FourVector hadron_p(pz, 0.0, 0.0, energy);
    FourVector hadron_x(0.0, 0.0, 0.0, 0.0);

    // create a JETSCAPE Hadron
    hadron_list.push_back(Hadron(hadron_label, hadron_id, hadron_status, 
                                hadron_p, hadron_x, hadron_mass));
  }

  std::unique_ptr<BulkMediaInfo> bulk_info_ptr;

  HadronicEMT hEMT(0.3,0.3,0.3);
  double dV = 0.027; // fm^3

  hEMT.GetBulkInfo(0.0,0.0,0.0,0.0,bulk_info_ptr,hadron_list);
  double energy_density = bulk_info_ptr->energy_density;
  EXPECT_NEAR(nparticles*std::sqrt(2.0)/dV,energy_density,1e-4);

  // add one more particle outside of the box and check that it is not taken
  // into account
  for (unsigned int ipart = 0; ipart < 1; ipart++) {
    const int hadron_label = 0;
    const int hadron_status = 11;
    const int hadron_id = 111;
    const double hadron_mass = 1.0;
    const double px = 1.0;
    const double energy = std::sqrt(hadron_mass*hadron_mass + px*px);
    FourVector hadron_p(px, 0.0, 0.0, energy);
    FourVector hadron_x(1.0, 0.0, 0.0, 0.0);

    // create a JETSCAPE Hadron
    hadron_list.push_back(Hadron(hadron_label, hadron_id, hadron_status, 
                                hadron_p, hadron_x, hadron_mass));
  }

  hEMT.GetBulkInfo(0.0,0.0,0.0,0.0,bulk_info_ptr,hadron_list);
  double energy_density1 = bulk_info_ptr->energy_density;
  EXPECT_NEAR(nparticles*std::sqrt(2.0)/dV,energy_density1,1e-4);

  // check the rest of the energy momentum tensor
  double T00 = bulk_info_ptr->tmn[0][0];
  EXPECT_NEAR(nparticles*std::sqrt(2.0)/dV,T00,1e-4);

  double T01 = bulk_info_ptr->tmn[0][1];
  EXPECT_NEAR(nparticles*1.0/dV,T01,1e-4);

  double T02 = bulk_info_ptr->tmn[0][2];
  EXPECT_NEAR(0.0,T02,1e-4);

  double T03 = bulk_info_ptr->tmn[0][3];
  EXPECT_NEAR(0.0,T03,1e-4);

  double T11 = bulk_info_ptr->tmn[1][1];
  EXPECT_NEAR(nparticles*1.0/std::sqrt(2.0)/dV,T11,1e-4);

  double T12 = bulk_info_ptr->tmn[1][2];
  EXPECT_NEAR(0.0,T12,1e-4);

  double T13 = bulk_info_ptr->tmn[1][3];
  EXPECT_NEAR(0.0,T13,1e-4);

  double T22 = bulk_info_ptr->tmn[2][2];
  EXPECT_NEAR(0.0,T22,1e-4);

  double T23 = bulk_info_ptr->tmn[2][3];
  EXPECT_NEAR(0.0,T23,1e-4);

  double T33 = bulk_info_ptr->tmn[3][3];
  EXPECT_NEAR(0.0,T33,1e-4);

  EXPECT_DOUBLE_EQ(T01,bulk_info_ptr->tmn[1][0]);
  EXPECT_DOUBLE_EQ(T02,bulk_info_ptr->tmn[2][0]);
  EXPECT_DOUBLE_EQ(T12,bulk_info_ptr->tmn[2][1]);
  EXPECT_DOUBLE_EQ(T03,bulk_info_ptr->tmn[3][0]);
  EXPECT_DOUBLE_EQ(T13,bulk_info_ptr->tmn[3][1]);
  EXPECT_DOUBLE_EQ(T23,bulk_info_ptr->tmn[3][2]);
}