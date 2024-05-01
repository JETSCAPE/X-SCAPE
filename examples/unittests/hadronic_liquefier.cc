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

#include "HadronicLiquefier.h"
#include "LiquefierBase.h"
#include "JetScapeParticles.h"
#include "FourVector.h"
#include "gtest/gtest.h"

using namespace Jetscape;


// Test the HadronicLiquefier constructor used for further testing
TEST(HadronicLiquefierTest, TEST_CONSTRUCTOR) {
    
    HadronicLiquefier lqf(true, 0.5, 0.5, 15., 15., 15., 100, 100, 100, false);
    
    EXPECT_EQ(true, lqf.get_covariant_smearing());
    EXPECT_EQ(0.5, lqf.get_sigma_transverse());
    EXPECT_EQ(0.5, lqf.get_sigma_longitudinal());
    EXPECT_EQ(15., lqf.get_xMax());
    EXPECT_EQ(15., lqf.get_yMax());
    EXPECT_EQ(15., lqf.get_zMax());
    EXPECT_EQ(100, lqf.get_Nx());
    EXPECT_EQ(100, lqf.get_Ny());
    EXPECT_EQ(100, lqf.get_Nz());
    EXPECT_EQ(false, lqf.get_hydro_Cartesian());
}

/** Test the creation of a Hadron and Droplet object for further testing and the 
 * accessability of its attributes
*/
TEST(HadronicLiquefierTest, CreateValidHadronTest) {

    int label = 1;
    int id = 211;
    int stat = 27;
    // Use the fastjet FourVector convention with the time as last element 
    FourVector p(0.0, 0.0, 0.0, 0.138);
    FourVector x(0.1, 0.2, 0.3, 1.0);
    double mass = 0.138;
    int charge = 0;
    int baryon_number = 0;
    int strangeness = 0;

    Hadron had = Hadron(label, id, stat, p, x, mass, charge, baryon_number, strangeness);

    EXPECT_EQ(1, had.plabel());
    EXPECT_EQ(id, had.pid());
    EXPECT_EQ(stat, had.pstat());
    EXPECT_EQ(true, p==had.p_in());
    EXPECT_EQ(true, x==had.x_in());
    EXPECT_EQ(mass, had.restmass());
    EXPECT_EQ(charge, had.charge());
    EXPECT_EQ(baryon_number, had.baryon_number());
    EXPECT_EQ(strangeness, had.strangeness());


    std::array<double, 4> js_real_array_x = {0.0, 0.0, 0.0, 1.0};
    std::array<double, 4> js_real_array_p = {0.0, 0.0, 0.0, 0.138};
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);

    // test each element of the x and p array
    for (int i = 0; i < 4; i++) {
        EXPECT_EQ(js_real_array_x[i], drop.get_xmu()[i]);
        EXPECT_EQ(js_real_array_p[i], drop.get_pmu()[i]);
    }
    EXPECT_EQ(baryon_number, drop.get_baryon_number());
    EXPECT_EQ(charge, drop.get_electric_charge());
    EXPECT_EQ(strangeness, drop.get_strangeness());
}

/** After testing the creation of the HadronicLiquefier object and the Hadron 
 * object, we can now test the smearing kernel functions of the 
 * HadronicLiquefier.
*/
TEST(HadronicLiquefierTest, SmearingKernelFunctionTest) {

    HadronicLiquefier lqf(true, 0.5, 0.5, 15., 15., 15., 100, 100, 100, false);

    // Test for a trivial point only the normalization factor
    double x_diff = 0.0;
    double y_diff = 0.0;
    double eta_diff = 0.0;
    double ux = 0.0;
    double uy = 0.0;
    double ueta = 0.0;
    double tau = 1.0;
    double gamma = 1.0;
    double expected_result_cov = 1. / (pow(M_PI, 1.5) 
                                * pow(lqf.get_sigma_transverse(), 3.0));
    double expected_result_gauss = 1. / (sqrt(M_PI) * lqf.get_sigma_longitudinal()
                                * M_PI * pow(lqf.get_sigma_transverse(), 2.0));

    double result1 = lqf.smearing_kernel_covariant_Milne(
        x_diff, y_diff, eta_diff, ux, uy, ueta, tau, gamma);

    ASSERT_DOUBLE_EQ(result1, expected_result_cov);

    // use the same arguments, also for the Cartesian smearing kernel every eta
    // is replaced by z
    double result2 = lqf.smearing_kernel_covariant_Cartesian(
        x_diff, y_diff, eta_diff, ux, uy, ueta, gamma);

    ASSERT_DOUBLE_EQ(result2, expected_result_cov);

    // test the Gaussian smearing kernel
    double result3 = lqf.smearing_kernel_gaussian(x_diff, y_diff, eta_diff);

    ASSERT_DOUBLE_EQ(result3, expected_result_gauss);

    // Test for a non-trivial point
    x_diff = 0.2;
    y_diff = 0.1;
    eta_diff = 0.05;
    ux = 0.2;
    uy = 0.3;
    ueta = 0.1;
    tau = 2.0;
    gamma = 1.5;
    double expected_result_cov1 = 3.304743965850829;
    double expected_result_cov2 = 1.7079807435145198;
    double expected_result_gauss1 = 1.1645639357902724;

    double result4 = lqf.smearing_kernel_covariant_Milne(
        x_diff, y_diff, eta_diff, ux, uy, ueta, tau, gamma);

    ASSERT_FLOAT_EQ(result4, expected_result_cov1);

    double result5 = lqf.smearing_kernel_covariant_Cartesian(
        x_diff, y_diff, eta_diff, ux, uy, ueta, gamma);

    ASSERT_FLOAT_EQ(result5, expected_result_cov2);

    double result6 = lqf.smearing_kernel_gaussian(x_diff, y_diff, eta_diff);

    ASSERT_FLOAT_EQ(result6, expected_result_gauss1);
}

TEST(HadronicLiquefierTest, TestComputeDropletKernelNormalizationCovariantMilne) {
    // First test covariant smearing, Milne hydro
    HadronicLiquefier lqf(true, 0.5, 0.5, 15., 15., 15., 128, 128, 128, false);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double t = 1.0;
    std::array<double, 4> js_real_array_x = {t, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 0;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);

    // Expectation is that the value is close to one for a hadron at rest and
    // a fine enough grid 
    double expected_result = 1.;

    double result = lqf.compute_drop_kernel_normalization(1.0, drop);
    double tau = t; // in this simplified case tau = t

    double tolerance = 1.e-6;

    ASSERT_NEAR(result, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestComputeDropletKernelNormalizationCovariantCartesian) {
    // First test covariant smearing, Cartesian hydro
    HadronicLiquefier lqf(true, 0.5, 0.5, 15., 15., 15., 128, 128, 128, true);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 0;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);

    // Expectation is that the value is close to one for a hadron at rest and
    // a fine enough grid 
    double expected_result = 1.;

    double result = lqf.compute_drop_kernel_normalization(0., drop);

    double tolerance = 1.e-6;

    ASSERT_NEAR(result, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestComputeDropletKernelNormalizationGaussMilne) {
    // First test Gaussian smearing, Milne hydro
    HadronicLiquefier lqf(false, 0.5, 0.5, 15., 15., 15., 128, 128, 128, false);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double t = 1.0;
    std::array<double, 4> js_real_array_x = {t, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 0;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);

    // Expectation is that the value is close to one for a hadron at rest and
    // a fine enough grid 
    double expected_result = 1.;

    double result = lqf.compute_drop_kernel_normalization(0., drop);
    double tau = t; // in this simplified case tau = t

    double tolerance = 1.e-6;

    ASSERT_NEAR(result, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestComputeDropletKernelNormalizationGaussCartesian) {
    // First test Gaussian smearing, Cartesian hydro
    HadronicLiquefier lqf(false, 0.5, 0.5, 15., 15., 15., 128, 128, 128, true);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 0;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);

    // Expectation is that the value is close to one for a hadron at rest and
    // a fine enough grid 
    double expected_result = 1.;

    double result = lqf.compute_drop_kernel_normalization(0., drop);

    double tolerance = 1.e-6;

    ASSERT_NEAR(result, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestGetSourceEnergyCovariantMilne) {
    // Create a HadronicLiquefier object for testing
    HadronicLiquefier lqf_cov_Milne(true, 0.5, 0.5, 15., 15., 15., 228, 228, 228, false);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double tau = 1.0;
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 0;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);
    double norm = lqf_cov_Milne.compute_drop_kernel_normalization(tau, drop);
    drop.set_normalization(norm);
    lqf_cov_Milne.add_a_hadronic_droplet(drop);

    // Test the function get_source_energy
    std::array<double, 4> jmu;

    // run over all grid points and integrate the j0 component
    double total_energy = 0.;
    for (int ieta = 0; ieta < lqf_cov_Milne.get_Nz(); ieta++) {
        for (int ix = 0; ix < lqf_cov_Milne.get_Nx(); ix++) {
            for (int iy = 0; iy < lqf_cov_Milne.get_Ny(); iy++) {
                const double eta = -lqf_cov_Milne.get_zMax() + ieta * lqf_cov_Milne.get_dz();
                const double x = -lqf_cov_Milne.get_xMax() + ix * lqf_cov_Milne.get_dx();
                const double y = -lqf_cov_Milne.get_yMax() + iy * lqf_cov_Milne.get_dy();
                
                lqf_cov_Milne.get_source_energy(tau, x, y, eta, jmu);
                total_energy += jmu[0];
            }
        }
    }
    // dz = tau * deta in this case
    total_energy *= (tau * lqf_cov_Milne.get_dx() * lqf_cov_Milne.get_dy() * lqf_cov_Milne.get_dz());

    double tolerance = 1.e-2;
    double expected_result = 0.138;
    ASSERT_NEAR(total_energy, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestGetSourceEnergyCovariantCartesian) {
    // Create a HadronicLiquefier object for testing
    HadronicLiquefier lqf_cov_Cartesian(true, 0.5, 0.5, 15., 15., 15., 128, 128, 128, true);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double t = 1.0;
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 0;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);
    double norm = lqf_cov_Cartesian.compute_drop_kernel_normalization(t, drop);
    drop.set_normalization(norm);
    lqf_cov_Cartesian.add_a_hadronic_droplet(drop);

    // Test the function get_source_energy
    std::array<double, 4> jmu;

    // run over all grid points and integrate the j0 component
    double total_energy = 0.;
    for (int ieta = 0; ieta < lqf_cov_Cartesian.get_Nz(); ieta++) {
        for (int ix = 0; ix < lqf_cov_Cartesian.get_Nx(); ix++) {
            for (int iy = 0; iy < lqf_cov_Cartesian.get_Ny(); iy++) {
                const double eta = -lqf_cov_Cartesian.get_zMax() + ieta * lqf_cov_Cartesian.get_dz();
                const double x = -lqf_cov_Cartesian.get_xMax() + ix * lqf_cov_Cartesian.get_dx();
                const double y = -lqf_cov_Cartesian.get_yMax() + iy * lqf_cov_Cartesian.get_dy();
                
                lqf_cov_Cartesian.get_source_energy(t, x, y, eta, jmu);
                total_energy += jmu[0];
            }
        }
    }
    total_energy *= (lqf_cov_Cartesian.get_dx() * lqf_cov_Cartesian.get_dy() * lqf_cov_Cartesian.get_dz());

    double tolerance = 1.e-2;
    double expected_result = 0.138;
    ASSERT_NEAR(total_energy, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestGetSourceEnergyGaussMilne) {
    // Create a HadronicLiquefier object for testing
    HadronicLiquefier lqf_Gauss_Milne(false, 0.5, 0.5, 15., 15., 15., 128, 128, 128, false);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double tau = 1.0;
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 0;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);
    double norm = lqf_Gauss_Milne.compute_drop_kernel_normalization(tau, drop);
    drop.set_normalization(norm);
    lqf_Gauss_Milne.add_a_hadronic_droplet(drop);

    // Test the function get_source_energy
    std::array<double, 4> jmu;

    // run over all grid points and integrate the j0 component
    double total_energy = 0.;
    for (int ieta = 0; ieta < lqf_Gauss_Milne.get_Nz(); ieta++) {
        for (int ix = 0; ix < lqf_Gauss_Milne.get_Nx(); ix++) {
            for (int iy = 0; iy < lqf_Gauss_Milne.get_Ny(); iy++) {
                const double eta = -lqf_Gauss_Milne.get_zMax() + ieta * lqf_Gauss_Milne.get_dz();
                const double x = -lqf_Gauss_Milne.get_xMax() + ix * lqf_Gauss_Milne.get_dx();
                const double y = -lqf_Gauss_Milne.get_yMax() + iy * lqf_Gauss_Milne.get_dy();
                
                lqf_Gauss_Milne.get_source_energy(tau, x, y, eta, jmu);
                total_energy += jmu[0];
            }
        }
    }
    // dz = tau * deta in this case
    total_energy *= (tau * lqf_Gauss_Milne.get_dx() * lqf_Gauss_Milne.get_dy() * lqf_Gauss_Milne.get_dz());

    double tolerance = 1.e-2;
    double expected_result = 0.138;
    ASSERT_NEAR(total_energy, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestGetSourceEnergyGaussCartesian) {
    // Create a HadronicLiquefier object for testing
    HadronicLiquefier lqf_Gauss_Cartesian(false, 0.5, 0.5, 15., 15., 15., 128, 128, 128, true);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double t = 1.0;
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 0;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);
    double norm = lqf_Gauss_Cartesian.compute_drop_kernel_normalization(t, drop);
    drop.set_normalization(norm);
    lqf_Gauss_Cartesian.add_a_hadronic_droplet(drop);

    // Test the function get_source_energy
    std::array<double, 4> jmu;

    // run over all grid points and integrate the j0 component
    double total_energy = 0.;
    for (int ieta = 0; ieta < lqf_Gauss_Cartesian.get_Nz(); ieta++) {
        for (int ix = 0; ix < lqf_Gauss_Cartesian.get_Nx(); ix++) {
            for (int iy = 0; iy < lqf_Gauss_Cartesian.get_Ny(); iy++) {
                const double eta = -lqf_Gauss_Cartesian.get_zMax() + ieta * lqf_Gauss_Cartesian.get_dz();
                const double x = -lqf_Gauss_Cartesian.get_xMax() + ix * lqf_Gauss_Cartesian.get_dx();
                const double y = -lqf_Gauss_Cartesian.get_yMax() + iy * lqf_Gauss_Cartesian.get_dy();
                
                lqf_Gauss_Cartesian.get_source_energy(t, x, y, eta, jmu);
                total_energy += jmu[0];
            }
        }
    }
    total_energy *= (lqf_Gauss_Cartesian.get_dx() * lqf_Gauss_Cartesian.get_dy() * lqf_Gauss_Cartesian.get_dz());

    double tolerance = 1.e-2;
    double expected_result = 0.138;
    ASSERT_NEAR(total_energy, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestAddHydroSourcesHadrons) {
    int label = 1;
    int id = 211;
    int stat = 27;
    // Use the fastjet FourVector convention with the time as last element 
    FourVector p(0.0, 0.0, 0.0, 0.138);
    FourVector x(0.1, 0.2, 0.3, 1.0);
    double mass = 0.138;
    int charge = 0;
    int baryon_number = 0;
    int strangeness = 0;

    Hadron had = Hadron(label, id, stat, p, x, mass, charge, baryon_number, strangeness);

    // add the same hadron 3 times and check the hadron_droplet_list size
    HadronicLiquefier lqf(true, 0.5, 0.5, 15., 15., 15., 128, 128, 128, false);
    std::vector<Hadron> hIn = {had, had, had};
    lqf.add_hydro_sources_hadrons(1.0, hIn);

    EXPECT_EQ(3, lqf.get_hadron_droplet_list_size());
}

TEST(HadronicLiquefierTest, TestGetSourceQuantityZero) {
    // Create a HadronicLiquefier object for testing
    HadronicLiquefier lqf_cov_Milne(true, 0.5, 0.5, 15., 15., 15., 228, 228, 228, false);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double tau = 1.0;
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 0;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);
    double norm = lqf_cov_Milne.compute_drop_kernel_normalization(tau, drop);
    drop.set_normalization(norm);
    lqf_cov_Milne.add_a_hadronic_droplet(drop);

    // run over all grid points and integrate
    double total_rhob = 0.;
    double total_rhoq = 0.;
    double total_rhos = 0.;
    for (int ieta = 0; ieta < lqf_cov_Milne.get_Nz(); ieta++) {
        for (int ix = 0; ix < lqf_cov_Milne.get_Nx(); ix++) {
            for (int iy = 0; iy < lqf_cov_Milne.get_Ny(); iy++) {
                const double eta = -lqf_cov_Milne.get_zMax() + ieta * lqf_cov_Milne.get_dz();
                const double x = -lqf_cov_Milne.get_xMax() + ix * lqf_cov_Milne.get_dx();
                const double y = -lqf_cov_Milne.get_yMax() + iy * lqf_cov_Milne.get_dy();

                total_rhob += lqf_cov_Milne.get_source_rhob(tau, x, y, eta);
                total_rhoq += lqf_cov_Milne.get_source_rhoq(tau, x, y, eta);
                total_rhos += lqf_cov_Milne.get_source_rhos(tau, x, y, eta);
            }
        }
    }
    // dz = tau * deta in this case
    total_rhob *= (tau * lqf_cov_Milne.get_dx() * lqf_cov_Milne.get_dy() * lqf_cov_Milne.get_dz());
    total_rhoq *= (tau * lqf_cov_Milne.get_dx() * lqf_cov_Milne.get_dy() * lqf_cov_Milne.get_dz());
    total_rhos *= (tau * lqf_cov_Milne.get_dx() * lqf_cov_Milne.get_dy() * lqf_cov_Milne.get_dz());

    double tolerance = 1.e-6;
    double expected_result = 0.0;
    ASSERT_NEAR(total_rhob, expected_result, tolerance);
    ASSERT_NEAR(total_rhoq, expected_result, tolerance);
    ASSERT_NEAR(total_rhos, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestGetSourceQuantityCovariantMilne) {
    // Create a HadronicLiquefier object for testing
    HadronicLiquefier lqf_cov_Milne(true, 0.5, 0.5, 15., 15., 15., 228, 228, 228, false);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double tau = 1.0;
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 1;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);
    double norm = lqf_cov_Milne.compute_drop_kernel_normalization(tau, drop);
    drop.set_normalization(norm);
    lqf_cov_Milne.add_a_hadronic_droplet(drop);

    // run over all grid points and integrate
    double total_rhob = 0.;
    for (int ieta = 0; ieta < lqf_cov_Milne.get_Nz(); ieta++) {
        for (int ix = 0; ix < lqf_cov_Milne.get_Nx(); ix++) {
            for (int iy = 0; iy < lqf_cov_Milne.get_Ny(); iy++) {
                const double eta = -lqf_cov_Milne.get_zMax() + ieta * lqf_cov_Milne.get_dz();
                const double x = -lqf_cov_Milne.get_xMax() + ix * lqf_cov_Milne.get_dx();
                const double y = -lqf_cov_Milne.get_yMax() + iy * lqf_cov_Milne.get_dy();

                total_rhob += lqf_cov_Milne.get_source_rhob(tau, x, y, eta);
            }
        }
    }
    // dz = tau * deta in this case
    total_rhob *= (tau * lqf_cov_Milne.get_dx() * lqf_cov_Milne.get_dy() * lqf_cov_Milne.get_dz());

    double tolerance = 1.e-6;
    double expected_result = 1.0;
    ASSERT_NEAR(total_rhob, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestGetSourceQuantityCovariantCartesian) {
    // Create a HadronicLiquefier object for testing
    HadronicLiquefier lqf_cov_Cartesian(true, 0.5, 0.5, 15., 15., 15., 228, 228, 228, false);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double t = 1.0;
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 1;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);
    double norm = lqf_cov_Cartesian.compute_drop_kernel_normalization(t, drop);
    drop.set_normalization(norm);
    lqf_cov_Cartesian.add_a_hadronic_droplet(drop);

    // run over all grid points and integrate
    double total_rhob = 0.;
    for (int ieta = 0; ieta < lqf_cov_Cartesian.get_Nz(); ieta++) {
        for (int ix = 0; ix < lqf_cov_Cartesian.get_Nx(); ix++) {
            for (int iy = 0; iy < lqf_cov_Cartesian.get_Ny(); iy++) {
                const double eta = -lqf_cov_Cartesian.get_zMax() + ieta * lqf_cov_Cartesian.get_dz();
                const double x = -lqf_cov_Cartesian.get_xMax() + ix * lqf_cov_Cartesian.get_dx();
                const double y = -lqf_cov_Cartesian.get_yMax() + iy * lqf_cov_Cartesian.get_dy();

                total_rhob += lqf_cov_Cartesian.get_source_rhob(t, x, y, eta);
            }
        }
    }
    total_rhob *= (lqf_cov_Cartesian.get_dx() * lqf_cov_Cartesian.get_dy() * lqf_cov_Cartesian.get_dz());

    double tolerance = 1.e-6;
    double expected_result = 1.0;
    ASSERT_NEAR(total_rhob, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestGetSourceQuantityGaussMilne) {
    // Create a HadronicLiquefier object for testing
    HadronicLiquefier lqf_Gauss_Milne(false, 0.5, 0.5, 15., 15., 15., 228, 228, 228, false);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double tau = 1.0;
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 1;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);
    double norm = lqf_Gauss_Milne.compute_drop_kernel_normalization(tau, drop);
    drop.set_normalization(norm);
    lqf_Gauss_Milne.add_a_hadronic_droplet(drop);

    // run over all grid points and integrate
    double total_rhob = 0.;
    for (int ieta = 0; ieta < lqf_Gauss_Milne.get_Nz(); ieta++) {
        for (int ix = 0; ix < lqf_Gauss_Milne.get_Nx(); ix++) {
            for (int iy = 0; iy < lqf_Gauss_Milne.get_Ny(); iy++) {
                const double eta = -lqf_Gauss_Milne.get_zMax() + ieta * lqf_Gauss_Milne.get_dz();
                const double x = -lqf_Gauss_Milne.get_xMax() + ix * lqf_Gauss_Milne.get_dx();
                const double y = -lqf_Gauss_Milne.get_yMax() + iy * lqf_Gauss_Milne.get_dy();

                total_rhob += lqf_Gauss_Milne.get_source_rhob(tau, x, y, eta);
            }
        }
    }
    // dz = tau * deta in this case
    total_rhob *= (tau * lqf_Gauss_Milne.get_dx() * lqf_Gauss_Milne.get_dy() * lqf_Gauss_Milne.get_dz());

    double tolerance = 1.e-6;
    double expected_result = 1.0;
    ASSERT_NEAR(total_rhob, expected_result, tolerance);
}

TEST(HadronicLiquefierTest, TestGetSourceQuantityGaussCartesian) {
    // Create a HadronicLiquefier object for testing
    HadronicLiquefier lqf_Gauss_Cartesian(true, 0.5, 0.5, 15., 15., 15., 228, 228, 228, false);

    // Create a 'trivial' droplet, use the usual FourVector convention, since
    // only arrays are used here
    const double t = 1.0;
    std::array<double, 4> js_real_array_x = {1.0, 0.5, 0.0, 0.0};
    std::array<double, 4> js_real_array_p = {0.138, 0.0, 0.0, 0.0};
    int baryon_number = 1;
    int charge = 0;
    int strangeness = 0;
    HadronDroplet drop = HadronDroplet(js_real_array_x, js_real_array_p, baryon_number, charge, strangeness);
    double norm = lqf_Gauss_Cartesian.compute_drop_kernel_normalization(t, drop);
    drop.set_normalization(norm);
    lqf_Gauss_Cartesian.add_a_hadronic_droplet(drop);

    // run over all grid points and integrate
    double total_rhob = 0.;
    for (int ieta = 0; ieta < lqf_Gauss_Cartesian.get_Nz(); ieta++) {
        for (int ix = 0; ix < lqf_Gauss_Cartesian.get_Nx(); ix++) {
            for (int iy = 0; iy < lqf_Gauss_Cartesian.get_Ny(); iy++) {
                const double eta = -lqf_Gauss_Cartesian.get_zMax() + ieta * lqf_Gauss_Cartesian.get_dz();
                const double x = -lqf_Gauss_Cartesian.get_xMax() + ix * lqf_Gauss_Cartesian.get_dx();
                const double y = -lqf_Gauss_Cartesian.get_yMax() + iy * lqf_Gauss_Cartesian.get_dy();

                total_rhob += lqf_Gauss_Cartesian.get_source_rhob(t, x, y, eta);
            }
        }
    }
    total_rhob *= (lqf_Gauss_Cartesian.get_dx() * lqf_Gauss_Cartesian.get_dy() * lqf_Gauss_Cartesian.get_dz());

    double tolerance = 1.e-6;
    double expected_result = 1.0;
    ASSERT_NEAR(total_rhob, expected_result, tolerance);
}