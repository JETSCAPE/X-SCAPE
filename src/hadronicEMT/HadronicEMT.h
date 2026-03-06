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

#ifndef HADRONICEMT_H
#define HADRONICEMT_H

#include "sigslot.h"
#include "BulkMediaInfo.h"
#include "JetScapeParticles.h"
#include "FluidDynamics.h"
#include <memory>
#include <array>
#include <gsl/gsl_eigen.h>

namespace Jetscape {

class HadronicEMT {
 private:
  double sigma_transverse_, sigma_longitudinal_;
  int smearing_covariant_;
  /// Initial state parameters (grid specifications)
  double xMax_, yMax_, zMax_, dx_, dy_, dz_;
  int Nx_, Ny_, Nz_;

  // EOS variables
  double e_lower_bound_, e_upper_bound_, e_spacing_;
  int e_length_;
  std::vector<double> T_table_;
  std::vector<double> s_table_;

  // define the (pseudo)-metric tensor g_{\mu\nu} as a 4x4 array
  const std::array<std::array<int, 4>, 4> g_ = {
      {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}}};

  // define the inverse (pseudo)-metric tensor g^{\mu\nu} as a 4x4 array
  const std::array<std::array<int, 4>, 4> g_inv_ = {
      {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}}};

  gsl_matrix *Tmualpha_gsl_;
  gsl_vector_complex *eigenvalues_;
  gsl_matrix_complex *eigenvectors_;
  gsl_eigen_nonsymmv_workspace *w_;

  // define a 4x4 array to store the energy-momentum tensor
  std::array<std::array<double, 4>, 4> Tmn_requested_point_;
  // define a length 4 array to store the flow velocity
  std::array<double, 4> umu_requested_point_;
  // define a variable to store the energy density
  double e_requested_point_;

 public:
  HadronicEMT();
  HadronicEMT(double sigma_transverse, double sigma_longitudinal,
              int smearing_covariant);
  ~HadronicEMT();

  void InitTask();

  /// Fill in bulk media info for (t,x,y,z) from current hadron list
  void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                   Jetscape::real z,
                   std::unique_ptr<BulkMediaInfo> &bulk_info_ptr,
                   std::vector<Hadron> &h_list);

  /// Determine particles in list to fluidize
  std::vector<bool> DetermineHadronsForFluidization(
      double T_critical, std::vector<Hadron> &current_hadrons,
      std::shared_ptr<FluidDynamics> fluid_dynamics_ptr = nullptr);

  void ComputeEnergyDensityAndFlowVelocity(
      const std::array<std::array<double, 4>, 4> &Tmn, double &e,
      std::array<double, 4> &umu);

  void ParseEigenSystem(gsl_matrix_complex *eigen_vectors,
                        gsl_vector_complex *eigen_values,
                        std::array<double, 4> &umu, double &e);

  void ResetEnergyMomentumTensor();
  std::array<std::array<double, 4>, 4> GetEnergyMomentumTensor() const;

  void AddParticleToEnergyMomentumTensor(Hadron &hadron, double scaling_factor);

  double smearing_kernel_covariant_Cartesian(
      const double x_diff, const double y_diff, const double z_diff,
      const double ux, const double uy, const double uz, const double gamma);

  double smearing_kernel_gaussian(const double x_diff, const double y_diff,
                                  const double z_diff);

  // setter functions for unit testing
  void set_lower_bound(double e_lower_bound) { e_lower_bound_ = e_lower_bound; }
  void set_upper_bound(double e_upper_bound) { e_upper_bound_ = e_upper_bound; }
  void set_spacing(double e_spacing) { e_spacing_ = e_spacing; }
  void set_length(int e_length) { e_length_ = e_length; }
  void set_T_table(const std::vector<double> &T_table) { T_table_ = T_table; }
  void set_s_table(const std::vector<double> &s_table) { s_table_ = s_table; }

  void read_EOS_from_file(const std::string &filename);
  double interpolate_1D_EOS(double e, const std::vector<double> &table) const;
  double get_T(double e) const;
  double get_s(double e) const;

  double get_t(double tau, double eta) const;
  double get_z(double tau, double eta) const;

  double get_ptau(double px, double pz, double eta) const;
  double get_peta(double px, double pz, double eta) const;
};

};  // namespace Jetscape

#endif  // HADRONICEMT_H