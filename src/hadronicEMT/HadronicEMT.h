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

#ifndef HADRONICEMT_H
#define HADRONICEMT_H

#include "sigslot.h"
#include "BulkMediaInfo.h"
#include "JetScapeParticles.h"
#include <memory>
#include <array>
#include <gsl/gsl_eigen.h>

namespace Jetscape {

class HadronicEMT {
private:
  double sigma_transverse_, sigma_longitudinal_;
  int smearing_covariant_;

  // define the (pseudo)-metric tensor g_{\mu\nu} as a 4x4 array
  const std::array<std::array<int, 4>, 4> g = {{{1, 0, 0, 0},
                                                {0, -1, 0, 0},
                                                {0, 0, -1, 0},
                                                {0, 0, 0, -1}}};

  // define the inverse (pseudo)-metric tensor g^{\mu\nu} as a 4x4 array
  const std::array<std::array<int, 4>, 4> g_inv = {{{1, 0, 0, 0},
                                                    {0, -1, 0, 0},
                                                    {0, 0, -1, 0},
                                                    {0, 0, 0, -1}}};

  gsl_matrix *Tmualpha_gsl;
  gsl_vector_complex *eigenvalues;
  gsl_matrix_complex *eigenvectors;
  gsl_eigen_nonsymmv_workspace *w;

  // define a 4x4 array to store the energy-momentum tensor
  std::array<std::array<double, 4>, 4> Tmn_requested_point;
  // define a length 4 array to store the flow velocity
  std::array<double, 4> umu_requested_point;
  // define a variable to store the energy density
  double e_requested_point;


public:
  HadronicEMT();
  HadronicEMT(double sigma_transverse, double sigma_longitudinal, 
              int smearing_covariant);
  ~HadronicEMT();

  void InitTask();

  /// Fill in bulk media info for (t,x,y,z) from current hadron list
  void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y, 
            Jetscape::real z, std::unique_ptr<BulkMediaInfo> &bulk_info_ptr,
            std::vector<Hadron> &h_list);

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
    const double ux, const double uy, const double uz,
    const double gamma);

  double smearing_kernel_gaussian(
    const double x_diff, const double y_diff, const double z_diff);

  double get_t(double tau, double eta) const;
  double get_z(double tau, double eta) const;

  double get_ptau(double px, double pz, double eta) const;
  double get_peta(double px, double pz, double eta) const;


};

}; // namespace Jetscape

#endif // HADRONICEMT_H