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
#include "JetScapeSignalManager.h"
#include "BulkMediaInfo.h"

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


namespace Jetscape {

/**
 * @brief Default constructor for the HadronicEMT class.
 * 
 * This constructor initializes the HadronicEMT object by calling the
 * InitTask() method to retrieve parameter values from the XML configuration
 * file. It allocates memory for GSL matrices and vectors used in the
 * calculations, and sets up the energy density and flow velocity
 * tables.
 */
HadronicEMT::HadronicEMT(){
  VERBOSE(8);
  InitTask(); // Get values of parameters from XML 

  Tmualpha_gsl_ = gsl_matrix_alloc(4, 4);
  eigenvalues_ = gsl_vector_complex_alloc(4);
  eigenvectors_ = gsl_matrix_complex_alloc(4, 4);
  w_ = gsl_eigen_nonsymmv_alloc(4);

  e_lower_bound_ = 0.0;
  e_upper_bound_ = 0.0;
  e_length_ = 100000;
  e_spacing_ = 0.0;
  T_table_.resize(e_length_);
  s_table_.resize(e_length_);

  read_EOS_from_file("./EOS/hotQCD/hrg_hotqcd_eos_SMASH_binary.dat");
}

/**
 * @brief Constructor for the HadronicEMT class used for unit testing.
 * 
 * This constructor initializes the HadronicEMT object with specified
 * transverse and longitudinal smearing parameters, and a covariance
 * flag. It also allocates memory for GSL matrices and vectors used
 * in the calculations.
 */
HadronicEMT::HadronicEMT(double sigma_transverse, double sigma_longitudinal, 
    int smearing_covariant)
    : sigma_transverse_(sigma_transverse),
    sigma_longitudinal_(sigma_longitudinal),
    smearing_covariant_(smearing_covariant){
  Tmualpha_gsl_ = gsl_matrix_alloc(4, 4);
  eigenvalues_ = gsl_vector_complex_alloc(4);
  eigenvectors_ = gsl_matrix_complex_alloc(4, 4);
  w_ = gsl_eigen_nonsymmv_alloc(4);

  e_lower_bound_ = 0.0;
  e_upper_bound_ = 0.0;
  e_length_ = 100000;
  e_spacing_ = 0.0;
  T_table_.resize(e_length_);
  s_table_.resize(e_length_);

  read_EOS_from_file("./EOS/hotQCD/hrg_hotqcd_eos_SMASH_binary.dat");

  // set the grid specifications from the initial state module by hand for this
  // test constructor
  xMax_ = 15.0;
  yMax_ = 15.0;
  zMax_ = 15.0;
  dx_ = 0.3;
  dy_ = 0.3;
  dz_ = 0.3;
  Nx_ = int(std::ceil(2 * xMax_ / dx_));
  Ny_ = int(std::ceil(2 * yMax_ / dy_));
  Nz_ = int(std::ceil(2 * zMax_ / dz_));
}

/**
 * @brief Destructor for the HadronicEMT class.
 * 
 * This destructor releases the memory allocated for the GSL matrices and
 * vectors used in the HadronicEMT calculations. It ensures that all
 * resources are properly freed to prevent memory leaks.
 */
HadronicEMT::~HadronicEMT() {
  gsl_matrix_free(Tmualpha_gsl_);
  gsl_vector_complex_free(eigenvalues_);
  gsl_matrix_complex_free(eigenvectors_);
  gsl_eigen_nonsymmv_free(w_);
}

/**
 * @brief Initializes the HadronicEMT object by reading parameters from XML.
 * 
 * This function retrieves the transverse and longitudinal smearing
 * parameters, the smearing covariance flag, and the grid specifications
 * from the XML configuration file. It sets up the necessary parameters
 * for the HadronicEMT object.
 */
void HadronicEMT::InitTask(){
  JSINFO << "Initialize HadronicEMT ...";

  sigma_transverse_ = JetScapeXML::Instance()->GetElementDouble({"HadronicEMT", "sigma_transverse"});
  sigma_longitudinal_ = JetScapeXML::Instance()->GetElementDouble({"HadronicEMT", "sigma_longitudinal"});
  smearing_covariant_ = JetScapeXML::Instance()->GetElementInt({"HadronicEMT", "smearing_covariant"});

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
}

/** 
 * @brief Computes the energy density and flow velocity from the stress-energy 
 * tensor T^{\mu\nu}.
 *
 * This function takes the stress-energy tensor T^{\mu\nu} for a single 
 * space-time point and computes the energy density (e) and flow velocity 
 * (u^\mu).
 * It solves the equation T^{\mu}_{\nu} u^{\nu} = e u^{\mu} using GSL for 
 * eigenvalue computation.
 * 
 * @param Tmn Reference to a 4x4 array representing the stress-energy tensor 
 * T^{\mu\nu}.
 * @param e Output parameter where the computed energy density will be stored.
 * @param umu Output parameter where the computed flow velocity u^\mu will be 
 * stored.
 * 
 * @details
 * The function multiplies T^{\mu\nu} with a metric tensor g_{\nu\alpha} to 
 * obtain T^{\mu}_{\alpha}, initializes T^{\mu}_{\alpha} as a 4x4 array of 
 * zeros, and then solves the eigenvalue problem using GSL. If the eigenvalue 
 * computation fails (returns GSL_FAILURE), it sets e to 0.0 and umu to 
 * {1., 0., 0., 0.}. Otherwise, it parses the eigenvalues and eigenvectors to 
 * determine the energy density and flow velocity.
 * 
 * @note The function assumes the existence of the following member variables:
 *   - gsl_matrix_complex *Tmualpha_gsl_: GSL matrix representing 
 *     T^{\mu}_{\alpha}.
 *   - gsl_vector_complex *eigenvalues_: GSL vector storing the eigenvalues.
 *   - gsl_matrix_complex *eigenvectors_: GSL matrix storing the eigenvectors.
 *   - gsl_eigen_nonsymmv_workspace *w_: GSL workspace for eigenvalue 
 *     computation.
 *   - g: Metric tensor g_{\nu\alpha} (assumed to be defined externally).
 * 
 * @see ParseEigenSystem for the function that parses eigenvalues and 
 * eigenvectors.
 */
void HadronicEMT::ComputeEnergyDensityAndFlowVelocity(
    const std::array<std::array<double, 4>, 4> &Tmn, double &e,
    std::array<double, 4> &umu) {

    // multiply T^{\mu\nu} with g_{\nu\alpha} to get T^{\mu}_{\alpha}
    // intialize T^{\mu}_{\alpha} as a 4x4 array of zeros
    std::array<std::array<double, 4>, 4> Tmualpha = {{{0., 0., 0., 0.},
                                                      {0., 0., 0., 0.},
                                                      {0., 0., 0., 0.},
                                                      {0., 0., 0., 0.}}};

    // solve T^{\mu}_{\nu} u^{\nu} = e u^{\mu} using GSL
    // fill the gsl matrix with the values of T^{\mu}_{\alpha}
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                Tmualpha[i][j] += Tmn[i][k] * g_[k][j];
            }
            gsl_matrix_set(Tmualpha_gsl_, i, j, Tmualpha[i][j]);
        }
    }

    // solve the eigenvalue problem
    int status = gsl_eigen_nonsymmv(Tmualpha_gsl_, eigenvalues_, 
                                    eigenvectors_, w_);

    // error handling
    if (status != GSL_SUCCESS) {
        VERBOSE(9) << "HadronicEMT: Eigenvalue problem could not be solved!";
        e = 0.0;
        umu = {1., 0., 0., 0.};
        return;
    }

    // parse the eigenvalues and eigenvectors to determine the eigenvalue with 
    // the a positive real part and the corresponding eigenvector
    // this eigenvector is the flow velocity u^{\mu} and the eigenvalue is 
    // the energy density e
    ParseEigenSystem(eigenvectors_, eigenvalues_, umu, e);
}

/** 
 * @brief Parses the eigenvalues and eigenvectors to determine the energy 
 * density and flow velocity.
 *
 * This function parses the eigenvalues and eigenvectors obtained from solving 
 * the eigenvalue problem T^{\mu}_{\nu} u^{\nu} = e u^{\mu} using GSL. It 
 * identifies the eigenvalue with a positive real part and a corresponding 
 * time-like eigenvector, which represents the flow velocity u^\mu.
 * 
 * @param eigen_vectors GSL matrix containing the eigenvectors of the eigenvalue 
 * problem.
 * @param eigen_values GSL vector containing the eigenvalues of the eigenvalue 
 * problem.
 * @param umu Output parameter where the computed flow velocity u^\mu will be 
 * stored.
 * @param e Output parameter where the computed energy density will be stored.
 * 
 * @details
 * The function iterates over the eigenvalues and eigenvectors to find a real 
 * eigenvalue that is positive and has a time-like eigenvector (where the 
 * Minkowski norm is positive). It checks for uniqueness of the solution and 
 * handles cases where multiple candidates exist. If no valid eigenvalue and 
 * eigenvector pair is found, it sets e to 0.0 and umu to {1., 0., 0., 0.}.
 * 
 * @note The function assumes that the eigenvalue problem has already been 
 * solved and the GSL matrices and vectors (eigenvectors and eigenvalues) are 
 * properly initialized and filled.
 */
void HadronicEMT::ParseEigenSystem(gsl_matrix_complex *eigen_vectors,
                        gsl_vector_complex *eigen_values,
                        std::array<double, 4> &umu, double &e) {
  // an exact solution corresponds to one real positive eigenvalue and a 
  // time-like eigenvector

  const double imag_tolerance = 1e-30;
  bool found_a_candidate = false;
  bool too_many_candidates = false;
  double eval_candidate = 0.0;
  std::array<double, 4> evec_candidate = {0.0, 0.0, 0.0, 0.0};
  double evec_candidate_norm_sq = 0.0;

  // loop over the 4 eigenvalues & eigenvectors
  for (int i = 0; i < 4; i++) {
    // get the real and imaginary part of the eigenvalue and store them in double
    double eval_real_part = GSL_REAL(gsl_vector_complex_get(eigen_values, i));
    double eval_imag_part = GSL_IMAG(gsl_vector_complex_get(eigen_values, i));

    // continue if the eigenvalue is real and positive
    if (eval_real_part > 0.0 && (fabs(eval_imag_part/eval_real_part) < imag_tolerance)) {
      // take the corresponding eigenvector to check it
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column(eigen_vectors, i);

      bool all_real = true;
      double norm_sq = 0.0;
      double evec[4];
      for (int j = 0; j < 4; j++) {
        // get real and imag elements of the eigenvector
        const double evec_real = GSL_REAL(gsl_matrix_complex_get(eigen_vectors, j, i));
        const double evec_imag = GSL_IMAG(gsl_matrix_complex_get(eigen_vectors, j, i));

        // check if the eigenvector is real
        if (evec_imag == 0. || (fabs(evec_imag/evec_real) < imag_tolerance)) {
          // compute Minkowski norm of the eigenvector
          if (j == 0) {
            norm_sq += evec_real * evec_real;
          } else {
            norm_sq -= evec_real * evec_real;
          }
          evec[j] = evec_real;
        } else {
          all_real = false;
          break;
        }
      }

      // for a time-like real eigenvector we have reached the solution
      if ((norm_sq > 0) && all_real) {
        // if there is another candidate, we might have a problem here..
        if (found_a_candidate) {
          const double similarity = 1e-3;
          // if the other candidate only differs in the imaginary part, then we
          // do not have a problem
          if ((fabs(1-eval_candidate/eval_real_part) > similarity) 
              && (eval_candidate != eval_real_part)
              && (fabs(evec_candidate[0]/evec[0]) > similarity)
              && (evec_candidate[0] != evec[0])
              && (fabs(evec_candidate[1]/evec[1]) > similarity)
              && (evec_candidate[1] != evec[1])
              && (fabs(evec_candidate[2]/evec[2]) > similarity)
              && (evec_candidate[2] != evec[2])
              && (fabs(evec_candidate[3]/evec[3]) > similarity)
              && (evec_candidate[3] != evec[3])){
            too_many_candidates = true;
          }
        }

        eval_candidate = eval_real_part;
        evec_candidate_norm_sq = norm_sq;
        evec_candidate = {evec[0], evec[1], evec[2], evec[3]};
        found_a_candidate = true;
      }
    }
  }

  if (found_a_candidate && !too_many_candidates) {
    e = eval_candidate;
    // assumption that time component is non-zero
    const double u0_sign = evec_candidate[0] / fabs(evec_candidate[0]);

    for (int i = 0; i < 4; i++) {
      umu[i] = u0_sign * evec_candidate[i] / sqrt(evec_candidate_norm_sq);
    }
  } else {
    e = 0.0;
    umu = {1., 0., 0., 0.};
  }
}

/**
 * @brief Resets the energy-momentum tensor to its initial state.
 *
 * This function initializes the `Tmn_requested_point_` tensor to a zero 4x4 
 * matrix.
 */
void HadronicEMT::ResetEnergyMomentumTensor() {
  Tmn_requested_point_ = {{{0., 0., 0., 0.},
                          {0., 0., 0., 0.},
                          {0., 0., 0., 0.},
                          {0., 0., 0., 0.}}};
}

/**
 * @brief Retrieves the current energy-momentum tensor.
 *
 * This function returns the current state of the energy-momentum tensor, 
 * represented as a 4x4 matrix of doubles. The tensor encapsulates the energy 
 * density, momentum density, and stress tensor components for the system.
 *
 * @return A 4x4 `std::array` representing the energy-momentum tensor.
 *
 * @note The returned tensor reflects the most recent state of 
 * `Tmn_requested_point_`.
 */
std::array<std::array<double, 4>, 4> HadronicEMT::GetEnergyMomentumTensor() const {
  return Tmn_requested_point_;
}

/**
 * @brief Adds the contribution of a particle to the energy-momentum tensor.
 *
 * This function updates the components of the energy-momentum tensor 
 * (`Tmn_requested_point_`) by incorporating the scaled four-momentum 
 * contributions of a given particle. The scaling factor is applied to the 
 * particle's momentum components before updating the tensor.
 *
 * @param hadron The particle whose four-momentum will contribute to the tensor. 
 * @param scaling_factor A factor to scale the particle's momentum and energy 
 *                       before adding it to the tensor.
 *
 * @note Ensure that the input `hadron` contains valid four-momentum data.
 */
void HadronicEMT::AddParticleToEnergyMomentumTensor(Hadron &hadron, 
                                                    double scaling_factor) {
  const FourVector p = hadron.p_in();
  const double e = p.t() * scaling_factor;
  const double px = p.x() * scaling_factor;
  const double py = p.y() * scaling_factor;
  const double pz = p.z() * scaling_factor;

  // add the scaled particle quantities to the T^{\mu\nu} tensor
  Tmn_requested_point_[0][0] += e;
  Tmn_requested_point_[0][1] += px;
  Tmn_requested_point_[0][2] += py;
  Tmn_requested_point_[0][3] += pz;
  Tmn_requested_point_[1][1] += px * px / e;
  Tmn_requested_point_[1][2] += px * py / e;
  Tmn_requested_point_[1][3] += px * pz / e;
  Tmn_requested_point_[2][2] += py * py / e;
  Tmn_requested_point_[2][3] += py * pz / e;
  Tmn_requested_point_[3][3] += pz * pz / e;
  
  Tmn_requested_point_[1][0] += px;
  Tmn_requested_point_[2][0] += py;
  Tmn_requested_point_[3][0] += pz;
  Tmn_requested_point_[2][1] += py * px / e;
  Tmn_requested_point_[3][1] += pz * px / e;
  Tmn_requested_point_[3][2] += pz * py / e;

}

/**
 * @brief Computes a covariant smearing kernel in Cartesian coordinates.
 *
 * This function evaluates a smearing kernel based on the relative position 
 * and four-velocity in Cartesian coordinates, using a Gaussian profile 
 * modified for Lorentz covariance.
 *
 * @param x_diff The x-component of the position difference.
 * @param y_diff The y-component of the position difference.
 * @param z_diff The z-component of the position difference.
 * @param ux The x-component of the four-velocity.
 * @param uy The y-component of the four-velocity.
 * @param uz The z-component of the four-velocity.
 * @param gamma The Lorentz factor (\( \gamma = \frac{1}{\sqrt{1 - v^2}} \)).
 *
 * @return The value of the smearing kernel at the given position and velocity.
 */
double HadronicEMT::smearing_kernel_covariant_Cartesian(
    const double x_diff, const double y_diff, const double z_diff,
    const double ux, const double uy, const double uz,
    const double gamma) {

  const double N = gamma / (pow(M_PI, 1.5) * sigma_transverse_ *
                   sigma_transverse_ * sigma_transverse_);
  // compute the squared distance and the scalar product r*u
  const double dr_squared =
      (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
  const double dr_dot_u = (x_diff * ux + y_diff * uy + z_diff * uz);

  return N * exp(-(dr_squared + dr_dot_u * dr_dot_u) /
                 (sigma_transverse_ * sigma_transverse_));
}

/**
 * @brief Computes a smearing kernel using a Gaussian profile.
 *
 * This function evaluates a smearing kernel in three dimensions, 
 * assuming Gaussian profiles for the transverse and longitudinal directions.
 *
 * @param x_diff The x-component of the position difference.
 * @param y_diff The y-component of the position difference.
 * @param z_diff The z-component of the position difference (longitudinal).
 *
 * @return The value of the smearing kernel at the given position.
 */
double HadronicEMT::smearing_kernel_gaussian(
    const double x_diff, const double y_diff, const double z_diff) {
  const double N_z = 1. / (sqrt(M_PI) * sigma_longitudinal_);
  const double N_trans = 1. / (M_PI * sigma_transverse_ * sigma_transverse_);

  const double r_trans_squared = x_diff * x_diff + y_diff * y_diff;
  const double deta_squared = z_diff * z_diff;

  const double exp_trans =
      N_trans * exp(-r_trans_squared / (sigma_transverse_ * sigma_transverse_));
  const double exp_z =
      N_z * exp(-deta_squared / (sigma_longitudinal_ * sigma_longitudinal_));
  return exp_trans * exp_z;
}

/**
 * @brief Calculates and fills the bulk medium properties at a given space-time point.
 *
 * This function computes the energy-momentum tensor \( T^{\mu\nu} \) and 
 * other bulk medium properties (e.g., energy density, flow velocity, 
 * temperature, entropy density) based on a list of hadrons and 
 * smearing kernels.
 *
 * @param t The time coordinate of the requested point (ignored).
 * @param x The x-coordinate of the requested point.
 * @param y The y-coordinate of the requested point.
 * @param z The z-coordinate of the requested point.
 * @param bulk_info_ptr A unique pointer to a `BulkMediaInfo` object to store the computed results.
 * @param h_list A vector of `Hadron` objects representing the hadrons in the vicinity of the point.
 *
 * @details
 * - Resets the energy-momentum tensor to zero using `ResetEnergyMomentumTensor()`.
 * - Loops over all hadrons in `h_list` to compute their contribution to the energy-momentum tensor 
 *   using either a covariant or Gaussian smearing kernel.
 * - Contributions are ignored for hadrons beyond a specified transverse or longitudinal skip distance.
 * - After computing the contributions, the energy-momentum tensor is stored in `bulk_info_ptr`.
 * - If the energy density is below a threshold (`rounding_error`), the bulk properties are set to zero.
 * - Otherwise, the energy density and flow velocity are computed using `ComputeEnergyDensityAndFlowVelocity`.
 * - The temperature and entropy density are obtained using the equation of state (EOS) functions `get_T` and `get_s`.
 *
 * @note The covariant smearing kernel is used if `smearing_covariant_` is set; otherwise, a Gaussian smearing kernel is applied.
 * 
 * @warning Assumes that all hadrons are at the same time step during the (half-)concurrent evolution.
 */
void HadronicEMT::GetBulkInfo(Jetscape::real t, Jetscape::real x,
                              Jetscape::real y, Jetscape::real z,
                              std::unique_ptr<BulkMediaInfo> &bulk_info_ptr,
                              std::vector<Hadron> &h_list) {

  bulk_info_ptr = std::make_unique<BulkMediaInfo>();

  ResetEnergyMomentumTensor();

  double skip_distance_transverse_ = 5.0 * sigma_transverse_;
  double skip_distance_longitudinal_ = 5.0 * sigma_longitudinal_;
  if (smearing_covariant_ == 1) {
    skip_distance_longitudinal_ = skip_distance_transverse_;
  }
  
  // loop over particles
  for (auto &ihad: h_list) {
    const FourVector r = ihad.x_in();
    // x,y,z of hadron
    // The time is ignored here, assumed to be at the same time step in the
    // (half-)concurrent evolution
    const double x_had = r.x();
    const double y_had = r.y();
    const double z_had = r.z();

    const double x_diff = x_had - x;
    if (std::abs(x_diff) > skip_distance_transverse_) {
      continue;
    }
    const double y_diff = y_had - y;
    if (std::abs(y_diff) > skip_distance_transverse_) {
      continue;
    }
    const double z_diff = z_had - z;
    if (std::abs(z_diff) > skip_distance_longitudinal_) {
      continue;
    }

    // get the momentum of the hadron
    const FourVector p = ihad.p_in();
    const double e = p.t();
    const double px = p.x();
    const double py = p.y();
    const double pz = p.z();
    // mass of the hadron
    const double m = sqrt(e * e - px * px - py * py - pz * pz);

    double value_kernel;
    if (smearing_covariant_) {
      const double ux = px / m;
      const double uy = py / m;
      const double uz = pz / m;
      const double gamma = e / m;
      value_kernel = smearing_kernel_covariant_Cartesian(x_diff, y_diff, z_diff,
                                                            ux, uy, uz, gamma);
    } else {
      value_kernel = smearing_kernel_gaussian(x_diff, y_diff, z_diff);
    }

    // give reference to the hadron and the value of the kernel
    AddParticleToEnergyMomentumTensor(ihad, value_kernel);
  }

  // fill T^{\mu\nu} in the bulk media info
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      bulk_info_ptr->tmn[i][j] = static_cast<Jetscape::real>(Tmn_requested_point_[i][j]);
    }
  }

  if (bulk_info_ptr->tmn[0][0] < rounding_error) {
    bulk_info_ptr->energy_density = 0.0;
    bulk_info_ptr->vx = 0.0;
    bulk_info_ptr->vy = 0.0;
    bulk_info_ptr->vz = 0.0;
    bulk_info_ptr->temperature = 0.0;
    bulk_info_ptr->entropy_density = 0.0;
    return;
  }

  // all particles added, then compute the energy density and flow velocity
  ComputeEnergyDensityAndFlowVelocity(Tmn_requested_point_, e_requested_point_, umu_requested_point_);

  // provide the energy density
  bulk_info_ptr->energy_density = e_requested_point_;
  // provide the flow velocity
  bulk_info_ptr->vx = umu_requested_point_[1];
  bulk_info_ptr->vy = umu_requested_point_[2];
  bulk_info_ptr->vz = umu_requested_point_[3];

  // provide the temperature and entropy density from the EOS
  bulk_info_ptr->temperature = get_T(e_requested_point_);
  bulk_info_ptr->entropy_density = get_s(e_requested_point_);
}

/**
 * @brief Determines which hadrons should be fluidized based on the local temperature of the medium.
 *
 * This function evaluates a 3D grid of spatial points and computes the bulk medium properties 
 * (e.g., temperature) at each point. If the temperature at a grid point exceeds the critical 
 * temperature \( T_{\text{critical}} \), hadrons near that point are flagged for fluidization.
 *
 * @param T_critical The critical temperature above which hadrons are flagged for fluidization.
 * @param current_hadrons A vector of `Hadron` objects representing the current hadron population.
 * @return A vector of booleans, where each entry corresponds to a hadron in `current_hadrons`. 
 *         `true` indicates that the corresponding hadron should be fluidized.
 *
 * @details
 * - The function initializes a boolean vector, `fluidize_hadrons`, with the same size as `current_hadrons`, 
 *   setting all entries to `false`.
 * - Loops over a 3D grid of spatial points with ranges defined by `Nx_`, `Ny_`, and `Nz_` and spacing 
 *   determined by `dx_`, `dy_`, and `dz_`.
 * - At each grid point, `GetBulkInfo` is called to compute the bulk medium properties.
 * - If the temperature at a grid point exceeds `T_critical`, all hadrons within a specified 
 *   smearing range around that point are flagged for fluidization.
 * - The smearing ranges are determined by `sigma_transverse_` and `sigma_longitudinal_`, with adjustments 
 *   for covariant smearing if `smearing_covariant_` is enabled.
 *
 * @note 
 * - Covariant smearing adjusts the longitudinal smearing range to match the transverse range.
 * - Assumes all hadrons and grid points are evaluated at the same time step (time coordinate fixed at \( t = 0 \)).
 *
 * @warning 
 * - The function loops over all grid points and hadrons, which may be computationally expensive for large grids 
 *   or hadron populations.
 */
std::vector<bool> HadronicEMT::DetermineHadronsForFluidization(double T_critical,
          std::vector<Hadron> &current_hadrons,
          std::shared_ptr<FluidDynamics> fluid_dynamics_ptr) {
  // vector of bools to determine which hadrons to fluidize
  std::vector<bool> fluidize_hadrons(current_hadrons.size(), false);

  double skip_distance_transverse_ = 5.0 * sigma_transverse_;
  double skip_distance_longitudinal_ = 5.0 * sigma_longitudinal_;
  if (smearing_covariant_ == 1) {
    skip_distance_longitudinal_ = skip_distance_transverse_;
  }

  std::unique_ptr<FluidCellInfo> fluid_cell_info_ptr;
  if (fluid_dynamics_ptr) {
    fluid_cell_info_ptr = std::make_unique<FluidCellInfo>();
  }

  // loop over the whole 3D grid
  for (int ix = 0; ix < Nx_; ix++) {
    const double x = -xMax_ + ix * dx_;
    for (int iy = 0; iy < Ny_; iy++) {
      const double y = -yMax_ + iy * dy_;
      for (int iz = 0; iz < Nz_; iz++) {
        const double z = -zMax_ + iz * dz_;

        bool above_Tc = false;

        if (fluid_dynamics_ptr) {
// TODO: This needs the hydro medium info at the current time step
          fluid_dynamics_ptr->GetHydroInfo(0.0, x, y, z, fluid_cell_info_ptr);
          above_Tc = (fluid_cell_info_ptr->temperature > T_critical);
        } else {
          std::unique_ptr<BulkMediaInfo> bulk_info_ptr;
          // get the bulk media info for the current space-time point
          // t = 0 is a dummy, because it is the current hadron list
          GetBulkInfo(0.0, x, y, z, bulk_info_ptr, current_hadrons);
          above_Tc = (bulk_info_ptr->temperature > T_critical);
        }

        if (!above_Tc) continue;

        int hadron_index = 0;
        // loop over particle positions and determine the ones to fluidize
        for (auto &ihad : current_hadrons) {
          const FourVector r = ihad.x_in();
          const double x_had = r.x();
          const double y_had = r.y();
          const double z_had = r.z();

          const double x_diff = x_had - x;
          const double y_diff = y_had - y;
          const double z_diff = z_had - z;

          if (std::abs(x_diff) <= skip_distance_transverse_ &&
              std::abs(y_diff) <= skip_distance_transverse_ &&
              std::abs(z_diff) <= skip_distance_longitudinal_) {
            fluidize_hadrons[hadron_index] = true;
          }

          hadron_index++;
        }
      }
    }
  }
  return fluidize_hadrons;
}

/**
 * @brief Reads the equation of state (EOS) data from a binary file.
 *
 * This function reads the EOS data from a specified binary file and populates 
 * the internal tables for energy density, temperature, and entropy density.
 *
 * @param filename The name of the binary file containing the EOS data.
 *
 * @throws std::runtime_error if the file cannot be opened or read.
 *
 * @details
 * - The function opens the specified binary file and reads the EOS data into 
 *   internal vectors `T_table_` and `s_table_`.
 * - The energy density is converted from units of GeV/fm^3 to 1/fm^4.
 * - The temperature is stored as \( T^5 \) for smoother interpolation.
 */
void HadronicEMT::read_EOS_from_file(const std::string &filename) {
  // Read the EOS table from file
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("HadronicEMT: Could not open EOS file " + filename);
  }

  double temp;
  for (int i = 0; i < e_length_; i++) {
    file.read(reinterpret_cast<char*>(&temp), sizeof(double));  // e
    temp /= hbarC;      // 1/fm^4
    if (i == 0) {
      e_lower_bound_ = temp;
    }
    if (i == 1) {
      e_spacing_ = temp - e_lower_bound_;
    }
    if (i == e_length_ - 1) {
      e_upper_bound_ = temp;
    }
    file.read(reinterpret_cast<char*>(&temp), sizeof(double));  // P, not used
    file.read(reinterpret_cast<char*>(&temp), sizeof(double));  // s
    s_table_[i] = temp;

    file.read(reinterpret_cast<char*>(&temp), sizeof(double));  // T
    temp /= hbarC;   // 1/fm
    T_table_[i] = std::pow(temp, 5);   // store T^5 for smooth curve
  }
}

/**
 * @brief Interpolates the equation of state (EOS) data for a given energy density.
 *
 * This function performs linear interpolation on the EOS data to obtain the 
 * temperature or entropy density at a specified energy density.
 *
 * @param e The energy density in 1/fm^4.
 * @param table The EOS table to interpolate (either temperature or entropy).
 * @return The interpolated value from the EOS table.
 *
 * @details
 * - The function computes the index for the energy density and clamps it within 
 *   the valid range of the table.
 * - It handles potential underflow and overflow cases by returning appropriate 
 *   values.
 */
double HadronicEMT::interpolate_1D_EOS(double e, const std::vector<double> &table) const {
  // This is a generic linear interpolation routine for EOS at zero mu_B
  // it assumes the input table has the following structure:
  //       e, P(e), T(e), s(e)
  // as one-dimensional arrays on an equally spacing lattice grid
  // units: e is in 1/fm^4

  // Compute the index and clamp it within the valid range
  int idx_e = std::clamp(static_cast<int>((e - e_lower_bound_) / e_spacing_),
                          0, e_length_ - 2);

  // Handle potential underflow
  if (e < e_lower_bound_) {
    return table[0] * e / e_lower_bound_;
  }

  // Handle potential overflow
  if (e > e_upper_bound_) {
    return table[e_length_ - 1];
  }

  // Linear interpolation
  const double frac_e = (e - (idx_e * e_spacing_ + e_lower_bound_)) 
                      / e_spacing_;
  const double temp1 = table[idx_e];
  const double temp2 = table[idx_e + 1];
  return temp1 * (1.0 - frac_e) + temp2 * frac_e;
}

/** 
 * @brief Computes the local temperature from the energy density.
 *
 * This function calculates the local temperature based on the energy density 
 * using the equation of state (EOS) data.
 *
 * @param e The energy density in GeV/fm^3.
 * @return The local temperature in GeV.
 */
double HadronicEMT::get_T(double e) const {
  const double e_fm4 = e / hbarC;  // 1/fm^4
  const double T5 = interpolate_1D_EOS(e_fm4, T_table_);  // returns e/T^5
  const double T = pow(T5, 0.2) * hbarC;  // GeV
  return T;
}

/**
 * @brief Computes the local entropy density from the energy density.
 * 
 * This function calculates the local entropy density based on the energy
 * density using the equation of state (EOS) data.
 * 
 * @param e The energy density in GeV/fm^3.
 * @return The local entropy density in GeV/fm^3.
 */
double HadronicEMT::get_s(double e) const {
  double e_fm4 = e / hbarC;  // 1/fm^4
  return interpolate_1D_EOS(e_fm4, s_table_);
}

/**
 * @brief Computes the Cartesian time coordinate from the proper time and pseudorapidity.
 * 
 * This function calculates the Cartesian time coordinate \( t \) based on the
 * proper time \( \tau \) and pseudorapidity \( \eta \).
 * 
 * @param tau The proper time in fm.
 * @param eta The pseudorapidity.
 * @return The Cartesian time coordinate \( t \) in fm.
 */
double HadronicEMT::get_t(double tau, double eta) const {
  return tau*cosh(eta);
}

/**
 * @brief Computes the Cartesian coordinate \( z \) from the proper time and pseudorapidity.
 * 
 * This function calculates the Cartesian coordinate \( z \) based on the
 * proper time \( \tau \) and pseudorapidity \( \eta \).
 * 
 * @param tau The proper time in fm.
 * @param eta The pseudorapidity.
 * @return The Cartesian coordinate \( z \) in fm.
 */
double HadronicEMT::get_z(double tau, double eta) const {
  return tau*sinh(eta);
}

/**
 * @brief Computes the proper time component of a four vector.
 * 
 * This function calculates the proper time component \( \tau \) of a four vector
 * based on the energy \( p^0 \), momentum \( p^3 \), and pseudorapidity \( \eta \).
 * 
 * @param p0 The energy component of the four vector.
 * @param p3 The momentum component of the four vector.
 * @param eta The pseudorapidity.
 * @return The proper time component \( \tau \) of the four vector.
 */
double HadronicEMT::get_ptau(double p0, double p3, double eta) const {
  return p0*cosh(eta) - p3*sinh(eta);
}

/**
 * @brief Computes the pseudorapidity component of a four vector.
 * 
 * This function calculates the pseudorapidity component \( \eta \) of a four vector
 * based on the energy \( p^0 \), momentum \( p^3 \), and proper time \( \tau \).
 * 
 * @param p0 The energy component of the four vector.
 * @param p3 The momentum component of the four vector.
 * @param tau The proper time.
 * @return The pseudorapidity component \( \eta \) of the four vector.
 */
double HadronicEMT::get_peta(double p0, double p3, double eta) const {
  return p3*cosh(eta) - p0*sinh(eta);
}

};