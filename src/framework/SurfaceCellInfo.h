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
// This is a general basic class for hydrodynamics

#ifndef SURFACECELLINFO_H
#define SURFACECELLINFO_H

#include <string>

#include "RealType.h"

namespace Jetscape {

/**
 * @class SurfaceCellInfo
 * @brief A class representing a single cell on a hyper-surface in the
 * simulation of heavy-ion collisions.
 *
 * This class holds data related to the hydrodynamic properties and conditions
 * of a surface cell in a heavy-ion collision simulation. The information
 * includes variables for energy density, entropy density, temperature, chemical
 * potentials, shear stress tensor, and other physical properties relevant to
 * the hydrodynamics of the system.
 */

class SurfaceCellInfo {
 public:
  // data structure for outputing hyper-surface information
  /**
   * @brief Tau (proper time).
   * @details The proper time of the cell in the surface.
   */
  Jetscape::real tau;

  /**
   * @brief X coordinate.
   * @details The x-coordinate of the cell in the surface.
   */
  Jetscape::real x;

  /**
   * @brief Y coordinate.
   * @details The y-coordinate of the cell in the surface.
   */
  Jetscape::real y;

  /**
   * @brief Eta (spatial rapidity).
   * @details The eta coordinate (spatial rapidity) of the cell.
   */
  Jetscape::real eta;

  /**
   * @brief Surface vector.
   * @details This is a 4-dimensional vector representing the surface element.
   *          It is used in the calculation of the flow dynamics.
   */
  Jetscape::real d3sigma_mu[4];

  /**
   * @brief Energy density.
   * @details Local energy density in the cell, measured in GeV/fm^3.
   */
  Jetscape::real energy_density;

  /**
   * @brief Entropy density.
   * @details Local entropy density in the cell, measured in 1/fm^3.
   */
  Jetscape::real entropy_density;

  /**
   * @brief Temperature.
   * @details Local temperature of the cell, measured in GeV.
   */
  Jetscape::real temperature;

  /**
   * @brief Thermal pressure.
   * @details Thermal pressure within the cell, measured in GeV/fm^3.
   */
  Jetscape::real pressure;

  /**
   * @brief Baryon density.
   * @details Net baryon density in the cell, measured in 1/fm^3.
   */
  Jetscape::real baryon_density;

  /**
   * @brief QGP fraction.
   * @details The fraction of the cell in the QGP+HRG phase, indicating the
   * presence of quark-gluon plasma.
   */
  Jetscape::real qgp_fraction;  //!< Fraction of quark gluon plasma assuming
                                //!< medium is in QGP+HRG phase.

  /**
   * @brief Electric charge density.
   * @details Net electric charge density [1/fm^3].
   */
  Jetscape::real electric_charge_density;

  /**
   * @brief Strangeness density.
   * @details Net strange density [1/fm^3].
   */
  Jetscape::real strangeness_density;

  /**
   * @brief Baryon chemical potential.
   * @details Net baryon chemical potential in the cell, measured in GeV.
   */
  Jetscape::real mu_B;

  /**
   * @brief Charge chemical potential.
   * @details Net charge chemical potential in the cell, measured in GeV.
   */
  Jetscape::real mu_Q;

  /**
   * @brief Strangeness chemical potential.
   * @details Net strangeness chemical potential in the cell, measured in GeV.
   */
  Jetscape::real mu_S;

  /**
   * @brief Flow velocity.
   * @details The flow velocity of the cell in 4 dimensions, used in
   * hydrodynamic simulations.
   */
  Jetscape::real umu[4];

  /**
   * @brief Shear stress tensor.
   * @details Shear stress tensor in the cell, representing the viscosity of the
   * medium. Measured in GeV/fm^3.
   */
  Jetscape::real pi[10];

  /**
   * @brief Bulk viscous pressure.
   * @details Bulk viscous pressure in the cell, measured in GeV/fm^3.
   */
  Jetscape::real bulk_Pi;

  /**
   * @brief Converts the SurfaceCellInfo object to a string representation.
   * @details This method formats the data contained in the SurfaceCellInfo
   * object into a human-readable string for easy output and logging.
   * @return A string representation of the SurfaceCellInfo object.
   */
  std::string sfi_to_string() const {
    std::string str =
        "tau = " + std::to_string(tau) + ", x = " + std::to_string(x) +
        ", y = " + std::to_string(y) + ", eta = " + std::to_string(eta) + "\n";
    str += "d3sigma_mu = " + std::to_string(d3sigma_mu[0]) + ", " +
           std::to_string(d3sigma_mu[1]) + ", " +
           std::to_string(d3sigma_mu[2]) + ", " +
           std::to_string(d3sigma_mu[3]) + "\n";
    str += "energy_density = " + std::to_string(energy_density) +
           ", entropy_density = " + std::to_string(entropy_density) +
           ", temperature = " + std::to_string(temperature) +
           ", pressure = " + std::to_string(pressure) +
           ", qgp_fraction = " + std::to_string(qgp_fraction) +
           ", mu_B = " + std::to_string(mu_B) +
           ", mu_Q = " + std::to_string(mu_Q) +
           ", mu_S = " + std::to_string(mu_S) + "\n";
    str += "umu = " + std::to_string(umu[0]) + ", " + std::to_string(umu[1]) +
           ", " + std::to_string(umu[2]) + ", " + std::to_string(umu[3]) + "\n";
    str += "pi = " + std::to_string(pi[0]) + ", " + std::to_string(pi[1]) +
           ", " + std::to_string(pi[2]) + ", " + std::to_string(pi[3]) + ", " +
           std::to_string(pi[4]) + ", " + std::to_string(pi[5]) + ", " +
           std::to_string(pi[6]) + ", " + std::to_string(pi[7]) + ", " +
           std::to_string(pi[8]) + ", " + std::to_string(pi[9]) + "\n";
    str += "bulk_Pi = " + std::to_string(bulk_Pi) + "\n";
    return str;
  }

  /**
   * @brief Equality operator for SurfaceCellInfo.
   * @details Compares two SurfaceCellInfo objects for equality based on their
   * member variables.
   * @param rhs The right-hand side SurfaceCellInfo object to compare with.
   * @return True if all member variables are equal, false otherwise.
   */
  bool operator==(const SurfaceCellInfo &rhs) const {
    return (tau == rhs.tau && x == rhs.x && y == rhs.y && eta == rhs.eta &&
            d3sigma_mu[0] == rhs.d3sigma_mu[0] &&
            d3sigma_mu[1] == rhs.d3sigma_mu[1] &&
            d3sigma_mu[2] == rhs.d3sigma_mu[2] &&
            d3sigma_mu[3] == rhs.d3sigma_mu[3] &&
            energy_density == rhs.energy_density &&
            entropy_density == rhs.entropy_density &&
            temperature == rhs.temperature && pressure == rhs.pressure &&
            qgp_fraction == rhs.qgp_fraction && mu_B == rhs.mu_B &&
            mu_Q == rhs.mu_Q && mu_S == rhs.mu_S && umu[0] == rhs.umu[0] &&
            umu[1] == rhs.umu[1] && umu[2] == rhs.umu[2] &&
            umu[3] == rhs.umu[3] && pi[0] == rhs.pi[0] && pi[1] == rhs.pi[1] &&
            pi[2] == rhs.pi[2] && pi[3] == rhs.pi[3] && pi[4] == rhs.pi[4] &&
            pi[5] == rhs.pi[5] && pi[6] == rhs.pi[6] && pi[7] == rhs.pi[7] &&
            pi[8] == rhs.pi[8] && pi[9] == rhs.pi[9] && bulk_Pi == rhs.bulk_Pi);
  }

  /**
   * @brief Default constructor.
   * @details Initializes the object with default values.
   */
  SurfaceCellInfo() = default;

  /**
   * @brief Destructor.
   * @details Cleans up any allocated resources (if any) when the object is
   * destroyed.
   */
  ~SurfaceCellInfo(){};
};

}  // namespace Jetscape

#endif  // SURFACECELLINFO_H
