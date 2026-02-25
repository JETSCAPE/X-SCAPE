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
// This is a general info basic class for xscape media modules

#ifndef BULKMEDIAINFO_H
#define BULKMEDIAINFO_H

#include "RealType.h"
#include <string>

namespace Jetscape {

class BulkMediaInfo {
 public:

 /**
   * @brief Container for local bulk-medium / hydrodynamic cell properties.
   *
   * This class stores thermodynamic and flow quantities for a single
   * hydrodynamic cell (or fluid cell) used across the Jetscape framework.
   * All values are initialized to sensible defaults by the constructor.
   */

   // data structure for outputing cell information
  Jetscape::real energy_density;   //!< Local energy density [GeV/fm^3].
  Jetscape::real entropy_density;  //!< Local entropy density [1/fm^3].
  Jetscape::real temperature;      //!< Local temperature [GeV].
  Jetscape::real pressure;         //!< Thermal pressure [GeV/fm^3].
  Jetscape::real qgp_fraction;     //!< Fraction of quark gluon plasma assuming
                                   //!< medium is in QGP+HRG phase.
  Jetscape::real mu_B;             //!< Net baryon chemical potential [GeV].
  Jetscape::real mu_C;             //!< Net charge chemical potential [GeV]
  Jetscape::real mu_S;        //!< Net strangeness chemical potential [GeV].
  Jetscape::real vx, vy, vz;  //!< Flow velocity.
  Jetscape::real pi[4][4];    //!< Shear stress tensor [GeV/fm^3].
  Jetscape::real bulk_Pi;     //!< Bulk viscous pressure [GeV/fm^3].
  Jetscape::real tmn[4][4];   //!< Energy momentum tensor [GeV/fm^3].

  std::string origin_id;  //!< string containing ID of originating bulk media.

  /** Default constructor.*/
  BulkMediaInfo();

  /**
   * @brief Scale this media info in-place by a scalar factor.
   *
   * Multiplies all numeric fields (densities, velocities, tensors, etc.)
   * by the scalar factor @p b.
   *
   * @param b Scalar factor to multiply all fields by.
   * @return Reference to the modified object (*this).
   */
  BulkMediaInfo inline operator*=(Jetscape::real b);

  /** Prints fluid cell properties to the screen. */
  // void Print();
};

// overload +-*/ for easier linear interpolation

/**
 * @brief Component-wise addition of two BulkMediaInfo objects.
 * adds \f$ c = a + b \f$
 *
 * Returns a new BulkMediaInfo equal to the component-wise sum of @p a and
 * @p b. The argument @p a is passed by value and modified, so this
 * implements c = a + b.
 *
 * @param a Left-hand operand (passed by value and modified).
 * @param b Right-hand operand (const reference).
 * @return BulkMediaInfo The component-wise sum.
 */
inline BulkMediaInfo operator+(BulkMediaInfo a, const BulkMediaInfo &b) {
  a.energy_density += b.energy_density;
  a.entropy_density += b.entropy_density;
  a.temperature += b.temperature;
  a.pressure += b.pressure;
  a.qgp_fraction += b.qgp_fraction;
  a.mu_B += b.mu_B;
  a.mu_C += b.mu_C;
  a.mu_S += b.mu_S;
  a.vx += b.vx;
  a.vy += b.vy;
  a.vz += b.vz;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      a.pi[i][j] += b.pi[i][j];
    }
  }
  a.bulk_Pi += b.bulk_Pi;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      a.tmn[i][j] += b.tmn[i][j];
    }
  }
  return a;
}

/**
 * @brief In-place scalar multiplication implementation.
 * Multiply the media info with a scalar factor
 *
 * This implements the body of the in-class declaration of
 * `operator*=`. All numeric members are multiplied by @p b.
 *
 * @param b Scalar factor to multiply by.
 * @return BulkMediaInfo& Reference to *this after scaling.
 */
BulkMediaInfo inline BulkMediaInfo::operator*=(Jetscape::real b) {
  this->energy_density *= b;
  this->entropy_density *= b;
  this->temperature *= b;
  this->pressure *= b;
  this->qgp_fraction *= b;
  this->mu_B *= b;
  this->mu_C *= b;
  this->mu_S *= b;
  this->vx *= b;
  this->vy *= b;
  this->vz *= b;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      this->pi[i][j] *= b;
    }
  }
  this->bulk_Pi *= b;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      this->tmn[i][j] *= b;
    }
  }
  return *this;
}

/**
 * @brief Multiply a BulkMediaInfo by a scalar (scalar on left).
 * multiply \f$ c = a * b \f$
 *
 * Returns a new object equal to @p b scaled by @p a.
 * @param a Scalar factor on the left.
 * @param b Media info to scale (passed by value).
 * @return BulkMediaInfo The scaled media info.
 */
inline BulkMediaInfo operator*(Jetscape::real a, BulkMediaInfo b) {
  b *= a;
  return b;
}

/**
 * @brief Multiply a BulkMediaInfo by a scalar (scalar on right).
 * multiply \f$ c = a * b \f$
 *
 * Returns a new object equal to @p a scaled by @p b.
 * @param a Media info to scale (passed by value).
 * @param b Scalar factor on the right.
 * @return BulkMediaInfo The scaled media info.
 */
inline BulkMediaInfo operator*(BulkMediaInfo a, Jetscape::real b) {
  a *= b;
  return a;
}

/**
 * @brief Divide a BulkMediaInfo by a scalar.
 * division \f$ c = a / b \f$
 *
 * Returns a new object equal to @p a multiplied by (1/b).
 * @param a Media info to scale (passed by value).
 * @param b Scalar divisor.
 * @return BulkMediaInfo The scaled media info.
 */
inline BulkMediaInfo operator/(BulkMediaInfo a, Jetscape::real b) {
  a *= 1.0 / b;
  return a;
}

}  // end namespace Jetscape

#endif  // BulkMediaInfo
