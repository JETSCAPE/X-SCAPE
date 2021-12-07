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
// This is a general info basic class for xscape media modules

#ifndef BULKMEDIAINFO_H
#define BULKMEDIAINFO_H

#include "RealType.h"

namespace Jetscape {

class BulkMediaInfo {
public:
  // data structure for outputing cell information
  Jetscape::real energy_density;  //!< Local energy density [GeV/fm^3].
  Jetscape::real entropy_density; //!< Local entropy density [1/fm^3].
  Jetscape::real temperature;     //!< Local temperature [GeV].
  Jetscape::real pressure;        //!< Thermal pressure [GeV/fm^3].
  Jetscape::real
      qgp_fraction; //!< Fraction of quark gluon plasma assuming medium is in QGP+HRG phase.
  Jetscape::real mu_B;       //!< Net baryon chemical potential [GeV].
  Jetscape::real mu_C;       //!< Net charge chemical potential [GeV]
  Jetscape::real mu_S;       //!< Net strangeness chemical potential [GeV].
  Jetscape::real vx, vy, vz; //!< Flow velocity.
  Jetscape::real pi[4][4];   //!< Shear stress tensor [GeV/fm^3].
  Jetscape::real bulk_Pi;    //!< Bulk viscous pressure [GeV/fm^3].

  /** Default constructor.*/
  BulkMediaInfo();

  /** @param b Multiply the fluid cell by scalar factor b. */
  BulkMediaInfo inline operator*=(Jetscape::real b);

  /** Prints fluid cell properties to the screen. */
  //void Print();
};

//overload +-*/ for easier linear interpolation
/// adds \f$ c = a + b \f$
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
  return a;
}

// Multiply the media info with a scalar factor
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
  return *this;
}

/// multiply \f$ c = a * b \f$
inline BulkMediaInfo operator*(Jetscape::real a, BulkMediaInfo b) {
  b *= a;
  return b;
}

/// multiply \f$ c = a * b \f$
inline BulkMediaInfo operator*(BulkMediaInfo a, Jetscape::real b) {
  a *= b;
  return a;
}

/// division \f$ c = a / b \f$
inline BulkMediaInfo operator/(BulkMediaInfo a, Jetscape::real b) {
  a *= 1.0 / b;
  return a;
}

} // end namespace Jetscape

#endif // BulkMediaInfo
