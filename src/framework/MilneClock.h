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

#ifndef MILNECLOCK_H
#define MILNECLOCK_H

#include "ModuleClock.h"
#include "MainClock.h"
#include "RealType.h"
#include <memory>

using Jetscape::real;

namespace Jetscape {

/**
 * @file MilneClock.h
 * @brief Module clock that uses Milne (tau-eta) coordinates.
 *
 * The MilneClock adapts timing information from the central `MainClock`
 * into proper-time (tau) and space-rapidity (eta) aware values for
 * modules that operate in Milne coordinates.
 */

class MilneClock : public ModuleClock {
 public:
  /**
   * @brief Construct a new MilneClock
   *
   * Initializes internal timing members to sensible defaults.
   */
  MilneClock();

  /**
   * @brief Virtual destructor
   */
  virtual ~MilneClock(){};

  /**
   * @brief Print basic information about this clock to the logger.
   *
   * This typically reports the clock id and configured tau range.
   */
  void Info();

  /**
   * @brief Transform time information from the `MainClock` to Milne time.
   *
   * This method locks the provided weak pointer to the main clock and
   * extracts the current time and delta time. It then computes the
   * Milne-specific tauMin_ and tauMax_ values based on the configured
   * `etaMax_`.
   *
   * @param mainClock Weak pointer to the `MainClock` instance providing
   *                  global timing information.
   */
  void Transform(std::weak_ptr<MainClock> mainClock);

  /**
   * @brief Set the maximum space-rapidity (eta) used for tau calculations.
   *
   * @param etaMax Maximum pseudorapidity (eta) for which tauMin/tauMax
   *               are computed.
   */
  void setEtaMax(const real etaMax) { etaMax_ = etaMax; }

  /**
   * @brief Get the minimum proper time (tau) available to this module.
   *
   * @return real The minimum proper time (tauMin_)
   */
  real getTMin() const { return (tauMin_); }

  /**
   * @brief Get the maximum proper time (tau) available to this module.
   *
   * @return real The maximum proper time (tauMax_)
   */
  real getTMax() const { return (tauMax_); }

  /**
   * @brief Get the current module time in the module's units (fm/c).
   *
   * @return double Current time propagated from the `MainClock`.
   */
  double GetCurrentTime() { return (currentModuleTime_); }

  /**
   * @brief Get the delta time (time step) used by the module.
   *
   * @return double The module time step propagated from `MainClock`.
   */
  double GetDeltaT() { return (moduleDeltaT_); }

 private:
  /**
   * @brief Current module time (proper time tau) taken from MainClock.
   *
   * Units: fm/c
   */
  real currentModuleTime_;

  /**
   * @brief Time step (delta tau) for the module.
   *
   * Units: fm/c
   */
  real moduleDeltaT_;

  /**
   * @brief Maximum space-rapidity (eta) used when computing tauMin_.
   */
  real etaMax_;

  /**
   * @brief Minimum proper time (tau) relevant for this module.
   *
   * Computed as currentMainTime / cosh(etaMax_).
   */
  real tauMin_;

  /**
   * @brief Maximum proper time (tau) relevant for this module.
   *
   * Typically set to the current main clock time when transformed.
   */
  real tauMax_;
};

}  // end namespace Jetscape

#endif  // MILNECLOCK_H
