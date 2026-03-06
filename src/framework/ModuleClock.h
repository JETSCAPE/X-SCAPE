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

// REMARK JP: Current module transform for testing, just *2 !!!

/**
 * @file ModuleClock.h
 * @brief Declaration of the ModuleClock class — a simple module-level clock.
 *
 * The ModuleClock provides a module-scoped clock that is transformed from
 * the main simulation clock (MainClock).  This file contains the class
 * declaration and inline accessors used by modules that require a
 * clock with a module-specific scaling/offset.
 */

#ifndef MODULECLOCK_H
#define MODULECLOCK_H

#include "ClockBase.h"
#include "MainClock.h"
#include "RealType.h"
#include <string>
#include <memory>

using Jetscape::real;
using std::string;

namespace Jetscape {

/**
 * @class ModuleClock
 * @brief Module-level clock derived from ClockBase.
 *
 * ModuleClock wraps a scaled view of the {@link MainClock} for use inside
 * modules.  In the current test implementation the module time and
 * delta-time are simply twice the corresponding values from the main
 * clock.  Transform() updates the module's internal times from a
 * supplied weak pointer to the MainClock.
 */
class ModuleClock : public ClockBase {
 public:
  /**
   * @brief Construct a new ModuleClock.
   *
   * Initializes the internal time values to sentinel values (typically
   * negative) so uninitialized usage can be detected.
   */
  ModuleClock();

  /**
   * @brief Virtual destructor.
   */
  virtual ~ModuleClock(){};

  // virtual void Transform(string mainClockRef, double mainClockCurrentTime);

  /**
   * @brief Update the module clock from the main clock.
   *
   * The function locks the provided weak pointer to obtain a shared
   * pointer to the MainClock.  If the MainClock exists, the module's
   * current time and delta-time are updated based on the main clock's
   * values (current implementation multiplies both by 2).  If the
   * MainClock has expired, a warning is logged and the program exits.
   *
   * @param mainClock weak pointer to the {@link MainClock}
   */
  virtual void Transform(std::weak_ptr<MainClock> mainClock);

  /**
   * @brief Print information about this clock to the logger.
   *
   * Overrides ClockBase::Info() to include module-specific time values.
   */
  virtual void Info();

  /**
   * @brief Get the current time of the module clock.
   * @return double current module time (units match MainClock)
   */
  inline double GetCurrentTime() { return currentModuleTime; }

  /**
   * @brief Get the delta-time (time step) for the module clock.
   * @return double module delta-time
   */
  inline double GetDeltaT() { return moduleDeltaT; }

  /**
   * @brief Return the maximum time considered by this clock.
   * @return real current module time
   */
  virtual real getTMax() const { return (currentModuleTime); }

  /**
   * @brief Return the minimum time considered by this clock.
   * @return real current module time
   */
  virtual real getTMin() const { return (currentModuleTime); }

 private:
  /**
   * @brief Current time for this module's clock.
   *
   * This value is derived from the MainClock (in the test transform it is
   * two times the MainClock current time).  Initialized to a sentinel
   * value in the constructor to indicate an uninitialized clock.
   */
  double currentModuleTime;

  /**
   * @brief Time step (delta t) for this module's clock.
   *
   * Derived from the MainClock delta-time and subject to the same
   * transformation factor as currentModuleTime.
   */
  double moduleDeltaT;
};

}  // end namespace Jetscape

#endif
