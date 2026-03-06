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

#ifndef TIMEMODULE_H
#define TIMEMODULE_H

#include "ModuleClock.h"
#include "MainClock.h"
//#include "JetScapeTask.h"

#include <string>
#include <memory>

using std::shared_ptr;
using std::string;

namespace Jetscape {

/**
 * @class TimeModule
 * @brief Helper base class providing access to module- and main-level clocks
 *
 * TimeModule centralizes common functionality for modules that operate
 * on a simulation time axis. It supports attaching an optional per-module
 * ModuleClock as well as a shared MainClock. When both are present the
 * module clock is transformed to the main clock before querying time or
 * delta-time information.
 */
class TimeModule  //: public JetScapeTask
{
 public:
  /**
   * @brief Default constructor
   *
   * Constructs a TimeModule without setting an explicit time range.
   * Module clocks are initially null until set with `AddModuleClock()`.
   */
  TimeModule();

  /**
   * @brief Constructor with explicit time range
   * @param t1 Module start time (inclusive)
   * @param t2 Module end time (exclusive)
   */
  TimeModule(double t1, double t2);

  /**
   * @brief Virtual destructor
   */
  virtual ~TimeModule(){};

  /**
   * @brief Print clock information to the logger
   *
   * Logs whether a main clock or module clock is used, prints clock
   * details and the current module time.
   */
  void ClockInfo();

  /**
   * @brief Attach a module-specific clock
   * @param m_mClock Shared pointer to the module's `ModuleClock`
   *
   * Once attached, the module clock will be used in preference to the
   * shared main clock for time queries.
   */
  void AddModuleClock(shared_ptr<ModuleClock> m_mClock) { mClock = m_mClock; }

  /**
   * @brief Return the attached module clock
   * @return Shared pointer to the module's `ModuleClock` or `nullptr` if none
   */
  shared_ptr<ModuleClock> GetModuleClock() const { return mClock; }

  /**
   * @brief Attach the shared main clock
   * @param m_mainClock Shared pointer to the application's `MainClock`
   *
   * Only the first call will set the main clock; subsequent calls are
   * ignored and a warning will be logged.
   */
  void AddMainClock(shared_ptr<MainClock> m_mainClock);

  /**
   * @brief Get the shared main clock
   * @return Shared pointer to the `MainClock` (may be `nullptr`)
   */
  static shared_ptr<MainClock> GetMainClock() { return mainClock; }

  /**
   * @brief Whether any clock (main or module) is in use
   * @return `true` if a clock has been registered, `false` otherwise
   */
  static bool ClockUsed() { return use_clock; }

  /**
   * @brief Whether this module has its own ModuleClock attached
   * @return `true` if a module clock is set, `false` otherwise
   */
  bool UseModuleClock() {
    if (mClock != nullptr)
      return true;
    else
      return false;
  }

  // static bool use_clock; //better in time based module base ...

  /**
   * @brief Get the module's current time in the main clock frame
   * @return Current time for the module
   *
   * If a module clock is attached it will be transformed to the main clock
   * prior to returning the time. If only the main clock is available its
   * current time is returned. If no clock is registered the function logs
   * a warning and terminates the program.
   */
  double GetModuleCurrentTime();

  /**
   * @brief Get the module's delta time (time step)
   * @return Delta time for the module
   *
   * Behavior mirrors `GetModuleCurrentTime()` with respect to module vs
   * main clock selection.
   */
  double GetModuleDeltaT();

  /**
   * @brief Check whether the current module time lies within the module range
   * @return `true` if current time t satisfies t0 <= t < tn
   */
  bool IsValidModuleTime() {
    if (GetModuleCurrentTime() >= t0 && GetModuleCurrentTime() < tn)
      return true;
    else
      return false;
  };

  /**
   * @brief Set the active time range for the module
   * @param t1 Start time (inclusive)
   * @param t2 End time (exclusive)
   */
  void SetTimeRange(double t1, double t2) {
    t0 = t1;
    tn = t2;
  };

  /**
   * @brief Get the module start time
   * @return Start time `t0`
   */
  const double GetTStart() const { return t0; };

  /**
   * @brief Get the module end time
   * @return End time `tn`
   */
  const double GetTEnd() const { return tn; };

 private:
  /**
   * @brief Optional module-specific clock. If non-null the module clock is
   *        transformed into the main clock frame for time queries.
   */
  shared_ptr<ModuleClock> mClock;

  /**
   * @brief Shared application main clock; set via `AddMainClock()`.
   */
  static shared_ptr<MainClock> mainClock;

  double t0;  // module start time; default is 0
  double tn;  // module end time; default is 100

  /**
   * @brief Whether any clock has been registered (module or main)
   */
  static bool use_clock;  // better in time based module base ...
};

}  // end namespace Jetscape

#endif
