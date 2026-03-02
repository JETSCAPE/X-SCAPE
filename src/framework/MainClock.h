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

// Remark JP: Just basic functionalities, has to be extended and potentially
// made a bit smarter ... Currently assuming uniform ... but should be able to
// be overwritten via Next and ++ operator Include reading parameters from XML
// ... (to be done ...)

#ifndef MAINCLOCK_H
#define MAINCLOCK_H

#include "ClockBase.h"
#include <string>
#include <memory>

using std::string;

namespace Jetscape {

/**
 * @class MainClock
 * @brief Simple main simulation clock derived from `ClockBase`.
 *
 * `MainClock` provides a uniform time stepping clock with a configurable
 * start time, end time and time step (deltaT). It supports advancing the
 * clock via `operator++()` and `Next()` and provides convenience accessors
 * and mutators for the time parameters.
 */
class MainClock : public ClockBase {
 public:
  /**
   * @brief Default constructor.
   *
   * Initializes the clock with default values:
   * - deltaT = 0.1
   * - startTime = 0.0
   * - endTime = 20.0 + deltaT
   */
  MainClock();

  /**
   * @brief Parameterized constructor.
   * @param m_id Identifier for the time reference frame (passed to base).
   * @param m_st Start time for the clock.
   * @param m_et End time for the clock (stored as m_et + deltaT internally).
   * @param m_dt Time step (delta t) between ticks.
   */
  MainClock(string m_id, double m_st, double m_et, double m_dt);

  /**
   * @brief Virtual destructor.
   */
  virtual ~MainClock(){};

  /**
   * @brief Reset clock to configured start time.
   *
   * Sets `currentTime` to `startTime`.
   */
  void Reset() { currentTime = startTime; }

  /**
   * @brief Reset clock to an explicit time.
   * @param m_resetTime Time to set as the current time.
   */
  void ResetToTime(double m_resetTime) { currentTime = m_resetTime; }

  /**
   * @brief Advance the clock by one time step (`deltaT`).
   * @note check usage ... maybe better to always use `Next()` function.
   * @return Reference to this `MainClock` after increment.
   */
  virtual MainClock& operator++();

  /**
   * @brief Advance the clock and indicate if the clock is still within
   * the configured end time.
   *
   * This increments `currentTime` by `deltaT` and returns `true` if
   * `currentTime` is still less than `endTime` after the increment,
   * otherwise returns `false`.
   *
   * @note check if max time check usage sufficient ...
   * @return `true` if more ticks are available, `false` otherwise.
   */
  virtual bool Next();

  /**
   * @brief Alias for `Next()`.
   * @return See `Next()`.
   */
  virtual bool Tick() { return Next(); }

  /**
   * @brief Print informational state of the clock.
   *
   * Calls `ClockBase::Info()` then logs start, end and current time.
   */
  void Info();

  /**
   * @brief Set the start time.
   * @param m_StartTime New start time.
   */
  void SetStartTime(double m_StartTime) { startTime = m_StartTime; }

  /**
   * @brief Set the end time.
   * @param m_EndTime New end time.
   */
  void SetEndTime(double m_EndTime) { endTime = m_EndTime; }

  /**
   * @brief Set the time step (delta t).
   * @param m_deltaT New time step.
   */
  void SetDeltaT(double m_deltaT) { deltaT = m_deltaT; }

  /** @brief Get the configured start time. */
  inline double GetStartTime() { return startTime; }

  /** @brief Get the configured end time. */
  inline double GetEndTime() { return endTime; }

  /** @brief Get the configured time step (deltaT). */
  inline double GetDeltaT() { return deltaT; }

  /** @brief Get the current time of the clock. */
  inline double GetCurrentTime() { return currentTime; }

 private:
  /** @brief Configured start time for the clock. */
  double startTime;

  /** @brief Configured end time for the clock. */
  double endTime;

  /** @brief Time step (delta t) used to advance the clock. */
  double deltaT;

  /** @brief Current time of the clock. */
  double currentTime;
};

}  // end namespace Jetscape

#endif