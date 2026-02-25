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

// Remark JP: Think about if this base class is truly necessary ...
// Anyways, keep for now in case changes are needed to current implementation
// idea ... Like putting more clock functions in to base class ...

/**
 * @file ClockBase.h
 * @brief Base class for clock objects used by the JETSCAPE framework.
 *
 * The ClockBase class provides a lightweight interface for clocks
 * that can be attached to framework components. It stores an identifier
 * and a reference frame identifier and exposes a minimal API used
 * throughout the framework for querying and reporting time-related
 * information.
 */

#ifndef CLOCKBASE_H
#define CLOCKBASE_H

#include <string>
#include <memory>

using std::string;

namespace Jetscape {

/**
 * @class ClockBase
 * @brief Lightweight base class for time-keeping objects.
 *
 * Derive from ClockBase to implement concrete clock behaviour
 * (e.g., event clocks, simulation clocks). The base class stores
 * an id string and a time reference frame id and provides common
 * accessors used by framework components.
 */
class ClockBase {
 public:

  /**
   * @brief Default constructor initializes id and time reference
   * frame id to empty strings.
   */
  ClockBase();

  /**
   * @brief Virtual destructor for proper cleanup of derived classes.
   */
  virtual ~ClockBase(){};

  /**
   * @brief Print basic information about this clock to the framework logger.
   *
   * The default implementation writes the clock `id` and its time
   * reference frame id using the framework logging facility. Derived
   * classes may override to provide additional details.
   */
  virtual void Info();

  /**
   * @brief Set the clock identifier and time reference frame id.
   */
  void SetId(string m_id) { id = m_id; }

  /**
   * @brief Set the time reference frame identifier for this clock.
   */
  void SetTimeRefFrameId(string m_time_id) { time_id = m_time_id; }
  // void SetCurrentTime(double m_CurrentTime) {currentTime = m_CurrentTime;}

  /**
   * @brief Get the clock identifier.
   * @return The id string previously set via `SetId()`.
   */
  const string GetId() const { return id; }

  /**
   * @brief Get the time reference frame identifier.
   * @return The time reference frame id previously set via
   *         `SetTimeRefFrameId()`.
   */
  const string GetTimeRefFrameId() const { return time_id; }

  /**
   * @brief Return the current time reported by the clock.
   *
   * The base implementation returns a sentinel value (-99.). Concrete
   * clocks should override this method to provide a meaningful time
   * value in the units appropriate for that clock.
   *
   * @return Current time as a double, or -99. if not implemented.
   */
  virtual double GetCurrentTime() { return -99.; }  // not clear if needed ...

  // static bool ClockUsed() { return use_clock; }

 private:
  /**
   * @brief Unique identifier for the clock instance.
   *
   * Typically set by framework configuration or by the owning component.
   */
  string id;

  /**
   * @brief Identifier for the frame of reference used by this clock.
   *
   * For example this may indicate the simulation frame, laboratory
   * frame, or another logical time frame used by modules.
   */
  string time_id;

  // static bool use_clock; //better in time based module base ...
};

}  // end namespace Jetscape

#endif