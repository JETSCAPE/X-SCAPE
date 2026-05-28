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

#ifndef JETSCAPEPEREVENT_H
#define JETSCAPEPEREVENT_H

#include "JetScape.h"

#include <memory>
#include <vector>

namespace Jetscape {

class JetScapeWriter;

/**
 * @class JetScapePerEvent
 * @brief JetScape driver that exposes the event loop one event at a time.
 *
 * JetScape::Exec() owns the full internal event loop and releases each module's
 * memory (via ClearTasks()) before the next event, so an external program never
 * gets to read the per-event results. JetScapePerEvent splits that loop so the
 * execution can be steered from an external loop:
 *
 *   - ExecInit()          : one-time pre-loop setup (collect writers, CheckExec,
 *                           snapshot active flags). Idempotent; also called
 *                           lazily by ExecPerEvent().
 *   - ExecPerEvent()      : run a single event, leaving all module data in memory.
 *   - ClearPerEvent()     : release that memory and advance the event counter.
 *
 * The current event number is taken from the (framework-global) event counter
 * JetScapeModuleBase::GetCurrentEvent(), which ClearPerEvent() advances. It
 * starts at 0 by default; SetStartEvent() can shift it (e.g. to resume a run or
 * run a sub-range) before the loop begins.
 *
 * Example usage:
 * @code
 * auto jetscape = make_shared<JetScapePerEvent>();
 * jetscape->Init();
 * // jetscape->SetStartEvent(100);  // optional, default 0
 * for (int i = 0; i < jetscape->GetNumberOfEvents(); i++) {
 *   jetscape->ExecPerEvent();
 *   // ... access module results for the current event here ...
 *   jetscape->ClearPerEvent();
 * }
 * jetscape->Finish();
 * @endcode
 */
class JetScapePerEvent : public JetScape {
 public:
  /** @brief Default constructor. */
  JetScapePerEvent();

  /** @brief Destructor. */
  virtual ~JetScapePerEvent();

  /**
   * @brief One-time setup before the external event loop.
   *
   * Collects active writers, runs CheckExec() on active modules and snapshots
   * the original active flags into taskOrgActiveMap. Idempotent: safe to call
   * once explicitly, or rely on the lazy call inside ExecPerEvent().
   */
  void ExecInit();

  /**
   * @brief Execute a single event.
   *
   * Mirrors the body of JetScape::Exec()'s loop EXCEPT it does NOT clear module
   * memory or increment the event counter, so the module results stay in memory
   * and can be read by an external program before ClearPerEvent() is called.
   *
   * The event number used for progress printout and hydro-reuse bookkeeping is
   * read from JetScapeModuleBase::GetCurrentEvent() (advanced by ClearPerEvent()).
   */
  void ExecPerEvent();

  /**
   * @brief Release per-event memory and advance the event counter.
   *
   * Performs the deferred ClearTasks() and IncrementCurrentEvent(). Call after
   * the external program has read the data for the event run by ExecPerEvent().
   */
  void ClearPerEvent();

  /**
   * @brief Disabled on this class.
   *
   * The inherited all-in-one Exec() loop and the per-event API both fill
   * taskOrgActiveMap independently, so mixing them is unsupported. This override
   * warns and exits to make that misuse impossible.
   */
  void Exec() override;

  /**
   * @brief Set the starting event number for the external loop.
   *
   * Shifts the framework-global event counter (GetCurrentEvent()) so the first
   * ExecPerEvent() runs as this event number. Defaults to 0, in which case
   * behavior is unchanged. Must be called before the first ExecPerEvent()
   * (i.e. before ExecInit() runs); it is applied once during ExecInit().
   *
   * @param m_start_event Starting event number (>= 0).
   */
  void SetStartEvent(int m_start_event) { start_event_ = m_start_event; }

  /**
   * @brief Get the configured starting event number.
   * @return The starting event number (default 0).
   */
  int GetStartEvent() const { return start_event_; }

 protected:
  /** @brief Active writers, collected once and reused across events. */
  std::vector<std::weak_ptr<JetScapeWriter>> vWriter;

  /** @brief Guard so ExecInit() runs exactly once. */
  bool exec_initialized_ = false;

  /** @brief Starting event number applied to the global counter in ExecInit(). */
  int start_event_ = 0;
};

}  // end namespace Jetscape

#endif
