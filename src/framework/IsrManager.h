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

#ifndef ISRMANAGER_H
#define ISRMANAGER_H

//#include "JetScapeTask.h"
#include "JetEnergyLossManager.h"
#include "JetClass.h"
#include "sigslot.h"

#include <vector>

namespace Jetscape {
/**
 * @class IsrManager
 * @brief Manager for Initial State Radiation (ISR) modules.
 *
 * Manages ISR modules derived from `JetEnergyLoss` and coordinates the
 * creation of ISR showers, integration with the `HardProcess`, and writing
 * ISR output via `JetScapeWriter` implementations.
 */
class IsrManager : public JetEnergyLossManager
// public std::enable_shared_from_this<IsrManager>
{
 public:
  /** Default constructor to create a jet energy loss manager. Sets task ID as
   * "JLossManager". Flag GetHardPartonListConnected is set to false.
   */
  /**
   * @brief Default constructor.
   *
   * Sets the task identifier to "IsrManager" and configures default
   * verbosity levels.
   */
  IsrManager();

  /** Destructor for the jet energy loss manager.
   */
  /**
   * @brief Virtual destructor.
   *
   * Performs cleanup of owned tasks and resources.
   */
  virtual ~IsrManager();

  /**
   * @brief Initialize the ISR manager and its registered ISR modules.
   *
   * This method initializes all registered ISR tasks/modules. It will
   * connect the manager to the `HardProcess` signal so ISR-produced
   * showers and partons can be updated into the `HardProcess`'s lists.
   *
   * @note If no ISR modules are registered the method will issue a warning
   * and exit with code -1.
   */
  virtual void Init() override;

  /**
   * @brief Execute ISR processing for the current event.
   *
   * Runs each registered ISR module (via `JetEnergyLossManager::Exec()`),
   * collects produced showers and final partons, and updates the
   * `HardProcess` with ISR showers and partons that should continue
   * in the simulation (those with non-negative `pstat`).
   */
  virtual void Exec() override;

  /*
  virtual void Clear();

  virtual void CalculateTime();

  virtual void ExecTime();

  virtual void InitPerEvent();

  virtual void FinishPerEvent();
  */

  /**
   * @brief Write ISR-related output using the provided writer.
   *
   * If the provided `JetScapeWriter` supports ISR ASCII/GZ output, the
   * manager will instruct each ISR module to serialize its shower data.
   *
   * @param w Weak pointer to a `JetScapeWriter` instance used for output.
   */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

 private:
};

}  // end namespace Jetscape

#endif
