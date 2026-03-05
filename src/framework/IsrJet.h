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

#ifndef ISRJET_H
#define ISRJET_H

/**
 * \file IsrJet.h
 * \brief Declaration of the IsrJet class which wraps ISR energy-loss modules.
 */

//#include "JetScapeTask.h"
#include "JetEnergyLoss.h"
#include "sigslot.h"

namespace Jetscape {
/**
 * @class IsrJet
 * @brief Manages initial-state radiation (ISR) jet energy-loss modules.
 *
 * `IsrJet` is a thin manager that groups ISR-related energy loss modules and
 * exposes initialization and lifecycle hooks used by the JETSCAPE framework.
 * It inherits from `JetEnergyLoss` and uses the framework's task infrastructure
 * to discover and initialize contained modules.
 */
class IsrJet : public JetEnergyLoss
// public std::enable_shared_from_this<IsrManager>
{
 public:
  /**
   * @brief Default constructor.
   *
   * Constructs an `IsrJet` manager and sets the task identifier to "IsrJet".
   * It calls the base `JetEnergyLoss` constructor to initialize common state.
   */
  IsrJet();

  /**
   * @brief Virtual destructor.
   *
   * Cleans up resources held by the manager. Most lifetime management is
   * handled by smart pointers in the framework; this destructor is virtual to
   * allow proper cleanup in derived classes.
   */
  virtual ~IsrJet();

  /**
   * @brief Initialize the ISR jet task and contained modules.
   *
   * This method logs basic startup information, verifies that at least one
   * energy-loss module is attached, and prints ISR timing configuration.
   * It is called by the framework when tasks are initialized.
   */
  virtual void InitTask();

 private:
};

}  // end namespace Jetscape

#endif
