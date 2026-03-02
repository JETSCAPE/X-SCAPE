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
// This is a general basic class for hadronic transport

#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "Afterburner.h"

namespace Jetscape {

/**
 * @brief Transport class inherits from the Afterburner class and adds X-SCAPE
 * specific functionalities.
 *
 * Ensures backward compatibility with the JETSCAPE framework.
 */
class Transport : public Afterburner {
 public:
  Transport() {
    VERBOSE(8);
    SetId("Transport");
  }

  ~Transport() {
    VERBOSE(8);
    disconnect_all();
  }

  /**
   * @brief Override Init; this function takes care of calling sub-tasks as well.
   */
  virtual void Init() override;

  /**
   * @brief Performs computations at the end of a time step.
   */
  virtual void ExecuteTask() override;

  /**
   * @brief Performs computations within one time step (evolution).
   */
  virtual void CalculateTime() override;

  /**
   * @brief Get the current list of hadrons in the transport as Jetscape hadrons.
   * @note Must be provided by all Transport implementations.
   * @return std::vector<Hadron> Current hadron list.
   */
  virtual std::vector<Hadron> GetCurrentHadronList() const = 0;
};

}  // end namespace Jetscape

#endif  // TRANSPORT_H