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

// This Transport class inherits from the Afterburner class and adds some new
// X-SCAPE specific functionalities. This ensures backward compatibility with
// the JETSCAPE framework
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

  /// Override Init here as function takes care of calling sub-tasks as well
  virtual void Init() override;

  /// Takes care of computations done at the end of a time step
  virtual void ExecuteTask() override;

  /// Takes care of the computations within one time step (evolution)
  virtual void CalculateTime() override;

  /// Get the current list of hadrons in the transport as Jetscape hadrons (has
  /// to be provided by all Transport implementations)
  virtual std::vector<Hadron> GetCurrentHadronList() const = 0;
};

}  // end namespace Jetscape

#endif  // TRANSPORT_H