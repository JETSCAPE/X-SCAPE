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
// This is a general basic class for hadronic afterburner

#ifndef BULKMEDIABASE_H
#define BULKMEDIABASE_H

#include "JetScapeModuleBase.h"
#include "SoftParticlization.h"
#include "RealType.h"
#include "BulkMediaInfo.h"
#include "sigslot.h"

namespace Jetscape {

/// Interface to hadronic afterburner
/**
 * @class BulkMediaBase
 * @brief General base class for bulk media modules (hadronic afterburner).
 *
 * Provides a common interface for bulk-media-like modules that run as
 * JetScape tasks. Subclasses should override task-related methods to
 * implement initialization, execution, and time-dependent calculations.
 */
class BulkMediaBase : public JetScapeModuleBase {
 public:
  /**
   * @brief Default constructor.
   *
   * Sets the default verbosity level and module id to "BulkMediaBase".
   */
  BulkMediaBase() {
    VERBOSE(8);
    SetId("BulkMediaBase");
  }

  /**
   * @brief Destructor.
   *
   * Performs cleanup and disconnects all signal/slot connections.
   */
  ~BulkMediaBase() {
    VERBOSE(8);
    disconnect_all();
  }

  /**
   * @brief Initialize the module and its sub-tasks.
   *
   * Overrides JetScapeModuleBase::InitTask to ensure that XML parameters are
   * loaded and that any registered sub-tasks are also initialized.
   */
  void Init() override;

  /**
   * @brief Execute the module's primary task.
   *
   * Subclasses should override this method to perform the bulk-media
   * execution behavior.
   */
  virtual void ExecuteTask();

  /**
   * @brief Perform time-based calculations for the module.
   *
   * This method is typically called each time-step; it triggers
   * CalculateTimeTask to perform the actual calculation.
   */
  virtual void CalculateTime();

  /**
   * @brief Query bulk-media information at a spacetime point.
   *
   * @param t Time coordinate at which to get bulk information.
   * @param x Spatial x coordinate.
   * @param y Spatial y coordinate.
   * @param z Spatial z coordinate.
   * @param bulk_info_ptr Output unique_ptr reference that should be filled by
   *                      the implementation with a newly created
   *                      BulkMediaInfo instance describing the medium.
   *
   * Default implementation does nothing; override to provide local medium
   * properties used by e.g. particlization or transport modules.
   */
  virtual void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                           Jetscape::real z,
                           std::unique_ptr<BulkMediaInfo> &bulk_info_ptr) {}

 protected:
};

}  // end namespace Jetscape

#endif  // BULKMEDIABASE_H
