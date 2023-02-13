/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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
class BulkMediaBase : public JetScapeModuleBase {
public:
  BulkMediaBase() {
    VERBOSE(8);
    SetId("BulkMediaBase");
  }

  ~BulkMediaBase() {
    VERBOSE(8);
    disconnect_all();
  }

  // Override Init here as function takes care of calling sub-tasks as well
  void Init() override;

  virtual void ExecuteTask();
  virtual void CalculateTime();
  virtual void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,Jetscape::real z,
			   std::unique_ptr<BulkMediaInfo> &bulk_info_ptr){}
protected:

};

} // end namespace Jetscape

#endif // BULKMEDIABASE_H
