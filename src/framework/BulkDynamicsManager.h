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

#ifndef BULKDYNAMICSMANAGER_H
#define BULKDYNAMICSMANAGER_H

//#include "JetScapeTask.h"
#include "JetScapeModuleBase.h"
#include "JetClass.h"
#include "sigslot.h"

#include <vector>

namespace Jetscape {
/** @class Bulk dynamics manager manager.
   */
class BulkDynamicsManager
    : public JetScapeModuleBase,
      public std::enable_shared_from_this<BulkDynamicsManager> {

public:
  /** Default constructor to create a bulk dynamics manager. Sets task ID as "BulkDynamicsManager". 
   */
  BulkDynamicsManager();

  /** Destructor for the bulk dynamics manager.
   */
  virtual ~BulkDynamicsManager();

  /** It initializes the tasks attached to the bulk dynamics manager. 
   */
  virtual void Init();

  /** 
  */
  virtual void Exec();

  /** It erases the tasks attached with the bulk dynamics manager. It can be overridden by other tasks.
   */
  virtual void Clear();

  virtual void CalculateTime();

  virtual void ExecTime();

  virtual void InitPerEvent();

  virtual void FinishPerEvent();

private:

  

};

} // end namespace Jetscape

#endif
