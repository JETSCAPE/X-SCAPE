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

#ifndef MCGlauberGenStringWrapper_H
#define MCGlauberGenStringWrapper_H

#include <memory>

#include "JetScapeModuleBase.h"
#include "InitialState.h"
#include "JetScapeLogger.h"
#include <MakeUniqueHelper.h>
#include "MCGlauberWrapper.h"
using namespace Jetscape;

class MCGlauberGenStringWrapper : public Jetscape::InitialState {
  // this is second 3DMCGlauber wrapper class to generate strings

public:
  MCGlauberGenStringWrapper();
  ~MCGlauberGenStringWrapper() {}


  /** Default Exec() function. It can be overridden by other tasks.
   */
  void Exec();

private:
  // Allows the registration of the module so that it is available to be
  // used by the Jetscape framework.
  static RegisterJetScapeModule<MCGlauberGenStringWrapper> reg;
  std::shared_ptr<InitialState> ini_MC; // temporary pointer to initial state
};

#endif  // MCGlauberGenStringWrapper_H
