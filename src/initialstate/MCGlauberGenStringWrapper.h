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
#include "JetScapeSignalManager.h"
#include "JetScapeXML.h"
#include <MakeUniqueHelper.h>
using namespace Jetscape;

class MCGlauberGenStringWrapper : public Jetscape::InitialState {
  // this is second 3DMCGlauber wrapper class to generate strings

public:
  MCGlauberGenStringWrapper();
  ~MCGlauberGenStringWrapper() {}

std::shared_ptr<InitialState> ini; // temporary pointer to initial state
  /** Default Exec() function. It can be overridden by other tasks.
   */
  void Exec();
  void Init();
  std::vector<double> Get_Proj_Remnant();
  std::vector<double> Get_Targ_Remnant();
private:
  // Allows the registration of the module so that it is available to be
  // used by the Jetscape framework.
  static RegisterJetScapeModule<MCGlauberGenStringWrapper> reg;
  std::vector<double> GenHardPartonPosAndMomProj_;
  std::vector<double> GenHardPartonPosAndMomTarg_;
};

#endif  // MCGlauberGenStringWrapper_H
