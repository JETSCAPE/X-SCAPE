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

#ifndef AFTERBURNER_H
#define AFTERBURNER_H

#include "JetScapeModuleBase.h"
#include "SoftParticlization.h"
#include "HadronizationManager.h"
#include "RealType.h"
#include "BulkMediaInfo.h"
#include "sigslot.h"

namespace Jetscape {

/// Interface to hadronic afterburner
class Afterburner : public JetScapeModuleBase {
public:
  Afterburner() {
    VERBOSE(8);
    SetId("Afterburner");
  }

  ~Afterburner() {
    VERBOSE(8);
    disconnect_all();
  }

  virtual void Init();
  virtual void Exec();
  virtual void CalculateTime();
  virtual void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,Jetscape::real z,
			   std::unique_ptr<BulkMediaInfo> &bulk_info_ptr){}
protected:
  /// Pointer to particlization sampler, which provides initial hadrons
  std::shared_ptr<SoftParticlization> soft_particlization_sampler_;
  /// Pointer to hadronization module, which provides fragmentation hadrons
  std::shared_ptr<HadronizationManager> hard_particlization_module_;

  // TODO(stdnmr) Move to .cc 
  std::vector<shared_ptr<Hadron>> GetFragmentationHadrons() {
    JSINFO << "Get fragmentation hadrons in Afterburner";
    vector<shared_ptr<Hadron>> tmp_list;
    hard_particlization_module_->GetHadrons(tmp_list);
    JSINFO << "Got " << tmp_list.size()
           << " fragmentation hadrons from HadronizationManager.";
    for (auto h : tmp_list) {
      std::cout << *h << std::endl;
    }
    // TODO(stdnmr) Check that hadrons have positions i.e. are from hybrid hadronization,
    // one could also get hadrons without position info here.
    return tmp_list;
  }
};

} // end namespace Jetscape

#endif // AFTERBURNER_H
