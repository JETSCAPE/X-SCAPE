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
  /// Gather all hadrons from soft particlization and fragmentation
  std::vector<std::vector<std::shared_ptr<Hadron>>> GatherAfterburnerHadrons();
  /// Get the events of soft particlization hadrons
  std::vector<std::vector<std::shared_ptr<Hadron>>> GetSoftParticlizationHadrons();
  /// Get the list of fragmentation hadrons
  std::vector<std::shared_ptr<Hadron>> GetFragmentationHadrons();
  /// Get the list of hadrons for the upcoming timestep from BulkDynamicsManager (will clear the list)
  std::vector<std::shared_ptr<Hadron>> GetTimetepParticlizationHadrons();
};

} // end namespace Jetscape

#endif // AFTERBURNER_H
