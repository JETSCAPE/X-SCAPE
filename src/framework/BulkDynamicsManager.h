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
#include "BulkMediaBase.h"
#include "JetClass.h"
#include "FluidCellInfo.h"
#include "BulkMediaInfo.h"
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

  void UpdateEnergyDeposit(int t, double edop){ UpdateEnergyDepositFromModules(t, edop); }

  void GetEnergyDensity(int t, double &edensity){ GetEnergyDensityFromModules(t, edensity); }

  void GetHydroCell(double t, double x, double y, double z,
                            std::unique_ptr<FluidCellInfo> &fCell) { GetHydroInfoFromModules(t, x, y, z, fCell); }

  void GetHydroStartTime(double &tau0){GetHydroStartTimeFromModules(tau0); }

  void UpdateEnergyDepositFromModules(int t, double edop);

  void GetEnergyDensityFromModules(int t, double &edensity);

  void GetHydroInfoFromModules(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
			    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);

  void GetHydroStartTimeFromModules(double &tau0);

  void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
                            std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);

  void InfoWrapper(std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr,std::unique_ptr<BulkMediaInfo> &bulk_info_ptr);

  /** Get the new hadrons for the upcoming timesteps and clear the vector for the next timestep
   */
  std::vector<shared_ptr<Hadron>> GetNewHadronsAndClear();

private:

  /** New hadrons for upcoming timestep of transport evolution,
   * to be filled at end of timestep by particlization routine.
   */
  std::vector<shared_ptr<Hadron>> new_hadrons_for_timestep;
  /** Switching temperature between media.
   */
  float Tc;

};

} // end namespace Jetscape

#endif
