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
// -----------------------------------------
// Fast ROOT writer for the MUSIC hydro evolution.
//
// Unlike RootBulkWriter, this module reads MUSIC's *native* in-memory evolution
// store directly (via MpiMusic::get_native_fluid_cell) and never materializes the
// JETSCAPE framework history (bulk_info.data). It is meant for hydro-only
// training-data runs and requires <Hydro><MUSIC><dump_hydro_only>1.
//
// Two grid modes (XML <grid_mode>):
//   * "native" : write MUSIC's own grid, optionally thinned in tau by <tau_stride>.
//                Pure grid-walk, zero interpolation, zero AoS.
//   * "grid"   : write a user-supplied grid; if it matches MUSIC's grid this is a
//                grid-walk, otherwise it interpolates (reusing EvolutionHistory::get
//                so the output matches RootBulkWriter bit-for-bit).
// -----------------------------------------

#ifndef FASTROOTBULKWRITER_H_
#define FASTROOTBULKWRITER_H_

#include <vector>
#include <string>

#include "JetScapeModuleBase.h"
#include "FluidEvolutionHistory.h"

#include "TTree.h"
#include "TFile.h"

// MpiMusic (the MUSIC wrapper) is defined in the GLOBAL namespace in
// MusicWrapper.h, not inside Jetscape. Forward-declare it there so the method
// signatures below bind to the real type rather than a stray Jetscape::MpiMusic.
class MpiMusic;

namespace Jetscape {

class FastRootBulkWriter : public JetScapeModuleBase {

public:
  FastRootBulkWriter();
  ~FastRootBulkWriter();

  virtual void Init();   // read xml configuration
  virtual void Exec();   // dump one event
  virtual void Clear() {}

private:
  void init_tree();
  // Fill v_data by walking MUSIC's native store directly (zero interpolation).
  void fill_native(::MpiMusic *music, const EvolutionHistory &g);
  // Fill v_data on the user grid, interpolating from the native store when the
  // grids differ (reuses EvolutionHistory::get).
  void fill_grid(::MpiMusic *music, const EvolutionHistory &g);

  TFile *f {nullptr};
  TTree *t {nullptr};

  std::vector<float> v_data;  // ragged [tau][x][y][eta][feature] for one event

  // xml input
  std::string out_file_name {"hydro_evo_fast.root"};
  std::string grid_mode {"native"};   // "native" | "grid"
  int tau_stride {1};                 // native mode: write every Nth stored step
  int include_preeq {0};              // reserved (preEq native store is separate)

  // user grid (grid mode only; 0 => default to MUSIC's value)
  float x_min{0}, dx{0}, y_min{0}, dy{0}, eta_min{0}, deta{0}, tau_min{0}, dtau{0};
  int ntau {0};                       // 0 => ragged to end of evolution

  // state / derived (filled at first Exec once the grid is known)
  bool isinit {false};
  int nx{0}, ny{0}, neta{0};
  int ntau_written {0};               // branch: number of tau steps in v_data
  float tau_freezeout {0};            // branch: tau one step past last stored step
  float eff_tau_min {0}, eff_dtau {0};// branch/metadata: actual tau origin & step
  const int nFeatures {4};            // energy_density, vx, vy, vz

  // MUSIC-native grid parameters (mirrors RootBulkWriter _MUSIC keys)
  int nX_MUSIC{0}, nY_MUSIC{0}, neta_MUSIC{0};
  float X_min_MUSIC{0}, dX_MUSIC{0}, Y_min_MUSIC{0}, dY_MUSIC{0};
  float eta_min_MUSIC{0}, deta_MUSIC{0}, tau_min_MUSIC{0}, dtau_MUSIC{0};

  // transient working history for grid-mode interpolation (no framework copy)
  EvolutionHistory work_hist_;

  static RegisterJetScapeModule<FastRootBulkWriter> reg;
};

} // end namespace Jetscape

#endif // FASTROOTBULKWRITER_H_
