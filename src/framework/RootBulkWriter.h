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
// JETSCAPE module for soft particlization
// This module will generate Monte-Carlo samples for soft hadrons
// -----------------------------------------
// #ifdef USE_ROOT FIXME

#ifndef ROOTBULKWRITER_H_
#define ROOTBULKWRITER_H_

#include <vector>

#include "JetScapeModuleBase.h"
#include "JetClass.h"
#include "JetScapeWriter.h"
#include "FluidEvolutionHistory.h"

#include "TTree.h"
#include "TFile.h"

namespace Jetscape {

class RootBulkWriter : public JetScapeModuleBase {

public:
  RootBulkWriter();
  ~RootBulkWriter();

  virtual void Init(); // Get values from xml file
  virtual void Show();
  virtual void Exec();
  virtual void Clear() {};

  void init_tree(const EvolutionHistory& bInfo);
  void FIXME_checkdata(std::string msg, const EvolutionHistory& bInfo);

  TFile *f {nullptr};
  TTree *t {nullptr};

  std::unique_ptr<float[]> data; // array of floats for the data, using fixed number of tau-steps -- using fixed size array for reading efficiency
  vector<float> v_data; // using a variable sized vector for variable number of tau-steps

  string in_music_name;

private:
  // state variables
  bool use_vec {false}; // use when maxTau <= 0
  bool isinit {false};

  // internal use
  // xml input:
  float x_min{0}, dx{0}, y_min{0}, dy{0}, tau_min{0}, dtau{0}, eta_min{0}, deta{0};

  string out_file_name {"_.root"};

  // derived:
  int nx{0}, ny{0}, neta{0}, ntau{0};
  int ntotal;

  // branch
  float tau_freezeout {0};
  int ntau_freezeout {0};

  float dtau_MUSIC {0};

  const int nFeatures {4}; // rho, vx, vy, vz

  static RegisterJetScapeModule<RootBulkWriter> reg;
};

} // end namespace Jetscape

#endif // ROOTBULKWRITER_H_
// #endif // USE_ROOT FIXME
