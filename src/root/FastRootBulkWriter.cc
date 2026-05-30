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
// Fast ROOT writer for the MUSIC hydro evolution (see FastRootBulkWriter.h).
// Reads MUSIC's native in-memory store directly, bypassing the framework AoS.
// -----------------------------------------

#include "FastRootBulkWriter.h"

#include <cmath>
#include <memory>
#include <string>

#include "JetScape.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include "MusicWrapper.h"

#include "TParameter.h"
#include "TString.h"

using namespace std;
using namespace Jetscape;

RegisterJetScapeModule<FastRootBulkWriter> FastRootBulkWriter::reg("FastRootBulkWriter");

FastRootBulkWriter::FastRootBulkWriter() { SetId("FastRootBulkWriter"); }

void FastRootBulkWriter::Init() {
  JSINFO << "Initializing FastRootBulkWriter ...";
  out_file_name = GetXMLElementText({"FastRootBulkWriter", "out_file_name"});

  std::string mode = GetXMLElementText({"FastRootBulkWriter", "grid_mode"}, false);
  if (mode == "native" || mode == "grid")
    grid_mode = mode;
  tau_stride = GetXMLElementInt({"FastRootBulkWriter", "tau_stride"}, false);
  if (tau_stride < 1)
    tau_stride = 1;
  include_preeq = GetXMLElementInt({"FastRootBulkWriter", "include_preequilibrium"}, false);
  if (include_preeq)
    JSWARN << "FastRootBulkWriter: include_preequilibrium is not yet implemented; "
              "dumping the hydro stage only.";

  // user grid (grid mode only)
  x_min = GetXMLElementDouble({"FastRootBulkWriter", "x_min"}, false);
  dx = GetXMLElementDouble({"FastRootBulkWriter", "dx"}, false);
  y_min = GetXMLElementDouble({"FastRootBulkWriter", "y_min"}, false);
  dy = GetXMLElementDouble({"FastRootBulkWriter", "dy"}, false);
  eta_min = GetXMLElementDouble({"FastRootBulkWriter", "eta_min"}, false);
  deta = GetXMLElementDouble({"FastRootBulkWriter", "deta"}, false);
  tau_min = GetXMLElementDouble({"FastRootBulkWriter", "tau_min"}, false);
  dtau = GetXMLElementDouble({"FastRootBulkWriter", "dtau"}, false);
  ntau = GetXMLElementInt({"FastRootBulkWriter", "ntau"}, false);

  JSINFO << "FastRootBulkWriter: out_file_name = " << out_file_name
         << ", grid_mode = " << grid_mode << ", tau_stride = " << tau_stride;
}

void FastRootBulkWriter::init_tree() {
  isinit = true;
  f = new TFile(out_file_name.c_str(), "RECREATE");
  t = new TTree("t", "Tree");
  // One basket per event so a reader can load a single event without reading all.
  t->SetAutoFlush(1);

  t->Branch("user_res", &v_data);
  t->Branch("ntau_freezeout", &ntau_written, "ntau_freezeout/I");

  // Self-describing metadata (mirrors RootBulkWriter param names so existing
  // uproot/FNO readers can reshape: [ntau_freezeout, nx, ny, neta, nFeatures]).
  f->cd();
  for (auto p : vector<std::tuple<bool, string, float>>{
           {true, "nFeatures", nFeatures}, {true, "nx", nx},
           {false, "x_min", x_min},        {false, "dx", dx},
           {true, "ny", ny},               {false, "y_min", y_min},
           {false, "dy", dy},              {true, "neta", neta},
           {false, "eta_min", eta_min},    {false, "deta", deta},
           {true, "ntau", ntau},           {false, "tau_min", eff_tau_min},
           {false, "dtau", eff_dtau},      {true, "tau_stride", tau_stride}}) {
    if (std::get<0>(p)) {
      TParameter<int> pi(std::get<1>(p).c_str(), (int)std::get<2>(p));
      pi.Write();
    } else {
      TParameter<float> pf(std::get<1>(p).c_str(), std::get<2>(p));
      pf.Write();
    }
  }
  TNamed gm("grid_mode", grid_mode.c_str());
  gm.Write();

  JSINFO << "FastRootBulkWriter tree initialized: grid_mode = " << grid_mode
         << ", nx = " << nx << ", ny = " << ny << ", neta = " << neta
         << ", tau_min = " << eff_tau_min << ", dtau = " << eff_dtau;
}

void FastRootBulkWriter::fill_native(::MpiMusic *music, const EvolutionHistory &g) {
  const int nxn = g.nx, nyn = g.ny, netan = g.neta;
  const long n_per_step = (long)nxn * nyn * netan;
  const int num_cells = music->get_number_of_fluid_cells();
  const int ntau_native = (n_per_step > 0) ? (int)(num_cells / n_per_step) : 0;

  nx = nxn;
  ny = nyn;
  neta = netan;
  eff_tau_min = g.tau_min;
  eff_dtau = g.dtau * tau_stride;

  v_data.clear();
  v_data.reserve((size_t)((ntau_native + tau_stride - 1) / tau_stride) * n_per_step *
                 nFeatures);

  FluidCellInfo cell;
  int written = 0;
  for (int it = 0; it < ntau_native; it += tau_stride) {
    const long base = (long)it * n_per_step;
    for (long ic = 0; ic < n_per_step; ic++) {
      music->get_native_fluid_cell((int)(base + ic), cell);
      v_data.push_back((float)cell.energy_density);
      v_data.push_back((float)cell.vx);
      v_data.push_back((float)cell.vy);
      v_data.push_back((float)cell.vz);
    }
    written++;
  }
  ntau_written = written;
}

void FastRootBulkWriter::fill_grid(::MpiMusic *music, const EvolutionHistory &g) {
  // Build a transient working history from the native store (one pass, freed with
  // the module). This avoids the *persistent* framework copy and reuses the exact
  // EvolutionHistory::get() interpolation, so output matches RootBulkWriter.
  const long n_per_step = (long)g.nx * g.ny * g.neta;
  const int num_cells = music->get_number_of_fluid_cells();
  const int ntau_native = (n_per_step > 0) ? (int)(num_cells / n_per_step) : 0;

  work_hist_.tau_min = g.tau_min;
  work_hist_.dtau = g.dtau;
  work_hist_.x_min = g.x_min;
  work_hist_.dx = g.dx;
  work_hist_.y_min = g.y_min;
  work_hist_.dy = g.dy;
  work_hist_.eta_min = g.eta_min;
  work_hist_.deta = g.deta;
  work_hist_.nx = g.nx;
  work_hist_.ny = g.ny;
  work_hist_.neta = g.neta;
  work_hist_.ntau = ntau_native;
  work_hist_.boost_invariant = g.boost_invariant;
  work_hist_.tau_eta_is_tz = false;
  work_hist_.data.resize(num_cells);
  for (int i = 0; i < num_cells; i++)
    music->get_native_fluid_cell(i, work_hist_.data[i]);

  // Output grid: user value, or MUSIC's where the user left it at 0.
  if (!dx) dx = g.dx;
  if (!dy) dy = g.dy;
  if (!deta) deta = g.deta;
  if (!dtau) dtau = g.dtau;
  if (!x_min) x_min = g.x_min;
  if (!y_min) y_min = g.y_min;
  if (!eta_min) eta_min = g.eta_min;
  if (!tau_min) tau_min = g.tau_min;

  nx = 2 * int(fabs(x_min) / dx) + 1;
  ny = 2 * int(fabs(y_min) / dy) + 1;
  neta = 2 * int(fabs(eta_min) / deta) + 1;
  eff_tau_min = tau_min;
  eff_dtau = dtau;

  const double tau_max = g.tau_min + (ntau_native - 1) * g.dtau;
  int n_out_tau = ntau;
  if (n_out_tau <= 0)
    n_out_tau = (dtau > 0) ? (int)((tau_max - tau_min) / dtau) + 1 : 0;

  v_data.clear();
  v_data.reserve((size_t)n_out_tau * nx * ny * neta * nFeatures);
  int written = 0;
  for (int k = 0; k < n_out_tau; k++) {
    double tau_In = tau_min + k * dtau;
    for (int ix = 0; ix < nx; ix++) {
      double x_In = x_min + ix * dx;
      for (int iy = 0; iy < ny; iy++) {
        double y_In = y_min + iy * dy;
        for (int ieta = 0; ieta < neta; ieta++) {
          double eta_In = eta_min + ieta * deta;
          auto c = work_hist_.get(tau_In, x_In, y_In, eta_In);
          v_data.push_back((float)c.energy_density);
          v_data.push_back((float)c.vx);
          v_data.push_back((float)c.vy);
          v_data.push_back((float)c.vz);
        }
      }
    }
    written++;
  }
  ntau_written = written;
  work_hist_.data.clear();
  work_hist_.data.shrink_to_fit();
}

void FastRootBulkWriter::Exec() {
  auto hydro_wp = JetScapeSignalManager::Instance()->GetHydroPointer();
  auto hydro = hydro_wp.lock();
  if (!hydro) {
    JSWARN << "FastRootBulkWriter: no hydro pointer found. Skipping.";
    return;
  }

  auto music = dynamic_pointer_cast<::MpiMusic>(hydro);
  if (!music) {
    JSWARN << "FastRootBulkWriter requires the MUSIC hydro module. Skipping.";
    return;
  }
  if (!music->get_dump_hydro_only()) {
    JSWARN << "FastRootBulkWriter requires <Hydro><MUSIC><dump_hydro_only>1 so "
              "MUSIC keeps its native evolution store. Skipping.";
    return;
  }

  const auto &g = hydro->get_bulk_info();  // grid metadata (data empty here)
  if (music->get_number_of_fluid_cells() <= 0 || g.nx <= 0 || g.neta <= 0) {
    JSWARN << "FastRootBulkWriter: MUSIC native store is empty "
              "(need output_evolution_to_memory=1). Skipping.";
    return;
  }

  if (grid_mode == "native")
    fill_native(music.get(), g);
  else
    fill_grid(music.get(), g);

  if (!isinit)
    init_tree();
  t->Fill();

  // Release MUSIC's native store now that this event has been written.
  music->clear_hydro_info_from_memory();
}

FastRootBulkWriter::~FastRootBulkWriter() {
  if (f) {
    f->cd();
    if (t)
      t->Write();
    f->Close();
  }
}
