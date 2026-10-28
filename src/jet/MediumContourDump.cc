/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#include "MediumContourDump.h"

#include <sys/stat.h>

#include <cmath>
#include <cstdio>
#include <vector>

#include "FluidDynamics.h"
#include "InitialState.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"

RegisterJetScapeModule<MediumContourDump> MediumContourDump::reg(
    "CustomModuleMediumDump");

bool MediumContourDump::done_ = false;

MediumContourDump::MediumContourDump() { SetId("MediumContourDump"); }

MediumContourDump::~MediumContourDump() {}

void MediumContourDump::InitTask() {
  auto text = [this](const char *name, const std::string &def) {
    auto *e = GetXMLElement({"Eloss", "CustomModuleMediumDump", name});
    const char *t = e ? e->GetText() : nullptr;
    return t ? std::string(t) : def;
  };
  auto num = [this](const char *name, double def) {
    auto *e = GetXMLElement({"Eloss", "CustomModuleMediumDump", name});
    double v = def;
    if (e) e->QueryDoubleText(&v);
    return v;
  };
  dir_ = text("dir", "medium_contours");
  dir_.erase(0, dir_.find_first_not_of(" \t\n"));
  dir_.erase(dir_.find_last_not_of(" \t\n") + 1);
  xmax_ = num("xmax", xmax_);
  dx_ = num("dx", dx_);
  tau0_ = num("tau0", tau0_);
  dtau_ = num("dtau", dtau_);
  taumax_ = num("taumax", taumax_);
  ncoll_samples_ = static_cast<long>(num("ncoll_samples", ncoll_samples_));
  contours_ = num("contours", 1.0) != 0.0;
  JSINFO << "MediumContourDump: writing to " << dir_ << ", |x|,|y| <= " << xmax_
         << " fm in steps of " << dx_ << " fm, dtau = " << dtau_ << " fm";
}

void MediumContourDump::DoEnergyLoss(double deltaT, double time, double Q2,
                                     vector<Parton> &pIn,
                                     vector<Parton> &pOut) {
  // the fluid exists once the first parton is propagated; dump it once per run
  if (done_) return;
  done_ = true;
  Dump();
}

void MediumContourDump::Dump() {
  mkdir(dir_.c_str(), 0755);
  const int n = static_cast<int>(std::lround(2.0 * xmax_ / dx_)) + 1;
  std::vector<double> xs(n);
  for (int i = 0; i < n; ++i) xs[i] = -xmax_ + i * dx_;

  if (contours_) {
    double tau0 = tau0_;
    if (tau0 < 0.0) GetHydroTau0Signal(tau0);
    if (tau0 <= 0.0) tau0 = dtau_;
    int nslices = 0;
    for (double tau = tau0; tau <= taumax_ + 1e-9; tau += dtau_) {
      char name[64];
      std::snprintf(name, sizeof(name), "%06.3f", tau);
      std::FILE *ft = std::fopen((dir_ + "/Tcontour" + name).c_str(), "w");
      std::FILE *fv = std::fopen((dir_ + "/Vcontour" + name).c_str(), "w");
      double tmax = 0.0;
      // x outer, y inner, as hydro_analysis writes the contours
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          std::unique_ptr<FluidCellInfo> cell;
          GetHydroCellSignal(tau, xs[i], xs[j], 0.0, cell);
          double temp = 0.0, vx = 0.0, vy = 0.0;
          if (cell) {
            temp = cell->temperature;
            vx = cell->vx;
            vy = cell->vy;
          }
          if (temp <= 0.05) temp = 0.0;  // as hydro_analysis
          tmax = std::max(tmax, temp);
          std::fprintf(ft, "%18.8f   %18.8f   %.8g\n", xs[i], xs[j], temp);
          std::fprintf(fv, "%18.8f   %18.8f   %.8g   %.8g\n", xs[i], xs[j], vx, vy);
        }
      }
      std::fclose(ft);
      std::fclose(fv);
      ++nslices;
      JSINFO << "MediumContourDump: tau = " << tau << " fm, T_max = " << tmax
             << " GeV";
      if (tmax == 0.0) break;  // the fluid has frozen out everywhere
    }
    JSINFO << "MediumContourDump: " << nslices << " time slices written";
  }

  // histogram of the binary-collision positions (NCOLLHISTO), bin centres xs
  auto ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  if (!ini) {
    JSWARN << "MediumContourDump: no initial state, no N_coll histogram";
    return;
  }
  std::vector<long> h(static_cast<size_t>(n) * n, 0);
  for (long s = 0; s < ncoll_samples_; ++s) {
    double t, x, y, z;
    ini->SampleABinaryCollisionPoint(t, x, y, z);
    const int i = static_cast<int>(std::floor((x + xmax_) / dx_ + 0.5));
    const int j = static_cast<int>(std::floor((y + xmax_) / dx_ + 0.5));
    if (i >= 0 && i < n && j >= 0 && j < n) ++h[static_cast<size_t>(i) * n + j];
  }
  std::FILE *fh = std::fopen((dir_ + "/ncoll_histo.dat").c_str(), "w");
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      std::fprintf(fh, "%18.8f   %18.8f   %ld\n", xs[i], xs[j],
                   h[static_cast<size_t>(i) * n + j]);
  std::fclose(fh);
  JSINFO << "MediumContourDump: N_coll histogram from " << ncoll_samples_
         << " samples written";
}
