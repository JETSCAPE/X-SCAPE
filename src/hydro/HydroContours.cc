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

#include "HydroContours.h"

#include <dirent.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <utility>

#include "JetScapeLogger.h"
#include "MakeUniqueHelper.h"

RegisterJetScapeModule<HydroContours> HydroContours::reg(
    "CustomModuleHydroContours");

namespace {

// files <prefix><tau> in dir, sorted by tau
std::vector<std::pair<double, std::string>> ListSlices(const std::string &dir,
                                                       const std::string &prefix) {
  std::vector<std::pair<double, std::string>> files;
  if (DIR *d = opendir(dir.c_str())) {
    while (dirent *e = readdir(d)) {
      const std::string name = e->d_name;
      if (name.compare(0, prefix.size(), prefix) != 0) continue;
      char *end = nullptr;
      const double tau = std::strtod(name.c_str() + prefix.size(), &end);
      if (end == name.c_str() + prefix.size() || *end != '\0') continue;
      files.emplace_back(tau, dir + "/" + name);
    }
    closedir(d);
  }
  std::sort(files.begin(), files.end());
  return files;
}

// rows of numbers (a first line that does not start with numbers is a header)
std::vector<std::vector<double>> ReadRows(const std::string &file, int ncol) {
  std::vector<std::vector<double>> rows;
  std::ifstream in(file);
  std::string line;
  while (std::getline(in, line)) {
    for (char &c : line)
      if (c == ',') c = ' ';
    std::istringstream is(line);
    std::vector<double> r(ncol);
    bool ok = true;
    for (int k = 0; k < ncol && ok; ++k) ok = static_cast<bool>(is >> r[k]);
    if (ok) rows.push_back(r);
  }
  return rows;
}

}  // namespace

HydroContours::HydroContours() : FluidDynamics() {
  hydro_status = NOT_START;
  SetId("HydroContours");
}

HydroContours::~HydroContours() {}

void HydroContours::Init() {
  JetScapeModuleBase::InitTask();
  JSINFO << "Initialize FluidDynamics : " << GetId() << " ...";
  InitializeHydro(GetParameterList());
  InitTask();
  InitTasks();
}

void HydroContours::InitTask() {
  auto *e = GetXMLElement({"Hydro", "CustomModuleHydroContours", "dir"});
  dir_ = (e && e->GetText()) ? e->GetText() : "medium_contours";
  dir_.erase(0, dir_.find_first_not_of(" \t\n"));
  dir_.erase(dir_.find_last_not_of(" \t\n") + 1);
  if (auto *t = GetXMLElement({"Hydro", "CustomModuleHydroContours", "Tc"}))
    t->QueryDoubleText(&tc_);
  if (auto *t = GetXMLElement({"Hydro", "CustomModuleHydroContours", "etamax"}))
    t->QueryDoubleText(&etamax_);

  const auto tfiles = ListSlices(dir_, "Tcontour");
  const auto vfiles = ListSlices(dir_, "Vcontour");
  if (tfiles.size() < 2 || tfiles.size() != vfiles.size()) {
    JSWARN << "HydroContours: need matching Tcontour/Vcontour files (at least 2) in "
           << dir_;
    exit(-1);
  }
  // grid from the first file: x outer, y inner, uniform and increasing
  const auto first = ReadRows(tfiles[0].second, 3);
  for (const auto &r : first) {
    if (xs_.empty() || r[0] != xs_.back()) xs_.push_back(r[0]);
    if (xs_.size() == 1) ys_.push_back(r[1]);
  }
  nx_ = static_cast<int>(xs_.size());
  ny_ = static_cast<int>(ys_.size());
  const int nt = static_cast<int>(tfiles.size());
  for (const auto &f : tfiles) taus_.push_back(f.first);
  temp_.assign(static_cast<size_t>(nt) * nx_ * ny_, 0.0);
  ux_.assign(temp_.size(), 0.0);
  uy_.assign(temp_.size(), 0.0);
  for (int it = 0; it < nt; ++it) {
    const auto rt = ReadRows(tfiles[it].second, 3);
    const auto rv = ReadRows(vfiles[it].second, 4);
    if (static_cast<int>(rt.size()) != nx_ * ny_ ||
        static_cast<int>(rv.size()) != nx_ * ny_) {
      JSWARN << "HydroContours: wrong number of rows in " << tfiles[it].second
             << " or " << vfiles[it].second;
      exit(-1);
    }
    for (size_t k = 0; k < rt.size(); ++k) {
      temp_[static_cast<size_t>(it) * nx_ * ny_ + k] = rt[k][2];
      ux_[static_cast<size_t>(it) * nx_ * ny_ + k] = rv[k][2];
      uy_[static_cast<size_t>(it) * nx_ * ny_ + k] = rv[k][3];
    }
  }
  // as medium-2D.f (readvelocities): no flow at grid points with T <= T_c
  double tmax = 0.0;
  for (int it = 0; it < nt; ++it)
    for (int ix = 0; ix < nx_; ++ix)
      for (int iy = 0; iy < ny_; ++iy) {
        double temp = Interpolate(temp_, taus_[it], xs_[ix], ys_[iy]);
        if (temp < tc_) temp = 0.0;
        tmax = std::max(tmax, temp);
        if (!(temp > tc_)) ux_[Index(it, ix, iy)] = uy_[Index(it, ix, iy)] = 0.0;
      }
  // JEWEL-type modules apply the pre-tau0 ramp themselves only when told;
  // here the medium starts at tau = 0 (ramp included), as in medium-2D.f
  hydro_tau_0 = 0.0;
  JSINFO << "HydroContours: " << nt << " time slices " << taus_.front() << " - "
         << taus_.back() << " fm, grid " << nx_ << " x " << ny_ << " ("
         << xs_.front() << " .. " << xs_.back() << " fm), T_max = " << tmax
         << " GeV, T_c = " << tc_ << " GeV, from " << dir_;
}

void HydroContours::InitializeHydro(Parameter parameter_list) {
  hydro_status = INITIALIZED;
}

void HydroContours::EvolveHydro() { hydro_status = FINISHED; }

// lower grid indices as in medium-2D.f (int truncates towards zero)
void HydroContours::Cell(int &it, int &ix, int &iy, double tau, double x,
                         double y) const {
  const double dtau = std::fabs(taus_[1] - taus_[0]);
  const double dx = std::fabs(xs_[1] - xs_[0]), dy = std::fabs(ys_[1] - ys_[0]);
  it = std::max(std::min(static_cast<int>((tau - taus_[0]) / dtau),
                         static_cast<int>(taus_.size()) - 2), 0);
  ix = std::max(std::min(static_cast<int>((x - xs_[0]) / dx), nx_ - 2), 0);
  iy = std::max(std::min(static_cast<int>((y - ys_[0]) / dy), ny_ - 2), 0);
}

double HydroContours::Bilinear(const std::vector<double> &g, int it, double x,
                               double y) const {
  int jt, ix, iy;
  Cell(jt, ix, iy, taus_[it], x, y);
  double yc[2];
  for (int j = 0; j < 2; ++j) {
    const double a = g[Index(it, ix, iy + j)], b = g[Index(it, ix + 1, iy + j)];
    yc[j] = a + (b - a) / (xs_[ix + 1] - xs_[ix]) * (x - xs_[ix]);
  }
  return yc[0] + (yc[1] - yc[0]) / (ys_[iy + 1] - ys_[iy]) * (y - ys_[iy]);
}

double HydroContours::Interpolate(const std::vector<double> &g, double tau,
                                  double x, double y) const {
  int it, ix, iy;
  Cell(it, ix, iy, tau, x, y);
  double yb[2][2];
  for (int j = 0; j < 2; ++j)
    for (int k = 0; k < 2; ++k) {
      const double a = g[Index(it, ix + j, iy + k)];
      const double b = g[Index(it + 1, ix + j, iy + k)];
      yb[j][k] = a + (b - a) / (taus_[it + 1] - taus_[it]) * (tau - taus_[it]);
    }
  double yc[2];
  for (int k = 0; k < 2; ++k)
    yc[k] = yb[0][k] + (yb[1][k] - yb[0][k]) / (xs_[ix + 1] - xs_[ix]) * (x - xs_[ix]);
  return yc[0] + (yc[1] - yc[0]) / (ys_[iy + 1] - ys_[iy]) * (y - ys_[iy]);
}

void HydroContours::GetHydroInfo(
    Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr) {
  fluid_cell_info_ptr = make_unique<FluidCellInfo>();
  FluidCellInfo &c = *fluid_cell_info_ptr;
  c.energy_density = c.entropy_density = c.pressure = 0.0;
  c.temperature = c.vx = c.vy = c.vz = 0.0;
  if (std::fabs(z) >= t) return;
  const double etas = 0.5 * std::log((t + z) / (t - z));
  if (std::fabs(etas) > etamax_) return;
  const double tau = std::sqrt(t * t - z * z);
  if (tau > taus_.back() || x > xs_.back() || y > ys_.back()) return;
  double temp, ux = 0.0, uy = 0.0;
  if (tau <= taus_[0]) {
    temp = Bilinear(temp_, 0, x, y) * tau / taus_[0];
  } else {
    temp = Interpolate(temp_, tau, x, y);
    ux = Interpolate(ux_, tau, x, y);
    uy = Interpolate(uy_, tau, x, y);
  }
  if (temp < tc_) temp = 0.0;
  c.temperature = temp;
  c.qgp_fraction = temp > 0.0 ? 1.0 : 0.0;
  const double ch = std::cosh(etas);
  c.vx = ux / ch;
  c.vy = uy / ch;
  c.vz = std::tanh(etas);
}
