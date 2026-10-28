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

#ifndef HYDROCONTOURS_H
#define HYDROCONTOURS_H

#include <string>
#include <vector>

#include "FluidDynamics.h"

using namespace Jetscape;

// Boost-invariant (2+1)D medium read from the contour files of the JEWEL
// (2+1)D medium interface (I. Kolbe, arXiv:2303.14166,
// github.com/isobelkolbe/jewel-2.4.0-2D): Tcontour<tau> ("x y T") and
// Vcontour<tau> ("x y v_x v_y" at mid-rapidity), e.g. written by
// CustomModuleMediumDump. Temperature and transverse velocity are interpolated
// exactly as in her medium-2D.f (trilinear in tau, x, y; before the first time
// slice T(tau0, x, y) tau/tau0 and no transverse flow; nothing beyond the
// largest tau, x or y; no medium below T_c or at |eta_s| > etamax; velocity
// zero at grid points below T_c), so that JETSCAPE modules and her JEWEL see
// the same medium. The lab-frame velocity adds Bjorken flow:
// (v_x / cosh eta_s, v_y / cosh eta_s, tanh eta_s).
class HydroContours : public FluidDynamics {
 public:
  HydroContours();
  ~HydroContours();

  // as FluidDynamics::Init, but without requiring a pre-equilibrium module (which
  // FluidDynamics::Init assumes for every module except MUSIC and Brick): the
  // contours already contain the whole evolution
  void Init() override;
  void InitTask();
  void InitializeHydro(Parameter parameter_list);
  void EvolveHydro();
  void GetHydroInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                    Jetscape::real z,
                    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);
  void GetHyperSurface(Jetscape::real T_cut,
                       SurfaceCellInfo *surface_list_ptr) {}

 private:
  // trilinear interpolation of grid g at (tau, x, y), indices as in medium-2D.f
  double Interpolate(const std::vector<double> &g, double tau, double x,
                     double y) const;
  double Bilinear(const std::vector<double> &g, int itau, double x,
                  double y) const;
  void Cell(int &it, int &ix, int &iy, double tau, double x, double y) const;
  size_t Index(int it, int ix, int iy) const {
    return (static_cast<size_t>(it) * nx_ + ix) * ny_ + iy;
  }

  std::string dir_;
  double tc_ = 0.17, etamax_ = 3.1;
  std::vector<double> taus_, xs_, ys_;
  std::vector<double> temp_, ux_, uy_;  // [tau][x][y]
  int nx_ = 0, ny_ = 0;

  static RegisterJetScapeModule<HydroContours> reg;
};

#endif  // HYDROCONTOURS_H
