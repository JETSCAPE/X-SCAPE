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
// -----------------------------------------
// This is a causal liquefier with the JETSCAPE framework
// -----------------------------------------

#ifndef CAUSALLIQUEFIER_H
#define CAUSALLIQUEFIER_H

#include "LiquefierBase.h"
#include "RealType.h"
#include <array>

#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <algorithm>
#include <string>

// <CausalLiquefier><frame>: get_bjorken_frame()/set_bjorken_frame() (bindings test this)
#define XSCAPE_CAUSAL_LIQUEFIER_FRAME 1
// <frame>fluid</frame>: get_frame()/set_frame(), fluid_e_min, fluid_v_max
#define XSCAPE_CAUSAL_LIQUEFIER_FLUID_FRAME 1

namespace Jetscape {

class CausalLiquefier : public Jetscape::LiquefierBase {
 private:
 public:
  // parameters (to be moved to xml)---------------------------
  double dtau, dx, dy, deta;  // in [fm]
  double tau_delay;           // in [fm]
  double time_relax;          // in [fm]
  double d_diff;              // in [fm]
  double width_delta;         // in [fm]
  //---------------------------
  double c_diff;
  double gamma_relax;
  //! The frame the causal diffusion is solved in (<CausalLiquefier><frame>):
  //!  - BJORKEN (bjorken, the default): the Bjorken flow's at the droplet's eta_s;
  //!  - LAB (lab): the lab frame, as in the paper and X-SCAPE before 2026-11;
  //!  - FLUID (fluid): the rest frame of the background fluid where the droplet
  //!    was made (Droplet::get_vfluid()); a droplet without a usable fluid
  //!    there (e < fluid_e_min_, or not looked up) takes the Bjorken frame.
  enum class Frame { LAB = 0, BJORKEN = 1, FLUID = 2 };
  Frame frame_ = Frame::BJORKEN;
  double fluid_e_min_ = 0.1;   //!< [GeV/fm^3] below it: the Bjorken frame
  double fluid_v_max_ = 0.99;  //!< the fluid speed is capped at it
  double dtau0_ = -1.;     // the initial hydro step (the XML's dtau)
  double dtau_late_ = 0.;  // the step from tau_switch_ on (0: constant)
  double tau_switch_ = -1.;

  CausalLiquefier();
  CausalLiquefier(double dtau, double dx, double dy, double deta);
  ~CausalLiquefier(){};

  void InitializeParameters();

  void smearing_kernel(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                       Jetscape::real eta, const Droplet drop_i,
                       std::array<Jetscape::real, 4> &jmu) const;

  bool droplet_may_contribute(const Droplet &drop_i, double tau_lo,
                              double tau_hi) const override;

  // The source is (n . kernel) x p^mu, deposited in the one step containing
  // tau_d + tau_delay: LiquefierBase can normalize it on the hydro grid.
  bool normalizable_on_grid() const override { return true; }
  double deposit_time(const Droplet &drop_i) const override {
    return drop_i.get_xmu()[0] + tau_delay;
  }
  //! the widest window of the schedule: the exact cut is deposit_window()
  double deposit_half_width() const override {
    return 0.5 * std::max(dtau, std::max(dtau0_, dtau_late_)) + 1e-3;
  }

  // The hydro's step schedule (MUSIC late_delta_tau): dtau is the step in
  // effect (the 1/dtau of the source), dtau0_ the initial one, dtau_late_
  // the one from tau_switch_ on (0: constant).  A query at tau deposits the
  // droplets with tau_d + tau_delay in [tau - step_before/2, tau + step_from/2),
  // so the windows of consecutive queries tile the axis also across the
  // switch, and both Runge-Kutta stages at the switch time see the same
  // droplets.
  void set_hydro_step(double dtau_in) override { dtau = dtau_in; }
  void set_late_hydro_step(double dtau_late, double tau_switch) override {
    dtau_late_ = dtau_late;
    tau_switch_ = tau_switch;
  }
  double query_time(double tau, double dtau_in, int k) const override;
  double step_before(double tau) const;
  double step_from(double tau) const;
  void deposit_window(double tau, double &lo, double &hi) const;

  double dumping(double t) const;

  double kernel_rho(double t, double r) const;
  double rho_smooth(double t, double r) const;
  double rho_delta(double t, double r) const;

  double kernel_j(double t, double r) const;
  double j_smooth(double t, double r) const;
  double j_delta(double t, double r) const;

  double get_t(double tau, double eta) const;
  double get_z(double tau, double eta) const;

  double get_ptau(double px, double pz, double eta) const;
  double get_peta(double px, double pz, double eta) const;

  // kept for code written for the two-frame option: true = BJORKEN, false = LAB
  bool get_bjorken_frame() const { return frame_ == Frame::BJORKEN; }
  void set_bjorken_frame(bool bjorken) {
    frame_ = bjorken ? Frame::BJORKEN : Frame::LAB;
  }
  //! "lab", "bjorken" or "fluid"
  std::string get_frame() const;
  void set_frame(const std::string &name);
  //! The velocity of drop_i's kernel frame relative to the Bjorken flow at the
  //! droplet: the background fluid's (fluid, its speed capped at
  //! fluid_v_max_), else 0 (lab, bjorken, fluid's fallback).
  std::array<double, 3> kernel_frame_velocity(const Droplet &drop_i) const;
  //! whether drop_i's kernel uses the fluid's velocity (fluid frame, usable cell)
  bool uses_fluid_velocity(const Droplet &drop_i) const;

  // For debug
  void set_t_delay(double new_tau_delay);
};

};  // namespace Jetscape

#endif  // CAUSALLIQUEFIER_H
