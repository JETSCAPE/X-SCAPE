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

#ifndef MUSICWRAPPER_H
#define MUSICWRAPPER_H

#include <memory>

#include "FluidDynamics.h"
#include "music.h"
#include "hydro_source_base.h"
#include "LiquefierBase.h"
#include "data_struct.h"
#include "JetScapeConstants.h"
#include "MakeUniqueHelper.h"
#include <string>

using namespace Jetscape;

class HydroSourceJETSCAPE : public HydroSourceBase {
 private:
  std::weak_ptr<LiquefierBase> liquefier_ptr;
  std::weak_ptr<HadronicLiquefier> hadronic_liquefier_ptr;

  double dtau;

  // Per-step droplet pruning, see prepare_list_for_current_tau_frame().
  double pruning_step_ = 0.1;  // [fm]; the hydro dtau once set_hydro_dtau() ran
  std::shared_ptr<LiquefierBase> liquefier_step_;  // locked once per step

 public:
  HydroSourceJETSCAPE() = default;
  ~HydroSourceJETSCAPE() {}


  // set the dtau of the hydro and if the hadronic source terms are present
  // add the value to the hadronic liquefier
  void set_hydro_dtau(double val) {
    dtau = val;
    pruning_step_ = val;
    if (!weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      hadronic_liquefier_ptr.lock()->set_hydro_dtau(dtau);
    }
  };

  void add_a_liquefier(std::shared_ptr<LiquefierBase> new_liquefier) {
    liquefier_ptr = new_liquefier;
    liquefier_step_.reset();
  }

  //! Called by MUSIC once per time step, before its Runge-Kutta substeps
  //! (query times tau and tau + dtau). A CausalLiquefier droplet deposits only
  //! in the step containing tau_drop + tau_delay, so the liquefier keeps just
  //! the droplets that can deposit at a query time in
  //! [tau - step, tau + 2 step]; the result is unchanged, because the others
  //! add exactly zero, and queries outside the window still see every droplet.
  //! Also locks the liquefier once per step instead of once per cell.
  void prepare_list_for_current_tau_frame(const double tau_local) {
    if (weak_ptr_is_uninitialized(liquefier_ptr)) return;
    liquefier_step_ = liquefier_ptr.lock();
    if (liquefier_step_) {
      liquefier_step_->prepare_active_droplets(tau_local - pruning_step_,
                                               tau_local + 2. * pruning_step_);
    }
  }

  void add_a_hadronic_liquefier(
      std::shared_ptr<HadronicLiquefier> new_liquefier) {
    hadronic_liquefier_ptr = new_liquefier;
  }

  int get_number_of_sources() const {
    int num_sources = 0;
    if (weak_ptr_is_uninitialized(liquefier_ptr)) {
      num_sources += 0;
    } else {
      num_sources += (liquefier_ptr.lock()->get_dropletlist_size());
    }

    if (weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      num_sources += 0;
    } else {
      num_sources += (hadronic_liquefier_ptr.lock()->get_dropletlist_size());
    }
    return num_sources;
  }

  double get_total_E_of_sources() const {
    double total_E = 0.0;
    if (weak_ptr_is_uninitialized(liquefier_ptr)) {
      total_E += 0.0;
    } else {
      total_E += (liquefier_ptr.lock()->get_dropletlist_total_energy());
    }

    if (weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      total_E += 0.0;
    } else {
      total_E +=
          (hadronic_liquefier_ptr.lock()->get_dropletlist_total_energy());
    }
    return total_E;
  }

  double get_net_baryon_number_of_sources() const {
    double net_baryon_number = 0.0;
    if (weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      net_baryon_number += 0.0;
    } else {
      net_baryon_number +=
          (hadronic_liquefier_ptr.lock()->get_dropletlist_net_baryon_number());
    }
    return net_baryon_number;
  }

  double get_net_electric_charge_of_sources() const {
    double net_electric_charge = 0.0;
    if (weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      net_electric_charge += 0.0;
    } else {
      net_electric_charge += (hadronic_liquefier_ptr.lock()
                                  ->get_dropletlist_net_electric_charge());
    }
    return net_electric_charge;
  }

  double get_net_strangeness_of_sources() const {
    double net_strangeness = 0.0;
    if (weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      net_strangeness += 0.0;
    } else {
      net_strangeness +=
          (hadronic_liquefier_ptr.lock()->get_dropletlist_net_strangeness());
    }
    return net_strangeness;
  }

  //! this function returns the energy source term J^\mu at a given point
  //! (tau, x, y, eta_s)
  void get_hydro_energy_source(const double tau, const double x, const double y,
                               const double eta_s, const FlowVec &u_mu,
                               EnergyFlowVec &j_mu) const {
    j_mu = {0.0};
    if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
      std::array<Jetscape::real, 4> jmu_tmp = {0.0};
      if (liquefier_step_) {
        liquefier_step_->get_source(tau, x, y, eta_s, jmu_tmp);
      } else {
        liquefier_ptr.lock()->get_source(tau, x, y, eta_s, jmu_tmp);
      }
      for (int i = 0; i < 4; i++) {
        j_mu[i] =
            jmu_tmp[i] / hbarC;  // convert the unit from GeV/fm^4 to 1/fm^5
      }
    }

    if (!weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      std::array<double, 4> jmu_tmp = {0.0};
      hadronic_liquefier_ptr.lock()->get_source_energy(tau, x, y, eta_s,
                                                       jmu_tmp);
      for (int i = 0; i < 4; i++) {
        // convert the unit from GeV/fm^4 to 1/fm^5
        j_mu[i] += jmu_tmp[i] / hbarC / dtau;
      }
    }
  }

  //! these functions return the B, Q, S density source terms at a given point
  //! (tau, x, y, eta_s)
  double get_hydro_rhob_source(const double tau, const double x, const double y,
                               const double eta_s, const FlowVec &u_mu) const {
    if (weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      return 0.0;
    }
    return hadronic_liquefier_ptr.lock()->get_source_rhob(tau, x, y, eta_s) /
           dtau;
  }

  double get_hydro_rhoq_source(const double tau, const double x, const double y,
                               const double eta_s, const FlowVec &u_mu) const {
    if (weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      return 0.0;
    }
    return hadronic_liquefier_ptr.lock()->get_source_rhoq(tau, x, y, eta_s) /
           dtau;
  }

  double get_hydro_rhos_source(const double tau, const double x, const double y,
                               const double eta_s, const FlowVec &u_mu) const {
    if (weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      return 0.0;
    }
    return hadronic_liquefier_ptr.lock()->get_source_rhos(tau, x, y, eta_s) /
           dtau;
  }
};

//! this is wrapper class for MUSIC so that it can be used as a external
//! library for the JETSCAPE integrated framework
class MpiMusic : public FluidDynamics {
 private:
  // int mode;            //!< records running mode
  std::unique_ptr<MUSIC> music_hydro_ptr;

  Jetscape::real freezeout_temperature;  //!< [GeV]
  int doCooperFrye;                      //!< flag to run Cooper-Frye freeze-out
                                         //!< for soft particles

  int flag_preEq_output_evo_to_memory;
  int flag_output_evo_to_file;
  int flag_output_evo_to_memory;
  int flag_surface_in_memory;
  bool flag_ensure_MusicWrapper_output;
  bool has_source_terms;
  std::shared_ptr<HydroSourceJETSCAPE> hydro_source_terms_ptr;

  int initialProfile_;

  // Allows the registration of the module so that it is available to be
  // used by the Jetscape framework.
  static RegisterJetScapeModule<MpiMusic> reg;

 public:
  MpiMusic();
  ~MpiMusic();

  void CalculateTime();
  void ExecTime();

  void InitializeHydro(Parameter parameter_list);
  int InitializeHydroEnergyProfile();

  void EvolveHydro();
  void EvolveHydroUpto(const double tauEnd);

  void GetHydroInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                    Jetscape::real z,
                    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);

  void GetHydroInfo_JETSCAPE(
      Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
      std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);
  void GetHydroInfo_MUSIC(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                          Jetscape::real z,
                          std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);

  void SetPreEqGridInfo();
  void SetHydroGridInfo();

  void PassPreEqEvolutionHistoryToFramework();
  void PassHydroEvolutionHistoryToFramework();
  void PassHydroSurfaceToFramework();

  void add_a_liquefier(std::shared_ptr<LiquefierBase> new_liquefier) {
    liquefier_ptr = new_liquefier;
    hydro_source_terms_ptr->add_a_liquefier(liquefier_ptr.lock());
  }

  void add_a_hadronic_liquefier(
      std::shared_ptr<HadronicLiquefier> new_liquefier) {
    hadronic_liquefier_ptr = new_liquefier;
    hydro_source_terms_ptr->add_a_hadronic_liquefier(
        hadronic_liquefier_ptr.lock());
  }

  void GetHyperSurface(Jetscape::real T_cut,
                       SurfaceCellInfo *surface_list_ptr){};
  void collect_freeze_out_surface();

  // ── Python-interop: preserve bulk_info across ClearTasks() ─────────────────
  // When set to true, Clear() skips clear_up_evolution_data() so that
  // bulk_info.data survives the end-of-event cleanup and can be inspected from
  // Python after Exec() returns.
  void set_preserve_bulk_info(bool v) { preserve_bulk_info_ = v; }
  bool get_preserve_bulk_info() const { return preserve_bulk_info_; }

  // ── Fast hydro-only ROOT dump support (used by FastRootBulkWriter) ──────────
  // When dump_hydro_only is set, EvolveHydro keeps MUSIC's native in-memory
  // evolution store and skips building the framework AoS (bulk_info.data).
  bool get_dump_hydro_only() const { return dump_hydro_only_; }
  void set_dump_hydro_only(bool v) { dump_hydro_only_ = v; }

  // Skip exporting MUSIC's freeze-out surface to the framework (and the
  // surface*.dat file collection). MUSIC still finds the surface internally —
  // it is the hydro stop condition — but in a hydro-only dump nothing consumes
  // the exported surface, so the hand-off is wasted work. See <skip_surface>.
  bool get_skip_surface() const { return skip_surface_; }
  void set_skip_surface(bool v) { skip_surface_ = v; }

  // Build MUSIC's freeze-out surface (1) or not (0). Without a surface MUSIC
  // stops on the equivalent max(e) < e_fo test (same stop step; music4gpu
  // only, CPU MUSIC ignores the setting and always builds it). Read from the
  // first <Hydro><MUSIC> block for every instance and overridden by the same
  // tag in the instance's own block; the setter overrides both and applies
  // from the next evolution on. See <freeze_out_surface>.
  bool get_freeze_out_surface() const { return freeze_out_surface_ != 0; }
  void set_freeze_out_surface(bool v) { freeze_out_surface_ = v ? 1 : 0; }

  // True if the last evolution stopped because the freeze-out surface
  // reached the transverse grid boundary (MUSIC's reRunHydro): the stored
  // evolution is then truncated. Reset at the start of every evolution.
  bool get_hit_grid_boundary() const { return hit_grid_boundary_; }

  // Thin pass-throughs to MUSIC's native in-memory store (music_hydro_ptr is
  // private, so the writer reaches it through these). music_hydro_ptr is only
  // created in InitializeHydro(), so guard against calls before that (Python).
  int get_number_of_fluid_cells() {
    return music_hydro_ptr ? music_hydro_ptr->get_number_of_fluid_cells() : 0;
  }
  void clear_hydro_info_from_memory() {
    if (music_hydro_ptr)
      music_hydro_ptr->clear_hydro_info_from_memory();
  }

  // Fill a caller-owned FluidCellInfo straight from MUSIC's native store at flat
  // index idx (tau-major x,y,eta order, same as EvolutionHistory::CellIndex).
  // Reuses MUSIC's u^mu->v conversion + hbarc scaling so values match the
  // ordinary framework path exactly.
  void get_native_fluid_cell(int idx, Jetscape::FluidCellInfo &out) {
    fluidCell fc;
    music_hydro_ptr->get_fluid_cell_with_index(idx, &fc);
    out.energy_density = fc.ed;
    out.entropy_density = fc.sd;
    out.temperature = fc.temperature;
    out.pressure = fc.pressure;
    out.vx = fc.vx;
    out.vy = fc.vy;
    out.vz = fc.vz;
    out.mu_B = 0.0;
    out.mu_C = 0.0;
    out.mu_S = 0.0;
    out.qgp_fraction = 0.0;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        out.pi[i][j] = fc.pi[i][j];
    out.bulk_Pi = fc.bulkPi;
  }

  //! Overrides FluidDynamics::Clear() to honour preserve_bulk_info_.
  void Clear();

  bool update_music_input_parameter(const std::string &filename,
                                    const std::string &key, int new_value);

  // Override the SetHydroStartTime function to set the initial time for MUSIC
  void SetHydroStartTime(double tau0) {
    VERBOSE(3) << "Setting hydro start time in MpiMusic to " << tau0;
    FluidDynamics::SetHydroStartTime(tau0);
    if (music_hydro_ptr) {
      VERBOSE(3) << "Setting Initial_time_tau_0 in MUSIC to " << tau0;
      music_hydro_ptr->set_parameter("Initial_time_tau_0", tau0);
    }
  }

 private:
  bool preserve_bulk_info_ = false;
  // Fast hydro-only ROOT dump (see set_dump_hydro_only); read from
  // <Hydro><MUSIC><dump_hydro_only> in MpiMusic::InitializeHydro.
  bool dump_hydro_only_ = false;
  // Skip the freeze-out surface hand-off / file collection (see
  // set_skip_surface); read from <Hydro><MUSIC><skip_surface>.
  bool skip_surface_ = false;
  bool hit_grid_boundary_ = false;
  int freeze_out_surface_ = 1;
  int reported_freeze_out_surface_ = -1;
  int ReadOwnMusicBlockInt(const char *tag) const;
  void ApplyFreezeOutSurface();
  void WarnIfGridBoundaryHit();
};

#endif  // MUSICWRAPPER_H
