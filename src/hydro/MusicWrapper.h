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

using namespace Jetscape;

class HydroSourceJETSCAPE : public HydroSourceBase {
 private:
  std::weak_ptr<LiquefierBase> liquefier_ptr;
  std::weak_ptr<HadronicLiquefier> hadronic_liquefier_ptr;

  double dtau;

 public:
  HydroSourceJETSCAPE() = default;
  ~HydroSourceJETSCAPE() {}

  // set the dtau of the hydro and if the hadronic source terms are present
  // add the value to the hadronic liquefier
  void set_hydro_dtau(double val) {
    dtau = val;
    if (!weak_ptr_is_uninitialized(hadronic_liquefier_ptr)) {
      hadronic_liquefier_ptr.lock()->set_hydro_dtau(dtau);
    }
  };

  void add_a_liquefier(std::shared_ptr<LiquefierBase> new_liquefier) {
    liquefier_ptr = new_liquefier;
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
      liquefier_ptr.lock()->get_source(tau, x, y, eta_s, jmu_tmp);
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
};

#endif  // MUSICWRAPPER_H
