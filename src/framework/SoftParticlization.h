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
// JETSCAPE module for soft particlization
// This module will generate Monte-Carlo samples for soft hadrons
// -----------------------------------------

#ifndef SOFTPARTICLIZATION_H_
#define SOFTPARTICLIZATION_H_

#include <vector>

#include "JetScapeModuleBase.h"
#include "JetClass.h"
#include "JetScapeWriter.h"
#include "SurfaceCellInfo.h"

namespace Jetscape {

/**
 * @brief JETSCAPE module for soft particlization
 *
 * This module will generate Monte-Carlo samples for soft hadrons from the
 * hydrodynamic output.
 */
class SoftParticlization : public JetScapeModuleBase {
 private:
  /// Flag for the connection status of the GetHydroHyperSurface signal
  bool HydroHyperSurfaceConnected_;

  /// Flag for the connection status of the ClearHydroHyperSurface signal
  bool ClearHydroHyperSurfaceConnected_;

  /// Parameters for building a surface from the stored evolution
  SurfaceFinderParams surface_params_;

  /// One-shot seed override (SetNextRandomSeed) and the seed last used
  long next_random_seed_ = 0;
  bool has_next_random_seed_ = false;
  long last_random_seed_ = 0;

 public:
  /**
   * @brief One event's hadrons as flat arrays, all samples concatenated in
   * order: what SetCompactHadronOutput() fills instead of Hadron_list_.
   */
  struct HadronArrays {
    std::vector<long long> sample_counts;  ///< hadrons per sample
    std::vector<int> pid, pstat;
    std::vector<float> p;     ///< (N, 4): E, px, py, pz [GeV]
    std::vector<float> x;     ///< (N, 4): t, x, y, z [fm]
    std::vector<float> mass;  ///< [GeV]

    /// Empty the arrays and free their memory.
    void clear() { HadronArrays().swap(*this); }
    void swap(HadronArrays &o) {
      sample_counts.swap(o.sample_counts);
      pid.swap(o.pid);
      pstat.swap(o.pstat);
      p.swap(o.p);
      x.swap(o.x);
      mass.swap(o.mass);
    }
  };

 protected:
  /// The last event's hadrons, if the module filled them in compact form
  HadronArrays compact_hadrons_;
  bool has_compact_hadrons_ = false;
  bool compact_hadrons_taken_ = false;

  /**
   * @brief The seed for this event: the SetNextRandomSeed override if one is
   * pending, else a draw from the module's generator. Recorded for
   * GetLastRandomSeed().
   */
  long NextRandomSeed() {
    last_random_seed_ = has_next_random_seed_
                            ? next_random_seed_
                            : static_cast<long>((*GetMt19937Generator())());
    has_next_random_seed_ = false;
    return last_random_seed_;
  }

 public:
  /**
   * @brief Construct a new SoftParticlization object
   */
  SoftParticlization();

  /**
   * @brief Destroy the SoftParticlization object
   */
  ~SoftParticlization();

  /**
   * @brief Initialize the SoftParticlization module
   * @note Override Init (not InitTask) here as sub-tasks are called as well.
   */
  void Init() override;

  /**
   * @brief Execute the SoftParticlization module
   */
  virtual void ExecuteTask();

  /**
   * @brief Clear the SoftParticlization module
   */
  virtual void ClearTask();

  /**
   * @brief Signal for getting the hydrodynamic hypersurface
   */
  sigslot::signal1<std::vector<SurfaceCellInfo> &, multi_threaded_local>
      GetHydroHyperSurface;

  /**
   * @brief Signal for building a surface from the hydro's stored evolution
   * (FluidDynamics::FindSurfaceFromEvolution). Used when the hydro did not
   * provide a surface of its own through GetHydroHyperSurface.
   */
  sigslot::signal2<SurfaceFinderParams, std::vector<SurfaceCellInfo> &,
                   multi_threaded_local>
      FindHydroHyperSurface;

  /**
   * @brief Signal for clearing the hydrodynamic hypersurface
   */
  sigslot::signal0<multi_threaded_local> ClearHydroHyperSurface;

  /**
   * @brief Surface parameters for FindHydroHyperSurface, read in Init() from
   * the optional XML keys SoftParticlization/{T_sw, surface_dtau,
   * surface_dx, surface_deta}.
   */
  const SurfaceFinderParams &GetSurfaceFinderParams() const {
    return surface_params_;
  }

  /**
   * @brief Use this seed for the next event's sampling instead of drawing one
   * from the module's generator (one-shot). This is how a stored surface is
   * re-sampled exactly: pass the seed GetLastRandomSeed() reported in the
   * original run.
   */
  void SetNextRandomSeed(long seed) {
    next_random_seed_ = seed;
    has_next_random_seed_ = true;
  }

  /**
   * @brief The seed the last event's sampling used.
   */
  long GetLastRandomSeed() const { return last_random_seed_; }

  /**
   * @brief Number of samples (oversamples) per event from the next event on.
   * Returns false if this module has no such setting. iSS reads it per event.
   */
  virtual bool SetNumberOfSamples(int n) { return false; }

  /**
   * @brief From the next event on, hand the hadrons over as flat arrays
   * (GetCompactHadrons) instead of one Hadron object per hadron in
   * Hadron_list_: ~8x less memory (44 instead of ~340 bytes per hadron), for
   * callers that only want the numbers (PyJetscape's soft_hadrons_numpy, which
   * takes them with TakeCompactHadrons). Hadron_list_ then stays empty, so
   * writers and afterburners see no bulk hadrons. Returns false if this module
   * cannot do it (it then fills Hadron_list_ as usual).
   */
  virtual bool SetCompactHadronOutput(bool on) { return false; }

  /**
   * @brief True if the last event's hadrons are in GetCompactHadrons() rather
   * than in Hadron_list_.
   */
  bool HasCompactHadrons() const { return has_compact_hadrons_; }

  /**
   * @brief The last event's hadrons in compact form (see HasCompactHadrons);
   * empty once TakeCompactHadrons() moved them out.
   */
  const HadronArrays &GetCompactHadrons() const { return compact_hadrons_; }

  /**
   * @brief Move the last event's compact hadrons into `out` without copying.
   * Returns false if they were already taken this event.
   */
  bool TakeCompactHadrons(HadronArrays &out) {
    if (compact_hadrons_taken_)
      return false;
    HadronArrays().swap(out);
    out.swap(compact_hadrons_);
    compact_hadrons_taken_ = true;
    return true;
  }

  /**
   * @brief Drop the compact hadrons (ClearTask does it every event).
   */
  void ClearCompactHadrons() {
    compact_hadrons_.clear();
    has_compact_hadrons_ = false;
    compact_hadrons_taken_ = false;
  }

  /**
   * @brief Set the GetHydroHyperSurfaceConnected flag
   *
   * @param m_GetHydroHyperSurfaceConnected Boolean flag
   */
  void SetGetHydroHyperSurfaceConnected(bool m_GetHydroHyperSurfaceConnected) {
    HydroHyperSurfaceConnected_ = m_GetHydroHyperSurfaceConnected;
  }

  /**
   * @brief Set the ClearHydroHyperSurfaceConnected flag
   *
   * @param m_ClearHydroHyperSurfaceConnected Boolean flag
   */
  void SetClearHydroHyperSurfaceConnected(
      bool m_ClearHydroHyperSurfaceConnected) {
    ClearHydroHyperSurfaceConnected_ = m_ClearHydroHyperSurfaceConnected;
  }

  /**
   * @brief Get the GetHydroHyperSurfaceConnected flag
   *
   * @return Boolean
   */
  bool GetGetHydroHyperSurfaceConnected() const {
    return HydroHyperSurfaceConnected_;
  }

  /**
   * @brief Get the ClearHydroHyperSurfaceConnected flag
   *
   * @return Boolean
   */
  bool GetClearHydroHyperSurfaceConnected() const {
    return ClearHydroHyperSurfaceConnected_;
  }

  /// List of hadrons
  std::vector<std::vector<shared_ptr<Hadron>>> Hadron_list_;

  /**
   * @brief Clear the hadron list
   */
  void ClearHadronList() { Hadron_list_.clear(); };

  /// Flag for boost invariance
  bool boost_invariance;

  /**
   * @brief Check the boost invariance
   *
   * @return Boolean
   */
  bool check_boost_invariance();
};

}  // end namespace Jetscape

#endif  // SOFTPARTICLIZATION_H_
