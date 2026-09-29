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

#ifndef LIQUEFIERBASE_H
#define LIQUEFIERBASE_H

#include "JetClass.h"
#include "sigslot.h"
#include "FluidCellInfo.h"

#include <array>
#include <string>
#include <vector>

/// Defined since LiquefierBase normalizes each droplet on the hydro grid
/// (set_hydro_grid(), normalize_active_droplets(), get_droplet_flux()), so
/// code built against older headers can test for it.
#define XSCAPE_LIQUEFIER_GRID_NORMALIZATION 1
#include "RealType.h"

namespace Jetscape {

/**
 * @brief Represents a localized energy-momentum contribution from a parton to
 * the fluid medium.
 *
 * In the JETSCAPE framework, a Droplet is a conceptual object used to bridge
 * the gap between discrete partonic information and continuous hydrodynamic
 * fields. Each droplet represents a localized "chunk" of energy and momentum to
 * be deposited into the hydrodynamic grid.
 */

class Droplet {
 private:
  std::array<Jetscape::real, 4> xmu;  ///< Position 4-vector
  std::array<Jetscape::real, 4> pmu;  ///< Momentum 4-vector

 public:
  /**
   * @brief Default constructor.
   *
   * Constructs an uninitialized droplet. The position and momentum vectors will
   * contain undefined values until explicitly set.
   */
  Droplet() = default;

  /**
   * @brief Construct a droplet from given position and momentum vectors.
   *
   * @param x_in The initial position 4-vector.
   * @param p_in The initial momentum 4-vector.
   */
  Droplet(std::array<Jetscape::real, 4> x_in,
          std::array<Jetscape::real, 4> p_in) {
    xmu = x_in;
    pmu = p_in;
  }

  /**
   * @brief Destructor.
   */
  ~Droplet(){};

  /**
   * @brief Get the position 4-vector of the droplet.
   *
   * @return A copy of the position vector.
   */
  std::array<Jetscape::real, 4> get_xmu() const { return (xmu); }

  /**
   * @brief Get the momentum 4-vector of the droplet.
   *
   * @return A copy of the momentum vector.
   */
  std::array<Jetscape::real, 4> get_pmu() const { return (pmu); }
};

/**
 * @brief Base class for converting partonic energy/momentum into hydrodynamic
 * sources ("liquefying").
 */
/**
 * @brief The hydro's computational grid: cell centres x_min + i dx (i < nx),
 * the same in y and eta.  Set by the hydro that samples get_source(), so the
 * liquefier can normalize each droplet on exactly the cells it is sampled at
 * (see LiquefierBase::normalize_active_droplets()).
 */
struct HydroGrid {
  int nx = 0, ny = 0, neta = 0;
  double x_min = 0., y_min = 0., eta_min = 0.;
  double dx = 0., dy = 0., deta = 0.;
  bool valid() const {
    return nx > 0 && ny > 0 && neta > 1 && dx > 0. && dy > 0. && deta > 0.;
  }
  bool operator==(const HydroGrid &o) const {
    return nx == o.nx && ny == o.ny && neta == o.neta && x_min == o.x_min &&
           y_min == o.y_min && eta_min == o.eta_min && dx == o.dx &&
           dy == o.dy && deta == o.deta;
  }
};

class LiquefierBase {
 private:
  std::vector<Droplet>
      dropletlist;  ///< List of droplets representing source contributions

  // ── normalization on the hydro grid ─────────────────────────────────────
  // The hydro samples get_source() once per cell centre and step, and adds
  // tau dtau dx dy deta J^mu per cell.  For a kernel that is not resolved by
  // the grid that sum is not the droplet's four-momentum (it can be 0 or
  // several times it).  Each droplet is therefore scaled by 1/flux, flux =
  // the sampled sum of its kernel (sampled_flux()), before it deposits.
  HydroGrid hydro_grid_;
  bool normalize_on_grid_ = true;  ///< <Liquefier><normalize_on_hydro_grid>
  std::vector<double> droplet_norm_;  ///< per droplet: 1/flux (1: not yet)
  std::vector<double> droplet_flux_;  ///< per droplet: flux (-1: not yet)
  /// A droplet whose kernel misses every cell centre (flux 0) is deposited
  /// whole into the cell nearest to it, at its deposit query time.
  struct PointDeposit {
    int idx;
    double tau_q, dtau, x, y, eta;
  };
  std::vector<PointDeposit> point_deposits_;
  int n_normalized_ = 0, n_point_ = 0, n_lost_ = 0;
  double E_droplets_ = 0., E_sampled_ = 0.;   // normalized droplets only
  double flux_min_ = 1e300, flux_max_ = 0.;
  bool in_light_cone(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                     Jetscape::real eta, const Droplet &drop_i) const;
  void reset_normalization();
  /// Indices into dropletlist of the droplets that can contribute to a query
  /// time in [active_tau_lo_, active_tau_hi_]; see prepare_active_droplets().
  /// An empty window (lo > hi) means not prepared: get_source() then loops
  /// over every droplet.
  std::vector<int> active_droplets_;
  double active_tau_lo_ = 1.;
  double active_tau_hi_ = 0.;
  bool GetHydroCellSignalConnected;  ///< Flag for whether signal connection to
                                     ///< hydro exists
  const int drop_stat;               ///< Droplet statistics
  const int miss_stat;               ///< Missed parton statistics
  const int neg_stat;                ///< Negative energy statistics
  const Jetscape::real
      hydro_source_abs_err;  ///< Error tolerance for hydro sources
  /// A vertex whose 4-momentum mismatch exceeds this fraction of its incoming
  /// energy is reported; smaller ones only get the p_missing parton.
  const double e_mom_rel_tol;
  bool
      threshold_energy_switch;  ///< Whether to apply energy threshold filtering
  double e_threshold;           ///< Energy threshold value

 public:
  /**
   * @brief Constructor.
   */
  LiquefierBase();

  /**
   * @brief Destructor that clears droplet list.
   */
  ~LiquefierBase() { ClearTask(); }

  /**
   * @brief Add a droplet to the internal list.
   * @param droplet_in Droplet to add
   */
  void add_a_droplet(Droplet droplet_in) {
    dropletlist.push_back(droplet_in);
    droplet_norm_.push_back(1.0);
    droplet_flux_.push_back(-1.0);
    invalidate_active_droplets();
  }

  /// The grid the hydro samples get_source() on; a different grid discards
  /// the normalizations computed so far.
  void set_hydro_grid(const HydroGrid &grid);
  const HydroGrid &get_hydro_grid() const { return hydro_grid_; }

  /// Switch the normalization on the hydro grid on or off (off: every
  /// droplet deposits what its point-sampled kernel sums to, as before).
  void set_normalize_on_hydro_grid(bool on) { normalize_on_grid_ = on; }
  bool get_normalize_on_hydro_grid() const { return normalize_on_grid_; }

  /// Whether this liquefier's source is kernel x p^mu, so that the sampled
  /// sum of the kernel is what a droplet deposits (CausalLiquefier: true).
  virtual bool normalizable_on_grid() const { return false; }

  /// The time a droplet deposits at (CausalLiquefier: tau_d + tau_delay).
  virtual double deposit_time(const Droplet &drop_i) const {
    return drop_i.get_xmu()[0];
  }
  /// smearing_kernel() is non-zero at a query time tau only for |tau -
  /// deposit_time()| within this (CausalLiquefier: its dtau / 2).
  virtual double deposit_half_width() const { return 0.0; }

  /**
   * @brief What the hydro receives from a droplet, as a fraction of its
   * four-momentum: sum over the query times tau_q = tau + k dtau near its
   * deposit time and over the cells of hydro_grid_ of tau_q dtau dx dy deta
   * (n . j), n . j from smearing_kernel() with a unit-energy droplet.
   * @param tau_q_out the query time nearest to the deposit (-1 if no query
   * time of this step reaches it yet).
   */
  double sampled_flux(const Droplet &drop_i, double tau, double dtau,
                      double *tau_q_out = nullptr) const;

  /**
   * @brief Normalize the droplets that deposit within reach of this step.
   *
   * Called once per hydro step, after prepare_active_droplets(), with the
   * step's tau and dtau (the query times are tau and tau + dtau).  Each
   * droplet is normalized once, before its first query: flux > 0 -> its
   * source is scaled by 1/flux; flux = 0 -> it is deposited whole into its
   * nearest cell; outside the grid -> it is lost (counted).
   */
  void normalize_active_droplets(double tau, double dtau);

  /// Sampled flux of droplet idx (-1 if not normalized yet).
  double get_droplet_flux(int idx) const { return droplet_flux_[idx]; }
  /// One-line summary of the normalization of this event's droplets.
  std::string normalization_summary() const;

  /**
   * @brief Restrict get_source() to the droplets that can contribute to a
   * query time in [tau_lo, tau_hi], e.g. one hydro time step with its
   * Runge-Kutta substeps.
   *
   * Uses droplet_may_contribute(). Queries outside the window, and any query
   * before the first call, still loop over every droplet, so a caller that
   * never prepares gets the unpruned result. Adding a droplet or clearing
   * the list drops the window.
   */
  void prepare_active_droplets(double tau_lo, double tau_hi);

  /// Drop the prepared window: get_source() loops over every droplet again.
  void invalidate_active_droplets() {
    active_droplets_.clear();
    active_tau_lo_ = 1.;
    active_tau_hi_ = 0.;
  }

  /// Droplets kept by the last prepare_active_droplets() call.
  int get_number_of_active_droplets() const {
    return static_cast<int>(active_droplets_.size());
  }

  /// Whether a window is prepared, i.e. get_source() at a query time inside
  /// [tau_lo, tau_hi] loops over the kept droplets only.  false after
  /// add_a_droplet() or ClearTask(), until the next prepare_active_droplets().
  bool active_droplets_prepared() const {
    return active_tau_lo_ <= active_tau_hi_;
  }

  /**
   * @brief Whether drop_i can give a non-zero smearing_kernel() for some
   * query time in [tau_lo, tau_hi]. Must never return false for a droplet
   * that contributes. The default keeps every droplet (no pruning).
   */
  virtual bool droplet_may_contribute(const Droplet &drop_i, double tau_lo,
                                      double tau_hi) const {
    return true;
  }

  /**
   * @brief Get number of droplet conversions performed.
   * @return Droplet statistic count
   */
  int get_drop_stat() const { return (drop_stat); }

  /**
   * @brief Get number of partons missed in processing.
   * @return Missed statistic count
   */
  int get_miss_stat() const { return (miss_stat); }

  /**
   * @brief Get number of partons with negative energy.
   * @return Negative statistic count
   */
  int get_neg_stat() const { return (neg_stat); }

  /**
   * @brief Get a specific droplet by index.
   * @param idx Index of the droplet
   * @return Droplet at that index
   */
  Droplet get_a_droplet(const int idx) const { return (dropletlist[idx]); }

  /**
   * @brief Check energy-momentum conservation between partons.
   * @param pIn Input partons
   * @param pOut Output partons
   */
  void check_energy_momentum_conservation(const std::vector<Parton> &pIn,
                                          std::vector<Parton> &pOut);

  /**
   * @brief Apply filtering to remove partons based on criteria.
   * @param pOut Partons to be filtered in-place
   */
  void filter_partons(std::vector<Parton> &pOut);

  /**
   * @brief Add hydrodynamic sources based on input/output partons.
   * @param pIn Input partons
   * @param pOut Output partons
   */
  void add_hydro_sources(std::vector<Parton> &pIn, std::vector<Parton> &pOut);

  // add hydro sources for hadrons is overriden in derived HadronicLiquefier
  // class
  /**
   * @brief Add hydrodynamic sources based on hadrons.
   * @note add hydro sources for hadrons is overriden in derived
   * HadronicLiquefier class
   * @param hIn Input hadrons
   */
  void add_hydro_sources_hadrons(std::vector<Hadron> &hIn){};

  /**
   * @brief Signal used to query hydro cell information.
   */
  sigslot::signal5<double, double, double, double,
                   std::unique_ptr<FluidCellInfo> &,
                   sigslot::multi_threaded_local>
      GetHydroCellSignal;

  /**
   * @brief Check if hydro cell signal is connected.
   * @return True if signal is connected
   */
  const bool get_GetHydroCellSignalConnected() {
    return GetHydroCellSignalConnected;
  }

  /**
   * @brief Set the signal connection flag.
   * @param connected Connection status
   */
  void set_GetHydroCellSignalConnected(bool m_GetHydroCellSignalConnected) {
    GetHydroCellSignalConnected = m_GetHydroCellSignalConnected;
  }

  /**
   * @brief Get number of droplets in the internal list.
   * @return Number of droplets
   */
  int get_dropletlist_size() const { return (dropletlist.size()); }

  /**
   * @brief Compute total energy from all droplets.
   * @return Total energy
   */
  Jetscape::real get_dropletlist_total_energy() const;

  /**
   * @brief Apply smearing kernel to a single droplet.
   * @param tau Proper time
   * @param x Transverse x
   * @param y Transverse y
   * @param eta Space-time rapidity
   * @param drop_i Droplet to smear
   * @param jmu Output source current 4-vector
   */
  virtual void smearing_kernel(Jetscape::real tau, Jetscape::real x,
                               Jetscape::real y, Jetscape::real eta,
                               const Droplet drop_i,
                               std::array<Jetscape::real, 4> &jmu) const {
    jmu = {0, 0, 0, 0};
  }

  /**
   * @brief Accumulate source term at a given space-time point.
   * @param tau Proper time
   * @param x Transverse x
   * @param y Transverse y
   * @param eta Space-time rapidity
   * @param jmu Output source current 4-vector
   */
  /// The source at one cell centre: every droplet's smearing_kernel(), each
  /// scaled by its normalization (see normalize_active_droplets()), plus the
  /// point deposits of droplets whose kernel misses every cell centre.
  void get_source(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                  Jetscape::real eta, std::array<Jetscape::real, 4> &jmu) const;

  /**
   * @brief Get hadronic source term at a given space-time point.
   * @note Functions for the hadronic droplet sources, overriden in derived
   * HadronicLiquefier class.
   * @param tau Proper time
   * @param x Transverse x
   * @param y Transverse y
   * @param eta Space-time rapidity
   * @param jmu Output source current 4-vector
   */
  void get_source_energy(const double tau, const double x, const double y,
                         const double eta, std::array<double, 4> &jmu) const {
    jmu = {0.0, 0.0, 0.0, 0.0};
  };

  /**
   * @brief Get baryon source term at a given space-time point.
   * @param tau Proper time
   * @param x Transverse x
   * @param y Transverse y
   * @param eta Space-time rapidity
   * @return Baryon source term
   */
  double get_source_rhob(const double tau, const double x, const double y,
                         const double eta) const {
    return 0.0;
  };

  /**
   * @brief Get charge source term at a given space-time point.
   * @param tau Proper time
   * @param x Transverse x
   * @param y Transverse y
   * @param eta Space-time rapidity
   * @return Charge source term
   */
  double get_source_rhoq(const double tau, const double x, const double y,
                         const double eta) const {
    return 0.0;
  };

  /**
   * @brief Get strangeness source term at a given space-time point.
   * @param tau Proper time
   * @param x Transverse x
   * @param y Transverse y
   * @param eta Space-time rapidity
   * @return Strangeness source term
   */
  double get_source_rhos(const double tau, const double x, const double y,
                         const double eta) const {
    return 0.0;
  };

  /**
   * @brief Clear all droplets from internal list.
   */
  virtual void ClearTask();
};

};  // namespace Jetscape

#endif  // LIQUEFIERBASE_H
