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
// ------------------------------------------------------------
// This is a hadronic liquefier for the JETSCAPE framework
// It implements a Gaussian smearing kernel and a covariant one
// ------------------------------------------------------------

#ifndef HADRONICLIQUEFIER_H
#define HADRONICLIQUEFIER_H

#include "LiquefierBase.h"
#include "InitialState.h"
#include "RealType.h"

namespace Jetscape {

class HadronDroplet {
private:
  std::array<double, 4> xmu;
  std::array<double, 4> pmu;
  int baryon_number;
  int electric_charge;
  int strangeness;
  double normalization;

public:
  HadronDroplet() = default;
  HadronDroplet(std::array<double, 4> x_in,
          std::array<double, 4> p_in) {
    xmu = x_in;
    pmu = p_in;
  }

  HadronDroplet(std::array<double, 4> x_in,
          std::array<double, 4> p_in,
          int baryon_number_in,
          int electric_charge_in,
          int strangeness_in) {
    xmu = x_in;
    pmu = p_in;
    baryon_number = baryon_number_in;
    electric_charge = electric_charge_in;
    strangeness = strangeness_in;
  }

  ~HadronDroplet(){};

  std::array<double, 4> get_xmu() const { return (xmu); }
  std::array<double, 4> get_pmu() const { return (pmu); }
  int get_baryon_number() const { return baryon_number; }
  int get_electric_charge() const { return electric_charge; }
  int get_strangeness() const { return strangeness; }

  void set_normalization(double norm) { normalization = norm; }
  double get_normalization() const { return normalization; }
};

class HadronicLiquefier : public LiquefierBase{
private:
  bool covariant_smearing_;
  bool hydro_Cartesian_;

  int Nx_, Ny_, Nz_;
  double dx_, dy_, dz_;
  double xMax_, yMax_, zMax_;
  double sigma_transverse_, sigma_longitudinal_;
  int skip_n_sigma_transverse_ = 5;
  int skip_n_sigma_longitudinal_ = 5;
  std::shared_ptr<InitialState> ini;

  std::vector<HadronDroplet> hadron_droplets_list;

  // this is the dtau of the hydro, set from the hydro module
  double dtau_ = -1.0; 


  enum QuantityType { BARYON_NUMBER, ELECTRIC_CHARGE, STRANGENESS };

public:
  /**
   * Default constructor. This reads the parameters from an XML file.
  */
  HadronicLiquefier();
  
  /**
   * Alternative constructor for unit testing. This sets all the needed parameters
   * manually instead of reading them from an XML file.
  */
  HadronicLiquefier(bool covariant, double sigma_transverse,
                    double sigma_longitudinal, double xMax,
                    double yMax, double zMax, int Nx, int Ny,
                    int Nz, bool hydro_Cartesian);

  ~HadronicLiquefier() { ClearTask(); }

  /**
   * Getter functions for the parameters.
  */
  bool get_covariant_smearing() { return covariant_smearing_; }
  bool get_hydro_Cartesian() { return hydro_Cartesian_; }
  int get_Nx() { return Nx_; }
  int get_Ny() { return Ny_; }
  int get_Nz() { return Nz_; }
  double get_dx() { return dx_; }
  double get_dy() { return dy_; }
  double get_dz() { return dz_; }
  double get_xMax() { return xMax_; }
  double get_yMax() { return yMax_; }
  double get_zMax() { return zMax_; }
  double get_sigma_transverse() { return sigma_transverse_; }
  double get_sigma_longitudinal() { return sigma_longitudinal_; }

  /**
   * Set the dtau of the hydro.
  */
  void set_hydro_dtau(double val) { dtau_ = val; };

  /**
   * Function to initialize all parameters from an XML file.
  */
  void InitializeParameters();

  /**
   * Covariant smearing kernel for Milne coordinates. Should be used in case of
   * a Milne evolution in the hydrodynamic model.
  */
  double smearing_kernel_covariant_Milne(const double x_diff,
                                         const double y_diff,
                                         const double eta_diff, const double ux,
                                         const double uy, const double ueta,
                                         const double tau,
                                         const double gamma) const;

  /**
   * Covariant smearing kernel for Cartesian coordinates. Should be used in case 
   * of a Cartesian evolution in the hydrodynamic model.
  */
  double smearing_kernel_covariant_Cartesian(const double x_diff,
                                             const double y_diff,
                                             const double z_diff,
                                             const double ux, const double uy,
                                             const double uz,
                                             const double gamma) const;

  /**
   * Gaussian smearing kernel for Cartesian coordinates.
  */
  double smearing_kernel_gaussian(const double x_diff, const double y_diff,
                                  const double eta_diff) const;

  /**
   * Function to compute the source energy at a given position in spacetime.
  */
  void get_source_energy(const double tau, const double x, const double y,
                         const double eta,
                         std::array<double, 4> &jmu) const;

  /**
   * Function to compute the source for a QuantityType at a given position in 
   * spacetime. The QuantityType can be BARYON_NUMBER, ELECTRIC_CHARGE, or
   * STRANGENESS.
  */
  double get_source_quantity(const double tau, const double x,
                             const double y, const double eta,
                             const QuantityType qtype) const;

  /**
   * Function to compute the source for the baryon density at a given position
   * in spacetime.
  */
  double get_source_rhob(const double tau, const double x, const double y,
                         const double eta) const;
  
  /**
   * Function to compute the source for the charge density at a given position
   * in spacetime.
  */
  double get_source_rhoq(const double tau, const double x, const double y,
                         const double eta) const;
  
  /**
   * Function to compute the source for the strangeness density at a given
   * position in spacetime.
  */
  double get_source_rhos(const double tau, const double x, const double y,
                         const double eta) const;

  /**
   * Function to compute the normalization factor for the smearing kernel and a
   * given Droplet. This is used to normalize the source due to grid effects.
   * In the Cartesian case, the tau parameter is not used. 
  */
  double compute_drop_kernel_normalization(
      const double tau, const HadronDroplet &drop_i) const;

  /**
   * Function to add the hadrons as sources for hydrodynamics by creating 
   * droplets and smear them with the selected smearing kernel.
   * 
   * In case the hydro runs in Cartesian coordinates, the position of the droplets
   * is given in the Cartesian coordinates. In case the hydro runs in Milne
   * coordinates, the position of the droplets is given in Milne coordinates.
  */
  void add_hydro_sources_hadrons(std::vector<Hadron> &hIn);

  /**
   * Function to add a droplet to the list of droplets.
  */
  void add_a_hadronic_droplet(HadronDroplet droplet_in) { 
    hadron_droplets_list.push_back(droplet_in);
  }

  int get_dropletlist_size() const { return (hadron_droplets_list.size()); }

  /**
   * Function to clear the list of droplets.
  */
  void clear_hadron_droplet_list() { hadron_droplets_list.clear(); }

  /**
   * Function to get the total energy of the droplets in the list.
   * This is used to check the energy conservation in the hydro.
  */
  Jetscape::real get_dropletlist_total_energy() const;

  /**
   * Function to get the net baryon number of the droplets in the list.
   * This is used to check the baryon number conservation in the hydro.
   */
  Jetscape::real get_dropletlist_net_baryon_number() const;

  /**
   * Function to get the net electric charge of the droplets in the list.
   * This is used to check the electric charge conservation in the hydro.
   */
  Jetscape::real get_dropletlist_net_electric_charge() const;

  /**
   * Function to get the net strangeness of the droplets in the list.
   * This is used to check the strangeness conservation in the hydro.
   */
  Jetscape::real get_dropletlist_net_strangeness() const;

  virtual void ClearTask();
};

}; // namespace Jetscape

#endif // HADRONICLIQUEFIER_H