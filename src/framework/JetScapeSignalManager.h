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

//  SignalManager instance class (meant as singelton)

#ifndef JETSCAPESIGNALMANAGER_H
#define JETSCAPESIGNALMANAGER_H

#include "Afterburner.h"
#include "Transport.h"
#ifdef USE_SMASH
#include "SMASHInitialStateWrapper.h"
#include "SMASHNucleusWrapper.h"
#endif
#include "InitialState.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "IsrManager.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "FluidDynamics.h"
#include "BulkDynamicsManager.h"
#include "HardProcess.h"
#include "JetScapeWriter.h"
#include "PreequilibriumDynamics.h"
#include "PartonPrinter.h"
#include "HadronPrinter.h"

#include <iostream>
#include <string>
#include <map>
#include "sigslot.h"

using namespace sigslot;

namespace Jetscape {

/**
 * @class JetScapeSignalManager
 * @brief Singleton class for managing signal-slot connections between
 *        JETSCAPE modules.
 *
 * The JetScapeSignalManager orchestrates the communication between different
 * modules of the JETSCAPE framework (e.g., initial state, hydro, energy loss,
 * hadronization, etc.) via the `sigslot` signal/slot mechanism. It maintains
 * and connects module interfaces such that partons, hydrodynamic cells, and
 * hadrons can be passed consistently between simulation stages.
 *
 * This class is designed as a singleton. Use JetScapeSignalManager::Instance()
 * to retrieve the global instance.
 *
 * Responsibilities:
 * - Store and provide weak pointers to all major physics modules.
 * - Establish signal-slot connections between modules (e.g., hydro ↔ jet energy
 *  loss, hard process ↔ hadronization).
 * - Maintain internal maps of connected signals for debugging and cleanup.
 * - Print and clean up signal maps.
 */

class
    JetScapeSignalManager  //: public
                           //: sigslot::has_slots<sigslot::multi_threaded_local>
{
 public:
  /**
   * @brief Retrieve the singleton instance of the SignalManager.
   *
   * If no instance exists, it is created. This ensures global coordination of
   * module signal connections.
   *
   * @return Pointer to the singleton JetScapeSignalManager.
   */
  static JetScapeSignalManager *Instance();

  /// @name Module pointer setters and getters
  /// Set and retrieve weak references to JETSCAPE modules.
  ///@{

  /// Set the InitialState module.
  void SetInitialStatePointer(shared_ptr<InitialState> m_initial) {
    initial_state = m_initial;
  }
  /// Get the InitialState module.
  weak_ptr<InitialState> GetInitialStatePointer() { return initial_state; }

#ifdef USE_SMASH
  /// Set the SMASH initial state wrapper module.
  void SetSMASHInitialStatePointer(
      shared_ptr<SmashInitialConditionWrapper> m_SMASH_initial) {
    transport_initial_state = m_SMASH_initial;
  }
  /// Get the SMASH initial state wrapper module.
  weak_ptr<SmashInitialConditionWrapper> GetSMASHInitialStatePointer() {
    return transport_initial_state;
  }
  /// Set the SMASH nucleus wrapper module.
  void SetSMASHNucleusPointer(shared_ptr<SMASHNucleusWrapper> m_SMASH_nucleus) {
    smash_nucleus = m_SMASH_nucleus;
  }
  /// Get the SMASH nucleus wrapper module.
  weak_ptr<SMASHNucleusWrapper> GetSMASHNucleusPointer() {
    return smash_nucleus;
  }
#endif

  /// Set the Pre-equilibrium dynamics module.
  void SetPreEquilibriumPointer(shared_ptr<PreequilibriumDynamics> m_pre_eq) {
    pre_equilibrium = m_pre_eq;
  }
  /// Get the Pre-equilibrium dynamics module.
  weak_ptr<PreequilibriumDynamics> GetPreEquilibriumPointer() {
    return pre_equilibrium;
  }

  /// Set the hydrodynamics module.
  void SetHydroPointer(shared_ptr<FluidDynamics> m_hydro) { hydro = m_hydro; }
  /// Get the hydrodynamics module.
  weak_ptr<FluidDynamics> GetHydroPointer() { return hydro; }

  /// Set the Bulk Dynamics Manager module.
  void SetBulkDynamicsManagerPointer(shared_ptr<BulkDynamicsManager> m_bulk) {
    bulk = m_bulk;
  }
  /// Get the Bulk Dynamics Manager module.
  weak_ptr<BulkDynamicsManager> GetBulkPointer() { return bulk; }

  /// Set the soft particlization module.
  void SetSoftParticlizationPointer(shared_ptr<SoftParticlization> m_soft) {
    softparticlization = m_soft;
  }
  /// Get the soft particlization module.
  weak_ptr<SoftParticlization> GetSoftParticlizationPointer() {
    return softparticlization;
  }

  /// Set the JetEnergyLossManager module.
  void SetJetEnergyLossManagerPointer(
      shared_ptr<JetEnergyLossManager> m_jloss) {
    jloss = m_jloss;
  }
  /// Get the JetEnergyLossManager module.
  weak_ptr<JetEnergyLossManager> GetJetEnergyLossManagerPointer() {
    return jloss;
  }

  /// Set the HardProcess module.
  void SetHardProcessPointer(shared_ptr<HardProcess> m_hardp) {
    hardp = m_hardp;
  }
  /// Get the HardProcess module.
  weak_ptr<HardProcess> GetHardProcessPointer() { return hardp; }

  /// Set the writer module (JetScapeWriter).
  void SetWriterPointer(shared_ptr<JetScapeWriter> m_writer) {
    writer = m_writer;
  }
  /// Get the writer module.
  weak_ptr<JetScapeWriter> GetWriterPointer() { return writer; }

  /// Set the HadronizationManager module.
  void SetHadronizationManagerPointer(
      shared_ptr<HadronizationManager> m_hadro) {
    hadro = m_hadro;
  }
  /// Get the HadronizationManager module.
  weak_ptr<HadronizationManager> GetHadronizationManagerPointer() {
    return hadro;
  }

  /// Set the PartonPrinter module.
  void SetPartonPrinterPointer(shared_ptr<PartonPrinter> m_pprinter) {
    pprinter = m_pprinter;
  }
  /// Get the PartonPrinter module.
  weak_ptr<PartonPrinter> GetPartonPrinterPointer() { return pprinter; }

  /// Set the HadronPrinter module.
  void SetHadronPrinterPointer(shared_ptr<HadronPrinter> m_hprinter) {
    hprinter = m_hprinter;
  }
  /// Get the HadronPrinter module.
  weak_ptr<HadronPrinter> GetHadronPrinterPointer() { return hprinter; }

  /// Set the JetEnergyLoss module.
  void SetEnergyLossPointer(shared_ptr<JetEnergyLoss> m_eloss) {
    eloss = m_eloss;
  }

  /// Get the JetEnergyLoss module.
  weak_ptr<JetEnergyLoss> GetEnergyLossPointer() { return eloss; }

  /// @name Signal connection methods
  /// Establish signal-slot connections between modules.
  ///@{

  void ConnectJetSignal(shared_ptr<JetEnergyLoss> j);
  void ConnectEdensitySignal(shared_ptr<JetEnergyLoss> j);
  void ConnectGetHydroTau0Signal(shared_ptr<JetEnergyLoss> j);
  void ConnectGetHydroCellSignal(shared_ptr<JetEnergyLoss> j);
  void ConnectGetHydroCellSignal(shared_ptr<LiquefierBase> l);
  void ConnectGetHydroCellSignal(shared_ptr<Hadronization> h);
  void ConnectGetHardPartonListSignal(shared_ptr<JetEnergyLossManager> jm);
  void ConnectGetHardPartonListSignal(shared_ptr<IsrManager> jm);
  void ConnectSentInPartonsSignal(shared_ptr<JetEnergyLoss> j,
                                  shared_ptr<JetEnergyLoss> j2);
  void ConnectGetFinalPartonListSignal(shared_ptr<HadronizationManager> hm);
  void ConnectTransformPartonsSignal(shared_ptr<Hadronization> h,
                                     shared_ptr<Hadronization> h2);
  void ConnectGetFinalHadronListSignal(shared_ptr<HadronPrinter> h);

  void ConnectGetHydroHyperSurfaceSignal(shared_ptr<Hadronization> h);
  void ConnectGetHydroHyperSurfaceSignal(shared_ptr<SoftParticlization> hSoft);
  void ConnectClearHydroHyperSurfaceSignal(
      shared_ptr<SoftParticlization> hSoft);

  /// Disconnect signals (not fully implemented).
  void DisconnectSignal(){};

  ///@}

  /**
   * @brief Clean up signal maps and counters.
   *
   * Removes obsolete signal connections and clears internal maps.
   * This is typically used at the end of an event or when resetting
   * the framework.
   */

  void CleanUp();

  /// @name Signal statistics
  /// Retrieve counters for established signal connections.
  ///@{
  int GetNumberOfJetSignals() { return num_jet_signals; }
  int GetNumberOfEdensitySignals() { return num_edensity_signals; }
  int GetNumberOfGetHydroCellSignals() { return num_GetHydroCellSignals; }
  ///@}

  /// @name Debugging utilities
  /// Print internal maps of connected signals.
  ///@{
  void PrintJetSignalMap();
  void PrintEdensitySignalMap();
  void PrintGetHydroCellSignalMap();
  void PrintSentInPartonsSignalMap();
  void PrintTransformPartonsSignalMap();
  ///@}

 private:
  /// Private constructor (singleton pattern).
  JetScapeSignalManager(){};
  /// Deleted copy constructor (singleton pattern).
  JetScapeSignalManager(JetScapeSignalManager const &){};

  /// Singleton instance pointer.
  static JetScapeSignalManager *m_pInstance;

  // Weak references to JETSCAPE modules
  weak_ptr<InitialState> initial_state;
  weak_ptr<PreequilibriumDynamics> pre_equilibrium;
#ifdef USE_SMASH
  weak_ptr<SmashInitialConditionWrapper> transport_initial_state;
  weak_ptr<SMASHNucleusWrapper> smash_nucleus;
#endif
  weak_ptr<FluidDynamics> hydro;
  weak_ptr<BulkDynamicsManager> bulk;
  weak_ptr<JetEnergyLossManager> jloss;
  weak_ptr<HardProcess> hardp;
  weak_ptr<JetScapeWriter> writer;
  weak_ptr<HadronizationManager> hadro;
  weak_ptr<Afterburner> afterburner;
  weak_ptr<Transport> transport;
  weak_ptr<PartonPrinter> pprinter;
  weak_ptr<HadronPrinter> hprinter;
  weak_ptr<JetEnergyLoss> eloss;
  weak_ptr<SoftParticlization> softparticlization;

  // Counters for signals
  int num_jet_signals = 0;
  int num_edensity_signals = 0;
  int num_GetHydroCellSignals = 0;
  int num_SentInPartons = 0;
  int num_TransformPartons = 0;

  // Maps for tracking connected signals
  map<int, weak_ptr<JetEnergyLoss>> jet_signal_map;
  map<int, weak_ptr<JetEnergyLoss>> edensity_signal_map;
  map<int, weak_ptr<JetEnergyLoss>> GetHydroCellSignal_map;

  map<int, weak_ptr<JetEnergyLoss>> SentInPartons_map;
  map<int, weak_ptr<Hadronization>> TransformPartons_map;
};

}  // end namespace Jetscape

#endif
