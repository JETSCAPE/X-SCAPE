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

#ifndef JETSCAPE_H
#define JETSCAPE_H

#include "JetScapeLogger.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeModuleBase.h"
#include "CausalLiquefier.h"
#include "HadronicLiquefier.h"
#include "HadronicEMT.h"
#include <unordered_map>

namespace Jetscape {

/**
 * @class JetScape
 * @brief Main driver class for the JetScape framework.
 *
 * The JetScape class orchestrates the execution of the JetScape framework.
 * It serves as the top-level task, responsible for initializing modules,
 * executing event loops, controlling event reuse (hydro), and coordinating
 * the writing of simulation outputs. The task list and writers can be
 * automatically determined from XML configuration files or set manually.
 *
 * Example usage:
 * @code
 * auto jetscape = make_shared<JetScape>();
 * jetscape->SetNumberOfEvents(1000);
 * jetscape->Init();
 * jetscape->Exec();
 * jetscape->Finish();
 * @endcode
 */

class JetScape : public JetScapeModuleBase,
                 public std::enable_shared_from_this<JetScape> {
 public:
  /**
   * @brief Default constructor.
   *
   * Creates the main JetScape task object.
   */
  JetScape();

  /**
   * @brief Destructor.
   *
   * Cleans up resources and finalizes framework execution.
   */
  virtual ~JetScape();

  /**
   * @brief Initialize the JetScape framework.
   *
   * Reads configuration from XML files and initializes all attached modules
   * and sub-tasks. Calls JetScapeTask::InitTasks() for recursive setup.
   */
  void Init() override;

  /**
   * @brief Execute the JetScape framework.
   *
   * Runs all tasks for the requested number of events. During execution, parton
   * showers may be printed via GetPartons(), and results are stored using
   * WriteTasks().
   */
  void Exec() override;

  /**
   * @brief Finalize the JetScape framework.
   *
   * Cleans up after execution, finalizes writers, and ensures proper shutdown.
   */
  void FinishTask() override;

  /**
   * @brief Set the total number of events to be generated.
   * @param m_n_events Number of events.
   */
  void SetNumberOfEvents(int m_n_events) { n_events = m_n_events; }

  /**
   * @brief Get the total number of events.
   * @return Number of events.
   */
  int GetNumberOfEvents() { return n_events; }

  /**
   * @brief Enable or disable hydro event reuse.
   * @param reuse_hydro Boolean flag to enable reuse.
   *
   * If enabled, hydro events can be reused multiple times
   * (controlled by SetNReuseHydro).
   */
  inline void SetReuseHydro(const bool reuse_hydro) {
    reuse_hydro_ = reuse_hydro;
  }

  /**
   * @brief Check if hydro event reuse is enabled.
   * @return True if hydro reuse is enabled, false otherwise.
   */
  inline bool GetReuseHydro() const { return reuse_hydro_; }

  /**
   * @brief Set the number of times a hydro event will be reused.
   * @param n_reuse_hydro Number of reuses.
   *
   * Note: Reuse must first be enabled with SetReuseHydro(true),
   * otherwise a warning will be issued.
   */
  inline void SetNReuseHydro(const unsigned int n_reuse_hydro) {
    if (!GetReuseHydro()) {
      JSWARN << "Number of hydro reusals set, but reusal not turned on.";
      JSWARN << "Try jetscape->SetReuseHydro (true);";
    }
    n_reuse_hydro_ = n_reuse_hydro;
  }

  /**
   * @brief Get the number of hydro event reuses.
   * @return Number of reuses.
   */
  inline unsigned int GetNReuseHydro() const { return n_reuse_hydro_; }

 protected:
  /**
   * @brief Compare XML configuration elements.
   *
   * Utility function for validating XML consistency.
   */
  void CompareElementsFromXML();

  /**
   * @brief Recursively build elements from XML.
   * @param elems Vector of XML element names.
   * @param mElement Pointer to XML element.
   */
  void recurseToBuild(std::vector<std::string> &elems,
                      tinyxml2::XMLElement *mElement);

  /**
   * @brief Recursively search XML elements.
   * @param elems Vector of XML element names.
   * @param uElement Pointer to XML element.
   */
  void recurseToSearch(std::vector<std::string> &elems,
                       tinyxml2::XMLElement *uElement);

  /**
   * @brief Read general parameters from XML configuration.
   */
  void ReadGeneralParametersFromXML();

  /**
   * @brief Determine task list from XML configuration.
   *
   * This may automatically set up module execution order.
   */
  void DetermineTaskListFromXML();

  /**
   * @brief Determine writer modules from XML configuration.
   */
  void DetermineWritersFromXML();

  /**
   * @brief Check for a specific writer in XML configuration.
   * @param writerName Name of the writer module.
   * @param outputFilename Output file name to associate with the writer.
   */
  void CheckForWriterFromXML(const char *writerName,
                             std::string outputFilename);

  /**
   * @brief Set the ID of a module from its XML configuration.
   * @param moduleElement Pointer to XML element for the module.
   * @param module Shared pointer to the module object.
   */
  void SetModuleId(tinyxml2::XMLElement *moduleElement,
                   shared_ptr<JetScapeModuleBase> module);

  /**
   * @brief Establish pointers to required submodules and tasks.
   */
  void SetPointers();

  /** Function to set the per event execution active flag so that
  if hadronization and Afterburner are attached and not per time step executed,
  that they will be automatically executed after the per time step modules are
  finished So currently possible workflow automatically executed correctly is:
  per event -> per timestep -> per event
   */

  /**
   * @brief Set per-event execution flags for modules.
   * 
   * Function to set the per event execution active flag so that
   * if hadronization and Afterburner are attached and not per time step executed,
   * that they will be automatically executed after the per time step modules are
   * finished So currently possible workflow automatically executed correctly is:
   * per event -> per timestep -> per event
   * 
   * @param start_of_event Boolean flag indicating whether this is the start
   * of an event.
   */
  void SetPerEventExecFlags(bool start_of_event);

  /**
   * @brief Reset per-event execution flags to original state.
   * 
   * Function to reset the per event execution active flags to its original
   * state.
   */
  void ResetPerEventExecFlags();

  /**
   * @brief Display JetScape banner.
   */
  void Show();

  /** @brief Number of events to run. */
  int n_events;

  /** @brief Interval for event progress printouts. */
  int n_events_printout;

  /** @brief Flag to enable/disable hydro reuse. */
  bool reuse_hydro_;

  /** @brief Number of hydro event reuses. */
  unsigned int n_reuse_hydro_;

  /** @brief Pointer to causal liquefaction module (optional). */
  std::shared_ptr<CausalLiquefier> liquefier;
  std::shared_ptr<HadronicLiquefier> hadronicLiquefier;
  std::shared_ptr<HadronicEMT> hadronicEMT;

  /**
   * @brief Flag for automatic task list determination from XML.
   *
   * If true, tasks are set up based on XML configuration instead of manual
   * calls to JetScapeTask::Add().
   */
  bool fEnableAutomaticTaskListDetermination;

  /** @brief List to store original SetActive flag settings. */
  std::unordered_multimap<int, bool> taskOrgActiveMap;
};

}  // end namespace Jetscape

#endif
