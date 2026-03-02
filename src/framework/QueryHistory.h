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

/**
 * @file QueryHistory.h
 * @brief Query history singleton used to find and query module histories.
 *
 * Provides a global access point to the module task map and helpers to
 * retrieve stored history data from modules derived from
 * `JetScapeModuleBase`.
 */

#ifndef QUERYHISTORY_H
#define QUERYHISTORY_H

#include "JetScapeModuleBase.h"
#include "cpp17/any.hpp"
#include "cpp17/variant.hpp"

#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
//#include <cstddef>

//#include "sigslot.h"
// using namespace sigslot;

// JP: Make sure to not introduce and memory leak here using an instance.
// Should be fine (see SignalManager) by using weak pointers. But make sure
// !!!!!

// Maybe change, namespaces for any and varaint ...
using namespace linb;
using namespace mpark;

namespace Jetscape {

/**
 * @class QueryHistory
 * @brief Singleton that stores a map of tasks and provides access to module
 *        histories.
 *
 * `QueryHistory` keeps a (multi)map of task id strings to weak pointers of
 * `JetScapeTask` instances. It is intended to be used as a global (singleton)
 * registry so other parts of the framework can lookup modules by name and
 * request their internal history (where applicable).
 */
class QueryHistory {
 public:
  /**
   * @brief Get the singleton instance of `QueryHistory`.
   *
   * If the instance does not yet exist it will be created.
   * @return pointer to the `QueryHistory` singleton.
   */
  static QueryHistory *Instance();

  /**
   * @brief Save the root (main) task used to build the internal task map.
   * @param m_main_task shared pointer to the main `JetScapeTask`.
   */
  void AddMainTask(std::shared_ptr<JetScapeTask> m_main_task) { main_task = m_main_task; }

  /**
   * @brief Rebuild the internal `taskMap` from the currently set main task.
   *
   * This clears the existing map and walks the main task's children (one
   * level deep) inserting task id -> weak_ptr pairs.
   */
  void UpdateTaskMap();

  /**
   * @brief Print the configured tasks hierarchy to the JetScape logger.
   */
  void PrintTasks();

  /**
   * @brief Print the contents of the internal `taskMap` with details.
   *
   * Each entry prints the task id, pointer value, active/multithread flags
   * and task number. If the task can be cast to `JetScapeModuleBase` the
   * `IsTimeStepped()` flag is also logged.
   */
  void PrintTaskMap();

  /**
   * @brief Return a copy of the internal task multimap.
   * @return unordered_multimap of task id to weak pointer of `JetScapeTask`.
   */
  std::unordered_multimap<std::string, std::weak_ptr<JetScapeTask>> GetTaskMap() { return taskMap; }

  /**
   * @brief Query a single module's history by module/ task name.
   * @param mName module/task id string to search for.
   * @return `any` containing the module history or an empty/zero value when
   *         not found or on error.
   */
  any GetHistoryFromModule(string mName);

  /**
   * @brief Query histories from all modules matching `mName`.
   *
   * This returns a `vector<any>` containing the `GetHistory()` results from
   * every module whose task id equals `mName`.
   *
   * @param mName module/task id string to search for.
   * @return vector of `any` objects with each module's history.
   */
  vector<any> GetHistoryFromModules(string mName);

 private:
  QueryHistory(){};
  QueryHistory(QueryHistory const &){};

  /**
   * @brief Pointer to the singleton instance.
   */
  static QueryHistory *m_pInstance;

  /**
   * @brief Map of task id -> weak pointer to the task instance.
   */
  std::unordered_multimap<std::string, std::weak_ptr<JetScapeTask>> taskMap;

  /**
   * @brief Weak pointer to the root/main `JetScapeTask` used to populate
   *        `taskMap`.
   */
  std::weak_ptr<JetScapeTask> main_task;
};

}  // end namespace Jetscape

#endif