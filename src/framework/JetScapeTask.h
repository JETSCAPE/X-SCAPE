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

/// \todo Implement truly recursively (not yet done)
/// (see https://root.cern.ch/doc/v608/TTask_8cxx_source.html l248)
/// (not really though ...)

#ifndef JETSCAPETASK_H
#define JETSCAPETASK_H

#include <list>
#include <vector>
#include <string>
#include <memory>

using std::shared_ptr;
using std::string;
using std::vector;
using std::weak_ptr;

namespace Jetscape {

// need forward declaration
class JetScapeWriter;
class JetScapeModuleMutex;
class Parton;

/**
 * @class JetScapeTask
 * @brief Base class for modular JetScape tasks.
 *
 * Provides a task-based interface for executing, initializing,
 * finalizing, and writing modules that make up a simulation.
 */
class JetScapeTask {
 public:
  /**
   * @brief Default constructor.
   *
   * Initializes the task with `active_exec` set to true and `id` to a
   * default value. `my_task_number_` is assigned using the `RegisterTask`
   * method.
   */
  JetScapeTask();

  /**
   * @brief Virtual destructor.
   */
  virtual ~JetScapeTask();

  /**
   * @brief Virtual initialization function.
   * Can be overridden by derived tasks.
   */
  virtual void Init();

  /** Init function for a JetScapeTask overridden by different modules/tasks.
   */
  virtual void InitTask(){};

  /** Recursive initialization of all the subtasks of the JetScapeTask i.e.
   * calling Init() for subtasks. Subtasks are of type JetScapeTask such as
   * Pythia Gun, Trento, Energy Loss Matter and Martini etc.
   */
  virtual void InitTasks();

  /**
   * @brief Virtual execution function.
   * Execute task itself (via ExecuteTask()) and call all subtasks. By default,
   * tasks only need to define an ExecuteTask() function as the subtasks will
   * automatically be called by this function. In special cases, this function
   * can also be overridden and the subtasks can be called in a specific way.
   */
  virtual void Exec();

  /**
   * @brief Virtual execution function.
   * Can be overridden by derived tasks.
   */
  virtual void ExecuteTask();

  /**
   * @brief Virtual execution function. 
   * Recursive execution of all subtasks of the JetScapeTask i.e. calling
   * Exec() for subtasks. If one wants to really make sure that no subtasks are
   * done, overwrite ExecuteTasks with empty function.
   */
  virtual void ExecuteTasks();

  /**
   * @brief Virtual clear function.
   * Clear task itself (via ClearTask()) and call the clear functions
   * for all subtasks. Usually not overriden by modules/tasks
   */
  virtual void Clear();

  /**
   * @brief Virtual clear function.
   * Clear function for a JetScapeTask overridden by different modules/tasks.
   */
  virtual void ClearTask(){};

  /** Recursively clearing of all subtasks of a JetScapeTask.
   */
  virtual void ClearTasks();

  /**
   * @brief Virtual finish function. 
   * Finish task itself (via FinishTask()) and call the finish functions
   * for all subtasks. Usually not overriden by modules/tasks.
   */
  virtual void Finish();

  /**
   * @brief Virtual finish function.
   * Finish function for a JetScapeTask overridden by different modules/tasks.
   */
  virtual void FinishTask(){};

  /**
   * @brief Virtual finish function.
   * Recursive finalization of all the subtasks of the JetScapeTask.
   */
  virtual void FinishTasks();

  /**
   * @brief Recursively write output from all subtasks.
   *
   * The active_exec flag is used to decide whether to write the output in
   * the file or not.
   *
   * @param w Weak pointer to the JetScapeWriter instance.
   */
  virtual void WriteTasks(weak_ptr<JetScapeWriter> w);

  /**
   * @brief Write output from an individual task.
   *
   * This is the default virtual WriteTask() function for a JetScapeTask.
   * It can be overridden by different modules/tasks.
   *
   * Current setup: Every task gets handed a pointer to the writer and can add
   * any information it likes using predefined functions like `WriteComment()`.
   *
   * This approach is flexible, but makes it difficult to properly store
   * information across multiple output formats. Ideally, writers should collect
   * and process this information intelligently.
   *
   * @param w Weak pointer to the JetScapeWriter instance.
   */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w){};

  /**
   * @brief Collect header for this task only.
   *
   * Should only be called by `CollectHeaders()`.
   *
   * @param w Weak pointer to the JetScapeWriter instance.
   */
  virtual void CollectHeader(weak_ptr<JetScapeWriter> w){};

  /**
   * @brief Recursively collect header information from all subtasks.
   *
   * Uses `active_exec` flag to decide whether to include output.
   *
   * @param w Weak pointer to the JetScapeWriter instance.
   */
  virtual void CollectHeaders(weak_ptr<JetScapeWriter> w);

  /**
   * @brief Add a task to the subtask list.
   *
   * @param m_tasks Shared pointer to the task to be added.
   */
  virtual void Add(shared_ptr<JetScapeTask> m_tasks);

  /**
   * @brief Get the task number of this task.
   *
   * @return The task number.
   */
  virtual const inline int GetMyTaskNumber() const { return my_task_number_; };

  /**
   * @brief Get the vector of subtasks.
   *
   * @return A vector of shared pointers to subtasks.
   */
  const vector<shared_ptr<JetScapeTask>> GetTaskList() const { return tasks; }

  /**
   * @brief Get a subtask at a specific index.
   *
   * @param i Index of the task.
   * @return Shared pointer to the task.
   */
  shared_ptr<JetScapeTask> GetTaskAt(int i) { return tasks.at(i); }

  /**
   * @brief Delete the last task in the subtask list.
   */
  void EraseTaskLast() { tasks.erase(tasks.begin() + (tasks.size() - 1)); }
  // funny syntax (why last() not working here!?)

  /**
   * @brief Delete the task at a specific index.
   *
   * @param i Index of the task to delete.
   */
  void EraseTaskAt(int i) { tasks.erase((tasks.begin() + i)); }

  /**
   * @brief Resize the subtask list.
   *
   * @param i New size of the task list.
   * If `i` is less than the current size, the list is truncated.
   */
  void ResizeTaskList(int i) { tasks.resize(i); }

  /**
   * @brief Remove all tasks from the subtask list.
   */
  void ClearTaskList() { tasks.clear(); }

  /**
   * @brief Get the number of subtasks.
   *
   * @return Number of tasks in the list.
   */
  int GetNumberOfTasks() { return (int)tasks.size(); }

  /**
   * @brief Check whether this task is active.
   *
   * @return True if active, false otherwise.
   */
  const bool GetActive() const { return active_exec; }

  /**
   * @brief Set the active execution flag.
   *
   * @param m_active_exec True to activate, false to deactivate.
   */
  void SetActive(bool m_active_exec) { active_exec = m_active_exec; }
  // needed to access tasks not recursively by default but individually ...
  // also usefull to prevent hydro if multiple read ins of the same event ...

  /**
   * @brief Set the task ID string.
   *
   * @param m_id String ID of the task.
   */
  void SetId(string m_id) { id = m_id; }

   /**
   * @brief Get the ID string of the task.
   *
   * @return The ID string.
   */
  const string GetId() const { return id; }

  /**
   * @brief Get the task's mutex.
   *
   * @return Shared pointer to the JetScapeModuleMutex.
   */
  const shared_ptr<JetScapeModuleMutex> GetMutex() const { return mutex; }

  /**
   * @brief Set the mutex object for this task.
   *
   * @param m_mutex Shared pointer to the mutex.
   */
  void SetMutex(shared_ptr<JetScapeModuleMutex> m_mutex) { mutex = m_mutex; }

  /**
   * @brief Get the multi-threading flag.
   * @return True if multi-threading is enabled, false otherwise.
   */
  bool GetMultiThread() { return multiThread; }

  /**
   * @brief Set the multi-threading flag for this task.
   * @param m_multiThread True to enable multi-threading, false to disable.
   */
  void SetMultiThread(bool m_multiThread) { multiThread = m_multiThread; }

  /**
   * @brief Get the task number of this task.
   * @return The task number.
   */
  int GetTaskNumber() { return my_task_number_; }

 private:
  /**
   * @brief Vector of subtasks.
   *
   * Can be made sortable to match predefined execution order.
   */
  vector<shared_ptr<JetScapeTask>> tasks;
  // list<shared_ptr<JetScapeTask>> tasks; // list vs vector any advantage of
  // list?

  /**
   * @brief Flag indicating whether the task is active.
   */
  bool active_exec;

  /**
   * @brief Flag indicating whether multi-threading is enabled for this task.
   */
  bool multiThread = false;

  /**
   * @brief Identifier string for this task.
   */
  string id;
  // if for example a search rather position ... (or always sort with predefined
  // order!?)

  /**
   * @brief Internal task number.
   */
  int my_task_number_;

  /**
   * @brief Mutex object used for thread safety.
   */
  shared_ptr<JetScapeModuleMutex> mutex;
};

}  // end namespace Jetscape

#endif
