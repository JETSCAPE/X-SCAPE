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

// JetScapeTask class
// \TODO: Implement truly recursively (not yet done)
// (see https://root.cern.ch/doc/v608/TTask_8cxx_source.html l248)
// (not really though ...)

#ifndef JETSCAPETASK_H
#define JETSCAPETASK_H

#include <list>
#include <vector>
#include <string>
#include <memory>

using std::vector;
using std::string;
using std::weak_ptr;
using std::shared_ptr;

namespace Jetscape {

// need forward declaration
class JetScapeWriter;
class JetScapeModuleMutex;
class Parton;

/**
 * A JetscapeTask has 4 central functions that it needs to provide: InitTask(),
 * ExecuteTask(), FinishTask() and ClearTask(). They are automatically called
 * recursively for all subtasks by calling Init(), Execute(), Finish() and
 * Clear() once from the top-level main (Jetscape) task to which all tasks are
 * added. Tasks usually only need to define there ...Task() functions.
 */
class JetScapeTask {

public:
  /** Default constructor to create a JetScapeTask. It sets the flag "active_exec" to true and  "id" to default string value.
   */
  JetScapeTask();

  /** Default destructor for a JetScapeTask.
   */
  virtual ~JetScapeTask();

  /** Initialize task itself (by calling InitTask()) and calling subtasks (via
   * InitTasks()). Usually not overriden by modules/tasks.
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

  /** Execute task itself (via ExecuteTask()) and call all subtasks. By default,
   * tasks only need to define an ExecuteTask() function as the subtasks will
   * automatically be called by this function. In special cases, this function
   * can also be overridden and the subtasks can be called in a specific way.
   */
  virtual void Exec();

  /** Main execution function for a JetScapeTask overridden by different modules/tasks.
   */
  virtual void ExecuteTask();

  /** Recursive execution of all subtasks of the JetScapeTask i.e. calling Exec()
   * for subtasks. If one wants to really make sure that no subtasks are done,
   * overwrite ExecuteTasks with empty function.
   */
  virtual void ExecuteTasks();

  /** Clear task itself (via ClearTask()) and call the clear functions
   * for all subtasks. Usually not overriden by modules/tasks
   */
  virtual void Clear();

  /** Clear function for a JetScapeTask overridden by different modules/tasks.
   */
  virtual void ClearTask(){};

  /** Recursively clearing of all subtasks of a JetScapeTask.
   */
  virtual void ClearTasks();

  /** Finish task itself (via FinishTask()) and call the finish functions
   * for all subtasks. Usually not overriden by modules/tasks.
   */
  virtual void Finish();

  /** Finish function for a JetScapeTask overridden by different modules/tasks.
   */
  virtual void FinishTask(){};

  /** Recursive finalization of all the subtasks of the JetScapeTask.
   */
  virtual void FinishTasks();

  /** Recursively write the output information of different tasks/subtasks of a JetScapeTask into a file.
      We use "active_exec" flag to decide whether to write the output in the file or not.
  */
  virtual void WriteTasks(weak_ptr<JetScapeWriter> w);

  //add here a write task (depending on if JetScapeWriter is initiallized and active) ...
  // Think about workflow ...
  /** A virtual function to define a default WriteTask() function for a JetScapeTask.
      It can be overridden by different modules/tasks.
      Current setup: Every task gets handed a pointer to the writer
      and can add any information it likes to it
      (using predefined functions like WriteComment())
      This is maximally flexible but makes it difficult to
      properly store information for a variety of outputs.
      E.g., sigmaGen: A HardProcess can easily write the xsec to any stream-type output
      using WriteComment. But to set it in a HepMC file, either HardProcess needs to
      make a case-by-case selection, meaning a new file format would need to percolate
      through multiple base classes, or the writer needs to know this information
      and implement WriteEvent appropriately. The latter is obviously better, but
      it's non-trivial to collect this information.
   */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w){};

  /** Should get called only by CollectHeaders. Maybe make protected?
      @param w is a pointer of type JetScapeWrite class.
  */
  virtual void CollectHeader(weak_ptr<JetScapeWriter> w){};

  /** Recursively collect the header information of different tasks/subtasks of a JetScapeTask into a writer.
      We use "active_exec" flag to decide whether to write the output in the file or not.
  */
  virtual void CollectHeaders(weak_ptr<JetScapeWriter> w);

  /** This function adds the module "m_tasks" into the vector of subtask of a JetScapeTask.
  */
  virtual void Add(shared_ptr<JetScapeTask> m_tasks);

  /** This function returns the current task number.
   */
  virtual const inline int GetMyTaskNumber() const { return my_task_number_; };

  /** This function returns the vector of tasks of a JetScapeTask.
   */
  const vector<shared_ptr<JetScapeTask>> GetTaskList() const { return tasks; }

  /** This function returns the task at  ith location in the vector of tasks of a JetScapeTask.*/
  shared_ptr<JetScapeTask> GetTaskAt(int i) { return tasks.at(i); }

  /** This function deletes the last task in the vector  of tasks of a JetScapeTask. */
  void EraseTaskLast() { tasks.erase(tasks.begin() + (tasks.size() - 1)); }
  //funny syntax (why last() not working here!?)

  /** This function deletes the task at ith location in the vector of tasks of a JetScapeTask.
   */
  void EraseTaskAt(int i) { tasks.erase((tasks.begin() + i)); }

  /** This function resizes the length of the vector of tasks to "i". If "i" is less than the current size, it will keep the first i elements of the vector of the tasks of a JetScapeTask.
   */
  void ResizeTaskList(int i) { tasks.resize(i); }

  /** This function removes all the tasks in the vector of tasks of a JetScapeTask and changes the size to 0.
   */
  void ClearTaskList() { tasks.clear(); }

  /** This function returns the number of tasks of a JetScapeTask stored in the vector array of tasks.
   */
  int GetNumberOfTasks() { return (int)tasks.size(); }

  /** This function tells whether the task is active or not i.e. whether its Exec() function is run.
   */
  const bool GetActive() const { return active_exec; }

  /** This functions sets the flag "active_exec" to true (active) or false(deactive).
   */
  void SetActive(bool m_active_exec) { active_exec = m_active_exec; }
  // needed to access tasks not recursively by default but individually ...
  // also usefull to prevent hydro if multiple read ins of the same event ...

  /** This function sets the string "id" of the task of a JetScapeTask.
   */
  void SetId(string m_id) { id = m_id; }

  /** This function returns the id of the task of a JetScapeTask.
   */
  const string GetId() const { return id; }

  /** This function returns the mutex of a JetScapeTask.
   */
  const shared_ptr<JetScapeModuleMutex> GetMutex() const { return mutex; }

  /** This function sets the "mutex" of a JetScapeTask.
   */
  void SetMutex(shared_ptr<JetScapeModuleMutex> m_mutex) { mutex = m_mutex; }

  bool GetMultiThread() {return multiThread;}

  void SetMultiThread(bool m_multiThread) {multiThread = m_multiThread;}

  int GetTaskNumber() {return my_task_number_;}

private:
  // can be made sortable to put in correct order or via xml file ...
  vector<shared_ptr<JetScapeTask>> tasks;
  //list<shared_ptr<JetScapeTask>> tasks; // list vs vector any advantage of list?

  /** Flag to set tasks to be inactive. When false, task execution is
   * deactivated i.e. their Exec() function is not run.
   */
  bool active_exec;
  bool multiThread = false;

  string id;
  // if for example a search rather position ... (or always sort with predefined order!?)

  int my_task_number_;
  shared_ptr<JetScapeModuleMutex> mutex;
};

} // end namespace Jetscape

#endif
