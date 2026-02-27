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

#ifndef JETSCAPEMODULEBASE_H
#define JETSCAPEMODULEBASE_H

#include <string>
#include <memory>
#include <random>
#include <map>

#include "JetScapeTask.h"
#include "JetScapeXML.h"
#include "TimeModule.h"
#include "sigslot.h"

#include "cpp17/any.hpp"
#include "cpp17/variant.hpp"

// Maybe change, namespaces for any and varaint ...
using namespace linb;
using namespace mpark;

namespace Jetscape {

class JetScapeWriter;

/**
 * @class JetScapeModuleBase
 * @brief Base class for modules in the Jetscape framework.
 *
 * This class is the base
 * class for all modules in the Jetscape framework. It provides a common
 * interface for all modules, including initialization, execution, and clearing.
 * It also provides a common interface for reading XML files and generating
 * random numbers.
 *
 * Inherits from JetScapeTask and
 * sigslot::has_slots<sigslot::multi_threaded_local>.
 */
class JetScapeModuleBase
    : public JetScapeTask,
      public sigslot::has_slots<sigslot::multi_threaded_local>,
      // public std::enable_shared_from_this<JetScapeModuleBase>,
      public TimeModule {
 public:
  /** @brief Default constructor for JetScapeModuleBase.
   *
   * Sets the XML file name to a default string value.
   */
  JetScapeModuleBase();

  /** @brief Constructor to create a JetScapeModuleBase with a name.
   *
   * @param m_name The name of the module.
   */
  JetScapeModuleBase(string m_name);

  /** @brief Destructor for the JetScapeModuleBase.
   */
  virtual ~JetScapeModuleBase();

  // virtual shared_ptr<JetScapeModuleBase> Clone() const {return nullptr;}

  /**
   * @brief Virtual function to initialize module.
   *
   * Checks if XML is loaded.
   */
  virtual void InitTask();

  /**
   * @brief Virtual function to execute module.
   *
   * Can be overriden by derived classes.
   */
  virtual void ExecuteTask(){};

  /**
   * @brief Virtual function to clear module.
   *
   * Can be overriden by derived classes.
   */
  virtual void ClearTask(){};

  /**
   * @brief Virtual function to execute subtasks of the module.
   *
   * Override execute tasks (at least for now) to keep is-time-stepped
   * handling on module level and not include in JetScapeTask
   */
  virtual void ExecuteTasks() override;

  /**
   * @brief Virtual function to clear module.
   *
   * Can be overriden by derived classes.
   */
  virtual void ClearTasks() override;

  // --------------

  //  Virtual functions to define what is done during the timestep. They can be
  //  overridden by different modules/tasks.

  virtual void CalculateTime(){};

  virtual void CalculateTimeTasks();

  virtual void CalculateTimeTask(){};

  //  Virtual functions to define what is done at the end of the timestep. They
  //  can be overridden by different modules/tasks.

  virtual void ExecTime(){};

  virtual void ExecTimeTasks();

  virtual void ExecTimeTask(){};

  // JP Might be necessary to initialize and finsih things per event if in per
  // time step execution mode ... (to be follow up) see also if one could reuse
  // the InitTask() etc from JetScapeTask since in principle it would be better
  // in the task class since exection is once per event ... well, with this
  // approach the lifetime per event per time step in class hierachry is not
  // very transparent and consistent ... Think about and try a more elegant
  // approach before pulishing the code ...

  virtual void InitPerEvent(){};

  virtual void InitPerEventTasks();

  virtual void FinishPerEvent(){};

  virtual void FinishPerEventTasks();

  // virtual void FinishPerEventTask() {}; // JP: see also in JetScapeTask ...
  // is it really used or would this be the per event way ...

  // --------------

  /**  A (virtual) function to check if module and attached submodules are
  consistently definded as be executed per time step or not. If for example a
  module is neither per event nor per timestep, like jet enerrgyloss modules,
  this check can be avoided by overwriting CheckExec() as an empty function.
   */
  virtual void CheckExec();

  /**  A virtual function to recursively call the CheckExec() function
   */
  virtual void CheckExecs();

  // --------------

  // JP: same can be done with variant if all datatypes are know
  // and put into the varaint definition --> elevated to framework like data
  // types maybe not ideal, to be discussed ...
  virtual any GetHistory() { return 0; }

  // --------------

  /**
   * @brief Sets main XML file name.
   *
   * @param m_name The name of the main XML file.
   */
  void SetXMLMainFileName(string m_name) { xml_main_file_name = m_name; }

  /**
   * @brief Gets main XML file name.
   *
   * @return xml_main_file_name The name of the main XML file.
   */
  string GetXMLMainFileName() { return xml_main_file_name; }

  /**
   * @brief Sets user XML file name.
   *
   * @param m_name The name of the user XML file.
   */
  void SetXMLUserFileName(string m_name) { xml_user_file_name = m_name; }

  /**
   * @brief Gets user XML file name.
   *
   * @return xml_user_file_name The name of the user XML file.
   */
  string GetXMLUserFileName() { return xml_user_file_name; }

  /**
   * @brief Gets the current event number.
   *
   * @return current_event The current event number.
   */
  static int GetCurrentEvent() { return current_event; }

  /**
   * @brief Increments the current event number.
   */
  static void IncrementCurrentEvent() { current_event++; }

  /**
   * @brief Gets shared pointer to random number based on Mersenne-Twister
   * algorithm.
   */
  shared_ptr<std::mt19937> GetMt19937Generator();

  /** Helper functions for XML parsing, wrapping functionality in JetScapeXML:
   */

  /**
   * @brief Retrieves an XML element from the configuration.
   *
   * @param path List of XML element names representing the path to the target
   * element.
   * @param isRequired If true, an error is raised if the element is missing.
   * @return Pointer to the retrieved XMLElement, or nullptr if not found and
   * isRequired is false.
   */
  tinyxml2::XMLElement *GetXMLElement(std::initializer_list<const char *> path,
                                      bool isRequired = true) {
    return JetScapeXML::Instance()->GetElement(path, isRequired);
  }

  /**
   * @brief Retrieves the text content of an XML element.
   *
   * @param path List of XML element names representing the path to the target
   * element.
   * @param isRequired If true, an error is raised if the element is missing.
   * @return The text content of the XML element as a std::string.
   */
  std::string GetXMLElementText(std::initializer_list<const char *> path,
                                bool isRequired = true) {
    return JetScapeXML::Instance()->GetElementText(path, isRequired);
  }

  /**
   * @brief Retrieves an integer value from an XML element.
   *
   * @param path List of XML element names representing the path to the target
   * element.
   * @param isRequired If true, an error is raised if the element is missing.
   * @return The integer value of the XML element.
   */
  int GetXMLElementInt(std::initializer_list<const char *> path,
                       bool isRequired = true) {
    return JetScapeXML::Instance()->GetElementInt(path, isRequired);
  }

  /**
   * @brief Retrieves a double value from an XML element.
   *
   * @param path List of XML element names representing the path to the target
   * element.
   * @param isRequired If true, an error is raised if the element is missing.
   * @return The double value of the XML element.
   */
  double GetXMLElementDouble(std::initializer_list<const char *> path,
                             bool isRequired = true) {
    return JetScapeXML::Instance()->GetElementDouble(path, isRequired);
  }

  // bool GetMultiThread() {return multiThread;}

  // void SetMultiThread(bool m_multiThread) {multiThread = m_multiThread;}

  /**
   * @brief Checks if the module is defined to evolve in time steps.
   * 
   * @return time_stepped True if the module evolves in time steps, false otherwise.
   */
  bool IsTimeStepped() const { return time_stepped; }

  /**
   * @brief Sets whether the module evolves in time steps.
   * 
   * @param m_time_stepped True if the module evolves in time steps, false otherwise.
   */
  void SetTimeStepped(bool m_time_stepped) { time_stepped = m_time_stepped; }

 protected:
  // template<typename T>
  // shared_ptr<T> GetThis() {return
  // std::static_pointer_cast<T>(shared_from_this());}

 private:
  std::string xml_main_file_name;
  std::string xml_user_file_name;
  static int current_event;
  shared_ptr<std::mt19937> mt19937_generator_;
  // bool multiThread = false;

  /// Decides whether the module evolves in time steps
  bool time_stepped;
};

/**
 * @class JetScapeModuleComponentFactory
 * @brief Factory for modules in the Jetscape framework.
 *
 * This class implements a static map (i.e. shared between all instances of
 * JetScapeModuleBase) consisting of a std::string of module names (i.e. name in
 * XML config) and a function that creates an instance of the module. This will
 * allow us to automatically add new modules to Jetscape without modifying the
 * framework classes.
 *
 * Based on: https://stackoverflow.com/a/582456 and
 * https://github.com/alisw/AliPhysics/blob/master/PWG/EMCAL/EMCALtasks/AliEmcalCorrectionComponent.h
 */

/// Template function for creating a new module. Used to register the module.
template <typename T>
shared_ptr<JetScapeModuleBase> createT() {
  return std::make_shared<T>();
}

// Factory to create and keep track of new modules
class JetScapeModuleFactory {
 public:
  virtual ~JetScapeModuleFactory() {}

  typedef std::map<std::string, shared_ptr<JetScapeModuleBase> (*)()> map_type;

  /// Creates an instance of an object based on the name if the name is
  /// registered in the map.
  static shared_ptr<JetScapeModuleBase> createInstance(std::string const &s) {
    map_type::iterator it = getMap()->find(s);
    if (it == getMap()->end()) {
      return 0;
    }
    return it->second();
  }

 protected:
  /// Creates and access the module map
  static map_type *getMap() {
    // We never delete the map (until program termination) because we cannot
    // guarantee correct destruction order
    if (!moduleMap) {
      moduleMap = new map_type;
    }
    return moduleMap;
  }

 private:
  /// Contains the map to all of the modules
  static map_type *moduleMap;
};

/**
 * @class RegisterJetScapeModule
 * @brief Registers Jetscape modules in the factory map
 */
template <typename T>
class RegisterJetScapeModule : public JetScapeModuleFactory {
 public:
  /// Registers the name of the module to map to a function that can create the
  /// module
  RegisterJetScapeModule(std::string const &s) {
    getMap()->insert(std::make_pair(s, &createT<T>));
  }
};

}  // end namespace Jetscape

#endif
