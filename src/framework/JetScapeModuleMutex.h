#ifndef JETSCAPEMODULEMUTEX_H
#define JETSCAPEMODULEMUTEX_H

#include <vector>
#include <memory>
#include "JetScapeTask.h"

using namespace std;
using std::shared_ptr;

namespace Jetscape {

/**
 * @class JetScapeModuleMutex
 * @brief Abstract interface for detecting mutual exclusions between JetScape
 * modules.
 *
 * The `JetScapeModuleMutex` class is designed to ensure that certain
 * `JetScapeTask` modules cannot be used together. For example, it may
 * prevent two modules that provide the same functionality from being
 * active simultaneously.
 *
 * Derived classes must implement the `CheckMutex` function, which examines
 * a collection of modules and returns whether they are mutually compatible.
 */
class JetScapeModuleMutex {
 public:
 /**
   * @brief Check if a set of modules violates mutex rules.
   *
   * This pure virtual method should be implemented by derived classes
   * to enforce mutual exclusion constraints between modules.
   *
   * @param modules A vector of shared pointers to `JetScapeTask` modules.
   * @return `true` if the modules are compatible, `false` if they violate
   *         any mutual exclusion rule.
   */
  virtual bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules) = 0;
};

}  // end namespace Jetscape

#endif
