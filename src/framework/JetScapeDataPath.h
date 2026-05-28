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

#ifndef JETSCAPEDATAPATH_H
#define JETSCAPEDATAPATH_H

#include <cstdlib>
#include <string>

namespace Jetscape {

// Root directory holding X-SCAPE read-only data assets (EOS tables, iSS tables,
// LBT tables, 3dMCGlauber inputs, ...). Honors the XSCAPE_DATA_DIR environment
// variable so a shared/installed asset tree can be located from an arbitrary
// per-run working directory. Falls back to "." (current working directory),
// preserving the historical build-directory behavior.
inline std::string GetXSCAPEDataDir() {
  const char *env = std::getenv("XSCAPE_DATA_DIR");
  if (env != nullptr && env[0] != '\0') {
    return std::string(env);
  }
  return std::string(".");
}

// Resolve a relative data-asset path against GetXSCAPEDataDir(). An absolute
// path (leading '/') is returned unchanged.
inline std::string XSCAPEDataPath(const std::string &relative_path) {
  if (!relative_path.empty() && relative_path[0] == '/') {
    return relative_path;
  }
  return GetXSCAPEDataDir() + "/" + relative_path;
}

}  // namespace Jetscape

#endif  // JETSCAPEDATAPATH_H
