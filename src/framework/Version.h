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

#ifndef VERSION_H
#define VERSION_H

#include <string>

namespace Jetscape {

/**
 * @file Version.h
 * @brief Version identifiers for the JETSCAPE framework and X-SCAPE layer.
 *
 * This header exposes two string constants used to identify the
 * X-SCAPE release and the corresponding JETSCAPE compatibility version.
 */

/// Version identifier for the JETSCAPE framework.
const std::string JetScapeVersion = "4.0.2";

/// Version identifier for the X-SCAPE framework.
const std::string XscapeVersion = "2.1";

}  // end namespace Jetscape

#endif
