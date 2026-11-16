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

#ifndef SMASHTABULATIONS_H
#define SMASHTABULATIONS_H

#include <cstdlib>
#include <string>

#include "smash/configuration.h"
#include "smash/library.h"

namespace Jetscape {

/**
 * SMASH's particle types, decay modes and tabulated resonance integrals are
 * process-wide: set up once per process, by whichever SMASH module comes
 * first (the initial state or the afterburner).
 *
 * The tabulations are kept in a directory between runs: SMASH stores each
 * table with the hash of its version and of the particle and decay lists,
 * re-tabulates one whose hash differs, and locks the directory while it writes
 * (so jobs can share it). The directory is $SMASH_TABULATIONS if set, else
 * ./smash_tabulations. Tabulating takes several seconds per process; before,
 * the initial-state module deleted the directory at every start.
 */
inline std::string SmashTabulationsPath() {
  const char *env = std::getenv("SMASH_TABULATIONS");
  return (env && env[0] != '\0') ? std::string(env)
                                 : std::string("./smash_tabulations");
}

//! sets up SMASH's particles, decays and tabulations unless done already in
//! this process
inline void InitializeSmashOnce(smash::Configuration &config,
                                const std::string &smash_version) {
  static bool done = false;
  if (done) return;
  smash::initialize_particles_decays_and_tabulations(config, smash_version,
                                                     SmashTabulationsPath());
  done = true;
}

}  // namespace Jetscape

#endif  // SMASHTABULATIONS_H
