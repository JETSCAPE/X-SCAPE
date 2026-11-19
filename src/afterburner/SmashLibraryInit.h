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
// -----------------------------------------------------------------------------
// Process-wide guard for the initialization of the SMASH particle and decay
// tables, shared by all SMASH-based modules (SMASH initial condition,
// SMASHNucleus and SMASH afterburner).
// -----------------------------------------------------------------------------

#ifndef SMASHLIBRARYINIT_H
#define SMASHLIBRARYINIT_H

#include <string>

#include "JetScapeLogger.h"
#include "smash/configuration.h"
#include "smash/input_keys.h"
#include "smash/library.h"

namespace Jetscape {

/**
 * SMASH keeps the ParticleType and DecayModes lists in process-wide globals,
 * which can only be built once (ParticleType::create_type_list throws
 * "Type list was already built!" otherwise). Several SMASH modules can live
 * in the same process, so the initialization is guarded in memory here.
 * The tabulations directory on disk is only a cache; SMASH itself checks its
 * validity via a hash of the SMASH version, particles and decay modes.
 *
 * If the tables were already built by another module, the particles and
 * decaymodes keys are removed from the given configuration and compared to
 * the ones used for the initialization.
 */
inline void InitializeSmashParticlesAndDecaysOnce(
    smash::Configuration &config, const std::string &smash_version,
    const std::string &tabulations_path) {
  static bool initialized = false;
  static std::string particles_used;
  static std::string decaymodes_used;

  if (!initialized) {
    particles_used = config.read(smash::InputKeys::particles);
    decaymodes_used = config.read(smash::InputKeys::decaymodes);
    smash::initialize_particles_decays_and_tabulations(config, smash_version,
                                                       tabulations_path);
    initialized = true;
    return;
  }

  JSINFO << "SMASH particles and decay modes already initialized, reusing them";
  const std::string particles = config.take(smash::InputKeys::particles);
  const std::string decaymodes = config.take(smash::InputKeys::decaymodes);
  if (particles != particles_used || decaymodes != decaymodes_used) {
    JSWARN << "SMASH particles/decaymodes files differ between SMASH modules "
              "in this run. SMASH can only hold one set per process, using the "
              "one of the first initialized SMASH module.";
  }
}

}  // namespace Jetscape

#endif  // SMASHLIBRARYINIT_H
