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

// REMARK: Old JetScape PSG w/o droplets etc ...
//  pretty much copy of the DoShower in JetEnergyLoss ....

#ifndef PARTONSHOWERGENERATORDEFAULT_H
#define PARTONSHOWERGENERATORDEFAULT_H

#include "PartonShowerGenerator.h"

namespace Jetscape {

class JetEnergyLoss;

/**
 * @class PartonShowerGeneratorDefault
 * @brief Default implementation of the PartonShowerGenerator interface.
 *
 * This class provides a basic implementation of the parton shower generation
 * algorithm. It evolves partons in discrete time steps, applying splitting
 * functions and energy-loss effects as defined by the JetEnergyLoss module.
 */
class PartonShowerGeneratorDefault : public PartonShowerGenerator {
 public:
  PartonShowerGeneratorDefault() : PartonShowerGenerator(){};
  virtual ~PartonShowerGeneratorDefault(){};

  /**
   * @brief Perform the parton shower evolution.
   *
   * This method takes a reference to a JetEnergyLoss module, which
   * provides the splitting functions and energy-loss dynamics.
   * The parton shower starts with the initiating parton provided by
   * the JetEnergyLoss object, and evolves recursively in discrete
   * time steps until the maximum evolution time is reached.
   *
   * The generated shower is recorded in the PartonShower object
   * associated with the JetEnergyLoss module.
   *
   * @param j Reference to the JetEnergyLoss module controlling
   *          the shower evolution.
   */
  virtual void DoShower(JetEnergyLoss &j);
};

}  // end namespace Jetscape

#endif
