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

#ifndef PARTONSHOWERGENERATOR_H
#define PARTONSHOWERGENERATOR_H

namespace Jetscape {

class JetEnergyLoss;

/**
 * @class PartonShowerGenerator
 * @brief Base class for generating parton showers.
 *
 * The PartonShowerGenerator class encapsulates the algorithm that
 * evolves partons into a full parton shower. The shower starts from
 * a single initiating hard parton and is recursively evolved in time,
 * with splittings and energy-loss effects modeled by a JetEnergyLoss module.
 *
 * The parton shower is represented as a graph consisting of
 * vertices (splitting points) and partons (edges).
 */
class PartonShowerGenerator {
 public:
  /**
   * @brief Default constructor.
   */
  PartonShowerGenerator(){};

  /**
   * @brief Virtual destructor.
   */
  virtual ~PartonShowerGenerator(){};

  virtual void DoShower(JetEnergyLoss &j){};
  virtual void DoCalculateTime(JetEnergyLoss &j){};

  virtual void DoExecTime(JetEnergyLoss &j){};
  virtual void DoInitPerEvent(JetEnergyLoss &j){};
  virtual void DoFinishPerEvent(JetEnergyLoss &j){};
};

}  // end namespace Jetscape

#endif
