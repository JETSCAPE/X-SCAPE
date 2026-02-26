// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef ISRSHOWERPSG_H
#define ISRSHOWERPSG_H

#include "PartonShowerGeneratorDefault.h"
#include "PartonShower.h"
#include <GTL/edge.h>

namespace Jetscape {

class JetEnergyLoss;
/**
 * @brief Initial-state radiation parton shower generator
 *
 * `IsrShowerPSG` implements hooks used by `JetEnergyLoss` to select
 * final-state partons or edges from a `PartonShower` at a specified
 * module time. It inherits the default behavior from
 * `PartonShowerGeneratorDefault` and overrides time-related callbacks.
 */
class IsrShowerPSG : public PartonShowerGeneratorDefault {
 public:
  IsrShowerPSG() : PartonShowerGeneratorDefault() {}
  virtual ~IsrShowerPSG(){};

    /**
     * @brief Update any internal time calculation for the provided module
     *
     * This method is called to allow the shower generator to update or
     * calculate time-related data for the `JetEnergyLoss` module. The
     * default ISR implementation currently only logs at verbose level 3.
     *
     * @param j Reference to the `JetEnergyLoss` module invoking the call.
     */
    virtual void DoCalculateTime(JetEnergyLoss &j);

    /**
     * @brief Execute the generator at the module time step
     *
     * This function inspects the attached `PartonShower` for final edges
     * that have finished before the module's current time and forwards the
     * corresponding partons to the `JetEnergyLoss` instance for processing.
     *
     * @param j Reference to the `JetEnergyLoss` module invoking the call.
     */
    virtual void DoExecTime(JetEnergyLoss &j);

    /**
     * @brief Initialize per-event state
     *
     * Called at the beginning of each event to reset or initialize any
     * per-event flags or state used by the generator.
     *
     * @param j Reference to the `JetEnergyLoss` module invoking the call.
     */
    virtual void DoInitPerEvent(JetEnergyLoss &j);
  // virtual void DoFinishPerEvent(JetEnergyLoss &j); //DEBUG only ...

 private:
    /**
     * @brief Collect final edges from a shower that finished before time t
     *
     * Iterates over the edges of the provided `PartonShower` and appends
     * any edges whose target vertex has no outgoing edges (final) and whose
     * time coordinate is less than `t` to the output vector `vE`.
     *
     * @param pS Shared pointer to the `PartonShower` to inspect.
     * @param t The time threshold; edges ending before this time are
     * considered final for processing.
     * @param[out] vE Vector to which qualifying edges will be appended.
     */
    void GetFinalEdgesForTime(shared_ptr<PartonShower> pS, double t,
                                                        vector<edge> &vE);

    /**
     * @brief Collect final partons from a shower that finished before time t
     *
     * Convenience wrapper that calls `GetFinalEdgesForTime` and then maps
     * the resulting edges to `Parton` instances using the `PartonShower`.
     *
     * @param pS Shared pointer to the `PartonShower` to inspect.
     * @param t The time threshold; partons whose finalizing edges end before
     * time `t` will be appended to `vP`.
     * @param[out] vP Vector to which qualifying `Parton` shared pointers will
     * be appended.
     */
    void GetFinalPartonsForTime(shared_ptr<PartonShower> pS, double t,
                                                            vector<std::shared_ptr<Parton>> &vP);
};

}  // end namespace Jetscape

#endif
