// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "IsrShowerPSG.h"
#include "PartonShower.h"
#include "JetEnergyLoss.h"
#include "JetScapeLogger.h"
#include <GTL/bfs.h>
#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>

#include <iostream>
#include <sstream>

using namespace std;

namespace Jetscape {

/**
 * @brief Collect final edges from a shower that finished before time t
 *
 * See `IsrShowerPSG::GetFinalEdgesForTime` declaration in the header for
 * full details. Iterates over all edges in the provided `PartonShower` and
 * appends edges whose target node has no outgoing edges (final nodes) and
 * whose end time is earlier than `t`.
 *
 * @param pS Shared pointer to the `PartonShower` to inspect.
 * @param t Time threshold; edges with end time < t are selected.
 * @param[out] vE Vector to append matching edges to.
 */
void IsrShowerPSG::GetFinalEdgesForTime(shared_ptr<PartonShower> pS, double t,
                                        vector<edge> &vE) {
  graph::edge_iterator eIt, eEnd;
  for (eIt = pS->edges_begin(), eEnd = pS->edges_end(); eIt != eEnd; ++eIt) {
    if (eIt->target().outdeg() < 1) {
      node nEnd = eIt->target();
      auto vEnd = pS->GetVertex(nEnd);
      double tEnd = vEnd->x_in().t();
      // DEBUG
      // cout<<tEnd<<endl;
      // cout<<nE<<endl;
      if (t > tEnd)
        vE.push_back(*eIt);
    }
  }
}


/**
 * @brief Collect final partons from a shower that finished before time t
 *
 * Wrapper around `GetFinalEdgesForTime` that converts the selected edges
 * into `Parton` shared pointers using the shower's `GetParton` method.
 *
 * @param pS Shared pointer to the `PartonShower` to inspect.
 * @param t Time threshold; partons whose finalizing edges end before `t`
 * will be appended to `vP`.
 * @param[out] vP Vector to append matching `Parton` shared pointers to.
 */
// not the most efficient way via GetFinalEdgesForTime ...
void IsrShowerPSG::GetFinalPartonsForTime(shared_ptr<PartonShower> pS, double t,
                                          vector<std::shared_ptr<Parton>> &vP) {
  vector<edge> vecE;
  GetFinalEdgesForTime(pS, t, vecE);

  for (auto e : vecE)
    vP.push_back(pS->GetParton(e));

  vecE.clear();
}

/**
 * @brief Update internal time calculations for the module
 *
 * Current implementation only writes a verbose log entry. Kept as a
 * separate override to allow future ISR-specific time computations.
 *
 * @param j Reference to the `JetEnergyLoss` module.
 */
void IsrShowerPSG::DoCalculateTime(JetEnergyLoss &j) { VERBOSE(3); }

// REMARK: Not the most elegant way to reuse the standard DoExecTime() in
// JetEnergyLoss ...
//         but seems to work. Think about how to make it more efficient and
//         avoid making things public ... !!!!

/**
 * @brief Execute the ISR shower generator for the current module time step
 *
 * This inspects the attached `PartonShower` for edges that finished prior
 * to the current module time. For each such edge the corresponding
 * `Parton` is queued into `j.pIn` and the start vertex is recorded in
 * `j.vStartVec` before invoking `JetEnergyLoss::DoExecTime` to process
 * them. The temporary containers are cleared afterwards.
 *
 * @param j Reference to the `JetEnergyLoss` module.
 */
void IsrShowerPSG::DoExecTime(JetEnergyLoss &j) {
  double currentTime = j.GetModuleCurrentTime() + j.GetModuleDeltaT();

  VERBOSE(2) << " t = " << currentTime;

  auto pS = j.GetShower();

  vector<edge> vecE;
  GetFinalEdgesForTime(pS, currentTime, vecE);

  if (vecE.size() > 0) {
    for (auto e : vecE) {
      // cout<<e<<" ";

      j.pIn.push_back(*pS->GetParton(e));
      j.vStartVec.push_back(e.target());
    }

    // cout<<endl;
  }

  // JP: Check if here not added an extra deltaT
  j.DoExecTime(j.GetModuleCurrentTime(), j.GetModuleDeltaT());

  j.pIn.clear();
  j.vStartVec.clear();

  vecE.clear();
}

/**
 * @brief Initialize per-event ISR generator state
 *
 * Sets per-event flags required by the `JetEnergyLoss` module. The ISR
 * implementation sets `foundchangedorig` to true and logs at verbose
 * level 2.
 *
 * @param j Reference to the `JetEnergyLoss` module.
 */
void IsrShowerPSG::DoInitPerEvent(JetEnergyLoss &j) {
  VERBOSE(2);

  j.foundchangedorig = true;
}

// for debug ...
/*
void IsrShowerPSG::DoFinishPerEvent(JetEnergyLoss &j)
{
  VERBOSE(2);

  auto pS=j.GetShower();
  pS->SaveAsGV("isr_fsr_"+std::to_string(j.GetMyTaskNumber())+".gv");

  //pS->PrintEdges(false);
  //cout<<&j<<endl;
}
*/

}  // end namespace Jetscape
