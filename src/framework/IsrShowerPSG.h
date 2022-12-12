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

class IsrShowerPSG : public PartonShowerGeneratorDefault
{
 public:

  IsrShowerPSG() : PartonShowerGeneratorDefault() {}
  virtual ~IsrShowerPSG() {};

  virtual void DoCalculateTime(JetEnergyLoss &j);
  virtual void DoExecTime(JetEnergyLoss &j);
  virtual void DoInitPerEvent(JetEnergyLoss &j);
  //virtual void DoFinishPerEvent(JetEnergyLoss &j); //DEBUG only ...

 private:

   void GetFinalEdgesForTime(shared_ptr<PartonShower> pS, double t, vector<edge> &vE);
   void GetFinalPartonsForTime(shared_ptr<PartonShower> pS, double t, vector<std::shared_ptr<Parton>> &vP);

};

} // end namespace Jetscape

#endif
