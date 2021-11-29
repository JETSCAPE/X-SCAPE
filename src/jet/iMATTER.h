//
//  iMATTER.h
//  
//
//  Created by Abhijit Majumder on 9/13/21.
//

#ifndef iMATTER_H
#define iMATTER_H

#include "JetEnergyLossModule.h"
#include "Pythia8/Pythia.h"
#include "Matter.h"

using namespace Jetscape;

class iMATTER : public JetEnergyLossModule<iMATTER> 
{
 public:
  
  iMATTER();
  virtual ~iMATTER();

  void Init();
  void DoEnergyLoss(double deltaT,double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w) {}; //funny, should not break if not not overriden !???
    void Dump_pIn_info(int i, vector<Parton> &pIn);

    
    
 protected:
  
  uniform_real_distribution<double> ZeroOneDistribution;
  
};


#endif /* iMATTER_hpp */
