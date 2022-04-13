//
//  iMATTER.h
//  
//
//  Created by Abhijit Majumder on 9/13/21.
//

#ifndef iMATTER_H
#define iMATTER_H

#include "InitialState.h"
#include "JetEnergyLossModule.h"
#include "Pythia8/Pythia.h"
#include "Matter.h"

using namespace Jetscape;

class Matter;

class iMATTER : public JetEnergyLossModule<iMATTER> 
{
 public:
  
  iMATTER();
  virtual ~iMATTER();

  void Init();
  void DoEnergyLoss(double deltaT,double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w) {}; //funny, should not break if not not overriden !???
   
    void printout_current();
    
    
    
    std::shared_ptr<InitialState> ini; // temporary pointer to initial state
    
    
    //double generate_L(double form_time);
    double sudakov_Pgg(double g0, double g1);
    double sud_val_GG(double h0, double h1, double h2, double loc_d, double E1);
    double sud_z_GG(double cg, double cg1, double loc_e, double l_fac, double E2);
    double P_z_gg_int(double cg, double cg1, double loc_e, double cg3, double l_fac, double E2);
    double sudakov_Pqg(double g0, double g1);
    double sud_val_QG(double h0, double h1, double h2, double loc_d, double E1);
    double sud_z_QG(double cg, double cg1, double loc_e, double l_fac, double E2);
    double P_z_qg_int(double cg, double cg1, double loc_e, double cg3, double l_fac, double E2);
    double sudakov_Pqq(double q0, double q1);

    double sud_val_QQ(double h0, double h1, double h2, double loc_d, double E1);
    double sud_z_QQ(double cg, double cg1, double loc_e, double l_fac, double E2);
    double P_z_qq_int(double cg, double cg1, double loc_e, double cg3, double l_fac, double E2);

    double generate_initial_virt(Parton p, FourVector location,double max_t);
    double generate_z( Parton p, FourVector r );
    double generate_L(double form_time);

    double invert_sudakov( double value , double min_t, double max_t);
    double Sudakov(double t_min, double t_max);
    
    
 private:
    
    Parton Parent,Sibling,Current;
    double P_A = 2510;
    double P_B = P_A; /// symmetric system should be overriden in init.
    Pythia8::PDF * pdf;
    double vir_factor=0.5;
    
    
 protected:
  
  uniform_real_distribution<double> ZeroOneDistribution;
  
};


#endif /* iMATTER_hpp */
