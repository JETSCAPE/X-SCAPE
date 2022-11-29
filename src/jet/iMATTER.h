/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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
//
//  iMATTER.h
//  
//
//  Created by Abhijit Majumder & Ismail Soudi on 9/13/21.
//

#ifndef iMATTER_H
#define iMATTER_H

#include "ISRRotation.h"
#include "JetEnergyLossModule.h"
#include "Pythia8/Pythia.h"
#include "Matter.h"
#include "HardProcess.h"
#include "gsl/gsl_rng.h"

using namespace Jetscape;

class Matter;

class iMATTER : public JetEnergyLossModule<iMATTER> 
{
 public:
  
   iMATTER();
   virtual ~iMATTER();

   double Q0;

   void Init();
   void DoEnergyLoss(double deltaT,double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
   void WriteTask(weak_ptr<JetScapeWriter> w) {}; //funny, should not break if not not overriden !???
   
   void printout_current();
    
   double alpha_s(double q2);
    
   // std::shared_ptr<InitialState> ini; // temporary pointer to initial state   
   // std::shared_ptr<HardProcess> Hard; // temporary pointer to hard   
   std::shared_ptr<ISRRotation> FinalRotation; 

   // Log of Sudakov Without Virtuality dependence part
   // g -> gg 
   inline double P_z_gg_int(double z);
   inline double P_z_gg(double z);
   inline double LogSud_Pgg(double z_min, double z_max); 
   inline double sudakov_Pgg_Integrand(double t, double y);
   inline double zDist_Pgg(double z, double t);
   double zDist_Pgg_int(double z_max, double t);

   // q -> qg 
   double LogSud_Pqq(double z_min, double z_max); 
   double sudakov_Pqq_Integrand(double t, double y);
   inline double zDist_Pqq(double z, double t);
   double zDist_Pqq_int(double z_max, double t);
   
   // g -> qqbar 
   double P_z_qg_int(double z);
   double P_z_qg(double z);
   double LogSud_Pqg(double z_min, double z_max); 
   double sudakov_Pqg_Integrand(double t, double y);
   inline double zDist_Pqg(double z, double t);
   double zDist_Pqg_int(double z_max, double t);
   
   // q -> gq 
   double P_z_qq_int(double z);
   double P_z_qq(double z);
   double LogSud_Pgq(double z_min, double z_max); 
   double sudakov_Pgq_Integrand(double t, double y);
   inline double zDist_Pgq(double z, double t);
   double zDist_Pgq_int(double z_max, double t);

    
   //double generate_L(double form_time);
   double sudakov_Pgg(double g0, double g1);
   double sud_val_GG(double h0, double h1, double h2, double loc_d, double E1);
   double sud_z_GG(double cg, double cg1, double loc_e, double l_fac, double E2);
   double sudakov_Pqg(double g0, double g1);
   double sud_val_QG(double h0, double h1, double h2, double loc_d, double E1);
   double sud_z_QG(double cg, double cg1, double loc_e, double l_fac, double E2);
   double sudakov_Pqq(double q0, double q1);

   double sud_val_QQ(double h0, double h1, double h2, double loc_d, double E1);
   double sud_z_QQ(double cg, double cg1, double loc_e, double l_fac, double E2);

   double generate_Forward_virt(Parton p, FourVector location,double max_t);
   double generate_initial_virt(Parton p, FourVector location,double max_t);
   double generate_z( Parton p, FourVector CurrentLocation, double tp);
   double generate_L(double form_time);

   double invert_Forward_sudakov( double value , double min_t, double max_t);
   double invert_Backward_sudakov( double value , double min_t, double max_t);

   // Takes function Dist(t,z_min, z_max)
   double invert_zDist( double value, std::function<double(double,double)> Dist, double t, double denom);

   // Probability of evolving backwards from t2 to t1 without branching //
   double Forward_Sudakov(double t1, double t2);
   double Backward_Sudakov(double t1, double t2);
   
    // Rotate Parton to the Axis
    void Rotate( FourVector &ToRotate, FourVector Axis, int icc);
    static void VERBOSE_OUTPUT(std::string st, Parton &pIn);
    
 private:
    
    int LabelOfTheShower, NPartonPerShower = 100000, MAX_COLOR;
    const double z_min_factor = 0.94; // this limits the parent momentum to be P_A/z_min_factor
    double TotalMomentumFraction, TotalMomentum;
    
    Parton Parent,Sibling,Current;
    int Current_Status = 1e8, Current_Label = -1;
    FourVector RotationVector;
    double P_A;
    double RealP_A = 2510;
    double P_B = P_A; /// symmetric system should be overriden in init.
    Pythia8::PDF * pdf;
    double vir_factor=0.25;
    double MomentumFractionCurrent, Maximum_z_frac,z_frac;
    int pid_Sib,pid_Par, Color_Sib = 0, Color_Par = 0, AntiColor_Sib = 0, AntiColor_Par = 0;
    std::string Fpath = "ISR-Partons.dat", Fpath1 = "ISR-Rotation.dat";
    std::ofstream *File, *File1, *File2, *File3, *File4, *File5;

    // Integration setup 
    static const int Nquadrature = 19; // !!! Only takes odd numbers !!! //
    static const int NLegendre = (Nquadrature - 1)/2;
    std::array<double, NLegendre> GaussLegendrePoints, GaussLegendreWeights;
    std::array<double, NLegendre * NLegendre> GaussLegendreDoubleWeights;
    std::array<double, Nquadrature> StieltjesPoints, StieltjesWeights;
    std::array<double, Nquadrature * Nquadrature> StieltjesDoubleWeights;
    void SetupIntegration();
    double DoubleIntegral(std::function<double(double,double)> & Integrand, double a, double b, double a1, double b1, double &Error, double epsabs);
    double SingleIntegral(std::function<double(double,double)> &Integrand, double t, double a, double b, double &Error, double epsabs);
    

    // Get pdf from pythia and check the bounds 
    double PDF(int pid, double z, double t);
    bool RotateShower(Parton& pIn);

    
 protected:
  
  uniform_real_distribution<double> ZeroOneDistribution;
  
};


#endif /* iMATTER_hpp */
