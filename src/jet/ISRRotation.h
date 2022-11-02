//
//  ISRRotation.h
//  
//
//  Created by Abhijit Majumder & Ismail Soudi on 9/13/21.
//

#ifndef ISRRotation_H
#define ISRRotation_H

#include "HadronizationManager.h"
#include "InitialState.h"
#include "JetEnergyLossModule.h"
#include "Pythia8/Pythia.h"
#include "Matter.h"
// #include "iMATTER.h"
#include "HardProcess.h"
#include "gsl/gsl_rng.h"
#include <memory.h>

using namespace Jetscape;


class ISRRotation : public JetEnergyLossModule<ISRRotation> 
{
 public:
  
   ISRRotation();
   virtual ~ISRRotation();

   void Init();
   void DoEnergyLoss(double deltaT,double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
   void WriteTask(weak_ptr<JetScapeWriter> w) {}; //funny, should not break if not not overriden !???
   void SetParameters(int LabelOfTheShower, int NPartonPerShower, int Current_label_);

   void printout_current();
   void SetLatestInitialParton(double px, double py, double pz, double E, int label);
   void ResetShower();
   void AddRemenant(Parton &Out, int label);
    
   
  //  std::shared_ptr<InitialState> ini; // temporary pointer to initial state   
  //  std::shared_ptr<HardProcess> Hard; // temporary pointer to Hard process   
    
 private:
    // std::shared_ptr<JetEnergyLossModule> *iMatterShower;
    int LabelOfTheShower, NPartonPerShower = 100000;
    bool AlreadyGeneratedPTForThisShower = false;
    std::array<double, 2> pT; double LatestPartonNewPz;
    double s0 = 0.2, sP0;
    double TotalMomentumFraction;

    std::string Fpath = "ISR-PT.dat";
    std::ofstream *File;
    
    int Current_Status = 1e8, Current_Label = -1;
    FourVector RotationVector;
    double P_A = 2510;
    double P_B = P_A; /// symmetric system should be overriden in init.

    // Integration setup 
    std::shared_ptr<gsl_rng>GSL_RNG;

    std::vector<double> Hotspots;
    int NumHotspots;
    std::array<double,2> Sample_PT();
    double DisP(double pxGeV, double pyGeV);
    std::array<double,2> GeneratPT();
    std::array<double,6> RotationMatrix;
    FourVector LatestInitialParton_Positive,LatestInitialParton_Negative;
    int LatestPartonLabel_Postive, LatestPartonLabel_Negative;
    void DefineRotationMatrix(double Dir);
    void RotateVector(FourVector &ToRotate);

    
 protected:
  
  uniform_real_distribution<double> ZeroOneDistribution;
  
};


#endif /* ISRRotation_hpp */
