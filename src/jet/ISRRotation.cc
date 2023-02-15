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
//  ISRRotation.cc
//  
//
//  Created by Abhijit Majumder & Ismail Soudi on 9/13/21.
//

#include "ISRRotation.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include "JetScapeXML.h"
#include "Matter.h"
#include "JetScapeConstants.h"
#include <string>
#include <thread>
#include <cmath>
#include "Pythia8/Pythia.h"
#include "tinyxml2.h"
#include <iostream>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include "helper.h"

// Needed for cubature integration //
// #include <boost/math/quadrature/gauss_kronrod.hpp>

#include "FluidDynamics.h"

#include <gsl/gsl_sf_dilog.h>

#define MAGENTA "\033[35m"

using namespace Jetscape;
#define INTRODUCE_PT 0


ISRRotation::ISRRotation():
TotalMomentumFraction(0.0)
{   
  SetId("ISRRotation");
  VERBOSE(8);
}

ISRRotation::~ISRRotation()
{
  VERBOSE(8);
}

void ISRRotation::InitTask()
{
  JSINFO<<"Intialize ISRRotation ...";

  P_A = GetXMLElementDouble({"Hard","PythiaGun","eCM"})/2.0;  /// Assuming symmetric system
  
  P_B = P_A ; /// assuming symmetric system, rewrite for non-symmetric collision.

  if (!P_A)
  {
      JSWARN << "Initial nucleon energy not found by iMATTER, assuming P_A = 2510 GeV" ;
      P_A = 2510 ;
      P_B = P_A; /// default symmetric assumption
  }
  sP0 = 0.5 / s0 / fmToGeVinv;
  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };

  GSL_RNG = std::shared_ptr<gsl_rng>(gsl_rng_alloc(gsl_rng_default), gsl_rng_free);

}

void ISRRotation::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{

  double blurb; // used all the time for testing

  auto ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  
  FourVector PlusZaxis(0.0,0.0,1.0,1.0);

  // if(  ) return;

  for (int in = 0; in < pIn.size();in++) /// we continue with the loop charade, even though the framework is just giving us one parton
  {
    // JSWARN << "Starting Rotation";

    // Getting The posiiton of the 3DGlauber Hotspots
    if (pIn[in].pz() >= 0) {
      Hotspots = ini->Get_quarks_pos_proj_lab();

    } else {
      Hotspots = ini->Get_quarks_pos_targ_lab();
    }
    // NumHotspots = Hotspots.size() / 3;
    NumHotspots = 3;
    auto Out = pIn[in];
    

    if (true) {
      //  if it's an initial parton we will generate pT 
      if(Out.pstat() < 1000 ){
      
      #if(INTRODUCE_PT == 1)
        // Generating final pT if not already done
        if(!AlreadyGeneratedPTForThisShower) {
          int Num =0;
          SampleAgain:
          Num++;
          pT = GeneratPT();
          // double pTX = std::abs(LatestInitialParton.z()) / 10.;
          // pT = std::array<double,2>{pTX,0.0};

          if(pT[0]*pT[0] + pT[1]*pT[1] >= LatestInitialParton.z() * LatestInitialParton.z() ){
            if(Num < 1000) goto SampleAgain;
            pT[0] = 0.0;
            pT[1] = 0.0;
          }

          (*File1) << "########################### " << std::endl;
          (*File1) << "# Latest Parton momentum " << LatestInitialParton.x() << " "
                  << LatestInitialParton.y() << " " << LatestInitialParton.z() << " "
                  << LatestInitialParton.t() << std::endl;
          DefineRotationMatrix( (Out.pz() >=0 ? 1.:-1.));
          AlreadyGeneratedPTForThisShower = true;
          (*File1) << "Generated new (px, py, new pz) : " << pT[0] << " " << pT[1] << " " << LatestPartonNewPz << std::endl;

          std::ofstream File6;
          File6.open("pTDist.dat",std::ofstream::app);
          File6 << pT[0] << " " << pT[1] << std::endl;
          File6.close();    
        }
      #endif
        
        // New momentum after rotating to get the pT
        FourVector p_Out(Out.px(), Out.py(),Out.pz(), Out.e());
        if(Out.plabel() == LatestPartonLabel_Postive ){
          AddRemenant(Out,LatestPartonLabel_Postive);
          VERBOSE(2) << "LatestPartonLabel_Postive used " << LatestPartonLabel_Postive;
        }
        if(Out.plabel() == LatestPartonLabel_Negative ){
          AddRemenant(Out,LatestPartonLabel_Negative);
          VERBOSE(2) << "LatestPartonLabel_Negative used " << LatestPartonLabel_Negative;
        }

      #if(INTRODUCE_PT == 1)
        RotateVector(p_Out);
        Out.reset_momentum(p_Out);
      #endif
      }

      int lab = Out.plabel();
      if (lab < 0 && lab > -NPartonPerShower){ 
        // Take only ISR partons to perform a shift in pT 
        // Chooses only the inital MPI partons from Pythia who undergo the scattering
          if (Out.pz() >= 0) {
              double OnShellEnergy =sqrt(Out.px() * Out.px() + Out.py() * Out.py() +Out.pz() * Out.pz() +Out.restmass() * Out.restmass());
              ini->CollisionPositiveRotatedMomentum[(-Out.plabel() - 1) / 2] = FourVector(Out.px(), Out.py(), Out.pz(),OnShellEnergy);

          }

          if (Out.pz() < 0) {
              // Out.set_t(0.0);
              double OnShellEnergy = sqrt(Out.px() * Out.px() + Out.py() * Out.py() +Out.pz() * Out.pz() +Out.restmass() * Out.restmass());
              ini->CollisionNegativeRotatedMomentum[(-Out.plabel() - 1) / 2] = FourVector(Out.px(), Out.py(), Out.pz(),OnShellEnergy);
          }
      }



      if ((Out.plabel() == Current_Label || std::abs(Out.pid()) == cid || std::abs(Out.pid()) == bid) && Out.pstat() < 0) {
        
        VERBOSE(1) << MAGENTA << " iMATTER Pushing particlelabel " << Out.plabel() << " status "
               << Out.pstat() << " pid " << Out.pid()
               << " e " << Out.e() << " px " << Out.px() << " py " << Out.py()<< " pz " << Out.pz()
               << " to MCGlauber for subtraction ";
        if (Out.e() < 0) {
          JSWARN << "Energy to subtract is negative !";
          exit(1);
        }
        ini->OutputHardPartonMomentum(Out.e(), Out.px(), Out.py(),
                                      Out.pz(),
                                      (Out.pz() >= 0.0 ? 1 : -1), P_A);
        
      }


      // Rotating Final state partons 
      if (pIn[in].pstat() == 1000) {

        FourVector p_Out(pIn[in].px(), pIn[in].py(), pIn[in].pz(), pIn[in].e());

        int Index = (pIn[in].plabel() - 1) / 2;

        auto CollisionPositive = ini->CollisionPositiveMomentum[Index];
        auto CollisionNegative = ini->CollisionNegativeMomentum[Index];

        auto CollisionPositive1 = ini->CollisionPositiveRotatedMomentum[Index];
        auto CollisionNegative1 = ini->CollisionNegativeRotatedMomentum[Index];

        if (CollisionPositive1.x() == 0 && CollisionNegative1.x() == 0 && CollisionPositive1.y() == 0 && CollisionNegative1.y() == 0 ) { 
          // Don't perform any boost
          return;
        }
        //  else if (CollisionPositive1.x() == CollisionPositive.x() && ) { // if only positive side doesn't need rotation
        //   CollisionPositive1 = ini->CollisionPositiveMomentum[Index];
        // } else if (CollisionNegative1.t() == 0) { // if only negative side doesn't need rotation
        //   CollisionNegative1 = ini->CollisionNegativeMomentum[Index];
        // }

        // double EnergySum = CollisionPositive.t() + CollisionNegative.t();
        double vx = (CollisionPositive.x() + CollisionNegative.x()); // EnergySum;
        double vy = (CollisionPositive.y() + CollisionNegative.y()); // EnergySum;
        double vz = (CollisionPositive.z() + CollisionNegative.z()); // EnergySum;

        // double EnergySum1 = CollisionPositive1.t() + CollisionNegative1.t();
        double vx1 = -(CollisionPositive1.x() + CollisionNegative1.x()); // EnergySum1;
        double vy1 = -(CollisionPositive1.y() + CollisionNegative1.y()); // EnergySum1;
        double vz1 = -(CollisionPositive1.z() + CollisionNegative1.z()); // EnergySum1;

        double DeltapxO2 = (vx + vx1) / 2.; 
        double DeltapyO2 = (vy + vy1) / 2.; 
        double DeltapzO2 = (vz + vz1) / 2.; 

        if (pIn[in].plabel() <= NPartonPerShower && (pIn[in].plabel() - 1) % 2 == 0) {
          ini->Olds = 2.0 * CollisionPositive.t() * CollisionNegative.t() -
                      2.0 * CollisionPositive.x() * CollisionNegative.x() -
                      2.0 * CollisionPositive.y() * CollisionNegative.y() -
                      2.0 * CollisionPositive.z() * CollisionNegative.z();
          ini->Oldt = -2.0 * CollisionPositive.t() * pIn[in].e() +
                      2.0 * CollisionPositive.x() * pIn[in].px() +
                      2.0 * CollisionPositive.y() * pIn[in].py() +
                      2.0 * CollisionPositive.z() * pIn[in].pz();
        }

        if (pIn[in].plabel() <= NPartonPerShower &&
            (pIn[in].plabel() - 1) % 2 == 1) {
          ini->Oldu = -2.0 * CollisionPositive.t() * pIn[in].e() +
                      2.0 * CollisionPositive.x() * pIn[in].px() +
                      2.0 * CollisionPositive.y() * pIn[in].py() +
                      2.0 * CollisionPositive.z() * pIn[in].pz();
        }

        if (!(std::abs(vx - vx1) < 1e-10 && std::abs(vy - vy1) < 1e-10 &&
              std::abs(vz - vz1) < 1e-10)) {
          // p_Out.boost(vx, vy, vz);
          // p_Out.boost(vx1, vy1, vz1);
          double NewPx = pIn[in].px() - DeltapxO2;
          double NewPy = pIn[in].py() - DeltapyO2;
          double NewPz = pIn[in].pz() - DeltapzO2;
          double NewE  = std::sqrt(NewPx * NewPx + NewPy * NewPy + NewPz * NewPz);
          p_Out.Set(NewPx, NewPy, NewPz, NewE);

          if (pIn[in].plabel() <= NPartonPerShower && (pIn[in].plabel() - 1) % 2 == 0) {
            ini->News = 2.0 * CollisionPositive1.t() * CollisionNegative1.t() -
                        2.0 * CollisionPositive1.x() * CollisionNegative1.x() -
                        2.0 * CollisionPositive1.y() * CollisionNegative1.y() -
                        2.0 * CollisionPositive1.z() * CollisionNegative1.z();
            ini->Newt = -2.0 * CollisionPositive1.t() * p_Out.t() +
                        2.0 * CollisionPositive1.x() * p_Out.x() +
                        2.0 * CollisionPositive1.y() * p_Out.y() +
                        2.0 * CollisionPositive1.z() * p_Out.z();
          }

          if (pIn[in].plabel() <= NPartonPerShower && (pIn[in].plabel() - 1) % 2 == 1) {
            ini->Newu = -2.0 * CollisionPositive1.t() * p_Out.t() +
                        2.0 * CollisionPositive1.x() * p_Out.x() +
                        2.0 * CollisionPositive1.y() * p_Out.y() +
                        2.0 * CollisionPositive1.z() * p_Out.z();
          }


          Parton Out = pIn[in];
          Out.reset_momentum(p_Out);

          pOut.push_back(Out);
          return;
        }
      }


      // Changing the status of the radiated stubs
      if(pIn[in].pstat() > 0 && pIn[in].pstat() != 1000){
        Out.set_stat(0);
      }
    }

    SkipRotation:

      pOut.push_back(Out);
      int iout = pOut.size() - 1;
      return;
  }

  // JSINFO << BOLDCYAN << " Moving to next time step " ;
  
  // std::cin >> blurb ;
}
// End of DoEnergyLoss

void ISRRotation::ResetShower(){
  AlreadyGeneratedPTForThisShower = false;
}
void ISRRotation::SetParameters(int LabelOfTheShower_, int NPartonPerShower_, int Current_label_){
  LabelOfTheShower = LabelOfTheShower_;
  NPartonPerShower = NPartonPerShower_;
  Current_Label = Current_label_;
}
double ISRRotation::DisP(double pxGeV, double pyGeV){
    double px = pxGeV * fmToGeVinv;
    double py = pyGeV * fmToGeVinv;


    // double x0=0.4,xy0=-0.4;
    // double x1=-0.1,xy1=1.0;
    // double x2=0.0,xy2=-0.2;
    // double x3=0.4,xy3=0.2;
    // return 2. * std::exp(-2. * s0 * s0 * (px * px + py * py)) * (2. + std::cos(px * (x0 - x1) + py * (xy0 - xy1)) 
    //           + std::cos(px * (x0 - x2) + py * (xy0 - xy2)) + std::cos(px * (x0 - x3) + py * (xy0 - xy3))
    //           + std::cos(px * (x1 - x2) + py * (xy1 - xy2)) + std::cos(px * (x1 - x3) + py * (xy1 - xy3))
    //           + std::cos(px * (x2 - x3) + py * (xy2 - xy3)));
    double Out = 2.0;
    for(int i=0; i<NumHotspots-1; i++){
      for(int j=i+1; j<NumHotspots; j++){
        Out += std::cos(px * (Hotspots[3 * i] - Hotspots[3 * j]) + py * (Hotspots[3 * i + 1] - Hotspots[3 * j + 1]));
      }
    }
    return 2. * std::exp(-2. * s0 * s0 * (px * px + py * py)) * Out;
}

std::array<double, 2> ISRRotation::GeneratPT() {
\
  int NSample = 1000;
  double px = gsl_ran_gaussian(GSL_RNG.get(), sP0);
  double py = gsl_ran_gaussian(GSL_RNG.get(), sP0);

  double Denom = DisP(px, py);

  double SigmaCandidate = 0.1;
  for (size_t i = 0; i < NSample; i++) {
    // //
    double pxCandidate = px + gsl_ran_gaussian(GSL_RNG.get(), SigmaCandidate);
    double pyCandidate = py + gsl_ran_gaussian(GSL_RNG.get(), SigmaCandidate);

    double u = ZeroOneDistribution(*GetMt19937Generator());

    double r = DisP(pxCandidate, pyCandidate) / Denom;
    if (u <= r) {
      px = pxCandidate;
      py = pyCandidate;
      Denom = DisP(px, py);
    }
  }

  return std::array<double, 2>{px, py};
}


void ISRRotation::DefineRotationMatrix(double Dir) {

    double Signb = Dir;
    double px,py,pz;
    if(Dir>=0){
      px = LatestInitialParton_Positive.x() + pT[0];
      py = LatestInitialParton_Positive.y() + pT[1];
      pz = Signb * std::sqrt(LatestInitialParton_Positive.z()*LatestInitialParton_Positive.z() - px*px - py*py);
    } else {
      px = LatestInitialParton_Negative.x() + pT[0];
      py = LatestInitialParton_Negative.y() + pT[1];
      pz = Signb * std::sqrt(LatestInitialParton_Negative.z()*LatestInitialParton_Negative.z() - px*px - py*py);
    }
    LatestPartonNewPz = pz;
    double E = sqrt(px * px + py * py + pz * pz);
    double pt = sqrt(px * px + py * py);

    double cosa, sina;
    if (pt == 0) {
        cosa = 1;
        sina = 0;
    } else {
        cosa = px / pt;
        sina = py / pt;
    }

    double cosb = pz / E;
    double abscosb = std::abs(cosb);
    double sinb = pt / E;

    double sinbHalf2 = 0.5*(1. - abscosb);// sin(b/2)^2
    double sin2a = 2.*cosa*sina; //sin(2a)
    double cosa2 = cosa*cosa;
    double sina2 = sina*sina;

    // RotationMatrix[0] = cosa;
    // RotationMatrix[1] = sina;
    // RotationMatrix[2] = cosb;
    // RotationMatrix[3] = sinb;


    RotationMatrix[0] = (cosa2 * abscosb + sina2); // R11
    RotationMatrix[1] = (sin2a * sinbHalf2); // -R12 & (- R21)
    RotationMatrix[2] = (Signb * cosa * sinb);// R13 & (- R31)

    RotationMatrix[3] = cosa2 + abscosb * sina2; // R22
    RotationMatrix[4] = (Signb * sina * sinb);// R23 & (- R32)
    RotationMatrix[5] = abscosb;// R33

  return;
}


void ISRRotation::RotateVector(FourVector &ToRotate) {
  //     input:  ToRotate, Axis=(px,py,pz) = (0,0,E)_{z}
  //     if i=-1, turn ToRotate in the direction (0,0,E)=>(px,py,pz)

  double &R11 = RotationMatrix[0];
  double &R12 = RotationMatrix[1];
  double &R13 = RotationMatrix[2];
  double &R22 = RotationMatrix[3];
  double &R23 = RotationMatrix[4];
  double &R33 = RotationMatrix[5];

  double wx = ToRotate.x();
  double wy = ToRotate.y();
  double wz = ToRotate.z();
  double e  = ToRotate.t();

  double wx1, wy1, wz1;
  wx1 =  wx * R11 - wy * R12 + wz * R13;
  wy1 = -wx * R12 + wy * R22 + wz * R23;
  wz1 = -wx * R13 - wy * R23 + wz * R33;

  
  ToRotate.Set(wx1,wy1,wz1,e);

  return;
}


void ISRRotation::SetLatestInitialParton(double px, double py, double pz, double E, int Label){
  if(pz >= 0){
    LatestInitialParton_Positive.Set(px,py,pz,E);
    LatestPartonLabel_Postive = Label;
  } else {
    LatestInitialParton_Negative.Set(px,py,pz,E);
    LatestPartonLabel_Negative = Label;
  }

}

void ISRRotation::AddRemenant(Parton &Out,int label){
  auto ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  auto Hard = JetScapeSignalManager::Instance()->GetHardProcessPointer().lock();
  Parton Rem = Out;
  Rem.set_label(label-1);
  double direction = (Rem.pz() >=0 ? 1.:-1.);
  int NHardScatterings = ini->pTHat.size();
  double Pz = (Rem.pz() >=0 ? 1.0:-1.0);

  Rem.reset_momentum(0.25 * Lambda_QCD, 0.25 * Lambda_QCD,Pz,0.0);
  Rem.set_color(Out.anti_color()); 
  Rem.set_anti_color(Out.color());
  if(Rem.pid() != 21 ) Rem.set_id(-Out.pid()); 
  Hard->PushRemnants(Rem);
  return;
}
