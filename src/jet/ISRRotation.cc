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
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "FluidDynamics.h"

#include <gsl/gsl_sf_dilog.h>

#define MAGENTA "\033[35m"

using namespace Jetscape;



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

void ISRRotation::Init()
{
  JSINFO<<"Intialize ISRRotation ...";
  File = new std::ofstream;
  File->open(Fpath.c_str(),std::ofstream::out);
  (*File) << "EventId z_frac x2 color anti_color max_color pid plabel status form_time t E Px Py Pz" << std::endl;
  File->close();    

  ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  if (!ini)
  {

    // If not vacuum case, give warning to add initial state module
    bool in_vac = GetXMLElementInt({"Eloss", "Matter", "in_vac"});
    if (!in_vac)
    {
      JSWARN << "No initial state module for ISRRotation, Please check whether you intend to "
                "add an initial state module.";
    }
  }
  else
  {
      JSINFO << BOLDCYAN << " Initial state module connected to i-MATTER";
  }
  Hard = JetScapeSignalManager::Instance()->GetHardProcessPointer().lock();
  if(Hard)
  {
      JSINFO << BOLDCYAN << " Hard process module connected to i-MATTER";
  }

  sP0 = 0.5 / s0 / fmToGeVinv;
  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };

  GSL_RNG = std::shared_ptr<gsl_rng>(gsl_rng_alloc(gsl_rng_default), gsl_rng_free);

}

void ISRRotation::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{

  double blurb; // used all the time for testing

  
  FourVector PlusZaxis(0.0,0.0,1.0,1.0);

  // if(  ) return;

  for (int in = 0; in < pIn.size();in++) /// we continue with the loop charade, even though the framework is just giving us one parton
  {
    JSWARN << "Starting Rotation";
    if (pIn[in].pz() >= 0) {
      Hotspots = ini->Get_quarks_pos_proj_lab();

    } else {
      Hotspots = ini->Get_quarks_pos_targ_lab();
    }

    NumHotspots = Hotspots.size() / 3;
    auto Out = pIn[in];

    if (true) {



      std::ofstream *File1 = new std::ofstream,*File2 = new std::ofstream,*File3 = new std::ofstream;
      File3->open("ISR-Col.dat", std::ofstream::app);
      File1->open("ISR-Rotation.dat", std::ofstream::app);
      (*File1) << "########################### " << std::endl;
      (*File1) << "# " << time << " " << Out.pid() << " "
               << Out.plabel() << " " << Out.pstat() << " "
               << Out.form_time() << " " << Out.t() << " "
               << Out.e() << " " << Out.px() << " " << Out.py()
               << " " << Out.pz() << " " << std::endl;
      if(Out.plabel() < 1000 ){
        // Generating final pT 
        if(!AlreadyGeneratedPTForThisShower) {
          pT = GeneratPT();

          if(pT[0]*pT[0] + pT[1]*pT[1] >= LatestInitialParton.z() * LatestInitialParton.z() ){
            pT[0] = 0.0;
            pT[1] = 0.0;
          }

          (*File1) << "########################### " << std::endl;
          (*File1) << "# Latest Parton momentum " << LatestInitialParton.x() << " "
                  << LatestInitialParton.y() << " " << LatestInitialParton.z() << " "
                  << LatestInitialParton.t() << std::endl;
          DefineRotationMatrix();
          AlreadyGeneratedPTForThisShower = true;
          (*File1) << "Generated new (px, py) : " << pT[0] << " " << pT[1] << std::endl;
        }

        FourVector p_Out(Out.px(), Out.py(),Out.pz(), Out.e());
        RotateVector(p_Out);
        Out.reset_momentum(p_Out);


        (*File1) << "# IncludingPT " << time << " " << Out.pid() << " "
                << Out.plabel() << " " << Out.pstat() << " "
                << Out.form_time() << " " << Out.t() << " "
                << Out.e() << " " << Out.px() << " " << Out.py()
                << " " << Out.pz() << " " << std::endl;
      }

      int lab = Out.plabel();
      if (lab < 0 && lab > -NPartonPerShower){ // Take only ISR partons to perform a shit in pT // Chooses only the inital MPI partons from Pythia who undergo the scattering
        (*File3) << GetCurrentEvent() << " " << Out.pid() << " "
                << Out.plabel() << " " << Out.form_time() << " "
                << Out.t() << " " << Out.e() << " " << Out.px()
                  << " " << Out.py() << " " << Out.pz() << std::endl;
          // (*File3) << "EventId pid plabel form_Time t E Px Py Pz" << std::endl;
          File3->close();
          if (Out.pz() >= 0) {
            (*File1) << " CollisionPositiveRotatedMomentum Pushing "
                    << Out.plabel() << std::endl;

            if (Out.pstat() <=-900) { // Set momentum to zero to signal that no rotation is needed from this side
              ini->CollisionPositiveRotatedMomentum[(-Out.plabel() - 1) / 2] = FourVector(0, 0, 0, 0);
            } else {
              double OnShellEnergy =sqrt(Out.px() * Out.px() + Out.py() * Out.py() +Out.pz() * Out.pz() +Out.restmass() * Out.restmass());
              ini->CollisionPositiveRotatedMomentum[(-Out.plabel() - 1) / 2] = FourVector(Out.px(), Out.py(), Out.pz(),OnShellEnergy);
            }

          }

          if (Out.pz() < 0) {
            (*File1) << " CollisionNegativeRotatedMomentum Pushing "
                    << Out.plabel() << std::endl;

            if (Out.pstat() <=-900) { // Set momentum to zero to signal that no rotation is needed from this side
              ini->CollisionNegativeRotatedMomentum[(-Out.plabel() - 1) / 2] = FourVector(0, 0, 0, 0);
            } else {
              double OnShellEnergy = sqrt(Out.px() * Out.px() + Out.py() * Out.py() +Out.pz() * Out.pz() +Out.restmass() * Out.restmass());
              ini->CollisionNegativeRotatedMomentum[(-Out.plabel() - 1) / 2] = FourVector(Out.px(), Out.py(), Out.pz(),OnShellEnergy);
            }
          }

          (*File1)
              << ini->CollisionPositiveRotatedMomentum[(-Out.plabel() - 1) /2].x()<< " "
              << ini->CollisionPositiveRotatedMomentum[(-Out.plabel() - 1) /2].y()<< " "
              << ini->CollisionPositiveRotatedMomentum[(-Out.plabel() - 1) /2].z()<< " "
              << ini->CollisionPositiveRotatedMomentum[(-Out.plabel() - 1) /2].t()<< " "
              << "\n ";
          (*File1)
              << ini->CollisionNegativeRotatedMomentum[(-Out.plabel() - 1) /2].x()<< " "
              << ini->CollisionNegativeRotatedMomentum[(-Out.plabel() - 1) /2].y()<< " "
              << ini->CollisionNegativeRotatedMomentum[(-Out.plabel() - 1) /2].z()<< " "
              << ini->CollisionNegativeRotatedMomentum[(-Out.plabel() - 1) /2].t()<< " "
              << "\n ";


            
      }



      if ((Out.plabel() == Current_Label || std::abs(Out.pid()) == cid || std::abs(Out.pid()) == bid) && Out.pstat() < 0) {
        JSWARN << MAGENTA << " Pushing label " << Out.plabel() << " status "
               << Out.pstat() << " pid " << Out.pid()
               << " to MCGlauber ";
        JSINFO << MAGENTA << " e " << Out.e();
        JSINFO << MAGENTA << " px " << Out.px();
        JSINFO << MAGENTA << " py " << Out.py();
        JSINFO << MAGENTA << " pz " << Out.pz();
        if (Out.e() < 0) {
          JSWARN << "Energy to subtract is negative !";
          exit(1);
        }
        ini->OutputHardPartonMomentum(Out.e(), Out.px(), Out.py(),
                                      Out.pz(),
                                      (Out.pz() >= 0.0 ? 1 : -1));
        
      }



      if (pIn[in].pstat() == 1000 ) {

        FourVector p_Out(pIn[in].px(), pIn[in].py(), pIn[in].pz(), pIn[in].e());

        int Index = (pIn[in].plabel() - 1) / 2;

        auto CollisionPositive = ini->CollisionPositiveMomentum[Index];
        auto CollisionNegative = ini->CollisionNegativeMomentum[Index];

        auto CollisionPositive1 = ini->CollisionPositiveRotatedMomentum[Index];
        auto CollisionNegative1 = ini->CollisionNegativeRotatedMomentum[Index];

        if (CollisionPositive1.t() == 0 && CollisionNegative1.t() == 0) { // Don't perform any boost
          File1->close();
          return;
        } else if (CollisionPositive1.t() == 0) { // if only positive side doesn't need rotation
          CollisionPositive1 = ini->CollisionPositiveMomentum[Index];
        } else if (CollisionNegative1.t() == 0) { // if only negative side doesn't need rotation
          CollisionNegative1 = ini->CollisionNegativeMomentum[Index];
        }

        double EnergySum = CollisionPositive.t() + CollisionNegative.t();
        double vx = (CollisionPositive.x() + CollisionNegative.x()) / EnergySum;
        double vy = (CollisionPositive.y() + CollisionNegative.y()) / EnergySum;
        double vz = (CollisionPositive.z() + CollisionNegative.z()) / EnergySum;

        double EnergySum1 = CollisionPositive1.t() + CollisionNegative1.t();
        double vx1 = -(CollisionPositive1.x() + CollisionNegative1.x()) / EnergySum1;
        double vy1 = -(CollisionPositive1.y() + CollisionNegative1.y()) / EnergySum1;
        double vz1 = -(CollisionPositive1.z() + CollisionNegative1.z()) / EnergySum1;

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

        (*File1) << "## " << vx << " " << vy << " " << vz << " " << std::endl;
        (*File1) << "## " << vx1 << " " << vy1 << " " << vz1 << " "
                 << std::endl;

        if (!(std::abs(vx - vx1) < 1e-10 && std::abs(vy - vy1) < 1e-10 &&
              std::abs(vz - vz1) < 1e-10)) {
          p_Out.boost(vx, vy, vz);

          if (pIn[in].plabel() <= NPartonPerShower &&
              (pIn[in].plabel() - 1) % 2 == 0) {
            ini->News = 2.0 * CollisionPositive1.t() * CollisionNegative1.t() -
                        2.0 * CollisionPositive1.x() * CollisionNegative1.x() -
                        2.0 * CollisionPositive1.y() * CollisionNegative1.y() -
                        2.0 * CollisionPositive1.z() * CollisionNegative1.z();
            ini->Newt = -2.0 * CollisionPositive1.t() * pIn[in].e() +
                        2.0 * CollisionPositive1.x() * pIn[in].px() +
                        2.0 * CollisionPositive1.y() * pIn[in].py() +
                        2.0 * CollisionPositive1.z() * pIn[in].pz();
          }

          if (pIn[in].plabel() <= NPartonPerShower &&
              (pIn[in].plabel() - 1) % 2 == 1) {
            ini->Newu = -2.0 * CollisionPositive1.t() * pIn[in].e() +
                        2.0 * CollisionPositive1.x() * pIn[in].px() +
                        2.0 * CollisionPositive1.y() * pIn[in].py() +
                        2.0 * CollisionPositive1.z() * pIn[in].pz();
            File2->open("Mandelstamn.dat", std::ofstream::app);
            (*File2) << GetCurrentEvent() << " " << ini->Olds << " "
                     << ini->Oldt << " " << ini->Oldu << " " << ini->News << " "
                     << ini->Newt << " " << ini->Newu << std::endl;
            File2->close();
          }

          p_Out.boost(vx1, vy1, vz1);

          Parton Out = pIn[in];
          Out.reset_momentum(p_Out);

          pOut.push_back(Out);
          auto out = pOut.size() - 1;
          // File1->precision(16);
          (*File1) << "# " << time << " " << pOut[out].pid() << " "
                   << pOut[out].plabel() << " " << pOut[out].pstat() << " "
                   << pOut[out].form_time() << " " << pOut[out].t() << " "
                   << pOut[out].e() << " " << pOut[out].px() << " "
                   << pOut[out].py() << " " << pOut[out].pz() << " "
                   << std::endl
                   << std::endl;
          return;
        }
      }

      File1->close();
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
void ISRRotation::SetParameters(int LabelOfTheShower_, int NPartonPerShower_){
  LabelOfTheShower = LabelOfTheShower_;
  NPartonPerShower = NPartonPerShower_;
}
double ISRRotation::DisP(double pxGeV, double pyGeV){
    double px = pxGeV * fmToGeVinv;
    double py = pyGeV * fmToGeVinv;
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


void ISRRotation::DefineRotationMatrix() {

    double px = LatestInitialParton.x() + pT[0];
    double py = LatestInitialParton.y() + pT[1];
    double pz = std::sqrt(LatestInitialParton.z()*LatestInitialParton.z() - px*px - py*py);

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
    double sinb = pt / E;


    RotationMatrix[0] = cosa;
    RotationMatrix[1] = sina;
    RotationMatrix[2] = cosb;
    RotationMatrix[3] = sinb;

  return;
}


void ISRRotation::RotateVector(FourVector &ToRotate) {
  //     input:  ToRotate, Axis=(px,py,pz) = (0,0,E)_{z}
  //     if i=-1, turn ToRotate in the direction (0,0,E)=>(px,py,pz)

  double &cosa = RotationMatrix[0];
  double &sina = RotationMatrix[1];
  double &cosb = RotationMatrix[2];
  double &sinb = RotationMatrix[3];

  double wx = ToRotate.x();
  double wy = ToRotate.y();
  double wz = ToRotate.z();
  double e  = ToRotate.t();

  double wx1, wy1, wz1;
  wx1 = wx * cosa * cosb - wy * sina + wz * cosa * sinb;
  wy1 = wx * sina * cosb + wy * cosa + wz * sina * sinb;
  wz1 = -wx * sinb + wz * cosb;

  
  ToRotate.Set(wx1,wy1,wz1,e);

  return;
}


void ISRRotation::SetLatestInitialParton(double px, double py, double pz, double E){
  LatestInitialParton.Set(px,py,pz,E);
}