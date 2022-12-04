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
//  iMATTER.cc
//  
//
//  Created by Abhijit Majumder & Ismail Soudi on 9/13/21.
//

#include "iMATTER.h"
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
#include <boost/filesystem.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include "helper.h"

// Needed for cubature integration //
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "FluidDynamics.h"

#include <gsl/gsl_sf_dilog.h>

#define MAGENTA "\033[35m"

using namespace Jetscape;


double y[] = {0.0,0.0,0.0,0.0};

const Parton initial(0,21,0,0.0,0.0,0.0,0.0,y) ; // a gluon with no momentum at the origin.

iMATTER::iMATTER(): Parent(initial) , Sibling(initial) , Current(initial)
{   
  SetId("iMATTER");
  VERBOSE(8);
}

iMATTER::~iMATTER()
{
  VERBOSE(8);
}

void iMATTER::Init()
{
    JSINFO<<"Intialize iMATTER ...";
    // alpha_s = 0.2;
    Q0 = GetXMLElementDouble({"Eloss", "Matter", "Q0"});
    vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});

    RealP_A = GetXMLElementDouble({"Hard","PythiaGun","eCM"})/2.0;  /// Assuming symmetric system
    
    P_A = RealP_A;

    P_B = P_A ; /// assuming symmetric system, rewrite for non-symmetric collision.
    
    if (!P_A)
    {
        JSWARN << "Initial nucleon energy not found by iMATTER, assuming P_A = 2510 GeV" ;
        P_A = 2510 ;
        P_B = P_A; /// default symmetric assumption
    }

    P_A *= z_min_factor;
    
    // Initialize random number distribution
    ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };

    Pythia8::Info info; 
    // Get Pythia data directory //
    std::stringstream pdfpath;
    pdfpath <<  getenv("PYTHIA8DIR") << "/share/Pythia8/pdfdata"; // usually PYTHIA8DATA leads to xmldoc but need pdfdata
    // boost::filesystem::exists();
    if(!boost::filesystem::exists(pdfpath.str().c_str()))
    {
        pdfpath.clear();
        pdfpath <<  getenv("PYTHIA8DIR") << "/share/pythia8/pdfdata"; // usually PYTHIA8DATA leads to xmldoc but need pdfdata
    }
    if( getenv("PYTHIA8DIR") == NULL){
        throw std::runtime_error("$PYTHIA8DIR environment variable not set, please make sure the variable points to the pythia8 installation folder found using (export PYTHIA8DIR=`pythia8-config --prefix`).");
    }
    if(!boost::filesystem::exists(pdfpath.str().c_str()))
    {
        throw std::runtime_error("iMatter cannot find the pdfdata folder of pythia. Tried looking in " + pdfpath.str());
    }
    pdf = new Pythia8::LHAGrid1( 2212, "20", pdfpath.str().c_str(), &info); /// Assuming its a proton
    

    // Setup the quadrature rules for calculating the z_distribution 
    SetupIntegration();


    FinalRotation = std::make_shared<ISRRotation>();
    FinalRotation->Init();

}

void iMATTER::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
    
    bool IsRotated = false;

    double blurb; // used all the time for testing
    
    FourVector PlusZaxis(0.0,0.0,1.0,1.0);
    auto ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    auto Hard = JetScapeSignalManager::Instance()->GetHardProcessPointer().lock();

    for (int in=0; in < pIn.size(); in++) /// we continue with the loop charade, even though the framework is just giving us one parton
    {  
        if( std::abs(time - GetMaxT()) <= 1e-10 )
        {

            if (  std::abs(pIn[in].pid()) >= cid && std::abs(pIn[in].pid()) != 21 &&  pIn[in].plabel() < 0) { //  if it's a parton not treated by iMatter 
                if( pIn[in].pz() >= 0) {
                    double OnShellEnergy = sqrt(pIn[in].px() * pIn[in].px() +
                                                pIn[in].py() * pIn[in].py() +
                                                pIn[in].pz() * pIn[in].pz() + pIn[in].restmass()*pIn[in].restmass()); 
                    ini->CollisionPositiveMomentum[(-pIn[in].plabel() - 1) / 2] = FourVector(pIn[in].px(),pIn[in].py(),pIn[in].pz(),OnShellEnergy);  
                    
                }

                if( pIn[in].pz() < 0)  {
                    double OnShellEnergy = sqrt(pIn[in].px() * pIn[in].px() +
                                                pIn[in].py() * pIn[in].py() +
                                                pIn[in].pz() * pIn[in].pz() + pIn[in].restmass()*pIn[in].restmass()); 
                    ini->CollisionNegativeMomentum[(-pIn[in].plabel() - 1) / 2] = FourVector(pIn[in].px(),pIn[in].py(),pIn[in].pz(),OnShellEnergy);
                }

                FinalRotation->SetLatestInitialParton(pIn[in].px(),pIn[in].py(),pIn[in].pz(),std::abs(pIn[in].pz()), pIn[in].plabel());

            }

            FinalRotation->SetParameters(LabelOfTheShower,NPartonPerShower, Current_Label);
            FinalRotation->DoEnergyLoss(deltaT,time,Q2,pIn,pOut);
            return;
        }

        if ( pIn[in].plabel()>0 ) return;
        // i-MATTER only deals with initial state (note the i -> in)
        
       
        if (std::isnan(pIn[in].e()) || std::isnan(pIn[in].px()) ||
        std::isnan(pIn[in].py()) || std::isnan(pIn[in].pz()) ||
        std::isnan(pIn[in].t()) || std::isnan(pIn[in].form_time()))
        {
            JSINFO << BOLDYELLOW << "Parton on entry busted on time step " << time;
            Matter::Dump_pIn_info(in, pIn);
            std::cin >> blurb;
        } /// usual dump routine for naned out partons
        
        
        if (  std::abs(pIn[in].pid()) >= cid && std::abs(pIn[in].pid()) != 21 ) { // neglect photons and heavy quarks
            // throw std::runtime_error(" Only light flavors allowed when using iMATTER for now");
            continue ; // neglect heavy quarks. 
        }
        


        if (pIn[in].plabel()>=0) { 
            continue ;
        } // Only accept initial state partons

        // Rotate the parton from the earlier times step 
        IsRotated = RotateShower(pIn[in]);

        // JSINFO << BOLDYELLOW <<" iMATTER::DoEnergyLoss at time = "<<time;
        VERBOSESHOWER(8)<< MAGENTA << " SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;
        double rNum = ZeroOneDistribution(*GetMt19937Generator());
    
        // Setting the Local assignments, position, momentum, velocity and velocity Mod
        double px = pIn[in].px();
        double py = pIn[in].py();
        double pz = pIn[in].pz();
        double e  = pIn[in].e();
        double pT2= px * px + py * py;
        double x_end = pIn[in].x_in().x();
        double y_end = pIn[in].x_in().y();
        double z_end = pIn[in].x_in().z();
        double t_end = pIn[in].time();

        // Set x2 the current particle's momentum fraction //
        double CurrentPlus =  (e + std::abs(pz)) / M_SQRT2;

        Maximum_z_frac = CurrentPlus / (CurrentPlus + Q0);

        double mass = pIn[in].restmass();
        
        double velocity[4];
        double OnShellEnergy = sqrt( px*px + py*py + pz*pz + mass*mass  );
        
        for (int j = 1; j <= 3; j++)
        {
            velocity[j] = pIn[in].p(j) / OnShellEnergy;
        }
        double velocityMod =
        std::sqrt(std::pow(velocity[1], 2) + std::pow(velocity[2], 2) +
                  std::pow(velocity[3], 2));

        if (velocityMod > 1.0 + rounding_error)
        {
            JSINFO << BOLDRED << " tachyonic propagation detected for parton passed from hard scattering, velocity mod = " << velocityMod;
            JSWARN << "velocityMod=" << std::setprecision(20) << velocityMod;
            Matter::Dump_pIn_info(0, pIn);
            //assert(velocityMod < 1.0 + rounding_error);
        }
        // JSINFO << BOLDYELLOW << " velocityMod = " << velocityMod << " vx = " << velocity[1] << " vy = " << velocity[2] << " vz = " << velocity[3];
        velocity[0] = velocityMod ;
        
        // JSINFO << BOLDYELLOW << " end location x = " << x_end << " y = " << y_end << " z = " << z_end << " t = " << t_end ;
        
        // parton position by propagating backward from t_end to time
        double x = x_end + velocity[1]*(time - t_end) ;
        double y = y_end + velocity[2]*(time - t_end) ;
        double z = z_end + velocity[3]*(time - t_end) ;
    
        
        
        // Find the density at the new location x,y,z, time
        double density_projectile=0.0;
        double density_target=0.0;
        
        if (!ini)
        {
            JSINFO << MAGENTA << "No initial state module, setting the starting location to "
                    "0. Make sure to add e.g. trento before PythiaGun.";
        }
        else
        {
            density_target = ini->Get_target_nucleon_density_lab(time, x, y, z);
            density_projectile = ini->Get_projectile_nucleon_density_lab(time, x, y, z);
        }

        //DEBUG:
        // std::cout<< " at time " << time << " x = " << x << " y= " << y << " and z = " << z << " target-density   = " << density_target << endl ;
        // std::cout<< " at time " << time << " x = " << x << " y= " << y << " and z = " << z << " projectile-density   = " << density_projectile << endl ;
       
        double t1 = pIn[in].t();
        
        // double z_frac = 0.5 ;
        
        FourVector Current_Location(x,y,z,time);
        
        
        // JSINFO << BOLDYELLOW << " formation time " << pIn[in].form_time() << " end time = " <<  t_end ;
        Current = pIn[in];
        
        if (pIn[in].pstat()==-1000) // parton stub from pythia, needs to be reset.
        {
            // Save Initial parton's momentum 
            {
                if( pIn[in].pz() >= 0) {
                    double OnShellEnergy = sqrt(pIn[in].px() * pIn[in].px() +
                                                pIn[in].py() * pIn[in].py() +
                                                pIn[in].pz() * pIn[in].pz() + pIn[in].restmass()*pIn[in].restmass()); 
                    ini->CollisionPositiveMomentum[(-pIn[in].plabel() - 1) / 2] = FourVector(pIn[in].px(),pIn[in].py(),pIn[in].pz(),OnShellEnergy);  
                    
                }

                if( pIn[in].pz() < 0)  {
                    double OnShellEnergy = sqrt(pIn[in].px() * pIn[in].px() +
                                                pIn[in].py() * pIn[in].py() +
                                                pIn[in].pz() * pIn[in].pz() + pIn[in].restmass()*pIn[in].restmass()); 
                    ini->CollisionNegativeMomentum[(-pIn[in].plabel() - 1) / 2] = FourVector(pIn[in].px(),pIn[in].py(),pIn[in].pz(),OnShellEnergy);
                }
            }

            pIn[in].set_jet_v(velocity);
            double pThat = ini->pTHat[(-pIn[in].plabel()-1)/2];
            double max_t = vir_factor * pThat * pThat;//*pIn[in].e()*pIn[in].e();



            // Set x2 the current particle's momentum fraction //
            CurrentPlus =  (e + std::abs(pz)) / M_SQRT2;
            MomentumFractionCurrent = CurrentPlus / ( M_SQRT2 * P_A );
            if(pIn[in].pz() >= 0.0) {
                TotalMomentum =  Hard->GetTotalMomentumPositive() - Current.e();
                TotalMomentumFraction =  Hard->GetTotalMomentumFractionPositive() - MomentumFractionCurrent;
                } 
            else {
                TotalMomentum =  Hard->GetTotalMomentumNegative() - Current.e();
                TotalMomentumFraction = Hard->GetTotalMomentumFractionNegative() - MomentumFractionCurrent;
                }
            MomentumFractionCurrent = MomentumFractionCurrent / (1.0 - TotalMomentumFraction);

            t1 = -generate_initial_virt(pIn[in], Current_Location, max_t);
        
            pIn[in].set_t(t1);
            
            pIn[in].set_mean_form_time();
        
            pIn[in].set_form_time( generate_L(std::abs( 2*e/t1 ) ) );
            
            pIn[in].set_stat(-900); // status for an unstable initial state parton moving backward in time.
        
            VERBOSE(1) << MAGENTA << " pSTAT = -1000 leads to new virt = " << t1 << " and New formation time = " << pIn[in].form_time() ;
            VERBOSE(1) << MAGENTA << " TotalMomentumFraction " << TotalMomentumFraction << " MomentumFractionCurrent " << MomentumFractionCurrent << " TotalEnergy " << TotalMomentum + Current.e() << " Current en "<< Current.e();



            // Set x2 the current particle's momentum fraction //
            CurrentPlus =  (e + std::abs(pz)) / M_SQRT2;
            MomentumFractionCurrent = CurrentPlus / ( M_SQRT2 * P_A );
            MomentumFractionCurrent = MomentumFractionCurrent / (1.0 - TotalMomentumFraction);
            
            Maximum_z_frac = CurrentPlus / (CurrentPlus + Q0);
            Current_Label = pIn[in].plabel();
            if(-pIn[in].plabel() >= NPartonPerShower ){
                JSWARN << "NHard partons: " << -pIn[in].plabel() << " allowed: " << NPartonPerShower;
                throw std::runtime_error(" More Hard partons than allowed");
            }

            LabelOfTheShower = (-Current_Label - 1) / 2;
            MAX_COLOR = -Current_Label * Hard->GetMax_ColorPerShower();
            Current_Status = NPartonPerShower;
            // JSINFO << MAGENTA << " MAX_COLOR " << MAX_COLOR;

            FinalRotation->SetLatestInitialParton(pIn[in].px(),pIn[in].py(),pIn[in].pz(),pIn[in].e(),pIn[in].plabel());
            FinalRotation->ResetShower();
            IsRotated = true;

        }
    
        
        Current = pIn[in];

               
        //printout_current();
        
        double split_time = t_end + Current.form_time();
        
        // JSINFO << BOLDYELLOW << "virtuality = " << Current.t() << " split time = " << split_time ;
        // JSINFO << BOLDYELLOW << " new px = " << Current.px() << " py =  " << Current.py() << " pz = " << Current.pz() << " energy = " << Current.e();
        

        // This logic should be worked out better for now
        // We have to do the last splitting two time step before the last step 
        // Because one time step is needed to rotate the shower and then one more for pT generation
        if (time < split_time && pIn[in].plabel() == Current_Label 
            && t1 < - Q0 * Q0 - error  && std::abs(time - GetMaxT()) > 1e-10 && (time + 2 * deltaT > GetMaxT() + rounding_error ) 
            )
        {
            if(abs(pIn[in].px()) >= rounding_error || abs(pIn[in].py()) >= rounding_error ){
                std::cout.precision(15);
                std::cerr.precision(15);
                JSWARN << " Current parton not along the z-axis ";
                JSWARN << " time = "<< time << " MaxT = "<< GetMaxT();
                std::cout << " deltaT = "<< deltaT << " " << (time - deltaT + error) << " MaxT = " << GetMaxT()<< std::endl;
                throw std::runtime_error("iMatter Crashed");
            }
            
            // JSINFO << BOLDRED << "Starting a Split " ;  

            // Direction of Current  //
            double Direction = (pIn[in].pz() >= 0.0 ? 1.0:-1.0);
            
            // FourVector p_Current(px,py,pz,OnShellEnergy);
            // Setting the Local assignments, position, momentum, velocity and velocity Mod
            // OnShellEnergy = sqrt( px*px + py*py + pz*pz);
            px = 0.0;
            py = 0.0;
            pz = Direction * pIn[in].pz();
            e  = pIn[in].e();
            pT2= 0.0;

            // Set x2 the current particle's momentum fraction //
            CurrentPlus =  (e + std::abs(pz)) / M_SQRT2;

            MomentumFractionCurrent = CurrentPlus / (M_SQRT2 * P_A);

            if (pIn[in].pz() >= 0.0) {
              TotalMomentum = Hard->GetTotalMomentumPositive();
              TotalMomentumFraction = Hard->GetTotalMomentumFractionPositive();
            } else {
              TotalMomentum = Hard->GetTotalMomentumNegative();
              TotalMomentumFraction = Hard->GetTotalMomentumFractionNegative();
            }
            
            VERBOSE(2) << BOLDYELLOW << " Current TotalMomentum = " << TotalMomentum;

            if(TotalMomentumFraction >= 1. || TotalMomentum >= P_A) {
                VERBOSE(2) << " Current TotalMomentumFraction = " << TotalMomentumFraction << " is larger than 1";
                VERBOSE(2) << " P_A = " << P_A;
                goto SkipSampling;
            }
            TotalMomentum         -= Current.e();
            TotalMomentumFraction -= MomentumFractionCurrent;
            MomentumFractionCurrent = MomentumFractionCurrent / (1.0 - TotalMomentumFraction);
            Maximum_z_frac = CurrentPlus / (CurrentPlus + Q0);


            // Virtuality of current //
            double max_t = std::abs(t1);
            
            // Virtuality of the parent //
            int NumbSampling =0;
            redovirt:
            NumbSampling++;
            double t = -generate_initial_virt(pIn[in], Current_Location, max_t );
            if(t > max_t){
                if(NumbSampling > 100)
                {
                    JSWARN << " Parent virtuality larger than current"; 
                    JSINFO << " t = " << t << " t1 = " << t1; 
                    JSINFO << " Current.status = " << Current.pstat(); 
                    return;
                    // exit(0);
                }
                goto redovirt;
            }
            
            NumbSampling =0;
            RedoSampling:
            NumbSampling++;
            z_frac = generate_z( pIn[in], Current_Location, std::abs(t)) ;
            double zb_frac = 1.0 - z_frac;

            double max_t2 = zb_frac * zb_frac * max_t / z_frac;

            if(max_t <= Q0 * Q0 || max_t2 <= Q0 * Q0){
                goto SkipSampling;
            }


            // Virtuality of the sibling //
            double t2 = generate_Forward_virt(pIn[in], Current_Location, max_t2 );
            
            double phi = 2*pi*ZeroOneDistribution(*GetMt19937Generator());
            
            double l_perp_x = std::cos(phi);
            
            double l_perp_y = std::sin(phi);
            
            // JSINFO << BOLDYELLOW << " z_frac = " << z_frac << " parent t = " << t << " current t1 = " << t1 << " sibling t2 = " << t2 << " resulting l_perp = " << l_perp_x ;
            
            // We take the proton going along the z-axis
            // The direction is determined from the Direction of the current 
            // double SiblingPT2  = -(zb_frac * (( t1 + zb_frac * pT2) / z_frac - t ) + t2) / z_frac;
            // double SiblingPT2  = zb_frac / z_frac * t - zb_frac * t1 / (z_frac * z_frac) - t2 / z_frac - zb_frac * zb_frac * pT2 / (z_frac * z_frac);
            double SiblingPT2  = -(zb_frac / z_frac / z_frac) * t1 + (zb_frac / z_frac) * t - t2 / z_frac;


            
            double SiblingPT   = 0.0;
            if (SiblingPT2 > 0) SiblingPT = std::sqrt(SiblingPT2);
            else{ 
                if(NumbSampling > 100)
                {
                    JSWARN << " SiblingPT2 negative = " << SiblingPT2 << " z = " << z_frac << " pT2 = " << pT2; 
                    JSINFO << " t = " << t << " t1 = " << t1 << " t2 = " << t2; 
                    JSINFO << " Current.pid() = " << Current.pid() << " pid_Par = " << pid_Par << " pid_Sib = " << pid_Sib; 
                    JSINFO << " E = " << Current.e() << " px = " << Current.px() << " py = " << Current.py() << " pz = " << Current.pz(); 
                    JSINFO << " Current.status = " << Current.pstat(); 
                    JSINFO << " Current.form_time = " << Current.form_time(); 
                    return;
                    // exit(0);
                }
                goto RedoSampling;
            }

            double ParentPlus  = CurrentPlus / z_frac;
            double SiblingPlus = zb_frac * ParentPlus;
            double SiblingMinu = 0.5 * (t2 + SiblingPT2) / SiblingPlus;
            double SiblingEn   = (SiblingPlus + SiblingMinu) / M_SQRT2;
            double SiblingPz   = Direction * (SiblingPlus - SiblingMinu) / M_SQRT2;
            double SiblingPx   = l_perp_x * SiblingPT;
            double SiblingPy   = l_perp_y * SiblingPT;
            double ParentEn    =  e + SiblingEn;
            double ParentPz    = Direction * pz + SiblingPz;
            double ParentPx    = SiblingPx;
            double ParentPy    = SiblingPy;


            FourVector p_Sibling(SiblingPx,SiblingPy,SiblingPz,SiblingEn);
            FourVector p_Parent ( ParentPx, ParentPy, ParentPz, ParentEn);
            
            Current_Status = Current.pstat() + 900;
            RotationVector = p_Parent;

            if( -pIn[in].plabel() < NPartonPerShower ){
                Current_Label = NPartonPerShower * pIn[in].plabel() - 1;  
            }
            else {
                Current_Label = pIn[in].plabel()-2;
            }
            


            double velocity_Sibling[4];
            double velocity_Parent[4];
            velocity_Sibling[1] = SiblingPx / SiblingEn;
            velocity_Sibling[2] = SiblingPy / SiblingEn;
            velocity_Sibling[3] = SiblingPz / SiblingEn;

            velocity_Sibling[0] = std::sqrt(std::pow(velocity_Sibling[1], 2) + std::pow(velocity_Sibling[2], 2) + std::pow(velocity_Sibling[3], 2));

            velocity_Parent[1] = ParentPx / ParentEn;
            velocity_Parent[2] = ParentPy / ParentEn;
            velocity_Parent[3] = ParentPz / ParentEn;

            velocity_Parent[0] = std::sqrt(std::pow(velocity_Parent[1], 2) + std::pow(velocity_Parent[2], 2) + std::pow(velocity_Parent[3], 2));
            Hard->SetNISRShower(Hard->GetNISRShower() + 1);
            // Hard->max_color += 1;
            MAX_COLOR += 1;
            Parton Sibling(Current_Label+1, pid_Sib , 0 + Current_Status, p_Sibling, Current_Location) ;
            Sibling.set_max_color(Hard->GetMax_Color() + Hard->GetNISRShower() * Hard->GetMax_ColorPerShower());
            Sibling.set_color(Color_Sib);
            Sibling.set_anti_color(AntiColor_Sib);
            Sibling.set_jet_v(velocity);
            double Position[4];
            Position[0] = Current_Location.t();
            Position[1] = Current_Location.x();
            Position[2] = Current_Location.y();
            Position[3] = Current_Location.z();
            Sibling.set_x(Position);
            VERBOSE_OUTPUT("Current",Current);
            VERBOSE_OUTPUT("Sibling",Sibling);


            Parton Parent(Current_Label, pid_Par , -900 + Current_Status, p_Parent, Current_Location) ;
            Parent.set_max_color(MAX_COLOR);
            Parent.set_color(Color_Par);
            Parent.set_anti_color(AntiColor_Par);
            Parent.set_jet_v(velocity);
            Parent.set_x(Position);
            VERBOSE_OUTPUT("Parent",Parent);

            double NewTotMomFract = TotalMomentumFraction + ParentPlus / ( M_SQRT2 * (P_A) );
            double NewTotMom      = TotalMomentum + ParentEn;
            VERBOSE(2) << BOLDYELLOW << "NewTotMomFract " << NewTotMomFract << " NewTotMom " << NewTotMom << " P_A "<< P_A;
            if(NewTotMomFract >= 1. || NewTotMom >= P_A ){
                VERBOSE(2) << BOLDYELLOW << "Skipped iMATTER split because no energy is available in the proton eCM";
                goto SkipSampling;
            }
            if(Parent.pz()>=0 ){
                Hard->SetTotalMomentumPositive(NewTotMom);
                Hard->SetTotalMomentumFractionPositive(NewTotMomFract);
            }
            else{
                Hard->SetTotalMomentumNegative(NewTotMom);
                Hard->SetTotalMomentumFractionNegative(NewTotMomFract);
            }

            if (std::isnan(Sibling.e()) || std::isnan(Sibling.px()) ||
            std::isnan(Sibling.py()) || std::isnan(Sibling.pz()) ||
            std::isnan(Sibling.t()) || std::isnan(Sibling.form_time()))
            {
                JSINFO << BOLDYELLOW << "Parton on entry busted on time step " << time;
                Matter::Dump_pIn_info(in, pOut);
                // std::cin >> blurb;
            } /// usual dump routine for naned out partons

            pOut.push_back(pIn[in]);
            pOut.push_back(Sibling);
            pOut.push_back(Parent);
            
            int iout = pOut.size()-1 ;
            
            // pOut[iout].set_t(t);
            
            pOut[iout].set_mean_form_time();
        
            pOut[iout].set_form_time( generate_L(std::abs( 2*ParentEn/t ) ) );

          
        }
        else if (IsRotated)
        {
            SkipSampling:

            pOut.push_back(pIn[in]);
            // int iout = pOut.size()-1;
            
            // int sign = pOut[iout].pz()>=0 ? 1:-1;
            
            // if (pIn[in].plabel() == Current_Label) ini->OutputHardPartonMomentum(pOut[iout].e(), pOut[iout].px() , pOut[iout].py() , pOut[iout].pz(), sign );
                
            
        }
        
        
        
        
    }
    
    // JSINFO << BOLDCYAN << " Moving to next time step " ;
    
    // std::cin >> blurb ;

    return;
}
// End of DoEnergyLoss



bool iMATTER::RotateShower(Parton& pIn){
    if( (pIn.pstat() + 900 == Current_Status || pIn.pstat() == Current_Status) && pIn.pstat() != -1000 ) {
        VERBOSE(1) << "Rotating parton pstat = " << pIn.pstat() 
                   << " plabel = " << pIn.plabel()
                   << " px = " << pIn.px()
                   << " py = " << pIn.py()
                   << " pz = " << pIn.pz();
        double Direction = (pIn.pz() >= 0.0 ? 1.0:-1.0);
        FourVector p_Parton(pIn.px(),pIn.py(),pIn.pz(),0.0);
        Rotate(p_Parton,RotationVector,1);
        pIn.reset_p(p_Parton.x(),p_Parton.y(),Direction * p_Parton.z());
        pIn.set_stat(pIn.pstat() + 1);
        double En = std::sqrt(p_Parton.x() * p_Parton.x() + p_Parton.y() * p_Parton.y() + p_Parton.z() * p_Parton.z());
        double velocity[4];
        velocity[1] =             p_Parton.x() / En;
        velocity[2] =             p_Parton.y() / En;
        velocity[3] = Direction * p_Parton.z() / En;

        velocity[0] = std::sqrt(std::pow(velocity[1], 2) + std::pow(velocity[2], 2) + std::pow(velocity[3], 2));
        pIn.set_jet_v(velocity);
        if(std::abs(p_Parton.x()) < error && std::abs(p_Parton.y()) < error && pIn.pstat() + 899 == Current_Status) 
        {
            FinalRotation->SetLatestInitialParton(pIn.px(),pIn.py(),pIn.pz(),std::abs(pIn.pz()), pIn.plabel());
        }
        VERBOSE(1) << "After Rotation parton pstat = " << pIn.pstat() 
                   << " plabel = " << pIn.plabel()
                   << " px = " << pIn.px()
                   << " py = " << pIn.py()
                   << " pz = " << pIn.pz();
        return true;
    }
    return false;
}

double iMATTER::generate_initial_virt(Parton p, FourVector location, double max_t)
{
    
    // double density_projectile=0.0;
    // double density_target=0.0;
    // auto ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    
    // if (!ini)
    // {
    //     JSINFO << MAGENTA << "No initial state module, setting the starting location to "
    //             "0. Make sure to add e.g. trento before PythiaGun.";
    // }
    // else
    // {
    //     density_target = ini->Get_target_nucleon_density_lab( location.t(), location.x(), location.y(), location.z() );
    //     density_projectile = ini->Get_projectile_nucleon_density_lab( location.t(), location.x(), location.y(), location.z() );
    // }
    
    double r = ZeroOneDistribution(*GetMt19937Generator());
    
   // double t = -1*r*std::abs(max_t);
    // definite negative virtuality
    
    double min_t = Q0 * Q0;
    
    
    double t = Q0 * Q0;
    if(max_t > Q0 * Q0) t = invert_Backward_sudakov(r, min_t , max_t );
    
    
    return(t);
}

double iMATTER::generate_Forward_virt(Parton p, FourVector location, double max_t)
{
    
    // double density_projectile=0.0;
    // double density_target=0.0;
    // auto ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    
    // if (!ini)
    // {
    //     JSINFO << MAGENTA << "No initial state module, setting the starting location to "
    //             "0. Make sure to add e.g. trento before PythiaGun.";
    // }
    // else
    // {
    //     density_target = ini->Get_target_nucleon_density_lab( location.t(), location.x(), location.y(), location.z() );
    //     density_projectile = ini->Get_projectile_nucleon_density_lab( location.t(), location.x(), location.y(), location.z() );
    // }
    
    double r = ZeroOneDistribution(*GetMt19937Generator());
    
   // double t = -1*r*std::abs(max_t);
    // definite negative virtuality
    
    double min_t = Q0 * Q0;
    
    
    double t = Q0 * Q0;
    if(max_t > Q0 * Q0) t = invert_Forward_sudakov(r, min_t , max_t );
    
    
    return(t);
}
    
double iMATTER::invert_Forward_sudakov( double value , double min_t, double max_t)
{
    
    
    if ( (value<=0)||(value>=1) )
    {
        JSINFO<< BOLDRED << " error in value passed to sudakov inverter  = " << value ;
        
        throw std::runtime_error(" value needs to be > 0 and < 1") ;
    }
    
    double abs_max_t = max_t;
    
    double abs_min_t = min_t;
    
    double denom = Forward_Sudakov(abs_min_t,abs_max_t);
    
    if (value <= 1./denom) {
        // Debug
        VERBOSE(1) << "Inverse_Forward_sudakov (value <= denom) value = "<< value << " denom = "<< denom
                    << " abs_min_t = " << abs_min_t << " abs_max_t = " << abs_max_t;
        return(min_t);
    }
    
    double lower_t = abs_min_t ;
    
    double upper_t = abs_max_t ;
    
    double abs_t = ( lower_t + upper_t )/2.0 ;
    
    double numer = Forward_Sudakov(abs_min_t, abs_t);
    
    double estimate = numer/denom;
    
    double diff = std::abs(value - estimate);
    
    double span = std::abs( upper_t - lower_t )/ abs_t ;
    
    int loop_counter = 1;
    
    while ( ( diff > approx )&&( span > error) )
    {
        
        if (loop_counter>1000) throw std::runtime_error(" stuck in infinite loop for finding virtuality ");
                    
        if (estimate<value)
        {
            lower_t = abs_t ;
        }
        else
        {
            upper_t = abs_t ;
        }
        abs_t = ( lower_t + upper_t )/2.0 ;
        
        numer = Forward_Sudakov(abs_min_t, abs_t);
        
        estimate = numer/denom;
        
        diff = std::abs(value - estimate);
        
        span = std::abs( upper_t - lower_t )/ abs_t ;
        
        loop_counter++;
    }

    // Debug
    // (*File) << "# Inverse_Forward_sudakov diff = "<< diff << " estimate = "<< estimate << " denom = "<< denom << " numer = "<< numer
    //         << " abs_min_t = " << abs_min_t << " abs_max_t = " << abs_max_t << " abs_t = " << abs_t << std::endl;
    
    return(abs_t); // returning a time like virtuality t
    
} 
double iMATTER::invert_Backward_sudakov( double value , double min_t, double max_t)
{
    
    
    if ( (value<=0)||(value>=1) )
    {
        JSINFO<< BOLDRED << " error in value passed to sudakov inverter  = " << value ;
        
        throw std::runtime_error(" value needs to be > 0 and < 1") ;
    }
    
    double abs_max_t = std::abs(max_t);
    
    double abs_min_t = std::abs(min_t);
    
    double numer = Backward_Sudakov(abs_min_t,abs_max_t);
    
    if (value <= numer) {
        // Debug
        VERBOSE(1) << "Invert_Backward_sudakov (value <= numer) value = "<< value << " numer = "<< numer
                    << " abs_min_t = " << abs_min_t << " abs_max_t = " << abs_max_t;
        return(min_t);
        }
    
    double lower_t = abs_min_t ;
    
    double upper_t = abs_max_t ;
    
    double abs_t = ( lower_t + upper_t )/2.0 ;
    
    double denom = Backward_Sudakov(abs_min_t, abs_t);
    
    double estimate = numer/denom;
    
    double diff = std::abs(value - estimate);
    
    double span = std::abs( upper_t - lower_t )/ abs_t ;
    
    int loop_counter = 1;
    
    while ( ( diff > approx )&&( span > error) )
    {
        
        if (loop_counter>1000) throw std::runtime_error(" stuck in infinite loop for finding virtuality ");
                    
        if (estimate<value)
        {
            lower_t = abs_t ;
        }
        else
        {
            upper_t = abs_t ;
        }
        abs_t = ( lower_t + upper_t )/2.0 ;
        
        denom = Backward_Sudakov(abs_min_t, abs_t);
        
        estimate = numer/denom;
        
        diff = std::abs(value - estimate);
        
        span = std::abs( upper_t - lower_t )/ abs_t ;
        
        loop_counter++;
    }

    // Debug
    // (*File) << "# Inverse_Backward_sudakov diff = "<< diff << " estimate = "<< estimate << " denom = "<< denom << " numer = "<< numer
    //         << " abs_min_t = " << abs_min_t << " abs_t = " << abs_t << std::endl;
    
    return(abs_t); // returning a time like virtuality t
    
}

double iMATTER::Forward_Sudakov(double t1, double t2)
{
    double sudakov = 0.0, logP;
    // double logt2 = (0.5 * alpha_s / pi ) * std::log( t2/t1);
    // double z_min = 1e-5;//MomentumFractionCurrent;
    // double z_max = 1. - z_min;//Maximum_z_frac;
    double x2 = MomentumFractionCurrent; 
    if ( std::abs(Current.pid()) <= sid  ) // A light quark
    {
        // sudakov = sudakov_Pgg(t_min , t_max)*pow( sudakov_Pqq(t_min , t_max) , nf );
        // double logP = LogSud_Pqg( z_min, z_max ) + LogSud_Pqq(z_min, z_max);
        logP = LogSud_Pqg(t1, t2) + LogSud_Pqq(t1, t2);
        // sudakov = std::exp(- logt2 * logP);
        sudakov = std::exp(-(0.5 * alpha_s((t1 + t2)/2.0) / pi ) * logP);
    }
    else if ( Current.pid()==gid )
    {
        // sudakov = sudakov_Pqg(t_min, t_max);
        // double logP = LogSud_Pgg( z_min, z_max ) + 2. * nf * LogSud_Pgq( z_min, z_max );
        // sudakov = std::exp(- logt2 * logP);
        logP = LogSud_Pgg(t1, t2) + 2. * nf * LogSud_Pqq(t1, t2);
        sudakov = std::exp(-(0.5 * alpha_s((t1 + t2)/2.0) / pi ) * logP);
        
    }
    else
    {
        JSWARN << " cannot generate sudakov for pid = " << Current.pid();
    }
    // (*File) << "# Sud = " << sudakov << " logP = " << logP << " t1 = " << t1 << " t2= " << t2  << " Current.pid()= " << Current.pid() 
    //         << " x2 = " << x2 << " LogSud_Pqq(t1, t2) = " << -(0.5 * alpha_s((t1 + t2)/2.0) / pi ) * LogSud_Pqq(t1, t2) << " LogSud_Pqg(t1, t2) = " << -(0.5 * alpha_s((t1 + t2)/2.0) / pi ) * LogSud_Pqg(t1, t2)
    //         << " LogSud_Pgg(t1, t2) = " << -(0.5 * alpha_s((t1 + t2)/2.0) / pi ) * LogSud_Pgg(t1, t2) << std::endl;
    return(sudakov);
}

double iMATTER::Backward_Sudakov(double t1, double t2)
{
    double x2 = MomentumFractionCurrent; 
    return Forward_Sudakov(t1,t2) / PDF(Current.pid(),x2,t2);
}

double iMATTER::invert_zDist( double value, std::function<double(double,double)> Dist, double t, double denom)
{
    if ( (value<=0)||(value>=1) )
    {
        JSINFO<< BOLDRED << " error in value passed to z Distribution inverter  = " << value ;
        
        throw std::runtime_error(" value needs to be > 0 and < 1") ;
    }
        
    double lower_z = MomentumFractionCurrent;
    
    double upper_z = Maximum_z_frac;
    
    double zVal = ( lower_z + upper_z )/2.0 ;
    
    double numer = Dist(zVal,t);
    
    double estimate = numer/denom;
    
    double diff = std::abs(value - estimate);
    
    double span = std::abs( upper_z - lower_z )/ zVal ;
    
    int loop_counter = 1;
    
    while ( ( diff > approx )&&( span > error) )
    {
        
        if (loop_counter>1000) throw std::runtime_error(" stuck in infinite loop for finding z ");
                    
        if (estimate<value)
        {
            lower_z = zVal ;
        }
        else
        {
            upper_z = zVal ;
        }
        zVal = ( lower_z + upper_z )/2.0 ;
        
        numer = Dist(zVal,t);
        
        estimate = numer/denom;
        
        diff = std::abs(value - estimate);
        
        span = std::abs( upper_z - lower_z ) / zVal ;
        
        loop_counter++;
    }
    
    return zVal; // returning the splitting fraction z
    

    
}

double iMATTER::generate_z( Parton p, FourVector CurrentLocation, double tp)
{
    double r  = ZeroOneDistribution(*GetMt19937Generator());
    double r1 = ZeroOneDistribution(*GetMt19937Generator());
    double zVal = 0.0;
    double ratio = 0.0;
    std::array<double,7> denomEach;
    if ((r > 1) || (r < 0) || (r1 > 1) || (r1 < 0))
    {
        throw std::runtime_error(" error in random number in z *GetMt19937Generator()");
    }

    if( Current.pid() == gid )
    {
        pid_Par = gid;
        denomEach[0] = zDist_Pgg_int(Maximum_z_frac,tp);
        double accum = denomEach[0];

        for (size_t id = 1; id <=3; id++)
        {
            /* code */
            pid_Par = id;
            denomEach[id]  = zDist_Pgq_int(Maximum_z_frac,tp);
            accum += denomEach[id];

            pid_Par = -id;
            denomEach[id+3]  = zDist_Pgq_int(Maximum_z_frac,tp);
            accum += denomEach[id+3];
        }

        ratio =  denomEach[0]/accum;
        if(r <= ratio ){ // g-> gg
            pid_Par = gid;
            pid_Sib = gid;
            if( ZeroOneDistribution(*GetMt19937Generator()) <= 0.5){
                Color_Sib = Current.anti_color();
                Color_Par = Current.color();
                AntiColor_Sib = MAX_COLOR + 1;
                AntiColor_Par = MAX_COLOR + 1;
            } else {
                Color_Sib = MAX_COLOR + 1;
                Color_Par = MAX_COLOR + 1;
                AntiColor_Sib = Current.color();
                AntiColor_Par = Current.anti_color();
            }
            std::function<double(double,double)> Fct =  [this](double z_max, double t) { return this->zDist_Pgg_int(z_max,t);};
            zVal = invert_zDist(r1,Fct,tp,denomEach[0]);
            // Stop searching for process //
            goto FoundzVal;
        } 
        
        for (size_t id = 1; id <=3; id++){ // q -> qg
            ratio +=  denomEach[id]/accum;
            if ( r <= ratio ){
                pid_Par = id;
                pid_Sib = id;
                Color_Sib = Current.anti_color();
                Color_Par = Current.color();
                AntiColor_Sib = 0;
                AntiColor_Par = 0;
                std::function<double(double,double)> Fct =  [this](double z_max, double t) { return this->zDist_Pgq_int(z_max,t);};
                zVal = invert_zDist(r1,Fct,tp,denomEach[id]);
                // Stop searching for process //
                goto FoundzVal;
            }
        }
        for (size_t id = 4; id <=6; id++){ // qbar -> qbar g
            ratio +=  denomEach[id]/accum;
            if ( r <= ratio ){
                pid_Par = -((id % 4) + 1);
                pid_Sib = pid_Par;
                Color_Sib = 0;
                Color_Par = 0;
                AntiColor_Sib = Current.color();
                AntiColor_Par = Current.anti_color();
                std::function<double(double,double)> Fct =  [this](double z_max, double t) { return this->zDist_Pgq_int(z_max,t);};
                zVal = invert_zDist(r1,Fct,tp,denomEach[id]);
                // Stop searching for process //
                goto FoundzVal;
            }
        }
    }
    else if( std::abs(Current.pid()) <= sid )
    {

        std::array<double,2> denomEach;
        pid_Par = Current.pid();
        double denomqq = zDist_Pqq_int(Maximum_z_frac,tp);
        pid_Par = gid;
        double denomqg = zDist_Pqg_int(Maximum_z_frac,tp);
        double accum = denomqq + denomqg;
        double ratio =  denomqq/accum;
        if(r <= ratio ){ // q->qg
            pid_Par = Current.pid();
            pid_Sib = gid;
            if(Current.pid() > 0 ){
                Color_Sib = MAX_COLOR + 1;
                Color_Par = MAX_COLOR + 1;
                AntiColor_Sib = Current.color();
                AntiColor_Par = 0;
            } else {
                Color_Sib = Current.anti_color();
                Color_Par = 0;
                AntiColor_Sib = MAX_COLOR + 1;
                AntiColor_Par = MAX_COLOR + 1;
            }
            std::function<double(double,double)> Fct =  [this](double z_max, double t) { return this->zDist_Pqq_int(z_max,t);};
            zVal = invert_zDist(r1,Fct,tp,denomqq);
        } 
        else{ // g->qqbar
            pid_Par = gid;
            pid_Sib = -Current.pid();
            if(Current.pid() > 0 ){
                Color_Sib = 0;
                Color_Par = Current.color();
                AntiColor_Sib = MAX_COLOR + 1;
                AntiColor_Par = MAX_COLOR + 1;
            } else {
                Color_Sib = MAX_COLOR + 1;
                Color_Par = MAX_COLOR + 1;
                AntiColor_Sib = 0;
                AntiColor_Par = Current.anti_color();
            }
            std::function<double(double,double)> Fct =  [this](double z_max, double t) { return this->zDist_Pqg_int(z_max,t);};
            zVal = invert_zDist(r1,Fct,tp,denomqg);
        }
    }

    FoundzVal:
    if( std::isnan(zVal) || zVal == 0.0 ){

      JSWARN << BOLDRED << "zVal is not a number zVal = " << zVal << " r = " << r << " r1 = " << r1;
      JSWARN << BOLDRED << "MomentumFractionCurrent = " << MomentumFractionCurrent;
      JSINFO << BOLDRED << "pid; Curent  = " << Current.pid() << " Parent = " << pid_Par << " Sibling = " << pid_Sib;
      JSWARN << BOLDRED << "ratio = " << ratio << " t = " << tp;
      JSWARN << BOLDRED << "denomEach[0] = " << denomEach[0] << " denomEach[1] = " << denomEach[1]
                        << " denomEach[2] = " << denomEach[2] << " denomEach[3] = " << denomEach[3] 
                        << " denomEach[4] = " << denomEach[4] << " denomEach[5] = " << denomEach[5]
                        << " denomEach[6] = " << denomEach[6];
      exit(0);
    }
    return zVal;
}

double iMATTER::generate_L(double form_time)
{
  double r, x_low, x_high, x, diff, span, val, arg, norm;

    if (form_time<0) form_time = std::abs(form_time) ;
    // we use a positve formation time to generate a splitting time and then make that negative before return.
    
    
  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*GetMt19937Generator());
  //    r = mtrand1();

  if ((r > 1) || (r < 0))
  {
    throw std::runtime_error(" error in random number in z *GetMt19937Generator()");
  }

  x_low = 0;

  x_high = 8.0 * form_time;
  // the value of x_high is slightly arbitrary, the erf function is more or less zero at this distance.
  // picking 10*form_time will not lead to any different results

  x = (x_low + x_high) / 2.0;

  span = (x_high - x_low) / x_high;

  arg = x / form_time / std::sqrt(pi);

  val = std::erf(arg);

  diff = std::abs(val - r);

  while ((diff > approx) && (span > error)) {
    if ((val - r) > 0.0) {
      x_high = x;
    } else {
      x_low = x;
    }

    x = (x_low + x_high) / 2.0;

    arg = x / form_time / std::sqrt(pi);

    val = std::erf(arg);

    diff = std::abs(val - r);

    span = (x_high - x_low) / x_high;
  }

  //    cout << " random number for dist = " << r << " distance generated = " << x << endl;

    x= -1.0*x;
    
    // For initial
    
  return (x);
}

inline double iMATTER::P_z_gg( double z )
{
    return 2. * Nc * (1. - z*(1.-z)*(1. - z*(1.-z))/(z*(1.-z)) );
}
inline double iMATTER::P_z_qq( double z )
{
    return Cf * ( 1. + z*z )/( 1. - z);
}
inline double iMATTER::P_z_qg( double z )
{
    return 0.5 * (z * z + (1. - z)*(1. - z));
}

inline double iMATTER::zDist_Pgg(double y, double t){

    double siny = std::sin(y / 2.0);
    double z = 1.0 - siny * siny;
    double Jac = std::sqrt(z * (1.0 - z));
    return Jac * alpha_s(t) / (2.0 * pi) * P_z_gg(z) * PDF(gid,MomentumFractionCurrent / z,t) / z;
}
inline double iMATTER::zDist_Pqq(double y, double t){

    double siny = std::sin(y / 2.0);
    double z = 1.0 - siny * siny;
    double Jac = std::sqrt(z * (1.0 - z));
    return Jac * alpha_s(t) / (2.0 * pi) * P_z_qq(z) * PDF(pid_Par,MomentumFractionCurrent / z,t) / z;
}
inline double iMATTER::zDist_Pqg(double y, double t){

    double siny = std::sin(y / 2.0);
    double z = 1.0 - siny * siny;
    double Jac = std::sqrt(z * (1.0 - z));
    return Jac * alpha_s(t) / (2.0 * pi) * P_z_qg(z) * PDF(gid,MomentumFractionCurrent / z,t) / z;
}
inline double iMATTER::zDist_Pgq(double y, double t){

    double siny = std::sin(y / 2.0);
    double z = 1.0 - siny * siny;
    double Jac = std::sqrt(z * (1.0 - z));
    return Jac * alpha_s(t) / (2.0 * pi) * P_z_qq(1.0 - z) * PDF(pid_Par,MomentumFractionCurrent / z,t) / z;
}

double iMATTER::zDist_Pgg_int(double z_max, double t){
    double Error;
    double z_min = MomentumFractionCurrent;
    // // Change of variables z -> y = 2 arcSin(\sqrt{1-z}) //
    double y_min= 2.0 * std::asin(std::sqrt(1.0 - z_max)), y_max = 2.0 * std::asin(std::sqrt(1.0 - z_min));

    std::function<double(double, double)> FctIntegrand = [this](double y, double t) { return this->zDist_Pgg(y,t);};
    double Sudakov = SingleIntegral(FctIntegrand, t, y_min, y_max,Error,1e0); // The tolerance is taken to be large for now //
    return Sudakov;
}
double iMATTER::zDist_Pqq_int(double z_max, double t){
    double Error;
    double z_min = MomentumFractionCurrent;
    // // Change of variables z -> y = 2 arcSin(\sqrt{1-z}) //
    double y_min= 2.0 * std::asin(std::sqrt(1.0 - z_max)), y_max = 2.0 * std::asin(std::sqrt(1.0 - z_min));

    std::function<double(double, double)> FctIntegrand = [this](double y, double t) { return this->zDist_Pqq(y,t);};
    double Sudakov = SingleIntegral(FctIntegrand, t, y_min, y_max,Error,1e0); // The tolerance is taken to be large for now //
    return Sudakov;
}
double iMATTER::zDist_Pgq_int(double z_max, double t){
    double Error;
    double z_min = MomentumFractionCurrent;
    // // Change of variables z -> y = 2 arcSin(\sqrt{1-z}) //
    double y_min= 2.0 * std::asin(std::sqrt(1.0 - z_max)), y_max = 2.0 * std::asin(std::sqrt(1.0 - z_min));

    std::function<double(double, double)> FctIntegrand = [this](double y, double t) { return this->zDist_Pgq(y,t);};
    double Sudakov = SingleIntegral(FctIntegrand, t, y_min, y_max,Error,1e0); // The tolerance is taken to be large for now //
    return Sudakov;
}
double iMATTER::zDist_Pqg_int(double z_max, double t){
    double Error;
    double z_min = MomentumFractionCurrent;
    // // Change of variables z -> y = 2 arcSin(\sqrt{1-z}) //
    double y_min= 2.0 * std::asin(std::sqrt(1.0 - z_max)), y_max = 2.0 * std::asin(std::sqrt(1.0 - z_min));

    std::function<double(double, double)> FctIntegrand = [this](double y, double t) { return this->zDist_Pqg(y,t);};
    double Sudakov = SingleIntegral(FctIntegrand, t, y_min, y_max,Error,1e0); // The tolerance is taken to be large for now //
    return Sudakov;
}

inline double iMATTER::P_z_gg_int( double z)
{  
    return  2.0 * Ca * (z * (2.0 - 0.5 * z + z * z / 3.0) + std::log((1.0 - z) / z) );
}
inline double iMATTER::P_z_qq_int( double z )
{
    return 0.5 * Cf * (z * (1. - 0.5 * z) - 4. * std::log(1. - z));
}
inline double iMATTER::P_z_qg_int( double z )
{
    return 0.5 * (z * (z - z * z - 1.));
}


// inline double iMATTER::LogSud_Pgg(double z_min, double z_max)
// {
//     return ( P_z_gg_int(z_max) - P_z_gg_int(z_min) );
// }
// inline double iMATTER::LogSud_Pqq(double z_min, double z_max)
// {
//     return ( P_z_qq_int(z_max) - P_z_qq_int(z_min) );
// }
// inline double iMATTER::LogSud_Pqg(double z_min, double z_max)
// {
//     return ( P_z_qg_int(z_max) - P_z_qg_int(z_min) );
// }
// inline double iMATTER::LogSud_Pgq(double z_min, double z_max)
// {
//     return ( P_z_qq_int(1. - z_min) - P_z_qq_int(1. - z_max) );
// }

inline double iMATTER::LogSud_Pgg(double t_min, double t_max)
{
    // double t_min2 = t_min * t_min;
    // double t_max2 = t_max * t_max;
    double R1 = t_min / t_max;
    double log1 = std::log(R1);
    // return Ca * (8.0 * std::log(t_max / t_min) + ((t_max - t_min) / 3.0) * (-23.0 - 4.0 * t_min / t_max + 2.0 * t_min2 / t_max2 + 12.0 * std::log(t_max / t_min - 1.0)));
    return Ca * (67. / 9. - 2. / 3. * pi * pi +
                 (R1 * R1) * (-4. / 9. * R1 + 1.0) - 8 * R1 +
                 log1 * (11. / 3. + 2.0 * log1) + 4. * gsl_sf_dilog(R1));
}
inline double iMATTER::LogSud_Pqq(double t_min, double t_max)
{
    // double t_min2 = t_min * t_min;
    // double t_max2 = t_max * t_max;
    double R1 = t_min / t_max;
    double log1 = std::log(R1);
    // return Cf * (0.5 * (t_max - t_min) *(-7.0 + 4.0 * std::log(t_max / t_min - 1.0)) + 3.0 * t_min * std::log(t_max/t_min));
    return Cf * ( 3.0 - pi * pi / 3.0 - 3.0 * R1 + 1.5 * log1 + log1 * log1 + 2.0 * gsl_sf_dilog(R1));
}
inline double iMATTER::LogSud_Pqg(double t_min, double t_max)
{
    // double t_min2 = t_min * t_min;
    // double t_max2 = t_max * t_max;
    double R1 = t_min / t_max;
    double log1 = std::log(R1);
    // return (t_max + t_min ) / 3.0 - t_min2 * (1.0 / t_max - t_min / t_max2) - t_min * std::log(t_max / t_min);
    return -13. / 18. + 2. / 9. * R1 * R1 * R1 - 0.5 * R1 * R1 + R1 - log1 / 3.0;
}
inline double iMATTER::LogSud_Pgq(double t_min, double t_max)
{
    // return Cf * (0.5 * (t_max - t_min) *(-7.0 + 4.0 * std::log(t_max / t_min - 1.0)) + 3.0 * t_min * std::log(t_max/t_min));
    return LogSud_Pqq(t_min, t_max);
}

void iMATTER::Rotate(FourVector &ToRotate, FourVector Axis, int icc)
{
    //     input:  ToRotate, Axis=(px,py,pz) = (0,0,E)_{z}
    //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
    //     if i=-1, turn ToRotate in the direction (0,0,E)=>(px,py,pz)

    double wx = ToRotate.x();
    double wy = ToRotate.y();
    double wz = ToRotate.z();
    double e  = ToRotate.t();

    double px = Axis.x();
    double py = Axis.y();
    double pz = Axis.z();

    double E = sqrt(px * px + py * py + pz * pz);
    double pt = sqrt(px * px + py * py);

    double w = sqrt(wx * wx + wy * wy + wz * wz);
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

    double sinbHalf2 = 0.5*(1. - cosb);// sin(b/2)^2
    double sin2a = 2.*cosa*sina; //sin(2a)
    double cosa2 = cosa*cosa;
    double sina2 = sina*sina;

    double R11 = (cosa2 * cosb + sina2);
    double R12 = (sin2a * sinbHalf2);
    double R13 = (cosa * sinb);

    double R22 = cosa2 + cosb * sina2;
    double R23 = (sina * sinb);
    
    double R33 = cosb;

    double wx1, wy1, wz1;
    if (icc == 1) {
        // wx1 = wx * cosb * cosa + wy * cosb * sina - wz * sinb;
        // wy1 = -wx * sina + wy * cosa;
        // wz1 = wx * sinb * cosa + wy * sinb * sina + wz * cosb;
        wx1 =  wx * R11 - wy * R12 - wz * R13;
        wy1 = -wx * R12 + wy * R22 - wz * R23;
        wz1 =  wx * R13 + wy * R23 + wz * R33;
    }

    else {
        // wx1 = wx * cosb * cosa - wy * sina + wz * cosa * sinb;
        // wy1 = wx * sina * cosb + wy * cosa + wz * sina * sinb;
        // wz1 = -wx * sinb + wz * cosb;


        wx1 =  wx * R11 - wy * R12 + wz * R13;
        wy1 = -wx * R12 + wy * R22 + wz * R23;
        wz1 = -wx * R13 - wy * R23 + wz * R33;
    }

    ToRotate.Set(wx1,wy1,wz1,e);

    return;
}


void iMATTER::printout_current()
{
    JSINFO << BOLDYELLOW << " testing Current pid = " << Current.pid() << " stat = " << Current.pstat() << " label = " << Current.plabel();
    
    JSINFO << BOLDYELLOW << " Current x = " << Current.x_in().x() << " y = " << Current.x_in().y() << " z = " << Current.x_in().z() << " t = " << Current.time() ;
    
    
    JSINFO << BOLDYELLOW << " Current px = " << Current.px() << " py = " << Current.py() << " pz = " << Current.pz() << " E = " << Current.e() ;
    
    
    JSINFO << BOLDYELLOW << " color = " << Current.color() << " acolor = " << Current.anti_color() << " virt = " << Current.t();
    
    JSINFO << BOLDYELLOW << " jet vx = " << Current.jet_v().x() << " vy = " << Current.jet_v().y() << " vz = " << Current.jet_v().z() << " v = " << Current.jet_v().t() ;
    
    JSINFO << BOLDYELLOW << " mean form time = " << Current.mean_form_time() << " form time = " << Current.form_time();
    
    return;
}

void iMATTER::SetupIntegration()
{
    auto helper  = boost::math::quadrature::gauss<double, NLegendre>::abscissa();
    auto helperW = boost::math::quadrature::gauss<double, NLegendre>::weights();
    // Populate the positive values of the Gauss-quadrature
    for (int i = NLegendre / 2; i < NLegendre; i++)
    {
      GaussLegendrePoints[i] = helper[i - NLegendre / 2];
      GaussLegendreWeights[i] = helperW[i - NLegendre / 2];
    }
    // Populate the negative values of the Gauss-quadrature
    for (int i = 0; i < NLegendre / 2; i++)
    {
      GaussLegendrePoints [NLegendre / 2 - 1 - i] = -GaussLegendrePoints[NLegendre / 2 + 1 + i];
      GaussLegendreWeights[NLegendre / 2 - 1 - i] = GaussLegendreWeights[NLegendre / 2 + 1 + i];
    }

    auto helper1  = boost::math::quadrature::gauss_kronrod<double, Nquadrature>::abscissa();
    auto helper1W = boost::math::quadrature::gauss_kronrod<double, Nquadrature>::weights();

    // Populate the positive values of the Stieltjes-quadrature
    for (int i = Nquadrature / 2; i < Nquadrature; i++)
    {
      StieltjesPoints[i]  = helper1[i - Nquadrature / 2];
      StieltjesWeights[i] = helper1W[i - Nquadrature / 2];
    }

    // Populate the positive values of the Stieltjes-quadrature
    for (int i = 0; i < Nquadrature / 2; i++)
    {
      StieltjesPoints [Nquadrature / 2 - 1 - i]  = -StieltjesPoints [Nquadrature / 2 + 1 + i];
      StieltjesWeights[Nquadrature / 2 - 1 - i]  =  StieltjesWeights[Nquadrature / 2 + 1 + i];
    }


    // Populate the Double Weights for faster cubature
    for (size_t i = 0; i < Nquadrature; i++)
    {
      for (size_t j = 0; j < Nquadrature; j++)
      {
        StieltjesDoubleWeights[j + i * Nquadrature] = StieltjesWeights[i] * StieltjesWeights[j];
      }
    }

    for (size_t i = 0; i < NLegendre; i++)
    {
      for (size_t j = 0; j < NLegendre; j++)
      {
        GaussLegendreDoubleWeights[j + i * NLegendre] = GaussLegendreWeights[i] * GaussLegendreWeights[j];
      }
    }

    return;
}

double iMATTER::DoubleIntegral(std::function<double(double,double)> & Integrand, double a, double b, double a1, double b1, double &Error, double epsabs)
{
    // Integrate a 2d std::function \int_{a}^{b} dx \int_{a1}^{b1}dy f(x,y)  //
    // The cubature is done using NLegendre points for Gauss-Legendre and (2*NLegendre + 1) points for Gauss-Stieltjes //
    // and the error is computed from the difference of the two //


    double kronrod_result = 0;
    double gauss_result = 0;
    double Jacobian = 0.25 * (b - a) * (b1 - a1);
    for (size_t i = 1; i < Nquadrature; i+=2)
    {
      double x = 0.5 * (b + a) + 0.5 * (b - a) * StieltjesPoints[i];
      for (size_t j = 1; j < Nquadrature; j+=2)
      {
        double y = 0.5 * (b1 + a1) + 0.5 * (b1 - a1) * StieltjesPoints[j];
        double w = StieltjesDoubleWeights[j + i * Nquadrature];
        double fct = Integrand(x, y);
        kronrod_result += w * fct;
        gauss_result += GaussLegendreDoubleWeights[(j / 2) + NLegendre * (i / 2)] * fct;
      }

      for (size_t j = 0; j < Nquadrature; j+=2)
      {
        double y = 0.5 * (b1 + a1) + 0.5 * (b1 - a1) * StieltjesPoints[j];
        double w = StieltjesDoubleWeights[j + i * Nquadrature];
        double fct = Integrand(x, y);
        kronrod_result += w * fct;
      }
    }

    for (size_t i = 0; i < Nquadrature; i+=2)
    {
      double x = 0.5 * (b + a) + 0.5 * (b - a) * StieltjesPoints[i];
      for (size_t j = 0; j < Nquadrature; j++)
      {
        double y = 0.5 * (b1 + a1) + 0.5 * (b1 - a1) * StieltjesPoints[j];
        double w = StieltjesDoubleWeights[j + i * Nquadrature];
        double fct = Integrand(x, y);
        kronrod_result += w * fct;

      }
    }

    kronrod_result *= Jacobian;
    gauss_result *= Jacobian;

    Error = std::abs(kronrod_result - gauss_result);
    if ( std::abs(kronrod_result) > 1e-15 ) Error /= std::abs(kronrod_result);

    if( Error > epsabs) throw std::runtime_error( "Estimate of the cubature error "+std::to_string(Error)+" is larger than the target "+ std::to_string(epsabs) + "\n Try larger target or using more Gauss points" );
    return kronrod_result;
}


double iMATTER::SingleIntegral(std::function<double(double,double)> &Integrand, double t, double a, double b, double &Error, double epsabs)
  {
    // Integrate a 1d with parameter t std::function \int_{a}^{b} dx f(x,y)  //
    // The cubature is done using NLegendre points for Gauss-Legendre and (2*NLegendre + 1) points for Gauss-Stieltjes //
    // and the error is computed from the difference of the two //

    double kronrod_result = 0;
    double gauss_result = 0;
    double Jacobian = 0.5 * (b - a);
    for (size_t i = 1; i < Nquadrature; i+=2)
    {
      double x = 0.5 * (b + a) + 0.5 * (b - a) * StieltjesPoints[i];
      double w = StieltjesWeights[i];
      double fct = Integrand(x,t);
      kronrod_result += w * fct;
      gauss_result += GaussLegendreWeights[(i / 2)] * fct;
    }

    for (size_t i = 0; i < Nquadrature; i+=2)
    {
      double x = 0.5 * (b + a) + 0.5 * (b - a) * StieltjesPoints[i];
      double w = StieltjesWeights[i];
      double fct = Integrand(x,t);
      kronrod_result += w * fct;
    }

    kronrod_result *= Jacobian;
    gauss_result *= Jacobian;

    Error = std::abs(kronrod_result - gauss_result);
    if ( std::abs(kronrod_result) > 1e-15 ) Error /= std::abs(kronrod_result);

    if( Error > epsabs) throw std::runtime_error( "Estimate of the quadrature error "+std::to_string(Error)+" is larger than the target "+ std::to_string(epsabs) + "\n Try larger target or using more Gauss points" );
    return kronrod_result;
  }


double iMATTER::PDF(int pid, double z, double t){
    if(pdf->insideBounds(z,t) == 0 ) throw std::runtime_error(" pdf out of bound z = " + std::to_string(z) + " t= " + std::to_string(t) + "\n");
    return pdf->xf(pid,z,t) / z;
}



double iMATTER::alpha_s(double q2) {
  double a, L2, q24, c_nf;

  L2 = std::pow(Lambda_QCD, 2);

  q24 = q2 / 4.0;

  c_nf = nf;

  if (q24 > 4.0) {
    c_nf = 4;
  }

  if (q24 > 64.0) {
    c_nf = 5;
  }

  if (q24 > L2) {
    a = 12.0 * pi / (11.0 * Nc - 2.0 * c_nf) / std::log(q24 / L2);
  } else {
    JSWARN << " alpha too large ";
    a = 0.6;
  }

  return a;
}


void iMATTER::VERBOSE_OUTPUT(std::string st, Parton &pIn){
    VERBOSE(2) << st 
               << " pid = " << pIn.pid()
               << " pstat = " << pIn.pstat()
               << " plabel = " << pIn.plabel()
               << " t = " << pIn.t()
               << " e = " << pIn.e()
               << " px = " << pIn.px()
               << " py = " << pIn.py()
               << " pz = " << pIn.pz()
               << " col = " << pIn.color()
               << " acol = " << pIn.anti_color()
               << " maxcol = " << pIn.max_color();
}
