//
//  iMATTER.cc
//  
//
//  Created by Abhijit Majumder on 9/13/21.
//

#include "iMATTER.h"
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
//#include "helper.h"

// Needed for cubature integration //
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "FluidDynamics.h"

#define MAGENTA "\033[35m"

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<iMATTER> iMATTER::reg("iMATTER");

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
    alpha_s = 0.2;
    Q0 = 2.0;

    P_A = GetXMLElementDouble({"Hard","PythiaGun","eCM"})/2.0;  /// Assuming symmetric system
    
    P_B = P_A ; /// assuming symmetric system, rewrite for non-symmetric collision.
    
    if (!P_A)
    {
        JSWARN << "Initial nucleon energy not found by iMATTER, assuming P_A = 2510 GeV" ;
        P_A = 2510 ;
        P_B = P_A; /// default symmetric assumption
    }
    
    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini)
    {

      // If not vacuum case, give warning to add initial state module
      bool in_vac = GetXMLElementInt({"Eloss", "Matter", "in_vac"});
      if (!in_vac)
      {
        JSWARN << "No initial state module for iMATTER, Please check whether you intend to "
                  "add an initial state module.";
      }
    }
    else
    {
        JSINFO << BOLDCYAN << " Initial state module connected to i-MATTER";
    }
    // Initialize random number distribution
    ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };
    
    Pythia8::Info info;
    // Get Pythia data directory //
    std::stringstream pdfpath;
    pdfpath <<  getenv("PYTHIA8DATA") << "/../pdfdata"; // usually PYTHIA8DATA leads to xmldoc but need pdfdata
    std::cerr << "Pythia path: " << pdfpath.str() << std::endl;
    pdf = new Pythia8::LHAGrid1( 2212, "20", pdfpath.str().c_str(), &info); /// Assuming its a proton
    
    vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});/// use the same virtuality factor as in the final state calculation

    SetupIntegration();

}

void iMATTER::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{

    double blurb; // used all the time for testing
    
    FourVector PlusZaxis(0.0,0.0,1.0,1.0);
    
    for (int in=0; in < pIn.size(); in++) /// we continue with the loop charade, even though the framework is just giving us one parton
    {
        if ( pIn[in].plabel()>0 ) return ;
        // i-MATTER only deals with initial state (note the i -> in)
        
       
        
        JSINFO << " " ;
        
        JSINFO << " ********************************************** " ;
        
        JSINFO << " pIn.plabel = " << pIn[in].plabel() << " pIn.pstat = " << pIn[in].pstat() << " pIn.px = " << pIn[in].px() << " pIn.py = " << pIn[in].py() << " pIn.pz = " << pIn[in].pz() << " pIn.e = " << pIn[in].e() ;
        
        if (std::isnan(pIn[in].e()) || std::isnan(pIn[in].px()) ||
        std::isnan(pIn[in].py()) || std::isnan(pIn[in].pz()) ||
        std::isnan(pIn[in].t()) || std::isnan(pIn[in].form_time()))
        {
            JSINFO << BOLDYELLOW << "Parton on entry busted on time step " << time;
            Matter::Dump_pIn_info(in, pIn);
            std::cin >> blurb;
        } /// usual dump routine for naned out partons
        
        if (pIn[in].pstat()>=0) continue ; // Only accept initial state partons
        
        if (pIn[in].pid()==22) continue ; // neglect photons. 
        
    
        //JSINFO << BOLDYELLOW << " pdfvalue = " << extPDF->xf(1,0.2,10);
    
        //std::cin >> blurb;
        JSINFO << BOLDYELLOW <<"iMATTER::DoEnergyLoss at time = "<<time;
        VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;
        double rNum = ZeroOneDistribution(*GetMt19937Generator());
    
        // Setting the Local assignments, position, momentum, velocity and velocity Mod
        double px = pIn[in].px();
        double py = pIn[in].py();
        double pz = pIn[in].pz();
        double e  = pIn[in].e();
        double x_end = pIn[in].x_in().x();
        double y_end = pIn[in].x_in().y();
        double z_end = pIn[in].x_in().z();
        double t_end = pIn[in].time();
        
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
        JSINFO << BOLDYELLOW << " velocityMod = " << velocityMod << " vx = " << velocity[1] << " vy = " << velocity[2] << " vz = " << velocity[3];
        velocity[0] = velocityMod ;
        
        JSINFO << BOLDYELLOW << " end location x = " << x_end << " y = " << y_end << " z = " << z_end << " t = " << t_end ;
        
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
        std::cout<< " at time " << time << " x = " << x << " y= " << y << " and z = " << z << " target-density   = " << density_target << endl ;
        
        std::cout<< " at time " << time << " x = " << x << " y= " << y << " and z = " << z << " projectile-density   = " << density_projectile << endl ;
       
        double t1 = pIn[in].t();
        
        double z_frac = 0.5 ;
        
        FourVector Current_Location(x,y,z,time);
        
        
        JSINFO << BOLDYELLOW << " formation time " << pIn[in].form_time() << " end time = " <<  t_end ;
        
        if (pIn[in].pstat()==-1000) // parton stub from pythia, needs to be reset.
        {
            pIn[in].set_jet_v(velocity);
            
            double max_t = -1*vir_factor*pIn[in].e()*pIn[in].e();
            
            t1 = generate_initial_virt(pIn[in], Current_Location, max_t);
        
            pIn[in].set_t(t1);
            
            pIn[in].set_mean_form_time();
        
            pIn[in].set_form_time( generate_L(std::abs( 2*e/t1 ) ) );
            
            pIn[in].set_stat(-900); // status for an unstable initial state parton moving backward in time.
        
            JSINFO << MAGENTA << " pSTAT = -1000 leads to new virt = " << t1 << " and New formation time = " << pIn[in].form_time() ;

        }
    
        
        Current = pIn[in] ;
               
        //printout_current();
        
        double split_time = t_end + Current.form_time();
        
        JSINFO << BOLDYELLOW << "virtuality = " << Current.t() << " split time = " << split_time ;
        JSINFO << BOLDYELLOW << " new px = " << Current.px() << " py =  " << Current.py() << " pz = " << Current.pz() << " energy = " << Current.e();
        

        
        if (time < split_time)
        {
            
            JSINFO << BOLDRED << "Starting a Split " ;
            
            // Set x2 the current particle's momentum fraction //
            MomentumFractionCurrent = Current.e() / P_A;
            
            z_frac = generate_z( pIn[in], Current_Location) ;
            
            double max_t = pIn[in].t() ;
            
            double t = generate_initial_virt(pIn[in], Current_Location, max_t );
            
            double t2 = -0.1;
            
            double l_perp_2 = (1-z_frac)*std::abs(t1) - std::abs(t)*z_frac*(1-z_frac) - std::abs(t2)*z_frac ;
            
            double l_perp = 0.0 ;
            
            if (l_perp_2 > 0) l_perp = std::sqrt(l_perp_2) ;
            
            double phi = 2*pi*ZeroOneDistribution(*GetMt19937Generator());
            
            double l_perp_x = l_perp*std::cos(phi);
            
            double l_perp_y = l_perp*std::sin(phi);
            
            JSINFO << BOLDYELLOW << " parent t = " << t << " sibling t2 = " << t2 << " resulting l_perp = " << l_perp ;
            
            int pid_a,pid_b;
            
            pid_a = pIn[in].pid();
            pid_b = gid;
            
            if (pIn[in].pid()==gid)
            {
                pid_a = gid;
                pid_b = gid;
            }
            
            if (pIn[in].pid()==qid)
            {
                pid_a = pIn[in].pid();
                pid_b = gid;
            }
            
            FourVector p_Sibling(pIn[in].px()+l_perp_x ,pIn[in].py()+ l_perp_y , pz*(1-z_frac)/z_frac ,  e*(1-z_frac)/z_frac  ) ;
            
            FourVector p_Parent(pIn[in].px()+l_perp_x, pIn[in].py()+l_perp_y , pz/z_frac, e/z_frac) ;
            
            
            // Rotate Back to the current direction 
            ReverseRotateParton(p_Sibling,Current.p_in());

            Parton Sibling(pIn[in].plabel()*10, pid_b , 0 , p_Sibling, Current_Location) ;
            
            Parton Parent(pIn[in].plabel()*10-1, pid_a , -900 , p_Parent, Current_Location) ;
            
            
            
            pOut.push_back(Sibling);
            
            pOut.push_back(Parent);
            
            int iout = pOut.size()-1 ;
            
            pOut[iout].set_t(t);
            
            pOut[iout].set_mean_form_time();
        
            pOut[iout].set_form_time( generate_L(std::abs( 2*e/t/z_frac ) ) );
            
          
        }
        else
        {

            pOut.push_back(pIn[in]);
            int iout = pOut.size()-1;
            
            int sign = 1 ;
            if ( pOut[iout].pz()<0 ) sign=-1;
            
            if (pOut[iout].pstat() == -900) ini->OutputHardPartonMomentum(pOut[iout].e(), pOut[iout].px() , pOut[iout].py() , pOut[iout].pz(), sign );
                
            
        }
        
        
        
        
    }
    
    JSINFO << BOLDCYAN << " Moving to next time step " ;
    
    std::cin >> blurb ;
}
// End of DoEnergyLoss



double iMATTER::generate_initial_virt(Parton p, FourVector location, double max_t)
{
    
    double density_projectile=0.0;
    double density_target=0.0;
    
    if (!ini)
    {
        JSINFO << MAGENTA << "No initial state module, setting the starting location to "
                "0. Make sure to add e.g. trento before PythiaGun.";
    }
    else
    {
        density_target = ini->Get_target_nucleon_density_lab( location.t(), location.x(), location.y(), location.z() );
        density_projectile = ini->Get_projectile_nucleon_density_lab( location.t(), location.x(), location.y(), location.z() );
    }
    
    double r = ZeroOneDistribution(*GetMt19937Generator());
    
   // double t = -1*r*std::abs(max_t);
    // definite negative virtuality
    
    double min_t = 0.5;
    
    
    double t = -1*invert_sudakov(r, min_t , max_t );
    
    
    return(t);
}
    
double iMATTER::invert_sudakov( double value , double min_t, double max_t)
{
    
    
    if ( (value<=0)||(value>=1) )
    {
        JSINFO<< BOLDRED << " error in value passed to sudakov inverter  = " << value ;
        
        throw std::runtime_error(" value needs to be > 0 and < 1") ;
    }
    
    double abs_max_t = std::abs(max_t);
    
    double abs_min_t = std::abs(min_t);
    
    double numer = Sudakov(abs_min_t,abs_max_t);
    
    if (value <= numer) return(min_t);
    
    double lower_t = abs_min_t ;
    
    double upper_t = abs_max_t ;
    
    double abs_t = ( lower_t + upper_t )/2.0 ;
    
    double denom = Sudakov(abs_min_t, abs_t);
    
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
        
        denom = Sudakov(abs_min_t, abs_t);
        
        estimate = numer/denom;
        
        diff = std::abs(value - estimate);
        
        span = std::abs( upper_t - lower_t )/ abs_t ;
        
        loop_counter++;
    }
    
    return(-1.0*abs_t); // returning a space like virtuality t
    
}

double iMATTER::Sudakov(double t1, double t2)
{
    double sudakov = 1;

    double logt1 = (0.5 * alpha_s / pi )*std::log( t1/Q0);
    double logt2 = (0.5 * alpha_s / pi )*std::log( t2/Q0);
    double z_min = 0.0;
    double z_max = 1.-z_min;
    double x2 = MomentumFractionCurrent; 
    
    if ( std::abs(Current.pid()) < 4  ) // A light quark
    {
        // sudakov = sudakov_Pgg(t_min , t_max)*pow( sudakov_Pqq(t_min , t_max) , nf );

        double logP = LogSud_Pqg( z_min, z_max ) + LogSud_Pqq(z_min, z_max);
        double sudakov1 = std::exp(- logt1 * logP) / PDF(Current.pid(),x2,t1);
        double sudakov2 = std::exp(- logt2 * logP) / PDF(Current.pid(),x2,t2);
        sudakov = sudakov2 / sudakov1;
    }
    else if ( Current.pid()==gid )
    {
        // sudakov = sudakov_Pqg(t_min, t_max);
        double logP = LogSud_Pgg( z_min, z_max ) + 2. * nf * LogSud_Pgq( z_min, z_max );
        double sudakov1 = std::exp(- logt1 * logP) / PDF(Current.pid(),x2,t1);
        double sudakov2 = std::exp(- logt2 * logP) / PDF(Current.pid(),x2,t2);
        sudakov = sudakov2 / sudakov1;
        
    }
    else
    {
        JSWARN << " cannot generate sudakov for pid = " << Current.pid();
    }
    
    return(sudakov);
}

double iMATTER::invert_zDist( double value, std::function<double(double,double)> Dist, double t, double denom)
{
    if ( (value<=0)||(value>=1) )
    {
        JSINFO<< BOLDRED << " error in value passed to z Distribution inverter  = " << value ;
        
        throw std::runtime_error(" value needs to be > 0 and < 1") ;
    }
        
    double lower_z = 0.0;
    
    double upper_z = 1.0;
    
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

double iMATTER::generate_z( Parton p, FourVector CurrentLocation)
{
    double r  = ZeroOneDistribution(*GetMt19937Generator());
    double r1 = ZeroOneDistribution(*GetMt19937Generator());
    double zVal = 0.0;

    if ((r > 1) || (r < 0))
    {
        throw std::runtime_error(" error in random number in z *GetMt19937Generator()");
    }

    if( Current.pid() == gid )
    {
        double denomPgg = zDist_Pgg_int(1.0,Current.t());
        // double ratio1 =  full;

        std::function<double(double,double)> FctPgg =  [this](double z_max, double t) { return this->zDist_Pgg_int(z_max,t);};
        
        zVal = invert_zDist(r1,FctPgg,Current.t(),denomPgg);
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
    return Jac * alpha_s / (2.0 * pi) * P_z_gg(z) * PDF(gid,MomentumFractionCurrent / z,t) / z;
}

double iMATTER::zDist_Pgg_int(double z_max, double t){
    double Error;

    // // Change of variables z -> y = 2 arcSin(\sqrt{1-z}) //
    double y_min= 2.0 * std::asin(std::sqrt(1.0 - z_max)), y_max = M_PI;

    std::function<double(double, double)> FctIntegrand = [this](double y, double t) { return this->zDist_Pgg(y,t);};
    double Sudakov = SingleIntegral(FctIntegrand, t, y_min, y_max,Error,1e-1); // The tolerance is taken to be large for now //
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


inline double iMATTER::LogSud_Pgg(double z_min, double z_max)
{
    return ( P_z_gg_int(z_max) - P_z_gg_int(z_min) );
}
inline double iMATTER::LogSud_Pqq(double z_min, double z_max)
{
    return ( P_z_qq_int(z_max) - P_z_qq_int(z_min) );
}
inline double iMATTER::LogSud_Pqg(double z_min, double z_max)
{
    return ( P_z_qg_int(z_max) - P_z_qg_int(z_min) );
}
inline double iMATTER::LogSud_Pgq(double z_min, double z_max)
{
    return ( P_z_qq_int(1. - z_min) - P_z_qq_int(1. - z_max) );
}


void iMATTER::ReverseRotateParton(FourVector &ToRotate, FourVector Axis )
{
    //     input:  ToRotate, Axis=(px,py,pz) = (0,0,E)_{z}
    //     if i=-1, turn ToRotate in the direction (0,0,E)=>(px,py,pz)

    double wx = ToRotate.x();
    double wy = ToRotate.y();
    double wz = ToRotate.z();

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

    double wx1 = wx * cosa * cosb - wy * sina + wz * cosa * sinb;
    double wy1 = wx * sina * cosb + wy * cosa + wz * sina * sinb;
    double wz1 = -wx * sinb + wz * cosb;

    ToRotate.Set(wx1,wy1,wz1,w);

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

    if( Error > epsabs) throw std::runtime_error( "Estimate of the cubature error "+std::to_string(Error)+" is larger than the target "+ std::to_string(epsabs) + "\n Try larger target or using more Gauss points" );
    return kronrod_result;
  }


double iMATTER::PDF(int pid, double z, double t){
    if(pdf->insideBounds(z,t) == 0 ) throw std::runtime_error(" pdf out of bound z = " + std::to_string(z) + " t= " + std::to_string(t) + "\n");
    return pdf->xf(pid,z,t) / z;
}