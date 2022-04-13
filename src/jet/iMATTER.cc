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

#include "FluidDynamics.h"

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
    pdf = new Pythia8::LHAGrid1( 2212, "20", "../../../../Dropbox/work/pythia8235/share/Pythia8/xmldoc/", &info); /// Assuming its a proton
    /// Locked to a specific PYTHIA directory
    /// should be updated to point to the default PYTHIA directory
    
    
    vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});/// use the same virtuality factor as in the final state calculation
    
}

void iMATTER::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{

    double blurb; // used all the time for testing
    
    
    
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

    double iMATTER::Sudakov(double t_min, double t_max)
    {
        double sudakov = 1;
        
        if ( Current.pid()==qid )
        {
            sudakov = sudakov_Pgg(t_min , t_max)*pow( sudakov_Pqq(t_min , t_max) , nf );
        }
        else if ( Current.pid()==gid )
        {
            sudakov = sudakov_Pqg(t_min, t_max);
        }
        else
        {
            JSWARN << " cannot generate sudakov for pid = " << Current.pid();
        }
        
        return(sudakov);
    }


    double iMATTER::generate_z( Parton p, FourVector r)
    {
        return(0.5);
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


double iMATTER::sudakov_Pgg(double g0, double g1)
{
    
    double Sudakov = 0.5;
    
    return(Sudakov) ;
}

double iMATTER::sudakov_Pqg(double g0, double g1)
{
    
    double Sudakov = 0.5;
    
    return(Sudakov) ;
}

double iMATTER::sudakov_Pqq(double g0, double g1)
{
    
    double Sudakov = 0.5;
    
    return(Sudakov) ;
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
