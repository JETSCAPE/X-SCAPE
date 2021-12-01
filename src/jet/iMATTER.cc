//
//  iMATTER.cc
//  
//
//  Created by Abhijit Majumder on 9/13/21.
//

#include "iMATTER.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "Matter.h"
#include <string>
#include <thread>
#include <cmath>
#include "Pythia8/Pythia.h"
#include "tinyxml2.h"
#include<iostream>
//#include "helper.h"

#include "FluidDynamics.h"

#define MAGENTA "\033[35m"

using namespace Jetscape;

iMATTER::iMATTER()
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

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };
}

void iMATTER::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{


    if ( pIn[0].plabel()>0 ) return ;
    // i-MATTER only deals with initial state (note the i)
    
    JSINFO << " test on pIn[1].plabel = " << pIn[1].plabel() << " pIn.pstat = " << pIn[1].pstat() << " pIn[1].px = " << pIn[1].px() << " pIn[1].py = " << pIn[1].py() << " pIn[1].pz() = " << pIn[1].pz() << " pIn[1].e = " << pIn[1].e() ;
    
    if (std::isnan(pIn[0].e()) || std::isnan(pIn[0].px()) ||
        std::isnan(pIn[0].py()) || std::isnan(pIn[0].pz()) ||
        std::isnan(pIn[0].t()) || std::isnan(pIn[0].form_time())) {
      JSINFO << BOLDYELLOW << "Parton on entry busted on time step " << time;
      Matter::Dump_pIn_info(0, pIn);
        double blurb;
        std::cin >> blurb;
    }
    
    double blurb;
    
    Pythia8::Info info;
    
    Pythia8::PDF* extPDF = new Pythia8::LHAGrid1( 2212, "20", "../../../../Dropbox/work/pythia8235/share/Pythia8/xmldoc/", &info);
    
    //JSINFO << BOLDYELLOW << " pdfvalue = " << extPDF->xf(1,0.2,10);
    
    //std::cin >> blurb;
    
  JSINFO << BOLDYELLOW <<"iMATTER::DoEnergyLoss at time = "<<time;
  VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;
 
  double rNum = ZeroOneDistribution(*GetMt19937Generator());

  //DEBUG:
    cout<< rNum << endl ;

  if (rNum<0.05 && std::abs(time)<3)
    {
      //DEBUG:
        JSINFO << MAGENTA <<"iMATTER module ..."<<time<<" "<<deltaT;
        
        //<<endl;
      //DEBUG<<"pIn : "<<pIn.front();//<<endl;

      //pOut.push_back(pIn.front());

      Parton pp=pIn[0];

      double x[4]={time,0,0,0};
      //Parton pNew(p.plabel(),p.pid(),p.pstat(),p.pt(),p.eta(),p.phi(),p.e(),x);
      double z=0.8;
      //double e=sqrt(p.pt()*z*p.pt()*z*cosh(p.eta())*cosh(p.eta())+p.m()*p.m());
      //double m=sqrt(p.e()*z*z*p.e()-p.pt()*z*p.pt()*z*cosh(p.eta())*cosh(p.eta()));
      //double pp=sqrt(p.e()*z*z*p.e()-p.m()*p.m());

      JSINFO << BOLDYELLOW << pp.plabel() << " pid = " << pp.pid() << " px = " << pp.px() << " py = " << pp.py() << " pz = " << pp.pz()  ;
      //cout<<p.m()<<" "<<m<<endl;
        
        JSINFO << BOLDRED << " x = " << pp.x_in().x() << " y = " << pp.x_in().y() << " z = " << pp.x_in().z() << " t = " << pp.x_in().t() ;

        double blurb;
        
        std::cin >> blurb;
        
      auto vp=pp.GetPseudoJet();
      fjcore::PseudoJet ve(0,0,0,0);
      ve.reset_PtYPhiM(pp.pt()*(1-z),pp.eta(),pp.phi(),0);

      auto p=vp-ve;

      Parton pNew(pp.plabel()+1,pp.pid(),pp.pstat(),p.pt(),p.eta(),p.phi(),p.e(),x);
      Parton pNew2(pp.plabel()+2,pp.pid(),pp.pstat(),ve.pt(),ve.eta(),ve.phi(),ve.e(),x);

      //cout<<pNew.pt()<<endl;
      //Parton (int label, int id, int stat, double pt, double eta, double phi, double e, double* x=0);
      //pOut.push_back(*std::make_shared<Parton>(p.plabel(),p.pid(),p.pstat(),p.pt(),p.eta(),p.phi(),p.e(),p.x_in()));

      pOut.push_back(pNew);
      pOut.push_back(pNew2);

      /*
      Parton pInNew(0,0,0,1,0,0,1,x);
      // Dummy test ... (works)
      pIn.push_back(pInNew);
      */
    }
}
