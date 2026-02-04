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

// Create a pythia collision at a specified point and return the two inital hard partons
#include "PythiaIsrGun.h"
#include <sstream>
#include <iostream>
#include <fstream>
#define MAGENTA "\033[35m"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<PythiaIsrGun> PythiaIsrGun::reg("PythiaIsrGun");

PythiaIsrGun::~PythiaIsrGun() { VERBOSE(8); }

void PythiaIsrGun::InitTask() {


  JSDEBUG << "Initialize PythiaIsrGun";
  VERBOSE(8);

  // Show initialization at INFO level
  readString("Init:showProcesses = off");
  readString("Init:showChangedSettings = off");
  readString("Init:showMultipartonInteractions = on");
  readString("Init:showChangedParticleData = off");
  if (JetScapeLogger::Instance()->GetInfo()) {
    readString("Init:showProcesses = on");
    readString("Init:showChangedSettings = on");
    readString("Init:showMultipartonInteractions = on");
    readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  readString("Next:numberShowInfo = 0");
  readString("Next:numberShowProcess = 0");
  readString("Next:numberShowEvent = 0");

  // Standard settings 
  readString("HardQCD:all = on"); // will repeat this line in the xml for demonstration
  // readString("HardQCD:gg2gg = on");
  // readString("HardQCD:gg2qqbar = on");
  // readString("HardQCD:qg2qg = on");
  // readString("HardQCD:qq2qq = on");
  // readString("HardQCD:qqbar2gg = on");
  // readString("HardQCD:qqbar2qqbarNew = on");
  readString("HardQCD:nQuarkNew = 3"); // Number Of Quark flavours
  readString("MultipartonInteractions:processLevel = 0"); 
  readString("MultipartonInteractions:nQuarkIn = 3"); // Number Of Quark flavours
 
  // readString("HardQCD:gg2ccbar = off");
  // readString("HardQCD:qqbar2ccbar = off");
  // readString("HardQCD:hardccbar = off");
  // readString("HardQCD:gg2bbbar = off");
  // readString("HardQCD:qqbar2bbbar = off");

  //  readString("HardQCD:gg2ccbar = on"); // switch on heavy quark channel
  //readString("HardQCD:qqbar2ccbar = on");
  readString("HadronLevel:Decay = off");
  readString("HadronLevel:all = on");
  readString("PartonLevel:ISR = off");
  readString("PartonLevel:MPI = on");
  //readString("PartonLevel:FSR = on");
  readString("PromptPhoton:all=on");
  readString("WeakSingleBoson:all=off");
  readString("WeakDoubleBoson:all=off");

  // For parsing text
  stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);
  numbf.setf(ios::showpoint);
  numbf.precision(1);
  stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

  std::string s = GetXMLElementText({"Hard", "PythiaGun", "name"});
  SetId(s);
  // cout << s << endl;

  // SC: read flag for FSR
  FSR_on = GetXMLElementInt({"Hard", "PythiaGun", "FSR_on"});
  if (FSR_on)
    readString("PartonLevel:FSR = on");
  else
    readString("PartonLevel:FSR = off");

  pTHatMin = GetXMLElementDouble({"Hard", "PythiaGun", "pTHatMin"});
  pTHatMax = GetXMLElementDouble({"Hard", "PythiaGun", "pTHatMax"});


  JSINFO << MAGENTA << "Pythia Gun with FSR_on: " << FSR_on;
  JSINFO << MAGENTA << "Pythia Gun with " << pTHatMin << " < pTHat < "
         << pTHatMax;

  numbf.str("PhaseSpace:pTHatMin = ");
  numbf << pTHatMin;
  readString(numbf.str());
  numbf.str("PhaseSpace:pTHatMax = ");
  numbf << pTHatMax;
  readString(numbf.str());

  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
  readString("Random:setSeed = on");
  numbi.str("Random:seed = ");
  unsigned int seed = 0;
  if (RandomXmlDescription) {
    tinyxml2::XMLElement *xmle =
        RandomXmlDescription->FirstChildElement("seed");
    if (!xmle)
      throw std::runtime_error("Cannot parse xml");
    xmle->QueryUnsignedText(&seed);
  } else {
    JSWARN << "No <Random> element found in xml, seeding to 0";
  }
  VERBOSE(7) << "Seeding pythia to " << seed;
  numbi << seed;
  readString(numbi.str());

  // Species
  readString("Beams:idA = 2212");
  readString("Beams:idB = 2212");

  // Energy
  eCM = GetXMLElementDouble({"Hard", "PythiaGun", "eCM"});
  numbf.str("Beams:eCM = ");
  numbf << eCM;
  readString(numbf.str());

  // Multiple nuceleon scattering and cross-section
  if (GetXMLElementInt({"Hard", "PythiaGun", "multi_scatter"}) == 1){
    multi_scatter = true;
    JSINFO << MAGENTA << "Initializing PythiaIsrGun with multiple nucleon scatterings enabled.";
  } else
  {
    multi_scatter = false;
  }
  cross_section = GetXMLElementDouble({"cross_section"});
  if ((cross_section < 0) && (multi_scatter == true)) {
    multi_scatter = false;
    JSWARN << "Cross section not set; cannot compute probability for multiple scatterings. Multiple nucleon scattering disabled.";
  }
  proj_A = GetXMLElementInt({"proj_A"});
  targ_A = GetXMLElementInt({"targ_A"});
  if (proj_A<0 || targ_A<0){
    JSWARN << "Projectile or target nuncleon numbers not set. Setting to d-Au values.";
    proj_A=2;
    targ_A=197;
  }

  std::stringstream lines;
  lines << GetXMLElementText({"Hard", "PythiaGun", "LinesToRead"}, false);
  int i = 0;
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue; // skip empty lines
    VERBOSE(7) << "Also reading in: " << s;
    readString(s);
  }

  // And initialize
  if (!init()) { // Pythia>8.1
    throw std::runtime_error("Pythia init() failed.");
  }

    std::ofstream sigma_printer;
    sigma_printer.open(printer, std::ios::trunc);

  // Check Pythia settings for multiple nucleon scatter compatibility
  if ((multi_scatter == true) && (!settings.flag("PhaseSpace:Bias2Selection"))) {
    JSWARN << "PhaseSpace:Bias2Selection is required for multiple nucleon scattering. Disabling multiple nucleon scattering.";
    multi_scatter = false;
  }

  //Check Pythia settings for if pTHat min is less than pTHat ref (can cause weird behaviors where almost everything will have multiple scatterings)
  if ((multi_scatter) && (settings.flag("PhaseSpace:Bias2Selection")) && (pTHatMin < settings.parm("PhaseSpace:pTHatRef"))) {
    JSWARN << "pTHatMin < pTHatRef can cause unexpected behavior with multiple nucleon scattering. Please check your settings.";
    throw std::runtime_error("pTHatMin < pTHatRef can cause unexpected behavior with multiple nucleon scattering. Please check your settings.");
  }
}

void PythiaIsrGun::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  JetScapeTask::WriteTasks(w);
}

void PythiaIsrGun::ExecuteTask() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();
  JSWARN << "Current Event #" << GetCurrentEvent() << "; PythiaIsrGun ExecuteTask called.";
  JSWARN << "The pTHat vector size is " << ini->pTHat.size() << "before clearing. Now clearing";
  ini->pTHat.clear();
  JSWARN << "Have cleared the pTHat vector. The size is now " << ini->pTHat.size();
  //Reading vir_factor from xml for MATTER
  double vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});

  vector<Pythia8::Particle> p62;

  // Initialze loop variables
  int NPP = 0;
  std::vector<int> IndexToSkip;
  std::vector<double> dummy_pTHat;

  //Initialize binary collision points
  std::vector<double> all_t;
  std::vector<double> all_x;
  std::vector<double> all_y;
  std::vector<double> all_z;
  ini->GetAllBinaryCollisionPoints(all_t, all_x, all_y, all_z);
  std::vector<std::vector<double>> all_projPos;
  ini->GetAllBinaryCollisionProjPos(all_projPos);
  std::vector<std::vector<double>> all_targPos;
  ini->GetAllBinaryCollisionTargPos(all_targPos);
  std::vector<int> AcceptedCollisionPoints; //Indices of accepted collision points
  int Ncoll = ini->GetNcoll();

  //Debug
  std::ofstream debug_file;
  debug_file.open("PIG_debug.txt", std::ios::out | std::ios::app);
  debug_file << "Event: " << GetCurrentEvent() << "\n";
  debug_file << "Ncoll: " << Ncoll << "\n";
  debug_file << "Index; \t Binary Collision PT; \t proj pos; \t targ pos \n"; 
  for (int i = 0; i < Ncoll; i++){
    debug_file << i << "; \t (" << all_t[i] << ", " 
    << all_x[i] << ", "
    << all_y[i] << ", "
    << all_z[i] << "); \t (";
    for (const auto& val : all_projPos[i]){
      debug_file << val << " ";
    }
    debug_file << "); \t (";
    for (const auto& val : all_targPos[i]){
      debug_file << val << " ";
    }
    debug_file << ")\n";
  } 

  std::ofstream debug_partons_file;
  debug_partons_file.open("PIG_partons.txt", std::ios::out | std::ios::app);
  debug_partons_file << "\n\nIn Event number " << GetCurrentEvent() << "\n";
  // Debug
  // JSWARN << "Size of all binary collision points: " << all_t.size();
  // JSWARN << "Number of binary collisions from MCGlauber: " << Ncoll;
  // JSWARN << "Size of all_projPos: " << all_projPos.size();
  // JSWARN << "Size of all_targPos: " << all_targPos.size();

  // sort by pt
  struct greater_than_pt {
    inline bool operator()(const Pythia8::Particle &p1,
                           const Pythia8::Particle &p2) {
      return (p1.pT() > p2.pT());
    }
  };

  // determine if two nucleons are at same location
  struct same_location {
    inline bool operator()(const std::vector<double> &p1,
                           const std::vector<double> &p2) {
      if (std::abs(p1[0]-p2[0]) < 1e-8
          && std::abs(p1[1]-p2[1]) < 1e-8
          && std::abs(p1[2]-p2[2]) < 1e-8 
          && std::abs(p1[3]-p2[3]) < 1e-8)
        return true;
      else
        return false;
    }
  };

    FourVector p_p;

  //Variables for checking duplicate collision points
  same_location same_location;
  std::vector<int> index_list; //Shuffled index list to pick collision points and avoid duplicates
  for (int idx =0; idx < Ncoll; idx++){
    index_list.push_back(idx);
  }
  std::shuffle(index_list.begin(), index_list.end(), *GetMt19937Generator());


  //Debug file to verify index matching
  std::ofstream index_match_file;
  index_match_file.open("PIG_index_matches.txt", std::ios::out | std::ios::app);
  index_match_file << "Event: " << GetCurrentEvent() << "\n";
  index_match_file << "Original Index \t Shuffled Index \n";
  for (int idx = 0; idx < index_list.size(); idx++){
    index_match_file << idx << " \t " << index_list[idx] << "\n";
  }
  index_match_file << "\n\n";
  index_match_file.close();

  
  // Loop over possible scatterings to select binary collision point. 
  // Max number of scatterings is min(proj_A, targ_A)
  FourVector x_p;
  for (int iscatt = 0; iscatt < std::min(proj_A, targ_A); iscatt++){
    int NSamplings = 0;
    p62.clear();
    bool flag62 = false; // reset for each scattering so interior loop runs
    bool accept_scatter = true;
    double ratio = 1.0;
    //Debug
    JSINFO << MAGENTA << "Entering new scattering loop: iscatt = " << iscatt << ". Ratio from last call= " << ratio;

    //Pick a collision point
    int icoll = -1;
    if (AcceptedCollisionPoints.size() == 0){
      icoll = index_list[iscatt];
      //Debug
      JSINFO << MAGENTA << "Selected collision point index " << icoll << " for scattering number " << iscatt;
    }
    else {
      //Debug
      JSWARN << "Entered else statement line 319";
      for (int accepted_idx : AcceptedCollisionPoints){
        //Debug
        JSWARN << "At top of for loop over accepted indices";
        JSINFO << MAGENTA << "Accepted index: " << accepted_idx;
        for (int shuffled_idx : index_list){
          bool same_proj, same_targ = true;
          //Debug
          JSWARN << "At top of for loop over shuffled indices";
          JSINFO<< MAGENTA << "Shuffled index: " << shuffled_idx;
          if (accepted_idx == shuffled_idx){
            //Debug
            JSINFO << "Skipped already accepted index: " << accepted_idx;
            continue; //skip already accepted points
          }
          same_proj = same_location(all_projPos[accepted_idx], all_projPos[shuffled_idx]);
          same_targ = same_location(all_targPos[accepted_idx], all_targPos[shuffled_idx]);
          if (!same_proj && !same_targ){
            //Debug
            JSINFO << "Found binary collision with no participant overlap from previous ones. Setting icoll";
            icoll = shuffled_idx;
            JSINFO << "icoll set to " << icoll << " for scatter number " << iscatt;
            break;
          }
        }
        //Debug
        JSWARN << "Out of for loop over shuffled indices";
        if (icoll != -1){//break if already found collision
          //Debug
          JSINFO << MAGENTA << "Already found a collision point, breaking loop";
          break;
        }
      }
      //Debug
      JSWARN << "Out of for loop over shuffled indices";
    }
    //Check icoll and pushback the accepted point
    if (icoll == -1){
      JSWARN << "All binary collision points have been used. Cannot add another scattering between nucleons";
      break; //exit the scattering loop
    }
    else { //Track accepted point and set binary collision position
      AcceptedCollisionPoints.push_back(icoll);
      x_p.Set(all_x[icoll], all_y[icoll], all_z[icoll], all_t[icoll]);
    }
    // Debug
    debug_file << "Selected collision point index " << icoll << " for scattering number " << iscatt << "\n";
    debug_file << "Length of AcceptedCollisionPoints: " << AcceptedCollisionPoints.size() << "\n";
    debug_file << "Collision Point Position (t,x,y,z): (" << x_p.t() << ", " << x_p.x() << ", " << x_p.y() << ", " << x_p.z() << ")\n";

    ReDoSampling:
    do { // loop over samplings in each scattering
      NSamplings++;
      p62.clear();
      IndexToSkip.clear();
      next();

      //Select indices to skip
      for (int parid = 0; parid < event.size(); parid++) {
        if (parid < 3)
          continue; // 0, 1, 2: total event and beams
        Pythia8::Particle &particle = event[parid];

        if (!(particle.isGluon() ||
              particle.isQuark())) { // Getting rid of diquark
          if (particle.status() == -31) {
            IndexToSkip.push_back(particle.daughter1());
            IndexToSkip.push_back(particle.daughter2());
            IndexToSkip.push_back(event[particle.daughter1()].mother1());
            IndexToSkip.push_back(event[particle.daughter1()].mother2());
          } else if (particle.status() == -33) {
            IndexToSkip.push_back(particle.mother1());
            IndexToSkip.push_back(particle.mother2());
            IndexToSkip.push_back(event[particle.mother1()].daughter1());
            IndexToSkip.push_back(event[particle.mother1()].daughter2());
          }
        }
      }

      //Warn for skipped indices
      for (auto &ToSkip : IndexToSkip) {
        JSWARN << " Skipping non-parton index " << ToSkip;
      }
  
      //Accept particles based on status and type
      for (int parid = 0; parid < event.size(); parid++) {
        if (parid < 3)
          continue; // 0, 1, 2: total event and beams
          
        Pythia8::Particle &particle = event[parid];
        if (!FSR_on) {
            
            if ( !( (particle.status() == -21) || (particle.status() == -23) || (particle.status() == -31) || (particle.status() == -33) )) continue ;
            for(auto &ToSkip: IndexToSkip) {
              if(ToSkip == parid){
                goto SkipParton;
              }
            }
            // if ( (particle.status() > -21)||(particle.status()<-23) ) continue ;

          // only accept particles after MPI
          //if (particle.status() != 62)
            //continue;
          // only accept gluons and quarks
          // Also accept Gammas to put into the hadron's list
          //if (fabs(particle.id()) > 5 &&
            //  (particle.id() != 21 && particle.id() != 22))
            //continue;

          // reject rare cases of very soft particles that don't have enough e to get
          // reasonable virtuality
          //if (particle.pT() < 1.0 / sqrt(vir_factor))
            //continue;

          //if(particle.id()==22) cout<<"########this is a photon!######" <<endl;
          // accept
        } else { // FSR_on true: use Pythia vacuum shower instead of MATTER
          if (!particle.isFinal())
            continue;
          // only accept gluons and quarks
          // Also accept Gammas to put into the hadron's list
          if (fabs(particle.id()) > 5 &&
              (particle.id() != 21 && particle.id() != 22))
            continue;
        }
        
          VERBOSE(1) << MAGENTA << " particle from pythiagun id = " << particle.id() << " pz = " << particle.pz() << " px = " << particle.px() << " py = " << particle.py() << " E =  " << particle.e() <<  " status = " << particle.status() << " idex "<< particle.index() << 
          " Color " << particle.col() << " " << particle.acol() << " Mothers " << particle.mother1() << " " << particle.mother2() << " daughter " << particle.daughter1() << " " << particle.daughter2();

          p62.push_back(particle);

          SkipParton:;
      }

      // if you want at least 2
      if (p62.size() < 2)
        continue;

      flag62 = true;

    } while (!flag62);

    //Decide whether to push the information to the framework by putting break here.
    double r = ZeroOneDistribution(*GetMt19937Generator());
    ratio = (GetEventWeight() * GetSigmaGen()) / cross_section;
    // if ( !multi_scatter || r > ratio ) accept_scatter = false;
    accept_scatter = true;

    //Debug
    JSINFO << MAGENTA << "At decision point of scattering loop for scatter # " << iscatt;
    JSINFO << MAGENTA << "weight = " << GetEventWeight() << " sigmaGen = " << GetSigmaGen() << " cross_section = " << cross_section;
    JSINFO << MAGENTA << "ratio = " << ratio << " r = " << r << " scatter_again = " << accept_scatter;
    JSINFO << MAGENTA << "Total number of partons so far NPP = " << NPP << ". p62 size = " << p62.size();
    debug_file << "Accept scatter number " << iscatt << " (0=T, 1=F): " << accept_scatter << "\n";

    if ((!accept_scatter) && (iscatt != 0)) {break;}
    debug_file << "Scatter progressed in scatter loop. Will now push information to framework. \n";

    //If scatter is accepted give to framework
    if (!ini)
    {
      JSINFO << BOLDYELLOW << "No initial state module, setting the starting location to "
                    "0. Make sure to add e.g. trento before PythiaIsrGun.";
    }
    else
    {
      bool pass = false;
      //x_p is already set from AcceptedCollisionPoints Loop
        ini->OutputHardCollisionPosition(x_p.t(), x_p.x(), x_p.y(), x_p.z());
    }

    // only update sigma printer for an event on first scatter (iscatt==0)
      if (!printer.empty() && iscatt==0){
            std::ofstream sigma_printer;
            sigma_printer.open(printer, std::ios::out | std::ios::app);

            sigma_printer << "sigma = " << GetSigmaGen() << " Err =  " << GetSigmaErr() << endl ;
            //sigma_printer.close();
            //JSINFO << BOLDYELLOW << " sigma = " << GetSigmaGen() << " sigma err = " << GetSigmaErr() << " printer = " << printer << " is " << sigma_printer.is_open() ;
      };

    // Loop through particles
    // Accept them all

    double initial_state_label = -1 ;
    double final_state_label = 1 ;
    // ini->pTHat.resize((p62.size())/4);
    dummy_pTHat.resize((p62.size())/4);
    int hCounter = 0;
    SetTotalMomentumPositive(0.0);
    SetTotalMomentumNegative(0.0);
    SetTotalMomentumFractionPositive(0.0);
    SetTotalMomentumFractionNegative(0.0);
    double TotalEnergyOfInitialStatePartons = 0.0;

    for (int np = 0; np < p62.size(); ++np) {
      Pythia8::Particle &particle = p62.at(np);

        if (particle.status()==-21 || particle.status()==-31 )
        {
            TotalEnergyOfInitialStatePartons += particle.e();
            if(particle.pz() >= 0.0) {
              SetTotalMomentumPositive(GetTotalMomentumPositive() + particle.e());
              SetTotalMomentumFractionPositive(GetTotalMomentumFractionPositive() + (particle.e() + particle.pz() ) / ( 0.94 * eCM));
              }
            else {
              SetTotalMomentumNegative(GetTotalMomentumNegative() +particle.e());
              SetTotalMomentumFractionNegative(GetTotalMomentumFractionNegative() + (particle.e() - particle.pz() ) / ( 0.94 * eCM));
              }
        }
    }

    VERBOSE(2) << "Negative Partons Momentum "<< GetTotalMomentumNegative()
            << " Positive Partons Momentum "<< GetTotalMomentumPositive()
            << " eCM " << eCM
            << " TotalEnergyOfInitialStatePartons = " << TotalEnergyOfInitialStatePartons;
    if(GetTotalMomentumFractionNegative() >= 1.0 || GetTotalMomentumFractionPositive() >= 1.){
      JSINFO << "Redoing Pythia Sampling Since MPI energy is larger than eCM/2.1 ";
      if(NSamplings < 1000){
        goto ReDoSampling;
      }
      JSWARN << "Negative Partons Momentum Fraction "<< GetTotalMomentumFractionNegative()
            << " Positive Partons Momentum Fraction "<< GetTotalMomentumFractionPositive()
            << " eCM " << eCM
            << " TotalEnergyOfInitialStatePartons = " << TotalEnergyOfInitialStatePartons;
      throw std::runtime_error("Pythia Isr Gun outputs more energy in the MPI partons than eCM");
    }

    // Decide which particles go to initial and final state modules
    // Update pTHat for each scattering; create and add partons 
    for (int np = 0; np < p62.size(); ++np) {
      Pythia8::Particle &particle = p62.at(np);

        int label = 0;
        int stat = 0;
        
        if (particle.status()==-21 || particle.status()==-31 )
        {
            label = initial_state_label;
            initial_state_label--;
            stat = -1000; // raw initial state status, must go to an initial state module
            debug_partons_file << "Found initial state particle to go to ISR. Label: " << label 
              << " Energy: " << particle.e() << "\n";  
        }
        if (particle.status()==-23 || particle.status()==-33)
        {
            label = final_state_label;
            final_state_label++;
            stat = 1000; // raw final state status, must go to a final state module with virtuality generation. 
            if( (label-1) % 2 == 0){
              dummy_pTHat[(label-1)/2] = particle.pT();
            }
        }

      FourVector p_p(particle.px(),particle.py(),particle.pz(),particle.e());
      
      auto ptn = make_shared<Parton>(label, particle.id(), stat, p_p, x_p);
      ptn->set_color(particle.col());
      ptn->set_anti_color(particle.acol()); 
      ptn->set_max_color(GetMax_ColorPerShower() * (np + 1));
      AddParton(ptn);
    }

    // Update the pTHat vector in initial state using dummy
    for (auto pT : dummy_pTHat){
      ini->pTHat.push_back(pT);
      debug_file << "pT pushed back for scatter " << iscatt << ": " << pT << "\n";
    }
    dummy_pTHat.clear();

    VERBOSE(4) << "PythiaIsrGun scattering number " << iscatt << "completed.";

    // Update NPP
    NPP += p62.size();
  }

  //Set max color for event and collision momenta vectors
  SetMax_Color(GetMax_ColorPerShower() * NPP); 
  FourVector Zeros(0,0,0,0);
  ini->CollisionNegativeMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->CollisionPositiveMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->CollisionNegativeRotatedMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->CollisionPositiveRotatedMomentum = std::vector<FourVector>(NPP/2,Zeros);
  // ini->ClearHardPartonMomentum();

  //Debug
  debug_file << "\n\n";
  debug_file.close();
  VERBOSE(8) << GetNHardPartons();
  debug_partons_file.close();
}
