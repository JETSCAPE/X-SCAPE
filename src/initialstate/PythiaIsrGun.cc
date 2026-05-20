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

  s = GetXMLElementText({"Hard", "PythiaGun", "name"});
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
    JSWARN << "Projectile or target nuncleon numbers not set. Setting to pp values.";
    proj_A=1;
    targ_A=1;
  }

  pythiaLines << GetXMLElementText({"Hard", "PythiaGun", "LinesToRead"}, false);
  int i = 0;
  while (std::getline(pythiaLines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue; // skip empty lines
    VERBOSE(7) << "Also reading in: " << s;
    readString(s);
  }

  outputFilename = GetXMLElementText({"outputFilename"});

  isFirstEvent = true;
  randState = rndm.getState(); //will need to move or delete probably to execute. Just here for test

  std::ofstream sigma_printer;
  sigma_printer.open(printer, std::ios::trunc);

  //Check for multi_scatter and Bias2Selection
  if ((multi_scatter) && (settings.flag("PhaseSpace:Bias2Selection"))) {
    JSWARN << "Multiple scatterings and Bias2Selection can lead to unintended behavior. Turning off multiple scatterings";
    multi_scatter = false;
  }
}

void PythiaIsrGun::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  JetScapeTask::WriteTasks(w);
}

void PythiaIsrGun::ExecuteTask() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();
  // JSWARN << "Current Event #" << GetCurrentEvent() << "; PythiaIsrGun ExecuteTask called.";
  // JSWARN << "The pTHat vector size is " << ini->pTHat.size() << "before clearing. Now clearing";
  ini->pTHat.clear();
  // JSWARN << "Have cleared the pTHat vector. The size is now " << ini->pTHat.size();
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
  std::vector<int> allProjIDs;
  ini->GetAllProjNucleonIDs(allProjIDs);
  std::vector<int> allTargIDs;
  ini->GetAllTargNucleonIDs(allTargIDs);
  // std::vector<int> targCharges;
  // ini->GetAllTargNucleonCharges(targCharges);
  std::vector<int> projCharges;
  ini->GetAllProjNucleonCharges(projCharges);
  std::vector<int> AcceptedCollisionPoints; //INDICES of accepted collision points
                                            //Used to ensure valid hard-scatt site
  int Ncoll = ini->GetNcoll();

  //Debug
  // std::ofstream debug_file;
  // debug_file.open("PIG_debug.txt", std::ios::out | std::ios::app);
  // debug_file << "Event: " << GetCurrentEvent() << "\n";
  // debug_file << "Ncoll: " << Ncoll << "\n";
  // debug_file << "Index; \t Binary Collision PT; \t proj pos; \t targ pos \n"; 
  // for (int i = 0; i < Ncoll; i++){
  //   debug_file << i << "; \t (" << all_t[i] << ", " 
  //   << all_x[i] << ", "
  //   << all_y[i] << ", "
  //   << all_z[i] << "); \t (";
  //   for (const auto& val : all_projPos[i]){
  //     debug_file << val << " ";
  //   }
  //   debug_file << "); \t (";
  //   for (const auto& val : all_targPos[i]){
  //     debug_file << val << " ";
  //   }
  //   debug_file << ")\n";
  // } 

  // std::ofstream debug_partons_file;
  // debug_partons_file.open("PIG_partons.txt", std::ios::out | std::ios::app);
  // debug_partons_file << "\n\nIn Event number " << GetCurrentEvent() << "\n";
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

    FourVector p_p;

  //Variables for checking duplicate collision points
  std::vector<int> index_list; //Shuffled index list to pick collision points and avoid duplicates
  for (int idx =0; idx < Ncoll; idx++){
    index_list.push_back(idx);
  }
  std::shuffle(index_list.begin(), index_list.end(), *GetMt19937Generator());

  // For outputting positions
  std::vector<int> passed_hard_position_idx;


  //Debug file to verify index matching
  // std::ofstream index_match_file;
  // index_match_file.open("PIG_index_matches.txt", std::ios::out | std::ios::app);
  // index_match_file << "Event: " << GetCurrentEvent() << "\n";
  // index_match_file << "Original Index \t Shuffled Index \n";
  // for (int idx = 0; idx < index_list.size(); idx++){
  //   index_match_file << idx << " \t " << index_list[idx] << "\n";
  // }
  // index_match_file << "\n\n";
  // index_match_file.close();

  // Loop over possible scatterings to select binary collision point. 
  // Max number of scatterings is min(proj_A, targ_A)
  FourVector x_p;
  for (int iscatt = 0; iscatt < std::min(targ_A, proj_A); iscatt++){
    int NSamplings = 0;
    p62.clear();
    bool flag62 = false; // reset for each scattering so interior loop runs
    bool doScatt = true;
    if (iscatt > 0 && !multi_scatter) break; //no multiple scatterings if off
  
    /*---Pick a collision point---*/
    int icoll = -1;
    if (AcceptedCollisionPoints.size() == 0){
      icoll = index_list[iscatt]; //Automatically select first point
    }
    else {
      for (int accepted_idx : AcceptedCollisionPoints){
        for (int shuffled_idx : index_list){
          bool same_proj, same_targ;
          if (accepted_idx == shuffled_idx){
            continue; //skip already accepted points
          }
          same_proj = (allProjIDs[accepted_idx] == allProjIDs[shuffled_idx]);
          same_targ = (allTargIDs[accepted_idx] == allTargIDs[shuffled_idx]);

          if (!same_proj && !same_targ){//Break internal check loop at index acceptable index
            icoll = shuffled_idx;
            break;
          }
        }
        if (icoll != -1){//break if already found collision to set icoll
          break;
        }
      }
    }
    //Check icoll and pushback the accepted point
    if (icoll == -1){//all binary points collisions already used
      break; //exit the scattering loop
    }
    else { //Track accepted point and set binary collision position
      AcceptedCollisionPoints.push_back(icoll);
      x_p.Set(all_x[icoll], all_y[icoll], all_z[icoll], all_t[icoll]); //passed to framework later
    }

    /*Set species for projectile beam*/
    std::string projSpecies;
    try{
      if (projCharges[icoll] == 1){//Set proton beam
        projSpecies = "Beams:idA = 2212";
      }
      else if (projCharges[icoll] == 0){//Set neutron beam
        projSpecies = "Beams:idA = 2112";
      }
      else {
        throw std::invalid_argument("Projectile species must have charge +1 or 0. Please check");
      }
    }
    catch (const std::invalid_argument &e) {
      std::cerr << "Caught exception: " << e.what() << std::endl;
      JSWARN << "Projectile species not set, setting to proton";
      projSpecies = "Beams:idA = 2212";
    }

    /*---Decide on proceeding with sampling and initialize---*/
    if (iscatt == 0) {//First scatter always happens
      DefaultInitializePythia(randState, doScatt, projSpecies);
    }
    else{
      JSINFO << MAGENTA << "Will attempt to initialize for totem scattering";
      TableInitializePythia(randState, doScatt, projSpecies);
    }
    if (!doScatt) {//Do not proceed with generation
      break;
    }

    // Debug
    // debug_file << "Selected collision point index " << icoll << " for scattering number " << iscatt << "\n";
    // debug_file << "Length of AcceptedCollisionPoints: " << AcceptedCollisionPoints.size() << "\n";
    // debug_file << "Collision Point Position (t,x,y,z): (" << x_p.t() << ", " << x_p.x() << ", " << x_p.y() << ", " << x_p.z() << ")\n";

    ReDoSampling:
    do { // loop over samplings in each scattering
      NSamplings++;
      flag62=false;
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

    //Update the state of pythia random number generator
    randState = rndm.getState();

    if (iscatt == 0){ // Set first scatter info for printer and getters
      first_sigmaGen = info.sigmaGen();
      first_sigmaErr = info.sigmaErr();
      first_ptHat = info.pTHat();
      first_weight = info.weight();
    }

    //If scatter was generated give to framework below

    //Resize vectors for hard process info
    ResizeTotalMomentumVectors(iscatt+1);
    ini->pTHat.resize(iscatt+1);

    if (!ini)
    {
      JSINFO << BOLDYELLOW << "No initial state module, setting the starting location to "
                    "0. Make sure to add e.g. trento before PythiaIsrGun.";
    }
    else
    {
      bool pass = false;
  
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
    SetTotalMomentumPositive(0.0, iscatt);
    SetTotalMomentumNegative(0.0, iscatt);
    SetTotalMomentumFractionPositive(0.0, iscatt);
    SetTotalMomentumFractionNegative(0.0, iscatt);
    double TotalEnergyOfInitialStatePartons = 0.0;

    for (int np = 0; np < p62.size(); ++np) {
      Pythia8::Particle &particle = p62.at(np);

        if (particle.status()==-21 || particle.status()==-31 )
        {
            TotalEnergyOfInitialStatePartons += particle.e();
            if(particle.pz() >= 0.0) {
              SetTotalMomentumPositive(GetTotalMomentumPositive(iscatt) + particle.e(), iscatt);
              SetTotalMomentumFractionPositive(GetTotalMomentumFractionPositive(iscatt) + (particle.e() + particle.pz() ) / ( 0.94 * eCM),iscatt);
              }
            else {
              SetTotalMomentumNegative(GetTotalMomentumNegative(iscatt) +particle.e(), iscatt);
              SetTotalMomentumFractionNegative(GetTotalMomentumFractionNegative(iscatt) + (particle.e() - particle.pz() ) / ( 0.94 * eCM), iscatt);
              }
        }
    }

    VERBOSE(2) << "Negative Partons Momentum "<< GetTotalMomentumNegative(iscatt)
            << " Positive Partons Momentum "<< GetTotalMomentumPositive(iscatt)
            << " eCM " << eCM
            << " TotalEnergyOfInitialStatePartons = " << TotalEnergyOfInitialStatePartons;
    if(GetTotalMomentumFractionNegative(iscatt) >= 1.0 || GetTotalMomentumFractionPositive(iscatt) >= 1.){
      JSINFO << "Redoing Pythia Sampling Since MPI energy is larger than eCM/2.1 ";
      if(NSamplings < 1000){
        goto ReDoSampling;
      }
      JSWARN << "Negative Partons Momentum Fraction "<< GetTotalMomentumFractionNegative(iscatt)
            << " Positive Partons Momentum Fraction "<< GetTotalMomentumFractionPositive(iscatt)
            << " eCM " << eCM
            << " TotalEnergyOfInitialStatePartons = " << TotalEnergyOfInitialStatePartons;
      throw std::runtime_error("Pythia Isr Gun outputs more energy in the MPI partons than eCM");
    }

    //If all sampling checks passed push back the collision position
    ini->OutputHardCollisionPosition(x_p.t(), x_p.x(), x_p.y(), x_p.z());
    passed_hard_position_idx.push_back(icoll);

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
            // debug_partons_file << "Found initial state particle to go to ISR. Label: " << label 
              // << " Energy: " << particle.e() << "\n";  
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
      ptn->set_hard_scattering(iscatt);
      AddParton(ptn);
    }

    // Update the pTHat vector in initial state using dummy
    for (auto pT : dummy_pTHat){
      ini->pTHat[iscatt].push_back(pT);
      // debug_file << "pT pushed back for scatter " << iscatt << ": " << pT << "\n";
    }
    dummy_pTHat.clear();

    VERBOSE(4) << "PythiaIsrGun scattering number " << iscatt << "completed.";

    // Update NPP
    NPP += p62.size();
  } // End of scattering loop

  //Set max color for event and collision momenta vectors
  SetMax_Color(GetMax_ColorPerShower() * NPP); 
  FourVector Zeros(0,0,0,0);
  ini->CollisionNegativeMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->CollisionPositiveMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->CollisionNegativeRotatedMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->CollisionPositiveRotatedMomentum = std::vector<FourVector>(NPP/2,Zeros);
  // ini->ClearHardPartonMomentum();

  // File for PIG summary of collision points and Ncoll
  std::ofstream PIG_summary;
  PIG_summary.open(outputFilename + std::string("_PIG_summary.dat"), std::ios::out | std::ios::app);
  PIG_summary << "# Event " << GetCurrentEvent() << ", Ncoll " << Ncoll << ", Nscatt: " 
    << passed_hard_position_idx.size() << ", Location (t,x,y,z,bool:hard_scatter)" << std::endl;
  for (int i=0; i<all_t.size(); i++){
    int hard_scatter = 0;
    for (int accepted_idx : passed_hard_position_idx){
      if (i == accepted_idx) {hard_scatter = 1;}
    }
    PIG_summary << all_t[i] << " " << all_x[i] << " " << all_y[i] << " " << all_z[i] << " " << hard_scatter << std::endl;
  }
  PIG_summary.close();
  //Debug
  // debug_file << "\n\n";
  // debug_file.close();
  VERBOSE(8) << GetNHardPartons();
  // debug_partons_file.close();
}

void PythiaIsrGun::DefaultInitializePythia(Pythia8::RndmState randState, bool &doScatt, std::string projSpecies){//Default initialization
  //For parsing text
  stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);
  numbf.setf(ios::showpoint);
  numbf.precision(1);
  stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

  /*Do all the defaults*/
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
  readString("Next:numberShowInfo = 0");
  readString("Next:numberShowProcess = 0");
  readString("Next:numberShowEvent = 0");
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

  numbf.str("PhaseSpace:pTHatMin = ");
  numbf << pTHatMin;
  readString(numbf.str());
  numbf.str("PhaseSpace:pTHatMax = ");
  numbf << pTHatMax;
  readString(numbf.str());

  /*Read in any additional lines to read*/
  pythiaLines.clear();
  pythiaLines.seekg(0, std::ios::beg); //start from beginning of stream
  while (std::getline(pythiaLines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue; // skip empty lines
    VERBOSE(7) << "Also reading in: " << s;
    readString(s);
  }

  /*Set the random state*/
  rndm.setState(randState);
  //set doScatt=True always

  /*Set the projectile species*/
  readString(projSpecies); //should be full line for Pythia

  // And initialize
  if (!init()) { // Pythia>8.1
    throw std::runtime_error("Pythia init() failed.");
  }
  doScatt=true;
}

void PythiaIsrGun::TableInitializePythia(Pythia8::RndmState randState, bool &doScatt, std::string projSpecies){//Initialize via table
  //Start by rolling random number
  double r = ZeroOneDistribution(*GetMt19937Generator());

  std::vector<tableRow> table = RetrieveTable();
  //Find the probability bin from r and probBin boundaries
  int ibin = -1;
  for (int i = table.size() - 1; i >= 0; i--) {//Go backwards since most will be in last bin
    if (r < table[i].probBin[1] && r >= table[i].probBin[0]) {
      ibin = i;
      break;
    }
  }
  if (ibin == -1) {
    throw std::runtime_error("Failed to identify bin for table initialization");
  }
  
  /*----Do the initialization----*/
  //Check if there should be no hard-scatter (first bin)
  if (ibin == table.size() - 1) {
    doScatt = false;
    return;
  }
  else {
    //For parsing text
    stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
    numbf.setf(ios::fixed, ios::floatfield);
    numbf.setf(ios::showpoint);
    numbf.precision(1);
    stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

    /*Do all the defaults*/
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
    readString("Next:numberShowInfo = 0");
    readString("Next:numberShowProcess = 0");
    readString("Next:numberShowEvent = 0");
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

    /*Read in any additional lines to read*/
    pythiaLines.clear();
    pythiaLines.seekg(0, std::ios::beg); //start from beginning of stream
    while (std::getline(pythiaLines, s, '\n')) {
      if (s.find_first_not_of(" \t\v\f\r") == s.npos)
        continue; // skip empty lines
      VERBOSE(7) << "Also reading in: " << s;
      readString(s);
    }

    /*Set the random state*/
    rndm.setState(randState);

    /*Set the projectile species*/
    readString(projSpecies); //should be full line for Pythia

    /*--Read in from the bin selected--*/
    tableRow binInfo = table[ibin];
    //pTHat bin
    numbf.str("PhaseSpace:pTHatMin = ");
    numbf << binInfo.pTHatBin[0];
    readString(numbf.str());
    numbf.str("PhaseSpace:pTHatMax = ");
    numbf << binInfo.pTHatBin[1];
    readString(numbf.str());
    readString("");
    //Process type
    readString(binInfo.processOn);
    readString("");
    readString(binInfo.processOff);
    readString("");
    
    // And initialize
    if (!init()) { // Pythia>8.1
      throw std::runtime_error("Pythia init() failed.");
    }
    doScatt = true;
    JSINFO << MAGENTA << "Successfully initialized Pythia for secondary scatter";
  }
}

std::vector<PythiaIsrGun::tableRow> PythiaIsrGun::RetrieveTable() {
  //scale to more recent pp cross-section STAR 2020. Table made with 42.07
  double x = 42.07/43.82;
  std::vector<tableRow> table = {
    { { 0.0000000000e+00, x* 1.0715277358e-02 }, { 5.0000000000e+00, 7.0000000000e+00 }, "HardQCD:all = on", "PromptPhoton:all = off" },
    { { x* 1.0715277358e-02, x* 1.2418494729e-02 }, { 7.0000000000e+00, 9.0000000000e+00 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2418494729e-02, x* 1.2810934567e-02 }, { 9.0000000000e+00, 1.1000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2810934567e-02, x* 1.2925635374e-02 }, { 1.1000000000e+01, 1.3000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2925635374e-02, x* 1.2964438558e-02 }, { 1.3000000000e+01, 1.5000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2925635374e-02, x* 1.2940426531e-02 }, { 1.5000000000e+01, 1.7000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2940426531e-02, x* 1.2948141745e-02 }, { 1.7000000000e+01, 2.0000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2948141745e-02, x* 1.2951087570e-02 }, { 2.0000000000e+01, 2.5000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2951087570e-02, x* 1.2951590239e-02 }, { 2.5000000000e+01, 3.0000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2951590239e-02, x* 1.2951688998e-02 }, { 3.0000000000e+01, 3.5000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2951688998e-02, x* 1.2951710180e-02 }, { 3.5000000000e+01, 4.0000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2951710180e-02, x* 1.2951714884e-02 }, { 4.0000000000e+01, 4.5000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2951714884e-02, x* 1.2951715943e-02 }, { 4.5000000000e+01, 5.0000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2951715943e-02, x* 1.2951716179e-02 }, { 5.0000000000e+01, 5.5000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.2951716179e-02, x* 1.2951716230e-02 }, { 5.5000000000e+01, 6.0000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.295171623000000e-02, x* 1.295171624243097e-02 }, { 6.000000000000000e+01, 7.000000000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.295171624243097e-02, x* 1.295171624281274e-02 }, { 7.000000000000000e+01, 8.000000000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.295171624281274e-02, x* 1.295171624281810e-02 }, { 8.000000000000000e+01, 9.000000000000000e+01 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.295171624281810e-02, x* 1.295171624281811e-02 }, { 9.000000000000000e+01, 1.000000000000000e+02 }, "HardQCD:all = on", "PromptPhoton:all = off" }, 
    { { x* 1.295171624281811e-02, x* 1.295356107133219e-02 }, { 5.000000000000000e+00, 7.000000000000000e+00 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295356107133219e-02, x* 1.295395966269295e-02 }, { 7.000000000000000e+00, 9.000000000000000e+00 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295395966269295e-02, x* 1.295407385491257e-02 }, { 9.000000000000000e+00, 1.100000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295407385491257e-02, x* 1.295411295149797e-02 }, { 1.100000000000000e+01, 1.300000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295411295149797e-02, x* 1.295412803167064e-02 }, { 1.300000000000000e+01, 1.500000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295412803167064e-02, x* 1.295413440185603e-02 }, { 1.500000000000000e+01, 1.700000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295413440185603e-02, x* 1.295413807381912e-02 }, { 1.700000000000000e+01, 2.000000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295413807381912e-02, x* 1.295413964598038e-02 }, { 2.000000000000000e+01, 2.500000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295413964598038e-02, x* 1.295413994786286e-02 }, { 2.500000000000000e+01, 3.000000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295413994786286e-02, x* 1.295414001193538e-02 }, { 3.000000000000000e+01, 3.500000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295414001193538e-02, x* 1.295414002615953e-02 }, { 3.500000000000000e+01, 4.000000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295414002615953e-02, x* 1.295414002929272e-02 }, { 4.000000000000000e+01, 4.500000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295414002929272e-02, x* 1.295414002996418e-02 }, { 4.500000000000000e+01, 5.000000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295414002996418e-02, x* 1.295414003010266e-02 }, { 5.000000000000000e+01, 5.500000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295414003010266e-02, x* 1.295414003012981e-02 }, { 5.500000000000000e+01, 6.000000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295414003012981e-02, x* 1.295414003013593e-02 }, { 6.000000000000000e+01, 7.000000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" }, 
    { { x* 1.295414003013593e-02, x* 1.295414003013613e-02 }, { 7.000000000000000e+01, 8.000000000000000e+01 }, "PromptPhoton:all = on", "HardQCD:all = off" },
    { {x* 1.295414003013613e-02, 1.0}, {0, 0}, "", ""} 
  };

  // //Debug
  //   std::vector<tableRow> table = { 
  //     {{0,.5}, {5,10}, "HardQCD:all = on", "PromptPhoton:all = off"},
  //     { {.5,1} , {20,30}, "HardQCD:all = off", "PromptPhoton:all = on"}
  //   };
  return table; 
}