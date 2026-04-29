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

#include "EAGun.h"
#include "Matter.h"
#include <sstream>
#include <iostream>
#include <fstream>
#define MAGENTA "\033[35m"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<EAGun> EAGun::reg("EAGun");

EAGun::~EAGun() { VERBOSE(8); }

void EAGun::InitTask() {
    JSDEBUG << "Initialize EAGun";
    VERBOSE(8);

    // Show initialization at INFO level
    readString("Init:showProcesses = off");
    readString("Init:showChangedSettings = off");
    readString("Init:showMultipartonInteractions = off");
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

    // For parsing text
    stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
    numbf.setf(ios::fixed, ios::floatfield);
    numbf.setf(ios::showpoint);
    numbf.precision(1);
    stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

    std::string s = GetXMLElementText({"Hard", "EAGun", "name"});
    SetId(s);

    // Initialize random number distribution
    ZeroOneDistribution = uniform_real_distribution<double>{0.0, 1.0};

    //for p/n sampling
    targZ = GetXMLElementDouble({"Hard", "EAGun", "targetZ"});
    targA = GetXMLElementDouble({"Hard", "EAGun", "targetA"});

    //intrinsic kt
    double ktsigma = 0.865;
    readString("BeamRemnants:primordialKT = on");
    settings.parm("BeamRemnants:primordialKTsoft", ktsigma);
    settings.parm("BeamRemnants:primordialKThard", ktsigma);
    settings.parm("BeamRemnants:halfMassForKT", 0);
    settings.parm("BeamRemnants:reducedKTatHighY", 0);
    settings.parm("BeamRemnants:primordialKTremnant", ktsigma);

    // initial kinematics
    eElectron = GetXMLElementDouble({"Hard", "EAGun", "electron_energy"});
    eProton = GetXMLElementDouble({"Hard", "EAGun", "proton_energy"});
    use_positron = GetXMLElementInt({"Hard", "EAGun", "use_positron"});
    photoproduction = GetXMLElementInt({"Hard", "EAGun", "photoproduction"});
    breitVir = GetXMLElementInt({"Hard", "EAGun", "breit_vir"});
    Q2pow = GetXMLElementDouble({"Hard", "EAGun", "Q2_pow"});
    Q2factor = GetXMLElementDouble({"Hard", "EAGun", "Q2_factor"});
    initial_virtuality_pT = GetXMLElementInt({"Eloss", "Matter", "initial_virtuality_pT"});

    // DIS parameters
    Q2 = GetXMLElementDouble({"Hard", "EAGun", "Q2"});
    nu = GetXMLElementDouble({"Hard", "EAGun", "nu"});

    // kinematic cuts
    // Q2min = GetXMLElementDouble({"Hard", "EAGun", "Q2min"});
    // Q2max = GetXMLElementDouble({"Hard", "EAGun", "Q2max"});
    Q2min = Q2-0.01;
    Q2max = Q2+0.01;
    numin = nu-0.1;
    numax = nu+0.1;
    W2min = GetXMLElementDouble({"Hard", "EAGun", "W2min"});
    W2max = GetXMLElementDouble({"Hard", "EAGun", "W2max"});
    xmin = GetXMLElementDouble({"Hard", "EAGun", "xmin"});
    xmax = GetXMLElementDouble({"Hard", "EAGun", "xmax"});
    ymin = GetXMLElementDouble({"Hard", "EAGun", "ymin"});
    ymax = GetXMLElementDouble({"Hard", "EAGun", "ymax"});

    // other Pythia settings
    readString("HadronLevel:Decay = off");
    readString("HadronLevel:all = off");
    readString("Print:quiet = on");
  
    // beam setup
    // readString("Beams:frameType = 2");
    // settings.parm("Beams:eA", 0.);
    // settings.parm("Beams:eB", eElectron);
    readString("Beams:frameType = 3");
    settings.parm("Beams:pzA", 0.);
    settings.parm("Beams:pzB", eElectron);
    // readString("Beams:allowMomentumSpread = off");
    // settings.parm("Beams:sigmaPxA", 0.);
    // settings.parm("Beams:sigmaPyA", 0.);
    // settings.parm("Beams:sigmaPzA", 0.);
    // settings.parm("Beams:maxDevA", 0.);

    // BeamA = p/n chosen randomly later

    // BeamB = electron
    if (use_positron) {
        readString("Beams:idB = -11");
        JSINFO << "Running with positron beam.";
    }
    else { readString("Beams:idB = 11"); }

    if (photoproduction) {
        readString("PDF:lepton2gamma = on");
        readString("PhotonParton:all = on");
        readString("Photon:Q2max = 1.0");
        readString("Photon:ProcessType = 0");
        readString("SoftQCD:nonDiffractive = on");
    }
    else {
        // Set up DIS process within some phase space.
        // Neutral current (with gamma/Z interference).
        readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
        // Uncomment to allow charged current.
        // readString("WeakBosonExchange:ff2ff(t:W) = on");

        // Phase-space cut: minimal Q2 of process.
        settings.parm("PhaseSpace:Q2Min", Q2min);
        settings.parm("PhaseSpace:Q2Max", Q2max);

        // Set dipole recoil on. Necessary for DIS + shower.
        readString("SpaceShower:dipoleRecoil = on");

        // Allow emissions up to the kinematical limit,
        // since rate known to match well to matrix elements everywhere.
        readString("SpaceShower:pTmaxMatch = 2");

        // QED radiation off lepton not handled yet by the new procedure.
        readString("TimeShower:QEDshowerByL = off");
        //readString("PartonShowers:model = 1");
        //readString("TimeShower:pTmaxMatch = 1");

        //special PDF
        readString("PDF:lepton = off");
        readString("PDF:useHard = on");
        //readString("PDF:pHardSet = LHAPDF6:PDF4LHC21_40"); //for special PDF setting
    }

    // SC: read flag for FSR
    FSR_on = GetXMLElementInt({"Hard", "EAGun", "FSR_on"});
    if (FSR_on) readString("PartonLevel:FSR = on");
    else readString("PartonLevel:FSR = off");

    JSINFO << MAGENTA << "EA Gun with FSR_on: " << FSR_on;

    // readString("Random:setSeed = on");
    // readString("Random:seed = 0");
    
    // random seed
    // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
    tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
    readString("Random:setSeed = on");
    numbi.str("Random:seed = ");
    unsigned int seed = 0;
    if (RandomXmlDescription) {
        tinyxml2::XMLElement *xmle = RandomXmlDescription->FirstChildElement("seed");
        if (!xmle) throw std::runtime_error("Cannot parse xml");
        xmle->QueryUnsignedText(&seed);
    } 
    else { JSWARN << "No <Random> element found in xml, seeding to 0"; }
    JSINFO << BOLDYELLOW << "Seeding pythia to " << seed;
    numbi << seed;
    readString(numbi.str());

    //Reading vir_factor from xml for MATTER
    vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});
    softMomentumCutoff = GetXMLElementDouble({"Hard", "EAGun", "softMomentumCutoff"});
    initial_virtuality_pT = GetXMLElementInt({"Eloss", "Matter", "initial_virtuality_pT"});
    if(vir_factor < rounding_error) {
        JSWARN << "vir_factor should not be zero or negative";
        exit(1);
    }

    std::stringstream lines;
    lines << GetXMLElementText({"Hard", "EAGun", "LinesToRead"}, false);
    int i = 0;
    while (std::getline(lines, s, '\n')) {
        if (s.find_first_not_of(" \t\v\f\r") == s.npos) continue; // skip empty lines
        VERBOSE(7) << "Also reading in: " << s;
        readString(s);
    }

    std::ofstream sigma_printer;
    sigma_printer.open(printer, std::ios::trunc);

    readString("PartonLevel:ISR = off");
    readString("PartonLevel:FSR = off");

    // And initialize
    // if (!init()) { // Pythia>8.1
    //     throw std::runtime_error("Pythia init() failed.");
    // }
    isFirstEvent = true;
}

void EAGun::ExecuteTask() {
    VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
    VERBOSE(8) << "Current Event #" << GetCurrentEvent();


    if (ZeroOneDistribution(*GetMt19937Generator()) < ((double) targZ)/((double) targA)) { //hit proton
        // cout << "PROTON!!!!" << endl;
        readString("Beams:idA = 2212");
    }
    else {
        // cout << "NEUTRON!!!!" << endl;
        readString("Beams:idA = 2112");
    }

    // settings.listChanged();
    if (!init()) { // Pythia>8.1
        throw std::runtime_error("Pythia init() failed.");
    }

    if (!isFirstEvent) { //load rng state
        rndm.setState(randState); 
    }

    bool flag62 = false;
    vector<Pythia8::Particle> p62;
    vector<Pythia8::Particle> p63;

    // sort by pt
    struct greater_than_pt {
        inline bool operator()(const Pythia8::Particle &p1,
                               const Pythia8::Particle &p2) {
            return (p1.pT() > p2.pT());
        }
    };

    Pythia8::Vec4 pProton, peIn, peOut, pPhoton, pProtonNoz, pStruck;
    Pythia8::RotBstMatrix fixedtargBoost, polarRot;
    // Pythia8::Vec4 tmp1,tmp2,tmp3,tmp4;
    double nu, Q2, W2, x, y;

    do {
        bool check = next();
        if (check==false) continue;

        // getting scattered electron index
        int elecID = 6;
        if (photoproduction) {
            for (int iElec=0; iElec<event.size(); iElec++) {
                if (abs(event[iElec].id()) == 11 and event[iElec].status() == 23) { elecID = iElec; }
            }
        }

        pProton = event[1].p();
        peIn    = event[2].p();
        peOut   = event[6].p();
        pPhoton = peIn - peOut;
        pStruck = event[3].p();

        // nu, Q2, W2, Bjorken x, y
        nu = pPhoton.eInFrame(pProton);
        Q2 = - pPhoton.m2Calc();
        W2 = (pProton + pPhoton).m2Calc();
        x  = Q2 / (2. * pProton * pPhoton);
        y  = (pProton * pPhoton) / (pProton * peIn);

        // kinematic cuts
        if(x < xmin or x > xmax) continue;
        if(y < ymin or y > ymax) continue;
        if(Q2 < Q2min or Q2 > Q2max) continue;
        if(W2 < W2min or W2 > W2max) continue;
        if(nu < numin or nu > numax) continue;

        // cout << endl;
        // cout << "BEFORE ROTATION" << endl;
        // cout << "incoming p: " << pProton[0] << " " << pProton[1] << " " << pProton[2] << " " << pProton[3] << " MAG " << pow( pow(pProton[1],2.) + pow(pProton[2],2.) + pow(pProton[3],2.), 0.5) << endl;
        // cout << "incoming e: " << peIn[0] << " " << peIn[1] << " " << peIn[2] << " " << peIn[3] << " tan " << sqrt(pow(peIn[1],2.)+pow(peIn[2],2.))/peIn[3] << " MAG " << pow( pow(peIn[1],2.) + pow(peIn[2],2.) + pow(peIn[3],2.), 0.5) << endl;
        // cout << "outgoing e: " << peOut[0] << " " << peOut[1] << " " << peOut[2] << " " << peOut[3] << " tan " << sqrt(pow(peIn[1],2.)+pow(peIn[2],2.))/peOut[3] << " MAG " << pow( pow(peOut[1],2.) + pow(peOut[2],2.) + pow(peOut[3],2.), 0.5) << endl;
        // cout << "photon    : " << pPhoton[0] << " " << pPhoton[1] << " " << pPhoton[2] << " " << pPhoton[3] << " MAG " << pow( pow(pPhoton[1],2.) + pow(pPhoton[2],2.) + pow(pPhoton[3],2.), 0.5) << endl;
        // cout << "incoming q: " << pStruck[0] << " " << pStruck[1] << " " << pStruck[2] << " " << pStruck[3] << " MAG " << pow( pow(pStruck[1],2.) + pow(pStruck[2],2.) + pow(pStruck[3],2.), 0.5) << endl;

        // tmp1 = pStruck;

        //boost to proton rest frame
        //...except pythia gives it small nonzero pz that we should keep
        pProtonNoz = event[1].p();
        pProtonNoz.pz(0.);

        fixedtargBoost = Pythia8::toCMframe(pProtonNoz, pPhoton, peIn);
        pProton.rotbst(fixedtargBoost);
        peIn.rotbst(fixedtargBoost);
        peOut.rotbst(fixedtargBoost);
        pPhoton.rotbst(fixedtargBoost);
        pStruck.rotbst(fixedtargBoost);

        // cout << "AFTER ROTATING" << endl;
        // cout << "incoming p: " << pProton[0] << " " << pProton[1] << " " << pProton[2] << " " << pProton[3] << " MAG " << pow( pow(pProton[1],2.) + pow(pProton[2],2.) + pow(pProton[3],2.), 0.5) << endl;
        // cout << "incoming e: " << peIn[0] << " " << peIn[1] << " " << peIn[2] << " " << peIn[3] << " tan " << sqrt(pow(peIn[1],2.)+pow(peIn[2],2.))/peIn[3] << " MAG " << pow( pow(peIn[1],2.) + pow(peIn[2],2.) + pow(peIn[3],2.), 0.5) << endl;
        // cout << "outgoing e: " << peOut[0] << " " << peOut[1] << " " << peOut[2] << " " << peOut[3] << " tan " << sqrt(pow(peIn[1],2.)+pow(peIn[2],2.))/peOut[3] << " MAG " << pow( pow(peOut[1],2.) + pow(peOut[2],2.) + pow(peOut[3],2.), 0.5) << endl;
        // cout << "photon    : " << pPhoton[0] << " " << pPhoton[1] << " " << pPhoton[2] << " " << pPhoton[3] << " MAG " << pow( pow(pPhoton[1],2.) + pow(pPhoton[2],2.) + pow(pPhoton[3],2.), 0.5) << endl;
        // cout << "incoming q: " << pStruck[0] << " " << pStruck[1] << " " << pStruck[2] << " " << pStruck[3] << " MAG " << pow( pow(pStruck[1],2.) + pow(pStruck[2],2.) + pow(pStruck[3],2.), 0.5) << endl;
        // cout << endl;
        
        // tmp2 = pStruck;

        // polarRot.reset();
        polarRot.rot(0., 2.*M_PI*ZeroOneDistribution(*GetMt19937Generator()));
        pProton.rotbst(polarRot);
        peIn.rotbst(polarRot);
        peOut.rotbst(polarRot);
        pPhoton.rotbst(polarRot);
        pStruck.rotbst(polarRot);

        // cout << "AFTER ROTATING" << endl;
        // cout << "incoming p: " << pProton[0] << " " << pProton[1] << " " << pProton[2] << " " << pProton[3] << " MAG " << pow( pow(pProton[1],2.) + pow(pProton[2],2.) + pow(pProton[3],2.), 0.5) << endl;
        // cout << "incoming e: " << peIn[0] << " " << peIn[1] << " " << peIn[2] << " " << peIn[3] << " tan " << sqrt(pow(peIn[1],2.)+pow(peIn[2],2.))/peIn[3] << " MAG " << pow( pow(peIn[1],2.) + pow(peIn[2],2.) + pow(peIn[3],2.), 0.5) << endl;
        // cout << "outgoing e: " << peOut[0] << " " << peOut[1] << " " << peOut[2] << " " << peOut[3] << " tan " << sqrt(pow(peIn[1],2.)+pow(peIn[2],2.))/peOut[3] << " MAG " << pow( pow(peOut[1],2.) + pow(peOut[2],2.) + pow(peOut[3],2.), 0.5) << endl;
        // cout << "photon    : " << pPhoton[0] << " " << pPhoton[1] << " " << pPhoton[2] << " " << pPhoton[3] << " MAG " << pow( pow(pPhoton[1],2.) + pow(pPhoton[2],2.) + pow(pPhoton[3],2.), 0.5) << endl;
        // cout << "incoming q: " << pStruck[0] << " " << pStruck[1] << " " << pStruck[2] << " " << pStruck[3] << " MAG " << pow( pow(pStruck[1],2.) + pow(pStruck[2],2.) + pow(pStruck[3],2.), 0.5) << endl;
        // cout << endl;

        // tmp3 = pStruck;

        p62.clear();
        p63.clear();

        if (!printer.empty()) {
            std::ofstream sigma_printer;
            sigma_printer.open(printer, std::ios::out | std::ios::app);

            sigma_printer << "sigma = " << GetSigmaGen() << " Err =  " << GetSigmaErr() << endl ;
            //sigma_printer.close();

            // JSINFO << BOLDYELLOW << " sigma = " << GetSigmaGen() << " sigma err = " << GetSigmaErr() << " printer = " << printer << " is " << sigma_printer.is_open() ;
        }

        // pTarr[0]=0.0; pTarr[1]=0.0;
        // pindexarr[0]=0; pindexarr[1]=0;

        //get the struck parton
        for (int parid = 0; parid < event.size(); parid++) {
            if (parid < 3) continue; // 0, 1, 2: total event and beams
            Pythia8::Particle &particle = event[parid];

            // if (particle.status() == 62) { cout << "62 is " << particle.id() << endl; }

            //skipping everything decayed
            if (!particle.isFinal()) continue;

            //catching scattered electron and beam remenants
            if (particle.isHadron() or particle.isLepton()) {
                AddHadron(EAGun::PythiaToJSHadron(particle));
                continue;
            }

            if (!FSR_on) {
                // only accept gluons and quarks
                // Also accept Gammas to put into the hadron's list
                if (fabs(particle.id()) > 5 && (particle.id() != 21 && particle.id() != 22)) continue;

                // reject rare cases of very soft particles that don't have enough e to get reasonable virtuality
                if (initial_virtuality_pT && (particle.pT() < softMomentumCutoff)) {
                    // this cutoff was 1.0/sqrt(vir_factor) in versions < 3.6
                    continue;
                } 
                else if (!initial_virtuality_pT && (particle.pAbs() < softMomentumCutoff)) continue;
            }
            else { // FSR_on true: use Pythia vacuum shower instead of MATTER
                // only accept gluons and quarks
                // Also accept Gammas to put into the hadron's list
                if (fabs(particle.id()) > 5 && (particle.id() != 21 && particle.id() != 22)) continue;
            }

            p62.push_back(particle);
        }

        // if you want at least 2
        // if (p62.size() < 2) continue;
        if (p62.size() < 1) continue;

        //now that we have a good struck parton
        //take all of the associated remnants
        for (int parid = 0; parid < event.size(); parid++) {
            if (parid < 3) continue; // 0, 1, 2: total event and beams
            Pythia8::Particle &particle = event[parid];

            if (particle.status() != 63) continue;

            p63.push_back(particle);
        }


        // Now have all candidates, sort them by pt
        // std::sort(p62.begin(), p62.end(), greater_than_pt());
        // check...
        // for (auto& p : p62 ) cout << p.pT() << endl;

        // event.list();
        // for (auto& p : p62 ) cout << "62 " << p.pz() << endl;
        // for (auto& p : p63 ) cout << "63 " << p.pz() << endl;

        flag62 = true;

        // int found61already = 0;
        // int numremn = 0;
        // int remni = 0;

        // for (int parid = 0; parid < event.size(); parid++) {
        //     Pythia8::Particle &particle = event[parid];
        //     if (particle.status() == -61) {
        //         if (found61already == 1) { cout << "ALREADY DID A -61 HERE" << endl; exit(-2); }
        //         Pythia8::Vec4 pmed = particle.p();

        //         std::ofstream foutx1;
        //         foutx1.open("fullrot-qmed1-v4-19.txt", std::ios_base::app);
        //         foutx1 << pmed[1] << " " << pmed[2] << " " << pmed[3] << endl;

        //         pmed.rotbst(fixedtargBoost);

        //         std::ofstream foutx2;
        //         foutx2.open("fullrot-qmed2-v4-19.txt", std::ios_base::app);
        //         foutx2 << pmed[1] << " " << pmed[2] << " " << pmed[3] << endl;

        //         pmed.rotbst(polarRot);

        //         std::ofstream foutx3;
        //         foutx3.open("fullrot-qmed3-v4-19.txt", std::ios_base::app);
        //         foutx3 << pmed[1] << " " << pmed[2] << " " << pmed[3] << endl;

        //         found61already = 1;
        //     }
        //     if (particle.status() == 63) { numremn++; remni=parid; }
        // }
        // if (numremn==1) {
        //     Pythia8::Particle &particle = event[remni];
        //     Pythia8::Vec4 pmed = particle.p();

        //     std::ofstream foutx1;
        //     foutx1.open("fullrot-remn1-v4-19.txt", std::ios_base::app);
        //     foutx1 << pmed[1] << " " << pmed[2] << " " << pmed[3] << endl;

        //     pmed.rotbst(fixedtargBoost);

        //     std::ofstream foutx2;
        //     foutx2.open("fullrot-remn2-v4-19.txt", std::ios_base::app);
        //     foutx2 << pmed[1] << " " << pmed[2] << " " << pmed[3] << endl;

        //     pmed.rotbst(polarRot);

        //     std::ofstream foutx3;
        //     foutx3.open("fullrot-remn3-v4-19.txt", std::ios_base::app);
        //     foutx3 << pmed[1] << " " << pmed[2] << " " << pmed[3] << endl;
        // }
    } while (!flag62);

    // std::ofstream fout1;
    // fout1.open("fullrot-qin1-v4-19.txt", std::ios_base::app);
    // fout1 << tmp1[1] << " " << tmp1[2] << " " << tmp1[3] << endl;

    // std::ofstream fout2;
    // fout2.open("fullrot-qin2-v4-19.txt", std::ios_base::app);
    // fout2 << tmp2[1] << " " << tmp2[2] << " " << tmp2[3] << endl;

    // std::ofstream fout3;
    // fout3.open("fullrot-qin3-v4-19.txt", std::ios_base::app);
    // fout3 << tmp3[1] << " " << tmp3[2] << " " << tmp3[3] << endl;

    // event passing kinematical cuts has been generated
    // identify where it occurred
    double xLoc[4];
    if (!ini) {
        VERBOSE(1) << "No initial state module, setting the starting location to "
                      "0. Make sure to add e.g. 3DGlauber before EAGun.";
        for (int i=0; i<4; i++) { xLoc[i] = 0.; }
    } 
    else {
        nucleonPositions = ini->GetTargetNucleonPositions();
        int randInd = (int)( ZeroOneDistribution(*GetMt19937Generator()) * nucleonPositions.size() );
        for (int i=0; i<4; i++) { xLoc[i] = nucleonPositions[randInd][i]; }
    }
    // cout << "INITIAL COLLISION IS AT: " << xLoc[0] << " " << xLoc[1] << " " << xLoc[2] << " " << xLoc[3] << endl;

    //debug: should only get one parton out
    if (p62.size()>1) { 
        cout << "SIZE IS " << p62.size() << endl;
        for (int np=0; np<p62.size(); np++) {
            Pythia8::Particle &particle = p62.at(np);
            std::cout << np << " Particle " << particle.id() << " E " << particle.e() << " MOM " << particle.px() << " " << particle.py() << " " << particle.pz() << endl;
        }
        // exit(-2);
    }

    Pythia8::Particle particle, remnant;
    Pythia8::Vec4 pPtn, pRmn;
    FourVector pParton, pRemnant;

    for (int np=0; np<p62.size(); np++) {
        particle = p62.at(np);

        pPtn = particle.p();
        pPtn.rotbst(fixedtargBoost);

        // std::ofstream fout4;
        // fout4.open("fullrot-qout1-v4-19.txt", std::ios_base::app);
        // fout4 << pPtn.px() << " " << pPtn.py() << " " << pPtn.pz() << endl;

        pPtn.rotbst(polarRot);
        pParton = FourVector(pPtn.px(), pPtn.py(), pPtn.pz(), pPtn.e());

        auto ptn = make_shared<Parton>(0, particle.id(), 0, pParton, xLoc);
        ptn->set_color(particle.col());
        ptn->set_anti_color(particle.acol());
        ptn->set_max_color(1000 * (np + 1));

        AddParton(ptn);

        std::cout << "Struck " << particle.id() << " " << particle.status() << " E " << ptn->e() << " MOM " << ptn->px() << " " << ptn->py() << " " << ptn->pz() << endl;

        //log the momentum
        // std::ofstream ptout;
        // ptout.open("fullrot-qout2-v4-19.txt", std::ios_base::app);
        // ptout << pPtn.px() << " " << pPtn.py() << " " << pPtn.pz() << endl;
    }


    for (int np=0; np<p63.size(); np++) {
        remnant = p63.at(np);

        pRmn = remnant.p();
        pRmn.rotbst(fixedtargBoost);

        // std::ofstream fout4;
        // fout4.open("fullrot-remnout1-v4-19.txt", std::ios_base::app);
        // fout4 << pRmn.px() << " " << pRmn.py() << " " << pRmn.pz() << endl;

        pRmn.rotbst(polarRot);
        pRemnant = FourVector(pRmn.px(), pRmn.py(), pRmn.pz(), pRmn.e());

        if (remnant.isHadron()) {
            auto rmn = make_shared<Hadron>(0, remnant.id(), 0, pRemnant, xLoc);

            AddHadron(rmn);
            std::cout << "Remnant Hadron " << remnant.id() << " " << remnant.status() << " E " << rmn->e() << " MOM " << rmn->px() << " " << rmn->py() << " " << rmn->pz() << endl;
        }
        else {
            auto rmn = make_shared<Parton>(0, remnant.id(), 0, pRemnant, xLoc);
            rmn->set_color(remnant.col());
            rmn->set_anti_color(remnant.acol());
            rmn->set_max_color(1000 * (np + 1));

            AddParton(rmn);
            std::cout << "Remnant Parton " << remnant.id() << " " << remnant.status() << " E " << rmn->e() << " MOM " << rmn->px() << " " << rmn->py() << " " << rmn->pz() << endl;
        }



        //log the momentum
        // std::ofstream ptout;
        // ptout.open("fullrot-remnout2-v4-19.txt", std::ios_base::app);
        // ptout << pRmn.px() << " " << pRmn.py() << " " << pRmn.pz() << endl;
    }

    randState = rndm.getState();
    isFirstEvent = false;

    VERBOSE(8) << GetNHardPartons();
}
