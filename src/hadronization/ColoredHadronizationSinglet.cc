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

#include "ColoredHadronizationSinglet.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include "tinyxml2.h"
#include "TMatrixD.h"
#include "TArrayD.h"

using namespace Jetscape;
using namespace Pythia8;

// Register the module with the base class
RegisterJetScapeModule<ColoredHadronizationSinglet>
    ColoredHadronizationSinglet::reg("ColoredHadronizationSinglet");

Pythia8::Pythia ColoredHadronizationSinglet::pythia("IntentionallyEmpty", false);

ColoredHadronizationSinglet::ColoredHadronizationSinglet() {
  SetId("MyHadroTest");
  VERBOSE(8);
}

ColoredHadronizationSinglet::~ColoredHadronizationSinglet() { VERBOSE(8); }

void ColoredHadronizationSinglet::InitTask() {

  std::string s = GetXMLElementText({"JetHadronization", "name"});
  JSDEBUG << s << " to be initializied ...";

  double p_read_xml =
      GetXMLElementDouble({"JetHadronization", "eCMforHadronization"});
  p_fake = p_read_xml / 6.;

  /*std::string weak_decays =
      GetXMLElementText({"JetHadronization", "weak_decays"});*/

  VERBOSE(2) << "Start Hadronizing using the PYTHIA module...";

  // Show initialization at DEBUG or high verbose level
  pythia.readString("Init:showProcesses = off");
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showMultipartonInteractions = off");
  pythia.readString("Init:showChangedParticleData = off");
  if (JetScapeLogger::Instance()->GetDebug() ||
      JetScapeLogger::Instance()->GetVerboseLevel() > 2) {
    pythia.readString("Init:showProcesses = on");
    pythia.readString("Init:showChangedSettings = on");
    pythia.readString("Init:showMultipartonInteractions = on");
    pythia.readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  if (JetScapeLogger::Instance()->GetDebug() ||
      JetScapeLogger::Instance()->GetVerboseLevel() > 2) {
    pythia.readString("Next:numberShowInfo = 1");
    pythia.readString("Next:numberShowProcess = 1");
    pythia.readString("Next:numberShowEvent = 1");
  }

  pythia.readString("ProcessLevel:all = off");
  pythia.readString("PartonLevel:FSR=off");

  // General settings for hadron decays
  std::string pythia_decays = GetXMLElementText({"JetHadronization", "pythia_decays"});
  double tau0Max = 10.0;
  double tau0Max_xml = GetXMLElementDouble({"JetHadronization", "tau0Max"});
	if(tau0Max_xml >= 0){tau0Max = tau0Max_xml;}
  else{JSWARN << "tau0Max should be larger than 0. Set it to 10.";}
  if(pythia_decays == "on"){
    JSINFO << "Pythia decays are turned on for tau0Max < " << tau0Max;
    pythia.readString("HadronLevel:Decay = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = " + std::to_string(tau0Max));
  } else {
    JSINFO << "Pythia decays are turned off";
    pythia.readString("HadronLevel:Decay = off");
  }

  // Settings for decays (old flag, will be depracted at some point)
  // This overwrites the previous settings if the user xml file contains the flag
  std::string weak_decays =
    GetXMLElementText({"JetHadronization", "weak_decays"});
  if (weak_decays == "off") {
    JSINFO << "Hadron decays are turned off.";
    JSWARN << "This parameter will be depracted at some point. Use 'pythia_decays' instead.\nOverwriting 'pythia_decays'.";
    pythia.readString("HadronLevel:Decay = off");
  } else if(weak_decays == "on") {
    JSINFO << "Hadron decays inside a range of 10 mm/c are turned on.";
    JSWARN << "This parameter will be depracted at some point. Use 'pythia_decays' and 'tau0Max' for more control on decays.\nOverwriting 'pythia_decays' and fix 'tau0Max' to 10.";
    pythia.readString("HadronLevel:Decay = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = 10.0");
  }

  std::stringstream lines;
  lines << GetXMLElementText({"JetHadronization", "LinesToRead"}, false);
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue; // skip empty lines
    JSINFO << "Also reading in: " << s;
    pythia.readString(s);
  }

  targZ = GetXMLElementDouble({"Hard", "EAGun", "targetZ"});
  targA = GetXMLElementDouble({"Hard", "EAGun", "targetA"});

  shiftHadrons = GetXMLElementInt({"JetHadronization", "shift_hadrons"});
  shiftType = GetXMLElementInt({"JetHadronization", "shift_type"});

  if (shiftHadrons) {
    if (!shiftType) {
      JSWARN << "Hadron shifting turned on but shifting algorithm not set, use 'shift_type' to choose. Defaulting to Monte Carlo.";
      shiftType = MC_SHIFT_TYPE;
    }
    if (shiftType>2) {
      JSWARN << "Invalid shift type. Options are 1=Voronoi and 2=Monte Carlo. Defaulting to Monte Carlo.";
      shiftType = MC_SHIFT_TYPE;
    }
  }
  if (!shiftHadrons && shiftType) {
    JSINFO << "Hadron shifting algorithm chosen even though shifting is turned off. Ignoring.";
  }

  pythia.init();

  //initial state pointer setting
  ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
}

void ColoredHadronizationSinglet::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  auto f = w.lock();
  if (!f)
    return;
  f->WriteComment("Hadronization Module : " + GetId());
  f->WriteComment("Hadronization to be implemented accordingly ...");
}


void ColoredHadronizationSinglet::GetCircumSphere(std::vector<double> A, std::vector<double> B, std::vector<double> C, std::vector<double> D, double &cx, double &cy, double &cz, double &cr) {
  double Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz, Dx,Dy,Dz;
  Ax=A[0]; Ay=A[1]; Az=A[2];
  Bx=B[0]; By=B[1]; Bz=B[2];
  Cx=C[0]; Cy=C[1]; Cz=C[2];
  Dx=D[0]; Dy=D[1]; Dz=D[2];
  // cout << Ax << " " << Ay << " " << Az << endl;
  // cout << Bx << " " << By << " " << Bz << endl;
  // cout << Cx << " " << Cy << " " << Cz << endl;
  // cout << Dx << " " << Dy << " " << Dz << endl;

  TMatrixD matA, matB, matX, matY, matZ;
  Double_t detA, detB, detX, detY, detZ;

  double Asq = Ax*Ax + Ay*Ay + Az*Az;
  double Bsq = Bx*Bx + By*By + Bz*Bz;
  double Csq = Cx*Cx + Cy*Cy + Cz*Cz;
  double Dsq = Dx*Dx + Dy*Dy + Dz*Dz;

  double mValA[] = {Ax, Ay, Az, 1,
                    Bx, By, Bz, 1,
                    Cx, Cy, Cz, 1,
                    Dx, Dy, Dz, 1};
  double mValB[] = {Asq, Ax, Ay, Az,
                    Bsq, Bx, By, Bz,
                    Csq, Cx, Cy, Cz,
                    Dsq, Dx, Dy, Dz};
  double mValX[] = {Asq, Ay, Az, 1,
                    Bsq, By, Bz, 1,
                    Csq, Cy, Cz, 1,
                    Dsq, Dy, Dz, 1};
  double mValY[] = {Asq, Ax, Az, 1,
                    Bsq, Bx, Bz, 1,
                    Csq, Cx, Cz, 1,
                    Dsq, Dx, Dz, 1};
  double mValZ[] = {Asq, Ax, Ay, 1,
                    Bsq, Bx, By, 1,
                    Csq, Cx, Cy, 1,
                    Dsq, Dx, Dy, 1};

  TArrayD mDataA(16, mValA);
  TArrayD mDataB(16, mValB);
  TArrayD mDataX(16, mValX);
  TArrayD mDataY(16, mValY);
  TArrayD mDataZ(16, mValZ);

  matA.Use(4,4,mDataA.GetArray());
  matB.Use(4,4,mDataB.GetArray());
  matX.Use(4,4,mDataX.GetArray());
  matY.Use(4,4,mDataY.GetArray());
  matZ.Use(4,4,mDataZ.GetArray());

  // cout << "A" << endl;
  detA = matA.Determinant();
  // cout << "B" << endl;
  detB = matB.Determinant();
  // cout << "X" << endl;
  detX = matX.Determinant();
  // cout << "Y" << endl;
  detY = matY.Determinant();
  // cout << "Z" << endl;
  detZ = matZ.Determinant();

  cx = detX/detA/2.;
  cy = detY/detA/-2.; //minus sign is not a typo
  cz = detZ/detA/2.;

  // cout << Ax << " " << Ay << " " << Az << endl;
  // cout << Bx << " " << By << " " << Bz << endl;
  // cout << Cx << " " << Cy << " " << Cz << endl;
  // cout << Dx << " " << Dy << " " << Dz << endl;
  // cout << "ANSWER " << cx << " " << cy << " " << cz << endl;

  cr = sqrt(cx*cx + cy*cy + cz*cz - detB/detA);

  // cout << "RADIUS " << cr << endl;
}


int ColoredHadronizationSinglet::IsInSphere(double cx, double cy, double cz, double cr, std::vector<double> P) {
    //is P inside of the sphere at (cx,cy,cz) with radius cr

    double dpx = P[0]-cx;
    double dpy = P[1]-cy;
    double dpz = P[2]-cz;

    double cr2 = cr*cr;
    double dp = dpx*dpx + dpy*dpy + dpz*dpz;

    if (dp>cr2) { return 1; }
    if (dp<cr2) { return -1; }
    return 0;
}

int ColoredHadronizationSinglet::IsInCircumsphere(std::vector<double> A, std::vector<double> B, std::vector<double> C, std::vector<double> D, std::vector<double> P) {
    //is P inside of the circumsphere defined by A,B,C,D

    double cx, cy, cz, cr;
    // cout << "\tCALLING CIRCUMSPHERE B" << endl;
    // cout << "\tA " << A[0] << " " << A[1] << " " << A[2] << endl;
    // cout << "\tB " << B[0] << " " << B[1] << " " << B[2] << endl;
    // cout << "\tC " << C[0] << " " << C[1] << " " << C[2] << endl;
    // cout << "\tD " << D[0] << " " << D[1] << " " << D[2] << endl;
    // cout << "\tP " << P[0] << " " << P[1] << " " << P[2] << endl;
    GetCircumSphere(A, B, C, D, cx, cy, cz, cr);
    // cout << "\tcenter is" << endl;
    // cout << "\t\t" << cx << " " << cy << " " << cz << " " << cr << endl;
    return IsInSphere(cx, cy, cz, cr, P);
}



ColoredHadronizationSinglet::Sphere ColoredHadronizationSinglet::InitSphere(double x, double y, double z, int pid) {
  Sphere retSphere;
  retSphere.x = x;
  retSphere.y = y;
  retSphere.z = z;
  if (pid == 2212 || pid == 2112) { retSphere.r = 1.; } //glauber nucleon radius
  else { retSphere.r = pow(2./3., 1./3.); } //mesons have 2/3 the quarks so just say 2/3 of volume for now
  return retSphere;
}

double ColoredHadronizationSinglet::d2Sphere(Sphere A, Sphere B) {
  return pow(A.x - B.x, 2) + pow(A.y - B.y, 2) + pow(A.z - B.z, 2);
}
      


std::vector<std::vector<int>> ColoredHadronizationSinglet::DelaunayTriangulate3D(std::vector<std::vector<double>> coords) {
  //each coords[i] is a x,y,z point
  //each tetrahedron is a list of 4 ints referring to the indices of coords[]
  std::vector<std::vector<int>> tetrahedra; //return val

  int numCoords = coords.size(); //anything past this index is part of the supertriangle

  //figure out supertriangle vertices and append coords
  //TODO: OPTIMIZE SUPERTRIANGLE COMPUTATION
  coords.push_back({50.,0,-25.});
  coords.push_back({50.*std::cos(pi*2./3.), 50.*std::sin(pi*2./3.), -25.1});
  coords.push_back({50.*std::cos(pi*4./3.), 50.*std::sin(pi*4./3.), -25.});
  coords.push_back({0,0,25.});
  tetrahedra.push_back({numCoords, numCoords+1, numCoords+2, numCoords+3});

  //add first point manually
  tetrahedra.push_back({numCoords, numCoords+1, numCoords+2, 0});
  tetrahedra.push_back({numCoords, numCoords+1, numCoords+3, 0});
  tetrahedra.push_back({numCoords, numCoords+2, numCoords+3, 0});
  tetrahedra.push_back({numCoords+1, numCoords+2, numCoords+3, 0});

  int isInside;
  for (int vi=1; vi<numCoords; vi++) { //add points one by one
    cout << "adding coord number " << vi << endl;
    //find all tetrahedra that are no longer valid
    std::vector<int> badTrianglesis;
    for (int trii=0; trii<tetrahedra.size(); trii++) { 
      std::vector<int> tetrai = tetrahedra[trii];
      //TODO: KEEP A LUT OR HASH TABLE OF CIRCUMCIRCLE INFO PER TETRAHEDRON
      // cout << "CALLING ISINCIRCUMSPHERE " << trii << " " << tetrai[0] << " " << tetrai[1] << " " << tetrai[2] << " " << tetrai[3] << endl;
      isInside = IsInCircumsphere(coords[tetrai[0]], coords[tetrai[1]], coords[tetrai[2]], coords[tetrai[3]], coords[vi]);
      // cout << "\t\t" << isInside << endl;
      if (isInside == -1) { badTrianglesis.push_back(trii); } //point is inside
    }
    cout << "bad triangles: " << badTrianglesis.size() << endl;

    // if (vi==3) { exit(1); }

    //find the boundary of the bad tetrahedra to connect the new point to
    std::vector<std::vector<int>> borderPlanes;
    bool b1, b2, b3, b4; //is the plane made by excluding the nth vertex on the boundary
    int cnt1, cnt2, cnt3, cnt4;
    for (int badtrii=0; badtrii<badTrianglesis.size(); badtrii++) { //for every bad triangle
      b1 = true; b2 = true; b3 = true; b4 = true;
      for (int othertrii=0; othertrii<badTrianglesis.size(); othertrii++) { //for every other triangle
        if (othertrii != badtrii) {
          cnt1 = count(tetrahedra[othertrii].begin(), tetrahedra[othertrii].end(), tetrahedra[badTrianglesis[badtrii]][0]);
          cnt2 = count(tetrahedra[othertrii].begin(), tetrahedra[othertrii].end(), tetrahedra[badTrianglesis[badtrii]][1]);
          cnt3 = count(tetrahedra[othertrii].begin(), tetrahedra[othertrii].end(), tetrahedra[badTrianglesis[badtrii]][2]);
          cnt4 = count(tetrahedra[othertrii].begin(), tetrahedra[othertrii].end(), tetrahedra[badTrianglesis[badtrii]][3]);
          if (cnt2>0 && cnt3>0 && cnt4>0) { b1=false; }
          if (cnt1>0 && cnt2>0 && cnt4>0) { b2=false; }
          if (cnt1>0 && cnt3>0 && cnt4>0) { b3=false; }
          if (cnt2>0 && cnt3>0 && cnt4>0) { b4=false; }
        }
      }
      if (b1) { borderPlanes.push_back({ tetrahedra[badTrianglesis[badtrii]][1], tetrahedra[badTrianglesis[badtrii]][2], tetrahedra[badTrianglesis[badtrii]][3] }); }
      if (b2) { borderPlanes.push_back({ tetrahedra[badTrianglesis[badtrii]][0], tetrahedra[badTrianglesis[badtrii]][2], tetrahedra[badTrianglesis[badtrii]][3] }); }
      if (b3) { borderPlanes.push_back({ tetrahedra[badTrianglesis[badtrii]][0], tetrahedra[badTrianglesis[badtrii]][1], tetrahedra[badTrianglesis[badtrii]][3] }); }
      if (b4) { borderPlanes.push_back({ tetrahedra[badTrianglesis[badtrii]][0], tetrahedra[badTrianglesis[badtrii]][1], tetrahedra[badTrianglesis[badtrii]][2] }); }
    }

    // cout << "removing" << endl;

    //remove the bad tetrahedra
    //easiest way to remove: swap with the last, then delete the last
    std::sort(badTrianglesis.begin(), badTrianglesis.end(), std::greater<>());
    for (int badtrii=0; badtrii<badTrianglesis.size(); badtrii++) {
      tetrahedra[badTrianglesis[badtrii]] = tetrahedra.back();
      tetrahedra.pop_back();
    }


    // cout << tetrahedra.size() << " adding new " << tetrahedra.size() + borderPlanes.size() << endl;

    //add the new tetrahedra formed using this point
    for (int planei=0; planei<borderPlanes.size(); planei++) {
      tetrahedra.push_back({borderPlanes[planei][0], borderPlanes[planei][1], borderPlanes[planei][2], vi});
    }
    cout << "NEW TETRA " << borderPlanes.size() << endl;
    
    //add back super
    tetrahedra.push_back({numCoords, numCoords+1, numCoords+2, numCoords+3});
    cout << "TOTAL TETRA " << tetrahedra.size() << endl;
  }

  //now remove anything involving the supertriangle
  for (unsigned tetrai = tetrahedra.size(); tetrai-- > 0; ) {
    for (int vi=0; vi<4; vi++) {
      if (tetrahedra[tetrai][vi] >= numCoords) {
        tetrahedra[tetrai] = tetrahedra.back();
        tetrahedra.pop_back();
        break;
      }
    }
  }

  cout << "FINAL NUM " << tetrahedra.size() << endl;
  return tetrahedra;
}

std::vector<ColoredHadronizationSinglet::Sphere> ColoredHadronizationSinglet::VoronoiHoles(std::vector<ColoredHadronizationSinglet::Sphere> inHadrons) {
  // cout << "Finding voronoi" << endl;

  std::vector<std::vector<double>> centers;
  for (int hi=0; hi<inHadrons.size(); hi++) {
    centers.push_back({inHadrons[hi].x, inHadrons[hi].y, inHadrons[hi].z});
  }

  // cout << "\tmaking Delaunay" << endl;
  std::vector<std::vector<int>> delaunay = DelaunayTriangulate3D(centers);
  // cout << "\tdone" << endl;
  cout << "delaunay size " << delaunay.size() << endl;

  //TODO: KEEP CIRCUMSPHERE DATA WITH TRIANGULATION
  std::vector<Sphere> holes;
  double cx, cy, cz, cr;
  int h1, h2, h3, h4;
  Sphere hole;
  for (int tetrai=0; tetrai<delaunay.size(); tetrai++) {
    h1 = delaunay[tetrai][0];
    h2 = delaunay[tetrai][1];
    h3 = delaunay[tetrai][2];
    h4 = delaunay[tetrai][3];

    // cout << "CALLING CIRCUMSPHERE A" << endl;
    GetCircumSphere(centers[h1], centers[h2], centers[h3], centers[h4], cx, cy, cz, cr);
    
    hole.x = cx;
    hole.y = cy;
    hole.z = cz;
    hole.r = cr - std::max({inHadrons[h1].r, inHadrons[h2].r, inHadrons[h3].r, inHadrons[h4].r});
    
    holes.push_back(hole);
  }

  //TODO: include the smaller midpoint holes and interpolation between
  return holes;
}

bool ColoredHadronizationSinglet::isOverlap(double *x, int pid, std::vector<ColoredHadronizationSinglet::Sphere> nucleonSpheres) {
  //check if we are inside a nucleon
  Sphere thisSphere = InitSphere(x[1], x[2], x[3], pid);

  Sphere nucSphere;
  for (int nuci=0; nuci<nucleonSpheres.size(); nuci++) {
    nucSphere = nucleonSpheres[nuci];
    if (d2Sphere(nucSphere, thisSphere) < pow(nucSphere.r + thisSphere.r, 2)) { return true; }
  }

  return false;
}

void ColoredHadronizationSinglet::randUnitR3(double &x, double &y, double &z) {
    // std::uniform_real_distribution<double> dist(0.0, 1.0);

    double r = sqrt(1.-z*z);
    double phi = 2.*M_PI*ZeroOneDistribution(*GetMt19937Generator());
    z = 2.*ZeroOneDistribution(*GetMt19937Generator()) - 1;

    x = r*cos(phi);
    y = r*sin(phi);
}


void ColoredHadronizationSinglet::DoHadronization(
    vector<vector<shared_ptr<Parton>>> &shower,
    vector<shared_ptr<Hadron>> &hOut, vector<shared_ptr<Parton>> &pOut) {

  cout << "HADRONIZING" << endl;

  // double centx=2.;
  // double centy=3.;
  // double centz=1.;
  // std::vector<double> pA = {1+centx,0+centy,0+centz};
  // std::vector<double> pB = {0+centx,1+centy,0+centz};
  // std::vector<double> pC = {0+centx,0+centy,1+centz};
  // // double pD[] = {pow(2,0.5)/2.,pow(2,0.5)/2.,0};
  // std::vector<double> pD = {-1+centx,0+centy,0+centz};

  // TRandom3 rng(0);

  // double rx,ry,rz,rr;
  // for (int ii=0; ii<10; ii++) {
  //   cout << "ON I " << ii << endl;
  //   rx = rng.Uniform(-1.0, 1.0);
  //   ry = rng.Uniform(-1.0, 1.0);
  //   rz = rng.Uniform(-1.0, 1.0);

  //   rr = rx*rx + ry*ry + rz*rz;

  //   std::vector<double> P = {rx+centx, ry+centy, rz+centz};
  //   int isin = IsInCircumsphere(pA, pB, pC, pD, P);

  //   cout << rr << " " << isin << endl;
  // }

  // exit(1);

  std::vector<std::array<double, 4>> nucleonPositions = ini->GetTargetNucleonPositions();
  std::vector<int> nucleonPids; //TODO: we currently just make these up, and keep all the nucleons...

  int pidi=0;
  for (; pidi<targZ; pidi++) { nucleonPids.push_back(2212); }
  for (; pidi<targA; pidi++) { nucleonPids.push_back(2112); }
  random_shuffle(nucleonPids.begin(), nucleonPids.end());

  std::vector<Sphere> nucleonSpheres;
  std::vector<Sphere> nuclearHoles;

  if (shiftHadrons) { 
    for (int i=0; i<nucleonPositions.size(); i++) {
      nucleonSpheres.push_back(InitSphere(nucleonPositions[i][1], nucleonPositions[i][2], nucleonPositions[i][3], nucleonPids[i]));
    }

    if (shiftType==VORONOI_SHIFT_TYPE) {
      nuclearHoles = VoronoiHoles(nucleonSpheres); //empty spaces inside
    }
  }
  // cout << "num holes " << nuclearHoles.size() << endl;

  Event &event = pythia.event;
  event.reset();
  double pz = p_fake;

  //initial hadron list
  //pre existing hadrons
  if(hOut.size() > 0){
    for(vector<shared_ptr<Hadron>>::iterator hadIter = hOut.begin(); hadIter<hOut.end();){
      if(hadIter->get()->pid() > 23){ //skipping leptons
        double massnow = hadIter->get()->e()*hadIter->get()->e() -
                        (hadIter->get()->px()*hadIter->get()->px() + hadIter->get()->py()*hadIter->get()->py() + hadIter->get()->pz()*hadIter->get()->pz());
        massnow = (massnow >= 0.) ? sqrt(massnow) : -sqrt(-massnow);
        event.append(hadIter->get()->pid(),0,0,0,hadIter->get()->px(),hadIter->get()->py(),hadIter->get()->pz(),hadIter->get()->e(),massnow);
        event[event.size()-1].vProd(0., 0., 0., 0.);
      }
      hadIter++;
    }
  }

  int pythiastat;
  JSDEBUG << "&&&&&&&&&&&&&&&&&&& the number of showers are: " << shower.size();
  for (unsigned int ishower = 0; ishower < shower.size(); ++ishower) {
    JSDEBUG << "&&&&&&&&&&&&&&&&&&& there are " << shower.at(ishower).size()
            << " partons in the shower number " << ishower;
    for (unsigned int ipart = 0; ipart < shower.at(ishower).size(); ++ipart) {
      double onshellE = pow(pow(shower.at(ishower).at(ipart)->px(), 2) +
                                pow(shower.at(ishower).at(ipart)->py(), 2) +
                                pow(shower.at(ishower).at(ipart)->pz(), 2),
                            0.5);

      if (shower.at(ishower).at(ipart)->pid() == 22) {

        VERBOSE(1) << BOLDYELLOW
                   << " photon found in colored hadronization with ";
        VERBOSE(1) << BOLDYELLOW
                   << "px = " << shower.at(ishower).at(ipart)->px();
        //cin >> blurb;
      }

      if (shower.at(ishower).at(ipart)->plabel() < 0) { pythiastat = 63; } //beam remnant
      else { pythiastat = 23; } //outgoing particle

      event.append(shower.at(ishower).at(ipart)->pid(), pythiastat,
                   shower.at(ishower).at(ipart)->color(),
                   shower.at(ishower).at(ipart)->anti_color(),
                   shower.at(ishower).at(ipart)->px(),
                   shower.at(ishower).at(ipart)->py(),
                   shower.at(ishower).at(ipart)->pz(), onshellE);
      event[event.size() - 1].vProd(shower.at(ishower).at(ipart)->x_in().comp(1)*Pythia8::FM2MM,
                                    shower.at(ishower).at(ipart)->x_in().comp(2)*Pythia8::FM2MM,
                                    shower.at(ishower).at(ipart)->x_in().comp(3)*Pythia8::FM2MM,
                                    shower.at(ishower).at(ipart)->x_in().comp(0)*Pythia8::FM2MM);
      //JSINFO << "parton " << shower.at(ishower).at(ipart)->pid() << " location from jetscape: " << event[event.size() - 1].vProd()*Pythia8::MM2FM;
    }
  }

  pythia.next();
  // event.list();

  unsigned int ip = hOut.size();
  Sphere thisSphere, nucSphere;
  double dnuc;
  for (unsigned int eventi = 0; eventi < event.size(); eventi++) {
    // cout << "ON " << eventi << " OF " << event.size() << endl;
    if (!event[eventi].isFinal())
      continue;
    //if ( !event[eventi].isHadron() )  continue;
    if (fabs(event[eventi].eta()) > 20)
      continue; //To prevent "nan" from propagating, very rare though

    double x[4] = {event[eventi].vProd().e()*Pythia8::MM2FM, 
            event[eventi].vProd().px()*Pythia8::MM2FM, 
            event[eventi].vProd().py()*Pythia8::MM2FM, 
            event[eventi].vProd().pz()*Pythia8::MM2FM};

    //JSINFO << "hadron " << event[eventi].id() << " location from pythia: " << event[eventi].vProd()*Pythia8::MM2FM;

    // cout << "ORIG POS " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;

    if (shiftHadrons) {

      //check if we are inside any nucleon
      if (!isOverlap(x, event[eventi].id(), nucleonSpheres)) continue;

      //it overlaps with a nucleon, have to move it
      // cout << "\toverlap" << endl;

      if (shiftType==VORONOI_SHIFT_TYPE) {
        //options are: the closest hole big enough to take us, or the outside of the nucleus
        //option 1:
        bool foundOpen;
        Sphere closestOpen = nuclearHoles[0];
        dnuc = 999999.;
        double thisdist;
        double x_o1[4] = {0.};

        for (int voroi=0; voroi<nuclearHoles.size(); voroi++) {
          thisdist = d2Sphere(nuclearHoles[voroi], thisSphere);
          if (thisdist < dnuc && nuclearHoles[voroi].r >= thisSphere.r) {
            foundOpen = true;
            dnuc = thisdist;
            closestOpen = nuclearHoles[voroi];
          }
        }
        //now find where we have to move it to -- for now just take the center
        //TODO: actually do the geometry problem
        if (foundOpen) {
          x_o1[1] = closestOpen.x;
          x_o1[2] = closestOpen.y;
          x_o1[3] = closestOpen.z;
          cout << "\tfound open space" << endl;
        }
        double d_o1[4] = {0.};
        for (int i=0; i<4; i++) { d_o1[i] = x[i] - x_o1[i]; }
        double d_o1sq = R3dot(d_o1, d_o1);

        // cout << "A" << endl;

        //TODO: THIS USES POSITION LOL NOT VELOCITY
        //option 2: find edge of nuclear medium by going way out and reeling it back in
        double x_o2[4] = {0.};

        double velx[4] = {0., event[eventi].px()/event[eventi].e(), 
                            event[eventi].py()/event[eventi].e(), 
                            event[eventi].pz()/event[eventi].e()};

        double vmag = pow(R3dot(velx,velx), 0.5);

        double ux = velx[1]/vmag;
        double uy = velx[2]/vmag;
        double uz = velx[3]/vmag;
        double mult = 100.;
        while (ini->Get_target_nucleon_density_lab(0., mult*ux, mult*uy, mult*uz) == 0) { mult -= 0.01; } //TODO: THIS ONLY SEES GLAUBER
        mult += (thisSphere.r + 0.01);
        x_o2[1] = x[1] + mult*ux;
        x_o2[2] = x[2] + mult*uy;
        x_o2[3] = x[3] + mult*uz;

        double d_o2[4] = {0.};
        for (int i=0; i<4; i++) { d_o2[i] = x[i] - x_o2[i]; }
        double d_o2sq = R3dot(d_o2, d_o2);

        // cout << "B" << endl;

        if (foundOpen) { //we have two options, see which one is closer
          if (d_o1sq<d_o2sq) {
            x[0] += pow(d_o1sq, 0.5);
            for (int i=1; i<4; i++) { x[i] = x_o1[i]; }
            cout << "\tgap closer" << endl;
          }
          else {
            x[0] += pow(d_o2sq, 0.5);
            for (int i=1; i<4; i++) { x[i] = x_o2[i]; }
            cout << "\toutside closer" << endl;
          }
        }
        else { //have to put it outside
          x[0] += pow(d_o2sq, 0.5);
          for (int i=1; i<4; i++) { x[i] = x_o2[i]; }
          cout << "\tonly have outside" << endl;
        }

        //TODO: add hadron to delaunay instead of recomputing the whole thing every time
        nucleonSpheres.push_back(thisSphere);
        nuclearHoles = VoronoiHoles(nucleonSpheres);
      }

      if (shiftType==MC_SHIFT_TYPE) {
        double rsearch = 1.; //how far out to go from x
        double currdot = -2.; //dot prod of two unit vectors
        bool foundHole = false;

        double velx[4] = {0., event[eventi].px()/event[eventi].e(), 
                            event[eventi].py()/event[eventi].e(), 
                            event[eventi].pz()/event[eventi].e()};

        double vmag = pow(R3dot(velx,velx), 0.5);
        double ux[4] = {velx[0], velx[1]/vmag, velx[2]/vmag, velx[3]/vmag}; //direction of motion
        double urand[4] = {0.};
        double tmpx[4] = {0.};
        double newt, newx, newy, newz;
        double tmpdot;

        while (!foundHole) {
          for (int mci=0; mci<100; mci++) { //number of directions to try
            randUnitR3(urand[1], urand[2], urand[3]);
            for (int i=1; i<4; i++) { tmpx[i] = x[i]+rsearch*urand[i]; } //randomly shifted point

            if (!isOverlap(tmpx, event[eventi].id(), nucleonSpheres)) { //is it empty
              // cout << "\t\tfound hole at r=" << rsearch << endl;

              //possibly reel it in a little bit
              do { 
                rsearch-=0.01;
                for (int i=1; i<4; i++) { tmpx[i] = x[i]+rsearch*urand[i]; }
              }
              while (!isOverlap(tmpx, event[eventi].id(), nucleonSpheres));
              rsearch+=0.01;

              // cout << "\t\tsettled on r=" << rsearch << endl;

              //see if it's closer to the direction of motion than what we have so far
              tmpdot = R3dot(ux,urand);
              if (tmpdot > currdot) {
                foundHole = true;
                // cout << "\t\t!!!better than others so far" << endl;
                currdot = tmpdot;
                newt = x[0] + rsearch/vmag;
                newx = x[1] + rsearch*urand[1];
                newy = x[2] + rsearch*urand[2];
                newz = x[3] + rsearch*urand[3];
              }
            }
          }

          rsearch += 1.; //go further out and try again
        }

        x[0] = newt;
        x[1] = newx;
        x[2] = newy;
        x[3] = newz;
      }
    }

    double hinf = 999999.;
    if (x[0]>hinf || x[1]>hinf || x[2]>hinf || x[3]>hinf) {
      // cout << "???????????????????????????????????????????????????????????????????????????????????????????????????" << endl;
      JSWARN << "Hadron formed at infinity. Setting it to the origin at t=1.";
      x[0] = 1.; x[1] = 0.; x[2] = 0.; x[3] = 0.;
    }

    hOut.push_back(make_shared<Hadron>(ip, event[eventi].id(), event[eventi].status(),
                                       event[eventi].pT(), event[eventi].eta(),
                                       event[eventi].phi(), event[eventi].e(), x));

    // cout << "END  POS " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;
    
    // cout << "PUSHED BACK" << event[eventi].status() << endl;
    ++ip;
  }


  // std::vector<std::array<double, 4>> nucleonPositions = ini->GetTargetNucleonPositions();
  // std::vector<int> nucleonPids; //TODO: we currently just make these up, and keep all the nucleons...

  //push the glauber nucleons
  double thisnuc_m;
  for (int nuci=0; nuci<nucleonPids.size(); nuci++) {
    if (nucleonPids[nuci]==2212) { thisnuc_m = 0.938; } //p
    else { thisnuc_m = 0.939; } //n
    hOut.push_back(make_shared<Hadron>(ip, nucleonPids[nuci], 91,
        FourVector(0., 0., 0., thisnuc_m), 
        FourVector(nucleonPositions[nuci][1], nucleonPositions[nuci][2], nucleonPositions[nuci][3], 0.), 
        thisnuc_m));
    ip++;
    // cout << "ADDED NUCLEON " << thisnuc_m << endl;
  }

  cout << "DONE HADRONIZING" << endl;
  // cout << hOut.size() << endl;
  shower.clear();
  // cout << hOut.size() << endl;
}