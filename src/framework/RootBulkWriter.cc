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
// -----------------------------------------
// JETSCAPE module for soft particlization
// This module will generate Monte-Carlo samples for soft hadrons
// -----------------------------------------
//#ifdef USE_ROOT

#include "RootBulkWriter.h"
#include <iostream>
#include <time.h>
#include <string>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#include "JetScapeSignalManager.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
//#include "JetScapeWriterRootHepMC.h"
#endif


// User modules derived from jetscape framework clasess
#include "TrentoInitial.h"
#include "AdSCFT.h"
#include "Matter.h"
#include "LBT.h"
#include "Martini.h"
#include "Brick.h"
#include "GubserHydro.h"
#include "MusicWrapper.h"
#include "PythiaGun.h"
#include "iSpectraSamplerWrapper.h"
#include "TrentoInitial.h"
#include "NullPreDynamics.h"
#include "PGun.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"
#include "ColorlessHadronization.h"
//#include "HydroFromFile.h"

#include <chrono>
#include <thread>

#include "TParameter.h"
#include <Riostream.h>
#include "TRandom.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

using namespace std;
using namespace Jetscape;

RegisterJetScapeModule<RootBulkWriter> RootBulkWriter::reg("RootBulkWriter");

// Forward declarations
// / -------------------------------------

void Show();

// -------------------------------------
void RootBulkWriter::Init() {

    JSINFO << "Initialzing RootBulkWriter ...";
    out_file_name = GetXMLElementText({"RootBulkWriter","out_file_name"});

    x_min = GetXMLElementDouble({"RootBulkWriter","x_min"});
    dx = GetXMLElementDouble  ({"RootBulkWriter","dx"});
    y_min = GetXMLElementDouble({"RootBulkWriter","y_min"});
    dy = GetXMLElementDouble  ({"RootBulkWriter","dy"});
    tau_min = GetXMLElementDouble({"RootBulkWriter","tau_min"});
    dtau = GetXMLElementDouble  ({"RootBulkWriter","dtau"});
    eta_min = GetXMLElementDouble({"RootBulkWriter","eta_min"});
    deta = GetXMLElementDouble  ({"RootBulkWriter","deta"});
    ntau = GetXMLElementInt  ({"RootBulkWriter","ntau"});

    bool ensure_MUSIC = (int)(GetXMLElementInt(
        {"RootBulkWriter", "ensure_MusicWrapper_output"}, false));
    if (!ensure_MUSIC) {
      JSWARN << " WARNING: RootBulkWriter: ensure_MusicWrapper_output not set "
                "to true. This will likely lead to erroneous output. Please "
                "fix your XML file.";
      JSINFO << " WARNING: RootBulkWriter: ensure_MusicWrapper_output not set "
                "to true. This will likely lead to erroneous output. Please "
                "fix your XML file.";
    }

    JSINFO << " RootBulkWriter initialized with output file name: " << out_file_name;
    JSINFO << " RootBulkWriter initialized with grid parameters: x_min = " << x_min << ", dx = " << dx << ", y_min = " << y_min << ", dy = " << dy << 
              ", tau_min = " << tau_min << ", dtau = " << dtau << ", eta_min = " << eta_min << ", deta = " << deta << ", ntau = " << ntau;
}

RootBulkWriter::RootBulkWriter() 
  // The information to init (the size of music, etc...) isn't present until
  // after the first Exec() call to MUSIC, so init is done in first Exec() instead.
{
    //JSINFO << " Adding RootBulkWriter ";
    SetId("RootBulkWriter");
} 

void RootBulkWriter::Show()
{
  INFO_NICE<<"------------------------------------------";
  INFO_NICE<<"| Bulk ROOT Writer JetScape Framework ... |";
  INFO_NICE<<"------------------------------------------";
  INFO_NICE;
}

void RootBulkWriter::init_tree(const EvolutionHistory& bInfo) {

  isinit=true;

  f=new TFile(out_file_name.c_str(),"RECREATE");
  t=new TTree("t","Tree");

  // save MUSIC parameters to the ROOT output file
  int nX_MUSIC = bInfo.nx;
  float dX_MUSIC = bInfo.dx;
  float X_min_MUSIC = bInfo.x_min;

  int nY_MUSIC = bInfo.nx;
  float dY_MUSIC = bInfo.dx;
  float Y_min_MUSIC = bInfo.x_min;

  int neta_MUSIC = bInfo.neta;
  float deta_MUSIC = bInfo.deta;
  float eta_min_MUSIC = bInfo.eta_min;;

  /* int ntau_MUSIC = bInfo.ntau; */
  dtau_MUSIC = bInfo.dtau;
  float tau_min_MUSIC = bInfo.tau_min;;

  // use parameters from xml
  if (!x_min) x_min = X_min_MUSIC;
  if (!dx) dx = dX_MUSIC;
  if (!y_min) y_min = Y_min_MUSIC;
  if (!dy) dy = dY_MUSIC;
  if (!tau_min) tau_min = tau_min_MUSIC;
  if (!dtau) dtau = dtau_MUSIC;
  /* if (!ntau) ntau = ntau_MUSIC; */
  if (!eta_min) eta_min = eta_min_MUSIC;
  if (!deta) deta = deta_MUSIC;\

  // assign to member variables
  nx = 2*int(fabs(x_min)/dx)+1;
  ny = 2*int(fabs(y_min)/dy)+1;
  neta = 2*int(fabs(eta_min)/deta)+1;

  //JP: Do not understand the need for the the +1 !??? And also the ntau !???? Follow up!!!
  //nx = 2*int(fabs(x_min)/dx);
  //ny = 2*int(fabs(y_min)/dy);
  //neta = 2*int(fabs(eta_min)/deta);

  use_vec = (ntau <= 0);
  if (use_vec) {
    JSINFO << " RootBulkWriter writing variable sized vector for tau steps ";
    t->Branch("user_res", &v_data);
  } else {
    JSINFO << " RootBulkWriter writing fixed sized array for " << ntau << " tau steps out to " << (tau_min+ntau*dtau) << " fm/c";
    ntotal = nx*ny*neta*ntau*nFeatures;
    data = std::make_unique<float[]>(ntotal); //new float[ntotal];
    t->Branch("user_res", data.get(), Form("user_res[%d]/F", ntotal));
  }

  JSINFO << " RootBulkWriter initialized with grid parameters: x_min = " << x_min << ", dx = " << dx << ", y_min = " << y_min << ", dy = " << dy << 
    ", tau_min = " << tau_min << ", dtau = " << dtau << ", eta_min = " << eta_min << ", deta = " << deta << ", ntau = " << ntau;
  JSINFO<<" neta = " << neta << " eta_min = " << eta_min << " deta = " << deta;
  JSINFO << " MUSIC grid parameters: x_min = " << X_min_MUSIC << ", dx = " << dX_MUSIC << ", y_min = " << Y_min_MUSIC << ", dy = " << dY_MUSIC << 
    ", tau_min = " << tau_min_MUSIC << ", dtau = " << dtau_MUSIC << ", eta_min = " << eta_min_MUSIC << ", deta = " << deta_MUSIC << ", ntau = " << bInfo.ntau;

  t->Branch("tau_freezeout", &tau_freezeout, "tau_freezeout/F");
  t->Branch("ntau_freezeout", &ntau_freezeout, "ntau_freezeout/I");

  TParameter<int> p_use_vec ("use_vec", (bool)(use_vec));
  p_use_vec.Write();

  // write the results as parameters to the output tree
  for (auto param : vector<std::tuple<bool,string,float>>{
    // the bool is for if it is an integer or not
    {true, "nFeatures", nFeatures},
    {true, "nx", nx},
    {false, "x_min", x_min},
    {false, "dx", dx},
    {true, "ny", ny},
    {false, "y_min", y_min},
    {false, "dy", dy},
    {true, "neta", neta},
    {false, "eta_min", eta_min},
    {false, "deta", deta},
    {true, "ntau", ntau},
    {false, "tau_min", tau_min},
    {false, "dtau", dtau},
    {true, "nX_MUSIC", nX_MUSIC},
    {false, "dX_MUSIC", dX_MUSIC},
    {false, "X_min_MUSIC", X_min_MUSIC},
    {true, "nY_MUSIC", nY_MUSIC},
    {false, "dY_MUSIC", dY_MUSIC},
    {false, "Y_min_MUSIC", Y_min_MUSIC},
    {true, "neta_MUSIC", neta_MUSIC},
    {false, "deta_MUSIC", deta_MUSIC},
    {false, "eta_min_MUSIC", eta_min_MUSIC},
    /* {true, "ntau_MUSIC", ntau_MUSIC}, */
    {false, "dtau_MUSIC", dtau_MUSIC},
    {false, "tau_min_MUSIC", tau_min_MUSIC}}) 
  {
    if (std::get<0>(param)) {
      TParameter<int> p_temp(std::get<1>(param).c_str(), (int)(std::get<2>(param)));
      p_temp.Write();
    } else {
      TParameter<float> p_temp(std::get<1>(param).c_str(), std::get<2>(param));
      p_temp.Write();
    }
  }
}

void RootBulkWriter::Exec() {
  auto hydro = JetScapeSignalManager::Instance()->GetHydroPointer();

  if (!hydro.lock()) {
    JSWARN << " No hydro pointer found for RootBulkWriter. "
           << " Skipping RootBulkWriter::Exec() logic.";
    return;
  }

  auto bInfo = hydro.lock()->get_bulk_info();
  if (!isinit) {
    init_tree(bInfo);
  }

  float _eta_min = bInfo.eta_min;
  float _eta_max = bInfo.EtaMax();
  float _tau_min = bInfo.tau_min;
  float _tau_max = bInfo.TauMax();

  //REMARK JP: Why are the actually different form the prequilibrium values !??? Follow up!
  //Because with strings tau0 = tau_min - string dtau (0.02) but first bin is nonsensical !!!
  //More puzzling when saving 2 or more events, first lower edensity for taubin=0 ... see dave tau0 = 0.58, so + 1 dtau!???

  JSINFO << " tau_min user( " << tau_min << " ) vs MUSIC ( " << _tau_min 
          << " )  tau_max MUSIC ( " << _tau_max << " ) "; 

  //REMAARK JP: Same issues as with init tree, #bins are not conistent ... not too big of a deal here since one 
  // can clean data for training, but still should be fixed for consistency and to avoid confusion ... Follow up!

  // NOTE: FIXME, below should use _tau_min instead of tau_min, this is why I am getting the empty steps past freezeout
  tau_freezeout = _tau_min + bInfo.ntau * dtau_MUSIC;
  ntau_freezeout = int((tau_freezeout - tau_min) / dtau) + 1; // this ends up being about 2 units too large (why?!?)
  
  // in practice, I am getting events with no energy distribution 
  if (use_vec) { // fill in vector of tau times
    v_data.clear();
    v_data.reserve(nx * ny * neta * nFeatures * ntau_freezeout);
    for (int itau = 0; itau < ntau_freezeout; itau++) {
      double tau_In = tau_min + itau * dtau;
      double FIXME_energy_sum = 0.;
      double FIXME_energy_sum_eta0 = 0.;
      /* std::cout << " FIXME tau_In: " << tau_In << std::endl; */
      for (int ix = 0; ix < nx; ix++) {
        double x_In = x_min + ix * dx;
        for (int iy = 0; iy < ny; iy++) {
          double y_In = y_min + iy * dy;
          for (int ieta = 0; ieta < neta; ieta++) {
            double eta_In = eta_min + ieta * deta;
            auto mCell = bInfo.get(tau_In, x_In, y_In, eta_In);
            /* if (x_In==0.0 && y_In==0.0 && eta_In==0.0) { */
              /* std::cout << " FIXME at center: e= " << mCell.energy_density << " -> "; */
                        /* << "  vx= " << mCell.vx */
                        /* << "  vy= " << mCell.vy */
                        /* << "  vz= " << mCell.vz << std::endl; */
            /* } */
            v_data.push_back((float)(mCell.energy_density));
            /* if (eta_In==0) FIXME_energy_sum_eta0 += mCell.energy_density; */
            /* FIXME_energy_sum += mCell.energy_density; */
            v_data.push_back((float)(mCell.vx));
            v_data.push_back((float)(mCell.vy));
            v_data.push_back((float)(mCell.vz));//(mCell.vz));
          } // loop ieta
        } // loop iy
      } // loop ix
      /* std::cout << " FIXME energy sum at tau " << tau_In << " = " << FIXME_energy_sum << "  and only at eta==0: " << FIXME_energy_sum_eta0 << std::endl; */
    } // loop itau
    // done with vector fill
    /* std::cout << " FIXME check: " << v_data.size() << " vs expected " */
              /* << nx * ny * neta * nFeatures * ntau_freezeout << std::endl; */
  } else { // use fixed-size for number of tau times
    int _ntau = std::min(ntau_freezeout, ntau);
    size_t index = 0;

    for (int k = 0; k < _ntau; k++) {
      double tau_In = tau_min + k * dtau;
      for (int ix = 0; ix < nx; ix++) {
        double x_In = x_min + ix * dx;
        for (int iy = 0; iy < ny; iy++) {
          double y_In = y_min + iy * dy;
          for (int ieta = 0; ieta < neta; ieta++) {
            double eta_In = eta_min + ieta * deta;
            auto mCell = bInfo.get(tau_In, x_In, y_In, eta_In);
            data[index++] = (float)(mCell.energy_density);
            data[index++] = (float)(mCell.vx);
            data[index++] = (float)(mCell.vy);
            data[index++] = (float)(mCell.vz);
          } // loop ieta
        } // loop iy
      } // loop ix
        VERBOSE(3) << " RootBulkWriter writing step " << k << " of " << ntau
                 << " rho_sum: " << " tau: " << tau_In;
    } // loop itau
    // pad with zeros for time past freezeout, if needed
    for (int i = index; i < ntotal; i++) {
      data[index++] = 0.;
    }
  } // end fixed arrays
  t->Fill();
}

RootBulkWriter::~RootBulkWriter() {
  f->cd();
  f->ls();
  t->Print();
  f->Write();
  f->Close();
}
//#endif // USE_ROOT
