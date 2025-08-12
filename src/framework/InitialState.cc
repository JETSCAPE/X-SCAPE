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

#include "InitialState.h"
#include "JetScapeWriter.h"
#include <iostream>

namespace Jetscape {

InitialState::~InitialState() {}

void InitialState::Init() {
  JetScapeModuleBase::InitTask();

  JSINFO << "Initialize InitialState ... " << GetId() << " ...";

  grid_max_x_ = GetXMLElementDouble({"IS", "grid_max_x"});
  grid_max_y_ = GetXMLElementDouble({"IS", "grid_max_y"});
  grid_max_z_ = GetXMLElementDouble({"IS", "grid_max_z"});
  grid_step_x_ = GetXMLElementDouble({"IS", "grid_step_x"});
  grid_step_y_ = GetXMLElementDouble({"IS", "grid_step_y"});
  grid_step_z_ = GetXMLElementDouble({"IS", "grid_step_z"});
  JSINFO << "x range for bulk evolution = [" << -grid_max_x_ << ", "
         << grid_max_x_ << "]";

  InitTask();
  InitTasks();
}

void InitialState::ExecuteTask() {
  // Do whatever is needed to figure out the internal temp...
}

void InitialState::ClearTask() {}

void InitialState::Write(weak_ptr<JetScapeWriter> w) {
  //Write out the original vertex so the writer can keep track of it...
  // auto f = w.lock();
  // if ( f ) f->Write(make_shared<Vertex>(initialVtx));
}

void InitialState::CollectHeader(weak_ptr<JetScapeWriter> w) {
  auto f = w.lock();
  if (f) {
    auto &header = f->GetHeader();
    header.SetNpart(GetNpart());
    header.SetNcoll(GetNcoll());
    header.SetEventCentrality(GetEventCentrality());
    header.SetTotalEntropy(GetTotalEntropy());
  }
}

std::tuple<double, double, double> InitialState::CoordFromIdx(int idx) {
  int nx = GetXSize();
  int ny = GetYSize();
  int nz = GetZSize();

  int ix = idx / (ny * nz);
  int iy = (idx - (ny * nz * ix))/ nz;
  int ieta = idx - (ny * nz * ix) - (nz * iy);

  return std::make_tuple(-grid_max_x_ + ix * grid_step_x_,
                         -grid_max_y_ + iy * grid_step_y_,
                         -grid_max_z_ + ieta * grid_step_z_);
}


void InitialState::SampleABinaryCollisionPoint(double &t, double &x,
                                               double &y, double &z) {
  if (num_of_binary_collisions_.size() == 0) {
    JSWARN << "num_of_binary_collisions is empty, setting the starting "
              "location to 0. Make sure to add e.g. trento before PythiaGun.";
  } else {
    std::discrete_distribution<> dist(
        begin(num_of_binary_collisions_),
        end(num_of_binary_collisions_)); // Create the distribution
    // Now generate values
    auto idx = dist(*GetMt19937Generator());
    auto coord = CoordFromIdx(idx);
    t = 0.0;
    x = std::get<0>(coord);
    y = std::get<1>(coord);
    z = 0.0;
  }
}

void InitialState::OutputHardCollisionPosition(double t, double x, double y, 
                                                                   double z) {}

void InitialState::OutputHardPartonMomentum(double E, double px, double py, double pz,
                                            int direction, double P_A) {}

void InitialState::ClearHardPartonMomentum() {}

  
void InitialState::GetHardPartonPosAndMomentumProj() {}

void InitialState::GetHardPartonPosAndMomentumTarg() {}

std::vector<double> InitialState::Get_projectile_nucleon_z_lab() {
    std::vector<double> Temp;
    for (int i = 0; i != 8; i++) {
        Temp.push_back(-1.);
    }
    return Temp;
}

std::vector<double> InitialState::Get_target_nucleon_z_lab() {
    std::vector<double> Temp;
    for (int i = 0; i != 8; i++) {
        Temp.push_back(-1.);
    }
    return Temp;
}

std::vector<double> InitialState::Get_quarks_pos_proj_lab() {
    std::vector<double> Temp;
    for (int i = 0; i != 9; i++) {
        Temp.push_back(-1.);
    }
    return Temp;
}

std::vector<double> InitialState::Get_quarks_pos_targ_lab() {
    std::vector<double> Temp;
    for (int i = 0; i != 9; i++) {
        Temp.push_back(-1.);
    }
    return Temp;
}

std::vector<double> InitialState::Get_remnant_proj() {
    std::vector<double> Temp;
    for (int i = 0; i != 4; i++) {
        Temp.push_back(-1.);
    }
    return Temp;
}

std::vector<double> InitialState::Get_remnant_targ() {
    std::vector<double> Temp;
    for (int i = 0; i != 4; i++) {
        Temp.push_back(-1.);
    }
    return Temp;
}

std::vector<double> InitialState::Get_Proj_Remnant() {
    std::vector<double> Temp;
    for (int i = 0; i != 4; i++) {
        Temp.push_back(-1.);
    }
    return Temp;
}

std::vector<double> InitialState::Get_Targ_Remnant() {
    std::vector<double> Temp;
    for (int i = 0; i != 4; i++) {
        Temp.push_back(-1.);
    }
    return Temp;
}

void InitialState::GenerateStrings() {
  // Do whatever is needed to figure out the internal temp...
  std::cout<<"Call the wrong GenerateStrings function..."<<std::endl;
}

} // end namespace Jetscape
