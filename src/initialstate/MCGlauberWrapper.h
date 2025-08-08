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

#ifndef MCGlauberWrapper_H
#define MCGlauberWrapper_H

#include <memory>

#include "JetScapeModuleBase.h"
#include "InitialState.h"
#include "JetScapeLogger.h"
#include <MakeUniqueHelper.h>
#include "InitialState.h"
#include "Parameters.h"
#include "EventGenerator.h"
using namespace Jetscape;

class MCGlauberWrapper : public Jetscape::InitialState {
  // this is wrapper class to read external files that
  // stores initial number of binary collisions and corresponding
  // configurations

public:
  MCGlauberWrapper();
  ~MCGlauberWrapper() {}

  /** Reads the input parameters from the XML file under the tag  <IS>. Calls InitTask(); This explicit call of InitTask() can be used for actual initialization of modules such as @a Trento if attached as a @a polymorphic class. It also initializes the tasks within the current module.
      @sa Read about @a polymorphism in C++.
   */
  //void Init();

  /** Default ExecuteTask() function. It can be overridden by other tasks.
   */
  void ExecuteTask();

  /** Default ClearTask() function. It can be overridden by other tasks.
   */
  void ClearTask();

  void InitTask();

  /** Default Write() function. It can be overridden by other tasks.
      @param w A pointer to the JetScapeWriter class.
   */
  virtual void Write(weak_ptr<JetScapeWriter> w) {}

  /** Generated number of binary collisions.
  */
  double GetNcoll() { return(static_cast<double>(ncoll_)); }

  void SampleABinaryCollisionPoint(double &t, double &x,
                                   double &y, double &z);
  double Get_total_nucleon_density_lab(double t, double x,
                                       double y, double z);
  double Get_target_nucleon_density_lab(double t, double x,
                                        double y, double z);
  double Get_projectile_nucleon_density_lab(double t, double x,
                                            double y, double z);
  void OutputHardCollisionPosition(double t, double x, double y, 
                                   double z);
  void OutputHardPartonMomentum(double E, double px, double py, double pz,
                                int direction, double P_A);
  void ClearHardPartonMomentum();
  void GetHardPartonPosAndMomentumProj();
  void GetHardPartonPosAndMomentumTarg();
  std::vector<double> Get_projectile_nucleon_z_lab();
  std::vector<double> Get_target_nucleon_z_lab();
  std::vector<double> Get_quarks_pos_proj_lab();
  std::vector<double> Get_quarks_pos_targ_lab();
  std::vector<double> Get_remnant_proj();
  std::vector<double> Get_remnant_targ();
  void GenerateStrings();
  std::vector< std::vector<double> > GetQCDStringList() {
    return(QCDStringList_);
}

double Get_total_hard_e() {return (proj_parton_e_ + targ_parton_e_);}
std::shared_ptr<InitialState> ini;

private:
  std::shared_ptr<MCGlb::EventGenerator> mc_gen_;
  std::vector<double> binary_collision_t_;
  std::vector<double> binary_collision_x_;
  std::vector<double> binary_collision_y_;
  std::vector<double> binary_collision_z_;
  std::vector< std::vector<double> > QCDStringList_;
  double hard_parton_x_, hard_parton_y_, hard_parton_z_, hard_parton_t_;
  double targ_parton_px_ = 0.0, targ_parton_py_ = 0.0;
  double targ_parton_pz_ = 0.0, targ_parton_e_ = 0.0;
  double proj_parton_px_ = 0.0, proj_parton_py_ = 0.0;
  double proj_parton_pz_ = 0.0, proj_parton_e_ = 0.0;
  std::shared_ptr<MCGlb::RandomUtil::Random> ran_gen_ptr_;
  std::shared_ptr<std::uniform_int_distribution<int>> rand_int_ptr_;
  int ncoll_ = -1;
  bool generateOnlyPositions_ = false;
  // Allows the registration of the module so that it is available to be
  // used by the Jetscape framework.
  static RegisterJetScapeModule<MCGlauberWrapper> reg;
};

#endif  // MCGlauberWrapper_H
