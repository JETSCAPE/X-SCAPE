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

#ifndef HARDPROCESS_H
#define HARDPROCESS_H

#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "JetClass.h"
#include <vector>

namespace Jetscape {

/**
     @class
     Interface for the hard process.
   */
class HardProcess : public JetScapeModuleBase {

public:
  /** Default constructor to create a Hard Process Physics task. Sets the task ID as "HardProcess".
  */
  HardProcess();

  /** Destructor for the Hard Process Physics task.
   */
  virtual ~HardProcess();

  /** It reads the input parameters relevant to the hard scattering from the XML file under the name tag <Hard>. Uses JetScapeSingnalManager Instance to retrieve the Initial State Physics information. Calls InitTask(); This explicit call can be used for actual initialization of modules such as @a PythiaGun if attached as a @a polymorphic class. It also initializes the tasks within the current module.
    @sa Read about @a polymorphism in C++. Override Init (not InitTask) here as sub-tasks are called as well.
  */
  void Init() override;

  /** Calls JetScapeTask::ExecuteTasks() for recursive execution of tasks attached to HardProcess module. It can be overridden by the attached module.
   */
  virtual void ExecuteTask();

  /** Erases the hard partons stored in the vector @a hp_list of the hard process module. It can be overridden by the attached module.
  */
  virtual void ClearTask();

  /** It writes the output information obtained from the HardProcess Task into a file.
      @param w is a pointer of type JetScapeWrite class.
  */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

  /** Collect header information for writer modules
      @param w is a pointer of type JetScapeWrite class.
  */
  virtual void CollectHeader(weak_ptr<JetScapeWriter> w);

  // connect the InitialState module with hard process
  /** A pointer of type InitialState class.
   */
  std::shared_ptr<InitialState> ini;

  /**
      @return The number of hard partons.
   */
  int GetNHardPartons() { return hp_list.size(); }

  /** @return A pointer to the Parton class for ith hard parton.
      @param i Index of a vector of the hard parton.
   */
  shared_ptr<Parton> GetPartonAt(int i) { return hp_list[i]; }

  /** @return A vector of the Parton class. These parton classes correspond to the hard partons.
   */
  vector<shared_ptr<Parton>> &GetPartonList() { return hp_list; }

  /** It adds a parton class pointer p into an existing vector of hard Parton class, and increases the vector size by 1.
      @param p Parton class pointer for a hard parton.
   */
  void AddParton(shared_ptr<Parton> p) { hp_list.push_back(p); }

  // Slots ...
  /** This function stores the vector of hard partons into a vector plist.
      @param plist A output vector of Parton class.
   */
  void GetHardPartonList(vector<shared_ptr<Parton>> &plist) { plist = hp_list; }

  int GetNPartonShowers() {return ps_list.size();}
  shared_ptr<PartonShower> GetPartonShowerAt(int i) {return ps_list[i];}
  vector<shared_ptr<PartonShower>>& GetPartonShowerList() {return ps_list;}

  void AddPartonShower(shared_ptr<PartonShower> p) {ps_list.push_back(p);}
  void GetPartonShowerList(vector<shared_ptr<PartonShower>> &pslist) {pslist=ps_list;}

  /** Generated cross section.
      To be overwritten by implementations that have such information.
  */
  virtual double GetSigmaGen() { return 1; };
  /** Generated cross section error.
      To be overwritten by implementations that have such information.
  */
  virtual double GetSigmaErr() { return 0; };

  /** Generated pt-hat
      To be overwritten by implementations that have such information.
  */
  virtual double GetPtHat() { return 0; };

  /** Generated weight.
      This is in addition to sigmaGen, e.g. coming from dynamic oversampling.
      To be overwritten by implementations that have such information.
  */
  virtual double GetEventWeight() { return 1; };

  /** It adds a Hadron class pointer h into an existing vector of Hadron class, and increases the vector size by 1.
      @param h Hadron class pointer for a hadron.
   */
  void AddHadron(shared_ptr<Hadron> h) { hd_list.push_back(h); }

  // Slots ...
  /** This function stores the vector of hadrons into a vector hlist.
      @param hlist an output vector of Hadron class.
   */
  void GetHadronList(vector<shared_ptr<Hadron>> &hlist) { hlist = hd_list; }

  /** @return A vector of the Hadron class.
   */
  vector<shared_ptr<Hadron>> &GetHadronList() { return hd_list; }

  /**
      @return The number of hadrons.
  */
  int GetNHadrons() { return hd_list.size(); }

  // Get max color of the current shower
  double GetMax_ColorPerShower() { return max_colorPerShower; }
  //   void SetMax_ColorPerShower(double col) { max_colorPerShower = col; }

  // Get max color the whole system  
  double GetMax_Color() { return max_color; }
  void SetMax_Color(double col) { max_color = col; }

  // Get Number of ISR showers
  double GetNISRShower() { return NISRShower; }
  void SetNISRShower(double NISR) { NISRShower = NISR; }


  // Get Total Momentum fraction for the Positive/Negative side 
  double GetTotalMomentumFractionPositive(){ return TotalMomentumFractionPositive;}
  double GetTotalMomentumFractionNegative(){ return TotalMomentumFractionNegative;}
  void SetTotalMomentumFractionPositive(double mom){ TotalMomentumFractionPositive = mom;}
  void SetTotalMomentumFractionNegative(double mom){ TotalMomentumFractionNegative = mom;}
  double GetTotalMomentumPositive(){ return TotalMomentumPositive;}
  double GetTotalMomentumNegative(){ return TotalMomentumNegative;}
  void SetTotalMomentumPositive(double mom){ TotalMomentumPositive = mom;}
  void SetTotalMomentumNegative(double mom){ TotalMomentumNegative = mom;}
  std::string printer;

  std::vector<Parton> GetRemnants() { return Remnants; }
  void PushRemnants(Parton par) { Remnants.push_back(par); }

private:
  // Think of always using unique_ptr for any vector in jetscape framework !???
  // To be discussed ...
  vector<shared_ptr<Parton>> hp_list;
  vector<shared_ptr<PartonShower>> ps_list;

  // A vector of Hadrons generated by Pythia
  vector<shared_ptr<Hadron>> hd_list;
  std::vector<Parton> Remnants;
    
  const int max_colorPerShower = 1000; 
  int max_color, NISRShower = 0;
  double TotalMomentumFractionPositive, TotalMomentumFractionNegative;
  double TotalMomentumPositive, TotalMomentumNegative;

};

} // end namespace Jetscape

#endif
