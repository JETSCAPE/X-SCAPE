/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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
 * @class HardProcess
 * @brief Interface for the hard scattering process in the JetScape framework.
 *
 * The HardProcess class defines the interface and base functionality
 * for simulating hard scattering events. It provides methods for
 * initialization, execution, clearing state, writing output, and handling event
 * metadata such as cross-sections and weights. Derived classes implement
 * specific physics modules such as Pythia-based event generators.
 */
class HardProcess : public JetScapeModuleBase {
 public:
  /**
   * @brief Default constructor.
   *
   * Creates a HardProcess physics task and sets the task ID as "HardProcess".
   */
  HardProcess();

  /**
   * @brief Destructor.
   *
   * Cleans up parton and hadron lists, disconnects signals, and finalizes
   * the task.
   */
  virtual ~HardProcess();

  /**
   * @brief Initialize the HardProcess module.
   *
   * It reads the input parameters relevant to the hard scattering from
   * the XML file under the name tag <Hard>. Uses JetScapeSingnalManager
   * Instance to retrieve the Initial State Physics information. Calls
   * InitTask(); This explicit call can be used for actual initialization
   * of modules such as @a PythiaGun if attached as a @a polymorphic class.
   * It also initializes the tasks within the current module.
   * @sa Read about @a polymorphism in C++. Override Init (not InitTask)
   * here as sub-tasks are called as well.
   */
  void Init() override;

  /**
   * @brief Execute the HardProcess module.
   *
   * Calls JetScapeTask::ExecuteTasks() for recursive execution of tasks
   * attached to HardProcess module. It can be overridden by the attached
   * module.
   */
  virtual void ExecuteTask();

  /**
   * @brief Clear the state of the HardProcess module.
   *
   * Erases the hard partons stored in the vector @a hp_list of the hard
   * process module. It can be overridden by the attached module.
   */
  virtual void ClearTask();

  /**
   * @brief Write hard process data to output.
   *
   * @param w Weak pointer to a JetScapeWriter module.
   *
   * Outputs hard parton list and related metadata to the writer.
   */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

  /**
   * @brief Collect header information for output.
   *
   * @param w Weak pointer to a JetScapeWriter module.
   *
   * Stores generated cross-sections, pt-hat, and event weights into
   * the output header.
   */
  virtual void CollectHeader(weak_ptr<JetScapeWriter> w);

  // connect the InitialState module with hard process

  /** @brief Pointer to the InitialState module. */
  std::shared_ptr<InitialState> ini;

  /**
      @return The number of hard partons.
   */
  int GetNHardPartons() { return hp_list.size(); }

  /**
   * @brief Get a pointer to a hard parton at a given index.
   * @param i Index of the hard parton.
   * @return Shared pointer to the Parton.
   */
  shared_ptr<Parton> GetPartonAt(int i) { return hp_list[i]; }

  /**
   * @brief Get the list of hard partons.
   *
   * @return A vector of the Parton class. These parton classes correspond to
   * the hard partons.
   */
  vector<shared_ptr<Parton>> &GetPartonList() { return hp_list; }

  /**
   * @brief Add a hard parton to the list.
   *
   * It adds a parton class pointer p into an existing vector of hard Parton
   * class, and increases the vector size by 1.
   *
   * @param p Parton class pointer for a hard parton.
   */
  void AddParton(shared_ptr<Parton> p) { hp_list.push_back(p); }

  // Slots ...

  /**
   * @brief Get the list of hard partons.
   *
   * This function stores the vector of hard partons into a vector plist.
   *
   * @param plist A output vector of Parton class.
   */
  void GetHardPartonList(vector<shared_ptr<Parton>> &plist) { plist = hp_list; }

  /**
   * @brief Get the number of parton showers.
   * @return The number of parton showers.
   */
  int GetNPartonShowers() { return ps_list.size(); }

  /**
   * @brief Get a pointer to a parton shower at a given index.
   * @param i Index of the parton shower.
   * @return Shared pointer to the PartonShower.
   */
  shared_ptr<PartonShower> GetPartonShowerAt(int i) { return ps_list[i]; }

  /**
   * @brief Get the list of parton showers.
   *
   * @return A vector of the PartonShower class. These parton showers are
   * generated by the hard process module.
   */
  vector<shared_ptr<PartonShower>> &GetPartonShowerList() { return ps_list; }

  /**
   * @brief Add a parton shower to the list.
   *
   * It adds a parton shower class pointer p into an existing vector of
   * PartonShower class, and increases the vector size by 1.
   *
   * @param p PartonShower class pointer for a parton shower.
   */
  void AddPartonShower(shared_ptr<PartonShower> p) { ps_list.push_back(p); }

  /**
   * @brief Get the list of parton showers.
   *
   * This function stores the vector of parton showers into a vector pslist.
   *
   * @param pslist An output vector of PartonShower class.
   */
  void GetPartonShowerList(vector<shared_ptr<PartonShower>> &pslist) {
    pslist = ps_list;
  }

  /**
   * @brief Get the generated cross section.
   *
   * Generated cross section. To be overwritten by implementations
   * that have such information.
   *
   * @return Generated cross section (default = 1).
   */
  virtual double GetSigmaGen() { return 1; };

  /**
   * @brief Get the generated cross section error.
   *
   * Generated cross section error.
   * To be overwritten by implementations that have such information.
   *
   * @return Generated cross section error (default = 0).
   */
  virtual double GetSigmaErr() { return 0; };

  /**
   * @brief Get the generated pt-hat.
   *
   * Generated pt-hat
   * To be overwritten by implementations that have such information.
   *
   * @return Generated pt-hat (default = 0).
   */
  virtual double GetPtHat() { return 0; };

  /**
   * @brief Get the generated event weight.
   *
   * Generated weight.
   * This is in addition to sigmaGen, e.g. coming from dynamic oversampling.
   * To be overwritten by implementations that have such information.
   *
   * @return Generated event weight (default = 1).
   */
  virtual double GetEventWeight() { return 1; };

  /**
   * @brief Get the list of hadrons.
   *
   * t adds a Hadron class pointer h into an existing vector of Hadron class,
   * and increases the vector size by 1.
   *
   * @param h Hadron class pointer for a hadron.
   */
  void AddHadron(shared_ptr<Hadron> h) { hd_list.push_back(h); }

  // Slots ...

  /**
   * @brief Get the list of hadrons.
   *
   * This function stores the vector of hadrons into a vector hlist.
   *  @param hlist an output vector of Hadron class.
   */
  void GetHadronList(vector<shared_ptr<Hadron>> &hlist) { hlist = hd_list; }

  /**
   * @brief Get the list of hadrons.
   * @return Reference to the internal list of Hadrons.
   */
  vector<shared_ptr<Hadron>> &GetHadronList() { return hd_list; }

  /**
   * @brief Get the number of hadrons.
   * @return Number of hadrons currently stored.
   */
  int GetNHadrons() { return hd_list.size(); }

  /**
   * @brief Get the maximum color of the current shower.
   * @return Maximum color of the current shower.
   */
  double GetMax_ColorPerShower() { return max_colorPerShower; }

  //   void SetMax_ColorPerShower(double col) { max_colorPerShower = col; }

  /**
   * @brief Get the maximum color of the whole system.
   * @return Maximum color of the whole system.
   */
  double GetMax_Color() { return max_color; }

  /**
   * @brief Set the maximum color of the whole system.
   * @param col Maximum color to set for the whole system.
   */
  void SetMax_Color(double col) { max_color = col; }

  /**
   * @brief Get the number of ISR showers.
   * @return Number of ISR showers.
   */
  double GetNISRShower() { return NISRShower; }

  /**
   * @brief Set the number of ISR showers.
   * @param NISR Number of ISR showers to set.
   */
  void SetNISRShower(double NISR) { NISRShower = NISR; }

  /**
   * @brief Get the total momentum fraction for the positive side.
   * @return Total momentum fraction for the positive side.
   */
  double GetTotalMomentumFractionPositive() {
    return TotalMomentumFractionPositive;
  }

  /**
   * @brief Get the total momentum fraction for the negative side.
   * @return Total momentum fraction for the negative side.
   */
  double GetTotalMomentumFractionNegative() {
    return TotalMomentumFractionNegative;
  }

  /**
   * @brief Set the total momentum fraction for the positive side.
   * @param mom Total momentum fraction to set for the positive side.
   */
  void SetTotalMomentumFractionPositive(double mom) {
    TotalMomentumFractionPositive = mom;
  }

  /**
   * @brief Set the total momentum fraction for the negative side.
   * @param mom Total momentum fraction to set for the negative side.
   */
  void SetTotalMomentumFractionNegative(double mom) {
    TotalMomentumFractionNegative = mom;
  }

  /**
   * @brief Get the total momentum for the positive side.
   * @return Total momentum for the positive side.
   */
  double GetTotalMomentumPositive() { return TotalMomentumPositive; }

  /**
   * @brief Get the total momentum for the negative side.
   * @return Total momentum for the negative side.
   */
  double GetTotalMomentumNegative() { return TotalMomentumNegative; }

  /**
   * @brief Set the total momentum for the positive side.
   * @param mom Total momentum to set for the positive side.
   */
  void SetTotalMomentumPositive(double mom) { TotalMomentumPositive = mom; }

  /**
   * @brief Set the total momentum for the negative side.
   * @param mom Total momentum to set for the negative side.
   */
  void SetTotalMomentumNegative(double mom) { TotalMomentumNegative = mom; }

  std::string printer;

  /**
   * @brief Get the list of remnants.
   * @return A vector of Parton objects representing the remnants.
   */
  std::vector<Parton> GetRemnants() { return Remnants; }

  /**
   * @brief Add a remnant parton to the list.
   * @param par Parton object representing the remnant to be added.
   */
  void PushRemnants(Parton par) { Remnants.push_back(par); }

 private:
  // Think of always using unique_ptr for any vector in jetscape framework !???
  // To be discussed ...

  /** @brief List of hard partons. */
  vector<shared_ptr<Parton>> hp_list;
  vector<shared_ptr<PartonShower>> ps_list;

  // A vector of Hadrons generated by Pythia

  /** @brief List of hadrons generated by the hard process. */
  vector<shared_ptr<Hadron>> hd_list;
  std::vector<Parton> Remnants;

  const int max_colorPerShower = 1000;
  int max_color, NISRShower = 0;
  double TotalMomentumFractionPositive, TotalMomentumFractionNegative;
  double TotalMomentumPositive, TotalMomentumNegative;
};

}  // end namespace Jetscape

#endif
