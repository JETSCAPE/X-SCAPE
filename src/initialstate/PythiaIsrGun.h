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

// Create a pythia collision at a specified point and return the two inital hard
// partons

#ifndef PYTHIAISRGUN_H
#define PYTHIAISRGUN_H

#include <memory>
#include <utility>
#include <vector>

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

//! Set when PythiaIsrGun understands <Hard><PythiaGun><pTHatBins>.
#define PYTHIAISRGUN_HAS_PTHAT_BINS 1

using namespace Jetscape;

/** Pythia8 hard process for the initial-state radiation (ISR) chain.

    Optional: several pTHat windows in one run, <Hard><PythiaGun><pTHatBins>
    ("min max min max ..."), as for PythiaGun: every window has its own Pythia
    instance (own seed and cross section), and event i uses window i mod K.
    Unlike PythiaGun, every event is a complete event: with ISR the strings, and
    so the medium, are made from what this event's hard scattering left, so the
    hydro is never reused (JetScape::Init refuses <setReuseHydro> with ISR).
    Without <pTHatBins> the only window is pTHatMin-pTHatMax, as before.

    Window 0 is this object itself, so the inherited Pythia members describe
    window 0 only; the getters below describe the window of the current event.
*/
class PythiaIsrGun : public HardProcess, public Pythia8::Pythia {
 private:
  double pTHatMin;  //!< of the current window
  double pTHatMax;  //!< of the current window
  double eCM;
  bool FSR_on;

  //! (pTHatMin, pTHatMax) per window; one entry without <pTHatBins>
  std::vector<std::pair<double, double>> pTHatBins_;
  std::vector<unsigned int> seedBin_;
  //! Pythia instances of windows 1..K-1 (window 0 is *this)
  std::vector<std::unique_ptr<Pythia8::Pythia>> extraPythia_;
  int activeBin_ = 0;

  void ReadPtHatBins();
  void ConfigurePythia(Pythia8::Pythia &py, int bin);
  Pythia8::Pythia &PythiaOf(int bin);
  const Pythia8::Info &InfoOf(int bin);

  // Allows the registration of the module so that it is available to be used by
  // the Jetscape framework.
  static RegisterJetScapeModule<PythiaIsrGun> reg;

 public:
  int NumberOfHardPartons;
  /** standard ctor
      @param xmlDir: Note that the environment variable PYTHIA8DATA takes
     precedence! So don't use it.
      @param printBanner: Suppress starting blurb. Should be set to true in
     production, credit where it's due
  */
  PythiaIsrGun(string xmlDir = "DONTUSETHIS", bool printBanner = false)
      : Pythia8::Pythia(xmlDir, printBanner), HardProcess() {
    SetId("UninitializedPythiaIsrGun");
  }

  ~PythiaIsrGun();

  void InitTask();
  void ExecuteTask();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

  // Getters (pTHat range of the current event's window)
  double GetpTHatMin() const { return pTHatMin; }
  double GetpTHatMax() const { return pTHatMax; }

  // Cross-section information in mb and event weight, of the current event's
  // window.
  double GetSigmaGen() { return InfoOf(activeBin_).sigmaGen(); };
  double GetSigmaErr() { return InfoOf(activeBin_).sigmaErr(); };
  double GetPtHat() { return InfoOf(activeBin_).pTHat(); };
  double GetEventWeight() { return InfoOf(activeBin_).weight(); };

  // pTHat windows (<pTHatBins>); a single window without it.
  int GetNPtHatBins() const { return static_cast<int>(pTHatBins_.size()); }
  int GetActivePtHatBin() const { return activeBin_; }
  double GetPtHatBinMin(int bin) const { return pTHatBins_.at(bin).first; }
  double GetPtHatBinMax(int bin) const { return pTHatBins_.at(bin).second; }
  unsigned int GetPythiaSeed(int bin) const { return seedBin_.at(bin); }
  //! one window's cross section so far [mb] and its error
  double GetSigmaGen(int bin) { return InfoOf(bin).sigmaGen(); }
  double GetSigmaErr(int bin) { return InfoOf(bin).sigmaErr(); }
  //! events Pythia accepted in one window so far
  long GetNAccepted(int bin) { return InfoOf(bin).nAccepted(); }
};

#endif  // PYTHIAISRGUN_H
