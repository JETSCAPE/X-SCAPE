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
#include <string>
#include <utility>
#include <vector>

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

//! Set when PythiaIsrGun understands <Hard><PythiaGun><pTHatBins>.
#define PYTHIAISRGUN_HAS_PTHAT_BINS 1
//! Set when PythiaIsrGun understands <partonYMax> / <partonYMode> (a rapidity
//! cut on the two outgoing partons of the hard scattering).
#define PYTHIAISRGUN_HAS_PARTON_Y_CUT 1

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

    Optional: a rapidity cut on the two outgoing partons of the hard scattering
    (Pythia status -23, before the ISR, which iMATTER adds later), <partonYMax>
    (e.g. 0.6) and <partonYMode>: "both" both have |y| < partonYMax, "any"
    either of them.  Their pT is equal at this stage, so PythiaGun's "leading"
    is not defined here and refused.  A rejected event is regenerated before
    anything else runs; the cross section is Pythia's sigmaGen times the
    acceptance kept/tried, per window.  Without <partonYMax> nothing is cut.
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

  double partonYMax_ = 0.;  //!< 0: no cut
  std::string partonYMode_ = "both";
  std::vector<long> yTried_, yKept_;  //!< per window: events tested / passing

  void ReadPtHatBins();
  void ReadPartonYCut();
  //! the two outgoing hard partons (status -23) against |y| < partonYMax_
  bool PassesPartonYCut(const std::vector<Pythia8::Particle> &partons) const;
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
  // window (with <partonYMax>: of the events passing the cut).
  double GetSigmaGen() { return GetSigmaGen(activeBin_); };
  double GetSigmaErr() { return GetSigmaErr(activeBin_); };
  double GetPtHat() { return InfoOf(activeBin_).pTHat(); };
  double GetEventWeight() { return InfoOf(activeBin_).weight(); };

  // pTHat windows (<pTHatBins>); a single window without it.
  int GetNPtHatBins() const { return static_cast<int>(pTHatBins_.size()); }
  int GetActivePtHatBin() const { return activeBin_; }
  double GetPtHatBinMin(int bin) const { return pTHatBins_.at(bin).first; }
  double GetPtHatBinMax(int bin) const { return pTHatBins_.at(bin).second; }
  unsigned int GetPythiaSeed(int bin) const { return seedBin_.at(bin); }
  //! one window's cross section so far [mb]: Pythia's estimate, times the
  //! acceptance of <partonYMax> if set (then the error includes the
  //! acceptance's binomial error)
  double GetSigmaGen(int bin) {
    const double s = InfoOf(bin).sigmaGen();
    return partonYMax_ > 0 ? s * GetYAcceptance(bin) : s;
  }
  double GetSigmaErr(int bin);
  //! Pythia's own estimate, including the events the cut rejected
  double GetSigmaGenRaw(int bin) { return InfoOf(bin).sigmaGen(); }
  double GetSigmaErrRaw(int bin) { return InfoOf(bin).sigmaErr(); }
  //! events Pythia accepted in one window so far (incl. those this gun
  //! rejected and regenerated, which also enter sigmaGen)
  long GetNAccepted(int bin) { return InfoOf(bin).nAccepted(); }

  // Rapidity cut on the outgoing hard partons (<partonYMax>; 0 = none).
  double GetPartonYMax() const { return partonYMax_; }
  std::string GetPartonYMode() const { return partonYMode_; }
  //! events that reached / passed the cut in one window (0 without it)
  long GetNYTried(int bin) const {
    return yTried_.empty() ? 0 : yTried_.at(bin);
  }
  long GetNYKept(int bin) const { return yKept_.empty() ? 0 : yKept_.at(bin); }
  //! kept / tried (1 without a cut or before the first event)
  double GetYAcceptance(int bin) const {
    const long n = GetNYTried(bin);
    return n > 0 ? double(GetNYKept(bin)) / n : 1.0;
  }
};

#endif  // PYTHIAISRGUN_H
