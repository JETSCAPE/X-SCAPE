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

#ifndef PYTHIAGUN_H
#define PYTHIAGUN_H

#include <memory>
#include <utility>
#include <vector>

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

//! Set when PythiaGun understands <Hard><PythiaGun><pTHatBins> (several pTHat
//! windows, one Pythia instance each), so code built against it can check.
#define PYTHIAGUN_HAS_PTHAT_BINS 1

using namespace Jetscape;

/** Pythia8 hard process.

    Optional: several pTHat windows in one run, <Hard><PythiaGun><pTHatBins>,
    e.g. "20 40 50 70" for the windows 20-40 and 50-70 GeV.  Every window has
    its own Pythia instance (settings as for pTHatMin/pTHatMax, its own seed and
    cross section), and event i uses window i mod K.  With <setReuseHydro> and
    nReuseHydro a multiple of K, every hydro event gets nReuseHydro / K jets in
    each window.  Without <pTHatBins> the only window is pTHatMin-pTHatMax and
    everything is as before.

    Window 0 is this object itself, so the inherited Pythia members (info,
    event, stat(), ...) describe window 0 only; the getters below describe the
    window of the current event.
*/
class PythiaGun : public HardProcess, public Pythia8::Pythia {
 private:
  double pTHatMin;  //!< of the current window
  double pTHatMax;  //!< of the current window
  double eCM;
  double vir_factor;
  bool initial_virtuality_pT;
  double softMomentumCutoff;
  bool FSR_on;
  bool softQCD;  //!< of the current window

  //! (pTHatMin, pTHatMax) per window; one entry without <pTHatBins>
  std::vector<std::pair<double, double>> pTHatBins_;
  std::vector<bool> softQCDBin_;
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
  static RegisterJetScapeModule<PythiaGun> reg;

 public:
  /** standard ctor
      @param xmlDir: Note that the environment variable PYTHIA8DATA takes
     precedence! So don't use it.
      @param printBanner: Suppress starting blurb. Should be set to true in
     production, credit where it's due
  */
  PythiaGun(string xmlDir = "DONTUSETHIS", bool printBanner = false)
      : Pythia8::Pythia(xmlDir, printBanner), HardProcess() {
    SetId("UninitializedPythiaGun");
  }

  ~PythiaGun();

  void InitTask();
  void ExecuteTask();

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
  //! Pythia's cross-section estimate of one window so far [mb]
  double GetSigmaGen(int bin) { return InfoOf(bin).sigmaGen(); }
  double GetSigmaErr(int bin) { return InfoOf(bin).sigmaErr(); }
  //! events Pythia accepted in one window so far (incl. those this gun
  //! rejected and regenerated, which also enter sigmaGen)
  long GetNAccepted(int bin) { return InfoOf(bin).nAccepted(); }
};

#endif  // PYTHIAGUN_H
