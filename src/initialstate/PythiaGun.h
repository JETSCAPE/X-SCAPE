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
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

//! Set when PythiaGun understands <Hard><PythiaGun><pTHatBins> (several pTHat
//! windows, one Pythia instance each), so code built against it can check.
#define PYTHIAGUN_HAS_PTHAT_BINS 1
//! Set when PythiaGun understands <partonYMax> / <partonYMode> (a rapidity cut on
//! the partons it hands to the framework).
#define PYTHIAGUN_HAS_PARTON_Y_CUT 1
//! Set when PythiaGun hands over final-state leptons next to the partons and
//! understands <minPartons> and <trigger>Z (with <triggerYMax>, <leptonPtMin>,
//! <leptonEtaMax>: the Z acceptance on its decay leptons).
#define PYTHIAGUN_HAS_Z_TRIGGER 1
//! Set when PythiaGun understands <extraHard>: additional hard scatterings at
//! other binary collisions of the event (and <extraHard><noTrigger>).
#define PYTHIAGUN_HAS_EXTRA_HARD 1

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

    Optional: a rapidity cut on the partons this gun hands to the framework
    (status 62 after ISR/MPI, or the final partons with FSR_on), <partonYMax>
    (e.g. 0.6) and <partonYMode> on their two hardest (by pT): "leading" the
    hardest has |y| < partonYMax, "both" the two hardest, "any" either of them.
    (At the process level the two hard partons have equal pT, so "leading" is
    only defined after ISR; with a single parton all three modes cut on it.)
    A rejected event is regenerated, like one with fewer than <minPartons>
    partons.  Pythia's sigmaGen still counts rejected events, so per window the
    events reaching the cut and those passing it are counted, and
    GetSigmaGen()/GetSigmaErr() return sigmaGen x kept/tried (GetSigmaGenRaw is
    Pythia's own).  An accepted event keeps all its partons.  Without
    <partonYMax> nothing is cut and everything is as before.

    Leptons: final-state electrons and muons (a Z's decay products, status 23
    after the resonance decay, not 62) are handed over next to the partons, and
    the energy loss and the liquefier pass them through like photons
    (Parton::isNonParton).  They are not partons: <softMomentumCutoff>,
    <minPartons> and <partonYMax> ignore them.  <minPartons> (default 2) is the
    number of quarks, gluons or photons an event needs; with fewer it is
    regenerated.

    Optional: <trigger>Z</trigger> (default none).  The event needs an
    opposite-charge, same-flavour lepton pair (the one with the largest
    invariant mass if there are several), the Z, which passes <triggerYMax>
    (|y| of the pair), <leptonPtMin> and <leptonEtaMax> (each lepton); each cut
    is off when empty.  A rejected event is regenerated and counted with
    <partonYMax>: an event is kept if it passes both, and kept/tried is the
    acceptance of both.  GetTriggerLeptons() gives the current event's pair.
    After 100000 rejected tries in one event PythiaGun throws (a process that
    never passes, e.g. <trigger>Z without a Z in <LinesToRead>).

    Optional: <extraHard><on>1</on> adds the event's other hard scatterings
    (js-contrib docs/Plans/PLAN_multi_hard.md).  After the triggering
    scattering (scattering 0, everything above unchanged) come N ~
    Poisson(N_coll sigma/sigma_NN) more, N_coll from the initial state, sigma
    the cross section above <pTHatMin> of a Pythia instance of their own
    (HardQCD:all, ISR on, MPI off, FSR as <FSR_on>, + <extraHard><LinesToRead>)
    from <warmup> events at Init, sigma_NN <sigmaNN> (empty: 3dMCGlauber's
    parametrization at <eCM>); <fixedN> >= 0 fixes N.  Each sits at its own
    binary collision, drawn uniformly (InitialState::
    SampleABinaryCollisionPointWith), and hands over its partons as the
    trigger does, after the trigger's, with plabel = its index (1, 2, ...) and
    colour tags offset by 100000 x index.  Their Pythia, the Poisson draw and
    the vertices use random streams of their own, and JetEnergyLossManager
    showers their partons with engines of their own (GetHardPartonScatter), so
    the trigger's partons and showers are the same with and without them.
    <noTrigger>1: no triggering scattering at all, only the extras (N may be
    0); the windows, <partonYMax> and <trigger> are then refused.
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

  double partonYMax_ = 0.;  //!< 0: no cut
  std::string partonYMode_ = "leading";
  //! per window: events tested / passing <partonYMax> and the trigger
  std::vector<long> yTried_, yKept_;

  int minPartons_ = 2;  //!< quarks, gluons, photons an event needs
  std::string trigger_ = "none";  //!< none | Z
  //! the Z acceptance (0: no cut)
  double triggerYMax_ = 0., leptonPtMin_ = 0., leptonEtaMax_ = 0.;
  //! this event's trigger leptons, hardest first (two with trigger Z)
  std::vector<Pythia8::Particle> triggerLeptons_;

 public:
  //! one hard scattering of the current event (0: the trigger)
  struct ScatterRecord {
    int scatter;     //!< 0 the trigger, 1, 2, ... the extras
    double pthat;    //!< Pythia's pTHat
    int code;        //!< Pythia's process code
    double x, y;     //!< the vertex [fm]
    int n_partons;   //!< partons (and leptons) handed over
    double e, pt;    //!< their summed energy and scalar pT [GeV]
  };

 private:
  // <extraHard>
  bool extraOn_ = false;
  bool extraNoTrigger_ = false;
  double extraPtHatMin_ = 15.;
  double extraSigmaNN_ = 0.;  //!< [mb]
  long extraWarmup_ = 20000;
  int extraFixedN_ = -1;      //!< < 0: Poisson
  unsigned int extraPythiaSeed_ = 0, extraRngSeed_ = 0;
  double extraSigma_ = 0., extraSigmaErr_ = 0.;  //!< from the warm-up [mb]
  std::unique_ptr<Pythia8::Pythia> extraHardPythia_;
  std::mt19937 extraRng_;  //!< the Poisson draw and the vertices
  // the current event
  int nExtra_ = 0;
  double extraMu_ = 0.;
  std::vector<int> partonScatter_;      //!< per hard parton
  std::vector<ScatterRecord> scatters_;

  void ReadExtraHard(unsigned int seed);
  void ConfigureExtraPythia();
  //! the particles of ``py``'s event to hand over (status 62 partons above
  //! softMomentumCutoff, leptons, or the final partons with FSR_on)
  void SelectPartons(Pythia8::Pythia &py, std::vector<Pythia8::Particle> &out);
  //! hands ``sorted`` over at (x, y) as scattering ``scatter`` and records it
  void AddScatterPartons(const std::vector<Pythia8::Particle> &sorted,
                         double x, double y, int scatter, double pthat,
                         int code);
  void AddExtraScatterings();

  void ReadPtHatBins();
  void ReadPartonYCut();
  void ReadTrigger();
  //! whether events are counted and cut (<partonYMax> or a trigger)
  bool CutsOn() const { return partonYMax_ > 0 || trigger_ == "Z"; }
  //! the rapidity cut on the pT-sorted particles to hand over (at least one
  //! parton; leptons are skipped)
  bool PassesPartonYCut(const std::vector<Pythia8::Particle> &sorted) const;
  //! the trigger's lepton pair among the particles to hand over ->
  //! triggerLeptons_; false if there is none
  bool FindTriggerPair(const std::vector<Pythia8::Particle> &sorted);
  //! triggerLeptons_ against <triggerYMax>, <leptonPtMin>, <leptonEtaMax>
  bool PassesTriggerCut() const;
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
  // window (with <partonYMax> or a trigger: of the events passing the cuts).
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
  //! acceptance of <partonYMax> and the trigger if set (then the error
  //! includes the acceptance's binomial error)
  double GetSigmaGen(int bin) {
    const double s = InfoOf(bin).sigmaGen();
    return CutsOn() ? s * GetYAcceptance(bin) : s;
  }
  double GetSigmaErr(int bin);
  //! Pythia's own estimate, including the events the cut rejected
  double GetSigmaGenRaw(int bin) { return InfoOf(bin).sigmaGen(); }
  double GetSigmaErrRaw(int bin) { return InfoOf(bin).sigmaErr(); }
  //! events Pythia accepted in one window so far (incl. those this gun
  //! rejected and regenerated, which also enter sigmaGen)
  long GetNAccepted(int bin) { return InfoOf(bin).nAccepted(); }

  // Rapidity cut on the handed-over partons (<partonYMax>; 0 = none).
  double GetPartonYMax() const { return partonYMax_; }
  std::string GetPartonYMode() const { return partonYMode_; }
  //! events that reached / passed the cuts (<partonYMax> and the trigger) in
  //! one window (0 without any)
  long GetNYTried(int bin) const {
    return yTried_.empty() ? 0 : yTried_.at(bin);
  }
  long GetNYKept(int bin) const { return yKept_.empty() ? 0 : yKept_.at(bin); }
  //! kept / tried (1 without a cut or before the first event)
  double GetYAcceptance(int bin) const {
    const long n = GetNYTried(bin);
    return n > 0 ? double(GetNYKept(bin)) / n : 1.0;
  }

  // Leptons and the trigger (<minPartons>, <trigger>; 0 = no cut).
  int GetMinPartons() const { return minPartons_; }
  std::string GetTrigger() const { return trigger_; }
  double GetTriggerYMax() const { return triggerYMax_; }
  double GetLeptonPtMin() const { return leptonPtMin_; }
  double GetLeptonEtaMax() const { return leptonEtaMax_; }
  //! the current event's trigger leptons, hardest first (two with <trigger>Z,
  //! else none)
  const std::vector<Pythia8::Particle> &GetTriggerLeptons() const {
    return triggerLeptons_;
  }

  // <extraHard>: the settings, the warm-up cross section and this event's
  // scatterings.
  bool GetExtraHardOn() const { return extraOn_; }
  bool GetExtraNoTrigger() const { return extraNoTrigger_; }
  double GetExtraPtHatMin() const { return extraPtHatMin_; }
  double GetExtraSigma() const { return extraSigma_; }
  double GetExtraSigmaErr() const { return extraSigmaErr_; }
  double GetExtraSigmaNN() const { return extraSigmaNN_; }
  long GetExtraWarmup() const { return extraWarmup_; }
  int GetExtraFixedN() const { return extraFixedN_; }
  unsigned int GetExtraPythiaSeed() const { return extraPythiaSeed_; }
  unsigned int GetExtraRngSeed() const { return extraRngSeed_; }
  //! the current event's number of extras and its Poisson mean
  int GetNExtra() const { return nExtra_; }
  double GetExtraMu() const { return extraMu_; }
  //! the current event's scatterings, the trigger's first (none with
  //! <noTrigger>)
  const std::vector<ScatterRecord> &GetScatterRecords() const {
    return scatters_;
  }
  int GetHardPartonScatter(int i) const override {
    return i >= 0 && i < static_cast<int>(partonScatter_.size())
               ? partonScatter_[i]
               : 0;
  }
};

#endif  // PYTHIAGUN_H
