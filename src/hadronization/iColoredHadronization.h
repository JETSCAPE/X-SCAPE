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

#ifndef iColoredHadronization_H
#define iColoredHadronization_H

#include "HadronizationModule.h"
#include "Pythia8/Pythia.h"

// Capability, for code built against several X-SCAPE versions (PyJetscape):
// per-event reseeding (<reseed_per_event>), the recorded Pythia input
// (<store_input>, GetLastInput) and HadronizeInput.
#define ICOLORED_HAS_INPUT_RECORD 1
// ... and the colour audit of the last event (GetLastColourAudit).
#define ICOLORED_HAS_COLOUR_AUDIT 1

using namespace Jetscape;

class iColoredHadronization
    : public HadronizationModule<iColoredHadronization> {
 public:
  iColoredHadronization();
  virtual ~iColoredHadronization();

  void InitTask();
  void DoHadronization(vector<vector<shared_ptr<Parton>>> &shower,
                       vector<shared_ptr<Hadron>> &hOut,
                       vector<shared_ptr<Parton>> &pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

  /// One entry of the parton list handed to Pythia (status 23): the final
  /// partons, the beam remnants and the colour-balancing parton, as appended.
  struct InputParton {
    int id, col, acol;
    double px, py, pz, e, m;
  };

  /**
   * Use this seed for the next DoHadronization/HadronizeInput (one-shot): the
   * Pythia generator is reseeded with it first.
   */
  void SetNextRandomSeed(unsigned int seed) {
    next_seed_ = seed;
    has_next_seed_ = true;
  }

  /// The seed the last hadronization used; 0 if it was not reseeded.
  unsigned int GetLastRandomSeed() const { return last_seed_; }

  /// With <store_input> 1: the last DoHadronization's Pythia input (else empty).
  const std::vector<InputParton> &GetLastInput() const { return last_input_; }

  /**
   * Hadronize a stored Pythia input (GetLastInput of an earlier run) the way
   * DoHadronization does after building it: no beam remnants or hydro needed.
   * With the seed of that run (SetNextRandomSeed) and the same Pythia version
   * and settings, the hadrons are the same.
   */
  void HadronizeInput(const std::vector<InputParton> &input,
                      vector<shared_ptr<Hadron>> &hOut);

  /// What the last DoHadronization found in the colour tags of its input.
  struct ColourAudit {
    int partons = 0;          ///< final shower partons handed to Pythia
    int remnants = 0;         ///< beam remnants handed to Pythia
    int holes = 0;            ///< pstat -1, not handed to Pythia
    double hole_energy = 0.;  ///< their energy [GeV], dropped
    int skipped = 0;          ///< other pstat < 0 (the liquefier's), skipped
    int unpaired_col = 0;     ///< colour tags without an anticolour partner
    int unpaired_acol = 0;    ///< anticolour tags without a colour partner
    int duplicate_tags = 0;   ///< extra uses of a tag as colour or anticolour
    int relabelled_tags = 0;  ///< duplicates relabelled (<colour_duplicates>)
  };
  const ColourAudit &GetLastColourAudit() const { return audit_; }

 private:
  double p_fake;

  // <JetHadronization><reseed_per_event> 1: every event draws a seed from the
  // module's generator and reseeds Pythia with it, so each event can be
  // reproduced on its own (GetLastRandomSeed). 0 (default): Pythia runs on
  // unbroken from its default seed.
  bool reseed_per_event_ = false;
  bool has_next_seed_ = false;
  unsigned int next_seed_ = 0;
  unsigned int last_seed_ = 0;
  // <JetHadronization><store_input> 1: keep the Pythia input (GetLastInput)
  bool store_input_ = false;
  std::vector<InputParton> last_input_;
  // <JetHadronization><colour_duplicates> relabel: give duplicated tags new
  // ones (a counted repair) instead of failing the event (fail, default)
  bool relabel_duplicates_ = false;
  ColourAudit audit_;

  void Reseed(unsigned int seed);
  void ReseedForEvent();
  // Count the colour tags of the record: returns the number of extra uses of
  // a tag (duplicates) and lists the unpaired tags, one entry per missing
  // partner.
  static int CountColourTags(const Pythia8::Event &event,
                             std::vector<int> &cols, std::vector<int> &acols);
  // Give every second, third, ... use of a tag a new tag; returns the number
  // of colour fields changed.
  static int RelabelDuplicateTags(Pythia8::Event &event);
  // the final hadrons of pythia.event after pythia.next()
  void CollectHadrons(vector<shared_ptr<Hadron>> &hOut);

  // Allows the registration of the module so that it is available to be used by
  // the Jetscape framework.
  static RegisterJetScapeModule<iColoredHadronization> reg;

 protected:
  static Pythia8::Pythia pythia;
};

#endif  // iColoredHadronization_H
