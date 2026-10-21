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

#include <map>
#include <string>
#include <vector>

// Capability, for code built against several X-SCAPE versions (PyJetscape):
// per-event reseeding (<reseed_per_event>), the recorded Pythia input
// (<store_input>, GetLastInput) and HadronizeInput.
#define ICOLORED_HAS_INPUT_RECORD 1
// ... and the colour audit of the last event (GetLastColourAudit).
#define ICOLORED_HAS_COLOUR_AUDIT 1
// ... and <colour_closure>, HadronizePartons (ShowerParton lists).
#define ICOLORED_HAS_COLOUR_CLOSURE 1

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

  /// A final parton as a shower hands it over (HadronizePartons).
  struct ShowerParton {
    int shower;  ///< the shower it belongs to; -1: a beam remnant (always
                 ///< handed to Pythia, whatever its pstat)
    int id, pstat, col, acol;
    double px, py, pz, e;
  };

  /**
   * Hadronize final partons (and, with shower -1, beam remnants) the way
   * DoHadronization does once it has collected them: holes and the
   * liquefier's partons left out, the duplicate check, the closure of
   * <colour_closure>, Pythia. For replays of stored partons and for tests.
   */
  void HadronizePartons(const std::vector<ShowerParton> &partons,
                        vector<shared_ptr<Hadron>> &hOut);

  /// What the last hadronization found in the colour tags of its input and
  /// how it closed the open lines (PLAN_lbt_colour_tags.md, C3).
  struct ColourAudit {
    int partons = 0;          ///< final shower partons handed to Pythia
    int remnants = 0;         ///< beam remnants handed to Pythia
    int holes = 0;            ///< pstat -1, not handed to Pythia
    double hole_energy = 0.;  ///< their energy [GeV], dropped
    int skipped = 0;          ///< other pstat < 0 (the liquefier's), skipped
    int unpaired_col = 0;     ///< open colour tags, before the closure
    int unpaired_acol = 0;    ///< open anticolour tags, before the closure
    int duplicate_tags = 0;   ///< extra uses of a tag as colour or anticolour
    int relabelled_tags = 0;  ///< duplicates relabelled (<colour_duplicates>)
    int reconnected = 0;      ///< open ends joined (B2), pairs
    int reconnected_cross = 0;  ///< joined across showers as a fallback
    int partners = 0;           ///< medium partners inserted for holes (B1)
    double partner_energy = 0.;  ///< their energy [GeV]
    int fakes = 0;               ///< soft end partons along the beam (B3)
    double fake_energy = 0.;     ///< their energy [GeV]
    int retries = 0;  ///< fragmentation failed, redone with soft end partons
  };
  const ColourAudit &GetLastColourAudit() const { return audit_; }

  /// The closure of open colour lines, as <colour_closure> (reconnect,
  /// reconnect_global, thermal, fake, none); InitTask sets it from the XML.
  void SetColourClosure(const std::string &mode);
  std::string GetColourClosure() const;

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
  // <JetHadronization><colour_closure>: reconnect (default), reconnect_global,
  // thermal, fake, none (PLAN_lbt_colour_tags.md, B)
  enum class Closure { Reconnect, ReconnectGlobal, Thermal, Fake, None };
  Closure closure_ = Closure::Reconnect;

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
  // everything after collecting the partons: pstat policy, duplicates,
  // closure, Pythia, the hadrons
  void HadronizeList(const std::vector<ShowerParton> &partons,
                     vector<shared_ptr<Hadron>> &hOut);
  // with <store_input> 1: keep the record as it goes to Pythia
  void RecordInput();
  // The record entries of the open colour (anticolour) tags, by tag.
  static void OpenEnds(const Pythia8::Event &event, std::map<int, int> &col_at,
                       std::map<int, int> &acol_at);
  // B2: join open colour and anticolour ends (within a shower, or across);
  // returns the number of pairs joined
  static int Reconnect(Pythia8::Event &event, const std::vector<int> &shower_of,
                       bool within_shower);
  // B1: the holes' medium partners for the lines still open
  void InsertPartners(Pythia8::Event &event, std::vector<int> &shower_of,
                      const std::vector<ShowerParton> &holes);
  // B3: soft end partons along the beam for every end still open
  void InsertFakes(Pythia8::Event &event, std::vector<int> &shower_of);
  void CloseColourLines(Pythia8::Event &event, std::vector<int> &shower_of,
                        const std::vector<ShowerParton> &holes, Closure mode);
  // the final hadrons of pythia.event after pythia.next()
  void CollectHadrons(vector<shared_ptr<Hadron>> &hOut);

  // Allows the registration of the module so that it is available to be used by
  // the Jetscape framework.
  static RegisterJetScapeModule<iColoredHadronization> reg;

 protected:
  static Pythia8::Pythia pythia;
};

#endif  // iColoredHadronization_H
