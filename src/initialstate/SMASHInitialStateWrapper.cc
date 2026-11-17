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
// -----------------------------------------------------------------------------
// This is a wrapper for SMASH hadronic transport with the JETSCAPE framework
// for the collider modus of SMASH to generate initial conditions
// -----------------------------------------------------------------------------

#include <chrono>
#include "SMASHInitialStateWrapper.h"
#include "Transport.h"
#include "SmashTabulations.h"

#include "smash/particles.h"
#include "smash/library.h"
#include "smash/particledata.h"
#include "smash/processbranch.h"

#include <math.h>
#include <string>
#include <map>
#include <algorithm>
#include <filesystem>

#include <boost/lexical_cast.hpp>

// Provide to_string overload for smash::PdgCode so SMASH's YAML conversion can
// encode keys
namespace smash {
inline std::string to_string(const PdgCode& code) {
  // Use the decimal representation as string, e.g., "2212"
  return std::to_string(code.get_decimal());
}
}  // namespace smash

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<SmashInitialConditionWrapper>
    SmashInitialConditionWrapper::reg("SMASHInitialState");

SmashInitialConditionWrapper::SmashInitialConditionWrapper() {
  SetId("SMASHInitialState");
}

void SmashInitialConditionWrapper::InitTask() {
  JSINFO << "SMASH: picking SMASH-specific configuration from xml file";
  std::string smash_config =
      GetXMLElementText({"IS", "SMASH", "SMASH_config_file"});
  std::string smash_hadron_list =
      GetXMLElementText({"IS", "SMASH", "SMASH_particles_file"});
  std::string smash_decays_list =
      GetXMLElementText({"IS", "SMASH", "SMASH_decaymodes_file"});
  // output path is just dummy here, because no output from SMASH is foreseen
  std::filesystem::path output_path("./");
  const std::string smash_version(SMASH_VERSION);

  auto config = smash::setup_config_and_logging(smash_config, smash_hadron_list,
                                                smash_decays_list);

  // Dynamic fluidization: with a Cartesian hydro next to SMASH, hadrons in
  // dense regions become hydro sources during the run (DynamicFluidization;
  // defaults after arXiv:2507.19389)
  {
    const std::string cartesian =
        GetXMLElementText({"Hydro", "CartesianHydro"}, false);
    const bool is_cartesian = cartesian.find("true") != std::string::npos;
    const int dyn = GetXMLElementInt({"BDM", "dynamic_fluidization"}, false);
    dynamic_fluidization_ = is_cartesian && dyn == 1;
    if (dyn == 1 && !is_cartesian) {
      JSWARN << "<BDM><dynamic_fluidization> needs <Hydro><CartesianHydro> "
                "true: off";
    }
    if (dynamic_fluidization_) {
      const double sigma =
          GetXMLElementDouble({"BDM", "fluidization_sigma"}, false);
      const double e_th =
          GetXMLElementDouble({"BDM", "fluidization_e_threshold"}, false);
      fluidization_string_time_ =
          GetXMLElementDouble({"BDM", "fluidization_string_time"}, false);
      fluidization_ = std::make_unique<DynamicFluidization>(sigma, e_th);
      JSINFO << "SMASH: dynamic fluidization, e_th = " << e_th
             << " GeV/fm^3, sigma = " << sigma << " fm, t_f = "
             << fluidization_string_time_ << " fm";
    }
  }

  // Take care of the random seed. This will make SMASH results reproducible.
  auto random_seed = (*GetMt19937Generator())();
  config.set_value(smash::InputKeys::gen_randomseed, random_seed);
  // Read in the rest of configuration
  if (IsTimeStepped()) {
    end_time_ = GetMainClock()->GetEndTime();
  } else {
    end_time_ = GetXMLElementDouble({"IS", "SMASH", "End_Time"});
  }
  config.set_value(smash::InputKeys::gen_endTime, end_time_);
  JSINFO << "End time until which SMASH initial condition propagates is "
         << end_time_ << " fm/c";

  int smash_projectile_protons =
      GetXMLElementInt({"IS", "SMASH", "ProjectileProtons"});
  int smash_projectile_neutrons =
      GetXMLElementInt({"IS", "SMASH", "ProjectileNeutrons"});
  int smash_target_protons = GetXMLElementInt({"IS", "SMASH", "TargetProtons"});
  int smash_target_neutrons =
      GetXMLElementInt({"IS", "SMASH", "TargetNeutrons"});

  auto to_pdg_map = [](const std::map<int, int>& input) {
    std::map<smash::PdgCode, int> output;
    for (const auto& [pdg_int, count] : input) {
      output.emplace(smash::PdgCode::from_decimal(pdg_int), count);
    }
    return output;
  };

  std::map<int, int> smash_projectile{{2212, smash_projectile_protons},
                                      {2112, smash_projectile_neutrons}};
  config.set_value(smash::InputKeys::modi_collider_projectile_particles,
                   to_pdg_map(smash_projectile));

  std::map<int, int> smash_target{{2212, smash_target_protons},
                                  {2112, smash_target_neutrons}};
  config.set_value(smash::InputKeys::modi_collider_target_particles,
                   to_pdg_map(smash_target));

  double smash_sqrtsnn = GetXMLElementDouble({"IS", "SMASH", "Sqrtsnn"});
  config.set_value(smash::InputKeys::modi_collider_sqrtSNN, smash_sqrtsnn);

  std::string smash_FermiMotion =
      GetXMLElementText({"IS", "SMASH", "FermiMotion"});
  smash::FermiMotion fm = ParseFermiMotion(smash_FermiMotion);
  config.set_value(smash::InputKeys::modi_collider_fermiMotion, fm);

  bool smash_CollisionsWithinNucleus =
      GetXMLElementInt({"IS", "SMASH", "CollisionsWithinNucleus"});
  config.set_value(smash::InputKeys::modi_collider_collisionWithinNucleus,
                   smash_CollisionsWithinNucleus);

  bool smash_impact_react_plane =
      GetXMLElementInt({"IS", "SMASH", "Impact", "RandomReactionPlane"});
  config.set_value(smash::InputKeys::modi_collider_impact_randomReactionPlane,
                   smash_impact_react_plane);

  int smash_impact_param_mode =
      GetXMLElementInt({"IS", "SMASH", "Impact", "Mode"});
  const double smash_impact_val =
      GetXMLElementDouble({"IS", "SMASH", "Impact", "Value"});
  std::string smash_impact_sample =
      GetXMLElementText({"IS", "SMASH", "Impact", "Sample"});
  const double smash_impact_valMin =
      GetXMLElementDouble({"IS", "SMASH", "Impact", "ImpactMin"});
  const double smash_impact_valMax =
      GetXMLElementDouble({"IS", "SMASH", "Impact", "ImpactMax"});

  if (smash_impact_param_mode == 0) {
    config.set_value(smash::InputKeys::modi_collider_impact_value,
                     smash_impact_val);
  } else if (smash_impact_param_mode == 1) {
    smash::Sampling sample = ParseImpactParameterSampling(smash_impact_sample);
    config.set_value(smash::InputKeys::modi_collider_impact_sample, sample);
    const std::array<double, 2> smash_impact_range = {smash_impact_valMin,
                                                      smash_impact_valMax};
    config.set_value(smash::InputKeys::modi_collider_impact_range,
                     smash_impact_range);
  } else {
    JSWARN << "This SMASH impact parameter mode does not exist in Jetscape";
    exit(1);
  }

  // particles, decays and tabulations, once per process (SmashTabulations.h)
  InitializeSmashOnce(config, smash_version);

  const double delta_t_sm = GetXMLElementDouble({"IS", "SMASH", "Delta_Time"});
  config.set_value(smash::InputKeys::gen_deltaTime, delta_t_sm);

  // Enforce timestep compatibility (temporarily)
  if (IsTimeStepped()) {
    const double delta_t_js = GetMainClock()->GetDeltaT();
    const double ts_rem = std::remainder(delta_t_js, delta_t_sm);
    const double ts_frac = delta_t_js / delta_t_sm;
    if (!(ts_rem < 1E-6 && ts_frac > 1.0)) {
      JSWARN << "Timesteps of SMASH (dt = " << delta_t_sm
             << ") and JETSCAPE (dt = " << delta_t_js
             << ") are incompatible."
                "SMASH IC timesteps should be a half, a third, etc. from "
                "JETSCAPE's";
    }
  }
  smash_collider_experiment_ =
      make_shared<smash::Experiment<smash::ColliderModus>>(config, output_path);
  config.clear();
  JSINFO << "Finish initializing SMASH initial condition";
}

smash::FermiMotion SmashInitialConditionWrapper::ParseFermiMotion(
    const std::string& s) {
  std::string lower = s;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

  if (lower == "off")
    return smash::FermiMotion::Off;
  if (lower == "on")
    return smash::FermiMotion::On;
  if (lower == "frozen")
    return smash::FermiMotion::Frozen;

  throw std::runtime_error("Invalid FermiMotion string: " + s);
}

smash::Sampling SmashInitialConditionWrapper::ParseImpactParameterSampling(
    const std::string& s) {
  std::string lower = s;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

  if (lower == "uniform")
    return smash::Sampling::Uniform;
  if (lower == "quadratic")
    return smash::Sampling::Quadratic;
  if (lower == "custom")
    return smash::Sampling::Custom;

  throw std::runtime_error("Invalid Sampling string: " + s);
}

void SmashInitialConditionWrapper::ExecuteTask() {
  VERBOSE(2) << "SMASH initial condition running: " << GetId() << "...";
  InitPerEvent();
  CalculateTimeTask();
  FinishPerEvent();
}

void SmashInitialConditionWrapper::InitPerEvent() {
  VERBOSE(3) << "Initializing SMASH initial condition event...";
  fluidize_ids_.clear();
  from_fluid_.clear();
  previous_ids_.clear();
  smash_collider_experiment_->initialize_new_event();
}

void SmashInitialConditionWrapper::CalculateTimeTask() {
  std::vector<shared_ptr<Hadron>> hadrons_to_add =
      Transport::GetTimestepParticlizationHadrons();
  std::vector<shared_ptr<Hadron>> hadrons_to_remove =
      Transport::GetTimestepHadronsToRemove();

  // Dynamic fluidization: the fluidized hadrons (status 29, label = SMASH
  // id) are removed by id, not by matching position and momentum
  std::unordered_set<int32_t> remove_ids;
  hadrons_to_remove.erase(
      std::remove_if(hadrons_to_remove.begin(), hadrons_to_remove.end(),
                     [&](const std::shared_ptr<Jetscape::Hadron>& h) {
                       if (h->pstat() != 29) return false;
                       remove_ids.insert(h->plabel());
                       return true;
                     }),
      hadrons_to_remove.end());

  VERBOSE(3) << "SMASH initial condition got " << hadrons_to_add.size()
             << " hadrons in this timestep.";
  VERBOSE(3) << "SMASH initial condition removed " << hadrons_to_remove.size()
             << " hadrons in this timestep.";

  // Check if absolute value of position and momentum vectors can be computed
  // without error in sqrt, otherwise remove hadrons from hadrons_to_add and
  // hadrons_to_remove, compute the absolute value of each jetscape hadron
  hadrons_to_add.erase(
      std::remove_if(hadrons_to_add.begin(), hadrons_to_add.end(),
                     [](const std::shared_ptr<Jetscape::Hadron>& h) {
                       return (!h->has_valid_momentum() || h->pid() == 22);
                     }),
      hadrons_to_add.end());
  hadrons_to_remove.erase(
      std::remove_if(hadrons_to_remove.begin(), hadrons_to_remove.end(),
                     [](const std::shared_ptr<Jetscape::Hadron>& h) {
                       return (!h->has_valid_momentum() || h->pid() == 22);
                     }),
      hadrons_to_remove.end());

  const double until_time =
      IsTimeStepped() ? GetMainClock()->GetCurrentTime() : end_time_;
  VERBOSE(3) << "Propagating SMASH IC until t = " << until_time;
  VERBOSE(3) << "End time until which SMASH initial condition propagates is "
             << end_time_ << " fm/c";

  smash::ParticleList add_list =
      get_smash_plist_from_JS_hadrons(hadrons_to_add);
  smash::ParticleList remove_list =
      find_smash_hadrons_and_get_exact_hadron_list(hadrons_to_remove);
  if (!remove_ids.empty()) {
    for (const auto& particle :
         *smash_collider_experiment_->first_ensemble()) {
      if (remove_ids.count(particle.id())) remove_list.push_back(particle);
    }
    VERBOSE(2) << "SMASH IC: " << remove_ids.size()
               << " fluidized hadrons to remove by id";
  }
  const auto prof_t0 = std::chrono::steady_clock::now();
  smash_collider_experiment_->run_time_evolution(
      until_time, std::move(add_list), std::move(remove_list));
  const auto prof_t1 = std::chrono::steady_clock::now();
  prof_smash_s_ += std::chrono::duration<double>(prof_t1 - prof_t0).count();
  if (dynamic_fluidization_) {
    DetermineFluidizingHadrons(until_time);
    prof_fluid_s_ += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - prof_t1).count();
  }
}

void SmashInitialConditionWrapper::UpdateFromFluidTags(double t) {
  const smash::Particles *particles =
      smash_collider_experiment_->first_ensemble();
  // the tagged particles that disappeared in this step (decayed, collided
  // inelastically, or were removed)
  std::unordered_set<int32_t> present;
  for (const auto &p : *particles) present.insert(p.id());
  std::vector<TaggedState> gone;
  for (const auto &kv : from_fluid_) {
    if (!present.count(kv.first)) gone.push_back(kv.second);
  }
  std::unordered_map<int32_t, TaggedState> tagged;
  const double match2 = from_fluid_match_ * from_fluid_match_;
  for (const auto &p : *particles) {
    const auto history = p.get_history();
    const smash::FourVector x = p.position(), mom = p.momentum();
    const double vx = mom.x1() / mom.x0(), vy = mom.x2() / mom.x0(),
                 vz = mom.x3() / mom.x0();
    bool tag = from_fluid_.count(p.id()) > 0 ||
               history.collisions_per_particle >= kFromFluid;
    if (!tag && !gone.empty() && !previous_ids_.count(p.id())) {
      // a new particle: its production vertex against the vanished tagged
      // particles, both at the time of its last collision
      const double tv = history.time_last_collision;
      const double xv = x.x1() - vx * (x.x0() - tv);
      const double yv = x.x2() - vy * (x.x0() - tv);
      const double zv = x.x3() - vz * (x.x0() - tv);
      for (const auto &g : gone) {
        const double dx = g.x + g.vx * (tv - g.t) - xv;
        const double dy = g.y + g.vy * (tv - g.t) - yv;
        const double dz = g.z + g.vz * (tv - g.t) - zv;
        if (dx * dx + dy * dy + dz * dz < match2) {
          tag = true;
          break;
        }
      }
    }
    if (tag) tagged[p.id()] = {x.x0(), x.x1(), x.x2(), x.x3(), vx, vy, vz};
  }
  from_fluid_.swap(tagged);
  previous_ids_.swap(present);
  (void)t;
}

void SmashInitialConditionWrapper::DetermineFluidizingHadrons(double t) {
  // The criterion of arXiv:2507.19389 (DynamicFluidization), with its
  // eligibility rules from the SMASH history: not unscattered nucleons (no
  // collision, or one elastic one), and string fragments only a time t_f
  // after their production in their rest frame.  Every particle counts in
  // the energy-momentum tensor.
  fluidize_ids_.clear();
  const auto prof_t0 = std::chrono::steady_clock::now();
  UpdateFromFluidTags(t);
  prof_tags_s_ += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - prof_t0).count();
  const smash::Particles *particles =
      smash_collider_experiment_->first_ensemble();
  std::vector<DynamicFluidization::Candidate> cand;
  std::vector<int32_t> ids;
  cand.reserve(particles->size());
  ids.reserve(particles->size());
  for (const auto &p : *particles) {
    const smash::FourVector mom = p.momentum(), pos = p.position();
    const double m = mom.abs();
    const auto history = p.get_history();
    bool eligible = p.type().is_hadron() && m > 0.;
    if (history.collisions_per_particle == 0 ||
        (history.collisions_per_particle == 1 &&
         history.process_type == smash::ProcessType::Elastic)) {
      eligible = false;
    }
    if (from_fluid_.count(p.id())) eligible = false;   // from the fluid
    if (eligible && (smash::is_string_soft_process(history.process_type) ||
                     history.process_type == smash::ProcessType::StringHard)) {
      const double gamma = mom.x0() / m;
      eligible = ((t - history.time_last_collision) / gamma >=
                  fluidization_string_time_);
    }
    cand.push_back({pos.x0(), pos.x1(), pos.x2(), pos.x3(), mom.x0(), mom.x1(),
                    mom.x2(), mom.x3(), m, eligible});
    ids.push_back(p.id());
  }
  const std::vector<bool> flags = fluidization_->Evaluate(cand);
  for (std::size_t i = 0; i < flags.size(); i++) {
    if (flags[i]) fluidize_ids_.insert(ids[i]);
  }
  VERBOSE(2) << "SMASH IC at t = " << t << ": " << fluidize_ids_.size()
             << " of " << cand.size() << " hadrons fluidize ("
             << from_fluid_.size() << " from the fluid)";
}

void SmashInitialConditionWrapper::FinishPerEvent() {
  JSINFO << "Finishing SMASH initial condition event...";
  if (dynamic_fluidization_) {
    JSINFO << "SMASH IC, dynamic fluidization: SMASH evolution "
           << prof_smash_s_ << " s; fluidization criterion " << prof_fluid_s_
           << " s (of which the from-fluid tags " << prof_tags_s_ << " s)";
  }
  prof_smash_s_ = prof_fluid_s_ = prof_tags_s_ = 0.;
  event_number_++;

  // SMASH within JETSCAPE only works with one (the first) ensemble
  smash::Particles* smash_particles =
      smash_collider_experiment_->first_ensemble();
  int ev_no = current_event_number();

  smash_collider_experiment_->do_final_decays();
  smash_collider_experiment_->final_output();

  if (ev_no > jetscape_hadrons_.size()) {
    jetscape_hadrons_.resize(ev_no);
  }

  fill_JS_hadrons_from_smash_particles(*smash_particles,
                                       jetscape_hadrons_[ev_no - 1]);

  smash_collider_experiment_
      ->increase_event_number();  // internal SMASH event counter
  JSINFO << jetscape_hadrons_[ev_no - 1].size()
         << " hadrons from SMASH initial condition.";
  JSINFO << "Finished SMASH collider event...";
}

void SmashInitialConditionWrapper::WriteTask(weak_ptr<JetScapeWriter> w) {
  JSINFO << "SMASH initial condition printout";
  auto f = w.lock();
  if (!f) {
    return;
  }
  f->WriteComment("JetScape module: " + GetId());
  for (const auto& event : jetscape_hadrons_) {
    int i = -1;
    for (const auto hadron : event) {
      f->WriteWhiteSpace("[" + to_string(++i) + "] H");
      f->Write(hadron);
    }
  }
}

std::vector<Hadron> SmashInitialConditionWrapper::GetCurrentHadronList() const {
  std::vector<Hadron> h_list;
  smash::Particles* smash_particles =
      smash_collider_experiment_->first_ensemble();

  for (const auto& particle : *smash_particles) {
    // 29: fluidizes now (dynamic fluidization; the label is then the SMASH
    // id, by which CalculateTimeTask removes it), else 28
    const bool fluidizes = fluidize_ids_.count(particle.id()) > 0;
    const int hadron_label = fluidizes ? particle.id() : 0;
    const int hadron_status = fluidizes ? 29 : 28;
    const int hadron_id = particle.pdgcode().get_decimal();
    smash::FourVector p = particle.momentum(), r = particle.position();
    const FourVector hadron_p(p.x1(), p.x2(), p.x3(), p.x0()),
        hadron_r(r.x1(), r.x2(), r.x3(), r.x0());
    const double hadron_mass = p.abs();
    const int charge = particle.type().charge();
    const int baryon_number = particle.type().baryon_number();
    const int strangeness = particle.type().strangeness();
    const auto history = particle.get_history();
    bool participant = false;
    if (history.collisions_per_particle > 0) {
      participant = true;
    }
    h_list.push_back(Hadron(hadron_label, hadron_id, hadron_status, hadron_p,
                            hadron_r, hadron_mass, charge, baryon_number,
                            strangeness, participant));
  }
  return h_list;
}

smash::ParticleList
SmashInitialConditionWrapper::get_smash_plist_from_JS_hadrons(
    const std::vector<shared_ptr<Hadron>>& JS_hadrons) {
  smash::ParticleList new_particles;
  for (const auto& JS_had : JS_hadrons) {
    const FourVector p = JS_had->p_in();
    const FourVector r = JS_had->x_in();
    smash::ParticleData new_p{
        smash::ParticleType::find(smash::PdgCode::from_decimal(JS_had->pid()))};
    new_p.set_4position(smash::FourVector(r.t(), r.x(), r.y(), r.z()));
    new_p.set_4momentum(smash::FourVector(p.t(), p.x(), p.y(), p.z()));
    if (dynamic_fluidization_) {
      // from the fluid: never fluidizes again (UpdateFromFluidTags)
      new_p.set_history(kFromFluid, 0, smash::ProcessType::None, r.t(), {});
    }
    new_particles.push_back(new_p);
  }
  return new_particles;
}

smash::ParticleList
SmashInitialConditionWrapper::find_smash_hadrons_and_get_exact_hadron_list(
    const std::vector<shared_ptr<Hadron>>& JS_hadrons) {
  // Get the current SMASH particle list
  smash::Particles* smash_particles =
      smash_collider_experiment_->first_ensemble();
  // New particle list to store the exact hadron properties
  // This will be filled with the exact hadron properties from SMASH
  // based on the Jetscape hadron list
  smash::ParticleList exact_hadron_list;

  // This function finds the hadrons in the SMASH particle list and returns a
  // new particle list with the exact hadron properties using the SMASH
  // particles.
  for (const auto& JS_had : JS_hadrons) {
    // Get the pdg code and momentum/position from the Jetscape hadron
    const int hadron_id = JS_had->pid();
    const FourVector p = JS_had->p_in();
    const FourVector r = JS_had->x_in();
    bool found = false;

    // Loop over the SMASH particles to find the matching hadron
    for (auto& smash_particle : *smash_particles) {
      // Check if the SMASH particle matches the Jetscape hadron
      if (smash_particle.pdgcode().get_decimal() != hadron_id) {
        continue;  // Skip if the pdg code does not match
      }

      // Check if the momentum and position match within
      // hadron_property_tolerance_
      if (std::abs(smash_particle.momentum().x0() - p.t()) >
              hadron_property_tolerance_ ||
          std::abs(smash_particle.momentum().x1() - p.x()) >
              hadron_property_tolerance_ ||
          std::abs(smash_particle.momentum().x2() - p.y()) >
              hadron_property_tolerance_ ||
          std::abs(smash_particle.momentum().x3() - p.z()) >
              hadron_property_tolerance_ ||
          std::abs(smash_particle.position().x0() - r.t()) >
              hadron_property_tolerance_ ||
          std::abs(smash_particle.position().x1() - r.x()) >
              hadron_property_tolerance_ ||
          std::abs(smash_particle.position().x2() - r.y()) >
              hadron_property_tolerance_ ||
          std::abs(smash_particle.position().x3() - r.z()) >
              hadron_property_tolerance_) {
        continue;  // Skip if the momentum or position does not match
      }
      // If we reach here, we have found a matching hadron, add this to the
      // exact hadron list
      double mass = smash_particle.effective_mass();
      auto is_particle_stable_and_with_invalid_mass =
          [&mass](const smash::ParticleData& p) {
            return p.type().is_stable() &&
                   std::abs(mass - p.pole_mass()) > rounding_error;
          };
      smash::ParticleData exact_hadron{smash::ParticleType::find(
          smash::PdgCode::from_decimal(JS_had->pid()))};
      smash::FourVector p = smash_particle.momentum();
      if (is_particle_stable_and_with_invalid_mass(smash_particle)) {
        // This modifies the momentum of the SMASH internal hadron to be
        // on-shell with the pole mass, but keeps the direction of the momentum
        // This is needed to ensure that the hadron can be found, ugly but
        // necessary
        smash_particle.set_4momentum(smash_particle.pole_mass(), p.threevec());
      }
      smash::FourVector r = smash_particle.position();
      exact_hadron.set_4position(r);
      exact_hadron.set_4momentum(p);
      exact_hadron_list.push_back(exact_hadron);
      found = true;
      break;  // Break the loop since we found the matching hadron
    }
    // If no matching hadron was found, we can print a warning
    if (!found) {
      JSWARN << "No matching hadron found for Jetscape hadron with ID: "
             << JS_had->pid();
    }
  }
  return exact_hadron_list;
}

void SmashInitialConditionWrapper::fill_JS_hadrons_from_smash_particles(
    const smash::Particles& smash_particles,
    std::vector<shared_ptr<Hadron>>& JS_hadrons) {
  JS_hadrons.clear();
  for (const auto& particle : smash_particles) {
    const int hadron_label = 0;
    const int hadron_status = 28;
    const int hadron_id = particle.pdgcode().get_decimal();
    smash::FourVector p = particle.momentum(), r = particle.position();
    const FourVector hadron_p(p.x1(), p.x2(), p.x3(), p.x0()),
        hadron_r(r.x1(), r.x2(), r.x3(), r.x0());
    const double hadron_mass = p.abs();
    const int charge = particle.type().charge();
    const int baryon_number = particle.type().baryon_number();
    const int strangeness = particle.type().strangeness();
    const auto history = particle.get_history();
    bool participant = false;
    if (history.collisions_per_particle > 0) {
      participant = true;
    }
    // Create a new Hadron object
    Hadron had(hadron_label, hadron_id, hadron_status, hadron_p, hadron_r,
               hadron_mass);

    // Set the properties using setter functions
    had.set_charge(charge);
    had.set_baryon_number(baryon_number);
    had.set_strangeness(strangeness);
    had.set_participant(participant);
    JS_hadrons.push_back(make_shared<Hadron>(had));
  }
}