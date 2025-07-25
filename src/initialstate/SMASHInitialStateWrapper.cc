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
// -----------------------------------------------------------------------------
// This is a wrapper for SMASH hadronic transport with the JETSCAPE framework
// for the collider modus of SMASH to generate initial conditions
// -----------------------------------------------------------------------------

#include "SMASHInitialStateWrapper.h"
#include "Transport.h"

#include "smash/particles.h"
#include "smash/library.h"

#include <math.h>
#include <string>
#include <map>
#include <filesystem>

#include <boost/lexical_cast.hpp>

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
  // store tabulation to make use of it if SMASH is used multiple times
  std::string tabulations_path("./smash_tabulations");
  // if tabulations directory exists, delete it
  if (std::filesystem::exists(tabulations_path)) {
    std::filesystem::remove_all(tabulations_path);
  }
  const std::string smash_version(SMASH_VERSION);

  auto config = smash::setup_config_and_logging(smash_config, 
                                                smash_hadron_list,
                                                smash_decays_list);

  // Take care of the random seed. This will make SMASH results reproducible.
  auto random_seed = (*GetMt19937Generator())();
  config.set_value({"General","Randomseed"}, random_seed);
  // Read in the rest of configuration
  if (IsTimeStepped()) {
    end_time_ = GetMainClock()->GetEndTime();
  } else {
    end_time_ = GetXMLElementDouble({"IS", "SMASH", "End_Time"});
  }
  config.set_value({"General","End_Time"}, end_time_);
  JSINFO << "End time until which SMASH initial condition propagates is " 
        << end_time_ << " fm/c";

  int smash_projectile_protons =
      GetXMLElementInt({"IS", "SMASH", "ProjectileProtons"});
  int smash_projectile_neutrons =
      GetXMLElementInt({"IS", "SMASH", "ProjectileNeutrons"});
  int smash_target_protons =
      GetXMLElementInt({"IS", "SMASH", "TargetProtons"});
  int smash_target_neutrons =
      GetXMLElementInt({"IS", "SMASH", "TargetNeutrons"});

  std::map<int, int> smash_projectile{{2212, smash_projectile_protons}, 
                                      {2112, smash_projectile_neutrons}};
  config.set_value({"Modi","Collider","Projectile","Particles"},
                    smash_projectile);

  std::map<int, int> smash_target{{2212, smash_target_protons}, 
                                  {2112, smash_target_neutrons}};
  config.set_value({"Modi","Collider","Target","Particles"},smash_target);

  double smash_sqrtsnn = GetXMLElementDouble({"IS", "SMASH", "Sqrtsnn"});
  config.set_value({"Modi","Collider","Sqrtsnn"}, smash_sqrtsnn);

  std::string smash_FermiMotion = 
        GetXMLElementText({"IS", "SMASH", "FermiMotion"});
  config.set_value({"Modi","Collider","Fermi_Motion"},smash_FermiMotion);

  bool smash_CollisionsWithinNucleus = 
        GetXMLElementInt({"IS", "SMASH", "CollisionsWithinNucleus"});
  config.set_value({"Modi","Collider","Collisions_Within_Nucleus"},
                    smash_CollisionsWithinNucleus);

  bool smash_impact_react_plane = 
          GetXMLElementInt({"IS", "SMASH", "Impact", "RandomReactionPlane"});
  config.set_value({"Modi","Collider","Impact","Random_Reaction_Plane"},
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
    config.set_value({"Modi","Collider","Impact","Value"}, smash_impact_val);
  } else if (smash_impact_param_mode == 1) {
    config.set_value({"Modi","Collider","Impact","Sample"},smash_impact_sample);
    const std::array<double, 2> smash_impact_range = {smash_impact_valMin,
                                                      smash_impact_valMax};
    config.set_value({"Modi","Collider","Impact","Range"}, smash_impact_range);
  } else {
    JSWARN << "This SMASH impact parameter mode does not exist in Jetscape";
    exit(1);
  }

  // Check if the tabulations directory exists, if not initialize particles and 
  // decays
  if (!std::filesystem::exists(tabulations_path)) {
    smash::initialize_particles_decays_and_tabulations(config, smash_version,
                                                     tabulations_path);
  }

  const double delta_t_sm = GetXMLElementDouble({"IS", "SMASH", "Delta_Time"});
  config.set_value({"General", "Delta_Time"}, delta_t_sm);

  // Enforce timestep compatibility (temporarily)
  if (IsTimeStepped()) {
    const double delta_t_js = GetMainClock()->GetDeltaT();
    const double ts_rem = std::remainder(delta_t_js, delta_t_sm);
    const double ts_frac = delta_t_js / delta_t_sm;
    if (!(ts_rem < 1E-6 && ts_frac > 1.0)) {
      JSWARN << "Timesteps of SMASH (dt = " << delta_t_sm
             << ") and JETSCAPE (dt = " << delta_t_js << ") are incompatible."
                "SMASH IC timesteps should be a half, a third, etc. from JETSCAPE's";
    }
  }
  smash_collider_experiment_ =
      make_shared<smash::Experiment<smash::ColliderModus>>(config, output_path);
  config.clear();
  JSINFO << "Finish initializing SMASH initial condition";
}

void SmashInitialConditionWrapper::ExecuteTask() {
  VERBOSE(2) << "SMASH initial condition running: " << GetId() << "...";
  InitPerEvent();
  CalculateTimeTask();
  FinishPerEvent();
}

void SmashInitialConditionWrapper::InitPerEvent() {
  VERBOSE(3) << "Initializing SMASH initial condition event...";
  smash_collider_experiment_->initialize_new_event();
}

void SmashInitialConditionWrapper::CalculateTimeTask() {
  std::vector<shared_ptr<Hadron>> hadrons_to_add = Transport::GetTimestepParticlizationHadrons();
  std::vector<shared_ptr<Hadron>> hadrons_to_remove = Transport::GetTimestepHadronsToRemove();

  //VERBOSE(2) << "SMASH initial condition got " << hadrons_to_add.size()  << " hadrons in this timestep.";
  //VERBOSE(2) << "SMASH initial condition removed " << hadrons_to_remove.size() << " hadrons in this timestep.";
  //JSINFO << "SMASH initial condition got " << hadrons_to_add.size()  << " hadrons in this timestep.";
  //JSINFO << "SMASH initial condition removed " << hadrons_to_remove.size() << " hadrons in this timestep.";

  // print the properties of the hadrons to be added
  /*for (const auto& hadron : hadrons_to_add) {
    FourVector p = hadron->p_in();
    FourVector r = hadron->x_in();
    JSINFO << "Adding hadron: " << hadron->pid() << " at (x, y, z, t) = ("
          << r.x() << ", " << r.y() << ", " << r.z() << ", " << r.t() << ") fm"
          << " with (px, py, pz, E) = (" << p.x() << ", " << p.y() << ", " << p.z() << ", " << p.t() << ") GeV";
  }*/

  const double until_time = IsTimeStepped() ? GetMainClock()->GetCurrentTime() : end_time_;
  //JSINFO << "Propagating SMASH IC until t = " << until_time;
  //JSINFO << "End time until which SMASH initial condition propagates is " << end_time_ << " fm/c";

  smash::ParticleList add_list = get_smash_plist_from_JS_hadrons(hadrons_to_add);
  smash::ParticleList remove_list = get_smash_plist_from_JS_hadrons(hadrons_to_remove);
  smash_collider_experiment_->run_time_evolution(until_time,std::move(add_list),std::move(remove_list));
}

void SmashInitialConditionWrapper::FinishPerEvent() {
  JSINFO << "Finishing SMASH initial condition event...";
  event_number_++;

  // SMASH within JETSCAPE only works with one (the first) ensemble
  smash::Particles *smash_particles = smash_collider_experiment_->first_ensemble();
  int ev_no = current_event_number();

  smash_collider_experiment_->do_final_decays();
  smash_collider_experiment_->final_output();

  if (ev_no > jetscape_hadrons_.size()) {
    jetscape_hadrons_.resize(ev_no);
  }

  fill_JS_hadrons_from_smash_particles(*smash_particles,
                                      jetscape_hadrons_[ev_no - 1]);

  smash_collider_experiment_->increase_event_number(); // internal SMASH event counter
  JSINFO << jetscape_hadrons_[ev_no - 1].size() << " hadrons from SMASH initial condition.";
  JSINFO << "Finished SMASH collider event...";
}


void SmashInitialConditionWrapper::WriteTask(weak_ptr<JetScapeWriter> w) {
  JSINFO << "SMASH initial condition printout";
  auto f = w.lock();
  if (!f) {
    return;
  }
  f->WriteComment("JetScape module: " + GetId());
  for (const auto &event : jetscape_hadrons_) {
    int i = -1;
    for (const auto hadron : event) {
      f->WriteWhiteSpace("[" + to_string(++i) + "] H");
      f->Write(hadron);
    }
  }
}

std::vector<Hadron> SmashInitialConditionWrapper::GetCurrentHadronList() const {
  std::vector<Hadron> h_list;
  smash::Particles* smash_particles = smash_collider_experiment_->first_ensemble();

  for (const auto &particle : *smash_particles) {
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
    h_list.push_back(Hadron(hadron_label, hadron_id, hadron_status, hadron_p, 
                            hadron_r, hadron_mass, charge, baryon_number, 
                            strangeness, participant));
  }
  return h_list;
}

smash::ParticleList SmashInitialConditionWrapper::get_smash_plist_from_JS_hadrons(
                    const std::vector<shared_ptr<Hadron>>& JS_hadrons) {
  smash::ParticleList new_particles;
  for (const auto& JS_had : JS_hadrons) {
    const FourVector p = JS_had->p_in();
    const FourVector r = JS_had->x_in();
    smash::ParticleData new_p{smash::ParticleType::find(smash::PdgCode::from_decimal(JS_had->pid()))};
    new_p.set_4position(smash::FourVector(r.t(), r.x(), r.y(), r.z()));
    new_p.set_4momentum(smash::FourVector(p.t(), p.x(), p.y(), p.z()));
    new_particles.push_back(new_p);
  }
  return new_particles;
}

void SmashInitialConditionWrapper::fill_JS_hadrons_from_smash_particles(
    const smash::Particles &smash_particles,
    std::vector<shared_ptr<Hadron>> &JS_hadrons) {
  JS_hadrons.clear();
  for (const auto &particle : smash_particles) {
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