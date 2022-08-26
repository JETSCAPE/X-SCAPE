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
// -----------------------------------------
// This is a wrapper for SMASH hadronic afterburner with the JETSCAPE framework
// -----------------------------------------

#include "SmashWrapper.h"

#include "smash/particles.h"
#include "smash/library.h"

#include <math.h>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<SmashWrapper> SmashWrapper::reg("SMASH");

SmashWrapper::SmashWrapper() {
  SetId("SMASH");
}

void SmashWrapper::InitTask() {
  JSINFO << "SMASH: picking SMASH-specific configuration from xml file";
  std::string smash_config =
      GetXMLElementText({"Afterburner", "SMASH", "SMASH_config_file"});
  std::string smash_hadron_list =
      GetXMLElementText({"Afterburner", "SMASH", "SMASH_particles_file"});
  std::string smash_decays_list =
      GetXMLElementText({"Afterburner", "SMASH", "SMASH_decaymodes_file"});
  // output path is just dummy here, because no output from SMASH is foreseen
  boost::filesystem::path output_path("./smash_output");
  // do not store tabulation, which is achieved by an empty tabulations path
  std::string tabulations_path("");
  const std::string smash_version(SMASH_VERSION);

  auto config = smash::setup_config_and_logging(smash_config, smash_hadron_list,
                                                smash_decays_list);

  // Take care of the random seed. This will make SMASH results reproducible.
  auto random_seed = (*GetMt19937Generator())();
  config["General"]["Randomseed"] = random_seed;
  // Read in the rest of configuration
  if (IsTimeStepped()) {
    end_time_ = GetMainClock()->GetEndTime();
  } else {
    end_time_ = GetXMLElementDouble({"Afterburner", "SMASH", "end_time"});
  }
  config["General"]["End_Time"] = end_time_;
  JSINFO << "End time until which SMASH propagates is " << end_time_ << " fm/c";
  only_final_decays_ =
      GetXMLElementInt({"Afterburner", "SMASH", "only_decays"});
  if (only_final_decays_) {
    JSINFO << "SMASH will only perform resonance decays, no propagation";
  }

  smash::initialize_particles_decays_and_tabulations(config, smash_version,
                                                     tabulations_path);

  // Enforce timestep compatibility (temporarily)
  if (IsTimeStepped()) {
    const double delta_t_js = GetMainClock()->GetDeltaT();
    const double delta_t_sm = config.read({"General", "Delta_Time"});
    const double ts_rem = std::remainder(delta_t_js, delta_t_sm);
    const double ts_frac = delta_t_js / delta_t_sm;
    if (!(ts_rem < 1E-6 && ts_frac > 1.0)) {
      JSWARN << "Timesteps of SMASH (dt = " << delta_t_sm
             << ") and JETSCAPE (dt = " << delta_t_js << ") are incompabitle."
                "SMASH timesteps should be a half, a third, etc. from Jetscape's";
    }
  }

  JSINFO << "Seting up SMASH Experiment object";
  smash_experiment_ =
      make_shared<smash::Experiment<AfterburnerModus>>(config, output_path);
  JSINFO << "Finish initializing SMASH";
}

std::vector<std::vector<shared_ptr<Hadron>>> TestHadronList() {
  std::vector<std::vector<shared_ptr<Hadron>>> hadron_list_list;
  std::vector<shared_ptr<Hadron>> hadron_list;
  unsigned int nparticles = 6;
  for (unsigned int ipart = 0; ipart < nparticles; ipart++) {
    const int hadron_label = 0;
    const int hadron_status = 11;
    const int hadron_id = 111; // current_hadron.pid;
    const double hadron_mass = 0.138;
    const double pz = 0.1  * ipart;
    const double energy = std::sqrt(hadron_mass*hadron_mass + pz*pz);
    FourVector hadron_p(pz, 0.0, 0.0, energy);
    FourVector hadron_x(ipart, 0.0, 0.0, ipart);


    // create a JETSCAPE Hadron
    hadron_list.push_back(make_shared<Hadron>(hadron_label, hadron_id,
                                          hadron_status, hadron_p, hadron_x,
                                          hadron_mass));
  }
  hadron_list_list.push_back(hadron_list);
  return hadron_list_list;
}

void SmashWrapper::ExecuteTask() {
  AfterburnerModus *modus = smash_experiment_->modus();
  // This is necessary to correctly handle indices of particle sets from hydro.
  // Every hydro event creates a new structure like jetscape_hadrons_
  // with as many events in it as one has samples per hydro
  modus->reset_event_numbering();
  modus->jetscape_hadrons_ = GatherAfterburnerHadrons();
  // modus->jetscape_hadrons_ = TestHadronList();
  const int n_events = modus->jetscape_hadrons_.size();
  JSINFO << "SMASH: obtained " << n_events << " events from particlization";
  for (unsigned int i = 0; i < n_events; i++) {
    JSINFO << "Event " << i << " SMASH starts with "
           << modus->jetscape_hadrons_[i].size() << " particles.";
    InitPerEvent();
    CalculateTimeTask();
    FinishPerEvent();
  }
}


void SmashWrapper::InitPerEvent() {
  if (IsTimeStepped()) {
    JSINFO << "Initalizing new time-stepped SMASH  event ...";
    AfterburnerModus *modus = smash_experiment_->modus();
    modus->reset_event_numbering();
    // modus->jetscape_hadrons_ = GatherAfterburnerHadrons();
    modus->jetscape_hadrons_ = TestHadronList();

    const int n_events = modus->jetscape_hadrons_.size();
    if (n_events > 1) {
      JSWARN << "In timestep mode SMASH only propagates one (= the first "
                "soft_particlization) event at the moment. No oversampling possible.";
    }
  }
  smash_experiment_->initialize_new_event();
}


smash::ParticleList SmashWrapper::convert_to_plist(std::vector<shared_ptr<Hadron>>& JS_hadrons) {
  // TODO Merge/generalize this to be also used in JS_hadrons_to_smash_particles()
  smash::ParticleList new_particles;
  for (const auto JS_had : JS_hadrons) {
    const FourVector p = JS_had->p_in();
    const FourVector r = JS_had->x_in();
    smash::ParticleData new_p{smash::ParticleType::find(smash::PdgCode::from_decimal(JS_had->pid()))};
    new_p.set_4position(smash::FourVector(r.t(), r.x(), r.y(), r.z()));
    new_p.set_4momentum(p.t(), p.x(), p.y(), p.z());

    std::cout << "particle to be put into smash: " << new_p << '\n';

    new_particles.push_back(new_p);
  }
  return new_particles;
}


void SmashWrapper::CalculateTimeTask() {

  // TODO(stdnmr) Get new hadrons from BDM here
  std::vector<shared_ptr<Hadron>> new_JS_hadrons = {}; 

  const double until_time = IsTimeStepped() ? GetMainClock()->GetCurrentTime() : end_time_;
  JSINFO << "Propgating SMASH until t = " << until_time;
  if (!only_final_decays_) {
    smash_experiment_->run_time_evolution(until_time, convert_to_plist(new_JS_hadrons));
  }
}

void SmashWrapper::FinishPerEvent() {
  JSINFO << "Finishing SMASH event...";

  AfterburnerModus *modus = smash_experiment_->modus();
  smash::Particles *smash_particles = smash_experiment_->first_ensemble();
  int ev_no = modus->current_event_number();

  smash_experiment_->do_final_decays();
  smash_experiment_->final_output();
  smash_particles_to_JS_hadrons(*smash_particles,
                                modus->jetscape_hadrons_[ev_no - 1]);
  JSINFO << modus->jetscape_hadrons_[ev_no - 1].size()
         << " hadrons from SMASH.";
}

void SmashWrapper::WriteTask(weak_ptr<JetScapeWriter> w) {
  JSINFO << "SMASH hadronic afterburner printout";
  auto f = w.lock();
  if (!f) {
    return;
  }
  AfterburnerModus *modus = smash_experiment_->modus();
  f->WriteComment("JetScape module: " + GetId());
  for (const auto &event : modus->jetscape_hadrons_) {
    int i = -1;
    for (const auto hadron : event) {
      f->WriteWhiteSpace("[" + to_string(++i) + "] H");
      f->Write(hadron);
    }
  }
}

void AfterburnerModus::JS_hadrons_to_smash_particles(
    const std::vector<shared_ptr<Hadron>> &JS_hadrons,
    smash::Particles &smash_particles) {
  smash_particles.reset();
  for (const auto JS_hadron : JS_hadrons) {
    const FourVector p = JS_hadron->p_in();
    const FourVector r = JS_hadron->x_in();
    const double mass = JS_hadron->restmass();
    smash::PdgCode pdgcode = smash::PdgCode::from_decimal(JS_hadron->pid());
    this->try_create_particle(smash_particles, pdgcode, r.t(), r.x(), r.y(),
                              r.z(), mass, JS_hadron->e(), JS_hadron->px(),
                              JS_hadron->py(), JS_hadron->pz());
  }
}

void SmashWrapper::smash_particles_to_JS_hadrons(
    const smash::Particles &smash_particles,
    std::vector<shared_ptr<Hadron>> &JS_hadrons) {
  JS_hadrons.clear();
  for (const auto &particle : smash_particles) {
    const int hadron_label = 0;
    const int hadron_status = -1;
    const int hadron_id = particle.pdgcode().get_decimal();
    smash::FourVector p = particle.momentum(), r = particle.position();
    const FourVector hadron_p(p.x1(), p.x2(), p.x3(), p.x0()),
        hadron_r(r.x1(), r.x2(), r.x3(), r.x0());
    const double hadron_mass = p.abs();
    JS_hadrons.push_back(make_shared<Hadron>(hadron_label, hadron_id,
                                             hadron_status, hadron_p, hadron_r,
                                             hadron_mass));
  }
}
