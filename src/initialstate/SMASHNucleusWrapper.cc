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
// for to create a nucleus at rest
// -----------------------------------------------------------------------------

#include "SMASHNucleusWrapper.h"
#include "smash/particles.h"
#include "smash/library.h"

#include <iterator>

namespace Jetscape {

// Register the module with the base class
RegisterJetScapeModule<SMASHNucleusWrapper>
    SMASHNucleusWrapper::reg("SMASHNucleusWrapper");

SMASHNucleusWrapper::SMASHNucleusWrapper() {
  SetId("SMASHNucleusWrapper");
  smash_nucleus_hadrons_ = new smash::Particles();
}

SMASHNucleusWrapper::~SMASHNucleusWrapper() {
  if (smash_nucleus_hadrons_) {
    delete smash_nucleus_hadrons_;
    smash_nucleus_hadrons_ = nullptr;
  }
}

void SMASHNucleusWrapper::InitTask() {
  JSINFO << "Initializing SmashNucleusWrapper ...";
  std::string smash_config =
      GetXMLElementText({"IS", "SMASHNucleus", "SMASH_config_file"});
  std::string smash_hadron_list =
      GetXMLElementText({"IS", "SMASHNucleus", "SMASH_particles_file"});
  std::string smash_decays_list =
      GetXMLElementText({"IS", "SMASHNucleus", "SMASH_decaymodes_file"});
  // store tabulation to make use of it if SMASH is used multiple times
  std::string tabulations_path("./smash_tabulations");
  // if tabulations directory exists, delete it
  if (std::filesystem::exists(tabulations_path)) {
    std::filesystem::remove_all(tabulations_path);
  }
  const std::string smash_version(SMASH_VERSION);

  auto config1 = smash::setup_config_and_logging(
      smash_config, smash_hadron_list, smash_decays_list);
  auto config2 = smash::setup_config_and_logging(
      smash_config, smash_hadron_list, smash_decays_list);

  // Take care of the random seed. This will make SMASH results reproducible.
  auto random_seed = (*GetMt19937Generator())();
  config1.set_value({"General", "Randomseed"}, random_seed);
  config2.set_value({"General", "Randomseed"}, random_seed);

  int smash_projectile_protons =
      GetXMLElementInt({"IS", "SMASHNucleus", "Protons"});
  int smash_projectile_neutrons =
      GetXMLElementInt({"IS", "SMASHNucleus", "Neutrons"});
  int Fermi_momenta_tmp =
      GetXMLElementInt({"IS", "SMASHNucleus", "FermiMomenta"});
  if (Fermi_momenta_tmp == 1) {
    Fermi_momenta_ = true;
  } else if (Fermi_momenta_tmp == 0) {
    Fermi_momenta_ = false;
  } else {
    JSWARN << "FermiMomenta is not set to 0 or 1";
    exit(1);
  }

  std::map<int, int> smash_projectile{{2212, smash_projectile_protons},
                                      {2112, smash_projectile_neutrons}};
  config1.set_value({"Modi", "Collider", "Projectile", "Particles"},
                    smash_projectile);
  config2.set_value({"Modi", "Collider", "Projectile", "Particles"},
                    smash_projectile);

  // Check if the tabulations directory exists, if not initialize particles and
  // decays
  if (!std::filesystem::exists(tabulations_path)) {
    smash::initialize_particles_decays_and_tabulations(config1, smash_version,
                                                       tabulations_path);
  }

  // Try to get the Collider configuration sub-configuration
  smash::Configuration modus_config1 =
      config1.extract_sub_configuration({"Modi"});
  smash::Configuration modus_config2 =
      config2.extract_sub_configuration({"Modi"});

  smash_nucleus_ = make_shared<NucleusModus>(std::move(modus_config1),
                                             std::move(modus_config2));
  config1.clear();
  config2.clear();
  JSINFO << "Finish initializing SMASH nucleus creation";
}

void SMASHNucleusWrapper::ExecuteTask() {
  VERBOSE(2) << "SMASH Nucleus generation running: " << GetId() << "...";
  smash_nucleus_hadrons_->reset();
  smash_nucleus_->initial_conditions(smash_nucleus_hadrons_);
}

std::vector<Hadron> SMASHNucleusWrapper::GetCurrentHadronList() const {
  std::vector<Hadron> h_list;
  for (auto it = smash_nucleus_hadrons_->begin();
       it != smash_nucleus_hadrons_->end(); ++it) {
    const auto &particle = *it;

    const int hadron_label = 0;
    const int hadron_status = 29;
    const int hadron_id = particle.pdgcode().get_decimal();
    smash::FourVector p = particle.momentum();
    FourVector hadron_p;
    if (Fermi_momenta_) {
      hadron_p.Set(p.x1(), p.x2(), p.x3(), p.x0());
    } else {
      hadron_p.Set(0.0, 0.0, 0.0, nucleon_mass_);
    }
    smash::FourVector r = particle.position();
    const FourVector hadron_r(r.x1(), r.x2(), r.x3(), r.x0());

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

} // end namespace Jetscape