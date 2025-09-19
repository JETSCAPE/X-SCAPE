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
#include "smash/input_keys.h"
#include "smash/particles.h"
#include "smash/library.h"

#include <iterator>

namespace YAML {
template <>
struct convert<smash::PdgCode> {
  static Node encode(const smash::PdgCode &pdg) {
    // encode as string, e.g., "pdg:221"
    // assuming you can convert it to string via operator<< or some method
    std::stringstream ss;
    ss << pdg;  // if operator<< is implemented
    return Node(ss.str());
  }

  static bool decode(const Node &node, smash::PdgCode &pdg) {
    if (!node.IsScalar()) return false;

    // parse integer from string
    int code = std::stoi(node.as<std::string>());
    pdg = smash::PdgCode(code);
    return true;
  }
};
}  // namespace YAML

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
  config1.set_value(smash::InputKeys::gen_randomseed, random_seed);
  config2.set_value(smash::InputKeys::gen_randomseed, random_seed);

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

  auto to_pdg_map = [](const std::map<int,int>& input) {
    std::map<smash::PdgCode,int> output;
    for (const auto& [pdg_int, count] : input) {
      output.emplace(smash::PdgCode(pdg_int), count);
    }
    return output;
  };
  std::map<int, int> smash_projectile{{2212, smash_projectile_protons},
                                      {2112, smash_projectile_neutrons}};
  std::map<smash::PdgCode,int> projectile_pdg = to_pdg_map(smash_projectile);
  config1.set_value(smash::InputKeys::modi_collider_projectile_particles,
                    projectile_pdg);
  config2.set_value(smash::InputKeys::modi_collider_projectile_particles,
                    projectile_pdg);

  // Check if the tabulations directory exists, if not initialize particles and
  // decays
  if (!std::filesystem::exists(tabulations_path)) {
    smash::initialize_particles_decays_and_tabulations(config1, smash_version,
                                                       tabulations_path);
  }

  // Try to get the Collider configuration sub-configuration
  smash::Configuration modus_config1 =
      config1.extract_complete_sub_configuration(
        smash::InputSections::m_collider);
  smash::Configuration modus_config2 =
      config2.extract_complete_sub_configuration(
        smash::InputSections::m_collider);

  smash_nucleus_ = make_shared<NucleusModus>(std::move(modus_config1),
                                             std::move(modus_config2));
  config1.clear();
  config2.clear();
  JSINFO << "Finish initializing SMASH nucleus creation";

  // Initialize the NucleonRadiusBlackDisk
  nucleon_radius_black_disk_ =
      GetXMLElementDouble({"IS", "SMASHNucleus", "NucleonRadiusBlackDisk"});
  if (nucleon_radius_black_disk_ <= 0.0) {
    JSWARN << "NucleonRadiusBlackDisk is not set to a positive value, using 0.6 fm as the default.";
    nucleon_radius_black_disk_ = 0.6;
  }
}

void SMASHNucleusWrapper::ExecuteTask() {
  VERBOSE(2) << "SMASH Nucleus generation running: " << GetId() << "...";
  smash_nucleus_hadrons_->reset();
  smash_nucleus_->initial_conditions(smash_nucleus_hadrons_);
  // Store the hadrons in the hadrons_ vector
  // This is done to have a copy of the hadrons in the SMASHNucleusWrapper
  StoreHadronsInWrapper();
}

void SMASHNucleusWrapper::StoreHadronsInWrapper() {
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
    hadrons_.push_back(Hadron(hadron_label, hadron_id, hadron_status, hadron_p,
                            hadron_r, hadron_mass, charge, baryon_number,
                            strangeness, participant));
  }
}

std::vector<Hadron> SMASHNucleusWrapper::GetCurrentHadronList() const {
  return hadrons_;
}

bool SMASHNucleusWrapper::IsHadronAtPosition(double t, double x,
                                        double y, double z) const {
  // Note: All nucleons are initialized at t=0; the 't' argument is ignored.
  // This function only checks if there is at least one hadron spatially
  // within nucleon_radius_black_disk_ of the given (x, y, z).
  for (const auto &hadron : GetCurrentHadronList()) {
    if (std::abs(hadron.x_in().t() - t) > rounding_error) {
      continue;
    }
    double dx = hadron.x_in().x() - x;
    double dy = hadron.x_in().y() - y;
    double dz = hadron.x_in().z() - z;
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (distance <= nucleon_radius_black_disk_) {
      return true;
    }
  }
  return false;
}

std::tuple<double, double, double, double> SMASHNucleusWrapper::BoostCoordinates(
  double x0, double x1, double x2, double x3, 
  double vx, double vy, double vz) const {
  const double beta2 = vx * vx + vy * vy + vz * vz;
  double gamma;
  if (std::sqrt(beta2) < 1.0) {
    gamma = 1.0 / std::sqrt(1.0 - beta2);
  } else {
    gamma = a_very_large_number;
    JSWARN << "Boost velocity is larger than 1, setting gamma to a very large number.";
  }
  // create array of the original coordinates
  double original_coordinates[4] = {x0, x1, x2, x3};
  // define 4x4 Lorentz transformation matrix
  double lorentz_matrix[4][4] = {
    {gamma, -gamma * vx, -gamma * vy, -gamma * vz},
    {-gamma * vx, (1 + ((gamma - 1) * vx * vx / beta2)), (gamma - 1) * vx * vy / beta2, (gamma - 1) * vx * vz / beta2},
    {-gamma * vy, (gamma - 1) * vy * vx / beta2, (1 + ((gamma - 1) * vy * vy / beta2)), (gamma - 1) * vy * vz / beta2},
    {-gamma * vz, (gamma - 1) * vz * vx / beta2, (gamma - 1) * vz * vy / beta2, (1 + ((gamma - 1) * vz * vz / beta2))}
  };
  // create an array to hold the boosted coordinates
  double boosted_coordinates[4] = {0.0, 0.0, 0.0, 0.0};
  // Perform the Lorentz transformation
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      boosted_coordinates[i] += lorentz_matrix[i][j] * original_coordinates[j];
    }
  }
  // Create the tuple with the boosted (primed) coordinates
  return std::make_tuple(boosted_coordinates[0], boosted_coordinates[1], 
                          boosted_coordinates[2], boosted_coordinates[3]);
}

std::vector<Hadron> SMASHNucleusWrapper::GetCurrentHadronListBoosted(double vx, double vy, double vz) {
  for (auto &hadron : hadrons_) {
    // Boost the hadron's position and momentum
    const FourVector r = hadron.x_in();
    const FourVector p = hadron.p_in();

    // Use the Lorentz transformation to boost the hadron's position
    auto boosted_positions = BoostCoordinates(r.t(), r.x(), r.y(), r.z(), vx, vy, vz);
    const double t_prime = std::get<0>(boosted_positions);
    const double x_prime = std::get<1>(boosted_positions);
    const double y_prime = std::get<2>(boosted_positions);
    const double z_prime = std::get<3>(boosted_positions);

    // Boost the hadron's momentum
    auto boosted_momentum = BoostCoordinates(p.t(), p.x(), p.y(), p.z(), vx, vy, vz);
    const double t_prime_mom = std::get<0>(boosted_momentum);
    const double x_prime_mom = std::get<1>(boosted_momentum);
    const double y_prime_mom = std::get<2>(boosted_momentum);
    const double z_prime_mom = std::get<3>(boosted_momentum);
    const FourVector boosted_momentum_vector(x_prime_mom, y_prime_mom, 
                                              z_prime_mom, t_prime_mom);
    // Update the hadron's position and momentum
    hadron.reset_momentum(boosted_momentum_vector);
    double new_x[4] = {t_prime, x_prime, y_prime, z_prime};
    hadron.set_x(new_x);
  }
  return hadrons_;
}



} // end namespace Jetscape