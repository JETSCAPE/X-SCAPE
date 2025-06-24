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

#include <string>
#include <fstream>
#include <stdio.h>
#include <sys/stat.h>


#include <vector>

#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include "MCGlauberWrapper.h"

// Register the module with the base class
RegisterJetScapeModule<MCGlauberWrapper> MCGlauberWrapper::reg("MCGlauber");

MCGlauberWrapper::MCGlauberWrapper() {
    SetId("MCGlauber");
    event_id_ = 0;
}


void MCGlauberWrapper::InitTask() {
    auto ran_seed = (*GetMt19937Generator())();
    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini) {
        JSWARN << "The MCGlauber pointer to the initial state is not set.";
    }

    int argc = 0;
    char* argv[1];
    mc_gen_ = std::shared_ptr<MCGlb::EventGenerator>(
        new MCGlb::EventGenerator("mcglauber.input", argc, argv, ran_seed));

    // overwrite input options
    int para_temp_int;
    double para_temp_double;
    std::string para_temp_string;

    para_temp_int = GetXMLElementInt(
        {"IS", "MCGlauber", "generateOnlyPositions"}, true);
    if (para_temp_int == 1) {
        generateOnlyPositions_ = true;
    } else {
        generateOnlyPositions_ = false;
    }

    para_temp_string = (
        GetXMLElementText({"IS", "MCGlauber", "projectile"}));
    mc_gen_->set_parameter("Projectile", para_temp_string);

    para_temp_string = (
        GetXMLElementText({"IS", "MCGlauber", "target"}));
    mc_gen_->set_parameter("Target", para_temp_string);

    para_temp_double = (
        GetXMLElementDouble({"IS", "MCGlauber", "sqrts"}));
    mc_gen_->set_parameter("roots", para_temp_double);

    para_temp_double = (
        GetXMLElementDouble({"IS", "MCGlauber", "b_min"}));
    mc_gen_->set_parameter("b_min", para_temp_double);

    para_temp_double = (
        GetXMLElementDouble({"IS", "MCGlauber", "b_max"}));
    mc_gen_->set_parameter("b_max", para_temp_double);

    para_temp_double = (
        GetXMLElementDouble({"IS", "MCGlauber", "cenMin"}));
    mc_gen_->set_parameter("cenMin", para_temp_double);

    para_temp_double = (
        GetXMLElementDouble({"IS", "MCGlauber", "cenMax"}));
    mc_gen_->set_parameter("cenMax", para_temp_double);

    para_temp_int = (
        GetXMLElementInt({"IS", "MCGlauber", 
            "nucleon_configuration_from_file"}));
    mc_gen_->set_parameter("nucleon_configuration_from_file", para_temp_int);

    para_temp_int = (GetXMLElementInt({"IS", "MCGlauber", "useQuarks"}));
    mc_gen_->set_parameter("useQuarks", para_temp_int);

    para_temp_double = (GetXMLElementDouble({"IS", "MCGlauber", "Q2"}));
    mc_gen_->set_parameter("Q2", para_temp_double);

    para_temp_double = (
        GetXMLElementDouble({"IS", "MCGlauber", "shadowing_factor"}));
    mc_gen_->set_parameter("shadowing_factor", para_temp_double);

    para_temp_int = (GetXMLElementInt({"IS", "MCGlauber", "baryon_junctions"}));
    mc_gen_->set_parameter("baryon_junctions", para_temp_int);

    para_temp_int = (
        GetXMLElementInt({"IS", "MCGlauber", "N_sea_partons"}));
    mc_gen_->set_parameter("N_sea_partons", para_temp_int);

    para_temp_double = (
        GetXMLElementDouble({"IS", "MCGlauber", "d_min"}));
    mc_gen_->set_parameter("d_min", para_temp_double);

    para_temp_double = (
        GetXMLElementDouble({"IS", "MCGlauber", "lambdaB"}));
    mc_gen_->set_parameter("lambdaB", para_temp_double);

    para_temp_double = (GetXMLElementDouble({"IS", "MCGlauber", "BG"}));
    mc_gen_->set_parameter("BG", para_temp_double);

    para_temp_int = (
        GetXMLElementInt({"IS", "MCGlauber", "Subtract_hard_momentum"}));
    mc_gen_->set_parameter("Subtract_hard_momentum", para_temp_int);

    para_temp_double = (GetXMLElementDouble({"IS", "MCGlauber", "lambdaBs"}));
    mc_gen_->set_parameter("lambdaBs", para_temp_double);

    para_temp_double = (
        GetXMLElementDouble({"IS", "MCGlauber", "baryonInStringProb"}));
    mc_gen_->set_parameter("baryonInStringProb", para_temp_double);

    para_temp_int = (GetXMLElementInt({"IS", "MCGlauber", 
        "fluct_Nstrings_per_NN_collision"}));
    mc_gen_->set_parameter("fluct_Nstrings_per_NN_collision", para_temp_int);

    para_temp_int = (GetXMLElementInt({"IS", "MCGlauber", 
        "QCD_string_production_mode"}));
    mc_gen_->set_parameter("QCD_string_production_mode", para_temp_int);

    para_temp_int = (
        GetXMLElementInt({"IS", "MCGlauber", "rapidity_loss_method"}));
    if (para_temp_int == 4) {
        mc_gen_->set_parameter("rapidity_loss_method", para_temp_int);

        para_temp_double = (
            GetXMLElementDouble({"IS", "MCGlauber", "ylossParam4At2"}));
        mc_gen_->set_parameter("ylossParam4At2", para_temp_double);

        para_temp_double = (
            GetXMLElementDouble({"IS", "MCGlauber", "ylossParam4At4"}));
        mc_gen_->set_parameter("ylossParam4At4", para_temp_double);

        para_temp_double = (
            GetXMLElementDouble({"IS", "MCGlauber", "ylossParam4At6"}));
        mc_gen_->set_parameter("ylossParam4At6", para_temp_double);

        para_temp_double = (
            GetXMLElementDouble({"IS", "MCGlauber", "ylossParam4At10"}));
        mc_gen_->set_parameter("ylossParam4At10", para_temp_double);

        para_temp_double = (
            GetXMLElementDouble({"IS", "MCGlauber", "ylossParam4var"}));
        mc_gen_->set_parameter("ylossParam4var", para_temp_double);
    } else {
        JSWARN << "The option " << para_temp_int
               << " for rapidity loss method is not supported. "
               << "Please specify these parameters in mcglauber.input!";
    }

    para_temp_double = (GetXMLElementDouble({"IS", "MCGlauber",
                                      "remnant_energy_loss_fraction"}));
    mc_gen_->set_parameter("remnant_energy_loss_fraction", para_temp_double);

    para_temp_int = (GetXMLElementInt({"IS", "MCGlauber", 
        "fluctuation_remnant_energy_loss_fraction"}));
    mc_gen_->set_parameter("fluctuation_remnant_energy_loss_fraction", 
                           para_temp_int);

    para_temp_double = (GetXMLElementDouble({"IS", "MCGlauber", 
        "remnant_energy_loss_fraction_var"}));
    mc_gen_->set_parameter("remnant_energy_loss_fraction_var", 
                           para_temp_double);

    para_temp_int = (GetXMLElementInt({"IS", "MCGlauber", 
        "evolve_QCD_string_mode"}));
    mc_gen_->set_parameter("evolve_QCD_string_mode", para_temp_int);

    para_temp_double = (GetXMLElementDouble({"IS", "MCGlauber", 
        "tau_form_mean"}));
    mc_gen_->set_parameter("tau_form_mean", para_temp_double);

    para_temp_double = (GetXMLElementDouble({"IS", "MCGlauber", 
        "tau_form_fluct_gamma_beta"}));
    mc_gen_->set_parameter("tau_form_fluct_gamma_beta", para_temp_double);

    mc_gen_->print_parameter_list();
    mc_gen_->initializeGlauberModel();
}


void MCGlauberWrapper::ClearTask() {
    VERBOSE(1) << "clear initial condition vectors";
    binary_collision_t_.clear();
    binary_collision_x_.clear();
    binary_collision_y_.clear();
    binary_collision_z_.clear();
    QCDStringList_.clear();
}


void MCGlauberWrapper::ExecuteTask() {
    ClearTask();
    VERBOSE(1) << "Run 3DMCGlauber to generate initial hard positions "
                     << "...";

    if (generateOnlyPositions_) {
        try {
            int iparticle=0;
            mc_gen_->generate_pre_events(); // generate one 3DGlauber event
            std::vector<MCGlb::CollisionEvent> collisionEvents = (
                mc_gen_->get_CollisionEventvector());
            ncoll_ = collisionEvents.size();
            rand_int_ptr_ = (
                std::make_shared<std::uniform_int_distribution<int>>(0, ncoll_-1));
            while (iparticle < ncoll_) {
                 auto xvec = (
                    collisionEvents[iparticle].get_collision_position());
                 binary_collision_t_.push_back(xvec[0]);
                 binary_collision_x_.push_back(xvec[1]);
                 binary_collision_y_.push_back(xvec[2]);
                 binary_collision_z_.push_back(xvec[3]);
                 iparticle++;
            }
            event_id_++;
        } catch (std::exception &err) {
            Jetscape::JSWARN << err.what();
            std::exit(-1);
        }
    } else {
        try {
            int iparticle=0;
            mc_gen_->generate_pre_events(); // TODO: change this function to generate_full_events
            std::vector<MCGlb::CollisionEvent> collisionEvents = (
                mc_gen_->get_CollisionEventvector());
            ncoll_ = collisionEvents.size();
            rand_int_ptr_ = (
                std::make_shared<std::uniform_int_distribution<int>>(0, ncoll_-1));
            while (iparticle < ncoll_) {
                 auto xvec = (
                    collisionEvents[iparticle].get_collision_position());
                 binary_collision_t_.push_back(xvec[0]);
                 binary_collision_x_.push_back(xvec[1]);
                 binary_collision_y_.push_back(xvec[2]);
                 binary_collision_z_.push_back(xvec[3]);
                 iparticle++;
            }
            event_id_++;
            ini->GenerateStrings();
        } catch (std::exception &err) {
            Jetscape::JSWARN << err.what();
            std::exit(-1);
        }
    }
}

void MCGlauberWrapper::SampleABinaryCollisionPoint(
        double &t, double &x, double &y, double &z) {
    const int rand_idx = (*rand_int_ptr_)(*GetMt19937Generator());
    t = binary_collision_t_[rand_idx];
    x = binary_collision_x_[rand_idx];
    y = binary_collision_y_[rand_idx];
    z = binary_collision_z_[rand_idx];
}


double MCGlauberWrapper::Get_total_nucleon_density_lab(
        double t, double x, double y, double z) {
    // get the summation of nucleon density over projectile and target
    // at the Lab frame. the unit is 1/fm^3
    return (mc_gen_->MCGlb_nucleon_density(t, x, y, z));
}


double MCGlauberWrapper::Get_target_nucleon_density_lab(
        double t, double x, double y, double z) {
    // get the target nucleon density at the Lab frame, 
    // target: moves to the -z direction. the unit is 1/fm^3
    return(mc_gen_->MCGlb_target_nucleon_density(t, x, y, z));
}


double MCGlauberWrapper::Get_projectile_nucleon_density_lab(
        double t, double x, double y, double z) {
    // get the projectile nucleon density at the Lab frame, 
    // projectile: moves to the +z direction. the unit is 1/fm^3
    return(mc_gen_->MCGlb_projectile_nucleon_density(t, x, y, z));
}

std::vector<double> MCGlauberWrapper::Get_projectile_nucleon_z_lab() {
    // get the z coordinate of all projectile nucleons at the Lab frame, 
    return(mc_gen_->MCGlb_projectile_nucleon_z());
}


std::vector<double> MCGlauberWrapper::Get_target_nucleon_z_lab() {
    // get the z coordinate of all target nucleons at the Lab frame, 
    return(mc_gen_->MCGlb_target_nucleon_z());
}

void MCGlauberWrapper::OutputHardCollisionPosition(double t, double x, double y, 
                                                                    double z) {
    hard_parton_t_ = t;
    hard_parton_x_ = x;
    hard_parton_y_ = y;
    hard_parton_z_ = z;
}

void MCGlauberWrapper::ClearHardPartonMomentum(){

    proj_parton_e_  = 0.0;
    proj_parton_px_ = 0.0;
    proj_parton_py_ = 0.0;
    proj_parton_pz_ = 0.0;
    targ_parton_e_  = 0.0;
    targ_parton_px_ = 0.0;
    targ_parton_py_ = 0.0;
    targ_parton_pz_ = 0.0;
}

void MCGlauberWrapper::OutputHardPartonMomentum(double E, double px, double py, double pz,
                                                int direction, double P_A) {
    // JSWARN <<  MAGENTA << " Pushing hard momentum to MCGlauber ";
    if (direction == 1) {
        proj_parton_e_ += E;
        proj_parton_px_ += px;
        proj_parton_py_ += py;
        proj_parton_pz_ += pz;
        // JSINFO <<  MAGENTA << " proj_parton_e_ " << proj_parton_e_; 
        // JSINFO <<  MAGENTA << " proj_parton_px_ " << proj_parton_px_;
        // JSINFO <<  MAGENTA << " proj_parton_py_ " << proj_parton_py_;
        // JSINFO <<  MAGENTA << " proj_parton_pz_ " << proj_parton_pz_;
    } else {
        targ_parton_e_ += E;
        targ_parton_px_ += px;
        targ_parton_py_ += py;
        targ_parton_pz_ += pz;
        // JSINFO <<  MAGENTA << " targ_parton_e_ " << targ_parton_e_;
        // JSINFO <<  MAGENTA << " targ_parton_px_ " << targ_parton_px_;
        // JSINFO <<  MAGENTA << " targ_parton_py_ " << targ_parton_py_;
        // JSINFO <<  MAGENTA << " targ_parton_pz_ " << targ_parton_pz_;
    }


    VERBOSE(2) << BOLDYELLOW << " proj_parton_e_ " << proj_parton_e_ 
                             << " proj_parton_pz_ " << proj_parton_pz_ 
                             << " targ_parton_e_ " << targ_parton_e_ 
                             << " targ_parton_pz_ " << targ_parton_pz_;

    if (targ_parton_e_ >= 0.95 * P_A || proj_parton_e_ >= 0.95 * P_A ) {
        throw std::runtime_error(
            "Energy to subtract from 3DMCGlauber >= 0.95 * P_A "
            + std::to_string(0.95 * P_A)+". Turn on Verbose for more info.");
    }

}


std::vector<double> MCGlauberWrapper::Get_quarks_pos_proj_lab() {
    // get the x, y, z of the three valence quarks of colliding projectile
    // The fourth parton is the soft ball
    // 3DGlauber attributes the remaining energy and momentum carried by the
    // sea quarks and gluons to a soft gluon cloud
    // Output formulation is (x,y,z, x,y,z, x,y,z, x,y,z)
    mc_gen_->GetHardPos(hard_parton_t_, hard_parton_x_, hard_parton_y_,
                        hard_parton_z_);
    return(mc_gen_->GetQuarkPosProj());
}


std::vector<double> MCGlauberWrapper::Get_quarks_pos_targ_lab() {
    // get the x, y, z of the three valence quarks of colliding target
    // The fourth parton is the soft ball
    // 3DGlauber attributes the remaining energy and momentum carried by the
    // sea quarks and gluons to a soft gluon cloud
    // Output formulation is (x,y,z, x,y,z, x,y,z, x,y,z)
    mc_gen_->GetHardPos(hard_parton_t_, hard_parton_x_, hard_parton_y_,
                        hard_parton_z_);
    return(mc_gen_->GetQuarkPosTarg());
}


std::vector<double> MCGlauberWrapper::Get_remnant_proj() {
    // get the fout-momentum (E, px, py, pz) of the remnant in projectile
    return(mc_gen_->GetRemMom_Proj());
}


std::vector<double> MCGlauberWrapper::Get_remnant_targ() {
    // get the fout-momentum (E, px, py, pz) of the remnant in target
    return(mc_gen_->GetRemMom_Targ());
}


void MCGlauberWrapper::GetHardPartonPosAndMomentumProj() {
    mc_gen_->GetMomandPos_Proj(hard_parton_t_, hard_parton_x_, hard_parton_y_,
                               hard_parton_z_, proj_parton_e_, proj_parton_px_,
                               proj_parton_py_, proj_parton_pz_);
}


void MCGlauberWrapper::GetHardPartonPosAndMomentumTarg() {
    mc_gen_->GetMomandPos_Targ(hard_parton_t_, hard_parton_x_, hard_parton_y_,
                               hard_parton_z_, targ_parton_e_, targ_parton_px_,
                               targ_parton_py_, targ_parton_pz_);
}

void MCGlauberWrapper::GenerateStrings() {
    // generate strings from 3D Glauber for MUSIC
    QCDStringList_.clear();
    auto stringList = mc_gen_->generate_strings();
    for (auto string_i: stringList) {
        std::vector<double> string_temp;
        for (auto ii: string_i) {
            string_temp.push_back(static_cast<double>(ii));
        }
        QCDStringList_.push_back(string_temp);
    }
    JSINFO << "Produced " << QCDStringList_.size() << " strings.";
}
