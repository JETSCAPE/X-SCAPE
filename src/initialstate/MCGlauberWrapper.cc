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
#include "MCGlauberWrapper.h"


// Register the module with the base class
RegisterJetScapeModule<MCGlauberWrapper> MCGlauberWrapper::reg("MCGlauber");

MCGlauberWrapper::MCGlauberWrapper() {
    SetId("MCGlauber");
    event_id_ = 0;
}


void MCGlauberWrapper::InitTask() {
        parameter_list_.read_in_parameters_from_file("mcgluaber.input");
        //int ran_seed = parameter_list_.get_seed();
        auto ran_seed = (*GetMt19937Generator())();
        auto gamma_beta = parameter_list_.get_tau_form_fluct_gamma_beta();

        mc_gen_ = std::shared_ptr<MCGlb::EventGenerator>(
                  new MCGlb::EventGenerator("mcgluaber.input", ran_seed));
}

void MCGlauberWrapper::Clear() {
    Jetscape::JSINFO << "clear initial condition vectors";
    binary_collision_t_.clear();
    binary_collision_x_.clear();
    binary_collision_y_.clear();
    binary_collision_z_.clear();
    HardPartonPosAndMomProj_.clear();
    HardPartonPosAndMomTarg_.clear();
}


void MCGlauberWrapper::Exec() {
    Clear();
    Jetscape::JSINFO << "Run 3DMCGlauber to generate initial hard positions "
                     << "...";
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

}


void MCGlauberWrapper::SampleABinaryCollisionPoint(
        double &t, double &x, double &y, double &z) {
    int rand_idx = (*rand_int_ptr_)(*GetMt19937Generator());
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

void MCGlauberWrapper::OutputHardPartonMomentum(double E, double px, double py, double pz,
                                                int direction) {
    if (direction == 1) {
        proj_parton_e_  = E;
        proj_parton_px_ = px;
        proj_parton_py_ = py;
        proj_parton_pz_ = pz;
    } else {
        targ_parton_e_  = E;
        targ_parton_px_ = px;
        targ_parton_py_ = py;
        targ_parton_pz_ = pz;
    }
}


std::vector<double> MCGlauberWrapper::GetHardPartonPosAndMomentumProj() {
    HardPartonPosAndMomProj_.push_back(hard_parton_t_);
    HardPartonPosAndMomProj_.push_back(hard_parton_x_);
    HardPartonPosAndMomProj_.push_back(hard_parton_y_);
    HardPartonPosAndMomProj_.push_back(hard_parton_z_);
    HardPartonPosAndMomProj_.push_back(proj_parton_e_);
    HardPartonPosAndMomProj_.push_back(proj_parton_px_);
    HardPartonPosAndMomProj_.push_back(proj_parton_py_);
    HardPartonPosAndMomProj_.push_back(proj_parton_pz_);
    return HardPartonPosAndMomProj_;
}

std::vector<double> MCGlauberWrapper::GetHardPartonPosAndMomentumTarg() {
    HardPartonPosAndMomTarg_.push_back(hard_parton_t_);
    HardPartonPosAndMomTarg_.push_back(hard_parton_x_);
    HardPartonPosAndMomTarg_.push_back(hard_parton_y_);
    HardPartonPosAndMomTarg_.push_back(hard_parton_z_);
    HardPartonPosAndMomTarg_.push_back(targ_parton_e_);
    HardPartonPosAndMomTarg_.push_back(targ_parton_px_);
    HardPartonPosAndMomTarg_.push_back(targ_parton_py_);
    HardPartonPosAndMomTarg_.push_back(targ_parton_pz_);
    return HardPartonPosAndMomTarg_;
}

void MCGlauberWrapper::GenerateStrings() {
    mc_gen_->generate_strings(); // generate 3D Glauber for MUSIC
}
