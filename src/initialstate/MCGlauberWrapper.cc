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


MCGlauberWrapper::~MCGlauberWrapper() {}

                                
void MCGlauberWrapper::InitTask() {
        parameter_list_.read_in_parameters_from_file("mcgluaber.input");
        //int ran_seed = parameter_list_.get_seed();
        auto ran_seed = (*GetMt19937Generator())();
        auto gamma_beta = parameter_list_.get_tau_form_fluct_gamma_beta();
        
        ran_gen_ptr_ = std::shared_ptr<RandomUtil::Random>(new RandomUtil::Random(ran_seed, 0.0, 1.0, gamma_beta));
                    
        mc_gen_=std::unique_ptr<MCGlb::EventGenerator>(new MCGlb::EventGenerator("mcgluaber.input", ran_seed));
        MCGlauber_ptr_ = std::unique_ptr<MCGlb::Glauber>(new MCGlb::Glauber(parameter_list_, ran_gen_ptr_));
}

void MCGlauberWrapper::Clear() {
    Jetscape::JSINFO << "clear initial condition vectors";
    binary_collision_t_.clear();
    binary_collision_x_.clear();
    binary_collision_y_.clear();
    binary_collision_z_.clear();
}


void MCGlauberWrapper::Exec() {
    Clear();
    Jetscape::JSINFO << "Run MCGlauber to generate initial hard positions ...";
    try {
        int iparticle=0;
        mc_gen_->generate_pre_events();
        
        std::vector<MCGlb::CollisionEvent> collisionEvents = mc_gen_->get_CollisionEventvector();
        ncoll_ = collisionEvents.size();
        rand_int_ptr_ = (
               std::make_shared<std::uniform_int_distribution<int>>(0, ncoll_-1));
        while(iparticle<ncoll_){
             // xvec[0],xvec[1],xvec[2] and xvec[3] are: t, x, y, z
             auto xvec = 
                       collisionEvents[iparticle].get_collision_position();
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



void MCGlauberWrapper::SampleABinaryCollisionPoint(double &t, double &x, double &y, double &z) {
    int rand_idx = (*rand_int_ptr_)(*GetMt19937Generator());
    t = binary_collision_t_[rand_idx];
    x = binary_collision_x_[rand_idx];
    y = binary_collision_y_[rand_idx];
    z = binary_collision_z_[rand_idx];
}



void MCGlauberWrapper::Write(weak_ptr<JetScapeWriter> w) {}
