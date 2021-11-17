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
        nev_ = GetXMLElementDouble({"IS", "MCGlauber","nEvents_MCGla"});
        parameter_list_.read_in_parameters_from_file("mcgluaber.input");
        int ran_seed = parameter_list_.get_seed();
        auto gamma_beta = parameter_list_.get_tau_form_fluct_gamma_beta();
        
        ran_gen_ptr_ = std::shared_ptr<RandomUtil::Random>(new RandomUtil::Random(ran_seed, 0.0, 1.0, gamma_beta));
                    
        mc_gen_=std::unique_ptr<MCGlb::EventGenerator>(new MCGlb::EventGenerator("mcgluaber.input", ran_seed));
        MCGlauber_ptr_ = std::unique_ptr<MCGlb::Glauber>(new MCGlb::Glauber(parameter_list_, ran_gen_ptr_));
}
void MCGlauberWrapper::Exec() {
    Clear();
    Jetscape::JSINFO << "Run MCGlauber to generate initial hard positions ...";
    try {
        int ievent=0;
        mc_gen_->generate_pre_events(nev_);
        
        std::vector<MCGlb::CollisionEvent> collisionEvents = mc_gen_->get_CollisionEventvector();
        ncoll_ = collisionEvents.size();
        
        while(ievent<ncoll_){
             // xvec[0],xvec[1],xvec[2] and xvec[3] are: t, x, y, z
             auto xvec = 
                    collisionEvents[ievent].get_collision_position();
             binary_collision_x_.push_back(xvec[1]);
             binary_collision_y_.push_back(xvec[2]);
             //std::cout<<"binary_collision_x_========= "<<binary_collision_x_[ievent]<<std::endl;
             //std::cout<<"binary_collision_y_========= "<<binary_collision_y_[ievent]<<std::endl;
             ievent++;
        }
        event_id_++;
    } catch (std::exception &err) {
        Jetscape::JSWARN << err.what();
        std::exit(-1);
    }

}


void MCGlauberWrapper::Clear() {
    Jetscape::JSINFO << "clear initial condition vectors";
}

void MCGlauberWrapper::Write(weak_ptr<JetScapeWriter> w) {}
