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
#include <iostream>
#include <stdio.h>
#include <sys/stat.h>
#include <vector>

#include "JetScapeLogger.h"
#include "MCGlauberGenStringWrapper.h"
#include "JetScapeSignalManager.h"
#include "JetScapeXML.h"

// Diagnostic counters defined in src/jet/iMATTER.cc (temporary).
extern int g_iMatter_nsplits;
extern int g_iMatter_skip_initial_budget;
extern int g_iMatter_skip_virtuality;
extern int g_iMatter_skip_final_budget;

// Register the module with the base class
RegisterJetScapeModule<MCGlauberGenStringWrapper> MCGlauberGenStringWrapper::reg("MCGlauberGenString");

MCGlauberGenStringWrapper::MCGlauberGenStringWrapper() {
    SetId("MCGlauberGenString");
    event_id_ = 0;
}

void MCGlauberGenStringWrapper::InitTask() {
    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini) {
        JSWARN << "No initialization module, "
               << "try: auto MCGlb = make_shared<MCGlauberWrapper>(); "
               << "jetscape->Add(MCGlb);";
    }
}


void MCGlauberGenStringWrapper::ExecuteTask() {
    VERBOSE(1) << "Run 3DMCGlauber second time to generate strings for MUSIC "
                     << "...";
    try {
        ini->GetHardPartonPosAndMomentumProj();
        ini->GetHardPartonPosAndMomentumTarg();

        // generate strings for MUSIC with wounded nucleons
        // for ISR setup (wound_nucleons=true) 
        ini->GenerateStrings(true, event_id_);

        // Unconditional (not gated by the JSINFO/info logger flag) so this
        // is always visible for cross-checking the average energy
        // subtraction between different JETSCAPE-framework builds.
        std::cout << "ENERGY_SUBTRACTED_EVENT "
                  << ini->GetHardPartonEnergySubtracted() << " GeV"
                  << std::endl;

        std::cout << "IMATTER_SHOWER_STATS_XSCAPE"
                  << " nsplits=" << g_iMatter_nsplits
                  << " skip_initial_budget=" << g_iMatter_skip_initial_budget
                  << " skip_virtuality=" << g_iMatter_skip_virtuality
                  << " skip_final_budget=" << g_iMatter_skip_final_budget
                  << std::endl;
        g_iMatter_nsplits = 0;
        g_iMatter_skip_initial_budget = 0;
        g_iMatter_skip_virtuality = 0;
        g_iMatter_skip_final_budget = 0;

        ini->ClearHardPartonMomentum();
    } catch (std::exception &err) {
        Jetscape::JSWARN << err.what();
        std::exit(-1);
    }
}

std::vector<std::vector<double>> MCGlauberGenStringWrapper::Get_remnant_proj() {
    return(ini->Get_remnant_proj());
}

std::vector<std::vector<double>> MCGlauberGenStringWrapper::Get_remnant_targ() {
    return(ini->Get_remnant_targ());
}

