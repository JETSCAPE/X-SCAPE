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
#include "MCGlauberGenStringWrapper.h"
#include "JetScapeSignalManager.h"
#include "JetScapeXML.h"

// Register the module with the base class
RegisterJetScapeModule<MCGlauberGenStringWrapper> MCGlauberGenStringWrapper::reg("MCGlauberGenString");

MCGlauberGenStringWrapper::MCGlauberGenStringWrapper() {
    SetId("MCGlauberGenString");
    event_id_ = 0;
}

void MCGlauberGenStringWrapper::Init() {
    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini) {
        JSWARN << "No initialization module, "
               << "try: auto MCGlb = make_shared<MCGlauberWrapper>(); "
               << "jetscape->Add(MCGlb);";
    }
}


void MCGlauberGenStringWrapper::Exec() {
    Jetscape::JSINFO << "Run 3DMCGlauber second time to generate strings for MUSIC "
                     << "...";
    try {
        ini->GetHardPartonPosAndMomentumProj();
        ini->GetHardPartonPosAndMomentumTarg();
        ini->GenerateStrings();// generate strings for the MUSIC
    } catch (std::exception &err) {
        Jetscape::JSWARN << err.what();
        std::exit(-1);
    }
}

