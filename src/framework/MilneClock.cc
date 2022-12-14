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

#include "MilneClock.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"

#include <iostream>

namespace Jetscape {

MilneClock::MilneClock() {
    SetId("MilneClock");
    etaMax_ = 0.;
    tauMin_ = 0.;
    tauMax_ = 0.;
}

void MilneClock::Info() {
    JSINFO << "Clock Name: " << GetId();
    JSINFO << "tauMin = " << tauMin_ << " fm/c";
    JSINFO << "tauMax = " << tauMax_ << " fm/c";
}

void MilneClock::Transform(std::weak_ptr<MainClock> mainClock) {
    auto p = mainClock.lock();
    if (p) {
        currentModuleTime_ = p->GetCurrentTime();
        moduleDeltaT_ = p->GetDeltaT();
        if (currentModuleTime_ > 0.) {
            tauMax_ = currentModuleTime_;
            tauMin_ = p->GetCurrentTime()/cosh(etaMax_);
        }
    } else {
        JSWARN << "Trying to transform module clock with no main clock ... ";
        exit(-1);
    }
}

}
