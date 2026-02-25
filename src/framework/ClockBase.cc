#include "ClockBase.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"

#include <iostream>

namespace Jetscape {

// bool ClockBase::use_clock = true;

/**
 * @file ClockBase.cc
 * @brief Implementation of the ClockBase methods.
 */

/**
 * @brief Default constructor initializes id and time reference 
 * frame id to empty strings.
 */
ClockBase::ClockBase() { id = ""; }

/**
 * @brief Log basic information about the clock.
 *
 * This implementation writes the clock id and the time reference
 * frame id to the framework logger at the info level. Derived
 * classes may override `Info()` to include additional fields.
 */
void ClockBase::Info() { JSINFO << GetId() << " " << GetTimeRefFrameId(); }

}  // namespace Jetscape