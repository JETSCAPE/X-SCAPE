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
// Test the generation of SMASH a nucleus initial condition
// -----------------------------------------------------------------------------

#include <iostream>
#include <time.h>
#include <chrono>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

#include "SMASHNucleusWrapper.h"

using namespace Jetscape;

// Forward declaration
void Show();

int main(int argc, char **argv) {
  clock_t t;
  t = clock();
  time_t start, end;
  time(&start);

  cout << endl;

  // DEBUG=true by default and REMARK=false
  // can be also set also via XML file (at least partially)
  JetScapeLogger::Instance()->SetInfo(true);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  //SetVerboseLevel (9 a lot of additional debug output ...)
  //If you want to suppress it: use SetVerboseLevel(0) or max  SetVerboseLevel(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);

  Show();

  auto jetscape = make_shared<JetScape>();
  const char *mainXMLName = "../config/jetscape_main.xml";
  const char *userXMLName = "../config/jetscape_user_SMASHNucleusTest.xml";

  jetscape->SetXMLMainFileName(mainXMLName);
  jetscape->SetXMLUserFileName(userXMLName);
  JetScapeXML::Instance()->OpenXMLMainFile(mainXMLName);
  JetScapeXML::Instance()->OpenXMLUserFile(userXMLName);

  auto smash_nucleus = make_shared<SMASHNucleusWrapper>();

  jetscape->Add(smash_nucleus);

  jetscape->Init();
  jetscape->Exec();

  // Get the hadrons from the last SMASH nucleus
  std::vector<Hadron> h_list = smash_nucleus->GetCurrentHadronList();

  // Write the last nucleus to file (the previous ones are overwritten if
  // there were any). The nucleus is usually extracted from the initial
  // condition event by event.
  std::ofstream outfile("SMASHNucleusTest.csv");
  outfile << "# pid,charge,t,x,y,z,E,p_x,p_y,p_z" << std::endl;
  for (const auto &hadron : h_list) {
    const FourVector hadron_r = hadron.x_in();
    const FourVector hadron_p = hadron.p_in();
    outfile << hadron.pid() << "," << hadron.charge() << "," << hadron_r.t()
            << "," << hadron_r.x() << "," << hadron_r.y() << "," << hadron_r.z()
            << "," << hadron_p.t() << "," << hadron_p.x() << "," << hadron_p.y()
            << "," << hadron_p.z() << std::endl;
  }
  outfile.close();

  // Check that the number of hadrons is 208
  if (h_list.size() != 208) {
    JSWARN << "Number of hadrons is not 208!";
  } else {
    JSINFO << "Number of hadrons is 208! This is the good result!";
  }

  jetscape->Finish();

  INFO_NICE << "Finished!";
  cout << endl;

  t = clock() - t;
  time(&end);
  printf("CPU time: %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
  printf("Real time: %f seconds.\n", difftime(end, start));
  return 0;
}

void Show() {
  INFO_NICE
      << "---------------------------------------------------------------";
  INFO_NICE
      << "| SMASH Nucleus Initial Condition Test X-SCAPE Framework ...  |";
  INFO_NICE
      << "---------------------------------------------------------------";
  INFO_NICE;
}