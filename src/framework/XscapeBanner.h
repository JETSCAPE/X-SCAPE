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

#ifndef JETSCAPE_BANNER_H
#define JETSCAPE_BANNER_H

#include "JetScapeLogger.h"
#include "Version.h"

namespace Jetscape {

void ShowXscapeBanner() {
  INFO_NICE << "*--------------------------------------------------------------*";
  INFO_NICE << "|                                                              |";
  INFO_NICE << "|                                                              |";
  INFO_NICE << "|                 %     /                                      |";
  INFO_NICE << "|                  %   /                                       |";
  INFO_NICE << "|                   % /             /"<< (char)92 << "                         |";
  INFO_NICE << "|                    X           /"<< (char)92 << "/ "<< " " << (char)92 << "                        |";
  INFO_NICE << "|                   / %       /" << (char)92 << "/"  << "   |  " << (char)92 << "/" << (char)92<< "                     |";
  INFO_NICE << "|                  /   %     /    % | %   " << (char)92 << "                    |";
  INFO_NICE << "|               __/     %___/"<< "      %|%     " << (char)92 << "/" << (char)92<< "__               |";
  INFO_NICE << "|                                                              |";
  INFO_NICE << "|                     XSCAPE by JETSCAPE                       |";
  INFO_NICE << "|                                                              |";
  INFO_NICE << "|           X-Ion Collisions with a Statistically              |";
  INFO_NICE << "|       and Computationally Advanced Program Envelope          |";
  INFO_NICE << "|                     http://jetscape.org                      |";
  INFO_NICE << "|                                                              |";
  INFO_NICE << "| Please cite arXiv:1903.07706 if you use this package for     |";
  INFO_NICE << "| scientific work.                                             |";
  INFO_NICE << "|                                                              |";
  INFO_NICE << "| JETSCAPE is provided without warranty under the terms        |";
  INFO_NICE << "| of the GNU GPLv3. It uses xxx code(s).                       |";
  INFO_NICE << "| See COPYING file for details.                                |";
  INFO_NICE << "|                                                              |";
  INFO_NICE << "*--------------------------------------------------------------*";
  INFO_NICE <<" XSCAPE version = "<<XscapeVersion<<" (includes JETSCAPE version = "<<JetScapeVersion<<")";
  INFO_NICE;
  //INFO_NICE << "*--------------------------------------------------------------*";
}

} // end namespace Jetscape

#endif
