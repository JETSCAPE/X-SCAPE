/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

/**
 * @brief Displays the X-SCAPE banner with ASCII art and license information.
 *
 * This function prints a stylized ASCII banner with the X-SCAPE logo and key
 * details about the framework, including citation information and the licensing
 * terms. The output is handled using the `INFO_NICE` logging macro.
 *
 * @note The function does not take any parameters or return values. It only
 * prints information to the standard output using the logger.
 */
void ShowXscapeBanner() {
  INFO_NICE
      << "*--------------------------------------------------------------*";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "|                 %     /                                      |";
  INFO_NICE
      << "|                  %   /                                       |";
  INFO_NICE << "|                   % /             /" << (char)92
            << "                         |";
  INFO_NICE << "|                    X           /" << (char)92 << "/ "
            << " " << (char)92 << "                        |";
  INFO_NICE << "|                   / %       /" << (char)92 << "/"
            << "   |  " << (char)92 << "/" << (char)92
            << "                     |";
  INFO_NICE << "|                  /   %     /    % | %   " << (char)92
            << "                    |";
  INFO_NICE << "|               __/     %___/"
            << "      %|%     " << (char)92 << "/" << (char)92
            << "__               |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "|                     XSCAPE by JETSCAPE                       |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "|           X-Ion Collisions with a Statistically              |";
  INFO_NICE
      << "|       and Computationally Advanced Program Envelope          |";
  INFO_NICE
      << "|                     http://jetscape.org                      |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "| Please cite arXiv:1903.07706 if you use this package for     |";
  INFO_NICE
      << "| scientific work.                                             |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "| JETSCAPE is provided without warranty under the terms        |";
  INFO_NICE
      << "| of the GNU GPLv3. It uses xxx code(s).                       |";
  INFO_NICE
      << "| See COPYING file for details.                                |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "*--------------------------------------------------------------*";
  INFO_NICE << " XSCAPE version = " << XscapeVersion
            << " (includes JETSCAPE version = " << JetScapeVersion << ")";
  INFO_NICE;
  // INFO_NICE <<
  // "*--------------------------------------------------------------*";
}

void ShowXscapeBanner2() {
  // ── box: ║ + 62 chars + ║ = 64 wide ──────────────────────────────────────
  // XSCAPE block layout per row (total 50 cols, centred with 6-col padding):
  //   X(9) + S(9) + C(8) + A(8) + P(8) + E(8) = 50
  INFO_NICE << "╔══════════════════════════════════════════════════════════════╗";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║      ██╗  ██╗  ██████╗  ██████╗ █████╗ ██████╗ ███████╗      ║";
  INFO_NICE << "║      ╚██╗██╔╝ ██╔════╝ ██╔════╝██╔══██╗██╔══██╗██╔════╝      ║";
  INFO_NICE << "║       ╚████╔╝ ╚█████╗  ██║     ███████║██████╔╝█████╗        ║";
  INFO_NICE << "║       ██╔═██╗  ╚════██╗██║     ██╔══██║██╔═══╝ ██╔══╝        ║";
  INFO_NICE << "║      ██╔╝  ██╗ ██████╔╝╚██████╗██║  ██║██║     ███████╗      ║";
  INFO_NICE << "║      ╚═╝   ╚═╝ ╚═════╝  ╚═════╝╚═╝  ╚═╝╚═╝     ╚══════╝      ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║                    X-SCAPE  by  JETSCAPE                     ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║         X-Ion Collisions with a Statistically and            ║";
  INFO_NICE << "║       Computationally Advanced Program Envelope              ║";
  INFO_NICE << "║                     http://jetscape.org                      ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║  Please cite arXiv:1903.07706 if you use this package for    ║";
  INFO_NICE << "║  scientific work.                                            ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║  JETSCAPE is provided without warranty under the terms of    ║";
  INFO_NICE << "║  the GNU GPLv3. See COPYING file for details.                ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "╚══════════════════════════════════════════════════════════════╝";
  INFO_NICE << " XSCAPE version = " << XscapeVersion
            << " (includes JETSCAPE version = " << JetScapeVersion << ")";
  INFO_NICE;
}
/**
 * @brief Modernized X-SCAPE banner.
 *
 * Same mountain/back-to-back-jet landscape as ShowXscapeBanner(), but the
 * percent-sign X arms are replaced by clean geometric \ / characters, and
 * the XSCAPE doom-font block logo is added below the art as a second visual.
 * The box uses Unicode double-line characters; inner width is 62 columns.
 */
void ShowXscapeBannerModern() {
  // ── landscape art: \ / X with right-side jet mountain ────────────────────
  INFO_NICE << "╔══════════════════════════════════════════════════════════════╗";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║                 \\     /                                      ║";
  INFO_NICE << "║                  \\   /                                       ║";
  INFO_NICE << "║                   \\ /             /" << (char)92 << "                         ║";
  INFO_NICE << "║                    X           /" << (char)92 << "/ "
            << " " << (char)92 << "                        ║";
  INFO_NICE << "║                   / \\       /" << (char)92 << "/"
            << "   |  " << (char)92 << "/" << (char)92
            << "                     ║";
  INFO_NICE << "║                  /   \\     /    % | %   " << (char)92
            << "                    ║";
  INFO_NICE << "║               __/     \\___/"
            << "      %|%     " << (char)92 << "/" << (char)92
            << "__               ║";
  INFO_NICE << "║                                                              ║";
  // ── XSCAPE doom-font block logo (X=9, S=9, C=8, A=8, P=8, E=8 → 50 cols) ─
  INFO_NICE << "║      ██╗  ██╗  ██████╗  ██████╗ █████╗ ██████╗ ███████╗      ║";
  INFO_NICE << "║      ╚██╗██╔╝ ██╔════╝ ██╔════╝██╔══██╗██╔══██╗██╔════╝      ║";
  INFO_NICE << "║       ╚████╔╝ ╚█████╗  ██║     ███████║██████╔╝█████╗        ║";
  INFO_NICE << "║       ██╔═██╗  ╚════██╗██║     ██╔══██║██╔═══╝ ██╔══╝        ║";
  INFO_NICE << "║      ██╔╝  ██╗ ██████╔╝╚██████╗██║  ██║██║     ███████╗      ║";
  INFO_NICE << "║      ╚═╝   ╚═╝ ╚═════╝  ╚═════╝╚═╝  ╚═╝╚═╝     ╚══════╝      ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║                         by  JETSCAPE                         ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║           X-Ion Collisions with a Statistically              ║";
  INFO_NICE << "║       and Computationally Advanced Program Envelope          ║";
  INFO_NICE << "║                     http://jetscape.org                      ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║         Please cite arXiv:1903.07706 if you use this         ║";
  INFO_NICE << "║                package for scientific work.                  ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "║        JETSCAPE is provided without warranty under the       ║";
  INFO_NICE << "║          terms of the GNU GPLv3. It uses xxx code(s).        ║";
  INFO_NICE << "║                 See COPYING file for details.                ║";
  INFO_NICE << "║                                                              ║";
  INFO_NICE << "╚══════════════════════════════════════════════════════════════╝";
  INFO_NICE << " XSCAPE version = " << XscapeVersion
            << " (includes JETSCAPE version = " << JetScapeVersion << ")";
  INFO_NICE;
}

}  // end namespace Jetscape

#endif
