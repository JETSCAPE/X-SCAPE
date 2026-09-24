#!/usr/bin/env bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
#
# Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
#
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

folderName="3dMCGlauber"
commitHash="96e659420ccadcfd721979a3dc4095c7289c0ed9"
# download the code package
rm -fr $folderName
git clone https://github.com/chunshen1987/3dMCGlauber.git --branch CSCAPE $folderName
cd $folderName
git checkout $commitHash
./get_LHAPDF.sh
