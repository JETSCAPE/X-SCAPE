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
commitHash="b74d26e399d27f8ad144375d204fdfbd36336608" # for xscape 1.2
# download the code package
rm -fr $folderName
git clone https://github.com/chunshen1987/3dMCGlauber.git --branch JETSCAPE $folderName
cd $folderName
git checkout $commitHash
./get_LHAPDF.sh
