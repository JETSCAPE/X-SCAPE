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
# commitHash="3eae883ac5f39959235881c7023cf5d26e216684" # for xscape 1.2
commitHash="09e0a7e06544fbd8f34f10c2de9cdeacf1a0976f" # commit with adjusted hotspot selction AND different algorithm
# commitHash="b927277f728d194221b56ebbbe6f92b347887bee" # commit with just algoritm change
# download the code package
rm -fr $folderName
git clone https://github.com/chunshen1987/3dMCGlauber.git --branch JETSCAPE $folderName
cd $folderName
git checkout $commitHash
./get_LHAPDF.sh
