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
# jhputschke/3dMCGlauber, a fork of chunshen1987/3dMCGlauber: the JETSCAPE commit
# X-SCAPE used before (71116fe) plus a cap on the draws for the target's hard-quark
# position, which could loop forever with ISR (~1 in 500-1300 O+O events);
# upstream PR chunshen1987/3dMCGlauber#7; branch zjet_LO_sub adds the strings with a
# hard process's energy taken out, for a jet leg's hydro (<LO_hard_subtraction>)
commitHash="3fdf20fb1ac00d9425b038f928fc10548a357161"
# download the code package
rm -fr $folderName
git clone https://github.com/jhputschke/3dMCGlauber.git --branch zjet_LO_sub $folderName
cd $folderName
git checkout $commitHash
./get_LHAPDF.sh
