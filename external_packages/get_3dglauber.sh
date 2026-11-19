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
# upstream PR chunshen1987/3dMCGlauber#7
commitHash="b99d680b2f2d44694053232804386d3d60474c04"
# download the code package
rm -fr $folderName
git clone https://github.com/jhputschke/3dMCGlauber.git --branch cap_hard_quark_draws $folderName
cd $folderName
git checkout $commitHash
./get_LHAPDF.sh
