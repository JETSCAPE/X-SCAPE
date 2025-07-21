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

# using a commit from the iSS repository that is compatible with the current X-SCAPE version
folderName="iSS"
commitHash="b00ee76357105030b8586c5e7a14d86a620bbe0c"

git clone https://github.com/chunshen1987/iSS -b XSCAPE iSS
cd $folderName
git checkout $commitHash

# Additional tables needed for 4D EoS (download only if necessary, large files)
cd iSS_tables/EOS_tables
bash download_HRG4D.sh

cd ../deltaf_tables/urqmd
bash download_NEoS4D_deltafCoeffs.sh