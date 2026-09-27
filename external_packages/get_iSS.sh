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
# (jhputschke/iSS common_seeds: chunshen1987/iSS XSCAPE d242555 + the FSSW yield-loop
# speed-up + correlated sampling, until they are merged upstream)
folderName="iSS"
commitHash="3192982552409ca8bfcc54aa90975571e0490db8"

git clone https://github.com/jhputschke/iSS -b common_seeds iSS
cd $folderName
git checkout $commitHash

# Additional tables needed for 4D EoS (download only if necessary, large files)
cd iSS_tables/EOS_tables
bash download_HRG4D.sh

cd ../deltaf_tables/urqmd
bash download_NEoS4D_deltafCoeffs.sh