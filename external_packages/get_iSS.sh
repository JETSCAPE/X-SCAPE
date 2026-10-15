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
# (jhputschke/iSS smash_ready: chunshen1987/iSS XSCAPE d242555 + the FSSW yield-loop
# speed-ups + correlated sampling (iSS_optim) + decays of the SMASH hadron list and the
# undecayed primaries (decay_smash_list, keep_primaries), until they are merged upstream)
folderName="iSS"
commitHash="6ac0333125f2a8521120eb4580e828cde428166e"

git clone https://github.com/jhputschke/iSS -b smash_ready iSS
cd $folderName
git checkout $commitHash

# Additional tables needed for 4D EoS (download only if necessary, large files)
cd iSS_tables/EOS_tables
bash download_HRG4D.sh

cd ../deltaf_tables/urqmd
bash download_NEoS4D_deltafCoeffs.sh