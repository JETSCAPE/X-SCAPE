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
commitHash="444410c34f8159fc126a2b3f68f8f78e3ee5cbad"

git clone https://github.com/chunshen1987/iSS -b JETSCAPE iSS
cd $folderName
git checkout $commitHash
