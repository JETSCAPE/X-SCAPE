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

# using a commit from the MUSIC repository that is compatible with the current X-SCAPE version
folderName="music4gpu"

git clone https://github.com/jhputschke/MUSIC4GPU.git -b XSCAPE $folderName
cd $folderName
cd EOS
bash download_hotQCD.sh SMASH_binary

# Finite muB not supprted on GPU yet!!!
# Download the 4D EoS tables (only needed for EOS 20, only download if necessary, large files)
# bash download_Neos4D.sh UrQMD