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

# 1) Download the SMASH code
# TODO: Change branch back to main SMASH version before public XSCAPE version
git clone --depth=1 https://github.com/smash-transport/smash.git --branch staudenmaier/mod_plist_ts smash/smash_code

# 2) Compile SMASH
cd smash/smash_code
mkdir build
cd build
cmake .. -DCMAKE_CXX_STANDARD=14 -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
number_of_cores=`nproc --all`
number_of_cores_to_compile=$(( ${number_of_cores} > 20 ? 20 : ${number_of_cores} ))
echo "Compiling SMASH using ${number_of_cores_to_compile} cores."
make -j${number_of_cores} smash_shared
