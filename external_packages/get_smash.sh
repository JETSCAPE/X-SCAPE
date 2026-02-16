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
git clone --depth=1 https://github.com/smash-transport/smash.git --branch SMASH-3.2.2 smash/smash_code

# 2) Compile SMASH
(
cd smash/smash_code
mkdir build
cd build
cmake_args=(
	-DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
)

if [[ -n "${Eigen3_DIR}" ]]; then
	cmake_args+=( -DEigen3_DIR=${Eigen3_DIR} )
fi

if [[ -n "${EIGEN3_INCLUDE_DIR}" ]]; then
	cmake_args+=( -DEIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR} )
fi

cmake .. "${cmake_args[@]}"
num_cores=${1:-1}
echo "Compiling SMASH using ${num_cores} cores."
make -j${num_cores} smash_shared
)
export SMASH_DIR="$(pwd)/smash/smash_code"
