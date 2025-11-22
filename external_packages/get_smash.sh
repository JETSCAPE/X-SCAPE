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

# If running on macOS, try to use Homebrew's Eigen installation for CMake
EIGEN_OPTION=""
if [[ "$(uname)" == "Darwin" ]]; then
    if command -v brew >/dev/null 2>&1; then
        PREFIX="$(brew --prefix eigen 2>/dev/null || true)"
        if [[ -n "${PREFIX}" ]]; then
            EIGEN_OPTION="-DEigen3_DIR=${PREFIX}/share/eigen3/cmake"
        fi
    fi
fi

cmake .. ${EIGEN_OPTION} -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
num_cores=${1:-1}
echo "Compiling SMASH using ${num_cores} cores."
make -j${num_cores} smash_shared
)
export SMASH_DIR="$(pwd)/smash/smash_code"
