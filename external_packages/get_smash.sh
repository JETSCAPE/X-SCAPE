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
mkdir -p build
cd build

# Ensure Eigen path is set
if [ -z "${EIGEN3_ROOT}" ]; then
	if command -v brew >/dev/null 2>&1; then
		export EIGEN3_ROOT="$(brew --prefix eigen)/include/eigen3"
	fi
fi

# Locate HepMC3 and Eigen. Add both if present.
prefix_path="${HEPMC3_DIR};${EIGEN3_ROOT}"

cmake .. \
	-DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config \
	-DCMAKE_PREFIX_PATH="${prefix_path}" \
	-DEIGEN3_ROOT="${EIGEN3_ROOT}" || exit 1

num_cores=${1:-1}
echo "Compiling SMASH using ${num_cores} cores."

if make -n | grep -q "smash_shared"; then
	make -j${num_cores} smash_shared || exit 1
else
	make -j${num_cores} || exit 1
fi
)
export SMASH_DIR="$(pwd)/smash/smash_code"
