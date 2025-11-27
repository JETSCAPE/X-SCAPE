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

# Inject patched FindEigen3.cmake to handle Eigen >=5 before configuring
cat > cmake/FindEigen3.cmake <<'EOF'
# Patched FindEigen3.cmake injected by JETSCAPE get_smash.sh
if(NOT Eigen3_FIND_VERSION)
	if(NOT Eigen3_FIND_VERSION_MAJOR)
		set(Eigen3_FIND_VERSION_MAJOR 2)
	endif()
	if(NOT Eigen3_FIND_VERSION_MINOR)
		set(Eigen3_FIND_VERSION_MINOR 91)
	endif()
	if(NOT Eigen3_FIND_VERSION_PATCH)
		set(Eigen3_FIND_VERSION_PATCH 0)
	endif()
	set(Eigen3_FIND_VERSION "${Eigen3_FIND_VERSION_MAJOR}.${Eigen3_FIND_VERSION_MINOR}.${Eigen3_FIND_VERSION_PATCH}")
endif()

macro(_eigen3_check_version)
	set(_eigen_macro_header "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h")
	if(NOT EXISTS "${_eigen_macro_header}")
		file(GLOB _eigen_macro_candidates "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/*Macros*.h")
		list(LENGTH _eigen_macro_candidates _cand_len)
		if(_cand_len GREATER 0)
			list(GET _eigen_macro_candidates 0 _eigen_macro_header)
		elseif(EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Version.h")
			set(_eigen_macro_header "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Version.h")
		endif()
	endif()
	if(NOT EXISTS "${_eigen_macro_header}")
		set(EIGEN3_VERSION "UNKNOWN")
		set(EIGEN3_VERSION_OK TRUE)
		message(STATUS "Eigen3 macro header not found; assuming version >= ${Eigen3_FIND_VERSION}.")
		return()
	endif()
	file(READ "${_eigen_macro_header}" _eigen3_version_header)
	string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _wv "${_eigen3_version_header}")
	set(EIGEN3_WORLD_VERSION "${CMAKE_MATCH_1}")
	string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _mv "${_eigen3_version_header}")
	set(EIGEN3_MAJOR_VERSION "${CMAKE_MATCH_1}")
	string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _nv "${_eigen3_version_header}")
	set(EIGEN3_MINOR_VERSION "${CMAKE_MATCH_1}")
	if(EIGEN3_WORLD_VERSION STREQUAL "" OR EIGEN3_MAJOR_VERSION STREQUAL "" OR EIGEN3_MINOR_VERSION STREQUAL "")
		string(REGEX MATCH "define[ \t]+EIGEN_VERSION_STRING[ \t]+\"([0-9]+\.[0-9]+\.[0-9]+)\"" _vs "${_eigen3_version_header}")
		if(CMAKE_MATCH_1)
			set(EIGEN3_VERSION "${CMAKE_MATCH_1}")
		else()
			set(EIGEN3_VERSION "UNKNOWN")
			set(EIGEN3_VERSION_OK TRUE)
			message(STATUS "Eigen3 version macros missing; assuming version >= ${Eigen3_FIND_VERSION}.")
			return()
		endif()
	else()
		set(EIGEN3_VERSION "${EIGEN3_WORLD_VERSION}.${EIGEN3_MAJOR_VERSION}.${EIGEN3_MINOR_VERSION}")
	endif()
	if(${EIGEN3_VERSION} VERSION_LESS ${Eigen3_FIND_VERSION})
		set(EIGEN3_VERSION_OK FALSE)
	else()
		set(EIGEN3_VERSION_OK TRUE)
	endif()
	if(NOT EIGEN3_VERSION_OK)
		message(STATUS "Eigen3 version ${EIGEN3_VERSION} found in ${EIGEN3_INCLUDE_DIR}, but at least version ${Eigen3_FIND_VERSION} is required")
	endif()
endmacro()

if(EIGEN3_INCLUDE_DIR)
	_eigen3_check_version()
	set(EIGEN3_FOUND ${EIGEN3_VERSION_OK})
else()
	find_path(EIGEN3_INCLUDE_DIR
		NAMES signature_of_eigen3_matrix_library
		HINTS ENV EIGEN3_ROOT ENV EIGEN3_ROOT_DIR
		PATHS ${CMAKE_INSTALL_PREFIX}/include
		PATH_SUFFIXES eigen3 eigen)
	if(EIGEN3_INCLUDE_DIR)
		_eigen3_check_version()
	endif()
	include(FindPackageHandleStandardArgs)
	find_package_handle_standard_args(Eigen3 DEFAULT_MSG EIGEN3_INCLUDE_DIR EIGEN3_VERSION_OK)
	mark_as_advanced(EIGEN3_INCLUDE_DIR)
endif()
EOF

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
	echo "Target 'smash_shared' not found; building default targets instead." 
	make -j${num_cores} || exit 1
fi
)
export SMASH_DIR="$(pwd)/smash/smash_code"
