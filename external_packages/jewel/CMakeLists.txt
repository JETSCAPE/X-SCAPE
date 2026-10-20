# JEWEL 2.6.0 (C++ port) as a library for the JETSCAPE module "Jewel".
#
# The JEWEL core is C++ (built here with C++2b); the modified PYTHIA 6.4.25 that
# JEWEL uses is Fortran. If a Fortran compiler is available it is built here,
# otherwise a prebuilt static library is used (default: lib/libpythia6jewel.a,
# made by build_pythia6.sh with any gfortran, e.g. outside the JETSCAPE docker
# image, which has no gfortran). The original stand-alone build of the port is
# CMakeLists.standalone.txt.

set(JEWEL_PYTHIA6_LIB "${CMAKE_CURRENT_SOURCE_DIR}/lib/libpythia6jewel.a"
    CACHE FILEPATH "prebuilt PYTHIA 6 for JEWEL (used without a Fortran compiler)")

include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER AND NOT JEWEL_USE_PREBUILT_PYTHIA6)
  enable_language(Fortran)
  add_library(pythia6jewel STATIC fortran/pythia6425mod-lhapdf6.f fortran/fio.f90)
  set_source_files_properties(fortran/pythia6425mod-lhapdf6.f PROPERTIES COMPILE_OPTIONS "-w")
  set_target_properties(pythia6jewel PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_compile_options(pythia6jewel PRIVATE -O2 -ffp-contract=off)
  set(JEWEL_PYTHIA6 pythia6jewel)
  message(STATUS "JEWEL: building PYTHIA 6 with ${CMAKE_Fortran_COMPILER}")
else()
  if(NOT EXISTS "${JEWEL_PYTHIA6_LIB}")
    message(FATAL_ERROR
      "JEWEL: no Fortran compiler and no prebuilt PYTHIA 6 at ${JEWEL_PYTHIA6_LIB}.\n"
      "Run external_packages/jewel/build_pythia6.sh where gfortran is available.")
  endif()
  set(JEWEL_PYTHIA6 "${JEWEL_PYTHIA6_LIB}")
  message(STATUS "JEWEL: using prebuilt PYTHIA 6 ${JEWEL_PYTHIA6_LIB}")
endif()

find_library(JEWEL_GFORTRAN_LIB NAMES gfortran libgfortran.so.5
             HINTS /usr/lib64 /usr/lib /usr/lib/x86_64-linux-gnu)
find_library(JEWEL_LHAPDF_LIB NAMES LHAPDF HINTS ENV LHAPDF_DIR PATH_SUFFIXES lib lib64)
if(NOT JEWEL_GFORTRAN_LIB)
  message(FATAL_ERROR "JEWEL: libgfortran not found")
endif()
if(NOT JEWEL_LHAPDF_LIB)
  message(FATAL_ERROR "JEWEL: LHAPDF not found (PYTHIA 6 is linked against its Fortran interface)")
endif()

add_library(jewel_jetscape STATIC
  src/io.cpp
  src/pythia6.cpp
  src/main_init.cpp
  src/event_setup.cpp
  src/cascade.cpp
  src/scattering.cpp
  src/kinematics.cpp
  src/sudakov_tables.cpp
  src/utilities_hepmc.cpp
  src/eix.cpp
  src/subtraction.cpp
  src/medium.cpp
  src/medium_simple.cpp
  src/medium_vac.cpp
  src/hfdiffusion.cpp
  src/jetscape_api.cpp)
target_include_directories(jewel_jetscape PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(jewel_jetscape PROPERTIES
  CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON LINKER_LANGUAGE CXX)
# same floating-point treatment as the stand-alone port (see PORTING.md); JETSCAPE's global
# warning flags are for its own code
target_compile_options(jewel_jetscape PRIVATE -O2 -ffp-contract=off -fno-builtin-sin
                       -fno-builtin-cos -w)
target_link_libraries(jewel_jetscape PUBLIC ${JEWEL_PYTHIA6} ${JEWEL_LHAPDF_LIB}
                      ${JEWEL_GFORTRAN_LIB})

# check of the interface without JETSCAPE (see the comment in the source)
add_executable(test_jetscape_brick tests/test_jetscape_brick.cpp)
target_link_libraries(test_jetscape_brick PRIVATE jewel_jetscape)
set_target_properties(test_jetscape_brick PROPERTIES CXX_STANDARD 17 LINKER_LANGUAGE CXX)
