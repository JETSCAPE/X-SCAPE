# X-SCAPE Installation Guide

## Contents

1. [Prerequisites](#1-prerequisites)
2. [Running with Docker / Singularity containers](#2-running-with-docker--singularity-containers)
3. [Obtaining the source](#3-obtaining-the-source)
4. [Downloading external packages](#4-downloading-external-packages)
5. [CMake options reference](#5-cmake-options-reference)
6. [Build configurations (quick recipes)](#6-build-configurations-quick-recipes)
7. [Installing to a prefix](#7-installing-to-a-prefix)
8. [Running from the build directory](#8-running-from-the-build-directory)
9. [Parallel runs with the launcher](#9-parallel-runs-with-the-launcher)
10. [Environment variables](#10-environment-variables)
11. [Data-asset search order](#11-data-asset-search-order)
12. [Interactive configurator](#12-interactive-configurator)
13. [Troubleshooting](#13-troubleshooting)
14. [Using X-SCAPE in a downstream CMake project](#14-using-x-scape-in-a-downstream-cmake-project)

---

## 1. Prerequisites

### Required (always)

| Package | Minimum version | Notes |
|---------|----------------|-------|
| CMake   | 3.5            | 3.14+ recommended |
| C++17 compiler | GCC 8 / Clang 9 / AppleClang 12 | |
| Boost   | 1.67           | components: `filesystem`, `program_options`, `system` |
| Pythia8 | 8.3            | set `PYTHIA8` or `PYTHIA8DIR` env var if not in PATH |
| HDF5    | 1.8            | C and C++ bindings required |
| zlib    | any            | usually present on all platforms |

### Optional (auto-detected)

| Package | Required by | Notes |
|---------|-------------|-------|
| GSL     | MUSIC, IP-Glasma, 3dMCGlauber, freestream-milne | Most physics options need it |
| HepMC3  | HepMC writer | auto-detected; set `HEPMC_DIR` if non-standard |
| ROOT    | ROOT output writer | set `ROOTSYS`; enable with `-DUSE_ROOT=ON` |
| MPI     | SMASH afterburner | auto-detected |
| OpenMP  | parallel loops | auto-detected; on macOS install via `brew install libomp` |
| OpenCL  | CLVisc GPU hydro | auto-detected |

### macOS notes

```bash
brew install boost pythia8 hdf5 gsl libomp cmake
```

`libomp` is found automatically via `brew --prefix libomp`; no manual path is needed.

---

## 2. Running with Docker / Singularity containers

The easiest way to get a complete, pre-built prerequisite environment without installing anything locally is to use an official JETSCAPE container image. Full step-by-step instructions, image tags, and version compatibility notes are maintained on the project wiki:

> **https://github.com/JETSCAPE/X-SCAPE/wiki/Doc.Installation**

The summary below covers the most common cases.

### GPU backends and container support

> **GPU hydro backends (`-DUSE_METAL=ON` / `-DUSE_CUDA=ON`) are not yet supported inside containers.**
>
> | Platform | Status |
> |----------|--------|
> | macOS — Metal | **Not supported.** Apple's Metal framework requires bare-metal access to the GPU; it is not available inside Docker Desktop for Mac or any container runtime on macOS. |
> | Linux — CUDA | **In preparation.** NVIDIA GPU pass-through support is being added to the container images. Linux x86\_64 is the primary target; Linux ARM (e.g. Graviton, Jetson, Apple Silicon Linux VM) will follow but may take considerably longer. |
>
> For GPU-accelerated runs, build from source on bare metal (sections [3](#3-obtaining-the-source)–[7](#7-installing-to-a-prefix)).

### Docker (macOS and Linux — CPU builds)

#### Install Docker

- **macOS**: [Docker Desktop for Mac](https://docs.docker.com/docker-for-mac/install/) — set CPUs and memory in Preferences → Advanced.
- **Linux**: [Docker Engine](https://docs.docker.com/install/) — add your user to the `docker` group:
  ```bash
  sudo groupadd docker && sudo usermod -aG docker $USER
  # log out and back in
  ```

#### Pull the image and start a container

```bash
# Create a shared directory on your host
mkdir -p ~/jetscape-docker
cd ~/jetscape-docker
git clone https://github.com/JETSCAPE/X-SCAPE.git
```

**macOS:**
```bash
docker run -it \
  -v ~/jetscape-docker:/home/jetscape-user \
  --name myXSCAPE \
  jetscape/base:stable
```

**Linux:**
```bash
docker run -it \
  -v ~/jetscape-docker:/home/jetscape-user \
  --name myXSCAPE \
  --user $(id -u):$(id -g) \
  jetscape/base:stable
```

For details on image tags and version compatibility see the [jetscape/base DockerHub page](https://hub.docker.com/r/jetscape/base).

#### Build inside the container

From **inside** the container:

```bash
cd /home/jetscape-user/X-SCAPE
mkdir build && cd build
cmake .. -DUSE_MUSIC=ON -DUSE_ISS=ON
make -j$(nproc)
```

Useful container management commands:

| Command | Effect |
|---------|--------|
| `exit` or `docker stop myXSCAPE` | Stop the container |
| `docker start -ai myXSCAPE` | Re-attach to a stopped container |
| `docker container ls -a` | List all containers (running and stopped) |
| `docker container rm myXSCAPE` | Delete the container |

### Singularity / Apptainer (HPC clusters)

Singularity (now [Apptainer](https://apptainer.org/)) can convert the same Docker image for use on HPC systems where Docker is not permitted:

```bash
# Pull and convert (run once; creates a .sif image file)
singularity pull xscape_base.sif docker://jetscape/base:stable

# Start an interactive shell with the X-SCAPE source directory bound in
singularity shell \
  --bind /path/to/X-SCAPE:/opt/xscape \
  xscape_base.sif

# Inside the Singularity shell:
cd /opt/xscape
mkdir -p build && cd build
cmake .. -DUSE_MUSIC=ON -DUSE_ISS=ON
make -j$(nproc)
```

For batch submission, wrap the build or run command in `singularity exec`:

```bash
singularity exec \
  --bind /path/to/X-SCAPE:/opt/xscape \
  xscape_base.sif \
  bash -c "cd /opt/xscape/build && ./runJetscape ../config/jetscape_user.xml"
```

---

## 3. Obtaining the source

```bash
git clone https://github.com/JETSCAPE/X-SCAPE.git
cd X-SCAPE
git submodule update --init --recursive   # GTL, Cornelius, trento, googletest, …
```

---

## 4. Downloading external packages

Each optional physics module has a dedicated download script in `external_packages/`:

| Script | Module | Notes |
|--------|--------|-------|
| `get_music.sh`         | MUSIC CPU hydro | required by `-DUSE_MUSIC=ON` (CPU) |
| `get_music4gpu.sh`     | MUSIC GPU hydro | required by `-DUSE_MUSIC=ON -DUSE_CUDA=ON` or `-DUSE_METAL=ON` |
| `get_iSS.sh`           | iSS particlization | required by `-DUSE_ISS=ON`; also downloads iSS tables |
| `get_3dglauber.sh`     | 3D-MCGlauber + LHAPDF | required by `-DUSE_3DGlauber=ON` |
| `get_lhapdf.sh`        | LHAPDF (standalone) | only needed if 3dMCGlauber LHAPDF_Lib is absent |
| `get_lbtTab.sh`        | LBT rate tables | required at runtime for LBT jet energy loss |
| `get_ipglasma.sh`      | IP-Glasma initial state | required by `-DUSE_IPGLASMA=ON` |
| `get_freestream-milne.sh` | freestream-milne pre-equilibrium | required by `-DUSE_FREESTREAM=ON` |
| `get_clvisc.sh`        | CLVisc GPU hydro | required by `-DUSE_CLVISC=ON` |
| `get_smash.sh`         | SMASH afterburner | required by `-DUSE_SMASH=ON` |
| `get_js_contrib.sh`    | js-contrib extensions | required by `-DUSE_JS_CONTRIB=ON` |

Run each script from the repository root, e.g.:

```bash
cd external_packages
./get_music.sh
./get_iSS.sh
./get_3dglauber.sh
./get_lbtTab.sh
cd ..
```

---

## 5. CMake options reference

### Physics modules

| Option | Default | Description |
|--------|---------|-------------|
| `-DUSE_MUSIC=ON` | OFF | Build MUSIC bulk hydro |
| `-DUSE_CUDA=ON`  | OFF | Use `music4gpu` with NVIDIA CUDA backend (implies `USE_MUSIC`) |
| `-DUSE_METAL=ON` | OFF | Use `music4gpu` with Apple Metal backend — macOS/Apple Silicon only (implies `USE_MUSIC`) |
| `-DUSE_ISS=ON`   | OFF | Build iSS particlization (requires `USE_MUSIC`) |
| `-DUSE_SMASH=ON` | OFF | Build SMASH afterburner (requires `USE_MUSIC`, `USE_ISS`, and `SMASH_DIR`) |
| `-DUSE_3DGlauber=ON` | OFF | Build 3D-MCGlauber initial state (requires GSL) |
| `-DUSE_IPGLASMA=ON`  | OFF | Build IP-Glasma initial state (requires GSL) |
| `-DUSE_FREESTREAM=ON` | OFF | Build freestream-milne pre-equilibrium (requires GSL) |
| `-DUSE_CLVISC=ON`    | OFF | Build CLVisc OpenCL GPU hydro (requires OpenCL) |

`USE_METAL` and `USE_CUDA` are mutually exclusive. `USE_METAL` is silently forced off on non-Apple platforms.

### Output and analysis

| Option | Default | Description |
|--------|---------|-------------|
| `-DUSE_ROOT=ON`  | OFF | Enable ROOT output writer (auto-locates ROOT via `ROOTSYS`) |

HepMC3 support is auto-detected; no flag is needed.

### Build utilities

| Option | Default | Description |
|--------|---------|-------------|
| `-Dunittests=ON`  | ON  | Build Google Test unit-test executables |
| `-DUSE_JS_CONTRIB=ON` | OFF | Enable js-contrib extension modules |
| `-DUSE_JS_FNO_HYDRO=ON` | OFF | Build FnoHydro contrib (requires libtorch ≈ 2 GB and ROOT; needs `USE_JS_CONTRIB`) |
| `-DUSE_JS_PYJETSCAPE=ON` | OFF | Build PyJetscape pybind11 bindings from source (needs `USE_JS_CONTRIB`) |

### Install layout

| Option | Default | Description |
|--------|---------|-------------|
| `-DCMAKE_INSTALL_PREFIX=<path>` | build dir | Root of the installed tree |
| `-DCMAKE_BUILD_TYPE=Release\|Debug` | Release | `Release` enables `-O3`; `Debug` enables `-g -O0` |

---

## 6. Build configurations (quick recipes)

All examples use an out-of-source build directory. Replace `/path/to/build` with your preferred location.

### Minimal (framework only, no physics modules)

```bash
cmake -S . -B /path/to/build
cmake --build /path/to/build -j$(nproc)
```

### Standard CPU physics chain (3DGlauber + MUSIC + iSS)

```bash
cmake -S . -B /path/to/build \
  -DUSE_MUSIC=ON \
  -DUSE_ISS=ON \
  -DUSE_3DGlauber=ON
cmake --build /path/to/build -j$(nproc)
```

### Full CPU chain including SMASH afterburner

```bash
cmake -S . -B /path/to/build \
  -DUSE_MUSIC=ON \
  -DUSE_ISS=ON \
  -DUSE_3DGlauber=ON \
  -DUSE_SMASH=ON \
  -DSMASH_DIR=/path/to/smash/build
cmake --build /path/to/build -j$(nproc)
```

### GPU hydro with Apple Metal (macOS / Apple Silicon)

```bash
cd external_packages && ./get_music4gpu.sh && cd ..
cmake -S . -B /path/to/build \
  -DUSE_MUSIC=ON \
  -DUSE_METAL=ON \
  -DUSE_ISS=ON \
  -DUSE_3DGlauber=ON
cmake --build /path/to/build -j$(nproc)
```

### GPU hydro with CUDA (Linux / NVIDIA)

```bash
cd external_packages && ./get_music4gpu.sh && cd ..
cmake -S . -B /path/to/build \
  -DUSE_MUSIC=ON \
  -DUSE_CUDA=ON \
  -DUSE_ISS=ON \
  -DUSE_3DGlauber=ON
cmake --build /path/to/build -j$(nproc)
```

### With ROOT output

```bash
cmake -S . -B /path/to/build \
  -DUSE_MUSIC=ON -DUSE_ISS=ON \
  -DUSE_ROOT=ON
cmake --build /path/to/build -j$(nproc)
```

---

## 7. Installing to a prefix

After building, install to a clean prefix to get the canonical `bin/`, `lib/`, `share/` layout:

```bash
cmake -S . -B /path/to/build \
  -DCMAKE_INSTALL_PREFIX=/path/to/prefix \
  -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_3DGlauber=ON
cmake --build /path/to/build -j$(nproc)
cmake --install /path/to/build
```

### Installed tree layout

```
<prefix>/
  bin/
    runJetscape          # main simulation driver
    run_in_workdir.sh    # parallel-run launcher
    readerTest
    FinalStateHadrons
    FinalStatePartons
    PythiaBrickTest
    JetScapePerEventTest
    MUSICMainClockTest   # (if USE_MUSIC + USE_ISS)
    MUSICTest            # (if USE_MUSIC + USE_ISS)
    PythiaIsrTest        # (if USE_MUSIC + USE_ISS + USE_3DGlauber)
    iSS.e                # (if USE_ISS)
    3dMCGlb.e            # (if USE_3DGlauber)
  lib/
    libJetScape.{so,dylib}
    libJetScapeReader.{so,dylib}
    libJetScapeThird.{so,dylib}
    libGTL.{so,dylib}
    libCornelius.{so,dylib}
    libhydroFromFile.{so,dylib}
    libtrento.a
    libmusic.{so,dylib}  # (if USE_MUSIC)
    libiSS.{so,dylib}    # (if USE_ISS)
    lib3dMCGlb.{so,dylib}# (if USE_3DGlauber)
  include/
    *.h                  # all framework headers
    XscapeInstallPaths.h # baked data-dir path (generated)
  share/xscape/
    config/              # jetscape_main.xml and all user XML templates
    EOS/                 # MUSIC equation-of-state tables  (if USE_MUSIC)
    iSS_tables/          # iSS particle tables              (if USE_ISS)
    iSS_parameters.dat                                      # (if USE_ISS)
    LBT-tables/          # LBT rate tables                  (if get_lbtTab.sh)
    tables/              # 3dMCGlauber valence-quark caches (if USE_3DGlauber)
    eps09/               # 3dMCGlauber nuclear PDF grids    (if USE_3DGlauber)
    LHAPDF_Lib/          # vendored LHAPDF                  (if USE_3DGlauber)
    nucleusConfigs/      # trento nucleus configuration files
    data_table/          # initial-state data tables
    mcglauber.input      # 3dMCGlauber parameter file       (if USE_3DGlauber)
    music_input          # MUSIC parameter file             (if USE_MUSIC)
```

RPATH is set at install time (`@loader_path/../lib` on macOS, `$ORIGIN/../lib` on Linux), so `runJetscape` finds its shared libraries without `LD_LIBRARY_PATH` / `DYLD_LIBRARY_PATH`.

### Running after install

Once installed, `runJetscape` resolves data assets via the path baked in at configure time, so it works from any CWD with no environment variables:

```bash
/path/to/prefix/bin/runJetscape \
  /path/to/prefix/share/xscape/config/jetscape_user.xml
```

Or with both XMLs explicit:

```bash
runJetscape /path/to/user.xml /path/to/main.xml
```

---

## 8. Running from the build directory

Traditional usage: `cd` into the build directory before running.

```bash
cd /path/to/build
./runJetscape ../config/jetscape_user.xml
```

The build directory contains symlinks / copies of all required data assets placed there by CMake at configure time (`EOS/`, `LBT-tables`, `iSS_tables/`, `tables/`, `nucleusConfigs/`, `mcglauber.input`, `music_input`, …).

---

## 9. Parallel runs with the launcher

`run_in_workdir.sh` lets you run many instances concurrently against **one** shared build or install tree. Each instance gets its own empty working directory for writable output; read-only assets are shared.

### Usage

```
run_in_workdir.sh -w RUN_DIR -u USER_XML [-d DATA_DIR] [-b RUNJETSCAPE] [-m MAIN_XML]

  -w RUN_DIR      Per-run working directory (created if absent).  REQUIRED.
  -u USER_XML     User XML configuration file.                    REQUIRED.
  -d DATA_DIR     Shared asset prefix (build dir or install prefix).
                  Default: $XSCAPE_DATA_DIR, or the directory of runJetscape.
  -b RUNJETSCAPE  Path to runJetscape executable. Default: DATA_DIR/runJetscape.
  -m MAIN_XML     Main XML file. Default: DATA_DIR/config/jetscape_main.xml.
```

### What the launcher does

1. Creates `RUN_DIR`.
2. Sets `XSCAPE_DATA_DIR`, `HYDROPROGRAMPATH`, `LBT_TABLES_PATH` to point at `DATA_DIR`.
3. Writes `RUN_DIR/jetscape_main.xml` with iSS table paths rewritten to absolute.
4. Copies `DATA_DIR/music_input` → `RUN_DIR/music_input` (MUSIC rewrites it in-place; the shared copy stays read-only).
5. Symlinks CWD-relative read-only asset dirs (`tables/`, `eps09/`, `LHAPDF_Lib/`, `nucleusConfigs/`, `data_table/`) into `RUN_DIR`.
6. Runs `runJetscape` from `RUN_DIR`; all outputs land there.

### Example — two concurrent runs from a build directory

```bash
BUILD=/path/to/build
USER_XML=$BUILD/config/jetscape_user_3DGlauber_MUSIC_iSS_test.xml

$BUILD/run_in_workdir.sh -d $BUILD -w /tmp/run1 -u $USER_XML &
$BUILD/run_in_workdir.sh -d $BUILD -w /tmp/run2 -u $USER_XML &
wait
```

### Example — two concurrent runs from an install prefix

```bash
PREFIX=/path/to/prefix
USER_XML=$PREFIX/share/xscape/config/jetscape_user_3DGlauber_MUSIC_iSS_test.xml

$PREFIX/bin/run_in_workdir.sh -w /tmp/run1 -u $USER_XML &
$PREFIX/bin/run_in_workdir.sh -w /tmp/run2 -u $USER_XML &
wait
```

After both finish, outputs are isolated in `/tmp/run1/` and `/tmp/run2/`.

### Notes on parallel safety

- `cache_tables` **must** remain `1` (the default) in `mcglauber.input`. Setting `cache_tables 0` regenerates the valence-quark sample files in `tables/`, which would race across concurrent runs sharing that directory.
- MUSIC writes fixed-name files (`surface*.dat`, `evolution_all_xyeta.dat`, …). Each run being in its own CWD is the only isolation; never point two concurrent instances at the same `RUN_DIR`.

---

## 10. Environment variables

| Variable | Used by | Description |
|----------|---------|-------------|
| `XSCAPE_DATA_DIR` | Framework (`JetScapeDataPath.h`) | Overrides the baked install path and the `./` default for all data asset lookups. Set automatically by `run_in_workdir.sh`. |
| `HYDROPROGRAMPATH` | MUSIC EOS (`eos_base.cpp`) | Directory containing `EOS/`. Falls back to `"."`. Set to `XSCAPE_DATA_DIR` by the launcher. |
| `LBT_TABLES_PATH` | LBT (`LBT.cc`) | Directory containing LBT rate tables (overrides `Eloss/Lbt/LBT_table_path` XML and `XSCAPE_DATA_DIR/LBT-tables`). Set to `XSCAPE_DATA_DIR/LBT-tables` by the launcher. |
| `SMASH_DIR` | CMake only | Points CMake at the SMASH build tree when `-DUSE_SMASH=ON`. |
| `PYTHIA8` / `PYTHIA8DIR` | CMake only | Override Pythia8 search path at configure time. |
| `ROOTSYS` | CMake only | Override ROOT search path when `-DUSE_ROOT=ON`. |

---

## 11. Data-asset search order

Asset lookups follow this priority chain (highest first):

1. **Explicit XML element** — e.g. `<LBT_table_path>`, `<mcglauber_input_file>`, `<iSS_table_path>`, `<MUSIC_input_file>`.
2. **Environment variable** — `LBT_TABLES_PATH` (LBT), `HYDROPROGRAMPATH` (MUSIC EOS).
3. **`XSCAPE_DATA_DIR`** — framework-level override for all other assets.
4. **Baked install path** — `<prefix>/share/xscape`, compiled in at configure time via `XscapeInstallPaths.h`. Only active when the install prefix is a real absolute directory (i.e. not the default build-dir value).
5. **Current working directory** — historical default; always the final fallback.

---

## 12. Interactive configurator

`configure.sh` provides a `dialog`-based TUI for choosing build options without writing CMake command lines. It requires the `dialog` utility:

```bash
# macOS
brew install dialog

# Ubuntu / Debian
sudo apt install dialog
```

Then from the repository root:

```bash
./configure.sh
```

---

## 13. Troubleshooting

**`runJetscape` cannot find `EOS/` tables**  
MUSIC EOS is located via `HYDROPROGRAMPATH`. If running outside the build dir or install prefix, set `HYDROPROGRAMPATH=/path/to/prefix/share/xscape` or use `run_in_workdir.sh`.

**`LBT-tables` not found**  
Download the tables first: `cd external_packages && ./get_lbtTab.sh`. Then set `LBT_TABLES_PATH` or `XSCAPE_DATA_DIR`, or add `<LBT_table_path>` to your user XML.

**`mcglauber.input` not found**  
Download 3D-MCGlauber first: `cd external_packages && ./get_3dglauber.sh`. At runtime the file is resolved as `$XSCAPE_DATA_DIR/mcglauber.input`. Override with `<mcglauber_input_file>` in your user XML.

**`iSS` tables not found**  
Run `./get_iSS.sh`. The tables are installed to `<prefix>/share/xscape/iSS_tables/`; `run_in_workdir.sh` rewrites the paths automatically.

**Shared-library not found at runtime (installed build)**  
The installed `runJetscape` has RPATH set to `<prefix>/lib`. If you moved the prefix after installation, set `DYLD_LIBRARY_PATH` (macOS) or `LD_LIBRARY_PATH` (Linux) to `<new_prefix>/lib`.

**Apple Silicon: `omp.h` not found**  
Install `libomp` via Homebrew (`brew install libomp`). CMake locates it automatically via `brew --prefix libomp`.

**Two runs writing the same output files**  
Use `run_in_workdir.sh` with distinct `-w` directories, or set a unique `<outputFilename>` in each user XML. Never run two instances in the same working directory.

---

## 14. Using X-SCAPE in a downstream CMake project

### Is `find_package(XSCAPE)` available?

Yes. After `cmake --install`, X-SCAPE installs a standard CMake package configuration to:

```
<prefix>/lib/cmake/XSCAPE/XSCAPEConfig.cmake
<prefix>/lib/cmake/XSCAPE/XSCAPETargets.cmake
```

Any downstream CMake project can consume these with `find_package(XSCAPE)` and link against the imported targets.

### Minimal downstream `CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.14)
project(MyAnalysis CXX)

# Tell CMake where the X-SCAPE install lives if it is not in a standard prefix.
# This can also be set on the command line:
#   cmake -DCMAKE_PREFIX_PATH=/path/to/prefix ...
list(APPEND CMAKE_PREFIX_PATH "/path/to/xscape/prefix")

find_package(XSCAPE REQUIRED)

add_executable(myAnalysis src/myAnalysis.cc)
target_link_libraries(myAnalysis PRIVATE XSCAPE::JetScape)
```

CMake will automatically:
- locate the X-SCAPE headers via `XSCAPE::JetScape`'s `INTERFACE_INCLUDE_DIRECTORIES`
- link the framework shared libraries
- bring in the Boost and HDF5 transitive dependencies

### Imported targets

| Target | Purpose |
|--------|--------|
| `XSCAPE::JetScape` | Core simulation framework (use this for custom modules and drivers) |
| `XSCAPE::JetScapeReader` | Standalone output reader (no simulation dependency) |

All other libraries (`JetScapeThird`, `GTL`, `libtrento`, `Cornelius`, `music`, `iSS`, `3dMCGlb`, …) are linked transitively through `XSCAPE::JetScape`; you do not need to list them explicitly.

### Convenience variables set by `find_package`

| Variable | Value |
|----------|-------|
| `XSCAPE_FOUND` | `TRUE` after a successful `find_package` |
| `XSCAPE_INCLUDE_DIR` | `<prefix>/include` — X-SCAPE framework headers |
| `XSCAPE_DATA_DIR` | `<prefix>/share/xscape` — installed data assets |

These paths are **relocatable**: they are resolved relative to the config-file location at use-time, so moving the install tree does not require a re-install.

### Writing a custom JETSCAPE module

A custom module typically inherits from a JETSCAPE base class. A complete minimal example:

```cmake
cmake_minimum_required(VERSION 3.14)
project(MyModule CXX)

find_package(XSCAPE REQUIRED)

add_library(MyJetModule SHARED src/MyJetModule.cc)
target_link_libraries(MyJetModule PRIVATE XSCAPE::JetScape)

# If you also need Pythia8 headers in your own code:
find_package(Pythia8 REQUIRED)
target_include_directories(MyJetModule PRIVATE ${PYTHIA8_INCLUDE_DIR})
target_link_libraries(MyJetModule PRIVATE ${PYTHIA8_LIBRARIES})
```

### Using output-reader utilities only

If your project only reads JETSCAPE output files (no simulation):

```cmake
find_package(XSCAPE REQUIRED)

add_executable(readOutput src/readOutput.cc)
target_link_libraries(readOutput PRIVATE XSCAPE::JetScapeReader)
```

### Locating the install prefix

If X-SCAPE is installed to a non-standard prefix, tell CMake at configure time:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/xscape/prefix
```

Alternatively, set the `XSCAPE_DIR` or `XSCAPE_ROOT` environment variable to the config directory:

```bash
export XSCAPE_DIR=/path/to/xscape/prefix/lib/cmake/XSCAPE
cmake -S . -B build
```

### Data assets in a downstream executable

`XSCAPE_DATA_DIR` is set after `find_package(XSCAPE)` and points at the installed `share/xscape/`. You can pass it as a compile definition:

```cmake
target_compile_definitions(myAnalysis PRIVATE
  XSCAPE_DATA_DIR="${XSCAPE_DATA_DIR}")
```

Alternatively, at runtime `runJetscape` and any driver linked against `XSCAPE::JetScape` will resolve data assets automatically through `JetScapeDataPath.h` (see [section 10](#10-environment-variables)).
