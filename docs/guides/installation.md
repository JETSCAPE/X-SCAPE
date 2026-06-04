# Installation

This page is a quick-start summary. The complete, maintained reference is
[`INSTALL.md`](https://github.com/JETSCAPE/X-SCAPE/blob/main/INSTALL.md) and the
[installation wiki](https://github.com/JETSCAPE/X-SCAPE/wiki/Doc.Installation).

## Easiest path: Docker

The official JETSCAPE container ships a complete pre-built prerequisite
environment (compilers, Boost, PYTHIA 8, HDF5, ROOT, FastJet, LHAPDF, …). Follow
the [Docker instructions](https://github.com/JETSCAPE/X-SCAPE/wiki/Doc.Installation.Docker.Linux)
and skip the prerequisite installation entirely.

## Prerequisites (native build)

| Package | Min version | Notes |
|---|---|---|
| CMake | 3.5 (3.14+ recommended) | |
| C++17 compiler | GCC 8 / Clang 9 / AppleClang 12 | |
| Boost | 1.67 | `filesystem`, `program_options`, `system` |
| PYTHIA 8 | 8.3 | set `$PYTHIA8` / `$PYTHIA8DIR` |
| HDF5 | 1.8 | C and C++ bindings |
| zlib | any | |

Optional (auto-detected): **GSL** (needed by most physics options), **HepMC3**,
**ROOT** (`-DUSE_ROOT=ON`), **MPI** (SMASH), **OpenMP**, **OpenCL** (CLVisc).

On macOS:

```bash
brew install boost pythia8 hdf5 gsl libomp cmake
```

## 1. Get the source and external packages

```bash
git clone https://github.com/JETSCAPE/X-SCAPE.git
cd X-SCAPE/external_packages
# fetch only the engines you intend to build:
./get_music.sh        # MUSIC hydro
./get_iSS.sh          # iSS particlization
./get_3dglauber.sh    # 3DGlauber initial state
./get_smash.sh        # SMASH afterburner
./get_freestream.sh   # freestream-milne pre-equilibrium
# GPU MUSIC drop-in (instead of get_music.sh):
./get_music4gpu.sh
```

Each external package is a separate engine wrapped by an X-SCAPE adapter in
`src/`; you only fetch and enable the ones you need.

## 2. CMake options

| Option | Builds |
|---|---|
| `-DUSE_MUSIC=ON` | MUSIC bulk hydro |
| `-DUSE_CUDA=ON` | `music4gpu` NVIDIA/CUDA backend (implies `USE_MUSIC`) |
| `-DUSE_METAL=ON` | `music4gpu` Apple Metal backend (implies `USE_MUSIC`) |
| `-DUSE_ISS=ON` | iSS particlization (needs `USE_MUSIC`) |
| `-DUSE_SMASH=ON` | SMASH afterburner (needs `USE_MUSIC`, `USE_ISS`, `SMASH_DIR`) |
| `-DUSE_3DGlauber=ON` | 3D-MCGlauber initial state (needs GSL) |
| `-DUSE_IPGLASMA=ON` | IP-Glasma initial state (needs GSL) |
| `-DUSE_FREESTREAM=ON` | freestream-milne pre-equilibrium (needs GSL) |
| `-DUSE_CLVISC=ON` | CLVisc OpenCL GPU hydro (needs OpenCL) |
| `-DUSE_ROOT=ON` | ROOT output / bulk writers |

`USE_CUDA` and `USE_METAL` are mutually exclusive.

## 3. Build recipes

=== "Framework only"

    ```bash
    cmake -S . -B build
    cmake --build build -j$(nproc)
    ```

=== "Standard CPU chain"

    ```bash
    cmake -S . -B build \
      -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_3DGlauber=ON
    cmake --build build -j$(nproc)
    ```

=== "Full chain + SMASH"

    ```bash
    cmake -S . -B build \
      -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_3DGlauber=ON \
      -DUSE_SMASH=ON -DSMASH_DIR=/path/to/smash/build
    cmake --build build -j$(nproc)
    ```

=== "GPU hydro (CUDA)"

    ```bash
    cmake -S . -B build_gpu -DUSE_MUSIC=ON -DUSE_CUDA=ON
    cmake --build build_gpu -j$(nproc)
    ```

!!! tip "OpenMP on many-core GPU hosts"
    GPU (`music4gpu`) builds default `OMP_WAIT_POLICY=passive` so idle OpenMP
    threads sleep instead of busy-spinning between GPU waits — important on
    many-core hosts. Override with `export OMP_WAIT_POLICY=active` or
    `export MUSIC_OMP_DEFAULTS=0`. See the project README and
    `external_packages/music4gpu/PORT_GPU.md` §9.10.

## 4. Verify the build

```bash
cd build
./runJetscape ../config/jetscape_user.xml   # a 2-event p+p test
ctest                                        # run the unit tests (if -Dunittests=ON)
```

## Interactive configurator

`./configure.sh` is an interactive helper that asks which modules you want and
emits the matching `cmake` command — handy for first-time setup. See
[`INSTALL.md` §12](https://github.com/JETSCAPE/X-SCAPE/blob/main/INSTALL.md).

## Community extensions (js-contrib)

[js-contrib](https://github.com/jhputschke/js-contrib) provides
community-maintained modules (e.g. FNO neural-network hydro, Python bindings)
compiled directly into X-SCAPE:

```bash
cd external_packages && ./get_js_contrib.sh
cd ../build && cmake .. -DUSE_JS_CONTRIB=ON -DUSE_JS_FNO_HYDRO=ON && make -j$(nproc)
```
