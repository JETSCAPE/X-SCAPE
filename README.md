# X-SCAPE 3.0

The X-ion collisions with a Statistically and Computationally Advanced Program Envelope (X-SCAPE) is the enhanced (and 2nd) project of the JETSCAPE
collaboration which extends the framework to include small systems created in p-A and p-p collisions, lower energy heavy-ion collisions and electron-Ion collisions.
The new framework allows for novel functionality such as the ability of the main simulation clock to go backwards and forwards, to deal systematically with initial
state and final state evolution. It allows for multiple bulk event generators to run concurrently while exchanging information via a new Bulk Dynamics Manager.
The X-SCAPE framework can be run using the new functionality or in JETSCAPE mode allowing for full backwards compatibility. New modules can also run in a hybrid fashion,
choosing to use or not use the new clock functionality.

The [JETSCAPE](http://jetscape.org) simulation framework is an overarching computational envelope for developing complete event generators for heavy-ion collisions.
It allows for modular incorporation of a wide variety of existing and future software that simulates different aspects of a heavy-ion collision.
For a full introduction to JETSCAPE, please see [The JETSCAPE framework](https://arxiv.org/abs/1903.07706).

Please cite [The JETSCAPE framework](https://arxiv.org/abs/1903.07706) if you use this package for scientific work.

## What is new in 3.0

The highlights below refer to X-SCAPE together with its Python interface, PyJetscape
(`import jetscape`, from [js-contrib](https://github.com/jhputschke/js-contrib)): some of them
are used through the Python interface, as the box below explains. The full list of changes
since X-SCAPE 2.1, including the physics bug fixes that change the output of existing setups,
is in [RELEASE_NOTES.md](RELEASE_NOTES.md).

> **ℹ️ Important: some of the new functionality needs PyJetscape, the Python interface from
> [js-contrib](https://github.com/jhputschke/js-contrib).** The new physics modules and settings
> run from the XML with `runJetscape`, and the two-stage hydro works there in part: see the
> example [`config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml`](config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml)
> ([below](#two-stage-hydro-and-music4gpu-from-xml)). For now the Python interface is required for:
>
> - the hadrons of both legs of an event (jet leg and background, to subtract one from the other),
> - several jets per background, with iSS on every jet leg,
> - oversampling per leg, and correlated oversampling of the two legs with common seeds,
> - explicit seeds per event and leg, to redo one event exactly,
> - the HDF5 production writers (hydro pair, particlize and hadron files) and the hadronization
>   of stored surfaces and partons later on CPUs (stage 2),
> - modules written in Python and the live PyTorch FNO hydro.
>
> Running the first five from the XML in C++ is **work in progress and will be provided in future
> releases**. PyJetscape comes compiled in the production container and the Mac wheel
> ([Install](https://jhputschke.github.io/js-contrib/install/)).

- **GPU hydro with MUSIC4GPU** (CUDA on NVIDIA, Metal on Apple Silicon), a drop-in replacement
  for CPU MUSIC, with its conservation fixes on by default and an opt-in active region:
  [docs/README_MUSIC4GPU.md](docs/README_MUSIC4GPU.md),
  [GPU acceleration](https://jhputschke.github.io/js-contrib/about/gpu/).
- **Two-stage hydro with restart**: the background runs once, and each jet leg (with the
  liquefier's droplets) restarts from the background's checkpoint:
  [docs/README_MUSIC4GPU.md](docs/README_MUSIC4GPU.md#two-stage-hydro-with-runjetscape),
  [Au+Au jets in two stages](https://jhputschke.github.io/js-contrib/run/auau-jet/).
  A worked example that runs from the XML alone is
  `config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml` (see
  [Two-stage hydro and MUSIC4GPU from XML](#two-stage-hydro-and-music4gpu-from-xml) below).
- **O+O with initial-state radiation**: 3dMCGlauber, PythiaIsrGun and iMATTER, MUSIC, Matter and
  LBT, and iColored hadronization through the jets' colour flow, all from the XML
  (`<Hard><ISR>`): [O+O with ISR](https://jhputschke.github.io/js-contrib/run/oo-isr/),
  [docs/README_ColourFlow.md](docs/README_ColourFlow.md).
- **Extra hard scatterings** (`<extraHard>`, with a minimum-bias mode):
  [Extra hard scatterings](https://jhputschke.github.io/js-contrib/run/extra-hard/).
- **Parallel showers**: the showers of one event on OpenMP threads (`<Eloss><nThreads>`), with
  an output that does not depend on the thread count:
  [release notes](RELEASE_NOTES.md#energy-loss-and-hadronization).
- **Faster iSS**: yields and delta f cached, correlated sampling, and sampling on OpenMP
  threads, hadrons bit-identical: [release notes](RELEASE_NOTES.md#particlization-iss).
- **HDF5 files and two-stage hadronization**: hydro surfaces and partons are stored on the GPU
  machine and hadronized later on CPUs, with reproducible seeds:
  [What a production writes](https://jhputschke.github.io/js-contrib/analysis/files/).
- **Install prefix and launcher**: `make install`, `find_package(XSCAPE)`, and many jobs from
  one build: [INSTALL.md](INSTALL.md) (sections 7 and 9),
  [docs/README_Launcher.md](docs/README_Launcher.md).
- **Python interface**: `import jetscape` drives X-SCAPE from Python and runs Python modules
  in the event loop (js-contrib's PyJetscape):
  [Python](https://jhputschke.github.io/js-contrib/python/).
- **Physics bug fixes** that change the output of existing setups (liquefier, CausalLiquefier,
  energy-loss start time, colour tags, MUSIC wrapper, SurfaceFinder, MUSIC4GPU fluxes):
  [release notes](RELEASE_NOTES.md#physics-bug-fixes). Do not combine files made before and
  after them.

## Documentation and installation

The main guide to installing and using X-SCAPE is the js-contrib website,
**<https://jhputschke.github.io/js-contrib/>**:

- [Install](https://jhputschke.github.io/js-contrib/install/): the production container, the
  xscape wheel for Apple Silicon Macs, an analysis-only environment, and the source build.
- [Run](https://jhputschke.github.io/js-contrib/run/): a first event, Au+Au jets in two stages,
  O+O with initial-state radiation, Z+jet, extra hard scatterings, clusters and GPUs.
- [Analyse](https://jhputschke.github.io/js-contrib/analysis/): what a production writes, and how
  to read hadrons, jets, hydro and bulk observables from it.
- [Python](https://jhputschke.github.io/js-contrib/python/): the Python interface and Python
  modules inside the event loop.

In this repository:

| Document | Content |
|---|---|
| [INSTALL.md](INSTALL.md) | Source build: prerequisites, external packages, CMake options, install prefix, launcher, environment variables |
| [RELEASE_NOTES.md](RELEASE_NOTES.md) | Changes since X-SCAPE 2.1 |
| [docs/README_MUSIC4GPU.md](docs/README_MUSIC4GPU.md) | The GPU hydro: build, OpenMP defaults, settings, two-stage hydro with `runJetscape` |
| [docs/README_ColourFlow.md](docs/README_ColourFlow.md) | Colour flow through Matter and LBT, and iColored hadronization |
| [docs/README_Launcher.md](docs/README_Launcher.md) | Many jobs in parallel from one build or install |
| [docs/README_BulkFast.md](docs/README_BulkFast.md) | `FastRootBulkWriter`: hydro evolution to ROOT |
| [CONTRIBUTING.md](CONTRIBUTING.md) | How to contribute |

### Which container is for what

| | For | Hardware |
|---|---|---|
| JETSCAPE images ([`jetscape/base`](https://hub.docker.com/r/jetscape/base)) | Development: build X-SCAPE from source inside them ([INSTALL.md](INSTALL.md), section 2) | CPU only |
| js-contrib production container ([`jhputschke/xscape-prod`](https://jhputschke.github.io/js-contrib/install/container/)) | Productions, your own scripts and notebooks, analysis; X-SCAPE comes compiled | NVIDIA GPU (Docker, Apptainer); CPU-only nodes for hadronization and analysis |
| xscape wheel ([Mac install](https://jhputschke.github.io/js-contrib/install/mac/)) | The same software as the production container, installed natively | Apple Silicon Mac (Metal GPU) |

GPU development containers are in progress. Until then, GPU development means a source build on
the machine itself ([INSTALL.md](INSTALL.md), sections 3 to 7, or the
[source build](https://jhputschke.github.io/js-contrib/install/source/) with conda).

## Building from source

The short version; [INSTALL.md](INSTALL.md) has every option and recipe.

```bash
cd external_packages
./get_music4gpu.sh      # or ./get_music.sh for CPU MUSIC
./get_iSS.sh
./get_3dglauber.sh
./get_js_contrib.sh     # optional: js-contrib (Python interface, FNO hydro)
cd ..
mkdir build && cd build
cmake .. -DUSE_MUSIC=ON -DUSE_CUDA=ON -DUSE_ISS=ON -DUSE_3DGlauber=ON   # -DUSE_METAL=ON on a Mac
make -j$(nproc)
```

`./configure.sh` offers the build options in a text menu (needs `dialog`).

## Running X-SCAPE/JETSCAPE

The main executable is `runJetscape` in the build directory. X-SCAPE is fully backwards
compatible: it contains the same functionalities and features as the latest JETSCAPE release.
Pass it an XML file with the settings you would like to run with:

```bash
./runJetscape ../config/jetscape_user.xml
```

The small-system chain with initial-state radiation runs from the XML as well (`<Hard><ISR>`),
for example `config/jetscape_user_OO_ISR_MatterLBT_iColored.xml`. Dedicated executables remain
for some X-SCAPE examples:

```bash
./PythiaIsrTest           # 3dMCGlauber + iMATTER initial-state shower
./SMASHInitialCondition   # SMASH initial condition
```

### Two-stage hydro and MUSIC4GPU from XML

[`config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml`](config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml)
runs one 0–10 % Au+Au event with a jet from `runJetscape` alone: 3dMCGlauber and PythiaGun, the
background hydro (MUSIC_1), Matter + LBT with the liquefier, the jet-leg hydro with the wake
(MUSIC_2), iSS on the jet leg's surface and Colorless hadronization, written as ASCII. It shows
the new per-instance MUSIC settings, MUSIC4GPU's late time step and iSS choosing the hydro it
samples (`<hydro_id>`). With a CUDA build on a GB10 the single event takes 22 s, of which about
7 s is LBT reading its tables at startup; `export LBT_TABLE_CACHE=/path/to/lbt_tables.cache`
removes that penalty after the first job (15 s; about 13 s per further event).

```bash
cd build
./runJetscape ../config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml
```

**Current limits of the XML compared with Python.** One job gives the hadrons of one leg only
(no paired jet leg minus background), cannot reuse a background for several jets, keeps no
boundary between oversamples, takes no explicit seeds per event and leg, and cannot hadronize
stored surfaces and partons later. The js-contrib productions do all of this from Python with
the same XML files ([Au+Au jets in two stages](https://jhputschke.github.io/js-contrib/run/auau-jet/)).
**XML support for these is in the works for future releases**: iSS following the hydro it
samples with background reuse, a second soft-particlization pass for the other leg, seed tags,
and oversampling per leg. The details, with a table of what runs from XML and what needs Python,
are in [docs/README_MUSIC4GPU.md](docs/README_MUSIC4GPU.md#example-one-two-stage-auau-event-from-xml).

### More examples

More test examples showcasing the new X-SCAPE framework functionalities, such as the main clock
and the Bulk Dynamics Manager, can be found in `./examples/custom_examples/` (for example
`PythiaBDMTest.cc` and `PythiaBrickTest.cc`). To run several jobs at once from one build, each
in a working directory of its own, use `examples/run_in_workdir.sh`
([docs/README_Launcher.md](docs/README_Launcher.md)).

### The XML Configuration

All of the JETSCAPE settings are specified by two XML files:
- Main XML file: *you don't modify this*
  - Contains default values for every possible module and parameter
- User XML file: *you provide this*
  - Contains a list of which modules to run, and which default parameter values to override

An example User XML file is provided at `config/jetscape_user.xml`.
You should adapt this as you like:
 - Set number of events to run
 - Set output format (`Writer` type) and filename
 - Set which modules to include (in order of execution)
 - Set any default parameter values (from Main XML file) to override

The Main XML file is located at `config/jetscape_main.xml`, and contains a list of
the default parameter settings which will be used for all activated modules (as specified by the User XML file),
if they are not overridden in the User XML file.

You can pass the path to your user XML file as a command-line argument to the `runJetscape` executable:
```bash
./runJetscape /path/to/my/user_config.xml
```

## JETSCAPE Output

JETSCAPE output can be generated in Ascii, gzipped Ascii, or HepMC format,
(for HepMC format in ROOT see `examples/custom_examples/PythiaBrickTestRoot.cc` and use cmake with `-DUSE_ROOT=ON`) and contains a full list of particles and the parton shower history.
You must specify which format you would like to activate in your User XML file.
With ROOT, `RootBulkWriter` and `FastRootBulkWriter` write the hydro evolution
([docs/README_BulkFast.md](docs/README_BulkFast.md)). The js-contrib productions write HDF5
files ([What a production writes](https://jhputschke.github.io/js-contrib/analysis/files/)).

### Analysis of JETSCAPE Output

Analysis of JETSCAPE output is generally beyond the scope of the JETSCAPE framework, and is the responsibility of the user.
The JETSCAPE docker container includes ROOT, python, fastjet, and several other tools that may be useful;
the js-contrib website has analysis guides and notebooks for the production files.

An example reading an ascii output file is provided:

```bash
./build/readerTest
```

which reads in the generated showers does some DFS search and shows the output. You can generate an output graph format which can be easily visualized using graphViz or other tools like Gephi (GUI for free for Mac) or more advanced, graph-tools (Python) and many more. Furthermore, as a "closure" test, the FastJet core package (compiled in our JetScape library) is used to perform a simple jetfinding (on the "final" partons, in graph language, incoming partons in a vertex with no outgoing partons/childs), and since the "shower" is perfectly collinear the jet pT is identical to the hard process parton pT (modulo some random new partons/roots in the final state, see above).

## JETSCAPE Tunes

There are several JETSCAPE tunes from publications available.
The XML files to reproduce the results are available in the [Default-tunes](https://github.com/JETSCAPE/Default-tunes) repository.
Please see [JETSCAPE Tunes](https://github.com/JETSCAPE/Default-tunes/blob/main/README.md) for more information.

Several example hydro profiles can be downloaded using `examples/get_hydroSample*`.

## Running JETSCAPE/X-SCAPE with LHAPDF

[LHAPDF](https://www.lhapdf.org/) is included in the JETSCAPE/X-SCAPE Docker image ([INSTALL.md](INSTALL.md), section 2). While inside the container, use [external_packages/get_lhapdf.sh](external_packages/get_lhapdf.sh) to download the desired set. For example, to download the JAM20-SIDIS_PDF_proton_nlo set, run:

```bash
cd ${JETSCAPE_DIR}/external_packages
./get_lhapdf.sh JAM20-SIDIS_PDF_proton_nlo
```

The above script sets an environment variable `LHAPATH='/home/jetscape-user/.local/share/LHAPDF'`. Adjust as needed.

Add the set to PythiaGun in your user XML file:

```xml
<LinesToRead>
  PDF:useHard = on
  PDF:pHardSet LHAPDF6:JAM20-SIDIS_PDF_proton_nlo
</LinesToRead>
```
Or run the example [config/jetscape_user_nPDF_test.xml](config/jetscape_user_nPDF_test.xml) user XML file.

```bash
cd ${JETSCAPE_DIR}/build
./runJetscape ../config/jetscape_user_nPDF_test.xml
```

## X-SCAPE modules

### 3DGlauber support

3DGlauber is a 3D initial state model. 3DGlauber generates the initial state for MUSIC and can be integrated into the JETSCAPE framework. To download the latest version of 3DGlauber, one can run the shell script under the external_packages folder,

```bash
    ./get_3dglauber.sh
```

To compile the 3DGlauber code together with the JETSCAPE framework, please turn on the 3DGlauber support option,

```bash
    mkdir build
    cd build
    cmake -DUSE_3DGlauber=ON ..
    make -j4
```

For more details see [3DGlauber](https://github.com/JETSCAPE/X-SCAPE/wiki/3DGlauber,-MUSIC,-iSS-and-Initial-State-Radiation).

### Initial State Shower using iMatter

To use the ISR shower of iMatter, please make sure that the environment variable `$PYTHIA8` is set and points to the directory where pythia8 is installed which can be found using `pythia8-config --prefix` (with conda, iMATTER finds the PDFs through `PYTHIA8DATA`). In the JETSCAPE docker container all environment variables, including the `$PYTHIA8` for iMatter, are set properly.

After using 3DGlauber support to compile JETSCAPE, one can use `./PythiaIsrTest` (in the build directory) to run iMatter and 3DGlauber which uses the xml user file `config/jetscape_user_iMATTERMCGlauber.xml`. For running 3DGlauber with Hydro (Music) please see [3DGlauber](https://github.com/JETSCAPE/X-SCAPE/wiki/3DGlauber,-MUSIC,-iSS-and-Initial-State-Radiation).

The ISR chain with jet energy loss in the medium (Matter + LBT) and the jets hadronized through their colour flow (`iColoredHadronization`) runs from the XML alone (`<Hard><ISR>`): `config/jetscape_user_OO_ISR_MatterLBT_iColored.xml`. How the colour tags of Matter and LBT and the closure of the colour lines that end in the medium (`<JetHadronization><colour_closure>`) work is in [docs/README_ColourFlow.md](docs/README_ColourFlow.md).

### Electron-Proton Gun

X-SCAPE includes an electron-proton gun capable of recreating both Deep-Inelastic Scattering (DIS) and photoproduction events. It should be run in photoproduction mode for when the virtuality of the exchanged photon (normally called Q^2) is under 1 GeV^2, and in DIS mode when it is larger. The photoproduction mode works well for all HERA energies and DIS can be run down to HERMES energies. It should not be run under this.

### SMASH Initial Condition

The SMASH Initial Condition module provides a hadronic initial state for heavy-ion collisions at low and intermediate energies. It can be used to provide the initial state for hydrodynamic evolution with MUSIC. The module is based on the SMASH hadronic transport approach, which is also used as a hadronic afterburner in JETSCAPE. For more details on how to install the SMASH Initial Condition module, see the [SMASH section](#smash-hadronic-afterburner) below.
For example XML files have a look at the [JETSCAPE tunes](#jetscape-tunes) section above.
For running with the SMASH initial condition, please use the example `./SMASHInitialCondition` in the build directory using example XML files provided in the aforementioned directory. The default XML file is `config/jetscape_user_SMASHInitialCondition.xml`.

## JETSCAPE modules

Since X-SCAPE is fully backwards compatible, all JETSCAPE modules can be used in X-SCAPE utilizing the JETSCAPE like per-event execution. Additional functions have to be implemented to extend towards per-time-step execution (normally achieved by refactoring the per-event code) using the new clock feature, allowing full concurrent running of all physics modules. More details are in [CONTRIBUTING](CONTRIBUTING.md).

### MUSIC support

MUSIC is a (3+1)D viscous hydrodynamical code developed at McGill university.
(Official website: http://www.physics.mcgill.ca/MUSIC)
MUSIC can be integrated into the JETSCAPE framework. To download the latest version of MUSIC, one can run the shell script under the external_packages folder,

```bash
    ./get_music.sh
```

This shell script will clone the latest version of MUSIC to external_packages folder.
It also setup the environment variables for MUSIC to run. Specifically, MUSIC
needs the folder path for the EOS tables. Please make sure the environment
variable HYDROPROGRAMPATH to be set to the path for MUSIC code package.

When compiling MUSIC with JETSCAPE, please turn on the MUSIC support option
when generating the cmake configuration file,

```bash
    mkdir build
    cd build
    cmake -DUSE_MUSIC=ON ..
    make
```

The GPU port MUSIC4GPU replaces CPU MUSIC (one or the other, never both): fetch it with
`./get_music4gpu.sh` and add `-DUSE_CUDA=ON` (NVIDIA) or `-DUSE_METAL=ON` (Apple Silicon).
Its settings, the OpenMP defaults of GPU builds and the two-stage hydro with `runJetscape` are
in [docs/README_MUSIC4GPU.md](docs/README_MUSIC4GPU.md).

To run JETSCAPE with MUSIC, one needs to use MPI commands,

```bash
    mpirun -np 1 ./MUSICTest
```
### iSS support

iSS is a Monte Carlo sampler code after the hydrodynamics and can be integrated
into the JETSCAPE framework. X-SCAPE fetches it from the fork
[jhputschke/iSS](https://github.com/jhputschke/iSS) (see [RELEASE_NOTES.md](RELEASE_NOTES.md#particlization-iss)).
To download it, one can run the shell script under the external_packages folder,

```bash
    ./get_iSS.sh
```

To compile the iSS code together with the JETSCAPE framework, please turn
on the iSS support option,

```bash
    mkdir build
    cd build
    cmake -DUSE_ISS=ON ..
    make
```

### Running JETSCAPE with CLVisc
In order to run clvisc in JETSCAPE, one has to download it in external\_packages/, using

```bash
sh get_clvisc.sh
```

Then compile and run the framework with a XML configuration file which turns clvisc on.
```bash
cd build/
cmake .. -DUSE_CLVISC=on
make
./runJetscape ../config/jetscape_clvisc.xml
```
If the cmake fails because OpenCL is not installed, please check it.
OpenCL is delivered in MacOS by default.
If you use linux machine with Nvidia GPUs, you will need to install CUDA,
which will provide OpenCL support.
If you use linux machine with AMD GPUs, Intel GPUs or any CPUs,
you will need to install AMD APP SDK.

## SMASH hadronic afterburner

SMASH [https://smash-transport.github.io] is a hadronic transport approach
developed at Frankfurt University and GSI by the group of
Prof. H. Elfner (nee Petersen). In JetScape SMASH can
serve as an afterburner (and in X-SCAPE as initial condition), useful to 
compute soft observables.

### Installing SMASH

SMASH is published on github at https://github.com/smash-transport/smash.
See SMASH Readme for libraries required by SMASH and how to install them.
When using the docker container, SMASH should be installed inside, since it 
relies on external libraries available in the container. 

```bash
  export EIGEN3_ROOT=<eigen install directory>/include/eigen3/
  export GSL_ROOT_DIR=$(gsl-config --prefix)
  export BOOST_ROOT=<boost install directory>
  export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8310
  export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8310

  export JETSCAPE_DIR=${HOME}/JETSCAPE-COMP
  export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code

  cd ${JETSCAPE_DIR}/external_packages
  ./get_smash.sh
```

To set the `SMASH_DIR` variable directly when using the docker container, one 
can run `source get_smash.sh`.

### Compiling JetScape with SMASH

The usage of SMASH in JETSCAPE as an afterburner requires hydro,
sampler and SMASH itself. Therefore, to use it in JETSCAPE,

```bash
    mkdir ${JETSCAPE_DIR}/build
    cd ${JETSCAPE_DIR}/build
    cmake -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_FREESTREAM=ON -DUSE_SMASH=ON ..
```

To run JetScape test with SMASH:

```bash
    cd build
    ./SMASHTest
```

With SMASH as afterburner, the iSS resonance decays must stay off
(`<Perform_resonance_decays>0`, the default), so that SMASH decays the resonances itself.

## js-contrib community extensions

[js-contrib](https://github.com/jhputschke/js-contrib) provides community-maintained physics modules
(analogous to [fastjet-contrib](https://fastjet.hepforge.org/contrib/)) that are built
together with X-SCAPE, and the production drivers and analysis tools behind the js-contrib
website. Current contribs:

| Contrib | Description | Extra deps |
|---------|-------------|------------|
| PyJetscape | pybind11 Python bindings, the Pythonic `jetscape` package, production drivers and analysis tools | pybind11 (PyTorch for the FNO hydro) |
| FnoHydro | Neural-network (FNO) hydrodynamics via LibTorch | ROOT, libtorch (~2 GB) |
| FastHydro | MC-Glauber and a fast pure-Python 3+1D hydro with X-SCAPE's energy-loss chain: paired background and jet evolutions | Python (PyJetscape) |
| Visualization | 3D visualization of the hydro medium evolution (PyVista) | Python (PyJetscape) |

Build them from source with:

```bash
cd external_packages
./get_js_contrib.sh        # clones js-contrib into external_packages/js-contrib

cd ../build
cmake .. -DUSE_JS_CONTRIB=ON -DUSE_JS_PYJETSCAPE=ON   # add -DUSE_JS_FNO_HYDRO=ON for the LibTorch FNO
make -j$(nproc)
```

The production container and the Mac wheel come with PyJetscape compiled
([Install](https://jhputschke.github.io/js-contrib/install/)).

## More information

More material on the physics behind JETSCAPE and how to use it can be found in the material of the JETSCAPE Summer Schools. The schools are a yearly event explaining the details of the approach. You can either sign-up for the next one or go through the material of the last school yourself. The material is found in the SummerSchool repositories under [the JETSCAPE organization](https://github.com/JETSCAPE).

## Troubleshooting

If you encounter a problem, please report the issue [here](https://github.com/JETSCAPE/X-SCAPE/issues).

# Contributing

Please see the [CONTRIBUTING](CONTRIBUTING.md) for instructions how to do contribute to the framework and development hints.
