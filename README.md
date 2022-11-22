# JETSCAPE 3.5.1

The [JETSCAPE](http://jetscape.org) simulation framework is an overarching computational envelope for developing complete event generators for heavy-ion collisions.
It allows for modular incorporation of a wide variety of existing and future software that simulates different aspects of a heavy-ion collision.
For a full introduction to JETSCAPE, please see [The JETSCAPE framework](https://arxiv.org/abs/1903.07706).

Please cite [The JETSCAPE framework](https://arxiv.org/abs/1903.07706) if you use this package for scientific work.

## Installation

Please see the [Installation Instructions](https://github.com/JETSCAPE/JETSCAPE/wiki/Doc.Installation).

## Running JETSCAPE

The main executable to generate JETSCAPE events is `runJetscape`, located in the `build/` directory.
To generate JETSCAPE events, you should pass an XML file specifying the settings with which you would like to run:

```
./runJetscape ../config/jetscape_user.xml
```

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
```
./runJetscape /path/to/my/user_config.xml
```

## JETSCAPE Output

JETSCAPE output can be generated in Ascii, gzipped Ascii, or HepMC format,
and contains a full list of particles and the parton shower history.
You must specify which format you would like to activate in your User XML file.

### Analysis of JETSCAPE Output

Analysis of JETSCAPE output is generally beyond the scope of the JETSCAPE framework, and is the responsibility of the user.
The JETSCAPE docker container includes ROOT, python, fastjet, and several other tools that may be useful.

An example reading an ascii output file is provided:

```bash
./build/readerTest
```

which reads in the generated showers does some DFS search and shows the output. You can generate an output graph format which can be easily visualized using graphViz or other tools like Gephi (GUI for free for Mac) or more advanced, graph-tools (Python) and many more. Furthermore, as a "closure" test, the FastJet core package (compiled in our JetScape library) is used to perform a simple jetfinding (on the "final" partons, in graph language, incoming partons in a vertex with no outgoing partons/childs), and since the "shower" is perfectly collinear the jet pT is identical to the hard process parton pT (modulo some random new partons/roots in the final state, see above).  

## JETSCAPE Tunes

Currently, there exists a pp tune [PP19](https://arxiv.org/abs/1910.05481), which can be run by:
```
./runJetscape ../config/jetscape_user_PP19.xml
```

Tuning of Pb-Pb is ongoing.
Several example hydro profiles can be downloaded using `examples/get_hydroSample*`.

## X-SCAPE modules
### 3DGlauber support 

3DGlauber is a 3D initial state model. 3DGlauber generates the initial state for
the MUSIC and can be integrated into the XETSCAPE framework. To download the lastest
version of 3DGlauber, one can run the shell script under the external_packages folder,

```bash
    ./get_3dglauber.sh
```

To compile the 3DGlauber code together with the XETSCAPE framework, please turn
 on the 3DGlauber support option,

```bash
    mkdir build
    cd build
    cmake -DUSE_3DGlauber=ON ..
    make
```

### MUSIC support

MUSIC is a (3+1)D viscous hydrodynamical code developed at McGill university.
(Official website: http://www.physics.mcgill.ca/MUSIC)
MUSIC can be integrated into the JETSCAPE framework. To download the lastest
version of MUSIC, one can run the shell script under the external_packages folder,

```bash
    ./get_music.sh
```

This shell script will clone the latest version of MUSIC to external_packages folder.
It also setup the enviroment variables for MUSIC to run. Specifically, MUSIC
needs the folder path for the EOS tables. Please make sure the enviroment
variable HYDROPROGRAMPATH to be set to the path for MUSIC code package.

When compiling MUSIC with JETSCAPE, please turn on the MUSIC support option
when generating the cmake configuration file,

```bash
    mkdir build
    cd build
    cmake -DUSE_MUSIC=ON ..
    make
```

To run JETSCAPE with MUSIC, one needs to use MPI commands,

```bash
    mpirun -np 1 ./MUSICTest
```
### iSS support 

iSS is a Monte Carlo sampler code after the hydrodynamics and can be integrated 
into the XETSCAPE framework. To download the lastest
version of iSS, one can run the shell script under the external_packages folder,

```bash
    ./get_iSS.sh
```

To compile the iSS code together with the XETSCAPE framework, please turn
on the iSS support option,

```bash
    mkdir build
    cd build
    cmake -DUSE_ISS=ON ..
    make
```

### Running JETSCAPE with CLVisc
In order to run clvisc in jetscape, one has to download it in external\_packages/, using 

```
sh get_clvisc.sh
```

Then compile and run the framework with a XML configuration file which turns clvisc on.
```
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
Prof. H. Elfner (nee Petersen).  In JetScape SMASH can
serve as an afterburner, useful to compute soft observables.

### Installing SMASH

SMASH is published on github at https://github.com/smash-transport/smash.
See SMASH Readme for libraries required by SMASH and how to install them.

```
  export EIGEN3_ROOT=<eigen install directory>/include/eigen3/
  export GSL_ROOT_DIR=$(gsl-config --prefix)
  export BOOST_ROOT=<boost install directory>
  export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8235
  export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8235

  export JETSCAPE_DIR=${HOME}/JETSCAPE-COMP
  export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code

  cd ${JETSCAPE_DIR}/external_packages
  ./get_smash.sh
```

### Compiling JetScape with SMASH

The usage of SMASH in JetScape as an afterburner requires hydro,
sampler and SMASH itself. Therefore, to use it in JetScape,

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

Currently the iSS sampler doesn't performs resonance decays after sampling.
For reasonable physics with SMASH these decays should be switched off.

## More information

More material on the physics behind JETSCAPE and how to use it can be found in the material of the JETSCAPE Summer Schools. The schools are a yearly event explaining the details of the approach. You can either sign-up for the next one or go through the material of the last school yourself. The material is found in the SummerSchool repositories under [the JETSCAPE organization](https://github.com/JETSCAPE).

## Troubleshooting

If you encounter a problem, please report the issue [here](https://github.com/JETSCAPE/JETSCAPE/issues).
Please be sure to include enough information so that we can reproduce your issue: your platform, JETSCAPE version, configuration file, and anything else that may be relevant.

# Contributing

Please see the [CONTRIBUTING](CONTRIBUTING.md) for instructions how to do contribute to the framework and development hints.
