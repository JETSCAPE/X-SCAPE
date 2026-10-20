# JEWEL 2.6.0 — C++23 port

A C++23 translation of the JEWEL 2.6.0 jet-quenching event generator
(K. C. Zapp, jewel.hepforge.org). The JEWEL core (`jewel-2.6.0.f`), both medium
models (`medium-simple.f`, `medium-vac.f`), the exponential integral (`meix.f`)
and the thermal-momentum subtraction (`subtraction.C`) are ported to C++. The
modified PYTHIA 6.4.25 distributed with JEWEL is **unchanged Fortran** and is
linked as a library; JEWEL and PYTHIA exchange data through the PYTHIA common
blocks exactly as before.

**This is a modified version of JEWEL** (MCnet guidelines §1.3): the program
logic is unchanged, only the implementation language of the JEWEL part differs.
Please cite JEWEL as described in the original README (manual, JHEP 1303 (2013)
080, EPJC 60 (2009) 617, and the references for V+jet, recoils, ISR quenching and
MPI) and refer to results as "JEWEL+PYTHIA". The copyright message is printed
unchanged, with one added line identifying the port. License: GPL v2 (`gpl-2.0.txt`).

## Validation

The port reproduces the original **event by event, bit for bit**: for the same
seed (`NJOB`) the HepMC output is byte-identical and the final random-number
state (`out.ran`) is identical, with lookup tables identical value by value.
This was checked for 5.02 TeV pp (200 000 events), Pb+Pb 0–10 % without recoils
(100 000) and with recoils + constituent subtraction (100 000), against the
original compiled with gfortran 16 `-O2 -ffp-contract=off` on macOS arm64.
`ctest` re-checks this on short runs.

Bit-identity depends on floating-point details (no FMA contraction, the same
libm, gfortran's handling of `sin`/`cos` pairs — see `PORTING.md`). Built with
other compilers or flags, results agree statistically but not event by event;
the same is true for the original Fortran code.

## Build

Requirements: a C++23 compiler (tested: Apple clang 21), gfortran (for PYTHIA 6),
CMake ≥ 3.25, LHAPDF 6 with the CT14nlo set (default `PDFSET 13100`).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/lhapdf
cmake --build build -j
cd build && ctest
```

This produces the same two executables as the Fortran distribution:
`jewel-2.6.0-vac` (vacuum / pp) and `jewel-2.6.0-simple` (simple medium model).
Options: `-DJEWEL_NO_FP_CONTRACT=OFF` allows FMA (faster on some CPUs, no longer
bit-identical); `-DJEWEL_BOUNDS_CHECK=ON` enables array bounds assertions.

## Running

Usage and parameter files are unchanged from the Fortran version:

```bash
export LHAPATH=/path/to/lhapdf/share/LHAPDF
./jewel-2.6.0-simple params.dat
```

`params.dat` has one `LABEL value` per line (lines starting with `#` are
comments) with the same labels and defaults as JEWEL 2.6.0 (`NEVENT`, `NJOB`,
`SQRTS`, `PTMIN`, `KEEPRECOILS`, `DOSUBTRACTION`, `MEDIUMPARAMS`, …; see the
manual). The medium parameter file named by `MEDIUMPARAMS` also keeps its
format (`TI`, `TAUI`, `CENTRMIN`, `CENTRMAX`, …). `params.example.dat` and
`medium.params.dat` from the distribution are included. Output files (log,
HepMC 2 IO_GenEvent, `splitint.dat`/`pdfs.dat`/`xsecs.dat` tables, `out.ran`)
have the same names and formats, and table files written by either version can
be read by the other. As in the Fortran version, parameter-file lines (and thus
file names) are limited to 80 characters.

Note (also true for the original): a run that reads existing table files is
statistically equivalent to, but not event-by-event identical with, a run that
has to compute them first.

## Layout

| path | content |
|---|---|
| `include/jewel/pythia6.hpp` | PYTHIA 6 common blocks and routines (`extern "C"`) |
| `include/jewel/commons.hpp` | JEWEL global state, one namespace per former COMMON block |
| `include/jewel/farray.hpp` | Fortran-layout arrays (1-based, column-major) |
| `include/jewel/medium.hpp` | medium model interface (`MediumSimple`, `MediumVac`) |
| `include/jewel/fortran_compat.hpp` | helpers reproducing Fortran arithmetic |
| `src/main_init.cpp` | program, parameters, PYTHIA initialisation, event loop |
| `src/event_setup.cpp` | interpretation of the PYTHIA event, colour, strings |
| `src/cascade.cpp` | medium-modified parton shower |
| `src/scattering.cpp` | scattering off medium partons, colour coherence |
| `src/kinematics.cpp` | scattering kinematics, momentum transfer |
| `src/sudakov_tables.cpp` | Sudakov factors, splitting functions, PDFs, cross sections, tables |
| `src/utilities_hepmc.cpp` | event-record utilities, HepMC output |
| `src/medium_simple.cpp`, `src/medium_vac.cpp` | medium models |
| `src/subtraction.cpp` | event-wide constituent subtraction (DOSUBTRACTION) |
| `fortran/` | modified PYTHIA 6.4.25 (unchanged) and a small Fortran I/O shim |
| `PORTING.md` | translation rules, including all bit-exactness details |

Known quirks of the original are preserved deliberately and marked
`NOTE(port)` in the source (for example the `lstart` logic in the
subtraction and a few reads just outside array bounds).
