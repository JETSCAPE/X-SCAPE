# JEWEL as a final-state shower in X-SCAPE

This directory holds the C++23 port of JEWEL 2.6.0 (with the modified PYTHIA
6.4.25 it uses) and an interface that lets JETSCAPE use JEWEL as its
final-state shower: the module `Jewel` (`src/jet/Jewel.{h,cc}`). JEWEL replaces
MATTER and LBT. The rest of the chain (initial state, PythiaGun, hydro,
hadronization, writers) is unchanged.

## Build

```bash
# where gfortran exists (not needed if the JETSCAPE build has a Fortran compiler)
external_packages/jewel/build_pythia6.sh
# then, e.g. in the JETSCAPE docker image
mkdir build && cd build
cmake .. -DUSE_JEWEL=ON [-DUSE_MUSIC=ON -DUSE_ISS=ON ...]
make -j8
```

The JEWEL library is compiled as C++23 (`-std=c++2b` with GCC 11) and linked
into `libJetScape`; `include/jewel/jfmt.hpp` stands in for `std::format` where
the standard library lacks it. PYTHIA 6 needs `libgfortran` and the Fortran
interface of LHAPDF at link time (both are in the docker image). No PDFs are
used at run time: JETSCAPE's PythiaGun makes the hard process.

## Use

```xml
<Eloss>
  <maxT> 20 </maxT>
  <Jewel>
    <param_file>my_jewel_params.dat</param_file>
    <Tmax> 0.6 </Tmax>
  </Jewel>
</Eloss>
```

All options and their defaults are in `config/jetscape_main.xml`; an example
is `config/jetscape_user_brick_Jewel.xml`. JEWEL parameters such as `Q0`,
`LPS`, `LQCD`, `ANGORD`, `RECMODE`, `SCATRECOIL`, `HEAVYFLAVOUR`, `HQDS` are
set in `param_file`, which has the format of the JEWEL parameter file. The
module appends `SQRTS`, `PTMAX`, `LOGFILE` and `KEEPRECOILS` from the XML and
writes the result to `run_param_file`. The lookup tables are computed at the
first start in the run directory (about 20 s) and are read in later runs. They
depend on the JEWEL parameters and on the medium bounds (`Tmax`, `Tc`,
`mdfactor`, `mdscalefac`, `nf`, `etamax`), so the module names them
`jewel_{splitint,pdfs,xsecs}_<hash>.dat` after a hash of all of these; changed
settings get new tables.

Parallel jobs: each job reads its own copy of the parameter file, but the JEWEL
log and the tables are shared within a directory. Run the first job alone (or
copy the tables in) before starting many jobs in one directory, so that no job
reads a table that another one is still writing; or give each job its own
directory.

## Heavy flavour

The JEWEL copy here contains the heavy-flavour extension of the port (branch `heavy-flavour` of
VandyHENP/jewel_cpp_port; physics and validation in its HEAVYFLAVOUR.md and VALIDATION.md): massive c and b in
the shower and the medium interactions (`HEAVYFLAVOUR T`), and optional Langevin transport of heavy quarks
(`HQDS` = 2πT·D_s, `HQP0` for a fall-off with momentum). Set them in `param_file`. In the module:

* the incoming heavy quark gets a mass-aware virtuality and z limits;
* the diffusion takes T and the full fluid velocity from `GetHydroCellSignal`; each step is in the local fluid
  rest frame, until T < T_c;
* a diffused heavy quark is handed to the framework at the point where it leaves the QGP (not at its last
  interaction), so that HybridHadronization can recombine it with thermal partons from the hypersurface;
* JEWEL's own recombination (`HQRECO`, through its strings) is switched off in library mode: hadronisation is
  the framework's.

Brick check (PGun b quark at y = 0, T = 0.3 GeV, 0 < t < 4 fm), energy lost by the b quark, JETSCAPE vs
stand-alone JEWEL: 1.02 vs 1.00 GeV (E = 10 GeV), 2.38 vs 2.18 GeV (20 GeV); with `HQDS 4`: 2.36 vs 2.22 and
4.40 vs 4.38 GeV. Pb+Pb 0–10 % with TRENTo + free streaming + MUSIC + JEWEL (`HQDS 4`) + HybridHadronization runs
and produces D and B mesons.

**Single returned parton.** `JetEnergyLoss` treats a module output of exactly one parton as free streaming and
puts the *input* parton into the shower graph. The module therefore copies a single returned parton into its
input (`pIn[i] = pOut.back()`), as MARTINI and LBT update the parton in place. Without this, a parton that JEWEL
returns alone lost its medium modification in the output: rare for light partons, about half of the b quarks
(dead cone).

## How it works

* **One-shot shower.** JEWEL does not evolve in time steps. On the first call
  for a hard parton the module runs JEWEL's complete cascade (`makecascade`)
  for it and returns all final partons at once. Every returned parton carries a
  `JewelUserInfo`, and later time steps leave it alone.
* **Starting virtuality.** The parton's maximum virtuality is `qmax_factor * pT`.
  JEWEL uses p_T of the hard scattering. The virtuality is sampled with JEWEL's
  Sudakov (`getmass`), as in `findvirtualities`. JEWEL rescales the two hard
  partons together in their CM frame, but JETSCAPE showers each parton on its
  own. So here the energy is kept and the three-momentum is reduced. There is
  also no colour partner, so the first emission has no angular constraint (as
  for JEWEL's MPI partons).
  The shower then conserves energy exactly but not the three-momentum: in the
  brick test (PythiaGun, p̂_T 50-60 GeV) |Δp|/E of a shower (with recoils and
  holes) has a median of 1 % and reaches about 65 % for soft (few GeV) ISR
  partons, the same in vacuum and in the medium.
* **Medium (`MediumJetscape`, in `src/jetscape_api.cpp`).** Temperature and
  flow come from `GetHydroCellSignal`, as in LBT. Scattering centres are thermal
  quarks or gluons with masses m_s (quarks) and m_D (gluons), sampled in the
  fluid rest frame and boosted with the cell velocity, as LBT does for its
  thermal partons. The Debye mass, the density, the quark/gluon fractions and
  the cut T < T_c follow `medium-simple.f`. The Bjorken flux factor of
  `medium-simple.f` (density / cosh η) is replaced by the general rate factor
  γ(1 − v·β_parton). The medium hook `Medium::getneffprobe` makes this possible;
  the media of the JEWEL distribution keep their results bit for bit. Before
  `tStart` there is no medium (`pre_tau0 0`, as in LBT); `pre_tau0 1` gives
  the linear rise T(τ0)·τ/τ0 of `medium-simple.f`. There are no interactions
  after `lifetime` (default: `Eloss/maxT`).
* **Output and status codes (as in LBT).** Shower partons have status 0. Medium
  recoils (`KEEPRECOILS T`) have status 1. The thermal momenta of the scattering
  centres before the scattering ("holes", for the 4-momentum subtraction) have
  status −1. A parton descended from a recoiled scattering centre counts as a
  recoil. JEWEL's colour indices are mapped into the colour range that
  PythiaGun reserved for the parton. The colour connection between the two
  hard partons is not known to JEWEL, so colourless hadronization is the
  natural choice.
* **Limitations.** JEWEL has global state (PYTHIA 6 common blocks), so the
  showers run one at a time (mutex). `Tmax` is a hard bound used by the veto
  algorithm. Cells above it are counted (written as a comment by the writers),
  and JEWEL logs "weight > 1". The soft-recoil cut `RECSOFTCUT` still boosts
  with the Bjorken rapidity (unused with the default `RECSOFTCUT 0`).
  Event-wide constituent subtraction (`DOSUBTRACTION`) is switched off; use the
  holes instead.
* **Checks (2026-10-20).** The JEWEL copy here, built stand-alone, reproduces the
  regression samples of the port (pp, Pb+Pb, Pb+Pb with recoils and
  subtraction) byte for byte. Brick (T = 0.3 GeV, L = 4 fm, 100 GeV quark): energy
  given to the medium 13.2 +- 0.6 GeV, as stand-alone JEWEL (13.16 +- 0.08).
  Pb+Pb 0-10 % with TRENTo + free streaming + MUSIC runs (about 30 recoils and
  holes per shower; shower + recoils - holes conserves the energy of every
  shower). With `fluxmax` 6 the rate veto is never exceeded in this test
  ("weight > 1" does not occur), although γ(1 − v·β) itself reaches about 19
  in a few cells (counted in the summary that the module prints at the end).
  MUSIC with the evolution in memory needs about 9 GB per job.
* **Fixed only in library mode** (the stand-alone program keeps the behaviour
  of JEWEL 2.6.0, so that it stays identical to the original):
  - when DOINSTATESCAT undoes a coherent scattering it resets the event record
    but not the stored scattering centres; they would become holes without
    recoils (stand-alone `DOSUBTRACTION` subtracts them);
  - COMPRESSEVENT clears the mother pointers, so the origin of each line
    (shower / recoil) is kept separately;
  - scatterings undone by REJECTSPLITTING (status -1) give neither recoil nor
    hole.

## (2+1)D contour media and the comparison with jewel-2.4.0-2D

Three small modules make it possible to run on exactly the medium of the JEWEL (2+1)D interface
(I. Kolbe, arXiv:2303.14166, github.com/isobelkolbe/jewel-2.4.0-2D):

* `CustomModuleMediumDump` (Eloss, `src/jet/MediumContourDump.*`) writes the fluid of an event at
  z = 0 as `Tcontour<tau>` ("x y T") and `Vcontour<tau>` ("x y v_x v_y") files, plus a histogram of
  the binary-collision positions (`ncoll_histo.dat`, her NCOLLHISTO format).
* `CustomModuleHydroContours` (Hydro, `src/hydro/HydroContours.*`) reads such files and interpolates
  them exactly as her medium-2D.f (trilinear in tau, x, y; linear rise before the first slice;
  T < T_c is no medium), with Bjorken flow added: v = (v_x / cosh eta_s, v_y / cosh eta_s, tanh eta_s).
* `CustomModuleNcollHisto` (IS, `src/initialstate/NcollHisto.*`) samples the hard vertices from the
  histogram (bin ~ N_coll, uniform within the bin).

`<flux_mode> 1 </flux_mode>` replaces the exact rate factor gamma (1 - v.beta) by medium-2D's
(cosh y_T - sinh y_T cos theta_T) / cosh eta_s. Example: `config/jetscape_user_Jewel_contours.xml`.

Comparison on one TRENTo+MUSIC Pb+Pb 5.02 TeV 0-10 % event (parton level, constituent-subtracted,
R = 0.4 jets, 25k/16k events): JETSCAPE + JEWEL agrees with stand-alone JEWEL 2.6 + medium-2D in
vacuum and, without recoils, in energy loss and all substructure observables. With recoils it loses
1-1.5 GeV more at 100-150 GeV (about 10 %) and its jet mass after subtraction grows by 0.9 +- 0.2
instead of 0.5 +- 0.1 GeV, both below 2 sigma. The flux formula changes the results by less than
their errors. Details: /data/rke_group/kunnar1/jewel-2D-comparison/README.md.
