# JEWEL in JETSCAPE: tests and cross-checks

This file collects every test and cross-check of the JEWEL final-state shower module (`Eloss/Jewel`,
branch `jewel` of X-SCAPE) done on ACCRE between 2026-10-19 and 2026-10-29: what was checked, how,
the result, and where the inputs and outputs live. The design is described in `README_JETSCAPE.md`.

Abbreviations: **A** = jewel-2.4.0 + medium-2D (I. Kolbe, arXiv:2303.14166,
github.com/isobelkolbe/jewel-2.4.0-2D); **B** = original Fortran JEWEL 2.6.0 + her medium-2D;
**C** = JETSCAPE + this module (C++ port of JEWEL 2.6.0). Paths: `$X` =
`/data/rke_group/kunnar1/jetscape-docker/X-SCAPE-jewel`, `$CMP` =
`/data/rke_group/kunnar1/jewel-2D-comparison`, `$PORT` = `/data/rke_group/Software/jewel-cpp-port`.

## 1. The vendored JEWEL is unchanged as a stand-alone program

| test | method | result |
|---|---|---|
| Regression of the vendored copy (after all library-mode changes) | `external_packages/jewel` built stand-alone like the port (gcc 14, no FMA), the three samples of `$PORT/runs/regression` rerun with their tables: pp vacuum, Pb+Pb, Pb+Pb with recoils + `DOSUBTRACTION`. Compare `events.hepmc` and `out.ran` | **byte-identical** for all three, both after the first round of changes (2026-10-20) and after the final ones (orphan centres, line origins; SLURM 14730165, 14733372); `$PORT/runs/regression-check-jetscape-10200929` |
| Library-mode changes are inert stand-alone | Every change is guarded by `cb::lib::active` or is a new virtual (`Medium::getneffprobe`) whose default returns the original `getneff` | follows from the regression above |
| gcc 11 compatibility (JETSCAPE image) | `std::format`/`std::println` replaced by `jewel/jfmt.hpp` (std when available, small formatter otherwise); all sources compile with `-std=c++2b` on gcc 11.5 | builds in `jetscape_base_v2.1.sif`; PYTHIA 6 built with the image's `f95` |

## 2. Physics of a single shower (library interface and module)

| test | method | result |
|---|---|---|
| Energy loss in a brick vs stand-alone JEWEL | `tests/test_jetscape_brick` (no JETSCAPE): 300 light quarks, E = 100 GeV, T = 0.3 GeV, L = 4 fm | **13.2 +- 0.6 GeV**, stand-alone JEWEL brick 13.16 +- 0.08 GeV (`HEAVYFLAVOUR.md` of the port) |
| Same through JETSCAPE | `config/jetscape_user_brick_Jewel.xml` (PythiaGun pT-hat 50-60 GeV, brick T = 0.3 GeV, 4 fm), 2000 events, parton level | energy in R < 0.4 around the leading parton: vacuum 0.731 +- 0.007, medium 0.613 +- 0.007 of its energy (gluons lose more than quarks) |
| Energy conservation per shower | shower partons + recoils - holes vs the initiating parton, every shower | exact (to the 6 digits of the ASCII writer, abs dE/E < 2e-5) in brick and Pb+Pb, after the fixes in section 3 |
| Three-momentum | same balance for px, py, pz | **not conserved**: the virtuality of a single parton is given at fixed energy (JEWEL distributes it between the two hard partons; JETSCAPE hands over one). Median 1 % of E; up to ~80 % for few-GeV ISR partons; identical in vacuum and medium |
| Recoils and holes | counts per shower; status codes 0 / 1 / -1 as in LBT | brick: 12.86 recoils = 12.86 holes per shower; Pb+Pb (MUSIC): 30.8 = 30.8 |
| Single returned parton (fix by the heavy-flavour session) | JetEnergyLoss keeps the *input* parton when a module returns exactly one parton, losing its medium modification (about half of the b quarks). Fix: update the input in place | b-quark brick (E = 10 / 20 GeV) now matches stand-alone JEWEL: 1.02/2.36 and 2.38/4.40 GeV vs 1.00/2.22 and 2.18/4.38 GeV |

## 3. Problems found and fixed in the module (library mode only)

| problem | how it showed up | fix and check |
|---|---|---|
| Undone coherent scatterings keep their stored scattering centres (**also in JEWEL 2.6.0 / 2.4.0**) | holes without recoils: in ~1 % of the showers shower + recoils - holes missed up to 90 % of the energy (debug switch `JEWEL_LIBDEBUG` showed the recorded recoil line was a shower line) | restore `nscatcen` when DOINSTATESCAT undoes the coherent scattering; afterwards no orphan holes and exact energy balance (brick, 100 events, same seeds: worst -89 % -> 1e-5) |
| COMPRESSEVENT clears all mother pointers | would label every shower parton as a recoil in events with > 20000 lines | origin (shower / recoil / deleted) kept per line and moved by COMPRESSEVENT; scattering-centre -> recoil line map updated |
| Scatterings undone by REJECTSPLITTING (status -1 lines) | possible holes without recoils | partons descending from status -1 lines and holes whose recoil is deleted are dropped |
| Rate veto bound in flowing media | "weight > 1" in GETDELTAT: 89 times in 77 Pb+Pb showers (up to 3.2) | flux factor gamma (1 - v.beta) also in the lower density bound: 15 times in 292 showers (max 1.78); `fluxmax` default 3 -> 6: **0** in 292 showers. gamma (1 - v.beta) can reach 22 = e^3.1 in a few cells; this produced no veto violation |
| Parameter file shared by parallel jobs | two jobs started together in one directory: one read a half-written file, ran with JEWEL defaults (no recoils) | each job reads a private copy (`run_param_file` keeps the last one); lookup tables named after a hash of all their inputs (stale tables can't be reused). Two jobs writing the same tables at once: byte-identical files (`cmp`) |
| Empty optional XML element crashed `InitTask` | `std::string(nullptr)` from tinyxml | `xml_text` helper |
| Module summary never written | `JetEnergyLoss::WriteTask` does not call sub-module writers | counters printed in `FinishTask` (showers, discarded, stored scatterings, holes dropped, cells above Tmax, rate queries above fluxmax, max flux) |
| Hydro queries outside the fluid | some fluid modules throw | exception -> no medium |

## 4. Full chain: TRENTo + free streaming + MUSIC + JEWEL

| test | result |
|---|---|
| Pb+Pb 5.02 TeV 0-10 %, 2 hydro events x 10 jet events (`build_music/test/pbpb_jewel.xml`) | runs (13 min, 16 cores); 30.8 recoils and holes per shower; energy balance per shower exact; memory 9.1 GB (MUSIC evolution kept in memory) |
| PythiaGun defaults | JETSCAPE's PythiaGun switches on **MPI and prompt photons** and drops partons below 2 GeV; stand-alone JEWEL (PPJJ) has neither MPI nor photons. For like-for-like comparisons: `PartonLevel:MPI = off`, `PromptPhoton:all = off` |
| Hadronisation of JEWEL output | colorless, `eCMforHadronization` 5020 (default in `jetscape_main.xml`): vacuum partons 1.10 TeV -> hadrons 1.74 TeV (beam partons added); with 0: 1.10 -> 1.10 TeV. Medium (with 0): positive partons 2.42 TeV -> hadrons 7.49 TeV, holes 1.33 -> 6.77 TeV. Colored (0): vacuum fine (1.10 -> 1.10), medium 2.42 TeV -> 3.54 TeV with 1.22 TeV left as partons. **Neither JETSCAPE hadroniser can be used with JEWEL's ~1000 soft recoils per event**; use parton-level output (hadron level needs JEWEL's own string building) |

## 5. Comparison with jewel-2.4.0-2D on one medium (`$CMP`, details in `$CMP/README.md`)

### Setup

* **Medium**: one TRENTo + free streaming + MUSIC Pb+Pb 5.02 TeV 0-10 % event written by
  `CustomModuleMediumDump` in her contour format (|x|, |y| <= 20 fm, 0.2 fm; tau = 0.2 ... 10 fm in
  0.2 fm steps) plus a 400k-sample binary-collision histogram (NCOLLHISTO). Checks: T(0,0) = 0.37 /
  0.25 / 0.175 GeV at tau = 1.2 / 4 / 8 fm; T > T_c until 8.6 fm; rms radius of the collisions 3.5 fm.
* **C reads the same files** through `CustomModuleHydroContours`, a line-by-line reproduction of
  medium-2D.f's interpolation (trilinear; linear rise before the first slice; zero flow at grid
  points below T_c; nothing beyond the largest x, y, tau; Bjorken flow added), and samples vertices
  from the same histogram (`CustomModuleNcollHisto`).
* **B** = original 2.6.0 Fortran with exactly her two core changes (longer file names; GETNEFF gets
  px, py), checked by `diff` against the release.
* **Identical settings**: pp at 5.02 TeV, CT14nlo (Pythia 8 through LHAPDF in C), QCD 2->2 with
  ISR, no MPI, pT-hat 50-500 GeV with a pT-hat^5 bias (event weights in all analyses), T_c 0.17,
  MDSCALEFAC 0.9, MDFACTOR 0.45, NF 3, ETAMAX 3.1, lifetime 13.4 fm (her GETLTIMEMAX), recoils kept,
  scattering centres written.
* **One analysis** (`$CMP/analysis/cmpana.cpp`): event-wide constituent subtraction of the thermal
  momenta (Milhano-Zapp, R_sub = 1), anti-kT R = 0.4, |y| < 2, dynamical grooming a = 0.1, 1, 2
  (z_G, theta_G, groomed mass), jet mass, R = 0.2 subjets (pT > 5 GeV), primary Lund plane.
* **Errors**: jackknife over seeds (`jackknife.py`); C medium and vacuum use the same seeds and hence
  the same hard events. **Equivalent energy loss** Delta from N_vac(>pT+Delta) = N_med(>pT), which
  removes the shape of the vacuum spectrum that a yield ratio mixes in.
* Statistics: A, B medium 25.6k events, vacuum 64k; C medium and vacuum 16k each.

### Results, parton level, subtracted

| | A (2.4.0) | B (2.6.0) | C (JETSCAPE) |
|---|---|---|---|
| thermal momenta per event | 582 (1253 GeV) | 593 (1113 GeV) | 629 (1118 GeV) |
| Delta at 100 / 150 / 250 GeV | 17.0+-0.7 / 18.4+-0.7 / 21.5+-1.3 | 14.4+-0.8 / 16.2+-0.7 / 16.4+-1.4 | 16.1+-0.6 / 17.6+-0.6 / 16.7+-1.0 |
| jet medium/vacuum 100-150 / 150-250 / 250-500 GeV (first 16 seeds) | 0.50 / 0.52 / 0.61 | 0.55 / 0.60 / 0.66 | 0.46 / 0.54 / 0.73 |
| mean M_jet med - vac (100-150 GeV) | 0.24+-0.19 | 0.49+-0.11 | 0.92+-0.23 |
| mean groomed M a = 0.1, med - vac | -0.89+-0.18 | -1.34+-0.13 | -0.98+-0.24 |
| mean z_G a = 0.1, vac -> med | 0.293 -> 0.205 | 0.246 -> 0.197 | 0.242 -> 0.199 |

### Cross-checks of the B-C comparison

| check | result |
|---|---|
| Vacuum | B and C agree (z_G 0.246 vs 0.242, jet mass 9.9 vs 9.8 GeV) |
| Without recoils (no subtraction) | agree within errors: Delta 20.7+-1.2 / 24.5+-1.3 / 28.4+-2.1 vs 22.7+-1.1 / 25.9+-1.3 / 26.0+-2.3 GeV; jet mass -4.25+-0.15 vs -4.00+-0.19 GeV; groomed mass -3.11+-0.15 vs -2.98+-0.17 GeV |
| Different pT dependence of the yield ratios | spectral shape: C's vacuum spectrum is slightly steeper at 100 GeV (n = 5.9 vs 5.5); in Delta B and C agree |
| Flux factor | C with medium-2D's formula (`flux_mode` 1): 613 instead of 636 scattering centres per event, Delta 16.0 / 17.0 / 15.6 GeV, jet mass +1.03+-0.19 GeV: unchanged within errors |
| B's orphan scattering centres | B balances energy exactly in vacuum but 17-18 % of the medium events are short by > 5 GeV (-2.9 GeV per event). "Bfix" (B with the undo fix of section 3): every event balances, jet observables unchanged (Delta 14.1 / 15.8 / 17.7 GeV, jet mass +0.31 GeV) |
| Remaining B-C difference | with recoils C loses 1-1.5 GeV more at 100-150 GeV (about 10 %, 1.5-1.8 sigma) and its jet mass grows by 0.4 GeV more (1.7 sigma); with half the statistics this looked like 2.5 sigma. **Not significant** |

### A vs B (JEWEL 2.4.0 vs 2.6.0)

| check | result |
|---|---|
| Energy loss | A loses more, parton level (Delta 17.0 / 18.4 / 21.5 vs 14.4 / 16.2 / 16.4 GeV) and hadron level (17.9 / 18.2 / 20.7 vs 14.9 / 15.8 / 19.7 GeV) |
| Vacuum partons | differ: 33 vs 29 final partons per event, larger z_G in A. Energy balance shows no double counting; JEWEL 2.4.0's parton-level output omits one beam remnant (11.5 % of the events miss up to 3 TeV at eta ~ 8; no effect on jets) |
| Vacuum hadrons | agree (z_G a = 0.1: 0.320 vs 0.317; jet mass 15.1 vs 14.7 GeV); each version has its own fragmentation tune (2.6: Professor LEP tune) |
| Hadron-level modification | same: jet mass -1.28+-0.13 vs -1.15+-0.10 GeV, z_G (a = 1) 0.208 -> 0.148 vs 0.193 -> 0.152, subjets 2.07 -> 2.53 vs 2.06 -> 2.55 |

## 6. Problems found outside this module

| where | problem | status |
|---|---|---|
| her medium-2D.f, PICKVTX (NCOLLHISTO) | `rndIndex = int(pyr(0)*ncollSum)` (all collisions) indexes `ncollNZ` (non-empty bins); `ncollpDist` (one row per collision) is built but unused. Emulated on our histogram: 1 % of the draws give real vertices, 70 % read past the array; 3 seeds segfaulted. With few collisions (He+Au): bins unweighted, index 0 read before the array in ~1/N_coll of the events | fixed in the comparison copies of A and B (first A/B medium runs kept as `$CMP/runs_bugvtx`); note to the author drafted (`$CMP/note_to_isobel.md`) |
| JEWEL 2.6.0 / 2.4.0, DOINSTATESCAT | undone coherent scattering keeps its scattering centres (section 3) | stand-alone port: switch `FIXCOHUNDO` (default F) drafted in `$CMP/port_switch/fixcohundo.patch`, tested in a scratch copy: F byte-identical to the regression samples, T balances every event (0 of 298 Pb+Pb events off by > 5 GeV, mean -0.003 GeV, vs 18 % without). Not applied to the port |
| JETSCAPE colorless hadronisation | `eCMforHadronization` default 5020 adds beam partons, also in AA | use 0 for AA (section 4) |
| JETSCAPE FluidDynamics::Init | requires a pre-equilibrium module for every hydro id except MUSIC and Brick | overridden in `CustomModuleHydroContours` |

## 7. Reproducing

* Stand-alone regression: `$PORT/runs/regression-check-jetscape-*/run.sbatch` (builds the vendored sources, compares byte for byte).
* Brick (library only): `test_jetscape_brick <params> <n> <E0> <T> <L> <seed>` built by `external_packages/jewel/CMakeLists.txt`.
* JETSCAPE brick: `config/jetscape_user_brick_Jewel.xml`; contour medium: `config/jetscape_user_Jewel_contours.xml`.
* Comparison: `$CMP/jetscape/dump_medium.xml` (medium), `$CMP/jewel/run_jewel.sbatch` (A, B, Bfix), `$CMP/jetscape/run_c.sbatch` (C; frozen binaries in `$CMP/bin_jetscape`), `$CMP/analysis/{cmpana.cpp,plot_cmp.py,jackknife.py}`; figures in `$CMP/figures`.
