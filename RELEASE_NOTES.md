# X-SCAPE contrib: changes since 2.1

Changes on the `contrib` branch ([`ef19d7f7`](https://github.com/JETSCAPE/X-SCAPE/commit/ef19d7f7), 2026-11-12) relative to the X-SCAPE 2.1 release (tag [`v2.1`](https://github.com/JETSCAPE/X-SCAPE/tree/v2.1), 2026-03-05). The third-party packages X-SCAPE fetches (iSS, 3dMCGlauber, MUSIC, trento, the LBT tables) are included. MUSIC4GPU's own changes are not listed here; see [jhputschke/MUSIC4GPU](https://github.com/jhputschke/MUSIC4GPU).

Fixes to code that was added after 2.1 are folded into the feature that added it.

- [Physics bug fixes](#physics-bug-fixes): change the physics output of existing setups.
- [Technical bug fixes](#technical-bug-fixes): crashes, hangs, memory, build and log noise. The output is unchanged.
- [Enhancements and new features](#enhancements-and-new-features): new or opt-in. Unless a row says otherwise, the defaults reproduce 2.1.

## Physics bug fixes

| Component | Fix | Effect in 2.1 | PR / commit |
| --- | --- | --- | --- |
| Liquefier | A free-streaming parton that the liquefier absorbs is marked as absorbed (status -11). Before, it kept status 0/1. | The parton was deposited in the medium and hadronized too. In two test jobs this double-counted 3% and 13% of the initial energy (issue [#156](https://github.com/JETSCAPE/X-SCAPE/issues/156)). | [#154](https://github.com/JETSCAPE/X-SCAPE/pull/154) |
| Liquefier | Each droplet is normalized on MUSIC's grid. Droplets at the hard vertex are deposited. | MUSIC received 0.004 to 6 times a droplet's 4-momentum. 40% of 0-10% Au+Au events were off by more than 20%. Vertex droplets got eta = 0/0 = NaN and were never deposited. On by default; `<normalize_on_hydro_grid>0` restores the old source. | [#155](https://github.com/JETSCAPE/X-SCAPE/pull/155) |
| CausalLiquefier | The causal diffusion is solved in the frame of the Bjorken flow at the droplet's eta_s (`<frame>bjorken`, the new default), not in the lab's Cartesian frame. | Away from midrapidity the deposit landed closer to midrapidity than the droplet (by about 0.4 units at \|eta_s\| = 1 and up to 0.9 units forward) and in a fraction of a cell. MUSIC reverted such cells, so deposited energy was lost: one test event kept 25% of it at tau = 2 fm/c. In Pb+Pb Z+jet the midrapidity wake came out 19% too large. **Behaviour change:** `<frame>lab` restores the old deposit; the two agree only at eta_s = 0. | [#189](https://github.com/JETSCAPE/X-SCAPE/pull/189) |
| Liquefier / MUSIC wrapper | The jet source's hydro dtau is set for every initial condition, not only without a pre-equilibrium module. | With string or pre-equilibrium initial conditions the droplet source was evaluated about 2.6 times too often. The hadronic liquefier's dtau stayed -1, so every droplet was deposited at every step. | [#146](https://github.com/JETSCAPE/X-SCAPE/pull/146) |
| Energy loss / MUSIC wrapper | `FluidDynamics::hydro_tau_0` and `hydro_tau_max` are initialized. MUSIC sets the energy-loss start time to its tau0, also when pre-equilibrium is kept in memory. | Matter, LBT and AdS/CFT read the start time as 0 in practice, so `<tStart>` had no effect. **Behaviour change:** medium-modified energy loss now starts at the hydro tau0. | [#174](https://github.com/JETSCAPE/X-SCAPE/pull/174) |
| Matter | No qhat table lookup at a negative index (ISR partons start at t < 0). | Garbage qhat values, results that changed between runs with the same seed, and segfaults (O+O with ISR). | [`8c975d37`](https://github.com/JETSCAPE/X-SCAPE/commit/8c975d37) ([#171](https://github.com/JETSCAPE/X-SCAPE/pull/171)) |
| Colour tags (Matter, iMATTER) | Colour tags come from one counter per event. | Matter's tags overlapped between showers: 85% of O+O events and 67-83% of pp Z+jet events had a tag used twice. Pythia joined colour lines across showers without a warning. Affects colour-aware hadronization only (iColored, Colored); kinematics are unchanged. | [#181](https://github.com/JETSCAPE/X-SCAPE/pull/181) |
| iColoredHadronization | Holes and liquefier partons are no longer handed to Pythia as ordinary partons. The tag-pairing loop no longer reads past the end of its vector. `pythia.next()` is checked. | A failed fragmentation wrote the input partons out as hadrons. Now the event throws. | [`ad1095be`](https://github.com/JETSCAPE/X-SCAPE/commit/ad1095be) ([#182](https://github.com/JETSCAPE/X-SCAPE/pull/182)) |
| MUSIC wrapper | MUSIC's `reRunHydro` flag is reset for every event, with a warning when the grid edge is hit. | After one event touched the transverse grid edge, every later event of that MUSIC instance stopped at its first freeze-out check, silently. | [`ca8dd84a`](https://github.com/JETSCAPE/X-SCAPE/commit/ca8dd84a) ([#141](https://github.com/JETSCAPE/X-SCAPE/pull/141)) |
| MUSIC wrapper | `EOS_to_use` is written before MUSIC is constructed. | The first MUSIC instance in a process ran with whatever EOS the input file still held. | [`a1753327`](https://github.com/JETSCAPE/X-SCAPE/commit/a1753327) |
| SurfaceFinder | The flow is built as u = gamma (1, v). Lattice cells are counted with `floor(L/d + 1e-4)`. | u.u was 1 - v^4 (13% short at v = 0.7). The last cell in tau, x, y and eta was dropped: on a 4x4x2 grid the surface covered 26 of 32 fm^2. HybridHadronization uses only ratios of u, so u did not change its output. | [`88099402`](https://github.com/JETSCAPE/X-SCAPE/commit/88099402) |
| MUSIC (CPU) | Bulk Pi in the eta-direction flux: `W_eta_p` / `W_eta_m` take `Pi_p` / `Pi_m` (they were swapped). | Wrong bulk-viscous contribution to the eta flux. | [MUSIC `9a304fb`](https://github.com/MUSIC-fluid/MUSIC/commit/9a304fb) ([#138](https://github.com/JETSCAPE/X-SCAPE/pull/138)) |
| MUSIC (CPU) | Charge source terms (`turn_on_QS`) use the coordinate factor instead of tau. `CoorType` can be set from X-SCAPE. | Wrong charge sources in Cartesian mode; Cartesian mode could not be selected from X-SCAPE. | [MUSIC `3bc760d`](https://github.com/MUSIC-fluid/MUSIC/commit/3bc760d) ([#138](https://github.com/JETSCAPE/X-SCAPE/pull/138)) |

## Technical bug fixes

| Component | Fix | Effect in 2.1 | PR / commit |
| --- | --- | --- | --- |
| 3dMCGlauber | The draws for the target's hard-quark position are capped. Pinned from the [jhputschke/3dMCGlauber](https://github.com/jhputschke/3dMCGlauber) fork; upstream PR [chunshen1987/3dMCGlauber#7](https://github.com/chunshen1987/3dMCGlauber/pull/7) was merged on 2026-11-04. | Endless loop in about 1 in 500-1300 O+O events with ISR (issue [#169](https://github.com/JETSCAPE/X-SCAPE/issues/169)). | [`21a06787`](https://github.com/JETSCAPE/X-SCAPE/commit/21a06787), [3dMCGlauber `b99d680`](https://github.com/jhputschke/3dMCGlauber/commit/b99d680) |
| EvolutionHistory | `get` on an empty history returns a vacuum cell. | `std::out_of_range` stopped the job, in about 5% of O+O events (events without a medium). | [`be68bc85`](https://github.com/JETSCAPE/X-SCAPE/commit/be68bc85) ([#165](https://github.com/JETSCAPE/X-SCAPE/pull/165)) |
| iSS wrapper | An event with a hydro but no freeze-out surface gives no hadrons and a warning. | `exit(1)`. | [`b8487c42`](https://github.com/JETSCAPE/X-SCAPE/commit/b8487c42) ([#164](https://github.com/JETSCAPE/X-SCAPE/pull/164)) |
| iSS | Decay channels and sampled-hadron vectors are freed. | About 1,900 leaked blocks per event. | [iSS `98deaff`](https://github.com/jhputschke/iSS/commit/98deaff) ([#171](https://github.com/JETSCAPE/X-SCAPE/pull/171)) |
| iSS | One hadron object is no longer leaked per sampled hadron (part of the yield speed-up). | Memory grew with every sampled hadron. | [iSS `01f7cf9`](https://github.com/jhputschke/iSS/commit/01f7cf9) ([#149](https://github.com/JETSCAPE/X-SCAPE/pull/149)) |
| MUSIC wrapper | The old framework copy of the evolution is freed before a larger one is allocated. | Peak memory of a 0-10% Au+Au jet job was 20 GiB; now 14 GiB. | [#160](https://github.com/JETSCAPE/X-SCAPE/pull/160) |
| PythiaIsrMUSIC | The example adds the pre-equilibrium module it needs. iMATTER finds the PDFs via `PYTHIA8DATA` (conda). | Crash in MUSIC's `InitArena`; iMATTER needed `$PYTHIA8DIR`. | [`5ae7d2a1`](https://github.com/JETSCAPE/X-SCAPE/commit/5ae7d2a1) ([#164](https://github.com/JETSCAPE/X-SCAPE/pull/164)) |
| Matter | No qhat lookup outside the forward light cone. | About 500 "Third instance" warnings per O+O event. | [`5900ebd0`](https://github.com/JETSCAPE/X-SCAPE/commit/5900ebd0) ([#164](https://github.com/JETSCAPE/X-SCAPE/pull/164)) |
| Liquefier | The energy-momentum check warns only above 5% of the incoming energy. | Warned at 1e-10 GeV: 19 warnings for 25 droplets in one Au+Au event. | [`d679b204`](https://github.com/JETSCAPE/X-SCAPE/commit/d679b204) ([#142](https://github.com/JETSCAPE/X-SCAPE/pull/142)) |
| Build (CMake) | OpenMP is enabled per language. | On macOS `OPENMP_FOUND` was false: the framework was built without OpenMP and SurfaceFinder and ThermPtnSampler ran on one thread. | [`3955cef6`](https://github.com/JETSCAPE/X-SCAPE/commit/3955cef6) |
| Build (CMake) | libiSS and lib3dMCGlb are linked with `-Bsymbolic`. | iSS and MUSIC4GPU both export a different `pretty_ostream`. Built together (CUDA, PyJetscape), iSS aborted in `malloc` on its first allocation. | [`f25ba028`](https://github.com/JETSCAPE/X-SCAPE/commit/f25ba028), [`73996c07`](https://github.com/JETSCAPE/X-SCAPE/commit/73996c07) |
| Build | `ThermPtnSampler.h` includes `<unordered_map>`. | Did not compile with some compilers. | [`16a658e9`](https://github.com/JETSCAPE/X-SCAPE/commit/16a658e9) |
| trento | `hdf5_utils.cxx` fixed. | Did not compile with newer compilers (Linux arm). | [`9980458c`](https://github.com/JETSCAPE/X-SCAPE/commit/9980458c) |
| LBT tables | `get_lbtTab.sh` runs with bash 3.2 and joins the parts in numerical order. | Failed under macOS's `/bin/bash` (`declare -A`). | [#186](https://github.com/JETSCAPE/X-SCAPE/pull/186) |
| JetScapeWriter | The destructor is defined out of line, so its type information is emitted once. | `dynamic_cast` to a writer failed across the Python module boundary. | [`9b3641b5`](https://github.com/JETSCAPE/X-SCAPE/commit/9b3641b5) |

## Enhancements and new features

### Event generation

| Component | Change | PR |
| --- | --- | --- |
| PythiaGun | `<pTHatBins>`: several pTHat windows per run, each with its own Pythia, seed and cross section. | [#152](https://github.com/JETSCAPE/X-SCAPE/pull/152) |
| PythiaGun | `<partonYMax>`: rapidity cut on the handed-over partons; the cross section is scaled by the acceptance. | [#153](https://github.com/JETSCAPE/X-SCAPE/pull/153) |
| PythiaGun | Z+jet: `<trigger>Z` with lepton acceptance cuts, `<minPartons>`. Final-state leptons are handed over and pass through energy loss like photons. | [#173](https://github.com/JETSCAPE/X-SCAPE/pull/173) |
| PythiaGun | `<extraHard>`: the event's other hard scatterings, each with its own stream; minimum-bias mode. | [#187](https://github.com/JETSCAPE/X-SCAPE/pull/187) |
| ISR | `<Hard><ISR>` in the XML task list builds the PythiaIsrMUSIC chain. `<JetHadronization><name>iColored`. PythiaIsrGun gets `<pTHatBins>`. The subtracted hard momentum is reported, with a warning for `Subtract_hard_momentum` other than 1. | [#164](https://github.com/JETSCAPE/X-SCAPE/pull/164) |
| PythiaIsrGun | `<partonYMax>` / `<partonYMode>`: rapidity cut on the two hard partons. | [#165](https://github.com/JETSCAPE/X-SCAPE/pull/165) |
| PGun | `<rapidity>` and `<phi>` of the parton (both 0 as before); `<points>`: several (pT, y) points in one run, event i takes point i mod K. | [#189](https://github.com/JETSCAPE/X-SCAPE/pull/189) |

### Energy loss and hadronization

| Component | Change | PR |
| --- | --- | --- |
| Matter, LBT, iColoredHadronization | Colour flow: Matter's recoils and holes keep their colour flow, LBT partons carry colour tags, and iColoredHadronization closes colour lines that end in the medium (`<colour_closure>`, default `reconnect`). See `README_ColourFlow.md`. | [#182](https://github.com/JETSCAPE/X-SCAPE/pull/182) |
| Energy loss (Matter, LBT) | `<Eloss><nThreads>`: the showers of an event run in parallel on OpenMP threads (N threads, or -1 for `OMP_NUM_THREADS`). Off by default (`0`, output unchanged). When on, every shower gets random streams and colour tags of its own and droplets are flushed in the serial order: the output does not depend on the thread count, but a seed gives other (equally valid) showers than with `0`. Runs serially, with a warning, unless all energy-loss modules declare support (Matter, LBT), showers run per event and there is no hadronic medium; ISR stays serial. O+O with ISR at 4 threads: energy loss 18 to 6.2 s per 40 events, jobs 13% faster; no gain in Au+Au jet jobs. | [#191](https://github.com/JETSCAPE/X-SCAPE/pull/191) |
| Matter | `<tStart_at_point>`: apply the start-time cut at the point where it is used (default off). | [#183](https://github.com/JETSCAPE/X-SCAPE/pull/183), [#184](https://github.com/JETSCAPE/X-SCAPE/pull/184) |
| LBT | `LBT_TABLES_PATH` / XML: table location configurable. `LBT_TABLE_CACHE`: tables in one block, shared read-only between jobs. | [#130](https://github.com/JETSCAPE/X-SCAPE/pull/130), [#163](https://github.com/JETSCAPE/X-SCAPE/pull/163) |
| Hadronization | Reproducible re-hadronization: `SetNextRandomSeed` / `GetLastRandomSeed` for iSS and Colorless, `<reseed_per_event>` for Colorless and iColored, and iColored's `<store_input>` / `HadronizeInput`. | [#148](https://github.com/JETSCAPE/X-SCAPE/pull/148), [#172](https://github.com/JETSCAPE/X-SCAPE/pull/172) |
| Liquefier | Only droplets that can deposit in the current step are evaluated. When none can, MUSIC skips the jet source (output unchanged). | [#141](https://github.com/JETSCAPE/X-SCAPE/pull/141), [#144](https://github.com/JETSCAPE/X-SCAPE/pull/144) |
| Liquefier | The jet source and CausalLiquefier follow MUSIC's late time step (`<late_delta_tau>`, `<late_tau_switch>`). | [#171](https://github.com/JETSCAPE/X-SCAPE/pull/171) |
| CausalLiquefier | `<frame>fluid`: the diffusion in the rest frame of the background fluid where the droplet was made (the droplet stores that cell's velocity and energy density). Opt-in; a droplet without a usable fluid takes the Bjorken frame. | [#189](https://github.com/JETSCAPE/X-SCAPE/pull/189) |

### Hydro (MUSIC)

| Component | Change | PR |
| --- | --- | --- |
| MUSIC wrapper | `<freeze_out_surface>`: switch off the surface export, globally or per instance. | [#142](https://github.com/JETSCAPE/X-SCAPE/pull/142) |
| MUSIC wrapper | `<slim_bulk_info>`: a 24-byte-per-cell framework copy for energy loss (peak memory 14.2 to 8.7 GiB in a 0-10% Au+Au jet job, output identical). | [#162](https://github.com/JETSCAPE/X-SCAPE/pull/162) |
| MUSIC wrapper | The framework copy is sized once and filled in parallel. | [#143](https://github.com/JETSCAPE/X-SCAPE/pull/143) |
| MUSIC wrapper | Two-stage runs (background, then a jet leg): `set_store_evolution`, `set_checkpoint_tau` / `set_restart_from`, per-instance XML and `<restart_from>` from `runJetscape`, `<output_momentum_anisotropy_vs_etas>`. The background's checkpoint is taken a deposit half width before the liquefier's `tau_delay`, so jet legs with a droplet from the hard vertex restart from it too (output unchanged). | [#166](https://github.com/JETSCAPE/X-SCAPE/pull/166), [#168](https://github.com/JETSCAPE/X-SCAPE/pull/168), [#180](https://github.com/JETSCAPE/X-SCAPE/pull/180), [#190](https://github.com/JETSCAPE/X-SCAPE/pull/190) |
| Liquefier / MUSIC wrapper | The jet source of a hydro step is filled from each droplet's light cone instead of asking at every cell (`LiquefierBase::get_source_on_grid`; used by MUSIC4GPU, bit-identical): 12 to 1.1 ms per fill. With the checkpoint change a Pb+Pb Z+jet job (1 background, 6 jets) takes 118 instead of 146 s. | [#190](https://github.com/JETSCAPE/X-SCAPE/pull/190) |
| Soft particlization | Particlization from stored evolutions (SurfaceFinder on the stored history). | [`88099402`](https://github.com/JETSCAPE/X-SCAPE/commit/88099402) |
| MUSIC (CPU) | Two source-term pointers; a flag for dynamic grid padding. Pin `cfcc264` to `cee9460`. | [#138](https://github.com/JETSCAPE/X-SCAPE/pull/138) |
| MUSIC4GPU | Optional GPU hydro build (`get_music4gpu.sh`, `USE_CUDA` / `USE_METAL`). MUSIC4GPU's own changes are not listed here. | [`51e58ce7`](https://github.com/JETSCAPE/X-SCAPE/commit/51e58ce7) |

### Particlization (iSS)

iSS is now pinned from the fork [jhputschke/iSS](https://github.com/jhputschke/iSS) (`sampling_omp` [`0be6daf`](https://github.com/jhputschke/iSS/commit/0be6daf), was chunshen1987/iSS `XSCAPE` `d242555`). The yield speed-up and the correlated sampling were merged upstream on 2026-11-04 ([chunshen1987/iSS#4](https://github.com/chunshen1987/iSS/pull/4)); the SMASH-ready hadrons and the threaded sampling are only in the fork.

| Component | Change | PR |
| --- | --- | --- |
| iSS | Speed-ups, hadrons bit-identical: yields once per (mu, T) and delta f coefficients once per cell, 8 species per pass over packed cell arrays, a tabulated momentum CDF. About 56 to 15 s per event. | [#149](https://github.com/JETSCAPE/X-SCAPE/pull/149), [#170](https://github.com/JETSCAPE/X-SCAPE/pull/170), [#171](https://github.com/JETSCAPE/X-SCAPE/pull/171) |
| iSS | Correlated sampling: random numbers addressed by the cell (Philox), switched in the XML. | [#151](https://github.com/JETSCAPE/X-SCAPE/pull/151) |
| iSS | Correlated sampling and the surface preparation run on OpenMP threads (`OMP_NUM_THREADS`), hadrons bit-identical at any thread count: 2.4 to 0.43 s per surface at 10 threads (1.1 million cells, 100 oversamples). Uncorrelated sampling stays on one thread. | [`cbc992e2`](https://github.com/JETSCAPE/X-SCAPE/commit/cbc992e2), [iSS `0be6daf`](https://github.com/jhputschke/iSS/commit/0be6daf) |
| SoftParticlization | Opt-in compact hadron output (iSS fills flat arrays); `SetNumberOfSamples`. | [#148](https://github.com/JETSCAPE/X-SCAPE/pull/148), [#150](https://github.com/JETSCAPE/X-SCAPE/pull/150) |
| iSS | SMASH-ready hadrons: `<decay_smash_list>`, `<keep_primaries>`; `SampleCells` for hadrons-only jobs. | [#172](https://github.com/JETSCAPE/X-SCAPE/pull/172) |

### Output, build and running

| Component | Change | PR |
| --- | --- | --- |
| Writers | `RootBulkWriter` and `FastRootBulkWriter` (hydro evolution to ROOT, e.g. for FNO training). Built only when ROOT is found. | [#129](https://github.com/JETSCAPE/X-SCAPE/pull/129), [#131](https://github.com/JETSCAPE/X-SCAPE/pull/131), [#157](https://github.com/JETSCAPE/X-SCAPE/pull/157), [#158](https://github.com/JETSCAPE/X-SCAPE/pull/158) |
| js-contrib | js-contrib built with X-SCAPE (`get_js_contrib.sh`, CMake), including the PyJetscape Python interface. `JetScapePerEvent`: the event loop driven from outside. | [`21ae1992`](https://github.com/JETSCAPE/X-SCAPE/commit/21ae1992), [`45eba887`](https://github.com/JETSCAPE/X-SCAPE/commit/45eba887) |
| Install | Full install (GNUInstallDirs, `find_package(XSCAPE)`). `XSCAPE_DATA_DIR` and `examples/run_in_workdir.sh` run many jobs from one build. | [#130](https://github.com/JETSCAPE/X-SCAPE/pull/130), [#132](https://github.com/JETSCAPE/X-SCAPE/pull/132), [#188](https://github.com/JETSCAPE/X-SCAPE/pull/188) |
| Install | macOS wheel: install component `xscape_runtime`, relative rpaths, Metal shaders. | [#185](https://github.com/JETSCAPE/X-SCAPE/pull/185) |
| Build | `configure.sh`: a text menu for the build options (needs `dialog`). | [`6fdccb63`](https://github.com/JETSCAPE/X-SCAPE/commit/6fdccb63) |
| 3dMCGlauber wrapper | `mcglauber.input` location configurable (`XSCAPE_DATA_DIR` / XML). | [#130](https://github.com/JETSCAPE/X-SCAPE/pull/130) |
| Configs | Finite-volume solver vs MUSIC comparison configs and run script. | [#136](https://github.com/JETSCAPE/X-SCAPE/pull/136) |
