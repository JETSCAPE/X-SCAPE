# MUSIC4GPU in X-SCAPE

MUSIC4GPU ([jhputschke/MUSIC4GPU](https://github.com/jhputschke/MUSIC4GPU), branch `XSCAPE`)
is a GPU port of MUSIC (CUDA on NVIDIA, Metal on Apple Silicon). This page covers the build,
the OpenMP defaults of GPU builds, the settings X-SCAPE's `MusicWrapper` reads for it, and the
two-stage hydro (background, then a jet leg restarted from it) with `runJetscape`. The physics
of the GPU port, its conservation fixes and the active region are described on the js-contrib
website: [GPU acceleration](https://jhputschke.github.io/js-contrib/about/gpu/).

## Build

A GPU-capable drop-in replacement, `music4gpu`, can be used instead of the CPU MUSIC
package (one or the other, never both). Download it with `./get_music4gpu.sh`, then
select it via the GPU backend flag (which both picks the package and enables that GPU
backend):

```bash
    cmake -DUSE_MUSIC=ON -DUSE_CUDA=ON ..   # NVIDIA / CUDA
    cmake -DUSE_MUSIC=ON -DUSE_METAL=ON ..  # macOS / Apple Silicon (Metal)
```

Plain `-DUSE_MUSIC=ON` (no GPU flag) keeps the original CPU `music` package.

## OpenMP defaults for GPU builds

The `music4gpu` backends run much of the
hydro on the GPU while the host does short OpenMP repack/output loops between
GPU waits. To stop idle OpenMP threads from busy-spinning between those loops —
which on a many-core host (e.g. a 20-core GB10 Grace) can pin every core and
inflate CPU time ~10–20× with no wall-clock benefit — GPU builds set
`OMP_WAIT_POLICY=passive` as a **default** at startup, so idle threads sleep
instead of spin. This is process-wide (it also affects 3DGlauber, iSS, …) and
only changes idle-thread behavior, so it is safe in practice. It is only a
default — anything you export yourself wins — and is fully revertible:

```bash
export OMP_WAIT_POLICY=active   # restore the OpenMP default (spinning)
export MUSIC_OMP_DEFAULTS=0     # disable the music4gpu OpenMP defaults entirely
export OMP_NUM_THREADS=8        # also worth capping on many-core hosts
```

(See MUSIC4GPU's [`docs/PORT_GPU.md`](https://github.com/jhputschke/MUSIC4GPU/blob/XSCAPE/docs/PORT_GPU.md)
§9.10 for the full rationale; locally it is
`external_packages/music4gpu/docs/PORT_GPU.md`.)

## music4gpu settings in X-SCAPE

Beyond the usual MUSIC parameters, `MusicWrapper` reads
or sets these. CPU MUSIC ignores all of them.

| setting | default | effect |
|---|---|---|
| `<Hydro><MUSIC><freeze_out_surface>` | 1 | 0: no freeze-out surface, MUSIC stops on max(e) < e_fo at the same step; per instance in its own `<Hydro><MUSIC>` block, or `MpiMusic::set_freeze_out_surface` |
| `<Hydro><MUSIC><late_delta_tau>`, `<late_tau_switch>` | 0, 0 (also when missing) | a larger MUSIC time step from `late_tau_switch` on (music4gpu `late_delta_tau`), with the freeze-out checks and the stored frames kept at their spacing; the jet source follows the step (`HydroSourceJETSCAPE::hydro_dtau_changed`, `CausalLiquefier`'s deposition windows tile the axis across the switch). The production XMLs in js-contrib set 0.05 fm/c from 3 fm/c |
| `<Hydro><MUSIC><output_momentum_anisotropy_vs_etas>` | 1 (also when missing) | 0: no `momentum_anisotropy_*` / `eccentricities_*` `.dat` files, which skips a serial pass over the grid at four steps of every run |
| `<Hydro><MUSIC><store_evolution>` | 1 (also when missing) | 0: the instance keeps its grid metadata and frame count (`bulk_info`) but stores no evolution frames (for a leg nobody reads); per instance in its own `<Hydro><MUSIC>` block, or `MpiMusic::set_store_evolution` |
| `<Hydro><MUSIC><restart_from>` | none | only in an instance's own block: the id (`<name>`) of another MUSIC instance; every evolution starts from that instance's latest checkpoint when it fits, and from tau0 otherwise, with the same output either way (the jet leg of the two-stage chain from the background's state; the log says `restart: started from the checkpoint at step …` or why not, `MUSIC_RESTART_VERIFY=1` checks it). Or `MpiMusic::set_restart_from(other)`, `get_restart_status()` |
| `<Hydro><MUSIC><checkpoint_tau>` | 0 | only in an instance's own block: keep a checkpoint of the state at the last freeze-out step by this tau [fm/c] in every evolution. 0: off, unless another instance names this one in `<restart_from>`; then just before `<Liquefier><CausalLiquefier><tau_delay>`, before which no droplet deposits: half the largest hydro step (`Delta_Tau`, `late_delta_tau`, the liquefier's `dtau`) and 1e-3 fm/c before it, so that a droplet made at the hard vertex (tau = 0) is not yet in reach of the jet source and the jet leg takes the checkpoint. A time between there and `tau_delay` is moved there (the log says so); other times are taken as given. Or `MpiMusic::set_checkpoint_tau(tau)` |

Two more MUSIC settings, which CPU MUSIC honours as well, can now also be set per instance: `<dump_hydro_only>` (1: keep MUSIC's native store, hand no copy to the framework's `bulk_info`) and `<skip_surface>` (1: do not export the freeze-out surface). As `freeze_out_surface` and `store_evolution`, the first `<Hydro><MUSIC>` block sets them for every instance and the same tag in an instance's own block overrides it (setters `set_dump_hydro_only`, `set_skip_surface`).

## Two-stage hydro with `runJetscape`

The jet leg (MUSIC_2, with the liquefier's droplets) starts from the background's (MUSIC_1) state and keeps no evolution of its own:

```xml
<Hydro>
  <AddLiquefier> true </AddLiquefier>
  <MUSIC>
    <name>MUSIC_2</name>
    <restart_from>MUSIC_1</restart_from>   <!-- MUSIC_1 keeps the checkpoint at tau_delay -->
    <dump_hydro_only>1</dump_hydro_only>   <!-- no framework copy of the jet leg -->
    <store_evolution>0</store_evolution>   <!-- and no evolution frames -->
    <skip_surface>0</skip_surface>         <!-- export its surface (iSS with hydro_id MUSIC_2) -->
  </MUSIC>
</Hydro>
```

Build the freeze-out surface on both legs or on neither: a checkpoint without the surface cannot serve a jet leg that builds one, and then every event runs from tau0 (with a warning). On the GB10 (0–10 % Au+Au, iSS + Colorless on MUSIC_2) this gives the same final-state hadrons as without the settings, at 6.3 instead of 15.9 GB peak memory and ~2 s less per event.

The js-contrib production drivers set these through the setters after `Init()`, which override the XML; see js-contrib's `docs/README_2stage.md`, `prod_AuAu_0_10_jet/README.md` and `docs/Plans/PLAN_music2_restart.md`.

## Example: one two-stage Au+Au event from XML

[`config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml`](../config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml)
runs one 0–10 % Au+Au event at 200 GeV with a jet, from `runJetscape` alone, without Python.
Its physics settings are those of js-contrib's production file
`AuAu_MCGlauber_MUSIC_0_10_jet.xml`. Each setting is explained next to it in the file.

```
IS (3dMCGlauber) -> Hard (PythiaGun) -> NullPreDynamics -> MUSIC_1 (background)
   -> Liquefier + Eloss (Matter + LBT)       droplets from the jet's energy loss
   -> MUSIC_2 (same strings + the droplets)  the jet leg: medium + wake
   -> iSS on MUSIC_2's freeze-out surface    bulk hadrons, status 11
   -> Colorless hadronization                jet fragments, status 0
```

Build with MUSIC4GPU, iSS and 3dMCGlauber, then run it from the build directory or in a
working directory of its own ([README_Launcher.md](README_Launcher.md)):

```bash
cmake .. -DUSE_MUSIC=ON -DUSE_CUDA=ON -DUSE_ISS=ON -DUSE_3DGlauber=ON   # -DUSE_METAL=ON on a Mac
make -j$(nproc)
./runJetscape ../config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml
# or, from the source directory:
examples/run_in_workdir.sh -d build -w RUN_DIR -u config/jetscape_user_AuAu_2stage_MUSIC2_hadrons.xml
```

The output, `AuAu_0_10_2stage_MUSIC2_hadrons.dat` (ASCII writer), holds the jet fragments
(status 0) and then iSS's hadrons of the jet leg (status 11).
On the GB10 (CUDA, `OMP_NUM_THREADS=10`, seed 1) the event takes 22 s and 9.5 GB peak memory, and
the file holds 32 jet hadrons and 7065 bulk hadrons (one iSS sample).

### What the file shows

| New in 3.0 | Where in the XML |
|---|---|
| Two MUSIC instances in one task list: the first `<Hydro><MUSIC>` block sets every instance, a block with `<name>` overrides it for that instance | `MUSIC_1`, `MUSIC_2` |
| Per instance: no freeze-out surface on the background, the surface only on the jet leg, no framework copy and no stored evolution of the jet leg | `<freeze_out_surface>`, `<skip_surface>`, `<dump_hydro_only>`, `<store_evolution>` |
| MUSIC4GPU's larger time step after 3 fm/c (0.02 to 0.05 fm/c), with the freeze-out checks and stored frames kept at their spacing | `<late_delta_tau>`, `<late_tau_switch>` |
| MUSIC4GPU's conservation fixes, on by default | `<symmetric_flux_fallback>`, `<quest_revert_as_cpu>` (main XML) |
| The liquefier's droplets deposited exactly on MUSIC's grid, diffused in the Bjorken frame | `<CausalLiquefier>` defaults (`normalize_on_hydro_grid`, `frame`) |
| Matter + LBT with recoils; energy loss starts at MUSIC's tau0 | `<Eloss>`, `<Matter><tStart_at_point>` |
| iSS samples a chosen hydro instance | `<SoftParticlization><hydro_id>MUSIC_2` |
| iSS decays the SMASH hadron list itself, for output that is not for SMASH | `<Perform_resonance_decays>1`, `<decay_smash_list>1` |
| 0–10 % centrality from the impact parameter (3dMCGlauber ignores `cenMin`/`cenMax` inside X-SCAPE) | `<b_max>4.7</b_max>` |

Settings to try on top of it:

- **Restart the jet leg from the background**: add `<restart_from>MUSIC_1</restart_from>`
  to the `MUSIC_2` block and set `<freeze_out_surface>1</freeze_out_surface>` in the first block
  (the part of the surface before the checkpoint comes from MUSIC_1; keep `<skip_surface>1`).
  The jet leg then skips the steps before the liquefier's `tau_delay`. In the event above it
  restarted at tau = 0.86 fm/c and wrote a bit-identical file, in 21 instead of 22 s; the
  saving grows with the number of jet legs per background, which needs Python today (below).
  See [Two-stage hydro with `runJetscape`](#two-stage-hydro-with-runjetscape) above.
- **Active region** (GPU only, opt-in): `<active_region_eps>1e-5</active_region_eps>` in the first
  block evolves only the blocks of cells near matter. Not bit-identical at the vacuum level.
- **Parallel showers**: `<Eloss><nThreads>-1</nThreads>`. Changes the events (one random
  stream per shower), so keep one setting per campaign.

Note that the two parts of the event are at different decay levels: Pythia decays π⁰, η and
Σ⁰ in the jet fragments, while iSS leaves them undecayed in the bulk (what SMASH calls stable).
Treat both alike before comparing with data.

### What needs Python today

The XML file gives one leg of one event. The js-contrib productions
([Au+Au jets in two stages](https://jhputschke.github.io/js-contrib/run/auau-jet/)) run the same
chain from the same two XML files, but step through the events from Python
(`JetScapePerEvent`) so that they can act between them:

| To get | XML today | Python drivers (js-contrib PyJetscape) |
|---|---|---|
| Hadrons of both legs of one event, for jet leg minus background | One `<SoftParticlization>` with one `<hydro_id>`: one leg per run | `hadronize.py` samples each leg's surface (`SampleCells`) |
| Several jets per background (`nReuseHydro` > 1) | The event loop switches iSS off on reused events together with MUSIC_1, although MUSIC_2 runs again: keep reuse off | `run_prod_jet.py --reuse N`, pT̂ windows sharing one background |
| Oversampling | `number_of_repeated_sampling` is fixed for the run, and the ASCII and HepMC writers keep no boundary between samples or legs | Samples per leg and per event (`--oversample-bg`); HDF5 hadron files keep every sample apart |
| Explicit seeds per event and leg (common seeds for both legs, redoing one event) | `<Random><seed>` and `<reseed_per_event>` only | `SetNextRandomSeed` / `GetLastRandomSeed` |
| Hadrons later, from stored surfaces and partons (stage 2 on CPUs) | No C++ reader: the hadrons are made in the job that ran the hydro | `hadronize.py`, `rehadronize_jets.py` |
| Python modules, the FNO hydro, hydro cells from numpy | `runJetscape` cannot load them | `@js.module`, PyJetscape |
| Work between events (diagnostics, restart status) | `runJetscape` calls `Exec()` once | per-event loop |

### XML improvements planned for future releases

Work is under way to close the gaps in the first four rows from the XML, so that a plain
`runJetscape` job gives paired hadrons of both legs. These are not part of 3.0:

1. **iSS follows the hydro it samples** when backgrounds are reused: with `<hydro_id>MUSIC_2`
   the jet leg is sampled on every event, also with `nReuseHydro` > 1.
2. **A second soft-particlization pass** for the other leg, so one job samples the background
   and the jet leg.
3. **Seed tags** for explicit seeds per event and leg.
4. **Oversampling per leg**, with the leg and the sample recorded in the output.
5. **A C++ reader for stored surfaces and partons**, for stage 2 without Python (if wanted).

Until then, use the Python drivers for productions, and this XML file to run and inspect one
event of the chain.
