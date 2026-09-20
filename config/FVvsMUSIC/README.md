# FV solver vs. MUSIC on one shared MC-Glauber IC

Runs a single Au+Au event through **MUSIC** and through the pure-Python finite-volume Milne
hydro in `FNO4d/loc_libs/fast_data/fv.py`, from the **same** initial condition, with **no jet
energy deposition**, and diffs the evolutions.

The IC comes from `fast_data`'s own tilted Bozek–Wyskieł MC-Glauber (`glauber.py`) — *not*
X-SCAPE's 3dMCGlauber or Trento — because the point is to test the two hydro solvers against
each other, not two different initial states.

Full design and the reasoning behind every setting: **`FNO4d/PLAN_fv_vs_music.md`**.
Measured results: **`FNO4d/RESULTS_fv_vs_music.md`** — the two codes agree to ~0.3 % in energy
density over a full central Au+Au lifetime in the ideal legs, and ~1–3 % with Israel–Stewart
shear.

## Files here

| file | what |
|---|---|
| `music_input_fv` | pristine MUSIC parameter file. `run_music_leg.py` copies it to `music_input_fv_<leg>` in the working directory, because `MpiMusic::InitializeHydro` rewrites `EOS_to_use` and `Include_Bulk_Visc` in place (`MusicWrapper.cc:124,243`). |
| `fv_vs_music_ideal_conformal.xml` | leg 1 — ideal, MUSIC `EOS_to_use 0` (ideal gas, dof 42.25) |
| `fv_vs_music_ideal_eos91.xml` | leg 2 — ideal, EOS 91 (hotQCD/SMASH lattice) |
| `fv_vs_music_is_eos91.xml` | leg 3 — Israel–Stewart shear, η/s = 0.08, no bulk, no second-order terms |
| `run_music_leg.py` | Python `InitialState` subclass + driver; exports the native store to `.npz` |

The FNO4d side lives in `FNO4d/workflow_fastdata/`: `fvmusic_eos_check.py`,
`fvmusic_make_ic.py`, `fvmusic_run_fv.py`, `fvmusic_compare.py`.

## How the IC gets in, without any new C++

`run_music_leg.py` subclasses `InitialState` through PyJetscape's `PyInitialState` trampoline
and builds the task list in Python:

```
FvMusicInitialState (Python)  ->  NullPreDynamics  ->  MUSIC
```

`NullPreDynamics` copies `ini->GetEntropyDensityDistribution()` straight into `e_` as **energy
density** — the name is a misnomer — with `u = (1,0,0,0)`, `π = 0`, `Π = 0`
(`NullPreDynamics.cc:36-56`). MUSIC's `Initial_profile 42` reader then indexes that vector as
`idx = (ny*neta)*ix + neta*iy + ieta` (`init.cpp:1239`), which is plain C order on
`(nx, ny, neta)` — exactly `e0[event].ravel()` from the HDF5.

Only `Exec` is overridden. `Init` is bound to `&JetScapeModuleBase::Init`
(`bind_framework.cc:162`) and dispatches virtually, so a Python `Init` calling `super().Init()`
would re-enter itself; leaving it to C++ also lets `InitialState::Init()` read the
`<IS><grid_*>` tags as usual. `InitialState::ExecuteTask()` is an empty stub, so nothing is
lost by not chaining.

## Running it

```bash
# 0. de-risk: how far apart are the two EoS implementations?  (sets the agreement floor)
conda activate fno_env
cd ~/FNO4d/workflow_fastdata
X=~/JetScape/X-SCAPE
python fvmusic_eos_check.py \
    --music-check $X/build_gpu/check_EoS_91_PST.dat \
    --eos-binary  $X/build_gpu/EOS/hotQCD/hrg_hotqcd_eos_SMASH_binary.dat \
    --plot out/fvmusic_eos_check.png
```

`check_EoS_91_PST.dat` is MUSIC's own dump of the EoS it actually uses:

```bash
cd $X/build_gpu
sed -e 's/^mode 2/mode 71/' -e 's/^EOS_to_use .*/EOS_to_use 91/' music_input > /tmp/mi_eos
./external_packages/music4gpu/src/MUSIChydro /tmp/mi_eos
```

```bash
# 1. the shared IC, plus an asymmetric probe for the gate
python fvmusic_make_ic.py --mode glauber --out out/ic_fvmusic_AuAu_b0.h5 \
    --eos-binary $X/build_gpu/EOS/hotQCD/hrg_hotqcd_eos_SMASH_binary.dat --overwrite
python fvmusic_make_ic.py --mode probe   --out out/ic_fvmusic_probe.h5 \
    --eos-binary $X/build_gpu/EOS/hotQCD/hrg_hotqcd_eos_SMASH_binary.dat --overwrite

# 2. THE GATE.  Nothing downstream means anything until this passes.
cd $X/build_gpu
MUSIC_FORCE_CPU=1 OMP_WAIT_POLICY=passive OMP_NUM_THREADS=8 \
python ../config/FVvsMUSIC/run_music_leg.py --leg ideal_conformal --probe \
    --ic ~/FNO4d/workflow_fastdata/out/ic_fvmusic_probe.h5 \
    --out out_fvcmp/music_probe.npz

# 3+4. one leg, both codes, then the diff
LEG=ideal_conformal
MUSIC_FORCE_CPU=1 OMP_WAIT_POLICY=passive OMP_NUM_THREADS=8 \
python ../config/FVvsMUSIC/run_music_leg.py --leg $LEG \
    --ic ~/FNO4d/workflow_fastdata/out/ic_fvmusic_AuAu_b0.h5 \
    --out out_fvcmp/music_$LEG.npz

cd ~/FNO4d/workflow_fastdata
python fvmusic_run_fv.py --leg $LEG --ic out/ic_fvmusic_AuAu_b0.h5 \
    --music-npz $X/build_gpu/out_fvcmp/music_$LEG.npz --out out/fv_$LEG.npz
python fvmusic_compare.py --fv out/fv_$LEG.npz \
    --music $X/build_gpu/out_fvcmp/music_$LEG.npz --plot out/fvmusic_$LEG.png
```

Repeat step 3+4 with `LEG=ideal_eos91` and `LEG=is_eos91`.

`MUSIC_FORCE_CPU=1` and `<beastMode>0` are deliberate: this is a numerics comparison, so it
must measure MUSIC's algorithm and not the fp32 GPU path. (For the ideal legs the GPU path is
off anyway — `advance.cpp:333` returns early when `viscosity_flag != 1`.)
`OMP_WAIT_POLICY=passive` with ~8 threads avoids idle OpenMP threads spinning.

## The gate, and what it catches

`--probe` injects one Gaussian blob at deliberately distinct, off-centre `x`, `y` and `η`,
runs a couple of τ steps, and asserts that MUSIC's frame 0 reproduces the injected array
cell for cell. Because no two axes are interchangeable, a transposed index, a flipped axis,
a wrong `s_factor` or a missing `hbarc` all move the arg-max cell or blow the tolerance.

Measured result:

```
arg-max cell : IC (40, 28, 19)  MUSIC (40, 28, 19)   OK
peak         : IC 30.000002  MUSIC 30.000002
max rel diff : all cells 7.535e-08, fluid cells 7.367e-08
eta axis     : fast_data [-5.00000, +5.00000]  MUSIC [-5.15625, +4.84375]  (offset -0.15625)
=> PASS
```

The run log should also show, and the driver prints enough to check:

```
Using Initial_profile = 42. Overwriting lattice dimensions:
neta = 33, nx = 65, ny = 65
Running ideal hydrodynamic simulations ...
Setting the initial viscous tensor to zero.      <- init.cpp:349-372, see below
Number of source terms: 0
Total E sources = 0 GeV.                         <- no jet deposition
```

## Settings that are load-bearing

Changing any of these silently invalidates the comparison.

- **`<MUSIC><EOS>` does NOT choose the EoS that MUSIC evolves with.** MUSIC builds its EoS
  in the constructor initialiser list, `eos(DATA.whichEOS)` (`music.cpp:24-26`), from
  `EOS_to_use` in the MUSIC *input file*. The wrapper's later
  `set_parameter("EOS", v)` only assigns `DATA.whichEOS`
  (`read_in_parameters.cpp:1012`) and never rebuilds that object; the wrapper then rewrites
  `EOS_to_use` in the file (`MusicWrapper.cc:124`) for iSS's benefit, so **the file ends up
  looking right while the run used the old EoS**. This is not specific to this study — it
  applies to any X-SCAPE MUSIC run whose XML `<EOS>` differs from the input file.
  `run_music_leg.py` therefore writes `EOS_to_use` into its per-leg copy before MUSIC is
  constructed, and the log line to check is `initialze EOS ideal gas` vs `reading EOS hotQCD`.
  *(Found the hard way: the first conformal leg silently ran hotQCD.)*
- **`<temperature_dependent_bulk_viscosity>` can only turn bulk viscosity ON.**
  `MusicWrapper.cc:236` is `if (flag_bulkvis != 0) set_parameter(...)`, so a `0` in the XML
  leaves whatever `Include_Bulk_Visc_Yes_1_No_0` the input file already had. The stock
  `build_gpu/music_input` ships `3`, so the value **must** be set to `0` in
  `music_input_fv` — which it is. *(Also found the hard way.)*
- **`<Preequilibrium><taus>` is MUSIC's τ₀**, not `<MUSIC><Initial_time_tau_0>`
  (`MusicWrapper.cc:391` → `PreequilibriumDynamics.cc:58`). `run_music_leg.py` refuses to run
  if it disagrees with the IC file's `tau0`.
- **`<T_dependent_Shear_to_S_ratio>` must be 0** on the ideal legs. Any positive value
  re-enables `Viscosity_Flag_Yes_1_No_0` at `MusicWrapper.cc:168` *regardless* of
  `shear_viscosity_eta_over_s`.
- **`<temperature_dependent_bulk_viscosity>` must be 0.** Nonzero sets
  `Include_Bulk_Visc = 1` at `MusicWrapper.cc:238` independently of the shear setting, and
  `fast_data`'s bulk sector is untested.
- **`viscosity_flag = 0` is what makes the ideal legs clean.** `NullPreDynamics` hardcodes
  `P_ = e/3`, and `initial_with_jetscape` turns that into
  `piBulk = P_in/hbarc − p_EoS(e)` — nonzero for a non-conformal EoS. `init.cpp:349-372`
  zeroes `piBulk_` and all 14 `Wmunu_` at τ₀ whenever `viscosity_flag == 0`, which is why the
  log line *"Setting the initial viscous tensor to zero"* matters. On leg 3 the stale value
  survives init and shows up in the τ₀ frame; it is dynamically inert (every consumer is
  gated on `turn_on_bulk`) but should be checked rather than assumed.
- **`<evolutionInMemory>0`.** With 1, `MusicWrapper.cc:505-508` overrides MUSIC's `Delta_Tau`
  with the pre-eq `dtau` and prepends pre-eq frames to the τ axis.
- **`<setReuseHydro>false`.** The `jetscape_main.xml` default reuses hydro across 10 events.
- **`beastMode` 0, never 2.** 2 grows `delta_tau` mid-run (`evolve.cpp:98-106`), so the
  recorded `dtau` stops describing the τ axis.
- **`output_evolution_every_N_x/y/eta` stay 1.** `HydroinfoMUSIC.cpp:365,377` and
  `grid_info.cpp:384` compute the skipped extent differently and agree only at 1.
- **Square transverse grid.** `Initial_profile 42` recovers `nx = sqrt(size/neta)` and sets
  `ny = nx` (`init.cpp:135-137`); separately, `MusicWrapper.cc:629-633` sets
  `bulk_info.ny = get_nx()` and takes `y_min`/`dy` from the x accessors, which is harmless
  only while `nx == ny`.

### Two harmless log lines

- `DeltaX = 0.317383 fm` in the startup banner is `X_grid_size_in_fm/(nx-1)` from the file
  (`read_in_parameters.cpp:607`), printed before `music.cpp:245` overwrites it with the IC's
  `dx = 0.3125`. The grid printed after `Using Initial_profile = 42` is the real one. This is
  also why `reset_dtau_use_CFL_condition` must stay **0** — `check_parameters()` runs before
  the overwrite, so a CFL-derived `Delta_Tau` would use the stale spacing.
- `non-zero eta/s = 0.08 is set with Include_Shear_Visc = 0` on the ideal legs: MUSIC checks
  the *file's* `Shear_to_S_ratio` during construction, before the wrapper applies the XML.
  The authoritative line is `Running ideal hydrodynamic simulations ...`.
