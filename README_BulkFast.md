# FastRootBulkWriter — fast hydro-evolution ROOT dumps

`FastRootBulkWriter` writes the MUSIC hydro evolution (energy density + flow
velocities) to a ROOT file **directly from MUSIC's native in-memory store**,
bypassing the JETSCAPE framework history (`bulk_info.data`) that the regular
`RootBulkWriter` reads back through per-cell interpolation. It is meant for
**hydro-only training-data generation** (e.g. the FNO pipeline).

See [BulkWriterImprovements.md](BulkWriterImprovements.md) for the design rationale.

---

## 1. When to use it

| | `RootBulkWriter` (legacy) | `FastRootBulkWriter` (new) |
|---|---|---|
| Reads from | framework AoS `bulk_info.data` | MUSIC native store (no AoS) |
| Interpolation | always | only in `grid` mode with a non-matching grid |
| Jet quenching in same run | yes | **no** (hydro-only) |
| Output grid | user grid (interpolated) | MUSIC grid, or user grid |

Use `FastRootBulkWriter` when you only want the hydro evolution ROOT file and are
**not** running energy loss / hadronization in the same job.

> ⚠️ **Hydro-only only.** `dump_hydro_only=1` leaves the framework medium empty, so
> any `Eloss` (`Matter`/`Martini`/`LBT`) / hadronization / afterburner module that
> queries the medium will throw (`EvolutionHistory is empty`). Remove `<Eloss>`
> (and any hadronization/afterburner) from the config. Derive a hydro-only config by
> deleting just those blocks, leaving the rest intact.

---

## 2. Quick start

### 2a. Tell MUSIC to keep its native store

Inside `<Hydro><MUSIC>` (alongside the existing evolution flags):

```xml
<output_evolution_to_memory>1</output_evolution_to_memory>  <!-- required -->
<dump_hydro_only>1</dump_hydro_only>                         <!-- skip the framework copy -->
<output_evolution_every_N_timesteps>1</output_evolution_every_N_timesteps>
```

`dump_hydro_only=1` makes `EvolveHydro` keep MUSIC's evolution in memory and **skip**
building `bulk_info.data`. It requires `output_evolution_to_memory=1`.

### 2b. Add the writer — `native` mode (fastest)

Top level, like the old `<RootBulkWriter>`:

```xml
<FastRootBulkWriter>
  <name> FastRootBulkWriter </name>
  <out_file_name> OO_test_fast.root </out_file_name>
  <grid_mode> native </grid_mode>     <!-- write MUSIC's own grid -->
  <tau_stride> 1 </tau_stride>        <!-- write every Nth stored timestep -->
</FastRootBulkWriter>
```

`native` mode does a pure grid-walk: no interpolation, no AoS. The spatial grid is
MUSIC's; control resolution at the source with `output_evolution_every_N_x/y/eta`, and
thin τ further here with `tau_stride`.

### 2c. Or `grid` mode (custom grid)

```xml
<FastRootBulkWriter>
  <name> FastRootBulkWriter </name>
  <out_file_name> OO_test_fastgrid.root </out_file_name>
  <grid_mode> grid </grid_mode>
  <tau_min> 0.524 </tau_min>  <dtau> 0.1 </dtau>  <ntau> 0 </ntau>  <!-- ntau=0 => to end -->
  <x_min> -10 </x_min>    <dx> 0.3125 </dx>
  <y_min> -10 </y_min>    <dy> 0.3125 </dy>
  <eta_min> -5 </eta_min> <deta> 0.3125 </deta>
</FastRootBulkWriter>
```

`grid` mode honours your grid; if it matches MUSIC's it grid-walks, otherwise it
interpolates (reusing `EvolutionHistory::get`). Any tag left at `0` defaults to MUSIC's
value.

### 2d. Run

Ready-made example configs live in [config/BulkFastTest/](config/BulkFastTest/). Run from
`build_gpu/` (so MUSIC's input/EOS paths resolve):

```bash
cd build_gpu
./runJetscape ../config/BulkFastTest/OO_one_event_fast.xml       # native -> OO_test_fast.root
./runJetscape ../config/BulkFastTest/OO_one_event_fastgrid.xml   # grid   -> OO_test_fastgrid.root
./runJetscape ../config/BulkFastTest/OO_one_event.xml            # legacy -> OO_test.root
```

You should see: `MUSIC dump_hydro_only: retaining native store (… cells), skipping
framework copy.`

---

## 3. XML reference (`<FastRootBulkWriter>`)

| tag | mode | default | meaning |
|---|---|---|---|
| `out_file_name` | both | `hydro_evo_fast.root` | output ROOT file |
| `grid_mode` | both | `native` | `native` or `grid` |
| `tau_stride` | native | `1` | write every Nth stored timestep |
| `include_preequilibrium` | both | `0` | reserved (not implemented) |
| `tau_min`,`dtau`,`ntau` | grid | MUSIC / `0` | τ origin, step, count (`ntau=0` ⇒ to end) |
| `x_min`,`dx` / `y_min`,`dy` / `eta_min`,`deta` | grid | MUSIC | output spatial grid |

Plus `<Hydro><MUSIC><dump_hydro_only>1` and `<output_evolution_to_memory>1`.

---

## 4. Output format & reading it

Branches per event:
- `user_res` — `vector<float>`, flattened `[tau][x][y][eta][feature]`, feature order
  `(energy_density, vx, vy, vz)`.
- `ntau_freezeout` — `int`, number of τ steps written this event.

Metadata (`TParameter`/`TNamed`, written once): `nFeatures, nx, ny, neta, x_min, dx,
y_min, dy, eta_min, deta, tau_min, dtau, ntau, tau_stride, grid_mode`. `tau_min`/`dtau`
are the **effective** values (in native mode `dtau` already includes `tau_stride`).

```python
import uproot, numpy as np
f = uproot.open("OO_test_fast.root")
nx, ny, neta, nF = (f[k].member("fVal") for k in ("nx","ny","neta","nFeatures"))
t = f["t"]
res  = t["user_res"].array(library="np")          # ragged: one vector per event
ntau = t["ntau_freezeout"].array(library="np")
arr  = np.asarray(res[0]).reshape(ntau[0], nx, ny, neta, nF)
energy, vx, vy, vz = arr[...,0], arr[...,1], arr[...,2], arr[...,3]
```

Same layout as the legacy `RootBulkWriter` (with `ntau=0`), so existing readers work
with at most a filename change.

---

## 5. What changed in the code (technical)

| file | change |
|---|---|
| `src/root/FastRootBulkWriter.{h,cc}` | **new** module; native grid-walk + grid-mode interpolation |
| `src/hydro/MusicWrapper.h` | `MpiMusic`: `get_number_of_fluid_cells()`, `clear_hydro_info_from_memory()`, `get_native_fluid_cell(idx, FluidCellInfo&)`, `set/get_dump_hydro_only`, `dump_hydro_only_` |
| `src/hydro/MusicWrapper.cc` | read `<dump_hydro_only>`; `EvolveHydro` keeps the native store + skips `PassHydroEvolutionHistoryToFramework` when set; defensive clear each event |
| `src/framework/JetScape.cc` | include + `DetermineTaskListFromXML` registration of `FastRootBulkWriter` |
| `src/root/RootBulkWriter.cc` | perf fix: `auto bInfo` → `const auto &bInfo` (was deep-copying the whole history every event) |
| `src/CMakeLists.txt` | compile `root/FastRootBulkWriter.cc` under `USE_ROOT` |
| `config/jetscape_main.xml` | **required** default entries for `<dump_hydro_only>` and `<FastRootBulkWriter>` (X-SCAPE rejects user-XML tags with no default — `recurseToSearch: tag unrecognized`) |

**Data path.** MUSIC stores each output timestep in `HydroinfoMUSIC::lattice_ideal`
(tau-major, then x, y, eta). `native` mode walks that store via
`MpiMusic::get_native_fluid_cell` (MUSIC's own u^μ→v + hbarc conversion) — no
interpolation, no `bulk_info.data`. `grid` mode builds a **transient**
`EvolutionHistory` from the native store once per event and calls `get()`; freed after
the event. The writer clears MUSIC's native store at the end of each `Exec`; `MpiMusic`
also clears it defensively at the start of each event.

> Build note: after editing `src/CMakeLists.txt`, **re-configure CMake**
> (`cmake build_gpu`) before `cmake --build`, or the new source isn't compiled and the
> module's self-registration never runs (writer silently absent → no file). Confirm with
> `nm -gU build_gpu/src/lib/libJetScape.dylib | grep FastRootBulk`.

---

## 6. Verification (O+O, 1 event, build_gpu)

Reproduce with the script [config/BulkFastTest/validate.py](config/BulkFastTest/validate.py)
**under a numpy+uproot Python** (conda `fno_env` — the system `python3` here has no numpy):

```bash
cd build_gpu
# (run the three configs from §2d first)
conda activate fno_env
python ../config/BulkFastTest/validate.py
```

Observed results:

- **Build + run:** clean build; all three configs run `EXIT=0`; the log confirms
  `MUSIC dump_hydro_only: retaining native store (… cells), skipping framework copy`, so
  `bulk_info.data` is never built.
- **`grid` mode vs legacy `RootBulkWriter`: bit-for-bit identical.** Same shape (both
  `ntau_freezeout=27`, `user_res` of 15,057,900 floats on the 65×65×33 grid) and
  `np.array_equal == True` (`max |abs diff| = 0`). Files ≈ 16.83 MB each. This confirms the
  fast path reproduces the legacy interpolation exactly when the output grid is the same.
- **`native` mode:** writes MUSIC's full grid (here 100×100×60, ntau=134 at MUSIC's
  `dtau≈0.02`), file ≈ 93 MB; no NaN/inf; peak energy density 7.48 vs legacy 7.36 — native
  is slightly higher because it keeps MUSIC's full spatial+τ resolution (the legacy
  coarser grid/τ slightly under-resolves the peak), which is expected, not an error.

> Run `validate.py` for exact current numbers — it prints shapes, file sizes, max/mean
> diff, and the native-vs-legacy peak.

---

## 7. Performance & caveats

- **native mode** removes the framework AoS build (`PassHydro…`), the per-event
  whole-history deep copy, and all interpolation. Biggest win on 3D grids, where the AoS
  dominates memory.
- **native mode writes MUSIC's full resolution** — it can be large (≈93 MB for one O+O
  event at 100×100×60×134; bigger on finer grids). Shrink it with `<tau_stride>`, with
  MUSIC's `output_evolution_every_N_*`, or use `grid` mode for a fixed coarse grid.
- **grid mode** removes the persistent AoS + the deep copy, but still builds one
  transient working copy and interpolates when the grid differs from MUSIC's.
- The fast dump starts at MUSIC's `hydroTau0` (pre-equilibrium steps are **not**
  prepended). A τ-origin offset vs a legacy file with a different `tau_min` is expected,
  not a bug.
- Hydro-only only: remove `<Eloss>`/hadronization/afterburner (see §1).
