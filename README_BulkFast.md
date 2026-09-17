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
<skip_surface>1</skip_surface>                               <!-- skip the freeze-out surface export -->
<output_evolution_every_N_timesteps>1</output_evolution_every_N_timesteps>
```

`dump_hydro_only=1` makes `EvolveHydro` keep MUSIC's evolution in memory and **skip**
building `bulk_info.data`. It requires `output_evolution_to_memory=1`.

`skip_surface=1` (optional) additionally skips exporting MUSIC's freeze-out surface to
the framework / `surface*.dat`. **MUSIC still finds the surface internally** — it is the
hydro stop condition, so it cannot be turned off without breaking termination — this only
drops the unused hand-off, which nothing consumes in a hydro-only dump. It is
output-neutral (the evolution ROOT file is byte-identical with it on or off).

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

In `native` mode on a 3D grid you may also see
`Error in <TBufferFile::WriteByteCount>: bytecount too large (more than 1073741822)`.
The data is still written; see §9.

### 2e. From Python (PyJetscape)

With js-contrib's PyJetscape built in-tree against this ROOT-enabled build, the writer is
the Python type `jetscape.FastRootBulkWriter`. Python can add it to a pipeline with
`create_module("FastRootBulkWriter")` or find it in `JetScape.GetTaskList()`. The file is
written and closed by `JetScape.Finish()`. `MpiMusic.get_native_evolution_numpy()` copies
the same native store straight into numpy, and `jetscape.fast_root_bulk.read_fast_root_bulk()`
reads the ROOT file back (§4). See the "C++ `FastRootBulkWriter` from Python"
section of `external_packages/js-contrib/contribs/PyJetscape/README.md` and
`example/python_fast_bulk_root_writer.py` there.

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

Plus, under `<Hydro><MUSIC>`: `<dump_hydro_only>1`, `<output_evolution_to_memory>1`, and
optionally `<skip_surface>1` (skip the unused freeze-out-surface export; output-neutral).

---

## 4. Output format & reading it

Tree `t`, one entry per event:
- `user_res` — `vector<float>`, flattened `[tau][x][y][eta][feature]`, feature order
  `(energy_density, vx, vy, vz)`.
- `ntau_freezeout` — `int`, number of τ steps written this event.
- `tau_freezeout` — `float`, τ one stored MUSIC step past the last stored step
  (`tau_min_MUSIC + N × dtau_MUSIC` for `N` stored steps, independent of `tau_stride`).

Metadata (`TParameter`, plus the `TNamed` `grid_mode`), written once at the first event:

| keys | meaning |
|---|---|
| `nFeatures`, `nx`, `ny`, `neta` | layout of `user_res` (`nFeatures` is always 4) |
| `tau_min`, `dtau` | **effective** τ axis of `user_res` (in `native` mode `dtau` already includes `tau_stride`) |
| `x_min`, `dx`, `y_min`, `dy`, `eta_min`, `deta`, `ntau` | the `<FastRootBulkWriter>` grid. In `grid` mode this is the output grid (tags left at 0 are filled with MUSIC's value, except `ntau`). In `native` mode these are just the XML values (0 by default) and **do not describe the data** |
| `nX_MUSIC`, `X_min_MUSIC`, `dX_MUSIC`, `nY_MUSIC`, `Y_min_MUSIC`, `dY_MUSIC`, `neta_MUSIC`, `eta_min_MUSIC`, `deta_MUSIC`, `tau_min_MUSIC`, `dtau_MUSIC` | MUSIC's stored grid; the coordinates of `user_res` in `native` mode |
| `tau_stride` | `native`-mode τ thinning |
| `use_vec` | always 1 (`user_res` is a `vector<float>`) |
| `grid_mode` (`TNamed`, value in the title) | `native` or `grid` |

```python
import uproot, numpy as np
f = uproot.open("OO_test_fast.root")
P = lambda k: f[k].member("fVal")
nx, ny, neta, nF = (P(k) for k in ("nx","ny","neta","nFeatures"))
t = f["t"]
res  = t["user_res"].array(library="np")          # ragged: one vector per event
ntau = t["ntau_freezeout"].array(library="np")
arr  = np.asarray(res[0]).reshape(ntau[0], nx, ny, neta, nF)
energy, vx, vy, vz = arr[...,0], arr[...,1], arr[...,2], arr[...,3]
x   = P("X_min_MUSIC") + P("dX_MUSIC") * np.arange(nx)   # native mode; grid mode: x_min, dx
tau = P("tau_min") + P("dtau") * np.arange(ntau[0])
```

A native-mode event is large (1.29 GB uncompressed for the O+O example), so load one
event at a time with `t["user_res"].array(entry_start=i, entry_stop=i+1, library="np")`.

From PyJetscape, `jetscape.fast_root_bulk.read_fast_root_bulk(path, entry_start,
entry_stop)` does all of the above: it returns one `(ntau, nx, ny, neta, 4)` array per
event and picks the right coordinate keys for the grid mode.

The branches and the shared metadata keys match the legacy `RootBulkWriter` with `ntau=0`
(which also writes `vector<float>`, `use_vec=1`), so existing readers work with at most a
filename change. The fast writer adds `tau_stride` and `grid_mode`.

Files from earlier versions of the fast writer have no `tau_freezeout` branch and no
`use_vec` or `*_MUSIC` keys, so native-mode coordinates can't be recovered from them;
`read_fast_root_bulk` returns `None` for those fields.

---

## 5. What changed in the code (technical)

| file | change |
|---|---|
| `src/root/FastRootBulkWriter.{h,cc}` | **new** module; native grid-walk + grid-mode interpolation; `tau_freezeout` branch and `use_vec` / `*_MUSIC` metadata matching the legacy writer; writes and closes the file in `FinishTask()` (called by `JetScape::Finish()`; the destructor calls it too, and repeated calls do nothing); read-only getters (`GetOutFileName()`, `GetNumberOfEventsWritten()`, `GetNtauWritten()`, …) for the Python bindings |
| `src/hydro/MusicWrapper.h` | `MpiMusic`: `get_number_of_fluid_cells()`, `clear_hydro_info_from_memory()` (both safe before `InitializeHydro()`), `get_native_fluid_cell(idx, FluidCellInfo&)`, `set/get_dump_hydro_only`, `set/get_skip_surface` |
| `src/hydro/MusicWrapper.cc` | read `<dump_hydro_only>` and `<skip_surface>`; `EvolveHydro` keeps the native store + skips `PassHydroEvolutionHistoryToFramework` when set; defensive clear each event; in dump mode sets `bulk_info.tau_min`/`dtau`/`ntau` from MUSIC's stored steps (see τ axis note below) |
| `src/framework/JetScape.cc` | include + `DetermineTaskListFromXML` registration of `FastRootBulkWriter` |
| `src/root/RootBulkWriter.cc` | perf fix: `auto bInfo` → `const auto &bInfo` (was deep-copying the whole history every event) |
| `src/CMakeLists.txt` | compile `root/FastRootBulkWriter.cc` under `USE_ROOT` |
| `config/jetscape_main.xml` | **required** default entries for `<dump_hydro_only>`, `<skip_surface>` and `<FastRootBulkWriter>` (X-SCAPE rejects user-XML tags with no default — `recurseToSearch: tag unrecognized`) |
| `external_packages/js-contrib/contribs/PyJetscape/` (separate repository) | `src/bind_root_bulk_writer.cc`: Python type `FastRootBulkWriter` (only with `USE_ROOT`; `jetscape.HAS_ROOT`). `src/bind_music.cc`: `get_native_evolution_numpy(tau_stride)` and the native-store accessors. `python/jetscape/fast_root_bulk.py`: reader. `example/python_fast_bulk_root_writer.py`: example |

**Data path.** MUSIC stores each output timestep in `HydroinfoMUSIC::lattice_ideal`
(tau-major, then x, y, eta). `native` mode walks that store via
`MpiMusic::get_native_fluid_cell` (MUSIC's own u^μ→v + hbarc conversion) — no
interpolation, no `bulk_info.data`. `grid` mode builds a **transient**
`EvolutionHistory` from the native store once per event and calls `get()`; freed after
the event. The writer clears MUSIC's native store at the end of each `Exec`; `MpiMusic`
also clears it defensively at the start of each event. The tree is written and the file
closed in `FinishTask()`, so the output is complete after `JetScape::Finish()` without
relying on destructor timing (which matters from Python).

**τ axis.** The writer takes MUSIC's grid (including `tau_min`, `dtau`) from `bulk_info`.
With `<Preequilibrium><evolutionInMemory>1` (the main-XML default) and initial profile 42,
`bulk_info` normally describes the *combined* pre-equilibrium + hydro history, so its
`tau_min` is the pre-equilibrium start. The native store holds hydro steps only, so in
dump mode `MpiMusic` now overrides `tau_min`, `dtau` and `ntau` with MUSIC's own stored
values (`get_hydro_tau0()`, `get_hydro_dtau()`, `get_ntau()`).

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

Later checks:

- **Python path** (PyJetscape, `example/python_fast_bulk_root_writer.py`): on the O+O
  native config, the file from `runJetscape`, the file from the Python-driven run, and
  `MpiMusic.get_native_evolution_numpy()` are **bit-for-bit identical**
  (`np.array_equal`, `max |diff| = 0`). A manual Python pipeline (Trento + MUSIC 2D,
  `tau_stride=5`) gives the same numpy-vs-ROOT match.
- **τ axis fix** (§5): on a Trento + NullPreDynamics 2D config with pre-equilibrium kept in
  memory (hydro starts at 0.5 fm/c), `native` mode now records `tau_min = 0.5` (was 0) and
  `tau_freezeout = 30.52` (was 30.02). With a fixed seed, `grid` mode with
  `<tau_min>1.0</tau_min>` matches the native rows at the same τ (energy density
  `max |diff| = 1.4e-6`); the old axis would have sampled τ + 0.5 fm/c. The O+O native
  output (`evolutionInMemory = 0`) is unchanged, bit for bit.
- **Large events:** ROOT C++ and uproot read identical data from native O+O events over
  1 GB; see §9.

---

## 7. Performance

Measured on the O+O 1-event config above, in `build_gpu`, with `/usr/bin/time -l`
(sequential runs; peak RSS was identical across two passes, so these are stable):

| run | wall time | peak RSS | output file |
|---|---|---|---|
| legacy `RootBulkWriter` | 13.3 s | **24.6 GB** | 16.1 MB |
| fast — `grid` mode | 11.2 s | **15.5 GB** | 16.1 MB |
| fast — `native` mode | 12.6 s | **11.0 GB** | 93.2 MB |

**Memory is the real win.** All three runs are identical until MUSIC finishes the hydro
solve at **~7.3 GB** (the native store: 80.4 M cells here). Then:

- **legacy** → **25.2 GB**: `PassHydroEvolutionHistoryToFramework` builds the framework
  AoS (`bulk_info.data`) *on top of* the still-live native store, so both coexist at peak.
- **grid** → **15.9 GB**: builds a *transient* `EvolutionHistory` (no second framework
  copy), freed after the event.
- **native** → **7.5 GB**: no AoS at all — just the native store + the output float vector.

So `native` mode cuts peak RSS **~55%** (24.6 → 11.0 GB) and `grid` mode **~37%**. The
≈18 GB framework-AoS materialization is what's eliminated; this saving scales with grid
size, so it grows on larger / 3+1D grids.

**`skip_surface` is a small, output-neutral extra.** A/B on the native config (2 passes):
peak RSS **unchanged** (10.96 GB both — the freeze-out surface is tiny next to the 80 M-cell
evolution store), wall time ~12.8–13.5 s → ~11.1–11.3 s (≈1.5–2 s, noisy). Output is
byte-identical (`np.array_equal`, same file size). So it trims a bit of wasted
surface-export work but is not a memory lever; the framework-AoS bypass above is the win.

**Runtime improvement is modest (~5–16%) and not the point.** The MUSIC hydro solve
dominates total time and is identical in all three runs — the writer/copy is a small
slice. `native` is slightly slower than `grid` here only because it writes 93 MB vs 16 MB.

> Numbers are for one boost-invariant O+O event on one machine — representative, not
> universal. The memory figures reproduced exactly across two passes; the wall-times are
> noisier at n=1. Re-measure with `/usr/bin/time -l ./runJetscape …` for your case.

## 8. Caveats

- **native mode** removes the framework AoS build (`PassHydro…`), the per-event
  whole-history deep copy, and all interpolation. Biggest win on 3D grids, where the AoS
  dominates memory.
- **native mode writes MUSIC's full resolution** — it can be large (≈93 MB for one O+O
  event at 100×100×60×134; bigger on finer grids). That is *disk*, not RAM. Shrink it with
  `<tau_stride>`, with MUSIC's `output_evolution_every_N_*`, or use `grid` mode for a
  fixed coarse grid. Uncompressed, that event is 1.29 GB, which is over ROOT's 1 GB
  per-object limit; see §9.
- **grid mode** removes the persistent AoS + the deep copy, but still builds one
  transient working copy and interpolates when the grid differs from MUSIC's.
- The fast dump starts at MUSIC's `hydroTau0` (pre-equilibrium steps are **not**
  prepended). A τ-origin offset vs a legacy file with a different `tau_min` is expected,
  not a bug.
- Files written by earlier versions with `<Preequilibrium><evolutionInMemory>1` and
  initial profile 42 record the **pre-equilibrium** start as `tau_min` (and a
  correspondingly shifted `tau_freezeout`), e.g. 0 instead of 0.5 fm/c. Their `native` data
  is correct but labelled with the wrong τ; `grid` mode with an explicit `<tau_min>`
  sampled the wrong times. Fixed in `MpiMusic::EvolveHydro` (§5, τ axis).
- Hydro-only only: remove `<Eloss>`/hadronization/afterburner (see §1).

---

## 9. Per-event size limit (ROOT's 1 GB object limit)

**Status: not fixed.** Large `native`-mode events currently work, but they log ROOT errors
and are close to a harder limit.

### Symptom

Writing prints, once per oversized event:

```
Error in <TBufferFile::WriteByteCount>: bytecount too large (more than 1073741822)
```

Reading the file with ROOT C++ prints, once per oversized entry:

```
Error in <TBufferFile::CheckByteCount>: object of class vector<float> read too many bytes: 1286400006 instead of 212658182
Warning in <TBufferFile::CheckByteCount>: vector<float>::Streamer() not in sync with data on file …, fix Streamer()
```

### Cause

When ROOT writes an object into a basket, it stores the object's length in a 32-bit field
in front of it. The two highest bits of that field are reserved for flags, so only 30 bits
hold the value: the largest length it can store is 2³⁰ − 2 = **1,073,741,822 bytes**.

`FastRootBulkWriter` writes each event as one `std::vector<float>` (`user_res`), so the
object size is

```
bytes per event = 16 × (ntau × nx × ny × neta) + 6      (4 features × 4 bytes, + header)
```

For the O+O example (§6):

| | O+O native event |
|---|---|
| cells | 134 × 100 × 100 × 60 = 80.4M |
| floats (× 4 features) | 321.6M |
| bytes | **1,286,400,006** |

That is over the limit. ROOT logs the error and keeps only the low 30 bits of the length:
1,286,400,006 − 2³⁰ = **212,658,182**, which is the value the read-side error reports.

### What still works

The vector also stores its own element count, and readers use that count to read the data.
ROOT then only complains that the header length doesn't match. Verified on a 2-event
native O+O file (seed 42) with ROOT 6.36.06 and uproot 5.6.9:

| entry | ntau | bytes | stored length | ROOT C++ vs uproot |
|---|---|---|---|---|
| 0 | 134 | 1,286,400,006 | 212,658,182 | identical bits |
| 1 | 199 | 1,910,400,006 | 836,658,182 | identical bits |

- `TTree::GetEntry` and uproot return bit-identical data for both entries (compared with a
  hash over every float's bit pattern).
- A bad entry doesn't disturb the next one. The writer calls `SetAutoFlush(1)`, so each
  event is in its own basket.
- uproot readers (§4, `validate.py`, PyJetscape's `read_fast_root_bulk`) are unaffected and
  print no errors.

### Where it is a real problem

1. **The ~2 GB limit.** ROOT's buffer sizes are signed 32-bit, so a single entry can't go
   much past ~2.15 GB. Entry 1 above is already **1.91 GB**. What happens past that has not
   been tested; expect the write to fail outright rather than only log an error.
2. **Other tools.** Anything that skips an object using its stored length instead of reading
   it could misread these entries. Only `TTree::GetEntry` and uproot have been checked, not
   `hadd`, `TTree::CopyTree`, RDataFrame or other readers.
3. **Noisy logs.** Every oversized event adds errors that can hide real ones.

### When it triggers

| limit | max cells per event (ntau × nx × ny × neta) | 100×100×60 grid (600k cells, 9.6 MB per τ step) |
|---|---|---|
| 1 GB (errors) | 67,108,863 | ntau ≥ **112** |
| ~2 GB (untested) | ≈ 134,000,000 | ntau ≳ **224** |

The number of steps is what `native` mode writes, after `tau_stride`. Boost-invariant 2D
runs are far below this (a 150×150 Trento event with 301 steps is 6.8M cells), and so is
`grid` mode on the legacy 65×65×33 grid (3.8M cells).

To stay under 1 GB, the largest allowed number of written τ steps is
`floor(67,108,863 / (nx × ny × neta))`.

### Options

1. **Configuration only (works now).** Keep each event under 1 GB with:
   - `<tau_stride>` (with `N` stored steps, `tau_stride ≥ ceil(N / max steps)`; stride 2
     brings both O+O events above under 1 GB at 67 and 100 steps),
   - MUSIC's `output_evolution_every_N_timesteps` or `output_evolution_every_N_x/y/eta`,
   - or `grid` mode on a coarser grid.
2. **One tree entry per τ step.** Add event-index and step-index branches and fill one
   entry per τ step (9.6 MB per entry on the O+O grid). Removes both limits for any grid
   size. Changes the file layout, so readers (§4, `validate.py`, PyJetscape's
   `read_fast_root_bulk`, FNO training readers) must regroup rows by event.
3. **Plain float array with a length branch.** Store `user_res` as a leaf array
   (`user_res[n]/F`) instead of a `vector<float>`. This should avoid the per-object length
   field and keep one entry per event, so uproot readers barely change. It does not remove
   the ~2 GB limit, which entry 1 above nearly reaches.
4. **Warn in the writer.** When an event exceeds 1 GB, log its size and the smallest
   `tau_stride` that would fit. Doesn't fix anything, but makes the problem obvious.

**Recommendation:** option 2 if full-resolution `native` dumps on 3D grids are the goal;
otherwise option 1 plus the warning from option 4.
