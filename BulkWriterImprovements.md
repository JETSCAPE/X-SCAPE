# RootBulkWriter — performance improvements & a "fast" framework-bypassing ROOT path

> Status: **implemented and building clean** in `build_gpu`; verified on O+O 1-event
> (see Verification). Usage: [README_BulkFast.md](README_BulkFast.md).

## Context

`RootBulkWriter` dumps the MUSIC hydro evolution (energy density + flow velocities) to a
ROOT file for ML/training use (e.g. the FNO pipeline). Two asks:

1. **Q1** — flesh out the §5.5 suggestions in `XSCAPE-CodeAnalysis.md` for speeding up the
   writer.
2. **Q2** — is a *fast* path possible that creates the ROOT hydro-evolution file
   **without copying the evolution into the JETSCAPE framework** first?

Reference workflow: the XML-steered `runJetscape` executable in `build_gpu/`, driven by
`build_gpu/OO_one_event.xml`, producing `OO_test.root`.

User decisions: a **new, separate writer module**; **hydro-only** scope (skip the
framework history); **two grid modes** — supply your own grid (interpolate only if it
doesn't match MUSIC's), or use MUSIC's grid and just thin the timesteps.

---

## How the data flows today (what "the copy into the framework" is)

```
MUSIC run_hydro()
  └─ per output timestep: OutputEvolutionDataXYEta_memory(arena, tau)   [grid_info.cpp:372]
        loops ix → iy → ieta → HydroinfoMUSIC::dump_ideal_info_to_memory  [HydroinfoMUSIC.cpp:373]
        → pushes fluidCell_ideal onto lattice_ideal (a std::vector)       [HydroinfoMUSIC.h:33]
        NATIVE STORE: tau-major, then x, then y, then eta.

MpiMusic::EvolveHydro()                                                  [MusicWrapper.cc:452]
  └─ PassHydroEvolutionHistoryToFramework()                             [MusicWrapper.cc:649]
        per native cell: new FluidCellInfo + copy 28 fields
        → bulk_info.data.push_back(*ptr)   (framework AoS, ~224 B/cell)
        then clear_hydro_info_from_memory()         ←—— "the copy into the framework"

RootBulkWriter::Exec()                                                   [RootBulkWriter.cc:246]
  const auto &bInfo = hydro.lock()->get_bulk_info();   // was `auto` → whole-history copy (FIXED)
  4-nested loop over the OUTPUT grid → bInfo.get(tau,x,y,eta)  // trilinear+linear interp
  t->Fill();
```

`EvolutionHistory::CellIndex` (`FluidEvolutionHistory.h:429`) =
`id_tau*nx*ny*neta + id_x*ny*neta + id_y*neta + id_eta` — the **same ordering** as MUSIC's
native `lattice_ideal`. So a grid-walk over either container is an identity map;
interpolation is only needed when the output grid differs from MUSIC's grid.

---

## Findings

- **A1 (fixed):** `RootBulkWriter::Exec` used `auto bInfo = …get_bulk_info()`, which
  deep-copied the **entire** history each event (`get_bulk_info()` returns
  `const EvolutionHistory&`). Changed to `const auto &bInfo`. One line, unconditional.
- **A2:** each `bInfo.get()` does 16 by-value `FluidCellInfo` copies (~3.5 KB) +
  interpolation per output cell; at a grid node this is an identity at full cost.
- **A4:** for a hydro-only dump, building `bulk_info.data` and reading it back is 100%
  overhead (and the dominant peak-RAM term for 3D grids — CodeAnalysis §3.4).

---

## Part A — existing `RootBulkWriter` (the §5.5 work)

- **Done:** A1 const-reference fix (`RootBulkWriter.cc:255`).
- **Superseded by Part B:** the §5.5 grid-walk and §5.2 out-param/τ-bracket `get()` are
  implemented in the new fast writer rather than retrofitted onto the legacy interpolation
  hot path (still used by jet quenching), to avoid destabilising it. The legacy writer
  keeps its exact behaviour apart from A1.

## Part B — new `FastRootBulkWriter` module (hydro-only, bypasses the copy)

Reads MUSIC's native in-memory store directly; never builds `bulk_info.data`.

- **B1 — module:** `src/root/FastRootBulkWriter.{h,cc}`, self-registered as
  `"FastRootBulkWriter"`, added to `JetScape::DetermineTaskListFromXML` and to the
  `USE_ROOT` sources in `src/CMakeLists.txt`. ROOT tree layout mirrors the legacy writer
  (`user_res` vector + `TParameter` metadata, `nFeatures=4`: energy_density, vx, vy, vz).
- **B2 — native access on MpiMusic** (`MusicWrapper.h`): inline
  `get_number_of_fluid_cells()`, `clear_hydro_info_from_memory()`, and
  `get_native_fluid_cell(idx, FluidCellInfo&)` (reuses MUSIC's u^μ→v + hbarc conversion so
  values match the framework path), plus `set/get_dump_hydro_only`. *Deviation from plan:*
  thin pass-throughs on `MpiMusic` (not new virtuals on `FluidDynamics`); the writer
  reaches them via `dynamic_pointer_cast<::MpiMusic>` (note: `MpiMusic` lives in the
  **global** namespace, so it must be forward-declared / cast as `::MpiMusic`).
- **B3 — hydro-only switch.** *Deviation from plan:* driven by XML
  `<Hydro><MUSIC><dump_hydro_only>1`, read in `MpiMusic::InitializeHydro` (robust to module
  init ordering). When set, `EvolveHydro` skips `PassHydroEvolutionHistoryToFramework`,
  calls `SetHydroGridInfo()` (grid scalars only), and retains the native store; it also
  clears the native store defensively at the top of each event. The writer clears it again
  at the end of its `Exec`.
- **B4 — two modes** (`<grid_mode>` = `native` | `grid`):
  - `native`: walk MUSIC's grid, optional `<tau_stride>` to thin τ. Zero interp, zero AoS.
  - `grid`: user grid; builds a **transient** local `EvolutionHistory` from the native
    store once per event and calls `EvolutionHistory::get()` so output matches the legacy
    writer. A transient working copy (freed each event), not the *persistent* framework AoS.
- **B5/B6 — XML + metadata:** see README_BulkFast.md. Metadata records effective
  `tau_min`/`dtau` (incl. `tau_stride`) and `grid_mode`.

### Decisions on the former open questions
- **Native-store lifetime:** writer clears at end of `Exec`; `MpiMusic` clears defensively
  at the top of `EvolveHydro`.
- **Pre-equilibrium:** excluded (hydro-only). `include_preequilibrium` is reserved (not
  wired). The fast dump starts at MUSIC's `hydroTau0`, which can differ from the legacy
  writer's configured `tau_min`.
- **boost-invariant / 2+1D:** native mode emits MUSIC's `neta`; grid mode reuses
  `EvolutionHistory` behaviour.

## Files changed
- `src/root/RootBulkWriter.cc` — A1 const-ref fix.
- **new** `src/root/FastRootBulkWriter.{h,cc}`.
- `src/hydro/MusicWrapper.h` / `.cc` — native accessors + `dump_hydro_only_`; `EvolveHydro`
  retains/clears native store, skips the AoS copy.
- `src/framework/JetScape.cc` — include + task-list registration.
- `src/CMakeLists.txt` — add `root/FastRootBulkWriter.cc` under `USE_ROOT`.
- `config/jetscape_main.xml` — default entries for `<dump_hydro_only>` and
  `<FastRootBulkWriter>` (required by the XML validator).
- `config/BulkFastTest/` — example configs, README, and `validate.py`.

## Verification — ✅ passing (O+O, 1 event, build_gpu)
Configs in `config/BulkFastTest/`. Run from `build_gpu/`, e.g.
`./runJetscape ../config/BulkFastTest/OO_one_event_fastgrid.xml`. Validate with
`config/BulkFastTest/validate.py` under a numpy+uproot Python (conda `fno_env`; the system
`python3` here has no numpy).

1. Build: `cmake build_gpu` (re-configure) then
   `cmake --build build_gpu --target runJetscape -j4` — clean; `FastRootBulkWriter.cc.o`
   compiled and symbols present in `libJetScape`.
2. All three configs run `EXIT=0`; log shows
   `MUSIC dump_hydro_only: retaining native store (… cells), skipping framework copy`, so
   `bulk_info.data` is never built (B3 works).
3. **`grid` mode vs legacy `RootBulkWriter`: bit-for-bit identical** (uproot+numpy):
   same shape (both `ntau=27`, 15,057,900 floats on 65×65×33) and `np.array_equal == True`
   (`max |abs diff| = 0`). Files ≈ 16.83 MB. Confirms the fast path reproduces the legacy
   interpolation exactly on a matching grid.
4. **`native` mode:** MUSIC's full grid (here 100×100×60, ntau=134 at `dtau≈0.02`); file
   ≈ 93 MB; no NaN/inf; peak energy 7.48 vs legacy 7.36 (native keeps full resolution, so
   it resolves the peak slightly higher — expected). Thin via `tau_stride` / MUSIC
   `output_evolution_every_N_*`.

### Two bugs found & fixed during verification
- **CMake source not added (root cause of an early "no file"):** the `src/CMakeLists.txt`
  edit to add `root/FastRootBulkWriter.cc` initially did not land, so the module wasn't
  compiled, `RegisterJetScapeModule` never ran, `createInstance` returned null, and the
  writer was silently absent from the task list. Fixed; re-configure CMake after editing
  CMakeLists and confirm the `.o` exists + symbol in `libJetScape`.
- **Main-XML defaults required:** every new user-XML tag needs a default entry in
  `config/jetscape_main.xml` or X-SCAPE aborts ("tag is unrecognized",
  `JetScape::recurseToSearch`).

> Also: remove `<Eloss>`/hadronization/afterburner from the user XML (they query the
> now-empty `bulk_info` and would throw).
