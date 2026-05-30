# X-SCAPE + music4gpu — Code Analysis & Engineering Roadmap

**Scope:** the X-SCAPE / JetScape C++ framework (`src/`, ~71k LOC) and the GPU hydro
port `music4gpu` (`external_packages/music4gpu/`, ~34k LOC, ~5.7k of it CUDA/Metal).

**Goal:** a prioritized engineering roadmap covering (1) memory issues, (2) design
improvements, (3) runtime improvements, and (4) GPU-acceleration opportunities beyond
the hydro stage. Every finding is grounded in the source with `file:line` references.

> This document is **analysis only** — it proposes no code changes. It is meant to
> drive future work. Each finding carries a severity, a rough effort, an impact, and a
> recommended next step; the master table in §7 sequences them.

**Legend**

| Field | Values |
|---|---|
| **Severity** | Low / Med / High (correctness or robustness risk) |
| **Effort** | S (< 1 day) · M (a few days) · L (weeks) |
| **Impact** | Low / Med / High (throughput, maintainability, or correctness payoff) |

**Method.** Findings come from direct reads of the files cited. The GPU-opportunity
section (§6) was verified against the actual hot loops. Memory claims (§3) are
static/code-evident — line numbers are the evidence; no runtime repro is implied.

---

## 1. Architecture at a glance

X-SCAPE is a **task-graph pipeline**: `JetScapeTask` nodes (`JetScapeTask.h`) are
executed recursively per event by `JetScape::Exec` (`JetScape.cc:1097-1243`). Modules
derive from `JetScapeModuleBase` and communicate through a sigslot
publish/subscribe layer (`JetScapeSignalManager`). The canonical stage order:

```
Initial state → Pre-equilibrium → Hydro (MUSIC / music4gpu) → Hard process
   → Jet energy loss (Matter/Martini/LBT/AdSCFT) → Hadronization → Afterburner → Writers
```

Two facts shape everything below:

1. **Execution is serial by default.** Module multithreading is hardcoded off
   (`multiTask = false`, `JetScape.cc:1162`, `JetEnergyLossManager.cc:214`). The only
   parallelism in `src/` is OpenMP in two spots (`ThermPtnSampler.cc:551,689`,
   `SurfaceFinder.cc:201,395`).
2. **The medium is the shared currency.** The hydro evolution is stored as a full 4D
   `EvolutionHistory` and read back by every downstream module through
   `GetHydroCellSignal` → `EvolutionHistory::get()`. Its storage and lookup cost
   (§3, §5) cut across the whole jet-quenching workload.

`music4gpu` plugs in at the hydro stage. It keeps three grid **snapshots**
(prev/curr/future) device-resident and runs seven hydro kernels per RK substep, with a
CPU fallback whenever the configuration is outside its support matrix
(`advance.cpp:347-370`). It has two memory back-ends chosen at runtime: **coherent**
(unified `cudaMallocManaged`, zero-copy — integrated / NVLink-C2C parts) and
**discrete** (device `cudaMalloc` + pinned staging + explicit DMA). It also has a Metal
back-end for Apple Silicon, mirroring the CUDA path one-to-one.

---

## 2. What's already done well (so we don't regress it)

Worth stating explicitly, because the roadmap should preserve these:

- **GPUGrid is RAII** (`~GPUGrid()` → `release()`, `GPUGrid.h:46`) with a clean
  free sweep that nulls every pointer (`GPUGrid_cuda.cu:198-244`). No leaks found on
  the nominal path.
- **Graceful degradation.** Device-absent, alloc-failure, EOS-upload-failure, and
  unsupported-config cases all fall back to CPU instead of crashing
  (`advance.cpp:46-57`, `:347-370`).
- **Transfer-skipping residency optimization.** When the GPU owns the authoritative
  state across a substep/step boundary, the host→device upload is skipped
  (`advance.cpp:387`), which *also* avoids re-truncating every cell fp64→fp32 each
  step (the rationale at `advance.cpp:383-386` notes eps_max drift drops from ~1e-2 to
  ~1e-4 over 100 steps). This is the single most important perf+accuracy design choice
  in the port — keep it.
- **No per-step allocation churn.** Snapshots and scratch are allocated once in
  `allocate()`; the only lazy (re)allocation is `evo_pack_out`, and it only grows
  (`CUDAPipelines.cu:266`).
- **Allocation failures are reported**, not silent — `alloc_*_buf` log the byte count
  and CUDA error string (`GPUGrid_cuda.cu:55-58, 69-72, 83-86`).

---

## 3. Memory

### 3.1 `buf_handles_[32]` is a fixed array with no bounds check — **Sev High (latent) · Effort S**

`GPUGrid` tracks every buffer to free in a raw C array:

```cpp
// GPUGrid.h:203-208
// 3 snapshots × 5 fields = 15
//   + dwmn + qi_out + uwrhs_out + theta_buf + a_buf + sigma_buf = 21
//   + eos_P + eos_dPde = 23
// Leave headroom for upcoming Phase-2 buffers (eos_T, eos_s, ...).
void* buf_handles_[32];
int   n_handles_ = 0;
```

The allocators append with **no capacity check**:

```cpp
// GPUGrid_cuda.cu:60 (and :74 for the device variant)
handles[n_handles++] = ptr;
```

The actual live count today is **29 of 32**:

| Group | Buffers | Count |
|---|---|---|
| Snapshots (prev/curr/future × {epsilon, rhob, u, Wmunu, pi_b}) | 3 × 5 | 15 |
| Scratch (dwmn, qi_out, uwrhs_out, uprhs_out, qi_source_buf, theta_buf, a_buf, sigma_buf) | | 8 |
| Reductions (reduce_eps_out, reduce_rhob_out) | | 2 |
| EOS (eos_P, eos_dPde, eos_s, eos_T) | | 4 |
| **Total** | | **29** |

The explanatory comment is **stale** (it stops at 23 and predates `uprhs_out`,
`qi_source_buf`, the two reduction scalars, and `eos_s`/`eos_T`). Only **3 slots**
remain. Adding ~3 more managed buffers — entirely plausible given the steady "Tier 3c
Phase N" growth visible in `gpu_types.h` — overflows the array and corrupts the
adjacent `GPUGrid` members (`n_handles_`, `Nx_`, the snapshot structs) with no
diagnostic.

**Recommended:** replace `void* buf_handles_[32]` with `std::vector<void*>` (push_back
in `alloc_*_buf`, sweep + `clear()` in `release()`). If a fixed array is preferred for
the Metal `__bridge` reverse-lookup, at minimum add
`assert(n_handles_ < kMaxHandles)` in both allocators and bump the comment.
*Test:* a unit allocate/release cycle; verify `n_handles_` after a full `allocate()`.

### 3.2 Grid index arithmetic is 32-bit — **Sev Low · Effort S**

```cpp
// GPUGrid_cuda.cu:128
Ncells_ = Nx * Ny * Neta;           // int * int * int

// GPUGrid_cuda.cu:39
static inline int cell_idx(int ix, int iy, int ieta, int Nx, int Ny) {
    return Nx * (Ny * ieta + iy) + ix;   // int result
}
```

`MUSICGridParams.Ncells` is likewise `int` (`gpu_types.h:23`). The product overflows
`INT_MAX` (2.1B) above ~1290³ cells (e.g. a 2048³ grid). Production MUSIC grids are
typically ≤ 512³ (134M), so this is **latent**, not active — but it's a silent
wraparound if someone scales up, and the per-component stride `m * Ncells_ + c`
(`GPUGrid_cuda.cu:272`, `music_kernels` SoA layout) overflows even sooner for the
14-component Wmunu field.

**Recommended:** compute `Ncells_` in `size_t`, return `size_t` from `cell_idx`, and
make the SoA stride `size_t`. *Test:* compile-time; optionally a (skipped) large-grid
allocation guard.

### 3.3 CUDA copies are unchecked — **Sev Med · Effort S**

The synchronous D2H copy-backs and the async H2D upload do not check their return
codes:

- D2H `cudaMemcpy`: `GPUGrid_cuda.cu:287, 290, 332-336, 407, 410, 436-440, 470`
- Async H2D `cudaMemcpyAsync`: `CUDAPipelines.cu:210`

A failed transfer (e.g. an unmapped pointer after a partial alloc failure, or an ECC
fault) propagates as **silently wrong physics** rather than an error. Kernel *launches*
are checked (`check_launch`, `CUDAPipelines.cu:46-52`), so the copies are the gap.

**Recommended:** a single `CUDA_CHECK(expr)` macro (log file/line + `cudaGetErrorString`,
set `gpu_ready_ = false` on failure) applied to every copy, `cudaEventRecord`,
`cudaMemsetAsync`, and stream/event create/destroy. This unifies with the design item
in §4.4. *Test:* fault injection via `MUSIC_CUDA_FORCE_DISCRETE=1` plus an artificially
small grid.

### 3.4 The full hydro history lives in host RAM — **Sev Med · Effort M**

`EvolutionHistory` stores the entire 4D medium as Array-of-Structs:

```cpp
// FluidEvolutionHistory.h:168
std::vector<FluidCellInfo> data;
```

`FluidCellInfo` (`FluidCellInfo.h:60-71`) is ~28 `Jetscape::real` fields
(energy/entropy/temperature/pressure/qgp_fraction, μ_B/μ_C/μ_S, vx/vy/vz, the full
`pi[4][4]` shear tensor = 16 reals, bulk_Pi). With `real = double` that's **~224 bytes
per cell**. Footprint:

```
bytes ≈ ntau · nx · ny · neta · 224
```

A 130×130×32 grid over ~60 τ-steps is ≈ 7 GB; a 261×261×64 grid over a few hundred
τ-steps runs into the tens-to-hundreds of GB. More than half of every cell
(`pi[4][4]`, 128 B) is dead weight for the many consumers that read only
`e, T, vx, vy, vz`. The history is built one heap `unique_ptr<FluidCellInfo>`
at a time in `MusicWrapper.cc:649-684`.

**Recommended:** (a) store in `float` unless double is needed downstream; (b) consider
SoA + a "with-shear / without-shear" split so jet modules that need only the flow
fields don't pay for `pi[4][4]`; (c) for music4gpu runs, the device already holds the
grid — a streaming or down-sampled history path would avoid the full materialization.
*Test:* peak-RSS comparison on a 3+1D event before/after.

### 3.5 `FluidCellInfo` is returned by value through the hottest function — **Sev Low (correctness) / see §5.2 for runtime**

`EvolutionHistory::get()` is **O(1)-indexed** (uniform-grid arithmetic, no search — good),
but each call makes **16 by-value `FluidCellInfo` copies**: two `GetAtTimeStep` calls
(tau, tau+1), each doing eight `GetFluidCell` returns-by-value, then trilinear + linear
interpolation (`FluidEvolutionHistory.cc:313-372`, `GetFluidCell` returns by value at
`FluidEvolutionHistory.h:452`). At ~224 B/struct that's ~3.5 KB copied per lookup. The
correctness side is fine; the cost side is in §5.2 because the function is called
per-parton per-microstep.

### 3.6 Parton/Hadron ownership — documented, no leak found

Energy-loss modules hold `vector<Parton>` by value and pass them by **reference**
(`JetEnergyLoss.h:278`, `Matter.h` `DoEnergyLoss(... vector<Parton>&, vector<Parton>&)`),
so there's no double-copy in the inner pipeline. The `PartonShower` is a GTL graph of
`shared_ptr<Parton>` kept alive for the event and cleared per event. No unbounded
growth was found. Worth documenting the model so future changes don't accidentally copy
shower nodes.

---

## 4. Design improvements

### 4.1 CUDA ⇄ Metal kernel duplication — **Sev Med (maintenance) · Effort L**

The per-cell hydro physics is implemented **twice**:

| CUDA | Metal | What |
|---|---|---|
| `music_kernels.cu` (~1,965 L, 77 KB) | `music_kernels.metal` (~2,191 L, 92 KB) | the seven hydro kernels |
| `CUDAPipelines.cu` | `MetalPipelines.mm` | launch/dispatch |
| `GPUGrid_cuda.cu` | `GPUGrid.mm` | buffer management |

The two kernel files are near-identical physics (same SoA layout, same stencils, same
Wmunu index map). Every physics fix — a stencil sign, a reconstruction tweak, a new
viscous term — must be applied in both, and the two will drift. `gpu_types.h` already
shows the seam: it shares the struct/constant definitions across C++/CUDA/Metal with
`#ifdef __METAL_VERSION__` / `__CUDACC__` shims (`gpu_types.h:151-175`).

**Recommended:** extend that shim approach to the kernel *bodies* — factor the per-cell
math into shared headers with thin `DEVICE_FN` / address-space macros, so CUDA and
Metal instantiate the same source. This is a large refactor; sequence it after the
quick wins. *Test:* GPU-vs-CPU parity (§8) must be byte-stable across the refactor.

### 4.2 Surface extraction shells out — **Sev Med · Effort M**

```cpp
// MusicWrapper.cc:547-568  (collect_freeze_out_surface)
system("rm surface.dat 2> /dev/null");
...
system_command << "cat surface_eps* >> " << surface_filename.str();
system(system_command.str().c_str());
...
system_command << "ln -s " << surface_filename.str() << " surface.dat";
system(...);
system("rm surface_eps* 2> /dev/null");
```

MUSIC writes per-slice `surface_eps*` files to disk, which are then concatenated,
symlinked, and deleted via four `system()` calls. Problems: non-portable (POSIX shell,
symlinks), **error-blind** (no check that `cat` succeeded or that any file existed),
fragile to the working directory, and **racy** for concurrent events sharing a CWD.
An in-memory path already exists — `PassHydroSurfaceToFramework` reads cells directly
via `get_surface_cell_with_index` (`MusicWrapper.cc:599-632`).

**Recommended:** assemble the surface in C++ (or route the sampler through the in-memory
surface) and drop the shell round-trip. *Test:* compare sampled-surface particle yields
before/after on a fixed seed.

### 4.3 `RootBulkWriter` carries acknowledged correctness debt — **Sev Med · Effort M**

The writer's own comments flag the issues (`RootBulkWriter.cc:243-255`):

- τ-bin count `ntau_freezeout` "ends up being about 2 units too large (why?!?)" (`:255`)
- empty energy-density steps past freeze-out (`:253`, `:257`)
- `tau_min` (user) vs `_tau_min` (MUSIC) mismatch with a `// NOTE: FIXME` to use the
  latter (`:253`)

These produce nonsensical first/last bins and inconsistent bin counts across events.

**Recommended:** pin down the intended (tau_min, dtau, ntau) contract between MUSIC's
grid and the writer, fix the binning, and delete the dead/commented debug scaffolding.
*Test:* a fixed event's bulk tree should have exactly `ntau` populated, physically
sensible steps.

### 4.4 Inconsistent GPU error handling — **Sev Med · Effort S**

Today: kernel launches are checked (`check_launch`), allocations are checked inline,
but copies/events are not (§3.3). Adopt one `CUDA_CHECK` convention across all
`*_cuda.cu` / `CUDAPipelines.cu` call sites. (Same work item as §3.3, listed here as the
design/consistency view.)

### 4.5 Land the working-tree cleanup — **Sev Low · Effort S**

`git status` shows `EPGun.{cc,h}` and `JetScapeWriterRootHepMC.{cc,h}` deleted but not
committed, with matching edits in `CMakeLists.txt` / `JetScape.cc` / `RootBulkWriter.cc`.
This looks like intentional dead-code removal mid-flight. **Recommended:** finish and
commit the removal (or revert) so the tree isn't left half-deleted; confirm no remaining
references in the build. *Test:* clean configure + build of both `build_gpu` and
`build_lite`.

---

## 5. Runtime improvements

### 5.1 The pipeline is serial; event-level parallelism is the biggest lever — **Impact High · Effort M**

`multiTask = false` is hardcoded (`JetScape.cc:1162`, `JetEnergyLossManager.cc:214`), so
each event runs every stage serially, and events run one at a time. Heavy-ion runs are
**embarrassingly parallel across events** — the cleanest throughput win is to run N
independent events concurrently (process- or thread-level), which needs no physics
changes. The dormant `multiTask` path (per-module `std::thread`, gated on
`GetMultiThread()`) is a secondary, finer-grained option but only parallelizes the
`CalculateTime()` phase and is currently untested.

**Recommended:** prioritize an event-parallel driver (or document the existing
launcher's parallelism). Treat `multiTask` as a separate, validated follow-up. *Test:*
wall-clock scaling vs. core count at fixed total events; identical per-event output.

### 5.2 `GetHydroCellSignal` copy overhead on the hottest path — **Impact Med · Effort S/M**

Every energy-loss module queries the medium **per parton, per microstep**:
`Matter.cc:772` (inside the `el_time` micro-step loop), `Martini.cc:201`, plus the
splitting path. Each query is a `get()` → 16 `FluidCellInfo` by-value copies (~3.5 KB)
+ interpolation (§3.5). For O(10²–10³) partons × O(10–10²) microsteps this dominates the
medium-interaction cost.

**Recommended:** add an out-param / const-ref `get()` overload that fills a caller-owned
`FluidCellInfo` (no return-by-value), reuse a single cell object across the microstep
loop, and cache the τ-bracket (`id_tau`, `tau0`, `tau1`) so consecutive microsteps at
the same τ skip the bracket search. *Test:* `MUSIC_PROFILE`/`perf` on a brick or full
event; bit-identical parton output.

### 5.3 `HadronicEMT` grid loops are scalar — **Impact Med · Effort S (now), High (GPU, see §6)**

`DetermineFluidization` is a triple-nested grid loop whose body loops over all hadrons —
O(N_grid · N_hadron) — and is fully scalar (`HadronicEMT.cc:680-727`). `GetBulkInfo`'s
smearing loop is likewise scalar (`:541-586`).

**Recommended (cheap, now):** OpenMP `collapse(3)` over the (ix, iy, iz) grid — the exact
pattern already used in `SurfaceFinder.cc:201` and `:395`. This is a near-free
multi-core speedup before any GPU work. **Also flags a correctness TODO:** the hydro is
queried at a hardcoded `t = 0.0` with a `// TODO: This needs the hydro medium info at
the current time step` (`HadronicEMT.cc:690-691`) — worth resolving. *Test:* identical
fluidization flags single- vs multi-threaded; resolve the t=0 TODO separately.

### 5.4 Residual GPU↔CPU sync points — **Impact Low/Med · Effort M (document first)**

- `reduce_max` issues a full `cudaStreamSynchronize` every timestep
  (`CUDAPipelines.cu:247`; called `evolve.cpp:298, 577`) to read the two scalars the CPU
  needs for the max-energy / freeze-out decision. It moves only 8 bytes but forces a
  per-step CPU↔GPU barrier.
- Freeze-out and output-evolution force a bulk D2H each step they fire, already mitigated
  by the `host_curr_fresh_` / `host_prev_fresh_` caching that collapses a same-step
  freeze-out + output into one copy (`advance.cpp:488-514`).

**Recommended:** document these as the intrinsic serialization of the current design;
explore overlapping the reduction with the next substep's kernels (the
`cudaEventDisableTiming` copy-stream machinery already exists) only if profiling shows
it matters. *Test:* per-step timeline (Nsight) to quantify the stall before investing.

### 5.5 `RootBulkWriter` interpolation loop — **Impact Low/Med · Effort S**

The 4-level loop calls `bInfo.get()` (full trilinear+linear interpolation, §3.5) for
every (τ, x, y, η) cell (`RootBulkWriter.cc:261-321`). For training-data dumps this is
the writer's whole cost. It's I/O+interpolation bound and only runs at output cadence,
so it's lower priority — but the per-cell `get()` shares the §5.2 win, and a direct
grid-walk (no re-interpolation when the output grid matches the hydro grid) would avoid
it entirely.

---

## 6. GPU acceleration beyond hydro

The hydro stage is ported. Below are the next candidates, **verified against the actual
loops**, ranked by impact ÷ effort. Each notes the parallelism pattern, data size, and a
kernel sketch.

### 6.1 ★ TOP TARGET — `HadronicEMT` particle→grid smearing — **Effort M · Impact High**

**Where:** `GetBulkInfo` smearing (`HadronicEMT.cc:541-586`) and `DetermineFluidization`
(`:680-727`). **Pattern:** classic particle-to-grid *scatter/deposit* — for each grid
cell, sum a Gaussian (or covariant) smearing kernel over hadrons within a 5σ box, into a
4×4 `T^{μν}`. **Why it's the best first target:** deterministic (no RNG), cache-friendly,
fixed small per-cell output, and the work today is O(N_grid · N_hadron) entirely on one
core. Data is tiny: O(10³–10⁴) hadrons × a few floats, uploaded once per timestep.

*Kernel sketch:*
```
upload hadron {x,y,z, p, m} array once per timestep
kernel<thread per grid cell>:
    accumulate T^{mu nu} by gathering hadrons in the 5σ neighborhood
    (gather variant → no atomics; or thread-per-hadron + atomicAdd into T^{mu nu})
fuse the fluidization-flag pass into the same kernel (it reuses the same neighborhood test)
```
Map the smearing kernel to a `__device__` function; the 5σ skip box becomes the
thread's loop bound. The same launch handles §5.3's grid loop. *Recommended as #1
non-hydro GPU target.* Start with the free OpenMP step (§5.3), then lift to CUDA.
*Test:* compare `T^{μν}` and fluidization flags vs. the CPU path on a fixed hadron list.

### 6.2 `Matter` elastic energy loss — **Effort L · Impact High**

**Where:** per-parton loop `Matter.cc:298`; inner micro-step trajectory `:715-…`;
per-step medium lookup `:772`; elastic collisions `colljet22`/`collHQ22` `:905-907`.
**Pattern:** embarrassingly parallel over partons, each integrating a micro-step
trajectory with RNG-driven rejection sampling, Lorentz boosts to/from the fluid rest
frame, and a per-step hydro interpolation. **Challenges:** per-parton RNG state
(cuRAND/Philox), branchy control flow (process-type and flavor branches), and the medium
lookup — which is why §6.6 (GPU-resident grid) is a prerequisite for real speedup rather
than PCIe thrash. **Impact** is high because jet quenching over many partons × microsteps
is a dominant cost.

*Kernel sketch:* one parton per thread (or warp), per-thread cuRAND state, medium read
from a device-resident hydro grid (§6.6), trajectory loop in-kernel; collect daughters/recoils into a compacted output buffer.

### 6.3 `Martini` radiative/elastic — **Effort M/L · Impact Med-High**

**Where:** per-parton loop `Martini.cc:160`; `DetermineProcess` `:293`;
`getNewMomentumRad`; rate tables via `use_table(p,k,…)` over NP×NK grids
(`:2926`+). **Pattern:** fully independent per parton; RNG + **2D bilinear rate-table
interpolation** that maps naturally to CUDA **constant/texture memory**. Lower
trajectory complexity than Matter (single-step decision per call), so a cleaner first
quenching kernel once the RNG/medium infrastructure from §6.2/§6.6 exists.

### 6.4 `LBT` — **Effort L · Impact Med**

Per-parton parallel at the outer level (`LBT.cc:206`+), but `LBT0`'s internals are
sequential rejection-sampling loops with heavy dependence on large rate tables. Suitable
only after Matter/Martini establish the per-parton RNG + table-in-device-memory pattern;
the sequential inner structure limits per-parton speedup.

### 6.5 Initial state & pre-equilibrium

- **IP-Glasma** color-field initialization is **FFT-based** → `cuFFT` is a natural,
  well-trodden GPU win (external package; needs work in `external_packages/ipglasma`).
- **MC-Glauber / TRENTo** are sequential Monte-Carlo samplers → Low GPU suitability.
- **CLVisc** is **already OpenCL/GPU** (`clvisc_wrapper/clvisc.cc:322-339`). The item
  here is *porting it to CUDA/HIP* so the pre-eq and hydro stages share one GPU stack
  with music4gpu (avoids OpenCL+CUDA coexistence and a CPU handoff).

### 6.6 ★ Cross-cutting enabler — a GPU-resident hydro grid — **Effort L · Impact High**

The recurring blocker for §6.2–§6.4 is that each parton microstep interpolates the
medium on the CPU (`GetHydroCellSignal`). music4gpu **already keeps the grid on the
device**. Exposing a read-only, device-resident hydro-cell sampler (texture-backed
trilinear interpolation) that jet kernels can call would let Matter/Martini/LBT/HadronicEMT
read the medium without per-microstep host round-trips. This is the highest-leverage GPU
*infrastructure* investment and a natural extension of the existing device grid — build
it alongside §6.1 and it unlocks the rest.

---

## 7. Prioritized roadmap

Sequenced by impact/severity ÷ effort. **Do the top band first** — all are S-effort with
real payoff.

| # | Item | § | Sev/Impact | Effort | Recommended action |
|---|---|---|---|---|---|
| 1 | `buf_handles_` bounds safety | 3.1 | High (latent) | S | `std::vector<void*>` or capacity assert |
| 2 | Unified `CUDA_CHECK` on copies/events | 3.3·4.4 | Med | S | one macro across `*_cuda.cu`/`CUDAPipelines.cu` |
| 3 | OpenMP `collapse(3)` on HadronicEMT grid | 5.3 | Med | S | reuse `SurfaceFinder.cc:201` pattern |
| 4 | `get()` by-ref + τ-bracket cache | 5.2·3.5 | Med | S/M | out-param overload; reuse cell object |
| 5 | 64-bit grid index arithmetic | 3.2 | Low | S | `size_t` for `Ncells`/`cell_idx`/stride |
| 6 | Land working-tree cleanup | 4.5 | Low | S | finish/commit EPGun + RootHepMC removal |
| 7 | De-shell surface extraction | 4.2 | Med | M | in-memory / C++ assembly |
| 8 | Fix `RootBulkWriter` τ-binning | 4.3 | Med | M | pin grid contract; fix off-by-~2 |
| 9 | Hydro-history footprint (float/SoA) | 3.4 | Med | M | float storage; shear split |
| 10 | Event-level parallelism | 5.1 | High | M | concurrent independent events |
| 11 | **GPU: HadronicEMT smearing kernel** | 6.1 | High | M | ★ top GPU target |
| 12 | **GPU: device-resident hydro sampler** | 6.6 | High | L | ★ unlocks jet-loss GPU work |
| 13 | GPU: Matter / Martini per-parton loss | 6.2·6.3 | High | L | after #12; cuRAND + tables in const/texture |
| 14 | Single-source CUDA/Metal kernels | 4.1 | Med | L | shared headers + macro shims |

**Top GPU target:** **HadronicEMT particle→grid smearing (#11)** — highest impact for
the effort, no RNG, and it pairs with the device-resident hydro sampler (#12) that
unlocks everything else.

---

## 8. Verification — how to validate findings and any fix derived from them

This document changes no code. To confirm a claim or test a future fix:

- **GPU↔CPU parity & profiling.** Run `runJetscape OO_one_event.xml` with
  `MUSIC_PROFILE=1`. Compare `MUSIC_FORCE_CPU` vs. GPU for eps_max drift (validates the
  residency rationale at `advance.cpp:383-386`). Exercise both memory back-ends with
  `MUSIC_CUDA_FORCE_DISCRETE=1` and `MUSIC_CUDA_FORCE_COHERENT=1`
  (`CUDAPipelines.cu:115-134`).
- **Builds.** Both `build_gpu` (USE_CUDA) and `build_lite` trees exist; any fix must
  compile-check against both, and the §4.5 cleanup against a clean configure.
- **Static memory items (§3.1, §3.2).** Code-evident — the cited lines are the evidence.
  For §3.1, print `n_handles_` after `allocate()` to confirm 29/32.
- **Runtime items.** Use `MUSIC_PROFILE` / `perf` / Nsight for §5.2 and §5.4; wall-clock
  scaling for §5.1; peak-RSS for §3.4.
- **Physics-touching changes (§4.2, §4.3, §6.x).** Require fixed-seed, fixed-event
  output comparison (parton/hadron yields, surface, bulk tree) before/after.

---

### Appendix — primary files referenced

- **GPU layer:** `external_packages/music4gpu/src/gpu/{GPUGrid.h, GPUGrid_cuda.cu, CUDAPipelines.cu, gpu_types.h, music_kernels.cu, music_kernels.metal, GPUGrid.mm, MetalPipelines.mm}`, `external_packages/music4gpu/src/{advance.cpp, evolve.cpp}`
- **Integration & writers:** `src/hydro/MusicWrapper.cc`, `src/root/RootBulkWriter.cc`, `src/framework/{JetScape.cc, JetEnergyLossManager.cc, FluidEvolutionHistory.cc, FluidEvolutionHistory.h, FluidCellInfo.h, FluidDynamics.h}`
- **GPU candidates / parallelism:** `src/hadronicEMT/HadronicEMT.cc`, `src/jet/{Matter.cc, Martini.cc, LBT.cc}`, `src/hadronization/{HybridHadronization.cc, ThermPtnSampler.cc}`, `src/framework/SurfaceFinder.cc`, `external_packages/clvisc_wrapper/clvisc.cc`
