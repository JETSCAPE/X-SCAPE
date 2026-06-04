# Hydrodynamics

## Physics motivation

The central discovery of the RHIC and LHC heavy-ion programs is that the QGP
behaves as a **nearly perfect fluid**: its shear viscosity to entropy density
ratio η/s is close to the conjectured quantum bound $1/4\pi$. Relativistic
**viscous hydrodynamics** is therefore the workhorse description of the bulk
evolution. It solves the conservation of the energy-momentum tensor (and, at
lower energies, the net-baryon current),

\[
\partial_\mu T^{\mu\nu} = 0, \qquad \partial_\mu N_B^{\mu} = 0,
\]

with $T^{\mu\nu} = \varepsilon\,u^\mu u^\nu - (p+\Pi)\Delta^{\mu\nu} +
\pi^{\mu\nu}$, closed by a lattice-QCD **equation of state** $p(\varepsilon)$ and
second-order (Israel–Stewart-type) relaxation equations for the shear
$\pi^{\mu\nu}$ and bulk $\Pi$ viscous stresses. Solving these on an
event-by-event fluctuating initial state turns initial **spatial** anisotropies
into final-state **momentum** anisotropies — the flow harmonics $v_n$ — whose
agreement with data is what pins down η/s.

For X-SCAPE the hydro stage plays a second role: it is the **medium that quenches
jets**. Every energy-loss query (`GetHydroCell`, see
[Signals](../framework/signals.md)) reads back the local temperature, flow, and
viscous tensor from the hydro evolution. The 4-D evolution is stored as an
`EvolutionHistory` of `FluidCellInfo` cells.

`FluidDynamics` (`src/framework/FluidDynamics.{h,cc}`) is the stage base class,
with the `HydroStatus` lifecycle (`NOT_START → INITIALIZED → EVOLVING →
FINISHED`), cell accessors, and the freeze-out surface interface.

## Modules

### MUSIC — (3+1)D viscous hydrodynamics

**XML:** `<Hydro><MUSIC>…</MUSIC></Hydro>` · **class:** `MpiMusic`
(`src/hydro/MusicWrapper.{h,cc}`, engine in `external_packages/music`, fetched by
`get_music.sh`, built with `-DUSE_MUSIC=ON`, run under MPI)

MUSIC is a (3+1)-dimensional relativistic **viscous** hydrodynamics code, the
default bulk evolution for full heavy-ion runs. It includes shear and bulk
viscosity, a lattice-QCD equation of state, optional net-baryon transport for
beam-energy-scan applications, and Cooper–Frye **freeze-out surface** extraction
(via Cornelius) for [particlization](particlization.md). It works in Milne
coordinates $(\tau, x, y, \eta_s)$.

A **GPU drop-in replacement, `music4gpu`**, runs the hydro kernels on
NVIDIA (CUDA) or Apple Silicon (Metal) devices; select it with
`-DUSE_MUSIC=ON -DUSE_CUDA=ON` (or `-DUSE_METAL=ON`). See
[Installation](../guides/installation.md) and `external_packages/music4gpu`.

*References:* Schenke, Jeon, Gale,
[arXiv:1004.1408](https://arxiv.org/abs/1004.1408); Paquet et al.,
[arXiv:1509.06738](https://arxiv.org/abs/1509.06738).

### CLVisc — GPU (OpenCL) (3+1)D viscous hydrodynamics

**XML:** `<Hydro><CLVisc>…` · **class:** `CLVisc`
(`src/hydro/CLViscWrapper.{h,cc}`, engine in `external_packages/clvisc_wrapper`,
built with `-DUSE_CLVISC=ON`)

CLVisc is a (3+1)D viscous hydro code written in **OpenCL** to run on GPUs (and
multi-core CPUs). It provides the same physics content as MUSIC with massive
parallel throughput, useful for large event ensembles.

*Reference:* Pang, Petersen, Wang,
[arXiv:1802.04449](https://arxiv.org/abs/1802.04449).

### Brick — static uniform medium

**XML:** `<Hydro><Brick>…` · **class:** `Brick`
(`src/hydro/Brick.{h,cc}`)

A **static, uniform "brick"** of QGP at a fixed temperature for a fixed
duration. Not a realistic fireball — its purpose is **controlled energy-loss
studies and validation**: every parton sees the same, known medium, so the
energy-loss formalism can be tested against analytic expectations. The standard
brick test (`brickTest`, `LBT_brickTest`) is the canonical quenching benchmark.

### Gubser — analytic flow

**XML:** `<Hydro><Gubser>…` · **class:** `GubserHydro`
(`src/hydro/GubserHydro.{h,cc}`)

An implementation of the **Gubser analytic solution** — a conformally symmetric,
azimuthally symmetric, radially expanding solution of ideal/viscous
hydrodynamics. Because it is exact, it is the standard **code-verification**
benchmark for hydro solvers and for the medium-interpolation machinery.

### HydroFromFile — replay a pre-computed evolution

**XML:** `<Hydro><HydroFromFile>…` · **class:** `HydroFromFile`
(`src/hydro/HydroFromFile.{h,cc}`)

Reads a tabulated hydro evolution (e.g. produced by an earlier MUSIC run, or by
an external code) from disk and replays it as the medium. This decouples the
(expensive) bulk evolution from the (cheap) jet evolution: compute a library of
hydro events once, then quench many jet samples through each. Several example
profiles can be downloaded with `examples/get_hydroSample*.sh`. Reads formats
including the freestream / IP-Glasma / MUSIC sample files in
`examples/test_*_files/`.

## Choosing a hydro module

| Goal | Use |
|---|---|
| realistic full event (default) | **MUSIC** (or **music4gpu** on GPU) |
| large ensembles on GPU | **CLVisc** |
| controlled jet-quenching test | **Brick** |
| hydro-code verification | **Gubser** |
| replay a stored medium for jet studies | **HydroFromFile** |

!!! tip "Reusing hydro"
    The bulk evolution dominates the cost of a full run while jets are cheap. The
    `<setReuseHydro>true</setReuseHydro>` / `<nReuseHydro>` settings (see
    [Configuration](../framework/configuration.md)) run many independent hard
    events on top of one hydro event.
