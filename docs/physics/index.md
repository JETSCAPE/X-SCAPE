# Physics Modules

X-SCAPE assembles a complete heavy-ion (or p–p / p–A / e–A) event from
interchangeable physics modules, one or more per stage. This section gives, for
each stage, the **physics motivation**, the **available modules**, their key
parameters, and the **primary literature**.

Start with **[The Multistage Picture](overview.md)** for the big picture of why
the simulation is staged the way it is.

## Stages and modules

| Stage | Modules (XML name) | Page |
|---|---|---|
| **Initial state** | `Trento`, `MCGlauber` (3DGlauber), `IPGlasma`, `InitialFromFile`, `NcollListFromFile`, SMASH-IC, particle guns (`PGun`, `PythiaGun`, `EPGun`, `epemGun`) | [Initial State](initial-state.md) |
| **Initial-state radiation** | `PythiaIsrGun` + iMATTER (`InitialStateRadiationTest`) | [ISR (iMATTER)](isr.md) |
| **Pre-equilibrium** | `FreestreamMilne`, `Glasma`, `NullPreDynamics` | [Pre-equilibrium](pre-equilibrium.md) |
| **Hydrodynamics** | `MUSIC`, `CLVisc`, `Brick`, `Gubser`, `HydroFromFile` | [Hydrodynamics](hydro.md) |
| **Hard process** | `PythiaGun`, `PGun`, `EPGun`, `epemGun` | [Hard Process](hard-process.md) |
| **Jet energy loss** | `Matter`, `Lbt`, `Martini`, `AdSCFT` (+ iMATTER) | [Jet Energy Loss](energy-loss.md) |
| **Hadronization** | `colorless` (Lund), `colored`, `hybrid` (recombination) | [Hadronization](hadronization.md) |
| **Particlization** | `iSS` (Cooper–Frye sampler) | [Particlization](particlization.md) |
| **Afterburner** | `SMASH` | [Afterburner](afterburner.md) |
| **Source terms** | Causal / hadronic liquefier | [Liquefier](liquefier.md) |

Many modules are **external packages** wrapped by a thin X-SCAPE adapter (e.g.
`MusicWrapper`, `SmashWrapper`, `iSpectraSamplerWrapper`). The wrapper lives in
`src/`; the physics engine is fetched into `external_packages/` by the
corresponding `get_*.sh` script and enabled with a CMake `-DUSE_*=ON` flag. See
[Installation](../guides/installation.md).

!!! info "References"
    Every module page cites its primary papers inline. A consolidated
    bibliography with arXiv links is on the [References](../references.md) page.
