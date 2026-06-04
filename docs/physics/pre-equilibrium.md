# Pre-equilibrium

## Physics motivation

Hydrodynamics assumes the matter is close enough to local thermal equilibrium
that an equation of state and transport coefficients make sense. But the
initial state is produced **far from equilibrium**, and it takes a finite time
$\tau_0 \sim 0.2$–$1.0$ fm/c for the system to relax. What happens in that gap
matters:

- The choice of hydro start time $\tau_0$ is otherwise arbitrary; a
  pre-equilibrium stage lets results become **insensitive to $\tau_0$** by
  modeling the approach to equilibrium explicitly.
- The pre-equilibrium dynamics builds up the early **transverse flow** and shapes
  the **shear-stress tensor** $\pi^{\mu\nu}$ that hydro inherits as its initial
  condition. Starting hydro from zero flow / zero shear is unphysical.
- It bridges a **CGC/glasma** or **free-streaming** initial state to a fluid by a
  matching (Landau matching) of the energy-momentum tensor.

`PreequilibriumDynamics` (`src/framework/PreequilibriumDynamics.{h,cc}`) is the
stage base class: it takes the initial-state $T^{\mu\nu}$, evolves it to
$\tau_0$, and provides the energy density, flow, and viscous tensor that
initialize [hydrodynamics](hydro.md).

## Modules

### Free-streaming (freestream-milne)

**XML:** `<Preequilibrium><FreestreamMilne>…` · **class:**
`FreestreamMilneWrapper` (`src/preequilibrium/FreestreamMilneWrapper.{h,cc}`)

The simplest physically motivated pre-equilibrium model: treat the deposited
matter as a collection of **non-interacting (free-streaming) massless partons**
that propagate ballistically from the initial proper time to $\tau_0$, then
perform a **Landau matching** of the streamed $T^{\mu\nu}$ to viscous
hydrodynamic variables (energy density, flow velocity, shear, bulk). Free
streaming is the opposite limit to ideal hydro (infinite vs. zero mean free
path); used over a short interval it provides a controlled, $\tau_0$-stabilizing
bridge and seeds the initial transverse flow and $\pi^{\mu\nu}$. The
`freestream-milne` engine works in Milne coordinates, matching the boost-invariant
hydro grid.

### Glasma

**XML:** `<Preequilibrium><Glasma>…` · **class:** `Glasma`
(`src/preequilibrium/Glasma.{h,cc}`)

When the initial state is [IP-Glasma](initial-state.md), the natural
pre-equilibrium stage is the **classical Yang–Mills evolution of the glasma
color fields** themselves. This module carries the early-time CGC field dynamics
forward to the hydro starting surface, consistently with a saturation-based
initial state.

### NullPreDynamics

**XML:** `<Preequilibrium><NullPreDynamics>…` · **class:** `NullPreDynamics`
(`src/preequilibrium/NullPreDynamics.{h,cc}`)

A **pass-through** placeholder: no pre-equilibrium evolution. Use it when the
initial-state model already hands a profile directly to hydro (e.g. TRENTo →
MUSIC with a chosen $\tau_0$), or for tests where you want to isolate the effect
of pre-equilibrium by switching it off. It satisfies the
`PreequilibriumDynamics` interface while doing nothing.

## Choosing a pre-equilibrium model

| Initial state | Natural pre-equilibrium |
|---|---|
| TRENTo (parametric) | **free-streaming** (or **none** with a chosen $\tau_0$) |
| IP-Glasma (CGC) | **Glasma** classical-field evolution |
| 3DGlauber (dynamical) | dynamical sourcing into hydro (often **none** as a separate stage) |
| controlled tests | **NullPreDynamics** |

## References

- Free streaming + Landau matching as a pre-equilibrium stage:
  Liu, Shen, Bass and collaborators (free-streaming approach);
  the `freestream-milne` implementation. See the
  [References](../references.md) page.
- Glasma / IP-Glasma classical Yang–Mills:
  Schenke, Tribedy, Venugopalan,
  [arXiv:1202.6646](https://arxiv.org/abs/1202.6646).
