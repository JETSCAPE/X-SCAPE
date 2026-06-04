# Initial State

## Physics motivation

The initial state sets the stage for everything downstream: it specifies *where*
and *how much* energy (and, at lower energies, net baryon number) is deposited
in the transverse plane at the moment of impact. Two features make this
non-trivial and physically rich:

- **Event-by-event fluctuations.** Nucleon positions, and the sub-nucleonic
  structure of each nucleon, fluctuate. These lumps in the initial geometry are
  what hydrodynamics converts into the anisotropic flow harmonics (v₂, v₃, …)
  measured in data. Getting the *fluctuation spectrum* right is essential.
- **Longitudinal / baseline structure.** At top RHIC and LHC energies a
  boost-invariant 2-D profile is a good approximation. At beam-energy-scan
  energies, and in small systems, the finite nuclear crossing time and net-baryon
  stopping demand a genuinely **3-D, dynamical** initial state.

X-SCAPE offers several initial-state models spanning these regimes, plus simple
"guns" for controlled tests.

`InitialState` (`src/framework/InitialState.{h,cc}`) is the stage base class; it
stores the sampled energy/entropy density on a grid and exposes it to
pre-equilibrium and hydro.

## Modules

### TRENTo — parametric initial condition

**XML:** `<IS><Trento>…</Trento></IS>` · **class:** `TrentoInitial`
(`src/initialstate/TrentoInitial.{h,cc}`, engine in `external_packages/trento`)

TRENTo ("**T**hickness-**R**educed **E**vent-by-event **N**uclear **To**pology)
is a *parametric* model that does not commit to a specific particle-production
mechanism. It builds the participant thickness functions $T_A$, $T_B$ of the two
nuclei and combines them through the **generalized mean**

\[
s(x,y) \;\propto\; \left(\frac{T_A^{\,p} + T_B^{\,p}}{2}\right)^{1/p},
\]

where the single parameter $p$ interpolates between known mechanisms
($p\to 0$ reproduces an EKRT/IP-Glasma-like geometric mean, $p=1$ a wounded-nucleon
arithmetic mean, $p=-1$ a harmonic mean). Because $p$ is continuous, TRENTo is
the workhorse for **Bayesian global analyses** of the QGP: it lets the data
choose the production scheme. Multiplicity fluctuations are tunable via a
Gamma-distributed nucleon fluctuation parameter, and nucleon width / structure
parameters control the granularity.

*Reference:* Moreland, Bernhard, Bass,
[arXiv:1412.4708](https://arxiv.org/abs/1412.4708).

### 3DGlauber / MC-Glauber — dynamical 3-D initial state

**XML:** `<IS>` with `MCGlauber` · **class:** `MCGlauberWrapper`,
`MCGlauberGenStringWrapper` (`src/initialstate/`, engine in
`external_packages/3dMCGlauber`, fetched by `get_3dglauber.sh`)

A fully **three-dimensional, dynamical** Glauber model. Rather than depositing
energy instantaneously on a 2-D plane, it forms **strings** between colliding
participants (down to the valence-quark level) and lets them **decelerate**,
depositing energy and **net baryon number** over a finite formation time and a
range of space-time rapidity. This is the appropriate initial state for:

- the **RHIC Beam Energy Scan**, where baryon stopping and finite-energy effects
  dominate;
- **small systems**, where the longitudinal structure and sub-nucleonic hotspots
  matter — it provides the dynamical initial state used by X-SCAPE's soft–hard
  framework (with the hard scattering emanating from colliding hotspots).

It feeds MUSIC directly and integrates with the
[Bulk Dynamics Manager](../framework/bulk-dynamics.md).

*References:* Shen & Schenke,
[arXiv:1710.00881](https://arxiv.org/abs/1710.00881);
[arXiv:1711.10544](https://arxiv.org/abs/1711.10544).

### IP-Glasma — saturation-based color fields

**XML:** `<IS>` with `IPGlasma` · **class:** `IPGlasmaWrapper`
(`src/initialstate/`, external `ipglasma`)

IP-Glasma combines the **IP-Sat** impact-parameter-dependent saturation model
for the nucleon's gluon distribution with a **classical Yang–Mills (glasma)**
evolution of the produced color fields. It captures both nucleonic *and*
sub-nucleonic color-charge fluctuations from first-principles-inspired CGC
dynamics, and naturally provides the early-time (pre-equilibrium) gluon fields.
Its field initialization is FFT-based.

*Reference:* Schenke, Tribedy, Venugopalan,
[arXiv:1202.6646](https://arxiv.org/abs/1202.6646).

### SMASH initial condition — hadronic initial state

**XML:** run via the `SMASHInitialCondition` executable · **class:**
`SMASHInitialStateWrapper`, `SMASHNucleusWrapper`
(`src/initialstate/`, see also [Afterburner](afterburner.md))

At **low and intermediate collision energies** the system does not start as a
deconfined plasma; the appropriate initial condition is a **hadronic transport**
description. The SMASH initial-condition module runs the SMASH cascade up to a
hypersurface of constant proper time and hands the resulting energy-momentum and
baryon distributions to MUSIC. It uses the same engine as the
[SMASH afterburner](afterburner.md).

*Reference:* SMASH — Weil et al.,
[arXiv:1606.06642](https://arxiv.org/abs/1606.06642).

### Reading a pre-computed initial state

| Module | XML / class | Use |
|---|---|---|
| Initial profile from file | `InitialFromFile` | read a tabulated energy/entropy density profile |
| Binary-collision list from file | `NcollListFromFile` | seed hard-scattering positions from an external $N_{\rm coll}$ list |

These let you decouple a (possibly expensive or external) initial-state
calculation from the rest of the run.

## Particle "guns" (controlled tests)

For testing energy loss, hadronization, or the framework itself without a full
collision, X-SCAPE provides parton/particle guns. These sit in the *hard
process* slot but define the initial particle content; see
[Hard Process](hard-process.md):

| Gun | XML | Produces |
|---|---|---|
| `PGun` | `<PGun>` | a single parton of chosen flavor/energy |
| `PythiaGun` | `<PythiaGun>` | a full PYTHIA hard-scattering event |
| `PythiaIsrGun` | `<PythiaIsrGun>` | PYTHIA hard process with initial-state radiation hooks (drives [iMATTER](isr.md)) |
| `EPGun` | `<EPGun>` | electron–proton DIS / photoproduction events |
| `epemGun` | `<epemGun>` | e⁺e⁻ annihilation events |

## Choosing an initial-state model

| You are studying… | Use |
|---|---|
| flow / Bayesian QGP extraction at RHIC-top / LHC | **TRENTo** |
| beam-energy scan, net-baryon dynamics, 3-D structure | **3DGlauber** |
| CGC / saturation physics, sub-nucleonic fluctuations | **IP-Glasma** |
| low-energy / hadronic initial state | **SMASH-IC** |
| controlled energy-loss / framework tests | a **gun** + `Brick`/`Gubser` hydro |
