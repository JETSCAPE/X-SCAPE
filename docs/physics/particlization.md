# Particlization & Sampling

## Physics motivation

Hydrodynamics describes the bulk as a **continuous fluid**, but detectors see
**discrete hadrons**, and the late hadronic stage is a particle cascade. At some
point the fluid description must be converted into particles — **particlization**.
This is done on a **freeze-out hypersurface** $\Sigma$, the surface of constant
temperature (or energy density) below which the fluid picture is abandoned.

The standard prescription is the **Cooper–Frye formula**, which gives the
invariant momentum spectrum of hadron species $i$ as a flux of the local thermal
distribution through the surface:

\[
E\frac{dN_i}{d^3p} \;=\; \frac{g_i}{(2\pi)^3}\!\int_\Sigma
   f_i\!\big(u_\mu p^\mu; T, \mu_i\big)\, p^\mu\, d\Sigma_\mu ,
\]

with viscous ($\delta f$) corrections from the shear $\pi^{\mu\nu}$ and bulk
$\Pi$ stresses so that particlization is consistent with the viscous fluid it
came from. To feed a hadronic **afterburner** (which needs individual particles,
not smooth spectra), the Cooper–Frye spectrum is **Monte-Carlo sampled** into a
concrete list of hadrons with positions and momenta on $\Sigma$, conserving the
relevant quantum numbers.

`SoftParticlization` (`src/framework/SoftParticlization.{h,cc}`) is the stage
base class; it obtains the freeze-out surface from hydro via the
`GetHydroHyperSurface` signal (see [Signals](../framework/signals.md)) and
produces the sampled hadron list for the [afterburner](afterburner.md).

## Modules

### iSS — iSpectraSampler (Cooper–Frye sampler)

**XML:** `<SoftParticlization><iSS>…` · **class:** `iSpectraSamplerWrapper`
(`src/hadronization/iSpectraSamplerWrapper.{h,cc}`, engine in
`external_packages/iSS`, fetched by `get_iSS.sh`, built with `-DUSE_ISS=ON`)

iSS performs the **Monte-Carlo Cooper–Frye sampling**: given the freeze-out
hypersurface produced by [MUSIC](hydro.md) (via the Cornelius surface finder),
it samples discrete hadrons species-by-species, including viscous $\delta f$
corrections, and outputs a particle list suitable for the
[SMASH afterburner](afterburner.md). It is the particlization step of the
established iEBE-VISHNU / MUSIC hybrid workflow.

*Reference:* iSS sampler, as used in iEBE-VISHNU — Shen, Qiu, Song, Bernhard,
Bass, Heinz, [arXiv:1409.8164](https://arxiv.org/abs/1409.8164).

### Freeze-out surface finding

The hypersurface itself is extracted during the hydro evolution. X-SCAPE uses
the **Cornelius** algorithm (`external_packages/Cornelius`,
`src/framework/SurfaceFinder.{h,cc}`) to construct the constant-temperature
surface element-by-element from the evolving energy-density field, producing the
$d\Sigma_\mu$ and the flow/viscous data each surface cell carries
(`SurfaceCellInfo`).

*Reference:* Cornelius surface finder — Huovinen & Petersen,
[arXiv:1206.3371](https://arxiv.org/abs/1206.3371).

!!! warning "Resonance decays and the afterburner"
    When iSS feeds the SMASH afterburner, iSS should **not** perform resonance
    decays itself — SMASH handles them as part of the hadronic cascade. Doing
    them twice double-counts. See the
    [afterburner](afterburner.md) page and the project README.

## Two routes to hadrons

It is worth keeping the two hadron-production mechanisms distinct:

| Source | Mechanism | Module |
|---|---|---|
| **jet partons** (hard sector) | fragmentation / recombination | [Hadronization](hadronization.md) |
| **bulk fluid** (soft sector) | Cooper–Frye particlization | **iSS** (this page) |

Both feed the final particle record and, optionally, the
[afterburner](afterburner.md).
