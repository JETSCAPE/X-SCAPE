# Hard Process

## Physics motivation

**Hard probes** — high transverse-momentum partons, heavy quarks, and
electroweak bosons — originate in the rare, large-momentum-transfer scatterings
that occur in the very first instant of the collision, *before* the medium
forms. Their production cross sections are calculable in perturbative QCD and are
well calibrated in p–p collisions, which makes them ideal **tomographic probes**:
whatever modification they show in A–A relative to p–p is attributable to the
QGP they traversed.

The hard-process stage produces these initial high-virtuality partons (and the
underlying event), which are then handed to the [energy-loss](energy-loss.md)
stage as the seeds of in-medium parton showers.

`HardProcess` (`src/framework/HardProcess.{h,cc}`) is the stage base class; it
stores the produced hard partons and exposes them through the
`GetHardPartonList` signal (see [Signals](../framework/signals.md)).

## Modules

### PythiaGun — PYTHIA 8 hard scattering

**XML:** `<Hard><PythiaGun>…</PythiaGun></Hard>` · **class:** `PythiaGun`
(`src/initialstate/PythiaGun.{h,cc}`)

Wraps **PYTHIA 8** to generate hard QCD (or other) scatterings. Key parameters:

| Tag | Meaning |
|---|---|
| `<pTHatMin>` / `<pTHatMax>` | the hard-scattering $\hat p_T$ window |
| `<eCM>` | center-of-mass energy (GeV) |
| `<LinesToRead>` | arbitrary PYTHIA settings (process selection, PDFs, …) |

PYTHIA produces the hard partons and (optionally) the full vacuum shower /
underlying event. With nuclear PDFs (via **LHAPDF**, see the
[README](https://github.com/JETSCAPE/X-SCAPE#running-jetscapex-scape-with-lhapdf))
it also models initial-state nuclear modification of the parton distributions.
`PythiaIsrGun` is the ISR-aware variant that drives [iMATTER](isr.md).

*Reference:* PYTHIA 8 — Sjöstrand et al.,
[arXiv:1410.3012](https://arxiv.org/abs/1410.3012).

### PGun — single-parton gun

**XML:** `<Hard><PGun>…` · **class:** `PGun`
(`src/initialstate/PGun.{h,cc}`)

Produces a **single parton** of chosen flavor, energy, and direction. The
controlled-input counterpart to PythiaGun: combined with a
[`Brick`](hydro.md) medium it isolates the energy-loss physics of one parton
with no underlying event. The standard quenching benchmarks use it.

### EPGun — electron–proton gun (DIS / photoproduction)

**XML:** `<Hard><EPGun>…` · **class:** `EPGun`
(`src/initialstate/EPGun.{h,cc}`)

X-SCAPE's entry point to **electron–ion physics**. It generates electron–proton
events in two regimes selected by the virtuality $Q^2$ of the exchanged photon:

- **DIS** (deep-inelastic scattering) for $Q^2 > 1\ \mathrm{GeV}^2$ — valid down
  to HERMES energies;
- **photoproduction** for $Q^2 < 1\ \mathrm{GeV}^2$ — works across all HERA
  energies.

Example configs: `config/jetscape_user_DIS.xml`,
`config/jetscape_user_photoproduction.xml`.

### epemGun — e⁺e⁻ annihilation

**XML:** `<Hard><epemGun>…` · **class:** `epemGun`
(`src/initialstate/epemGun.{h,cc}`)

Generates $e^+e^-$ annihilation events — the cleanest QCD environment, with no
hadronic initial state — useful for validating the vacuum parton shower and
hadronization against the wealth of LEP data.

## Choosing a hard-process module

| You want… | Use |
|---|---|
| realistic jets in p–p / A–A | **PythiaGun** |
| a single, controlled parton (energy-loss tests) | **PGun** |
| ISR + small-system soft–hard | **PythiaIsrGun** (→ [iMATTER](isr.md)) |
| electron–ion (DIS / photoproduction) | **EPGun** |
| e⁺e⁻ validation | **epemGun** |
