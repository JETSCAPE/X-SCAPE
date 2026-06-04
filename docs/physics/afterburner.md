# Hadronic Afterburner

## Physics motivation

Particlization produces a hadron gas on the freeze-out surface, but those
hadrons are **not** the final state. As the fireball continues to expand the
hadrons keep **scattering** (elastically and resonantly) and **decaying** until
the gas is so dilute that interactions cease — **kinetic freeze-out**. This late
hadronic phase:

- reshapes **soft observables** — spectra, particle ratios, and especially the
  differential flow of identified particles, which is sensitive to the
  species-dependent hadronic cross sections;
- handles **resonance regeneration and decay** consistently;
- matters most for **low-energy** collisions and for the **dilute late stage** of
  high-energy ones, where a fluid description is no longer justified.

A microscopic **hadronic transport** model is the right tool. X-SCAPE uses
**SMASH**.

`Afterburner` (`src/framework/Afterburner.{h,cc}`) is the stage base class; it
receives the sampled hadron list from [particlization](particlization.md) and
runs the cascade.

## Module: SMASH

**XML:** `<Afterburner><SMASH>…` · **class:** `SmashWrapper`
(`src/afterburner/SmashWrapper.{h,cc}`, engine in `external_packages/smash`,
fetched by `get_smash.sh`, built with `-DUSE_SMASH=ON`)

**SMASH** ("**S**imulating **M**any **A**ccelerated **S**trongly-interacting
**H**adrons") is a relativistic hadronic transport approach. It propagates the
particlized hadrons and solves the relativistic Boltzmann equation for the hadron
gas, including binary scatterings, resonance formation and decay, and string
excitation at higher energies, until interactions cease. In X-SCAPE it serves
two roles:

- as the **afterburner** following [iSS](particlization.md) particlization of a
  MUSIC hydro evolution — the standard route to soft observables;
- as a **hadronic initial condition** at low/intermediate energies (the
  [SMASH initial-condition](initial-state.md) module shares the same engine).

A full afterburner run therefore needs hydro + sampler + SMASH together:

```bash
cmake -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_FREESTREAM=ON -DUSE_SMASH=ON ..
make
./SMASHTest
```

*Reference:* SMASH — Weil et al.,
[arXiv:1606.06642](https://arxiv.org/abs/1606.06642);
[arXiv:1808.06832](https://arxiv.org/abs/1808.06832).

!!! warning "Switch off resonance decays in the sampler"
    SMASH performs resonance decays as part of the cascade. The
    [iSS](particlization.md) sampler must therefore have its own resonance
    decays **disabled**, or they are double-counted. See the project README.

## Where it fits

```mermaid
flowchart LR
    HY["Hydro (MUSIC)"] -->|"freeze-out surface"| iSS["iSS<br/>Cooper-Frye sampling"]
    iSS -->|"hadron list"| SMASH["SMASH<br/>hadronic cascade"]
    SMASH -->|"final hadrons"| OUT["Writers"]
```

`CascadeTest` (`src/afterburner/CascadeTest.{h,cc}`) is a standalone test of the
afterburner stage.
