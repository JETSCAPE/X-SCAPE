# Hadronization

## Physics motivation

Quarks and gluons are never observed in isolation: **confinement** forces them
to combine into color-singlet hadrons before they reach the detector. The
energy-loss stage produces a shower of (colored) partons; hadronization is the
non-perturbative step that turns them into the pions, kaons, protons, and other
hadrons that make up a jet and the bulk.

In a heavy-ion collision two distinct mechanisms operate:

- **Fragmentation** (string / cluster) dominates at high $p_T$, where a fast
  parton fragments much as it would in vacuum. This is the Lund-string picture
  of PYTHIA.
- **Recombination / coalescence** becomes important at intermediate $p_T$ in the
  presence of a medium: a shower parton can pick up a thermal partner from the
  QGP and **coalesce** into a hadron. Recombination naturally explains the
  observed **baryon-to-meson enhancement** and the constituent-quark scaling of
  flow at intermediate $p_T$ — features fragmentation alone cannot reproduce.

X-SCAPE provides modules for both, and a **hybrid** that combines them.

`Hadronization` (`src/framework/Hadronization.{h,cc}`) and its
`HadronizationManager` are the stage base/manager; modules derive from
`HadronizationModule<Derived>` (CRTP) and implement `DoHadronization(...)`.
The choice is made with `<JetHadronization><name>…</name>`.

## Modules

### Colorless hadronization (Lund string)

**XML:** `<JetHadronization><name>colorless</name>` · **class:**
`ColorlessHadronization` (`src/hadronization/ColorlessHadronization.{h,cc}`)

Hands the final shower partons to **PYTHIA's Lund string fragmentation**,
ignoring the detailed color flow of the in-medium shower (hence "colorless" — the
strings are assigned in a simplified way). This is the standard, robust choice
for **high-$p_T$ jet** studies, where fragmentation dominates and the medium's
effect on color flow is sub-leading.

*Reference:* Lund string model in PYTHIA 8 — Sjöstrand et al.,
[arXiv:1410.3012](https://arxiv.org/abs/1410.3012).

### Colored hadronization

**XML:** `<JetHadronization><name>colored</name>` · **classes:**
`ColoredHadronization`, `iColoredHadronization`
(`src/hadronization/`)

Retains the **color-tag information** of the in-medium shower when forming
strings, so that the color connections built up during energy loss (including
recoils) are respected at fragmentation. `iColoredHadronization` is the variant
used in the time-stepped / ISR-aware (`i`-prefixed) workflow.

### Hybrid hadronization (recombination + string)

**XML:** `<JetHadronization><name>hybrid</name>` · **class:**
`HybridHadronization` (`src/hadronization/HybridHadronization.{h,cc}`)

The most complete option: it combines **recombination** at intermediate $p_T$
with **string fragmentation** of the leftover partons at high $p_T$. Each shower
parton may either coalesce with a near-by (shower or thermal) parton into a
hadron, or, failing that, be connected into a Lund string. Thermal partners are
drawn from the local medium by the **thermal parton sampler**
(`ThermPtnSampler`, `src/hadronization/ThermPtnSampler.{h,cc}`), which samples a
thermal distribution consistent with the hydro freeze-out. Hybrid hadronization
is what lets X-SCAPE describe the full $p_T$ range and baryon/meson chemistry in
one consistent picture.

*References:* recombination/coalescence — Han, Fries, Cao
([arXiv:1601.01583](https://arxiv.org/abs/1601.01583)); the JETSCAPE hybrid
hadronization module (see the [References](../references.md) page).

## Choosing a hadronization module

| You are studying… | Use |
|---|---|
| high-$p_T$ jet fragmentation | **colorless** (Lund) |
| color-flow effects from in-medium showers | **colored** |
| full $p_T$ range, baryon/meson ratios, coalescence | **hybrid** |

The bulk fluid is converted to hadrons by a *different* mechanism — Cooper–Frye
sampling on the freeze-out surface — documented under
[Particlization](particlization.md).
