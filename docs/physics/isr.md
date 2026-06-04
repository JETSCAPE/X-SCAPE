# Initial-State Radiation (iMATTER)

## Physics motivation

Before two partons scatter at high momentum transfer, each is part of an
incoming hadron and carries some space-like virtuality. As it is pulled toward
the hard vertex it radiates — **initial-state radiation (ISR)** — building up a
*space-like* (time-reversed) parton shower. In the vacuum this is standard
DGLAP backward evolution, handled inside PYTHIA. In a nuclear collision the
incoming partons traverse cold (and, in the overlap region, forming) nuclear
matter, so the ISR can itself be **medium-modified**, and the *longitudinal*
position of each initial splitting matters for where the hard scattering sits
relative to the developing fireball.

For small systems and e–A collisions, X-SCAPE needs ISR that:

- conserves four-momentum **exactly** against the hard process and the bulk
  (so the soft–hard energy budget closes — see the
  [Bulk Dynamics Manager](../framework/bulk-dynamics.md));
- tracks the **space-time** location of each space-like split, not just its
  momentum, so the shower can be embedded in a dynamical initial state.

**iMATTER** ("initial-state MATTER") is the module that does this. It is the
time-reversed analogue of the [MATTER](energy-loss.md) final-state shower: the
same virtuality-ordered formalism, run **backward** in time from the hard vertex
toward the incoming hadrons, with the longitudinal location of the initial
splits included. This is precisely the kind of physics the X-SCAPE
[reversible clock](../framework/clock.md) was designed to support — the main
clock runs backward through the ISR and forward through the final-state
evolution.

*References:* the X-SCAPE small-system framework,
[arXiv:2308.02650](https://arxiv.org/abs/2308.02650); soft–hard framework with
exact four-momentum conservation,
[arXiv:2407.17443](https://arxiv.org/abs/2407.17443). The underlying MATTER
formalism is in [arXiv:1301.5323](https://arxiv.org/abs/1301.5323).

## Modules and how they fit together

| Component | Class / file | Role |
|---|---|---|
| ISR-aware hard process | `PythiaIsrGun` (`src/initialstate/PythiaIsrGun.{h,cc}`) | runs the PYTHIA hard scattering and exposes the incoming partons for ISR |
| ISR manager | `IsrManager` (`src/framework/IsrManager.{h,cc}`) | the stage manager that drives the space-like shower (analogue of `JetEnergyLossManager`) |
| ISR shower generator | `IsrShowerPSG` (`src/framework/IsrShowerPSG.{h,cc}`) | generates the backward shower (analogue of the final-state `PartonShowerGenerator`) |
| ISR jet object | `IsrJet` (`src/framework/IsrJet.{h,cc}`) | the space-like parton/shower bookkeeping |
| ISR rotation | `ISRRotation` (`src/jet/ISRRotation.{h,cc}`) | kinematic alignment of the ISR shower |
| iMATTER engine | `iMATTER` (`src/jet/iMATTER.{h,cc}`) | the virtuality-ordered backward-evolution physics |
| Test driver | `InitialStateRadiationTest` (`src/initialstate/`) | standalone test harness |
| ISR output | `JetScapeWriterIsrStream` (`src/framework/`) | writes the ISR shower record |

The ISR machinery mirrors the final-state energy-loss machinery one-to-one — an
`IsrManager`/`IsrShowerPSG` pair playing the roles of
`JetEnergyLossManager`/`PartonShowerGenerator` — which is what makes "the same
shower, run backwards" a literal statement in the code.

## Running it

iMATTER is exercised through its own executable rather than the generic
`runJetscape` path:

```bash
# requires 3DGlauber support compiled in and $PYTHIA8 set
./PythiaIsrTest          # uses config/jetscape_user_iMATTERMCGlauber.xml
```

This runs iMATTER together with the [3DGlauber](initial-state.md) dynamical
initial state. To couple it to hydro (MUSIC) as well, see
`config/jetscape_user_iMATTERMCGlauberMUSIC.xml` and the
[3DGlauber wiki page](https://github.com/JETSCAPE/X-SCAPE/wiki/3DGlauber,-MUSIC,-iSS-and-Initial-State-Radiation).

!!! note "Prerequisites"
    Set `$PYTHIA8` to your PYTHIA 8 install (`pythia8-config --prefix`) and
    build with `-DUSE_3DGlauber=ON`. The JETSCAPE Docker container sets these
    up for you.
