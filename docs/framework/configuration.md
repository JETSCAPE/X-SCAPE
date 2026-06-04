# XML Configuration

An X-SCAPE run is fully specified by **two XML files**. This two-file design
separates the (large, rarely edited) catalogue of every parameter from the
(small, run-specific) choice of what to actually do.

| File | You edit it? | Contains |
|---|---|---|
| **Main** — `config/jetscape_main.xml` | No | the default value of *every* parameter of *every* module |
| **User** — e.g. `config/jetscape_user.xml` | Yes | which modules to run, in what order, and which defaults to override |

At startup `JetScapeXML` (`src/framework/JetScapeXML.{h,cc}`, a singleton built
on the bundled **tinyxml2**) loads the main file, then overlays the user file:
any tag present in the user file replaces the main-file default; anything absent
falls back to the main default. Modules then read their settings through the
typed helpers on `JetScapeModuleBase`:

```cpp
double Q0  = GetXMLElementDouble({"Eloss", "Matter", "Q0"});
int in_vac = GetXMLElementInt   ({"Eloss", "Matter", "in_vac"});
string had = GetXMLElementText  ({"JetHadronization", "name"});
```

The path argument is the chain of nested tag names; pass `isRequired = false`
to make a missing element non-fatal.

## A minimal user file

`config/jetscape_user.xml` — a proton–proton jet run with vacuum MATTER and
colorless (Lund-string) hadronization:

```xml
<?xml version="1.0"?>
<jetscape>

  <nEvents> 2 </nEvents>
  <JetScapeWriterAscii> on </JetScapeWriterAscii>

  <!-- Hard process: PYTHIA hard scattering -->
  <Hard>
    <PythiaGun>
      <pTHatMin>50</pTHatMin>
      <pTHatMax>70</pTHatMax>
      <eCM>5020</eCM>
    </PythiaGun>
  </Hard>

  <!-- Jet energy loss: MATTER in vacuum -->
  <Eloss>
    <Matter>
      <in_vac> 1 </in_vac>
    </Matter>
  </Eloss>

  <!-- Hadronization -->
  <JetHadronization>
    <name>colorless</name>
  </JetHadronization>

</jetscape>
```

Run it with:

```bash
./runJetscape ../config/jetscape_user.xml
```

## Top-level / global settings

These live directly under `<jetscape>` (defaults shown from the main file):

| Tag | Meaning |
|---|---|
| `<nEvents>` | number of events to generate |
| `<setReuseHydro>` / `<nReuseHydro>` | reuse one hydro event for *N* hard events (the bulk is expensive; jets are cheap) |
| `<Random><seed>` | RNG seed; `0` means draw a fresh seed |
| `<outputFilename>` | base name for output files |
| `<JetScapeWriterAscii>` etc. | enable/disable each [writer](io.md) |
| `<nEvents_printout>` | progress-printout cadence |

## Stage blocks

Modules are grouped under stage tags that mirror the
[task graph](architecture.md). The **presence** of a module block activates it;
its execution order follows the canonical pipeline:

| Stage tag | Holds | Page |
|---|---|---|
| `<IS>` | initial-state module (`<Trento>`, `<initial_profile>`, …) | [Initial State](../physics/initial-state.md) |
| `<Preequilibrium>` | pre-equilibrium module | [Pre-equilibrium](../physics/pre-equilibrium.md) |
| `<Hydro>` | hydro module (`<MUSIC>`, `<Brick>`, `<CLVisc>`, …) | [Hydrodynamics](../physics/hydro.md) |
| `<Hard>` | hard process (`<PythiaGun>`, `<PGun>`, `<EPGun>`, …) | [Hard Process](../physics/hard-process.md) |
| `<Eloss>` | one or more energy-loss modules (`<Matter>`, `<Lbt>`, `<Martini>`, `<AdSCFT>`) | [Jet Energy Loss](../physics/energy-loss.md) |
| `<JetHadronization>` | hadronization choice (`name = colorless / colored / hybrid`) | [Hadronization](../physics/hadronization.md) |
| `<SoftParticlization>` | Cooper–Frye sampler (`<iSS>`) | [Particlization](../physics/particlization.md) |
| `<Afterburner>` | hadronic transport (`<SMASH>`) | [Afterburner](../physics/afterburner.md) |

Multiple energy-loss modules inside `<Eloss>` are applied as a **multistage**
shower (e.g. MATTER at high virtuality handing off to LBT at low virtuality);
see the [energy-loss](../physics/energy-loss.md) page.

## `LinesToRead` — passing raw settings to a backend

Some external packages (notably PYTHIA) take their own configuration strings.
These are passed verbatim inside a `<LinesToRead>` block, e.g. to load an
LHAPDF set:

```xml
<LinesToRead>
  PDF:useHard = on
  PDF:pHardSet = LHAPDF6:JAM20-SIDIS_PDF_proton_nlo
</LinesToRead>
```

## Tunes

Published parameter sets ("tunes") are distributed as ready-to-run user XML
files in the [Default-tunes](https://github.com/JETSCAPE/Default-tunes)
repository — start from a tune rather than from scratch when reproducing a
result. The bundled `config/jetscape_user_*.xml` files cover the common module
combinations (MUSIC, CLVisc, 3DGlauber, SMASH, DIS/photoproduction, nPDF, …).
