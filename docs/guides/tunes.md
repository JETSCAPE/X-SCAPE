# Published Tunes

## What a "tune" is

A physical prediction from X-SCAPE depends on dozens of parameters — the
switching virtuality $Q_0$, the $\hat q$ normalization, η/s, the freeze-out
temperature, the TRENTo geometry parameters, and many more. A **tune** is a
*complete, frozen parameter set* — distributed as a ready-to-run user XML file —
that reproduces the configuration used in a specific publication, typically the
result of a global/Bayesian calibration against experimental data.

Tunes matter because they make published results **reproducible** and give you a
**physically sensible starting point**: rather than guessing parameter values,
begin from the tune closest to your system and energy and override only what you
need.

## The Default-tunes repository

The JETSCAPE Collaboration maintains its published tunes in a dedicated
repository:

> **[github.com/JETSCAPE/Default-tunes](https://github.com/JETSCAPE/Default-tunes)**

It contains the user XML files (and any auxiliary tables) needed to reproduce
the results of the corresponding papers, organized by analysis. See the
repository's
[README](https://github.com/JETSCAPE/Default-tunes/blob/main/README.md) for the
current list of tunes, the publication each corresponds to, and which X-SCAPE
modules and external packages each one requires.

## How to run a tune

A tune is just a user XML file, so it runs through the standard
[`runJetscape`](running.md) path. Clone the tunes repository alongside X-SCAPE
and point the executable at the desired tune file:

```bash
# 1. Get the tunes (anywhere on disk)
git clone https://github.com/JETSCAPE/Default-tunes.git

# 2. Build X-SCAPE with the modules the tune requires
#    (check the tune's README — e.g. MUSIC + iSS + 3DGlauber + SMASH)
cmake -S . -B build -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_3DGlauber=ON
cmake --build build -j$(nproc)

# 3. Run the tune's user XML through runJetscape
cd build
./runJetscape /path/to/Default-tunes/<analysis>/<tune_user>.xml
```

The tune's user XML still overlays the same `config/jetscape_main.xml` defaults
described in [XML Configuration](../framework/configuration.md); it simply sets
the modules and the overridden parameter values that define the publication.

!!! tip "Match the build to the tune"
    A tune that uses MUSIC, iSS, 3DGlauber, and SMASH will fail at startup if
    those modules were not compiled in. Read the tune's README first and enable
    the matching `-DUSE_*` options (see [Installation](installation.md)).

!!! warning "Version compatibility"
    Tunes are calibrated against a specific framework/module version. Reproducing
    a published result exactly may require the X-SCAPE/JETSCAPE version noted in
    the tune's README; newer code can shift results. Check the tune's
    documentation for the intended version.

## Citing a tune

When you use a tune, cite **both** the JETSCAPE framework paper
([arXiv:1903.07706](https://arxiv.org/abs/1903.07706)) **and** the specific
publication the tune reproduces, as listed in the
[Default-tunes README](https://github.com/JETSCAPE/Default-tunes/blob/main/README.md).
See the [References](../references.md) page for the framework and module
citations.

## Related sample configurations

Beyond the published tunes, the X-SCAPE repository bundles many example user
files under `config/jetscape_user_*.xml` covering common module combinations
(MUSIC, CLVisc, 3DGlauber, SMASH, DIS/photoproduction, nPDF, two-stage hydro,
…). These are good starting templates even when they are not formal tunes.
Sample hydro profiles for `HydroFromFile` can be fetched with
`examples/get_hydroSample_*.sh`.
