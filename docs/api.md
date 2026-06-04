# API Reference

The complete C++ API is generated from the in-source documentation comments by
**Doxygen** and published alongside this site.

<div class="grid cards" markdown>

- :material-code-tags: **[Browse the Doxygen API Reference →](/X-SCAPE/api/cpp/index.html)**

    Every framework and physics-module class, with inheritance and
    collaboration diagrams, member documentation, and source cross-references.

</div>

The API reference covers all of X-SCAPE's own code under `src/`:

| Area | Namespace / key classes |
|---|---|
| **Orchestration** | `Jetscape::JetScape`, `JetScapeTask`, `JetScapeModuleBase`, `JetScapeModuleFactory` |
| **Signals** | `JetScapeSignalManager` |
| **Clocks** | `ClockBase`, `MainClock`, `ModuleClock`, `MilneClock`, `TimeModule` |
| **Bulk coordination** | `BulkDynamicsManager`, `BulkMediaBase`, `BulkMediaInfo`, `LiquefierBase`, `QueryHistory` |
| **Stage base classes** | `InitialState`, `PreequilibriumDynamics`, `FluidDynamics`, `HardProcess`, `JetEnergyLoss`, `JetEnergyLossModule`, `Hadronization`, `HadronizationModule`, `SoftParticlization`, `Afterburner` |
| **Event record** | `Parton`, `Hadron`, `PartonShower`, `JetScapeParticleBase`, `JetScapeEvent` |
| **Medium** | `FluidCellInfo`, `FluidEvolutionHistory`, `SurfaceCellInfo`, `SurfaceFinder` |
| **Physics modules** | `Matter`, `LBT`, `Martini`, `AdSCFT`, `iMATTER`, `MpiMusic`, `CLVisc`, `Brick`, `GubserHydro`, `TrentoInitial`, `MCGlauberWrapper`, `IPGlasmaWrapper`, `PythiaGun`, `EPGun`, `HybridHadronization`, `ColorlessHadronization`, `iSpectraSamplerWrapper`, `SmashWrapper`, … |
| **I/O** | `JetScapeWriter*`, `JetScapeReader*`, `RootBulkWriter`, `FastRootBulkWriter` |
| **Configuration / logging** | `JetScapeXML`, `JetScapeLogger` |

!!! note "How the reference is built"
    The CI workflow runs `doxygen docs/Doxyfile` and places the HTML under
    `site/api/cpp/`, so it is served at
    [`/X-SCAPE/api/cpp/`](/X-SCAPE/api/cpp/index.html) (this MkDocs page itself
    is served at `/X-SCAPE/api/`). To regenerate it locally:

    ```bash
    pip install -r docs/requirements.txt   # for the MkDocs site
    sudo apt-get install doxygen graphviz  # for the API reference
    mkdocs build --site-dir site           # builds this site into ./site
    doxygen docs/Doxyfile                  # adds the API under ./site/api/cpp
    mkdocs serve                           # preview the site at localhost:8000
    ```

    The legacy whole-tree config `JetScapeDoxy.conf` in the repository root is
    still available for a standalone (non-Pages) Doxygen build.
