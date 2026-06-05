# API Reference

The complete C++ API is generated from the in-source documentation comments by
**Doxygen** and rendered as native pages in this site by the
[mkdoxy](https://github.com/JakubAndrysek/MkDoxy) plugin. It is built
automatically by `mkdocs build` / `mkdocs serve` — no separate step.

<div class="grid cards" markdown>

- :material-format-list-bulleted-type: **[Class List →](api/annotated.md)**

    Every framework and physics-module class, with member documentation,
    inheritance, and source cross-references.

- :material-sitemap: **[Class Hierarchy →](api/hierarchy.md)**

    The inheritance tree across all stages and modules.

- :material-folder-multiple: **[Files →](api/files.md)**

    Browse by source file under `src/`.

- :material-tag-multiple: **[Namespaces →](api/namespaces.md)**

    `Jetscape` and related namespaces.

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
    The [mkdoxy](https://github.com/JakubAndrysek/MkDoxy) plugin runs Doxygen
    over `src/` during the normal site build and converts the result into the
    Material-themed pages linked above. Requirements:

    ```bash
    pip install -r docs/requirements.txt   # includes mkdoxy
    # Doxygen must be on PATH (e.g. `brew install doxygen graphviz`
    # or `sudo apt-get install doxygen graphviz`)
    mkdocs serve                           # API is generated + served live
    ```

    The legacy whole-tree config `JetScapeDoxy.conf` in the repository root is
    still available for a standalone (non-Pages) Doxygen build.
