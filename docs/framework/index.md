# The X-SCAPE Framework

X-SCAPE is built as a **modular, task-based envelope**: the physics lives in
interchangeable *modules*, and the framework's job is to assemble them into a
pipeline, run them, let them exchange information, and write the result. Nothing
in the core depends on the details of any particular physics model — modules are
discovered by name at run time and wired together according to an XML
configuration.

This section documents the machinery that makes that possible.

<div class="grid cards" markdown>

- **[Architecture & Execution Model](architecture.md)** — the task graph, the
  per-event execution loop, and how X-SCAPE adds concurrent, time-stepped
  evolution on top of JETSCAPE's per-event model.

- **[Tasks, Modules & the Factory](tasks-modules.md)** — `JetScapeTask`,
  `JetScapeModuleBase`, the self-registration factory, and the base classes a
  new physics module derives from.

- **[Signals & Data Flow](signals.md)** — the `sigslot` publish/subscribe layer
  through which modules request partons, medium information, and surfaces from
  one another.

- **[The Dynamical Clock](clock.md)** — `MainClock`, `ModuleClock`,
  `MilneClock`, and the time-stepped execution that can advance *and rewind*.

- **[Bulk Dynamics Manager](bulk-dynamics.md)** — the coordinator that lets
  several bulk generators run concurrently and feeds a unified medium to the
  hard sector.

- **[XML Configuration](configuration.md)** — the two-file (main + user) XML
  system, how parameters are resolved, and how modules read them.

- **[Output, Writers & Readers](io.md)** — Ascii, gzip, HepMC and ROOT writers,
  the bulk writers, and the reader/`PartonShower` graph utilities.

</div>

## Key source locations

| Concern | Files (`src/framework/`) |
|---|---|
| Top-level orchestration | `JetScape.{h,cc}` |
| Task / module base | `JetScapeTask.{h,cc}`, `JetScapeModuleBase.{h,cc}` |
| Signals | `JetScapeSignalManager.{h,cc}`, `sigslot.h` |
| Clocks | `ClockBase`, `MainClock`, `ModuleClock`, `MilneClock`, `TimeModule` |
| Bulk coordination | `BulkDynamicsManager.{h,cc}`, `BulkMediaBase`, `BulkMediaInfo` |
| Stage base classes | `InitialState`, `PreequilibriumDynamics`, `FluidDynamics`, `HardProcess`, `JetEnergyLoss`, `Hadronization`, `SoftParticlization`, `Afterburner`, `LiquefierBase` |
| Event record | `JetScapeParticles`, `JetClass`, `PartonShower`, `JetScapeEvent` |
| Configuration | `JetScapeXML.{h,cc}` |
| Logging | `JetScapeLogger.{h,cc}` |
| I/O | `JetScapeWriter*`, `JetScapeReader*` |
