# Tasks, Modules & the Factory

## `JetScapeTask` — the node type

Every runnable object in X-SCAPE derives from `JetScapeTask`
(`src/framework/JetScapeTask.{h,cc}`). A task is a node in the execution tree:
it owns a `vector<shared_ptr<JetScapeTask>>` of children and exposes the
recursive lifecycle described in the
[execution model](architecture.md#the-task-graph).

Key members:

```cpp
class JetScapeTask {
public:
  virtual void Init();                 // calls InitTask() then InitTasks()
  virtual void InitTask() {}           // override: this task's own init
  virtual void Exec();                 // calls ExecuteTask() then ExecuteTasks()
  virtual void ExecuteTask() {}        // override: this task's own work
  virtual void Clear();                // per-event reset
  virtual void Finish();               // shutdown
  virtual void WriteTask(weak_ptr<JetScapeWriter> w) {}

  virtual void Add(shared_ptr<JetScapeTask> t);   // attach a child
  const vector<shared_ptr<JetScapeTask>> GetTaskList() const;
  shared_ptr<JetScapeTask> GetTaskAt(int i);

  void SetActive(bool);   bool GetActive() const;     // skip without removing
  void SetId(string);     const string GetId() const; // human-readable id
};
```

The `active_exec` flag (`SetActive`/`GetActive`) lets the framework skip a task
without detaching it from the tree — used, for instance, to suppress a second
hydro run when re-reading the same event.

## `JetScapeModuleBase` — the module type

Physics modules derive from `JetScapeModuleBase`
(`src/framework/JetScapeModuleBase.{h,cc}`), which is a `JetScapeTask` plus three
mix-ins:

```cpp
class JetScapeModuleBase
    : public JetScapeTask,
      public sigslot::has_slots<sigslot::multi_threaded_local>, // signals
      public TimeModule {                                       // clocks
  ...
};
```

It provides the services every module needs:

- **XML access** — typed helpers that read from the merged main+user
  configuration:
  ```cpp
  std::string GetXMLElementText({"Eloss", "Matter", "name"});
  int         GetXMLElementInt ({"Eloss", "Matter", "in_vac"});
  double      GetXMLElementDouble({"Eloss", "Matter", "Q0"});
  tinyxml2::XMLElement* GetXMLElement({"Eloss", "Matter"});
  ```
  Each takes a path (a list of nested tag names) and an optional `isRequired`
  flag. See [XML Configuration](configuration.md).
- **Random numbers** — `GetMt19937Generator()` returns a seeded
  `std::mt19937` so module-level randomness is reproducible and centrally
  controlled.
- **Per-event bookkeeping** — `GetCurrentEvent()` / `IncrementCurrentEvent()`.
- **Time-stepped hooks** — `CalculateTime()`, `ExecTime()`, `IsTimeStepped()`,
  inherited from the clock-aware design (see [The Dynamical Clock](clock.md)).

## The self-registration factory

X-SCAPE never hard-codes the list of modules. Instead each module
**registers itself** under the name it will be referenced by in XML, using a
static `RegisterJetScapeModule` instance:

```cpp
// in Matter.cc
RegisterJetScapeModule<Matter> Matter::reg("Matter");
```

This inserts `("Matter", &createT<Matter>)` into a global map at static-init
time. When the framework parses the configuration and encounters
`<Matter>...</Matter>`, it calls:

```cpp
auto module = JetScapeModuleFactory::createInstance("Matter");
```

which returns a freshly constructed `shared_ptr<Matter>` — no framework header
needs to know `Matter` exists. **Adding a new module therefore requires no
changes to the framework**: write the class, register it, and reference it in
XML. The full set of registered names is listed under
[Physics Modules](../physics/index.md).

## Stage base classes

Between `JetScapeModuleBase` and a concrete module sits an **abstract stage
base class** that defines the interface for that stage and wires up the
appropriate signals. A new module derives from the stage base, not from
`JetScapeModuleBase` directly:

| Stage | Base class | A module implements… |
|---|---|---|
| Initial state | `InitialState` | `InitTask()`, `Exec()` to fill the entropy/energy density profile |
| Initial-state radiation | `JetEnergyLoss` (+ `IsrManager`) | space-like shower hooks |
| Pre-equilibrium | `PreequilibriumDynamics` | `InitializePreequilibrium()`, `EvolvePreequilibrium()` |
| Hydrodynamics | `FluidDynamics` | `InitializeHydro()`, `EvolveHydro()`, cell access |
| Hard process | `HardProcess` | `InitTask()`, `Exec()` to produce hard partons |
| Jet energy loss | `JetEnergyLossModule<T>` (CRTP) | `DoEnergyLoss(deltaT, time, Q2, in, out)` |
| Hadronization | `HadronizationModule<T>` (CRTP) | `DoHadronization(...)` |
| Particlization | `SoftParticlization` | sample the freeze-out surface |
| Afterburner | `Afterburner` | run hadronic transport on the sampled hadrons |
| Source terms | `LiquefierBase` | `add_hydro_sources(...)` |

The energy-loss and hadronization bases use the **Curiously Recurring Template
Pattern** (`JetEnergyLossModule<Derived>`, `HadronizationModule<Derived>`) so
the base can clone the derived type for per-shower instances without virtual
construction.

## The event record

Three data types carry the physics content of an event:

- **`Parton`** (`JetClass.{h,cc}`) — a four-momentum, PDG id, color tags,
  virtuality / form-time bookkeeping, and a position. Energy-loss modules pass
  `vector<Parton>&` in and out of `DoEnergyLoss`.
- **`PartonShower`** (`PartonShower.{h,cc}`) — a directed graph (built on the
  bundled **GTL** graph library) of `shared_ptr<Parton>` nodes connected by
  splitting vertices. This is the full shower history, kept alive for the event
  and cleared afterwards. The [readers](io.md) walk this graph.
- **`Hadron`** (`JetClass.{h,cc}`) and **`JetScapeParticles`** — the final-state
  particles produced by hadronization and the afterburner.

`Parton` and `Hadron` both derive from a common `JetScapeParticleBase`
(four-momentum + PDG identity), with PDG metadata provided by `pdgcode.{h,cc}`.

!!! tip "Writing a new module — checklist"
    1. Pick the stage and derive from its base class.
    2. Implement the stage's pure-virtual hooks plus `InitTask()`.
    3. Read parameters with the `GetXMLElement*` helpers.
    4. Register: `RegisterJetScapeModule<MyModule> MyModule::reg("MyName");`
    5. Add `<MyName>...</MyName>` defaults to `config/jetscape_main.xml`.
    6. (Optional) implement `CalculateTime()`/`ExecTime()` to run time-stepped.
    See [CONTRIBUTING](https://github.com/JETSCAPE/X-SCAPE/blob/main/CONTRIBUTING.md).
