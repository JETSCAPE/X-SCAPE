# Examples

X-SCAPE ships two kinds of examples:

- **`examples/`** — standard driver programs and analysis helpers.
- **`examples/custom_examples/`** — programs that build the task tree
  *explicitly in C++* instead of from XML. These demonstrate the framework API
  directly and are the place to learn the new X-SCAPE clock / Bulk Dynamics
  Manager features.

## Standard drivers (`examples/`)

| File | Purpose |
|---|---|
| `runJetscape.cc` | the main XML-driven executable (`runJetscape`) |
| `readerTest.cc` | read Ascii output, reconstruct the shower graph, run a FastJet closure test |
| `FinalStatePartons.cc`, `FinalStateHadrons.cc` | minimal final-state readers |
| `PythiaIsrTest.cc` | iMATTER ISR driver |
| `graph_fancy.py` | export a shower to a Graphviz/Gephi graph |
| `convert_binary_to_ASCII_output.py` | convert binary dumps to Ascii |
| `get_hydroSample_*.sh` | download sample hydro profiles for `HydroFromFile` |

## Custom (programmatic) examples (`examples/custom_examples/`)

These require disabling automatic task-list construction:

```xml
<enableAutomaticTaskListDetermination> false </enableAutomaticTaskListDetermination>
```

and editing `CMakeLists.txt` to build the chosen executable.

| File | Demonstrates |
|---|---|
| `PythiaBrickTest.cc` | PYTHIA + Brick + MATTER, built by hand; clock objects |
| `PythiaBDMTest.cc` | the **Bulk Dynamics Manager** in time-stepped mode |
| `PythiaNoBDMTest.cc` | the same chain *without* the BDM, for comparison |
| `MUSICMainClockTest.cc` | MUSIC driven by the **main clock** |
| `TwoStagesHydro.cc`, `TwoStagesHydroFromFile.cc` | two-stage (pre-eq → hydro) evolution |
| `PythiaIsrMUSIC.cc` | iMATTER ISR coupled to MUSIC |
| `SMASHInitialCondition.cc`, `SMASHNucleusTest.cc` | SMASH initial state |
| `LBT_brickTest.cc`, `brickTest.cc` | energy loss in a static brick |
| `CLViscTest.cc`, `MUSICTest.cc` | GPU / MPI hydro drivers |
| `hydroJetTest.cc`, `hydroFileTest.cc` | jets on a stored hydro medium |

## The programmatic API in one screen

The custom examples follow this pattern (condensed from
`PythiaBrickTest.cc`) — it mirrors exactly the
[task graph](../framework/architecture.md) you would otherwise specify in XML:

```cpp
using namespace Jetscape;

// 1. Create the root task and point it at the XML defaults.
auto jetscape = make_shared<JetScape>();
jetscape->SetXMLMainFileName("../config/jetscape_main.xml");
jetscape->SetXMLUserFileName("../config/jetscape_user_test.xml");
jetscape->SetId("primary");

// 2. (X-SCAPE) Attach clocks for time-stepped, reversible evolution.
auto mClock = make_shared<MainClock>("SpaceTime", -0.1, 0.1, 0.1); // t0, t_end, dt
jetscape->AddMainClock(mClock);

// 3. Build the pipeline by adding modules as sub-tasks, in order.
auto pythiaGun = make_shared<PythiaGun>();
auto hydro     = make_shared<Brick>();
jetscape->Add(pythiaGun);
jetscape->Add(hydro);

// 4. Energy loss is a manager -> per-shower loss -> modules sub-tree.
auto jlossmanager = make_shared<JetEnergyLossManager>();
auto jloss        = make_shared<JetEnergyLoss>();
auto matter       = make_shared<Matter>();
jloss->Add(matter);                 // add LBT / Martini / AdSCFT here too
jlossmanager->Add(jloss);
jetscape->Add(jlossmanager);

// 5. Run: Init -> Exec (per event) -> Finish.
jetscape->Init();
jetscape->Exec();
jetscape->Finish();
```

The same five steps drive *every* run; the XML path just builds the tree for
you. See [Tasks, Modules & the Factory](../framework/tasks-modules.md) and
[The Dynamical Clock](../framework/clock.md) for the API used here, and the
[API Reference](../api.md) for the full class documentation.

## Unit tests

`examples/unittests/` holds the Google Test suite (built with
`-Dunittests=ON`). Run them from the build directory with `ctest`. They are a
good, compiling reference for how individual classes are used.
