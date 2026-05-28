# Running multiple X-SCAPE instances in parallel

This document describes how to run several `runJetscape` simulations
**concurrently from a single build (or install) tree**, each in its own working
directory. Previously this required creating a separate full build directory per
concurrent run, because every module looked up its data files and wrote its
output relative to the current working directory (CWD), so two runs in the same
directory collided.

Phase 1 of the install refactor makes the read-only data assets locatable from a
shared prefix and provides a launcher that gives each run an isolated working
directory. See `PlanInstall.md` for the full roadmap (Phase 2 adds a proper
`make install` with a `bin/ lib/ share/` layout).

## The model

* **Read-only assets are shared.** EOS tables, iSS tables, LBT tables,
  3dMCGlauber inputs/tables, the config XML, etc. live once in a shared "data
  directory" (for an in-build run this is simply the build directory).
* **Each run gets its own working directory.** All writable output and the few
  files a module rewrites in place stay in that per-run directory, so concurrent
  runs never collide.

## How assets are located

| Asset | Resolution (first match wins) |
|-------|-------------------------------|
| 3dMCGlauber `mcglauber.input` | XML `IS/MCGlauber/mcglauber_input_file` → `$XSCAPE_DATA_DIR/mcglauber.input` → `./mcglauber.input` |
| MUSIC EOS tables | env `HYDROPROGRAMPATH` (set to the data dir) → `./EOS` |
| LBT tables | XML `Eloss/Lbt/LBT_table_path` → env `LBT_TABLES_PATH` → `$XSCAPE_DATA_DIR/LBT-tables` → `./LBT-tables` |
| iSS tables / parameters | XML `SoftParticlization/iSS/iSS_table_path`, `iSS_input_file`, … (the launcher rewrites these to absolute paths) |
| MUSIC input file | XML `Hydro/MUSIC/MUSIC_input_file` (CWD-relative `music_input`; the launcher copies a writable per-run copy) |

`XSCAPE_DATA_DIR` (default `.`) is the single framework-wide knob introduced for
this. When it is unset, every fallback reduces to the historical
build-directory behavior, so existing workflows are unaffected.

### New optional XML knobs

Two documented-but-inactive knobs were added to `config/jetscape_main.xml`. Add
them with an explicit path to override; **do not leave them empty** (an empty XML
element is unsupported and will crash the parser):

```xml
<IS>
  <MCGlauber>
    <mcglauber_input_file>/abs/path/mcglauber.input</mcglauber_input_file>
    ...
  </MCGlauber>
</IS>

<Eloss>
  <Lbt>
    <LBT_table_path>/abs/path/LBT-tables</LBT_table_path>
    ...
  </Lbt>
</Eloss>
```

## The launcher: `examples/run_in_workdir.sh`

The launcher wires the model together for one run. For each invocation it:

1. creates the per-run working directory (`-w`);
2. writes a per-run copy of the main XML with the iSS paths rewritten to
   absolute locations under the data dir;
3. copies `music_input` into the run directory (MUSIC rewrites it in place, so
   it must be per-run);
4. symlinks the read-only directories the vendored code still reads CWD-relative
   (`tables`, `eps09`, `LHAPDF_Lib`, `nucleusConfigs`, `data_table`);
5. exports `XSCAPE_DATA_DIR`, `HYDROPROGRAMPATH`, and `LBT_TABLES_PATH`;
6. runs `runJetscape <user.xml> <per-run main.xml>` with the run directory as CWD.

### Usage

```
examples/run_in_workdir.sh -w RUN_DIR -u USER_XML [-d DATA_DIR] [-b RUNJETSCAPE] [-m MAIN_XML]

  -w RUN_DIR      Per-run working directory (created if missing). REQUIRED.
  -u USER_XML     User XML config. REQUIRED.
  -d DATA_DIR     Shared asset prefix. Default: $XSCAPE_DATA_DIR if set,
                  else the directory containing the runJetscape executable.
  -b RUNJETSCAPE  Path to runJetscape. Default: DATA_DIR/runJetscape.
  -m MAIN_XML     Main XML config. Default: DATA_DIR/config/jetscape_main.xml,
                  else DATA_DIR/../config/jetscape_main.xml (in-build layout).
```

### Single run

```bash
# Build tree at ./build, run in /tmp/run1
examples/run_in_workdir.sh \
  -d "$PWD/build" \
  -w /tmp/run1 \
  -u "$PWD/config/jetscape_user_3DGlauber_MUSIC_iSS_SMASH_test.xml"
```

### Many runs in parallel (one build tree)

```bash
BUILD="$PWD/build"
USER_XML="$PWD/config/jetscape_user_3DGlauber_MUSIC_iSS_SMASH_test.xml"
for i in $(seq 1 8); do
  examples/run_in_workdir.sh -d "$BUILD" -w "/tmp/runs/run_$i" -u "$USER_XML" \
    > "/tmp/runs/run_$i.log" 2>&1 &
done
wait
```

Each `/tmp/runs/run_$i` ends up with its own `event.dat`, `events_summary.dat`,
`strings_event_*.dat`, `surface_eps_*`, freeze-out/eccentricity files, etc. The
shared assets in `$BUILD` are never modified. Give each run a distinct
`<Random><seed>` and/or `<outputFilename>` in the user XML if you want
statistically independent events.

## Notes & limitations

* **`config/` location for in-build runs.** The build dir does *not* contain
  `config/` (it lives at the source root). The launcher falls back to
  `DATA_DIR/../config/jetscape_main.xml`, which works when the build directory is
  inside the source tree (e.g. `X-SCAPE/build`). For an out-of-source build, pass
  `-m <source>/config/jetscape_main.xml` explicitly. Phase 2 installs `config/`
  under `share/xscape` so the default resolves automatically.
* **3dMCGlauber `tables/` must stay read-only.** Keep `cache_tables 1` (the
  shipped default) so the valence-quark sample tables are only read, never
  rewritten/regenerated. `cache_tables 0` is incompatible with shared assets and
  parallel runs.
* **SMASH afterburner** needs its own data files; if your config enables SMASH,
  make sure those are available from the run directory as well.
* **Never run two instances in the same directory.** MUSIC issues
  `rm surface*.dat` and writes fixed-name outputs; isolation comes from each run
  owning its CWD.
