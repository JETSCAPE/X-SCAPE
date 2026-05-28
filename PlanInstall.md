# X-SCAPE: genuine install + parallel runs without per-run build dirs

## Context

Today X-SCAPE is run directly from inside a `build/` directory. To run several
instances at once, users currently create several full build directories,
because every module resolves its files relative to the current working
directory (CWD = the build dir):

- **Read-only assets** (executable, libs, EOS tables, iSS tables, LBT tables,
  3dMCGlauber tables, config XML, input files) are materialized into the build
  dir at *configure* time via `file(COPY)` / `configure_file` / symlink, and
  looked up CWD-relative at runtime.
- **Writable outputs** use fixed names in CWD (`surface.dat`,
  `events_summary.dat`, `OSCAR.DAT`, `evolution_all_xyeta.dat`, …), and MUSIC
  even **rewrites its input file in place** (`update_music_input_parameter`).

So two concurrent runs in the same dir collide. The fix is to (1) relocate
read-only assets into a shared install prefix discoverable via one env var, and
(2) let each run execute in its own *empty, lightweight* working directory for
writable output. Then "multiple instances at once" = multiple cheap run dirs
against one install — no rebuilds, no duplicate build trees.

**Decisions (confirmed with user):** staged delivery (MVP first, full install
second); make the *whole* framework relocatable; single data-root env var
`XSCAPE_DATA_DIR` (extending the existing `HYDROPROGRAMPATH` pattern), set by a
launcher.

## What already works (no code change)

- **MUSIC EOS** already honors env var `HYDROPROGRAMPATH`
  (`external_packages/music/src/eos_base.cpp` `get_hydro_env_path()`, default
  `"."`). Map `HYDROPROGRAMPATH=$XSCAPE_DATA_DIR`.
- **iSS** is fully XML-redirectable: `iSS_table_path`,
  `iSS_particle_table_path`, `iSS_input_file`, `iSS_working_path`
  (`src/hadronization/iSpectraSamplerWrapper.cc` → `external_packages/iSS/src/iSS.cpp`).
- **3dMCGlauber sample caches** (`tables/{proton,neutron}_valence_quark_samples*.dat`)
  are READ-only with the shipped default `cache_tables 1`
  (`external_packages/3dMCGlauber/src/Parameters.cpp` `get_cached_tabels()`); the
  `rm`/regen path (`Glauber.cpp`, `Nucleus.cpp` `system("./Metropolis.e …")`)
  only fires when `cache_tables 0` or files are missing. Keep `cache_tables 1`
  and ship the caches → runs never write `tables/`.
- **Primary JETSCAPE outputs** (ascii/hepmc/root) are already XML-configurable
  via `<outputFilename>` (`src/framework/JetScape.cc`).

## True blockers (need code) for full relocation

- **`mcglauber.input` hard-coded** — `src/initialstate/MCGlauberWrapper.cc:45-46`
  (`new MCGlb::EventGenerator("mcglauber.input", …)`).
- **`LBT-tables/` hard-coded** — `src/jet/LBT.cc:95-96` and the `ifstream("LBT-tables/…")`
  reads in `read_tables()` (ratedata, ratedata-HQ, dNg_over_dt_*, distB/F).
- **3dMCGlauber CWD-relative dirs** `tables/`, `eps09/`, `LHAPDF_Lib` (vendored
  sources). Handled by launcher symlink (preferred — avoids editing vendored
  code); deeper code change optional/deferred.
- **MUSIC input in-place rewrite** — `MusicWrapper.cc:105,224` + defn `727-778`.
  Handled by giving each run a *local copy* of `music_input` (see Phase 1.3),
  which also preserves the MUSIC↔iSS shared-file coupling.

---

## Phase 1 — MVP: parallel runs against the existing build dir

Goal: run unlimited concurrent instances against ONE existing build dir, each in
its own empty CWD. No bin/lib/share refactor yet. The build dir acts as
`$XSCAPE_DATA_DIR`.

1. **Data-root helper (framework).** Add a tiny helper that returns the data
   root: `getenv("XSCAPE_DATA_DIR")`, else a baked default (Phase 2), else `"."`.
   Put it where both initial-state and jet modules can use it (e.g. a small
   function in `src/framework`). Reuse for the two code blockers below.

2. **3dMCGlauber input path configurable (CODE).**
   `src/initialstate/MCGlauberWrapper.cc:45-46` — resolve the input file via a
   new XML element `{"IS","MCGlauber","mcglauber_input_file"}` (fallback
   `"mcglauber.input"`, optionally prefixed by the data-root helper); pass to the
   `EventGenerator` ctor. Add `<mcglauber_input_file>` under `IS/MCGlauber` in
   `config/jetscape_main.xml`.

3. **LBT-tables dir configurable (CODE).** `src/jet/LBT.cc` — in `InitTask()`
   resolve a base dir once: env `LBT_TABLES_PATH` → else data-root helper +
   `/LBT-tables` → else `"LBT-tables"`. Prefix `setParameter(dir+"/LBT.input")`
   (L95-96) and every `ifstream("LBT-tables/…")` in `read_tables()`. Add
   `<LBT_table_path>` under `Eloss/Lbt` as an alternative knob. (Leave the
   standalone-only `../hydroProfile/...` path untouched.)

4. **MUSIC input → per-run local copy (LAUNCHER, no code).** The launcher copies
   the shared `music_input` into the run CWD and points `<MUSIC_input_file>` at
   that local copy, with `<iSS_working_path>` = the run CWD. Since MUSIC rewrites
   whatever `MUSIC_input_file` points at, and iSS symlinks that same path into
   `iSS_working_path/music_input`, the per-run local copy keeps the shared copy
   read-only while preserving the coupling.

5. **Launcher script (NEW FILE)** e.g. `examples/run_in_workdir.sh`. Per run:
   - `mkdir -p $RUNDIR && cd $RUNDIR`
   - `export XSCAPE_DATA_DIR=<build-or-prefix>`; `export HYDROPROGRAMPATH=$XSCAPE_DATA_DIR`;
     `export LBT_TABLES_PATH=$XSCAPE_DATA_DIR/LBT-tables`
   - symlink read-only CWD-relative dirs: `tables`, `eps09`, `LHAPDF_Lib`,
     `nucleusConfigs`, `data_table`
   - `cp $XSCAPE_DATA_DIR/music_input ./music_input` (writable per-run — NOT a symlink)
   - run `runJetscape <abs user.xml> <abs main.xml>`; outputs land in `$RUNDIR`.

After Phase 1: parallel instances work from one build dir; outputs isolated;
shared assets read-only.

---

## Phase 2 — Full clean `bin/ lib/ share/` install

1. **Install prefix + GNUInstallDirs** (`CMakeLists.txt`). Remove/relax the block
   (~L229-234) that FORCES `CMAKE_INSTALL_PREFIX` to the build dir. Add
   `include(GNUInstallDirs)`; define `XSCAPE_DATADIR = ${CMAKE_INSTALL_DATADIR}/xscape`.

2. **Install executables via TARGETS** (`CMakeLists.txt`).
   `install(TARGETS runJetscape … RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})`;
   gate optional MUSIC/iSS/Isr test executables behind their existing `if()`
   blocks. Build-tree `EXECUTABLE_OUTPUT_PATH` (L545) may stay for dev.

3. **Install libraries via TARGETS — replace the DIRECTORY copy hack**
   (`CMakeLists.txt:307-323`). Use `install(TARGETS …)` for the export list at
   L732-751 (`JetScape JetScapeThird GTL libtrento Cornelius` + conditional
   `hydroFromFile ipglasma_lib 3dMCGlb music iSS clviscwrapper`) plus
   `JetScapeReader`, into `${CMAKE_INSTALL_LIBDIR}`. Fix sub-package CMakes that
   install into `${CMAKE_HOME_DIRECTORY}`:
   `external_packages/iSS/src/CMakeLists.txt`,
   `external_packages/3dMCGlauber/src/CMakeLists.txt` (incl. `Metropolis.e`),
   `external_packages/music/src/CMakeLists.txt` → bindir/libdir.

4. **Install config + data to share/xscape** (alongside the existing configure-time
   `file(COPY)` at L646-725; both coexist — copies for in-build dev, `install()`
   for the relocatable tree):
   - `config/` → `${XSCAPE_DATADIR}/config`
   - `external_packages/music/EOS/` → `.../EOS`; `examples/test_music_files/music_input` → `.../music_input`
   - `external_packages/iSS/iSS_tables/` → `.../iSS_tables`; `iSS_parameters.dat` → `.../`
   - `external_packages/3dMCGlauber/{tables,eps09,LHAPDF_Lib}/` → `.../…`;
     `external_packages/3dMCGlauber/input` → `.../mcglauber.input` (rename)
   - `external_packages/trento/nucleusConfigs/`, `src/initialstate/data_table/`,
     `external_packages/LBT-tables/` (guard `EXISTS`) → `.../…`

5. **Bake the data root + env override** (single mechanism). `configure_file` a
   generated header (`XscapeInstallPaths.h.in`) with
   `#define XSCAPE_INSTALL_DATADIR "${CMAKE_INSTALL_FULL_DATADIR}/xscape"`. The
   data-root helper (Phase 1.1) prefers `getenv("XSCAPE_DATA_DIR")`, else this
   baked default, else `"."`. `examples/runJetscape.cc:58-59`: fall back to
   `XSCAPE_INSTALL_DATADIR "/config/…"` when `../config/…` is absent (keep argv
   override). Installed `config/jetscape_main.xml` uses absolute
   `$XSCAPE_DATA_DIR/...` paths for iSS / MUSIC input / mcglauber input.

6. **RPATH so installed `runJetscape` finds `<prefix>/lib`** (`CMakeLists.txt`,
   before targets):
   ```
   set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
   set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
   if(APPLE)
     set(CMAKE_MACOSX_RPATH ON)
     set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
   else()
     set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
   endif()
   ```
   Reconcile the 3dMCGlauber sub-CMake `INSTALL_RPATH` (references
   `${CMAKE_HOME_DIRECTORY}/LHAPDF_Lib/lib/`) so LHAPDF resolves at the prefix.
   Ensure all shared libs install to one `${CMAKE_INSTALL_LIBDIR}`.

7. **Install the launcher** to `${CMAKE_INSTALL_BINDIR}`, defaulting
   `XSCAPE_DATA_DIR` to the baked `${CMAKE_INSTALL_FULL_DATADIR}/xscape`.

---

## Critical files

- `/Users/du8478/JetScape/X-SCAPE/CMakeLists.txt` (install layout, RPATH, data install)
- `/Users/du8478/JetScape/X-SCAPE/src/CMakeLists.txt` (lib install, data_table)
- `/Users/du8478/JetScape/X-SCAPE/src/initialstate/MCGlauberWrapper.cc` (mcglauber.input path)
- `/Users/du8478/JetScape/X-SCAPE/src/jet/LBT.cc` (LBT-tables path)
- `/Users/du8478/JetScape/X-SCAPE/src/hydro/MusicWrapper.cc` (in-place rewrite; per-run copy strategy)
- `/Users/du8478/JetScape/X-SCAPE/examples/runJetscape.cc` (config XML fallback)
- `/Users/du8478/JetScape/X-SCAPE/config/jetscape_main.xml` (new path knobs)
- new: `examples/run_in_workdir.sh`, `XscapeInstallPaths.h.in`

## Verification

1. **Install tree:** `cmake -S . -B /tmp/xbuild -DUSE_MUSIC=ON -DUSE_ISS=ON
   -DUSE_3DGlauber=ON -DCMAKE_INSTALL_PREFIX=/tmp/xprefix && cmake --build
   /tmp/xbuild -j && cmake --install /tmp/xbuild`. Confirm
   `/tmp/xprefix/{bin/runJetscape, lib/libJetScape.*, share/xscape/{config,EOS,
   iSS_tables,LBT-tables,tables,eps09,LHAPDF_Lib,nucleusConfigs,data_table,
   mcglauber.input,music_input}}`. `otool -L` / `readelf -d` shows rpath →
   `<prefix>/lib`; run with empty `DYLD_/LD_LIBRARY_PATH`.
2. **Single run, arbitrary empty CWD:** launcher + a 3DGlauber+MUSIC+iSS user XML
   (e.g. `config/jetscape_user_3DGlauber_MUSIC_iSS_*.xml`). Outputs only in the
   run dir; assert shared `music_input` and `tables/*samples*.dat` mtimes
   unchanged.
3. **Two concurrent runs:** two run dirs, distinct seeds / `<outputFilename>`.
   Both finish; outputs differ; NO shared-asset mtime changed; diff the dirs.

## Risks / flags

- **LBT** unusable from a shared tree without Phase 1.3; static table reads done
  once.
- **`cache_tables 0`** is incompatible with shared read-only assets / parallel
  runs (it `rm`s + regenerates `tables/` and would race). Document: keep
  `cache_tables 1`.
- **3dMCGlauber `Metropolis.e`** uses a literal `./` (`Nucleus.cpp`); only fires
  if sample files missing — shipped caches + symlinked `tables/` avoid it.
- **MUSIC `music.cpp`** does `rm surface*.dat` + fixed-name outputs — safe only
  because each run owns its CWD; never run two instances in one dir.
- **Vendored `external_packages/*`**: prefer launcher symlinks + framework-side
  wrappers over editing vendored sources, to ease future re-syncs.
- Sub-package CMakes install into `${CMAKE_HOME_DIRECTORY}` and trento creates a
  `bin/`; apply the Phase 2.3 redirects consistently or the layout will be
  inconsistent.
