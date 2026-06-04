# Running X-SCAPE

## The main executable

A standard JETSCAPE-mode run uses `runJetscape`, built into your build
directory, and takes a [user XML file](../framework/configuration.md):

```bash
cd build
./runJetscape ../config/jetscape_user.xml
```

This reads the main XML (`config/jetscape_main.xml`, all defaults) overlaid with
your user XML (modules to run + overrides), generates the requested number of
events, and writes them in the enabled [output formats](../framework/io.md).

## X-SCAPE / specialized executables

Some X-SCAPE capabilities and external-package combinations have their own
executables (built when the corresponding `-DUSE_*` option is on). They each
default to a matching `config/jetscape_user_*.xml`:

| Executable | What it runs | Default config |
|---|---|---|
| `runJetscape` | generic per-event pipeline | argument (e.g. `jetscape_user.xml`) |
| `MUSICTest` | MUSIC hydro chain (run under MPI) | `jetscape_user_MUSIC.xml` |
| `PythiaIsrTest` | iMATTER ISR + 3DGlauber | `jetscape_user_iMATTERMCGlauber.xml` |
| `SMASHInitialCondition` | SMASH hadronic initial state | `jetscape_user_SMASHInitialCondition.xml` |
| `SMASHTest` | hydro + iSS + SMASH afterburner | — |
| `brickTest`, `LBT_brickTest` | energy loss in a static brick | — |
| `readerTest` | read Ascii output, FastJet closure | — |

MUSIC-based runs use MPI:

```bash
mpirun -np 1 ./MUSICTest
```

## Configuring a run

Everything is in the XML. See [XML Configuration](../framework/configuration.md)
for the full scheme. The essentials:

- set `<nEvents>`;
- enable one or more [writers](../framework/io.md) (`<JetScapeWriterAscii>on…`);
- add a module block per stage you want
  (`<Hard>`, `<Eloss>`, `<Hydro>`, `<JetHadronization>`, …);
- override any default parameter by repeating its tag inside the module block.

Start from a [published tune](https://github.com/JETSCAPE/Default-tunes) or one
of the bundled `config/jetscape_user_*.xml` files rather than from scratch.

## Output and analysis

Analysis is left to the user (it is outside the framework's scope). The
container ships ROOT, Python, and FastJet. To inspect Ascii output:

```bash
./readerTest                                  # walk showers, run FastJet
python ../examples/graph_fancy.py             # export shower to Graphviz/Gephi
python ../examples/convert_binary_to_ASCII_output.py
```

See [Output, Writers & Readers](../framework/io.md).

## Parallel runs

Execution is serial per event, so production scales by running **many
independent events as separate jobs**. A launcher script for batched parallel
runs is documented in
[`README_Launcher.md`](https://github.com/JETSCAPE/X-SCAPE/blob/main/README_Launcher.md)
and [`INSTALL.md` §9](https://github.com/JETSCAPE/X-SCAPE/blob/main/INSTALL.md).

## Nuclear PDFs (LHAPDF)

Download a set and reference it from `PythiaGun`:

```bash
cd external_packages && ./get_lhapdf.sh JAM20-SIDIS_PDF_proton_nlo
```

```xml
<LinesToRead>
  PDF:useHard = on
  PDF:pHardSet = LHAPDF6:JAM20-SIDIS_PDF_proton_nlo
</LinesToRead>
```

Or run the ready-made `config/jetscape_user_nPDF_test.xml`.
