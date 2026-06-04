# Output, Writers & Readers

## Writers

Output is handled by **writer modules** that subscribe to the task tree and are
handed each task during the per-event `WriteTasks()` sweep (see the
[execution model](architecture.md)). Any number of writers can be active at
once; each is toggled by a tag in the [configuration](configuration.md).

| Writer | XML tag | Format | Source |
|---|---|---|---|
| Ascii | `JetScapeWriterAscii` | human-readable text: full parton-shower history + final particles | `JetScapeWriterStream` |
| Ascii (gzip) | `JetScapeWriterAsciiGZ` | gzip-compressed Ascii | `JetScapeWriterStream` |
| HepMC | `JetScapeWriterHepMC` | HepMC3 event record (conforms to [arXiv:1912.08005](https://arxiv.org/abs/1912.08005) App. A) | `JetScapeWriterHepMC` |
| ROOT HepMC | `JetScapeWriterRootHepMC` | HepMC in a ROOT file (requires `-DUSE_ROOT=ON`) | `JetScapeWriterRootHepMC` |
| Final-state partons | `JetScapeWriterFinalStatePartonsAscii` | only the final partons | `JetScapeWriterFinalStateStream` |
| Final-state hadrons | `JetScapeWriterFinalStateHadronsAscii` | only the final hadrons | `JetScapeWriterFinalStateStream` |
| Qₙ vector | `JetScapeWriterQnVectorAscii` | flow Qₙ vectors | `JetScapeWriterQnCalculator` |

### Bulk writers (medium dumps)

For studies that need the **medium itself** — e.g. training surrogate
(neural-network) hydro models or visualizing the fireball — X-SCAPE provides
writers that dump the 4-D evolution rather than the particle record. These live
in `src/root/` and require `-DUSE_ROOT=ON`:

| Writer | XML tag | Notes |
|---|---|---|
| `RootBulkWriter` | `RootBulkWriter` | the (τ, x, y, η) bulk fields written to a ROOT tree |
| `FastRootBulkWriter` | `FastRootBulkWriter` | a faster, streamlined variant (see `README_BulkFast.md`) |

The bulk-writer design and its trade-offs are documented in
`BulkWriterImprovements.md` and `README_BulkFast.md` in the repository root.

### Stream filtering

`JetScapeWriterStreamFilter` lets a writer subset what it emits (e.g. only
final-state particles, or only a particular status code) without changing the
producing modules.

## The event record content

The Ascii/HepMC writers emit, per event:

- an **event header** (`JetScapeEventHeader`) — cross section, weight, impact
  parameter, event-plane angle, and any module-contributed header fields,
  gathered by the `CollectHeaders()` sweep;
- the **parton-shower graph** — every parton and every splitting vertex, so the
  full DAG of the shower can be reconstructed offline;
- the **final hadrons** from hadronization and the afterburner.

## Readers

To analyse Ascii output, X-SCAPE ships reader classes
(`src/reader/JetScapeReader*.{h,cc}`) that rebuild the in-memory objects from a
file:

- `JetScapeReader` — reconstructs the full `PartonShower` graph (using the
  bundled **GTL** graph library), so you can walk the shower, run a
  depth-first search, or feed the final partons to FastJet.
- `JetScapeReaderFinalStateHadrons` — a lightweight reader for the
  final-state-hadron format.

The example `examples/readerTest.cc` reads a shower, walks the graph, and runs a
simple FastJet jet-finding "closure" check; `examples/graph_fancy.py` exports
the shower to a graph format viewable in Graphviz/Gephi.

```bash
./build/readerTest          # read showers, DFS, FastJet closure test
python examples/convert_binary_to_ASCII_output.py   # convert binary dumps
```

## Choosing an output format

| You want… | Use |
|---|---|
| quick inspection / debugging | `JetScapeWriterAscii` |
| compact storage of many events | `JetScapeWriterAsciiGZ` |
| interoperability with HEP analysis tools | `JetScapeWriterHepMC` |
| only final particles for an analysis | final-state partons/hadrons writers |
| the medium fields (ML training, viz) | `RootBulkWriter` / `FastRootBulkWriter` |
| flow-observable Qₙ vectors | `JetScapeWriterQnVectorAscii` |
