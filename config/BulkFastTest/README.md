# BulkFastTest — example configs for FastRootBulkWriter

Example/test XML configs for the hydro-only fast ROOT dump (`FastRootBulkWriter`).
See [../../README_BulkFast.md](../../README_BulkFast.md) for the full usage guide and
[../../BulkWriterImprovements.md](../../BulkWriterImprovements.md) for the design notes.

## Files

| file | what it is |
|---|---|
| `OO_one_event.xml` | baseline O+O 1-event config using the **legacy** `RootBulkWriter` (for comparison) |
| `OO_one_event_fast.xml` | hydro-only, `FastRootBulkWriter` in **`native`** mode (writes MUSIC's own grid) |
| `OO_one_event_fastgrid.xml` | hydro-only, `FastRootBulkWriter` in **`grid`** mode on the legacy 65×65×33 grid |
| `validate.py` | compares the three ROOT outputs (numpy + uproot) |

The two `*_fast*` configs differ from the baseline only by:
- removing the `<Eloss>` block (hydro-only — nothing may query the now-empty medium),
- adding `<dump_hydro_only>1` inside `<Hydro><MUSIC>`,
- replacing `<RootBulkWriter>` with `<FastRootBulkWriter>`.

## How to run

`runJetscape` resolves MUSIC's input/EOS paths relative to the current directory, so run
from `build_gpu/` and pass the config by path:

```bash
cd build_gpu
./runJetscape ../config/BulkFastTest/OO_one_event_fast.xml       # native -> OO_test_fast.root
./runJetscape ../config/BulkFastTest/OO_one_event_fastgrid.xml   # grid   -> OO_test_fastgrid.root
./runJetscape ../config/BulkFastTest/OO_one_event.xml            # legacy -> OO_test.root
```

## Validate

Use a Python with numpy + uproot (the system `python3` here has neither; use conda
`fno_env`):

```bash
cd build_gpu
conda activate fno_env
python ../config/BulkFastTest/validate.py
```

## Verified result (O+O, 1 event)

- **`grid` mode vs legacy `RootBulkWriter`: bit-for-bit identical** — same shape (both
  `ntau_freezeout=27`, 15,057,900 floats on 65×65×33) and `np.array_equal == True`
  (`max |abs diff| = 0`). Files ≈ 16.83 MB.
- **`native` mode:** MUSIC's full grid (here 100×100×60, ntau=134 at `dtau≈0.02`); file
  ≈ 93 MB; no NaN/inf; peak energy 7.48 vs legacy 7.36 (native keeps full resolution).
  Thin with `<tau_stride>` or MUSIC's `output_evolution_every_N_*`.

> ⚠️ Requires `config/jetscape_main.xml` to contain default entries for `<dump_hydro_only>`
> and `<FastRootBulkWriter>` (already added). Without them X-SCAPE aborts with
> "tag is unrecognized".
