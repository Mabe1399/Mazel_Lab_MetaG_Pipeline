# Pipeline smoke-test data

Small, synthetic, fully reproducible input data for validating pipeline changes on the
cluster before trusting them against real sequencing data. Intended for `dev`-branch
development — not a substitute for a real pilot run.

## What's here

- `reads/` — 2 samples (`sampleA`, `sampleB`), synthetic paired-end FASTQ, already
  "cleaned" (`{sample}_R{1,2}_cleaned.fastq.gz`). Simulated from 3 mock reference
  "species" genomes (80kb/60kb/40kb random sequence) at **different relative
  abundances per sample**, so binning tools have a real differential-coverage signal
  to work with — not just two identical samples.
- `generate_test_reads.py` — the generator that produced `reads/` (fixed seed —
  rerun it to regenerate from scratch rather than hand-editing the FASTQ files).
- `test_metadata.tsv` — minimal 2-row metadata file for the iTOL-annotation step
  (module 04), matching `config.test.yaml`'s `DATA_VARIABLES` column names.
- `config.test.yaml` — a ready-to-run config wired to this test data.

## What's deliberately NOT here

- **No host genome / Kraken2 database.** `config.test.yaml` sets `steps.QC: false`,
  so the pipeline goes straight from these pre-cleaned reads to assembly — host
  filtering and taxonomic profiling are skipped entirely. This is the "core path"
  (assembly → binning → CheckM2 → GTDB-Tk/dRep → final report), which is what every
  change in this repo's recent history has actually touched.
- **No CheckM2 or GTDB-Tk database.** Both are multi-GB, cluster-specific, and
  already need to exist on whatever cluster you're testing on — they can't be
  committed to the repo. `config.test.yaml` leaves `CheckM.db` and `GTDB_Tk.db`
  blank with a comment marking them; **you must fill both in before running**, or
  Snakemake will fail immediately (a blank `GTDB_Tk.db` fails at DAG-build time with
  "Input and output files have to be specified as strings", not just at runtime).

## Running it

```bash
# one-time: these are run outputs, not test data, so they're gitignored and must
# exist before the pipeline's onstart validation will pass
mkdir -p test_data/_run/{output,scratch,logs,benchmarks}

cd workflow
snakemake -n --configfile ../test_data/config.test.yaml          # dry run first
snakemake --use-conda --profile ../profile/slurm --configfile ../test_data/config.test.yaml
```

All paths in `config.test.yaml` are relative to `workflow/` (matching the pipeline's
own convention), so this works unmodified on any clone — only the two database paths
above need a per-cluster edit.

## Known limitations

- Assembly/binning results are not biologically meaningful (synthetic random
  sequence) — this validates that the pipeline *runs end-to-end without crashing*
  and that new rule/script changes wire together correctly, not assembly/binning
  *quality*.
- `Bin_filtering`/`dRep` completeness/contamination thresholds are lowered from
  production (50/10 instead of 90/5) since small synthetic contigs don't score as
  cleanly as a real assembled genome — don't carry these over to a real config.
- The generated reads were validated for FASTQ structural correctness (read/quality
  length matching, valid record structure) but **not test-assembled locally** in the
  session that created them (no SPAdes available) — the first real cluster run is
  the actual end-to-end validation.
