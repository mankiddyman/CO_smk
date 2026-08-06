# Session kickoff: build the CO_smk QC battery

*Paste this as the first message of a new conversation.*

---

I want to build an automated QC/diagnostics battery for my Snakemake pipeline
`CO_smk`, which infers meiotic crossover landscapes from single-cell pollen
RNA-seq. Last week I audited a species by hand across many terminal sessions and
want that work turned into something that runs automatically for every future
species.

## Read these first

Both are committed in the repo:

```bash
cd /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk
cat docs/qc_battery_brief.md    # what to build: stages, plots, numbers, flags
cat docs/spondias_audit.md      # the hand audit; source of truth for expected numbers
cat docs/USAGE.md               # pipeline overview, inputs, gotchas
```

The brief is the spec. The audit is how you verify the build is correct — if a
number the new code computes disagrees with the audit, **the new code is wrong
until proven otherwise**; each of those numbers was derived and cross-checked.

## Environment

- HPC, interactive on `hpc09`. Driver env:
  `micromamba activate /netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk`
- Repo: `/netscratch/dep_mercier/grp_marques/Aaryan/CO_smk`
- Per-rule conda envs under `.snakemake/conda/` (Snakemake-managed, `--use-conda`)
- **`python` is aliased to ipython** — always `python3` for scripts, or argparse
  flags get eaten
- Always `snakemake <target> -n --cores N --use-conda --rerun-triggers mtime params`.
  Put the target **before** the flags; `--rerun-triggers` takes multiple values
  and will swallow a trailing positional.
- **Never `--forceall`** on `Spondias_hap1` — its GTF is hand-staged at
  `results/annotation/Spondias_hap1.gtf` and `gff_to_gtf` would clobber it.

## Repo conventions

- Per-species parameters resolve through `param_for(sample_id, section, key)` in
  the Snakefile, reading `config/config.yaml` blocks: a `default:` plus optional
  species blocks keyed on the `species` column of `samples.csv`. Cuscuta's values
  are **pinned** so its published figures can't drift. New code must read
  parameters this way, not hardcode them.
- Every parameter in config carries its sweep numbers and reasoning as inline
  comments. Keep that up.
- `results/`, `logs/`, `benchmarks/`, `diag/`, `resources/` are gitignored.
  Code lives in `workflow/`, docs in `docs/`.
- Commit style: deliberate staging (never `git add .`), messages that record the
  evidence and not just the change.

## Reference numbers — Spondias_hap1

Verify against these. Full detail in the audit.

| | |
|---|---|
| markers | 1,471,619 · median spacing 33 bp · 2 gaps > 1 Mb genome-wide |
| het SNP depth | q1 79, median 91, q3 104, p90 117 (`max_dp: 150`) |
| filter attrition | AF window 18.5%, QUAL 2.4%, max_dp 2.2%, min_dp 0.1% |
| barcodes | 8,497 called · median 1,127 UMI · 863 genes |
| genotyping depth | 6,001,202 cell×marker obs · **mean DP 1.20 · 89.5% single-read** |
| contamination | 0.57% at well-behaved markers; ~3.5% ambient per droplet |
| cells selected | 153 |
| crossovers | 1,740 · 11.37/cell · 1.42 per bivalent |
| map | **1137 cM over 284 Mb (4.0 cM/Mb)** · 0 chromosomes below the 50 cM floor |
| detection | median 132 markers/chromosome; ~148 needed for 90%; **not saturated** |

## Bugs we already hit — do not repeat

- **Contamination must be measured per *marker*, not per observation.** A fixed
  ratio window (e.g. "clean = ≤0.10 or ≥0.90") makes deeper observations look
  dirtier, because at DP=10 one contaminating read gives exactly 0.10. We got a
  bogus 13%→32% "rising with depth" result this way. Classify markers by their
  aggregate behaviour across cells instead.
- **AD and DP sparse matrices have different non-zero entry sets** — a marker
  with 2 REF reads has a DP entry and no AD entry. Never `paste` the files or zip
  their value columns; join on (row, col) with AD defaulting to 0.
- **`ast.parse` does not validate a Snakefile** (`rule x:` isn't Python). Use
  `snakemake -n`.
- **Every parameter sweep needs an independent freshness check.** Two of mine
  silently reported stale results. Confirm the thing actually changed (e.g.
  `wc -l cellSNP.samples.tsv` per iteration), and never hide stderr with
  `> /dev/null 2>&1` without checking the exit code.
- **Allele-ratio plots at low DP show a discreteness lattice**, not populations —
  at DP=3 only 0, ⅓, ⅔, 1 are achievable. Don't misread the spikes.

## Working style

Bullet points, dense. Practical steps over theory. Scripted edits via heredoc or
a small Python patch script rather than hand-editing — long sessions have mangled
files otherwise. Verify after every edit. Anti-overengineering: build in small
vertical slices, test each against Spondias before moving on.

## Task

Build Tier 1 from the brief: three staged report rules producing
`qc/{sample}/{01_markers,02_cells,03_cos}/` with plots, a machine-readable
`summary.tsv`, and `flags.txt` listing detected anomalies.

Start with stage 01 (markers) end-to-end and verify it against the reference
numbers before touching stage 02. Suggested checklist is at the bottom of the
brief.

Where do you want to start?
