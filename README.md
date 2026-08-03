# SDSPNN

`SDSPNN` is the self-contained Julia numerical package accompanying the paper on a state-dependent metric projection neural network for variational inequalities. It studies the flow

\[
\dot{x}=\lambda\left[P_{S,M(x)^{-1}}\!\left(x-\alpha M(x)F(x)\right)-x\right],
\]

where `F` is the variational-inequality operator, `S` is the feasible set, and the projection metric depends on the current state. The repository reproduces the numerical tables, stored figure data, audit records, and four figures used in the numerical section. It is a reproducibility package, not a general-purpose solver library.

All runtime inputs are encoded in the repository. No external dataset or parent-directory file is required. Every generated file stays below `results/`, which is intentionally excluded from Git.

## Requirements

- Julia 1.12.6 is recommended. It is the version recorded in `Manifest.toml`.
- `Project.toml` declares compatibility with Julia 1.10 or newer.
- Git is needed only to clone the repository.
The first package installation can take several minutes because the environment includes the plotting stack.

## Clone and instantiate

```powershell
git clone https://github.com/mmogib/SDSPNN.git
cd SDSPNN
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The committed `Project.toml` and `Manifest.toml` define the complete Julia environment. Run the commands below from the repository root.

For the closest match to the recorded environment, use one Julia thread. In PowerShell:

```powershell
$env:JULIA_NUM_THREADS = '1'
```

In Bash or another POSIX shell:

```bash
export JULIA_NUM_THREADS=1
```

## Reproduce the numerical section

The stages are ordered. Compute scripts have the prefix `sNN_`; plot-only consumers have the prefix `pNN_`.

### 1. Validate the source and numerical primitives

```powershell
julia --project=. scripts/s00_parse_checks.jl
julia --project=. test/runtests.jl
julia --project=. scripts/s10_record_gates.jl
```

The parser checks every Julia script without executing it. The test suite checks the analytic constants, metric bounds, projections, integrator verification, integrability residuals, trajectory classifier, artifact accounting, and repository path containment. The gate script records the numerical gate evidence in the artifact database.

### 2. Compute the table and figure artifacts

```powershell
julia --project=. scripts/s20_t1_constants.jl
julia --project=. scripts/s30_f1_annulus.jl
julia --project=. scripts/s40_f2_certified_radius.jl
julia --project=. scripts/s50_f4_integrability.jl
julia --project=. scripts/s60_f5_variation_band.jl
julia --project=. scripts/s70_active_projection.jl
```

These scripts create pending artifacts in `results/experiments.db`. They do not draw figures.

### 3. Validate, finalize, export, and audit

```powershell
julia --project=. scripts/s90_g6_accounting.jl
julia --project=. scripts/s95_final_audit.jl
```

`s90_g6_accounting.jl` checks run accounting and row counts before it atomically promotes every artifact from `pending` to `final`. Only then does it write the CSV rows and JSON manifests below `results/exports/`. `s95_final_audit.jl` is read-only. It revalidates the final manifests, selection hashes, gates, and the values used in the numerical discussion.

### 4. Generate the figures

Run the plot consumers only after `s90_g6_accounting.jl` has finalized the artifacts:

```powershell
julia --project=. scripts/p30_f1_annulus.jl
julia --project=. scripts/p40_f2_certified_radius.jl
julia --project=. scripts/p50_f4_integrability.jl
julia --project=. scripts/p60_f5_variation_band.jl
```

Each plot script reads a final database artifact and writes one PDF plus a provenance sidecar under `results/figures/`. Plot scripts do not solve an ODE, project a point, differentiate a residual, or write to the artifact database.

## Output map

Artifact identifiers are internal provenance keys. They do not need to match the final numbering in the manuscript.

| Numerical output | Compute script | Artifact ID | Plot consumer | Generated output |
|---|---|---|---|---|
| Hypothesis and constants table | `scripts/s20_t1_constants.jl` | `T1_R2_20260801` | — | `results/exports/t1_r2_20260801_*` |
| Annulus trajectories and phase portrait | `scripts/s30_f1_annulus.jl` | `F1_R2_20260801` | `scripts/p30_f1_annulus.jl` | `results/figures/f1_annulus_r2.pdf` |
| Certified radius and theorem envelope | `scripts/s40_f2_certified_radius.jl` | `F2_R2_20260801` | `scripts/p40_f2_certified_radius.jl` | `results/figures/f2_certified_radius_r2.pdf` |
| Integrability residual and controls | `scripts/s50_f4_integrability.jl` | `F4_R2_20260801` | `scripts/p50_f4_integrability.jl` | `results/figures/f4_integrability_residual_r2.pdf` |
| Metric-variation certification band | `scripts/s60_f5_variation_band.jl` | `F5_R2_20260801` | `scripts/p60_f5_variation_band.jl` | `results/figures/f5_variation_band_r2.pdf` |
| Active weighted-projection table | `scripts/s70_active_projection.jl` | `AP1_R2_20260802` | — | `results/exports/ap1_r2_20260802_*` |

Every export prefix expands to a row file ending in `_rows.csv` and a manifest ending in `_manifest.json`. Every figure PDF has a neighboring `_provenance.txt` sidecar that records the artifact manifest hash and selection-query hash.

## Generated directory structure

After a complete run, the ignored `results/` tree has this form:

```text
results/
├── experiments.db
├── exports/
│   ├── *_rows.csv
│   └── *_manifest.json
├── figures/
│   ├── *.pdf
│   └── *_provenance.txt
└── logs/
    └── *.log
```

Only `results/.gitkeep` is tracked. The database, exports, logs, figures, and provenance sidecars are generated locally and never committed.

## Artifact lifecycle and reruns

Compute scripts use skip-by-default resume. If an artifact already exists, its compute script leaves it unchanged. To replace one artifact deliberately, pass `--force` to that compute script, then rerun finalization and the relevant plot consumer. For example:

```powershell
julia --project=. scripts/s30_f1_annulus.jl --force
julia --project=. scripts/s90_g6_accounting.jl
julia --project=. scripts/p30_f1_annulus.jl
```

The database is the single source for plots. Plot consumers refuse missing or non-final artifacts. This prevents a figure from mixing rows from an incomplete run.

## Reproducibility and comparison

- Production code resolves all paths from the repository root exported as `JCODE_ROOT`.
- The compute pipeline is deterministic and uses no random input.
- Compute and plot stages are deliberately separate. Every plotted number must already exist in a finalized artifact.
- The exported CSV files and JSON manifest hashes are the preferred comparison targets.
- Raw PDF hashes can differ across operating systems because fonts, LaTeX installations, or plotting backends may embed environment-specific metadata.

For a concise machine-readable review of the completed run, use:

```powershell
julia --project=. scripts/s95_final_audit.jl
```

## Troubleshooting

**A plot reports a missing or non-final artifact.** Run all compute scripts, then `scripts/s90_g6_accounting.jl`, before running any plot consumer.

**A compute script says the artifact already exists.** This is the normal resume behavior. Use `--force` only when you intend to replace that artifact.

**The manifest was created by a newer Julia version.** Use Julia 1.12.6 for the recorded environment. If you resolve dependencies with another Julia version, the environment is no longer an exact manifest reproduction.

**The first plot command is slow.** Julia may be compiling `Plots`, `PGFPlotsX`, and their dependencies. This does not rerun the numerical experiment.

**You launched a script from another directory.** The implementation still writes under this repository's `results/`, but running from the repository root keeps the commands and environment selection unambiguous.

## Source layout

```text
src/       SDSPNN module, problems, metrics, projections, solver, residuals, and artifact I/O
scripts/   ordered compute, audit, and plot entry points
test/      numerical, accounting, active-projection, and path-containment tests
results/   generated local artifacts; ignored except for .gitkeep
```

## License

This code is available under the [MIT License](LICENSE).
