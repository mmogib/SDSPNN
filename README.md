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

`s90_g6_accounting.jl` checks run accounting and row counts for all six artifacts before finalizing them and writing CSV rows and JSON manifests below `results/exports/`. `s95_final_audit.jl` checks final status, accounting, stored selection hashes, and agreement between the table and figure values of the certificate threshold. It prints the recorded gates and numerical evidence without changing the artifacts.

### 4. Generate the figures

Run the plot consumers only after `s90_g6_accounting.jl` has finalized the artifacts:

```powershell
julia --project=. scripts/p30_f1_annulus.jl
julia --project=. scripts/p40_f2_certified_radius.jl
julia --project=. scripts/p50_f4_integrability.jl
julia --project=. scripts/p60_f5_variation_band.jl
```

Each plot script reads a final database artifact and writes one PDF plus a provenance sidecar under `results/figures/`. Plot scripts do not solve an ODE, project a point, differentiate a residual, or write to the artifact database.

## Paper-to-script output map

Artifact identifiers are internal provenance keys. They do not need to match the final numbering in the manuscript.

| Paper output (§7) | Compute script | Artifact ID | Plot consumer | Generated output |
|---|---|---|---|---|
| Table 1 — metric, residual and certificate quantities | `scripts/s20_t1_constants.jl` | `T1_R2_20260801` | — | `results/exports/t1_r2_20260801_rows.csv` |
| Figure 1 — annulus | `scripts/s30_f1_annulus.jl` | `F1_R2_20260801` | `scripts/p30_f1_annulus.jl` | `results/figures/f1_annulus_r2.pdf` |
| Figure 2 — integrability residual, four metrics | `scripts/s50_f4_integrability.jl` | `F4_R2_20260801` | `scripts/p50_f4_integrability.jl` | `results/figures/f4_integrability_residual_r2.pdf` |
| Figure 3 — local certificate at τ = 1 | `scripts/s40_f2_certified_radius.jl` | `F2_R2_20260801` | `scripts/p40_f2_certified_radius.jl` | `results/figures/f2_certified_radius_r2.pdf` |
| Figure 4 — global certificate across τ | `scripts/s60_f5_variation_band.jl` | `F5_R2_20260801` | `scripts/p60_f5_variation_band.jl` | `results/figures/f5_variation_band_r2.pdf` |
| Table 2 — projection at an active constraint | `scripts/s70_active_projection.jl` | `AP1_R2_20260802` | — | `results/exports/ap1_r2_20260802_rows.csv` |

Every artifact has a `_rows.csv` export and a `_manifest.json` file with the lowercase artifact ID as prefix. Table 2's CSV contains per-point data; `s70_active_projection.jl` prints its aggregated table to the console and `results/logs/`. Every figure PDF has a neighboring `_provenance.txt` sidecar recording its manifest and selection-query hashes.

For the numbers quoted in §7, run `scripts/s95_final_audit.jl`: it checks the shared certificate threshold and prints the constants, trajectory errors, radius gap, global-certificate ratio, and active-projection evidence. The table and figure computations above supply those values.

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
- Compare exported CSV data and manifest fields other than `created_utc`, `finalized_utc`, and `final_manifest_hash`. Creation time enters the manifest hash, so it changes on a fresh run even when the data agree; use it to trace a figure to its own run. Floating-point results can also vary across platforms.
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
