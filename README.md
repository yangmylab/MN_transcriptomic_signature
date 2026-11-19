# MN transcriptomic signature

Analysis code, metadata and result tables for a compartment-resolved glomerular transcriptomic
signature and a four-gene classifier for membranous nephropathy (MN).

The study assembles 13 public GEO accessions (2,165 samples; 453 MN samples from 178 unique
patients) into seven non-redundant cohort sets, after patient-level de-duplication. Discovery uses
NEPTUNE glomerular microarrays, independent external validation uses the European Renal cDNA Bank
(ERCB), and cross-platform replication uses two RNA-seq cohorts plus a laser-capture cohort. A
236,107-cell kidney atlas is used for cell-type attribution and NNLS deconvolution.

> Research use only. The four-gene panel (SLAMF8, TRPC6, MYOZ2, PPP1R3C) is a discrimination
> instrument, not a validated diagnostic or prognostic test, and it requires per-platform
> recalibration.

## Repository contents

| Path | Contents |
|---|---|
| `R/` | Numbered analysis scripts (`00_` … `99_verify_*`), one step per script |
| `tools/` | Figure compositor, cohort-design builder, environment record, provenance test |
| `metadata/` | Cohort definitions, Track B panel register and exclusions, PLA2R table, symbol map |
| `Manuscript/` | Manuscript (DOCX), `references.tsv`, figure copies |
| `Docs/` | Consolidated results, findings ledger, audit, MYOZ2 provenance |
| `output/tables/` | Result tables cited by the manuscript (CSV) |
| `output/figures/` | Main and supplementary figure plates (PDF + 300-dpi JPEG) |
| `environment/` | `renv.lock`, `sessionInfo.txt`, package inventories |
| `reference/` | Pandoc reference document used for the manuscript build |

Not included, because it is large and regenerable: `work/` (normalised matrices, networks, atlas,
caches), `.Rlib/` (local package library), `scrna/extracted/`, and per-step intermediate figures.

## Data availability

All primary data are public GEO accessions, listed in `metadata/` and in the manuscript Methods.
Raw and processed GEO archives are not redistributed here. Harmonised matrices are regenerable from
the accessions with the scripts below, given the shared data root described under *Requirements*.

## Requirements

- R 4.6.1 with Bioconductor 3.23
- Python 3.14 (figure composition only)
- Pandoc 3.8.3, TeX Live and Ghostscript (manuscript and figure builds)
- A local package library containing the Brainarray CDF packages, `makecdfenv`, `oligo`, WGCNA
  1.74, `ConsensusClusterPlus`, `caret`, `glmnet`, `randomForest`, `xgboost`, `Boruta`, `metafor`,
  `msigdbr`, `rms` and `babelgene`; the remaining packages are recorded in `environment/renv.lock`

## Reproducing the analysis

```r
install.packages("renv")
renv::restore(lockfile = "environment/renv.lock")
```

Rebuild the local Brainarray CDF library from the GEO archives (no network needed):

```bash
Rscript R/00d_install_cdf.R
Rscript R/02a_build_cdf.R
for p in work/cdf/*; do R CMD INSTALL -l .Rlib "$p"; done
```

Point the scripts at the shared data root, then run the numbered steps in order:

```bash
export PAPER1_DATA=/path/to/MN/Datasets
```

| Steps | Role |
|---|---|
| `03d`, `04a`–`04d` | QC policy, ComBat, limma, random-effects meta-analysis with ERCB replication, enrichment |
| `05a`–`05c` | Separability gate, robustness, partition stability |
| `06a` | Exploratory subtyping (negative result) |
| `07a`–`07z` | Single-cell manifest, QC, atlas, pseudobulk, bridge, deconvolution |
| `08a`–`08f` | WGCNA, traits, stability, preservation, hubs, enrichment |
| `09a`–`09f` | Panel selection, validation, robustness, mechanism, sensitivity, PLA2R |
| `10b`–`10d` | Track B harmonised matrices, attrition ladder, panel-size frontier |
| `99_verify_*` | Verification scripts (see below) |

Steps marked in the manuscript as long-running should be run in the background; the full pipeline
takes hours and requires about 30 GB of working space.

## Verification

Each reported number is asserted by a verifier that fails loudly on any change:

| Verifier | Assertions |
|---|---:|
| `99_verify_phase_a.R` | 54 |
| `99_verify_phase_b.R` | 86 |
| `99_verify_downstream.R` | 32 |
| `99_verify_scrna.R` | 42 |
| `99_verify_phase_d.R` | 117 |
| `99_verify_phase_e.R` | 106 |
| `99_verify_trackb.R` | 37 (Track B benchmark) |

## Citation

If you use this code, please cite the manuscript and the archived release:

- Manuscript: *A compartment-resolved glomerular transcriptomic signature and a four-gene
  classifier for membranous nephropathy across three independent cohorts* (in preparation).
- Code archive DOI: to be issued by Zenodo at the first tagged release.

## Licence

MIT. See `LICENSE`.
