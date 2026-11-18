# Paper 1 — scRNA-seq arm: data request list

**Status:** waiting for user download (standing rule: the agent never downloads datasets).
**Target directory:** `scrna/raw/<GSE>/` (project root moved 2026-10-19 to `/run/media/meng/MY_Data/MN/Paper_1`; `R/_lib.R` resolves it automatically)
**Purpose:** cell-type resolution of the bulk signature, a deconvolution reference matrix (unblocks Phase C without CIBERSORTx), and single-cell validation of the Phase E panel.

All accessions below were checked against NCBI GEO on **2026-10-18** (supplementary-file listing only — sizes are the FTP-reported archive sizes).

## Tier 1 — required for the Paper-1 integration

| # | Accession | Layer | Design | Files to download | Size |
|---|---|---|---|---|---|
| 1 | **GSE271188** | Kidney **snRNA-seq** | 6 MN, 13 IgAN, 7 MCD, 7 control, 4 DKD (n = 37) | `GSE271188_RAW.tar` (+ `filelist.txt` if available) | **2.9 GB** |
| 2 | **GSE171458** | Kidney scRNA-seq | 6 IMN (anti-PLA2R+) + 2 healthy | `GSE171458_RAW.tar` | **71 MB** |
| 3 | **GSE241302** | Kidney scRNA-seq | 3 IMN + 1 healthy | `GSE241302_RAW.tar` | **171 MB** |

Tier-1 total ≈ **3.1 GB**. GSE271188 is the backbone: it is the only cohort with MN *and* IgAN/MCD/control in one study, which is exactly the comparator structure Paper 1 needs.

## Tier 2 — optional, strengthens specific claims

| # | Accession | Layer | Why | Files | Size |
|---|---|---|---|---|---|
| 4 | GSE286184 | Xenium spatial | localise panel genes / circuits in situ | `GSE286184_RAW.tar` | 320 MB |
| 5 | GSE283826 | PBMC scRNA (6 MN + 3 healthy) | circulating counterpart of the tissue signature | `GSE283826_PBMC_with_annotation.rds.gz` | 366 MB |
| 6 | GSE341102 | PBMC scRNA (5 treatment-naïve MN + 3 healthy) | independent PBMC replication | `GSE341102_RAW.tar` | 502 MB |
| 7 | GSE199850 | Kidney scRNA (1 HBV-MN) | secondary-MN contrast | `GSE199850_RAW.tar` | 6.3 MB |
| 8 | GSE193512 | Urinary cells (2 MN, remission vs no remission) | treatment-response signal | `GSE193512_RAW.tar` | 9.2 MB |

## How to place the files

```
Paper_1/scrna/raw/
├── GSE271188/GSE271188_RAW.tar
├── GSE171458/GSE171458_RAW.tar
└── GSE241302/GSE241302_RAW.tar
```

Leaving the archives as `.tar` is fine — `R/07a_scrna_manifest.R` extracts them into `scrna/extracted/<GSE>/` on first run (and never touches the archives). If you prefer to extract manually, do it there yourself; the manifest checker accepts either state.

## What I do once the files are in place

1. `Rscript R/07a_scrna_manifest.R` — validates the archives, extracts, writes `scrna/samples.csv` with one row per 10x sample and its disease label.
2. Import → QC → integration → annotation → pseudobulk → bulk mapping → deconvolution, per `INTEGRATION_PLAN.md`.
