# Handoff — Project Handoff Document

**Project:** Paper 1 — Machine-learning biomarker discovery & molecular subtyping of membranous nephropathy (MN) from multi-cohort bulk transcriptomics
**Workspace (project root):** `/home/yang/Documents/MN/Paper_1`
**Last updated:** 2026-10-18 (analysis work performed 2026-10-15 → 2026-10-16)
**Status:** Phases 0, 1 and A complete and verified. Phase B not started.
**Version control:** none — this directory is **not** a git repository.

> **Formatting note:** the standing instruction is `.md` **plus** a `.docx` (仿宋 for CJK, Times New Roman for Latin, red code). `pandoc` is not installed on this machine, so only this `.md` exists. See §3.

---

## 1. Objective (What we are doing)

Membranous nephropathy is an antibody-mediated glomerular disease with heterogeneous outcomes, and no tissue-molecular biomarker is used clinically. This project asks whether the glomerular and tubulointerstitial transcriptomes of MN contain (a) a compact, compartment-resolved gene signature that separates MN from other proteinuric glomerulopathies (FSGS, MCD, IgAN, lupus nephritis, ANCA-GN) and healthy tissue, and (b) reproducible molecular subtypes with distinct prognostic trajectories.

The plan is a 14-week, multi-cohort bulk-transcriptomic study: harmonise every public MN cohort, discover a machine-learning biomarker panel (LASSO / SVM-RFE / Boruta), and derive molecular subtypes with drug-repositioning follow-up. Target journals are *Kidney International*, *JASN*, *CJASN*, *EBioMedicine*, *Briefings in Bioinformatics*.

**Explicit deliverable of this block:** a verified, non-redundant, normalised expression resource plus QC — i.e. Phases 0, 1 and A of the plan in `MN_Paper1_Bulk_ML_Biomarkers_Subtypes.md`.

**Success criterion met:** `Rscript R/99_verify_phase_a.R` exits 0 with 51/51 checks passing.

---

## 2. Completed Items (What has been done)

### 2.1 Phase 0 — Data consolidation and staging ✅

- All **13 GEO accessions** supplied manually by the user in `/home/yang/Documents/MN/Datasets` (6.6 GB) were inventoried and verified.
- The 10 `*_RAW.tar` archives were extracted into `Paper_1/Datasets/<GSE>/raw/`; the source folder was left **untouched as the immutable master** (user's explicit choice).
- CEL counts per cohort were confirmed against the `tar -tf` listings: 196 / 195 / 111 / 169 / 244 / 193 / 153 / 0 (Agilent) / 175 / 187.
- `R/00c_repair_counts.R` fixed a **corruption in GSE197307**: the header read `ENSEMBL_ene_ID` and every Ensembl ID had lost its `G` (`ENS00000000003`). Repair is deterministic (`sub("^ENS","ENSG")`), guarded by `stopifnot(all(grepl("^ENS[0-9]+$", ids)))`. Also standardised GSE182380 (symbol-keyed) and GSE216841 (dropped its duplicate `ensembl_gene_id` column), and harvested `metadata/ensg2symbol.tsv` (19,847 entries) from GSE216841 — no download needed.

### 2.2 Phase 1 — Harmonisation, de-duplication and design ✅

- `R/01a_build_sample_metadata.R` produced `output/tables/sample_metadata.csv`: **2165 samples** across 13 accessions with accession, GSM, sample title, patient id, compartment, canonical disease, raw disease label, batch, assay, platform, cohort.
  - Canonical disease vocabulary collapses cohort-specific spellings (ERCB uses **`MGN`** for MN, `LD`/`TN` for controls, plus `HT`, `TMD`, `CD`, `FSGS_MCD`, `AAV`).
  - Living-donor ids are numbered 1..n independently of patients, so donors were namespaced to `LD-<n>` to stop them colliding with real patient ids.
- `R/01b_dedup_patients.R` produced `patient_map.csv`, `dedup_summary.csv`, `release_subset_relations.csv`, `recommended_releases.csv`.
- `R/01c_design_options.R` produced `proposed_tiers.csv` and the discovery-vs-validation overlap table.

**Headline intermediate results (these changed the study design):**

| Quantity | Value |
|---|---|
| Total samples across 13 releases | **2165** |
| MN samples | **453** |
| **Unique MN patients** | **178** |
| Redundancy, NEPTUNE-numbered | 355 samples / 119 patients = **2.98×** |
| Redundancy, ERCB-numbered | 78 samples / 39 patients = **2.00×** |

**Four redundancy groups (6 accessions, 7 containment relations):**

- `GSE108112` (TI) ⊂ `GSE133288` — 43/43 MN patients, 100 %
- `GSE108109` (Glom) ⊂ `GSE200828` — 44/44, 100 %
- `GSE104948` (Glom) ≡ `GSE99339` (Glom) — 21/21, mutually identical
- `GSE104954` (TI) ≡ `GSE99325` (TI) — 18/18, mutually identical

**Decisions taken and their rationale (user-approved 2026-10-16):**

| Decision | Rationale |
|---|---|
| Tiering: NEPTUNE microarray = discovery; ERCB microarray = independent external validation; NEPTUNE RNA-seq = cross-platform replication | ERCB shares **0 %** of MN patients with discovery; the RNA-seq tiers do not |
| Drop `GSE108109` and `GSE108112` entirely | 100 % redundant re-releases |
| Keep `GSE99339`/`GSE99325` only as a **normalisation-robustness layer** | Same patients as GSE104948/GSE104954 but a different Brainarray revision |
| RNA-seq tiers keep **all** patients and are labelled **"cross-platform replication"** | NEPTUNE biopsied each patient in both compartments and on both platforms, so these tiers are the same people measured twice, not independent validation |

### 2.3 Phase A — Preprocessing and QC ✅

- `R/00d_install_cdf.R` installed `makecdfenv` into the isolated workspace library.
- `R/02a_build_cdf.R` built **5 Brainarray CDF packages** from the `.cdf.gz` files that ship *inside* the GEO tars — no external download.

| Package | Array | Source |
|---|---|---|
| `hugene21sthsentrezgcdf` | HuGene-2.1-ST | `GPL19983_hugene21st_Hs_ENTREZG.cdf.gz` |
| `hgu133plus2hsentrezgcdf` | HG-U133 Plus 2 | GPL22945 |
| `hgu133plus2hsentrezg2cdf` | HG-U133 Plus 2 (other revision) | GPL19109 |
| `hgu133ahsentrezgcdf` | HG-U133A | GPL24120 |
| `hgu133ahsentrezg2cdf` | HG-U133A (other revision) | GPL19184 |

- `R/02b_rma_arrays.R` normalised **7 array accessions** (GSE133288, GSE200818, GSE200828, GSE104948, GSE104954, GSE99339, GSE99325) → `work/rma/`.
- `R/02c_rnaseq_normalize.R` normalised **3 RNA-seq cohorts** → `work/rnaseq/`.
- `R/03b_build_cohort_sets.R` assembled the **7 cohort sets** → `work/sets/` (native-gene and symbol-keyed versions plus a samples table each).
- `R/03c_qc.R` produced per-cohort PCA + outlier screen → `output/qc/` (7 PDF + 7 JPEG, per the standing preference of PDF + JPEG, no grid lines, step-prefixed filenames).
- `R/99_verify_phase_a.R` — **51 assertions, all passing, exit 0.**

**Key parameters and thresholds actually used:**

| Step | Parameter |
|---|---|
| Array normalisation (ST) | RMA background (`oligo::backgroundCorrect(method="rma")`) → quantile (`oligo::normalize`) → log2 → **median polish per Brainarray gene** (`stats::medpolish`, the `overall + col` term) |
| Array normalisation (U133) | `affy::rma()` with the matching Brainarray CDF |
| Chunking | `CHUNK = 40` samples per chunk, `NCPU = min(6, detectCores()-1)`; two-pass quantile reference across the whole cohort |
| RNA-seq | `edgeR::filterByExpr` → TMM → `limma::voom` (log2-CPM `E`) |
| GSE216841 | already-normalised decimals → `log2(x+1)`; **not** re-counted |
| QC PCA | top 2000 variable genes, `prcomp(scale.=TRUE)` |
| QC outlier flag | median sample–sample correlation < `median − 3 × MAD` |

**The 7 cohort sets (final):**

| Set | Releases | Samples | Genes | MN | FSGS | MCD | Other | Ctrl | Role |
|---|---|---|---|---|---|---|---|---|---|
| `disc_glom_array` | GSE200828 | 153 | 25,582 | 51 | 31 | 19 | 46 | 6 | Discovery (glomerular) |
| `disc_ti_array` | GSE133288 + GSE200818 | 269 | 25,582 | 51 | 78 | 63 | 72 | 5 | Discovery (tubulointerstitial) |
| `val_glom_array` | GSE104948 | 196 | 19,769 | 21 | 18 | 14 | 0 | 26 | **Independent** external validation |
| `val_ti_array` | GSE104954 | 195 | 19,769 | 18 | 13 | 13 | 0 | 26 | **Independent** external validation |
| `cross_glom_rnaseq` | GSE197307 | 282 | 15,739 | 62 | 93 | 89 | 29 | 8 | Cross-platform replication |
| `cross_ti_rnaseq` | GSE182380 | 199 | 12,605 | 59 | 70 | 32 | 37 | 0 | Cross-platform replication |
| `pla2r_lcm` | GSE216841 | 34 | 20,321 | 12 | 0 | 14 | 0 | 8 | PLA2R+/− sub-analysis |

**Leakage accounting (verified):** `val_glom_array` and `val_ti_array` have **0** patients overlapping discovery; `cross_glom_rnaseq` has **208** overlapping samples and `cross_ti_rnaseq` **168** — recorded in the `role` / `MN_overlapping_discovery` columns of `output/tables/cohort_sets_summary.csv`, not hidden.

**QC result (`output/tables/qc_summary.csv`):**

| Set | PC1 % | PC2 % | min median-corr | flagged |
|---|---|---|---|---|
| `disc_glom_array` | 12.0 | 9.0 | 0.458 | 3 |
| `disc_ti_array` | 17.5 | 11.4 | 0.281 | 29 |
| `val_glom_array` | 59.8 | 8.4 | 0.459 | 71 |
| `val_ti_array` | 68.5 | 10.3 | 0.528 | 64 |
| `cross_glom_rnaseq` | 14.5 | 7.3 | 0.188 | 7 |
| `cross_ti_rnaseq` | 22.7 | 7.3 | 0.183 | 26 |
| `pla2r_lcm` | 18.2 | 10.8 | 0.247 | 0 |

**Critical validation result:** the HuGene-2.1-ST gene-level RMA reimplementation was checked against the authors' own Brainarray matrix (`GSE200828_processed_data.txt.gz`) — **Pearson r = 0.985–0.987 per sample, overall 0.9859**, across all 25,582 genes.

---

## 3. Current Blockers (Where we are stuck)

### B1 — `.docx` cannot be generated (standing formatting rule unmet)
- **Symptom:** `pandoc` is not installed (`which pandoc` → not found).
- **Impact:** the standing rule "every generated `.md` also gets a `.docx` with 仿宋 + Times New Roman + red code" is unmet for this and all previous documents.
- **Tried:** nothing — the tool is simply absent. The reference template from memory lived on a Windows machine (`D:\Wechat_Miniprogram\...\tnr_fangsong_redcode.docx`) and was already noted as needing regeneration on this machine.
- **Needs:** install `pandoc`, then regenerate the reference doc.

### B2 — Clinical variables not acquired
- **Symptom:** NEPTUNE baseline eGFR, UPCR, and longitudinal remission/relapse status are not in any downloaded file.
- **Impact:** Phase C step 11 (correlate infiltrates with eGFR/UPCR) and Phase F step 21 (subtype–outcome Cox models, eGFR slope) cannot run.
- **Needs:** NEPTUNE data-sharing application, or per-sample annotations from published supplements. 【To be confirmed with the user which route he prefers】

### B3 — Deconvolution packages not installed
- **Symptom:** CIBERSORTx is not an R package (needs the web/Linux binary and a licence token); `xCell`, `MCP-counter`/`immunedeconv`, `quanTIseq` are **not** present in either library.
- **Impact:** Phase C blocked.
- **Tried:** package inventory only.

### B4 — PCA outlier decisions pending
- **Symptom:** the median-correlation screen flags 200 samples across the 7 sets (dominant: `val_glom_array` 71, `val_ti_array` 64).
- **Impact:** whether any are excluded changes every downstream DEG/model result.
- **Status:** these are **candidates, not exclusions**. Awaiting the user's decision. 【Requires user】

### B5 — GSE73953 (PBMC arm) not normalised
- **Symptom:** it is an Agilent one-channel array (`GPL4133`, plain `.txt`, no CEL) and was never added to `02b_rma_arrays.R`.
- **Impact:** the optional peripheral cross-check is unavailable. Low priority (8 MN samples).

### B6 — Two series matrices are incomplete
- **Symptom:** `GSE108112_series_matrix.txt.gz.part` (18 MB) and `GSE133288_series_matrix.txt.gz.part` (22 MB) fail `gzip -t`.
- **Impact:** **none** — both cohorts' phenotypes were recovered from CEL filenames (`B1_TI_MN_94`), and GSE108112 is dropped by design anyway.

---

## 4. Next Steps (What to do next)

### P0 — Decide the PCA outlier policy 【owner: user】
- **Resume from:** `output/tables/qc_summary.csv` and `output/qc/03c_<set>_sample_qc.csv`.
- **Action:** review the flagged samples; decide none / correlation-only / manual review. Anything excluded must be recorded in a table for the methods section.
- **Blocks:** Phase B and everything after.

### P1 — Phase B: batch integration + differential expression 【owner: agent】
- **Resume from:** `work/sets/*_expr.rds` (or `*_expr_symbol.rds` for cross-platform views) and `output/tables/cohort_sets_summary.csv`.
- **Action:** ComBat on log-intensities per platform family (study as batch, disease as biological covariate) → limma DEGs per compartment → random-effects meta-analysis (`metafor`) → GO/KEGG/Reactome (`clusterProfiler`) and GSEA (Hallmark, C7).
- **Thresholds from the plan:** |log2FC| > 1, BH-adjusted P < 0.05.
- **Must NOT** merge GSE133288 with GSE200818 blindly — `disc_ti_array` is already patient-deduplicated, so do not re-introduce GSE108112.
- **Blocks:** Phase D (WGCNA) and Phase E (ML panel).

### P2 — Phase C: deconvolution 【owner: agent, blocked by B3】
- **Resume from:** the same `work/sets/` matrices.
- **Action:** install `immunedeconv`/`xCell`/`MCP-counter`/`quanTIseq`; CIBERSORTx requires the user to obtain the binary + token.

### P3 — Acquire clinical variables 【owner: user】
- See B2. Blocks the outcome analyses in Phases C and F.

### P4 — Optional: normalise GSE73953 【owner: agent】
- Add an Agilent branch to `R/02b_rma_arrays.R` or a small sibling script reading the `.txt` files in `Datasets/GSE73953/raw/`.

### P5 — Reproducibility package 【owner: agent, deferred to Phase G】
- No `renv.lock` exists. The plan requires one. Generate it early, because the environment is split across two libraries (see §6).

---

## 5. Pitfalls & Anti‑patterns (Never repeat these)

### P-1 — Do not trust the plan's cohort arithmetic
The plan asserted ">500 MN samples" and instructed merging releases that are subsets of one another. **Always compute patient-level overlap across *all* releases before building any merged matrix.** Concretely: 13 releases → 178 patients, and 4 releases are 100 % contained in another.

### P-2 — Never estimate cross-cohort overlap from one disease class or one compartment
I first reported "13 shared patients" between the microarray discovery set and the RNA-seq layer. That was wrong: it compared MN patients within the *same* compartment only. The correct figure is **50 / 51 MN patients**, and **208 / 282 samples (74 %)** overall, because NEPTUNE biopsied each patient in both compartments and on both platforms. **Compute overlap over every disease class and every accession before claiming independence.**

### P-3 — `affy` cannot process HuGene-2.1-ST. Do not try.
Both `ReadAffy(cdfname=…)` and `justRMA(cdfname=…)` fail with *"The affy package is not designed for this array type."* Use `oligo` at probe level plus the Brainarray CDF probe indices (see the working implementation in `R/02b_rma_arrays.R`).

### P-4 — `preprocessCore::rcModelPLM()$Estimates` is NOT the sample effect vector
It has length `nrow + ncol` (probe effects then sample effects). A length check against `ncol` silently fails and drops you into a slow fallback. The validated path in this project uses `stats::medpolish()`'s `overall + col` instead — keep it, it is what produced r = 0.986.

### P-5 — Never call `match()` against the 1.02 M-row probe index per gene
It rebuilds a hash table on every call: **436 s** for 25,582 genes. Build a lookup vector once (`pos_of[row_names] <- seq_along(row_names)`).

### P-6 — Never process a 244-array ST cohort in one gulp on this box
It reached **12.1 GB RSS of 14 GB** and had to be killed. Use the chunked two-pass quantile design already implemented (`CHUNK = 40`).

### P-7 — CDF package names must end in `cdf`
`affy`'s `cleancdfname()` appends `cdf` to any name that does not already end in it. A package called `…gcdfb` becomes a lookup for `…gcdfbcdf` and fails.

### P-8 — `make.cdf.package()` has two traps
It only *builds* a source package (you must then `R CMD INSTALL -l .Rlib <dir>`), and it resolves its input as `file.path(getwd(), filename)` — passing an absolute path doubles it. Pass a bare filename after `setwd()`.

### P-9 — Do not derive the array platform from the series matrix
For two-platform releases (GSE104948/104954/99339/99325) the two series matrices list **overlapping** GSM ids, so a `match()` picks the wrong one and `affy` dies with *"does not seem to have the correct dimensions"*. Read each CEL's own header instead: `affyio::read.celfile.header(f, info="basic")$cdfName`.

### P-10 — Do not `rbind()` two platforms of one release
They share samples but not genes. Align on the **union of genes** and fill `NA` where absent (implemented in `R/02b_rma_arrays.R`), then drop incomplete genes before PCA/correlation (implemented in `R/03c_qc.R`).

### P-11 — Do not redirect a script's stdout to the same file its `writeLines()` writes
The script's own log write truncates the captured stdout. Use `.out` for the shell redirect and let the script write its `.log`.

### P-12 — Do not assume count-style filenames contain counts
`GSE216841_MNd_ncounts_annot.txt.gz` holds **already-normalised decimals**; it also carries `hgnc_GRCh38p12` and `description` columns beside the gene id, so a naive "first column = id, rest = samples" reader silently treats text columns as samples (this produced `NA`s and 36 "samples" instead of 34).

### P-13 — Do not use raw integer donor ids as patient ids
Living donors are numbered 1..n independently of patients, so `1_Living_donor` collides with patient `1`. Namespace them (`LD-<n>`). Symptom before the fix: controls silently deduplicated away and comparator counts were wrong.

### P-14 — Data acquisition is manual in this project
Per a standing user rule, **never download datasets**. Present a concrete request list and let the user download. (Package installation is *not* covered by that rule and was performed directly.)

---

## 6. Key Resources & Environment Info

### 6.1 Compute
| Item | Value |
|---|---|
| OS / shell | Linux amd64, bash |
| CPU | 16 cores |
| RAM | **14 GB total** (~11 GB usable) — this is the binding constraint |
| Disk | 470 GB volume, **360 GB free** |
| Project footprint | 25 GB total: `Datasets/` 6.8 GB, `work/` 18 GB |

### 6.2 Languages and libraries
- **R 4.6.1 (2026-06-24) "Happy Hop"**
- **Two libraries are in play** — this is essential to reproduce anything:
  1. Workspace library `/home/yang/Documents/MN/Paper_1/.Rlib` — **48 packages**, prepended by `R/_lib.R` (`.libPaths(c(WS_LIB, .libPaths()))`). Holds `makecdfenv`, the 5 Brainarray CDF packages, `oligo`, `oligoClasses`, `affyPLM`, `arrayQualityMetrics`, `WGCNA`, `ConsensusClusterPlus`, `caret`, `glmnet`, `randomForest`, `xgboost`, `Boruta`, `metafor`, `msigdbr`, `rms`, `babelgene`, `hugene21sttranscriptcluster.db`, `hgu133b.db`, and since this session `pd.hugene.2.1.st`.
  2. User library `/home/yang/R/x86_64-pc-linux-gnu-library/4.6` — **794 packages**. Holds `limma`, `edgeR`, `DESeq2`, `sva`, `GEOquery`, `Biobase`, `annotate`, `preprocessCore`, `affy`, `affyio`, `BiocManager`, `GSVA`, `clusterProfiler`, `org.Hs.eg.db`, `hgu133plus2.db`, `hgu133a.db`, `ggplot2`, `pheatmap`, `ComplexHeatmap`, `e1071`, `biomaRt`, `ensembldb`, `reshape2`.
  3. System library `/usr/lib/R/site-library` is **empty**; `/usr/lib/R/library` holds only base packages.
- **R package mirrors used:** `https://mirrors.tuna.tsinghua.edu.cn/CRAN/` and `https://mirrors.tuna.tsinghua.edu.cn/bioconductor` (set in `R/00d_install_cdf.R`).
- **No `renv.lock`, no `sessionInfo()` capture, no git.** 【Gap — see P5】

### 6.3 Data locations
| Path | Contents |
|---|---|
| `/home/yang/Documents/MN/Datasets` | **Immutable master**: the user's manual download, 6.6 GB, 13 accessions (10 `*_RAW.tar`, 3 count matrices, 2 processed matrices) |
| `Paper_1/Datasets/<GSE>/raw/` | Extracted CEL files + the Brainarray `.cdf.gz` / mapping files shipped in the tars |
| `Paper_1/Datasets/<GSE>/processed/` | Repaired count matrices (GSE197307/GSE182380/GSE216841) and the authors' `GSE200818/200828_processed_data.txt.gz` |
| `Paper_1/Datasets/<GSE>/*_series_matrix.txt.gz` | GEO metadata (11 of 13 intact; see B6) |
| `Paper_1/metadata/ensg2symbol.tsv` | 19,847-entry Ensembl→symbol map harvested from GSE216841 |
| `Paper_1/work/rma/` | 7 normalised array cohorts (`.rds` + `.txt.gz`) |
| `Paper_1/work/rnaseq/` | 3 normalised RNA-seq cohorts (+ `DGEList` objects for the two count cohorts) |
| `Paper_1/work/sets/` | The 7 cohort sets: `*_expr.rds`, `*_expr_symbol.rds`, `*_samples.csv` |
| `Paper_1/work/cel/` | **16 GB** of decompressed CELs — deletable once `work/rma/` is final |
| `Paper_1/work/rma_cache/` | 1.1 GB chunk caches — deletable once `work/rma/` is final |
| `Paper_1/work/cdf/` | Built CDF source packages (provenance for the `.Rlib` installs) |
| `Paper_1/output/tables/` | 10 result tables |
| `Paper_1/output/qc/` | 7 PCA PDF + 7 JPEG + 7 per-sample QC tables |
| `Paper_1/logs/` | Run logs (`rma_arrays.log`, `rnaseq_normalize.log`, `qc.log`, `build_cohort_sets.log`, `repair_counts.log`, plus the earlier download logs) |

**Reference/annotation:** Brainarray Human ENTREZG custom CDFs, version 19 (shipped inside the GEO tars); HGNC symbols via `org.Hs.eg.db`. No genome build is used — all data are pre-quantified expression matrices.

### 6.4 Cross-machine resume
1. Copy the whole project directory (at minimum `R/`, `Datasets/`, `metadata/`, `output/`, and the `*.rds` under `work/`; `work/cel` and `work/rma_cache` are regenerable).
2. Install R 4.6.x, then restore the user library packages listed in §6.2(2) from CRAN/Bioconductor.
3. Rebuild the workspace library: `Rscript R/00d_install_cdf.R` then `Rscript R/02a_build_cdf.R` followed by `R CMD INSTALL -l .Rlib work/cdf/<pkg>` for each of the 5 CDF packages (they are built from the `.cdf.gz` files inside `Datasets/*/raw/`, so no network needed).
4. Sanity check: `Rscript R/99_verify_phase_a.R` must exit 0.

### 6.5 API keys / credentials
None used. CIBERSORTx (needed later, see B3) requires a **licence token** obtained from the CIBERSORTx site — it does not exist on this machine.

---

## 7. Startup Instruction for the Next Session

> **When starting a new session, the very first message must be: "Read handoff.md first."**

The new session should then restate the current state from this document to confirm understanding — specifically: 13 accessions → 2165 samples → 453 MN samples → **178 unique MN patients** → 7 non-redundant cohort sets → Phase A verified (51/51) → Phase B not started.

If the previous session was interrupted, first record what actually completed before the interruption. The last verified state is: **`Rscript R/99_verify_phase_a.R` exits 0.** If that command now fails, the interruption damaged something and must be diagnosed before any new analysis.
