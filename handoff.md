# Handoff — Project Handoff Document

**Project:** Paper 1 — Machine-learning biomarker discovery & molecular subtyping of membranous nephropathy (MN) from multi-cohort bulk transcriptomics
**Workspace (project root):** `/run/media/meng/MY_Data/MN/Paper_1` (moved 2026-10-19 from the old `/home/yang/Documents/MN/Paper_1`; every script resolves the root from `R/_lib.R`, so no path edits were needed)
**Data root (shared):** `/run/media/meng/MY_Data/MN/Datasets` — the GEO inputs live **outside** the project here, shared with Paper_2/Paper_3. `R/_lib.R` sets `DATA` to it and `scrna/manifest.tsv` points its archives at `../Datasets/<GSE>_RAW.tar`.
**Last updated:** 2026-11-18 12:15 — **checkpoint refresh.** Three blocks of work today, all verified: (i) **F62** corrected the gate's cross-platform transfer (§2.20) — 05a/05b now quantile-map each validation cohort onto the discovery reference (label-blind) before applying training statistics, giving ERCB **0.861–0.922** with no per-diagnosis inversions; (ii) a **`bioanalyst` BAP audit** (§2.21) — verdict *reliable with important limitations*, one Major metadata defect fixed (`gsm` mixed 1,623 real accessions with 473 local indices, 134 colliding; provenance columns + 3 assertions), plus **P0 gate-stability** (§2.23: 10 partitions, ERCB 0.912–0.928 at K=all, range 0.016) and its own Fig S6; (iii) the **figure set re-laid-out** (§2.22) — Fig 3 is now the separability gate, functional profiling moved to S4, new S5/S6, and four two-page plates fixed. Six verifiers pass **54 / 86 / 42 / 117 / 106 / 32**. Deliverables: `Docs/{BAP_AUDIT_2026-11-18.md, Findings_Ledger.csv, CONSOLIDATED_RESULTS.md}`, `Manuscript/SUBMISSION_CHECKLIST.md`, **12 plates**.
**Status:** Phases 0, 1, A, B, C, D and E complete, **re-derived under the audited pipeline and fully re-verified**; the figure set and the evidence ledger are submission-ready from the computational side. Locked panel = **4 genes** with ERCB AUC **0.816 (0.721–0.910)**, RNA-seq 0.901, LCM 0.929, other-CDF 0.848, tubulointerstitial 0.449/0.539 (compartment-specific); label-permutation p = 0.005, refit null p = 0.11 (**not** significant — reported as such). The gate's external ERCB AUC is **0.861–0.922** (harmonised transfer, F07 + F62) and is **partition-stable** (0.912–0.928 across 10 draws); the whole-transcriptome model beats the four-gene panel in ERCB (0.922 vs 0.816), so the panel's claim is interpretability and transportability, not superiority. Manuscript is in EBioMedicine format with 44 verified references, a **377-word** Summary and a **4,336-word** main text (which now reports the pathway result), 12 plates (Fig 1–6, S1–S6), and rebuilt `.docx` files. Remaining work is external: clinical variables (B2), participant map (B7), and IF validation of the panel proteins (E4).
**Version control:** **git, initialised 2026-11-18** (**326** tracked files, 24.5 MB of content). `git log` is authoritative (do not trust a commit count written here): `6f8cebe` initial import, `c7ba326` ignore `Rplots.pdf`, then handoff corrections. `core.fileMode` and `core.autocrlf` are both `false` (exfat). `work/`, `.Rlib/`, `scrna/extracted/`, `logs/` and five bulk CSVs >3 MB are ignored as regenerable; see `.gitignore` for the exact rules.

---

## 1. Objective (What we are doing)

Membranous nephropathy is an antibody-mediated glomerular disease with heterogeneous outcomes, and no tissue-molecular biomarker is used clinically. This project asks whether the glomerular and tubulointerstitial transcriptomes of MN contain (a) a compact, compartment-resolved gene signature that separates MN from other proteinuric glomerulopathies (FSGS, MCD, IgAN, lupus nephritis, ANCA-GN) and healthy tissue, and (b) reproducible molecular subtypes with distinct prognostic trajectories.

The plan is a 14-week, multi-cohort bulk-transcriptomic study: harmonise every public MN cohort, discover a machine-learning biomarker panel (LASSO / SVM-RFE / Boruta), and derive molecular subtypes with drug-repositioning follow-up. Target journals are *Kidney International*, *JASN*, *CJASN*, *EBioMedicine*, *Briefings in Bioinformatics*.

**Deliverables completed so far:** (i) a verified, non-redundant, normalised expression resource plus QC (Phases 0, 1, A); (ii) a batch-integrated discovery matrix, per-cohort limma DEGs, a random-effects meta-DEG signature with I² and a held-out ERCB validation of every meta-DEG, and GO/Reactome/KEGG + GSEA profiling (Phase B); (iii) a separability gate showing the ML signature is learnable (external AUC 0.87–0.92) with exploratory subtypes reported as negative; (iv) a kidney single-cell atlas arm delivering **Phase C: cell-type attribution, a signature matrix and NNLS deconvolution of all 7 bulk cohorts**; (v) **Phase D: signed co-expression networks per compartment, module–trait mapping that survives composition adjustment, a reproducibility audit, ERCB preservation of every module, and the hub-gene Phase E candidate space**; (vi) **Phase E: a locked 5-gene panel selected by three orthogonal methods inside nested CV and evaluated once in ERCB and the independent RNA-seq cohort**.

**Success criteria met:** `Rscript R/99_verify_phase_a.R` exits 0 (51/51), `99_verify_phase_b.R` exits 0 (86/86), `99_verify_downstream.R` exits 0 (27/27), `99_verify_scrna.R` exits 0 (42/42), `99_verify_phase_d.R` exits 0 (117/117) and **`99_verify_phase_e.R` exits 0 (102/102)**.

---

## 2. Completed Items (What has been done)

### 2.1 Phase 0 — Data consolidation and staging ✅

- All **13 GEO accessions** supplied manually by the user (6.6 GB) were inventoried and verified; their extracted copies now live in the **shared** `/run/media/meng/MY_Data/MN/Datasets` (18 GB with the archives), reachable as `DATA` — see §6.3.
- The 10 `*_RAW.tar` archives were extracted into `Datasets/<GSE>/raw/`; the source folder was left **untouched as the immutable master** (user's explicit choice). That source folder is now gone — the in-project `Datasets/` copy is the master (§6.3).
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
| Tiering: NEPTUNE microarray = discovery; ERCB microarray = independent external validation; NEPTUNE RNA-seq = separate RNA-seq cohort (independent) | ERCB shares **0 %** of MN patients with discovery; the RNA-seq series carry **no participant id**, so they cannot be linked to the array patients at all (see §2.3 and P-2) |
| Drop `GSE108109` and `GSE108112` entirely | 100 % redundant re-releases |
| Keep `GSE99339`/`GSE99325` only as a **normalisation-robustness layer** | Same patients as GSE104948/GSE104954 but a different Brainarray revision |
| RNA-seq tiers keep **all** patients and are labelled **"independent RNA-seq cohort"**; `MN_overlapping_discovery` is left `NA` | GSE197307/GSE182380 sample titles are local sample indices and their series matrices carry only disease + tissue — patient-level linkage is impossible. An earlier "50/51 shared MN patients, 208/282 shared samples" figure was an identifier collision and has been **removed** |

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
- `R/03d_qc_policy.R` — **QC outlier policy frozen** (user decision 2026-10-18): all **200 flagged samples are retained**; flags are annotation only. `output/tables/qc_outlier_policy.csv` (per-set counts + decision) and `output/tables/qc_flagged_samples_retained.csv` (one row per flagged sample) are the methods-section record. A sensitivity analysis can drop them later without recomputing anything.
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
| QC outlier flag | median sample–sample correlation < `median − 3 × MAD` (flag only, no exclusion) |

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

**Leakage accounting (CORRECTED 2026-10-18):** `val_glom_array` and `val_ti_array` have **0** participants overlapping discovery (ERCB is a different cohort). For the RNA-seq series the question cannot be answered: their titles are local sample indices (`Sample_1`…`Sample_199` in GSE182380 — exactly 1..199) and the series matrices have **no participant field**, so any key-based overlap with the array participant ids is an identifier collision. Measured collisions are small and non-informative (MN 11, FSGS 9, MCD 0, control 6), which is consistent with chance. `R/03b_build_cohort_sets.R` now writes `role = "independent RNA-seq cohort"` and `MN_overlapping_discovery = NA` for both cross sets, and `99_verify_phase_a.R` asserts that.

> **Amended claim for the manuscript:** the discovery array, the two RNA-seq cohorts and ERCB are three *independent* sample sources; the earlier statement that RNA-seq measured "the same people twice" was wrong.

**QC result (`output/tables/qc_summary.csv`):**

| Set | PC1 % | PC2 % | min median-corr | flagged (retained) |
|---|---|---|---|---|
| `disc_glom_array` | 12.0 | 9.0 | 0.458 | 3 |
| `disc_ti_array` | 17.5 | 11.4 | 0.281 | 29 |
| `val_glom_array` | 59.8 | 8.4 | 0.459 | 71 |
| `val_ti_array` | 68.5 | 10.3 | 0.528 | 64 |
| `cross_glom_rnaseq` | 14.5 | 7.3 | 0.188 | 7 |
| `cross_ti_rnaseq` | 22.7 | 7.3 | 0.183 | 26 |
| `pla2r_lcm` | 18.2 | 10.8 | 0.247 | 0 |

**Critical validation result:** the HuGene-2.1-ST gene-level RMA reimplementation was checked against the authors' own Brainarray matrix (`GSE200828_processed_data.txt.gz`) — **Pearson r = 0.985–0.987 per sample, overall 0.9859**, across all 25,582 genes.

### 2.4 Phase B — Batch integration, DEGs, meta-analysis, functional profiling ✅ (2026-10-18)

Every script follows the standing R conventions: blank line between single-line statements, plots printed to the device then saved as **PDF + JPEG**, **no grid lines**, and every figure prefixed with its generating step id.

#### Step 6 — `R/04a_combat_integrate.R`

- Patient de-duplication left exactly one multi-batch discovery matrix: `disc_ti_array` = **GSE133288 (244) + GSE200818 (25)**, both HuGene-2.1-ST (GPL19983) but normalised as separate series → `sva::ComBat(batch = accession, mod = ~group)` on both the Entrez- and symbol-keyed matrices.
- `disc_glom_array` is a single series → **no ComBat**, recorded as such (not silently skipped).

| Metric (disc_ti_array) | Before | After |
|---|---|---|
| Batch silhouette (top-2000 genes, Euclidean) | 0.097 | **−0.033** |
| Batch R² on PC1–5 | 0.168 | **0.0003** |
| Disease silhouette | −0.035 | −0.033 |
| PC1 variance | 16.9 % | 17.5 % |

Batch structure is removed while the biological covariate is untouched. Outputs: `work/combat/disc_ti_array_expr_combat.rds`, `..._expr_combat_symbol.rds`, `output/tables/combat_summary.csv`, `output/qc/04a_disc_ti_array_pca_{before,after,before_after}.{pdf,jpeg}`.

#### Step 7 — `R/04b_deg_limma.R`

- Group vocabulary: `Control` (LD + TN), `FSGS`, `MCD`, `Other` (IgAN, SLE, DN, HT, RPGN, TMD, FSGS_MCD, NEPTUNE "Other"); `UNKNOWN` samples dropped.
- Gene filter per cohort: detected in ≥ 25 % of samples (> 5 log2 intensity for arrays, > 1 log2 CPM for RNA-seq) **and** complete (no NA) **and** non-constant.
- `lmFit(~0 + group)` → `contrasts.fit` → `eBayes(trend = TRUE, robust = TRUE)`. Contrasts: MN vs Control, FSGS, MCD, Other, plus a binary MN-vs-nonMN model. Thresholds **|log2FC| > 1 and BH < 0.05**.
- Discovery TI uses the ComBat matrix; the un-adjusted matrix is refitted for a sensitivity table.

| Set | Tested genes | MN vs Ctrl | vs FSGS | vs MCD | vs Other | vs nonMN |
|---|---|---|---|---|---|---|
| `disc_glom_array` | 7,882 | **904** | 6 | 32 | 17 | 14 |
| `disc_ti_array` | 7,616 | **323** | 0 | 14 | 22 | 1 |
| `val_glom_array` | 11,180 | **2,202** | 1,011 | 289 | 210 | 393 |
| `val_ti_array` | 9,493 | **1,802** | 0 | 69 | 48 | 89 |
| `cross_glom_rnaseq` | 13,487 | **1,537** | 134 | 687 | 57 | 190 |
| `cross_ti_rnaseq` | 11,863 | n/a (0 donors) | 2 | 7 | 0 | 0 |
| `pla2r_lcm` | 19,591 | **2,048** | n/a | 1,521 | n/a | 1,467 |

**ComBat sensitivity:** logFC Pearson **0.983–0.999** (Spearman 0.968–0.998), DEG Jaccard 0.63–0.87. ComBat changes 323 vs 272 DEGs for MN vs Control — it is doing real work but does not overturn the biology. Tables: `output/tables/deg/<set>__<contrast>.csv` (32 files), `deg_summary.csv`, `deg_all_significant.csv` (15,096 rows), `deg_combat_sensitivity.csv`.

#### Step 8 — `R/04c_meta_analysis.R`

Tiered so that independence is preserved:

| Tier | Glom | TI |
|---|---|---|
| **Primary meta** (discovery + cross-platform replication) | `disc_glom_array` + `cross_glom_rnaseq` | `disc_ti_array` + `cross_ti_rnaseq` |
| **Validation (held out of the meta)** | `val_glom_array` (ERCB) | `val_ti_array` (ERCB) |
| **Sensitivity meta** (ERCB added) | 3 cohorts | 2–3 cohorts |

Random effects via `metafor::rma(method = "REML")`, with a **DerSimonian–Laird fallback** for genes where Fisher scoring does not converge (724 such genes at k = 3 were recovered this way). BH correction on the pooled P. **meta-DEG = |pooled log2FC| > 1, pooled BH < 0.05, same direction in every contributing cohort.**

| Compartment | Contrast | k | Pooled genes | meta-DEGs (up/down) | median I² |
|---|---|---|---|---|---|
| Glom | MN vs Control | 2 | 7,522 | **584** (353/231) | 0.0 |
| Glom | MN vs FSGS | 2 | 7,522 | **18** (16/2) | 50.4 |
| Glom | MN vs MCD | 2 | 7,522 | **60** (46/14) | 57.1 |
| Glom | MN vs Other | 2 | 7,522 | **16** (13/3) | 0.0 |
| Glom | MN vs nonMN | 2 | 7,522 | **18** (18/0) | 75.6 |
| TI | MN vs Control | **1** | 7,616 | **323** (242/81) | n/a |
| TI | MN vs FSGS | 2 | 6,438 | 2 (0/2) | 0.0 |
| TI | MN vs MCD | 2 | 6,438 | 4 (3/1) | 0.0 |
| TI | MN vs Other | 2 | 6,438 | 7 (1/6) | 0.0 |
| TI | MN vs nonMN | 2 | 6,438 | 1 (0/1) | 0.0 |

TI MN-vs-Control can only be k = 1 because GSE182380 (TI RNA-seq) has **no healthy donors** — stated in the table rather than papered over.

**Top glomerular meta-DEGs (MN vs Control):** CCNL1 (−1.53), TET2 (−1.84), FOXK1 (+1.37), HBA2/HBA1 (+3.5/+3.4), LUC7L3 (−1.71), CLK4 (−2.07), TARID (−3.18), ECM1 (+3.34), BGN (+1.51), CD151 (+1.46), GADD45A (+1.49), ALB (−2.92). I² is 0 for most of them, i.e. the array and RNA-seq cohorts agree.

**Independent ERCB validation of every meta-DEG** (`output/tables/meta_validation.csv`):

| Compartment | Contrast | meta-DEGs testable | replicated (BH<0.05, same sign) | % | logFC Pearson (all genes) |
|---|---|---|---|---|---|
| Glom | MN vs Control | 456 | 252 | 55.3 | 0.30 |
| Glom | MN vs FSGS | 16 | 6 | 37.5 | 0.19 |
| Glom | MN vs MCD | 46 | 18 | 39.1 | 0.27 |
| Glom | MN vs Other | 12 | 10 | **83.3** | 0.33 |
| Glom | MN vs nonMN | 17 | 14 | **82.4** | 0.26 |
| TI | MN vs Control | 216 | 70 | 32.4 | 0.24 |
| TI | MN vs FSGS | 2 | 0 | 0 | 0.26 |
| TI | MN vs MCD | 2 | 1 | 50.0 | 0.14 |
| TI | MN vs Other | 7 | 4 | 57.1 | 0.34 |
| TI | MN vs nonMN | 1 | 0 | 0 | 0.13 |

**Honest caveat to carry into the manuscript:** MN-vs-Control meta-DEGs replicate at 55 % (glom) / 32 % (TI) in ERCB, and TI MN-vs-FSGS is essentially null in every cohort (array 0, RNA-seq 2). The *glomerular* compartment carries the discriminating signal; the tubulointerstitial compartment mainly separates MN from healthy donors. Do not describe the TI signature as validated against other glomerulopathies.

**All-cohorts sensitivity:** of the 584 glomerular control meta-DEGs, **207 stay significant** when ERCB is pooled in (269 meta-DEGs at k = 3). Full tables in `output/tables/meta/meta_allcohorts_*.csv`, summary in `meta_sensitivity_allcohorts.csv`.

#### Step 9 — `R/04d_enrichment.R`

- Gene sets: MSigDB **v2025.1.Hs** symbol GMT files in `work/msigdb/` (Hallmark 50, C7 5,219, Reactome 1,787, KEGG-legacy 186). These are reference annotation fetched as GMTs; `msigdbr` 26.1.1 needs a download on first use and was therefore not used at run time.
- **ORA** on the meta-DEGs of every contrast with ≥ 10 meta-DEGs (6 contrasts × GO BP/CC/MF + Reactome + KEGG = **30 tables**), universe = genes actually tested in that meta.
- **GSEA** (`fgsea`, epsilon 0, seed 20261018) on the pooled z-score of all 10 metas, Hallmark and C7 (**20 tables**).

| Result | Value |
|---|---|
| Glom MN vs Control, Hallmark | EMT NES **+2.38** (padj 4e-10), apical junction +2.20, myogenesis +2.23, fatty-acid metabolism −2.08 |
| Glom MN vs nonMN, Hallmark | E2F targets −1.95, MYC targets V1 −1.76, fatty-acid metabolism −1.74, EMT +1.73 |
| TI MN vs Control, Hallmark | TNFA signalling via NFKB −2.30, DNA repair +1.87, androgen response −1.95 |
| TI MN vs FSGS, Hallmark | MYC targets V1 −2.28, IFN-γ response −2.47, bile-acid metabolism +2.24, allograft rejection −2.36 |
| C7 significant sets | 829/5,007 (Glom control), 732/4,998 (TI control), 211–427 in the comparator contrasts |

**Expected-signal audit** (`output/tables/enrichment_expected_signals.csv`, 125 keyword × contrast × collection rows, **52 significant**): ECM/collagen organization is the most robust glomerular signal (GO_CC `extracellular matrix` padj 2.1e-5, GO_BP `extracellular matrix organization` 6.1e-4, GO_MF `extracellular matrix structural constituent` 5.4e-4). Complement, FcγR-mediated phagocytosis and antigen presentation are significant in **some** contrasts (KEGG complement & coagulation padj 0.18 glom control; Reactome complement 0.46) but not universally — they are transcript-level weaker than ECM remodelling. Report exactly that, not a stronger claim.

#### Verification

`R/99_verify_phase_b.R` — **86 assertions, 0 failures, exit 0**: QC policy completeness, ComBat matrix integrity and effect size, 32 DEG tables with threshold compliance, 10 metas with table↔summary agreement and I² bounds, ERCB validation scoping, k = 3 sensitivity, 30 ORA + 20 GSEA runs, expected-signal audit, and PDF/JPEG pairing for every figure (53 DEG/meta + 37 enrichment figures). Phase A still passes 51/51.

### 2.5 Separability gate — can MN be told apart from its mimics? ✅ (2026-10-18)

> ⚠️ **SUPERSEDED — see §2.19.** The ERCB AUCs below (0.866–0.915) came from cohort-internal
> standardisation. With training-set statistics (audit F07) the external range is 0.674–0.837.

Phase B's thresholded DEGs suggested MN-vs-other-GN was thin (18 / 60 / 16 meta-DEGs). The gate tests the claim the paper actually needs, before any WGCNA/LASSO investment: regularised elastic-net classification (alpha 0.5, `glmnet`, lambda by inner 5-fold CV), features reselected **inside every outer fold**, genes z-scored within cohort, trained on `disc_glom_array` and tested blind in ERCB (0 shared participants).

| Check (MN vs FSGS/MCD/Other) | Value |
|---|---|
| Nested CV AUC on discovery (K = 50…all) | **0.915 – 0.932** |
| External AUC in ERCB | **0.866 – 0.915** (K = all: 0.915, DeLong CI 0.853–0.977) |
| Independent RNA-seq cohort AUC | 0.926 – 0.942 |
| Secondary endpoint (incl. healthy donors) | 0.911 — no better than the GN-only endpoint |
| Reverse direction (ERCB → discovery) | 0.849 |

**The gate passes.** The honest reading: the 584-vs-18 meta-DEG counts were a *threshold* artefact — discrimination lives in many genes of moderate effect, exactly what an ML panel exploits. `MN vs all non-MN including donors` (0.911) is no stronger than `MN vs other GN` (0.915), so the signature is not merely a disease-vs-normal axis.

**Robustness (`R/05b_gate_robustness.R`, `output/tables/gate_robustness.csv`)**

| Check | Result |
|---|---|
| Permutation null (20 label shuffles) | median AUC 0.500, 95th pct 0.703 → **no shuffle matched the observed AUC** (empirical p = 1/21 = 0.048) |
| ERCB QC-flagged samples removed (71 of 170) | AUC **0.941** (n = 120) — better, so the signal is not a QC artefact |
| Removing flagged samples in both cohorts | AUC 0.938 |
| Reverse direction | AUC 0.849 |

**Per-diagnosis AUC in ERCB** (`gate_auc_by_comparator.csv`): HT 0.98, TMD 0.98, IgA 0.96, SLE 0.92, RPGN 0.91, MCD 0.90, DN 0.85, **FSGS 0.80**, **FSGS_MCD 0.67**. FSGS is the hard comparator, as expected clinically.

**Caveat:** the multiclass (MN/FSGS/MCD/Other) model is *not* good at naming the correct alternative — balanced accuracy 0.44, MN recall 20/21 but heavy leakage of other classes into "MN". Claim a **MN-detection** classifier with a single decision boundary, not a diagnostic multi-classifier.

### 2.6 Exploratory molecular subtypes — negative on reproducibility ⚠️ (2026-10-18)

`R/06a_subtypes.R`: ConsensusClusterPlus (ward.D2, Euclidean, top 2000 variable genes, 1000 resamples, pItem 0.9 / pFeature 0.8), run **separately per platform** because the cohorts cannot be patient-linked.

| Platform | Stable k (PAC < 0.1) | PAC at best k | Cluster sizes |
|---|---|---|---|
| `disc_glom_array` (51 MN) | **k = 2 and k = 3** | 0.0015 (k=2), 0.035 (k=3) | 27 / 24 at k=2 |
| `cross_glom_rnaseq` (62 MN) | **none** | 0.215 (k=5) | 41 / 10 / 6 / 3 / 2 — unstable |

- **Cross-platform concordance is poor:** subtype-defining genes Jaccard **0.007** (6 shared of 889), maximum centroid correlation **0.50**.
- **The array split is technically confounded:** C1 vs C2 differ in median sample–sample correlation (0.686 vs 0.638, Wilcoxon p = 4e-4) and mean expression (3.52 vs 3.40, p ≈ 0). With no clinical or demographic covariates available, this cannot be separated from biology.
- **ERCB transportability exists but is weak:** the array framework assigns all 21 ERCB MN samples (10 C1 / 11 C2, centroid r reported per sample), yet the RNA-seq framework agrees with it only at ARI 0.17.

> **Verdict for the manuscript:** do **not** publish "molecular subtypes of MN" from these data yet. Report the analysis as exploratory and negative — it is a genuine finding that supports waiting for clinical outcomes (B2). Revisit only with eGFR/UPCR/remission plus, ideally, the PLA2R strata.

### 2.7 scRNA-seq arm — COMPLETE, Phase C delivered ✅ (2026-10-19)

Paper 2 already holds the verified single-cell/spatial accession list; Paper 1 **consumes** it rather than building a second atlas. The Tier-1 archives the user downloaded (3.2 GB) live in the **shared** `MN/Datasets` (not in the project) and `scrna/readiness.csv` reports **READY 3/3**.

**Step 07a / 07a0 — manifest & phenotype** ✅ `scrna/samples.csv` (49 samples: `path`, `prefix`, disease label), `scrna/sample_metadata.tsv` (parsed from the GEO series matrices: the `*_RAW.tar` archives carry no phenotype table). The archives themselves are **not** in the project: `scrna/manifest.tsv` points them at `../Datasets/<GSE>_RAW.tar` (the shared MN collection), and 10x samples are read **in place** from the untarred directory by filename prefix — no per-sample link folders, because this volume is exfat (P-34).

**Step 07b — per-sample QC** ✅ **49/49 samples, 432,378 cells kept**

| Accession | Layer | Samples | Cells kept | Doublets |
|---|---|---|---|---|
| GSE271188 | kidney snRNA | 37 | 407,096 | not applied (single-nucleus) |
| GSE171458 | kidney scRNA | 8 | 9,412 | scDblFinder dbr 0.04 |
| GSE241302 | kidney scRNA | 4 | 15,870 | scDblFinder dbr 0.04 |

Filtering: MAD ±3 on log1p(nCount)/log1p(nFeature) per sample plus a hard percent.mt < 20 % ceiling. The 07b run was the step the previous session left **interrupted** (see §3 B10): the resume guard treats an existing file as done, so the two 0-byte objects had to be deleted before re-running.

**Step 07c — integrated atlas** ✅ 236,107 cells × 2,019 genes (2,000 HVG + 54 markers), Harmony on `accession_donor`, 35 clusters, **16** marker-module cell types

| Cell type | Cells | Cell type | Cells |
|---|---|---|---|
| TAL | 71,400 | Podocyte | 6,617 |
| PT | 57,637 | Mesangium | 5,461 |
| DCT | 24,488 | Fibroblast | 4,387 |
| PEC | 19,809 | Mast | 3,910 |
| CD_PC | 18,608 | T_NK | 3,070 |
| Endothelium | 15,435 | Plasma | 1,881 |
| B | 1,079 | Myeloid | 1,693 |
| Neutrophil | 535 | pDC | 97 |

Cells per disease group: IgAN 73,822 / MN 53,206 / Control 49,349 / MCD 37,101 / DKD 22,629. 16,619 cells (7 %) have a module-score margin < 0.05 (low confidence).

**Design change (user decision 2026-10-19):** `MAX_CELLS` 3,000 → **6,000** and `MAX_GENES` 1,500 → **2,000**. The old caps were forced by the 14 GB laptop; this workstation has 61 GB. Re-running with the same seed reproduced the identical 236,107 cells.

**Step 07d — pseudobulk DE & composition** ✅ 16 pseudobulk matrices (`work/scrna/pseudobulk/<celltype>.rds`, genes × donors), **534,128 DE rows** (13 cell types × 4 contrasts: DKD/IgAN/MCD/MN vs control), limma-voom with donor as the unit of replication.

- Composition: **15 MN vs 10 control donors**; no cell type differs after BH (smallest p_adj 0.218 for CD_PC and Myeloid at raw p 0.02-0.03) — report as a negative.

**Step 07e — bulk ↔ scRNA bridge** ✅

| Product | Result |
|---|---|
| `signature_matrix.csv` | **14,751 genes × 16 cell types** written — aligned by gene name from the 28,396-gene union, then restricted to genes measurable in the bulk cohorts |
| `bulk_deg_celltype_attribution.csv` | 946 meta-DEG rows (top cell type + tau specificity) |
| `bulk_vs_pseudobulk_concordance.csv` | best Pearson **Endothelium 0.293**, TAL 0.274, PEC 0.230, CD_PC 0.224, … Podocyte 0.147; **B −0.031, T_NK −0.064** |

Bulk-vs-pseudobulk concordance is weak everywhere (~55-60 % direction agreement). That is a real limit of comparing whole-tissue bulk logFC with one cell type's pseudobulk logFC, not a bug — say so in the manuscript.

**Step 07f — NNLS deconvolution (Phase C)** ✅ 1,328 samples × 16 cell types, **0 degenerate fits**

- **Primary signature = specificity-filtered** (tau > 0.5 on the CPM scale → 2,450 of 14,751 genes) with CIBERSORT-style quantile normalisation of mixtures together with the signature columns. `deconvolution_fractions.csv`, `deconvolution_tests.csv` (**55 of 321 testable comparisons significant**; 79 all-zero comparisons flagged `testable = FALSE`).
- **Sensitivity = unfiltered all-gene signature** → `deconvolution_fractions_allgenes.csv` + `deconvolution_signature_sensitivity.csv`. The two disagree on the dominant cell type in **1,325 of 1,328 samples** — that instability is the honest headline, and it is why the filtered signature is the primary.
- Biology check that passes: glomerular discovery samples are **Podocyte-dominant** (MN 0.661 vs other-GN 0.609, padj 3.3e-4; MN vs Control median Podocyte 0.64 vs Neutrophil 0.03).
- **Caveats to carry:** pDC (97 atlas cells) and Neutrophil (535) signature columns rest on very few cells and carry shared/ambient signal (the Neutrophil profile has the highest S100A9 and LYZ in the atlas). Treat those two fractions as unreliable. Tubulointerstitial fits are less coherent than glomerular ones (B 0.17, Fibroblast 0.28 in MN).

**Verification:** `Rscript R/99_verify_scrna.R` → **42/42 pass, exit 0** (infrastructure + data assertions, including the deconvolution regression guards and the 10x-prefix/link-folder guards).

**Three deliverables this arm adds to Paper 1:** cell-type attribution of the meta-DEGs and gate features; a kidney-derived signature matrix plus NNLS deconvolution (Phase C without a licence token); pseudobulk MN-vs-control per cell type for panel justification.

**Design constraints still standing:** only 6 MN donors in the reference cohort (pseudobulk DE is hypothesis-grade); snRNA under-represents cytoplasmic transcripts so attribution must be assay-aware; **no patient-level linkage to the bulk cohorts** (B7), so all comparisons are signature-level.

### 2.8 Data consolidation into the shared MN collection ✅ (2026-10-19)

The three MN projects share accessions, so the GEO inputs were lifted **out of Paper 1** into one flat, accession-keyed collection.

| Step | Result |
|---|---|
| Move | `Paper_1/Datasets` (13 extracted accessions, 7.9 GB) → **`/run/media/meng/MY_Data/MN/Datasets`**, now **18 GB / 34 items**: `<GSE>_RAW.tar` (archive) beside `<GSE>/{raw,processed}/` (extraction). No collisions; `mv` within the volume was instant |
| Index | `Datasets/manifest.tsv` (accession → layer, projects, status, contents) + `README.md`/`.docx` (layout, house rules, per-project accession lists) |
| Code | `R/_lib.R` resolves `DATA` → `<MN>/Datasets` (env `PAPER1_DATA` overrides; falls back to a project-local `Datasets/`); the three `00*_download*.sh` helpers were repointed off the dead `/home/yang/Documents/MN/Paper_1` path; `scrna/manifest.tsv` archives → `../Datasets/<GSE>_RAW.tar` |
| Dedup (6 GB) | 3.3 GB: the three scRNA archives existed twice (`scrna/raw/` copies deleted). 3.5 GB: `extracted/<GSE>/samples/<s>/` were **symlinks that a copy onto exfat dereferenced into real files** — 07a now records dir + filename prefix, 07b reads the triple with `ReadMtx(feature.column = 2)`, and the folders were deleted |
| Proof | one snRNA sample re-processed end-to-end reproduced exactly (6,477 → 6,463 cells); `GSE241302 R1_MN` reads 36,601 × 17,213 → 11,047 cells; totals unchanged (49 samples / **432,378 cells**) |
| Guard | `99_verify_scrna.R` grew to **42** assertions: 10x prefix → exactly one barcode/feature/matrix triple, no link folders, tier-1 archives resolve |
| Register | all four verifiers exit 0 (51 / 86 / 27 / 42); 0-byte scan clean; free space 55 GB → 49 GB used |


### 2.9 Phase D — WGCNA, module–trait, reproducibility, preservation, hub genes ✅ (2026-10-19)

Five scripts, all following the house conventions (PDF + JPEG, no grid lines, step-prefixed figures, `nThreads = detectCores() − 1`). **Glomerular compartment first**, because Phase B, the separability gate and the deconvolution all put the signal there.

#### Step 08a — `R/08a_wgcna_network.R` (signed networks)

- Gene space: complete, detected in ≥ 25 % of samples (> 5 log2, the same filter as `04b`), non-constant, then the **top 8,000 by MAD** (single block, TOM ≈ 0.5 GB). The MAD cap is a no-op for both compartments (7,882 and 7,616 genes pass); it is there to keep the block bounded.
- `blockwiseModules(networkType = "signed", TOMType = "signed", corType = "bicor", maxPOutliers = 0.05, minModuleSize = 30, deepSplit = 2, mergeCutHeight = 0.25, pamRespectsDendro = FALSE, saveTOMs = FALSE, randomSeed = 20261019)`.
- Soft power = **smallest power with scale-free R² ≥ 0.85**; the fallback (≥ 0.80, then best available) is written into the table and never silently applied.

| Compartment | Set | Samples | Genes used | Soft power | Scale-free R² | Modules (excl. grey) | Genes assigned | Grey |
|---|---|---|---|---|---|---|---|---|
| Glomerular | `disc_glom_array` | 153 | 7,882 | **18** | 0.872 | **9** | 4,057 (51.5 %) | 3,825 |
| Tubulointerstitial | `disc_ti_array` (ComBat) | 269 | 7,616 | **10** | 0.872 | **7** | 6,466 (84.9 %) | 1,150 |

- **Traits**: MN status + a non-MN-GN indicator + the **14 reliable deconvolution fractions**. `Neutrophil` and `pDC` are excluded from the trait set (B8: 535 and 97 atlas cells, ambient-contaminated), and equally from the attribution evidence line in 08e.
- The 3 (glom) / 29 (TI) QC-flagged samples are **flagged, not dropped** (frozen policy, B4); an outlier-removal refit is reported in `wgcna_outlier_sensitivity.csv`.

#### Step 08b — `R/08b_module_trait.R` (modules vs traits)

Two layers, kept separate on purpose: **primary** = `bicor(eigengene, trait)` with Student p and BH over the whole module × trait matrix; **adjusted** = `limma` `lmFit(MEs ~ MN + Podocyte + Fibroblast + Myeloid)`, which is the guard against a module that tracks MN only because MN samples are more podocyte-rich. Gene-level module significance (mean within-module bicor) is reported per trait.

| Compartment | Module | n genes | MN bicor | BH padj | MN effect after composition adjustment | Preserved |
|---|---|---|---|---|---|---|
| Glom | **pink** | 84 | **+0.553** | 4.7e-12 | **padj 3.6e-6** | yes (Z 21.4) |
| Glom | **red** | 144 | **+0.456** | 7.4e-8 | **padj 4.1e-4** | yes (Z 27.7) |
| Glom | black | 113 | +0.384 | 1.0e-5 | padj 0.090 (n.s.) | yes (Z 29.0) |
| Glom | green | 233 | −0.204 | 0.034 | padj 0.090 (n.s.) | yes (Z 36.5) |
| Glom | turquoise | 2,304 | +0.034 | n.s. | n.s. | yes (Z 84.6) |
| TI | — | — | none | — | — | — |

- **The glomerular MN program is two small modules (pink 84, red 144), and both survive composition adjustment.** The 2,304-gene turquoise module is a generic glomerular/structural program with no MN association — a useful negative control.
- Magenta is strongly **podocyte-negative** (bicor −0.675, padj 1.4e-19): a tubule/stromal axis inside glomerular samples.
- **TI: no module tracks MN** — consistent with Phase B's caveat that the tubulointerstitial compartment mainly separates MN from healthy donors (and there are only 5 TI donors). Do not present a TI MN module.
- Zero-MAD traits (MN itself at 33 % prevalence, plus every tubule-specific fraction that is > 50 % exact zeros in glomeruli) fall back to Pearson; the `cor_method` column records which pairs did.

#### Step 08c — `R/08c_wgcna_stability.R` (is the structure reproducible at all?)

Three diagnosis-stratified split halves plus the flagged-sample removal, each refit **re-selecting its own soft power** (a full-cohort power over-sparsifies a half-size network — see P-38). Scored two ways: gene-membership Jaccard and eigengene bicor against the best-matching module of the other fit.

| Compartment | Split Jaccard | Split eigengene bicor | Flagged-removal Jaccard | Flagged-removal eigengene bicor | Verdict |
|---|---|---|---|---|---|
| Glomerular | 0.292 | **0.961** | 0.552 | **0.994** | axes reproducible, membership not |
| Tubulointerstitial | 0.316 | **0.919** | 0.553 | **0.983** | reproducible |

The honest reading: **the co-expression axes are highly reproducible, but which gene lands in which module is partly a resampling artefact** (typical for WGCNA at n ≈ 150). That is exactly why a hub gene is only trusted after independent ERCB preservation, and why Phase E must not select on membership alone.

#### Step 08d — `R/08d_module_preservation.R` (ERCB preservation)

`modulePreservation` with our own module assignment as the reference, 100 permutations for the primary tests and 50 for the secondary cross-platform one, on the 5,911 / 5,626 / 7,574 genes measurable in both. Verdict = **Zsummary.qual > 10 AND qsummary.qual < 0.05** (the `.qual` Z scores run much larger than the legacy `Zsummary.pres`, so a bare threshold would call everything preserved).

| Test | Modules | Preserved | Median Z | Max Z |
|---|---|---|---|---|
| Glomerular → `val_glom_array` (ERCB) | 11 (9 + grey + gold) | **9/9 real modules** | 29.0 | 84.6 |
| TI → `val_ti_array` (ERCB) | 9 (7 + grey + gold) | **7/7 real modules** | 26.6 | 61.7 |
| Glomerular → `cross_glom_rnaseq` (secondary, cross-platform) | 11 | 9 | 26.6 | 90.1 |

- **Grey and gold fail in every test** (Z −27…+1.2, q 0.08–0.80) — the test discriminates rather than rubber-stamping.
- The cross-platform result is reported but **never used to qualify a module**; only the ERCB primary tests gate the hub stage.

#### Step 08e — `R/08e_hub_genes.R` (hubs and the Phase E candidate space)

- Hub = **|kME| > 0.8 inside an ERCB-preserved module**: **708 hub genes** (209 glomerular across 9 modules, 499 TI across 7). Hubs of unpreserved modules are excluded rather than quietly carried forward.
- Candidate space = union of four independent evidence lines — hub, meta-DEG (best contrast per gene), elastic-net gate feature at a **reduced** K (B9: the gate's `K = "all"` row is the whole 13,379-gene model and does not count), and cell-type attribution with **tau > 0.5 to a reliable cell type**.
- `phaseE_candidate_space.csv` (3,279 genes with provenance), **`phaseE_candidate_core.csv` (515 genes with ≥ 2 evidence lines)**. Five genes carry all four lines.
- Top candidates: **ECM1** (pink hub, kME 0.890, MN_vs_Control logFC +3.34, podocyte tau 0.84, gate K = 100), IGFBP2, ATP10A, MXRA8, and TI **C1R** (fibroblast).
- **Caveat to carry:** the glomerular **blue** module is a proximal-tubule module (hubs DDC, FBP1, SLC22A6, ACE2, ALDOB, SLC27A2) — PT transcripts in glomerular samples reflect tubular contamination, and blue's MN correlation is null (−0.060). Treat it as a composition axis, not an MN program, when Phase E ranks candidates.
- `R/99_verify_phase_d.R` — **117 assertions, exit 0**.

### 2.10 Module functional identity — what each module actually is ✅ (2026-10-19)

`R/08f_module_enrichment.R` runs ORA on every non-grey module (universe = that compartment's network genes) against GO BP/CC/MF, Hallmark, Reactome and KEGG-legacy, then assigns each module a programme label from its strongest term. 16 modules labelled, 70,117 module × collection × term rows. The script is **resumable** (it reuses `08f_module_enrichment_all.csv` unless `--force` is passed) because the ORA pass takes ~12 min.

| Compartment | Module | n | Programme label (strongest term) | Supporting terms | MN bicor (adj padj) |
|---|---|---|---|---|---|
| Glom | **pink** | 84 | cell surface (GO_CC 2.9e-05) | calcium ion binding, external encapsulating structure, morphogenesis | **+0.553 (3.6e-06)** |
| Glom | **red** | 144 | metal ion transmembrane transporter activity (GO_MF 6.0e-03) | lysosome, glycosaminoglycan degradation, N-glycosylation | **+0.456 (4.1e-04)** |
| Glom | black | 113 | anchoring junction (GO_CC 4.8e-04) | **Reactome ECM organization**, EMT, Z disc, actin binding | +0.384 (0.090) |
| Glom | green | 233 | transmembrane signaling receptor activity (2.2e-17) | blood-vessel development | −0.204 (0.090) |
| Glom | magenta | 81 | inflammatory response (2.1e-14) | defense response, immune regulation | −0.075 (n.s.) |
| Glom | blue | 405 | organic acid metabolic process (2.4e-55) | oxoacid/carboxylic acid metabolism | −0.060 (n.s.) |
| Glom | brown / yellow / turquoise | 384 / 309 / 2,304 | mitochondrial inner membrane / RNP complex / tRNA processing | — | n.s. |
| TI | turquoise / blue / brown / yellow / green / red / black | 1,829 / 1,743 / 1,279 / 829 / 445 / 194 / 147 | translation / chromatin organization / small-molecule catabolism / signaling receptor binding / vesicle-mediated transport / cell adhesion / carboxylic acid metabolism | — | all n.s. |

**Reading:** the two MN modules are (i) a **podocyte surface / calcium-signalling / matrix-interface** programme (pink: cell surface, calcium binding, external encapsulating structure) and (ii) a **lysosomal / glycosylation / ion-transport** programme (red). The third MN-tilted module (black) is **ECM organization + EMT**, and it loses significance once cell composition is adjusted — an ECM-remodelling readout rather than a cell-intrinsic MN signature. Blue is the proximal-tubule contamination axis (organic-acid metabolism, MN-null), consistent with the 08e caveat. Magenta is an inflammatory infiltrate module (podocyte-negative).

### 2.11 Phase E — the locked 5-gene biomarker panel ✅ (2026-10-19)

> ⚠️ **SUPERSEDED — see §2.19.** The locked panel is now FOUR genes (ICAM1 dropped) with ERCB AUC
> 0.816 (0.721–0.910). The numbers below are the pre-audit record.

`R/09a_panel_selection.R` → `R/09b_panel_validation.R`. Selection runs **inside** the cross-validation folds (10 repeats × 5 folds), so the reported AUC is not selection-optimistic.

**Candidate space:** 353 glomerular core genes → **287 after an ERCB measurability constraint** (a panel gene must be present in ≥ 95 % of ERCB samples; 68 genes dropped, including **GPC6**, which is NA in 125/196 ERCB samples because that release is a two-platform union). Endpoint: MN vs other GN, controls excluded (147 discovery samples: 51 MN / 96 other).

| Selector | Mean nested-CV AUC | Median genes selected per fold |
|---|---|---|
| LASSO (λ.1se) | 0.847 | 11.5 |
| SVM-RFE (linear, size by inner CV) | **0.859** | 20 |
| Boruta (Confirmed, 100 runs) | 0.807 | 24 |

**Consensus rule:** selected by ≥ 2 methods at frequency ≥ 0.60 (frequency measured over all 50 folds).

| Gene | LASSO | SVM-RFE | Boruta | Methods | Module | OR (95 % CI) | p | Meta-DEG | Cell type (τ) |
|---|---|---|---|---|---|---|---|---|---|
| **SLAMF8** | 0.86 | 0.82 | 1.00 | 3 | — | **3.26 (1.33–7.97)** | 0.0097 | MN vs nonMN +1.02 | B (0.76) |
| **TRPC6** | 0.84 | 0.48 | 1.00 | 2 | — | 1.90 (0.56–6.45) | 0.30 | MN vs nonMN +1.17 | **Podocyte (0.83)** |
| **MYOZ2** | 0.62 | 0.46 | 1.00 | 2 | **pink** | **3.97 (1.30–12.14)** | 0.016 | MN vs MCD +2.21 | — |
| **ICAM1** | 0.70 | 0.20 | 0.64 | 2 | — | **0.39 (0.22–0.70)** | 0.0017 | MN vs Control +1.10 | Myeloid (0.70) |
| **PPP1R3C** | 0.70 | 0.02 | 0.80 | 2 | — | **2.67 (1.28–5.55)** | 0.0087 | MN vs Control +1.35 | Mesangium (0.44) |

**Locked-model performance** (logistic regression, genes z-scored with discovery mean/sd, no tuning on validation data):

| Cohort | n (MN/other) | AUC | 95 % CI (DeLong) | Sens / Spec at Youden |
|---|---|---|---|---|
| Discovery (apparent) | 147 (51/96) | 0.935 | 0.894–0.976 | 0.882 / 0.885 |
| **ERCB external** | 170 (21/149) | **0.881** | **0.812–0.951** | 0.905 / 0.758 |
| **RNA-seq independent** | 274 (62/212) | **0.894** | **0.856–0.932** | 0.887 / 0.750 |
| ERCB incl. donors (secondary) | 196 (21/175) | 0.897 | 0.835–0.959 | 0.905 / 0.794 |
| *Reference: gate elastic net (13,379 genes)* | 170 (21/149) | *0.915* | *0.853–0.977* | — |

**Permutation nulls (200 each):** ERCB label shuffles → null median 0.502, max 0.660, **empirical p = 0.005**. Discovery-label refits (refit the panel on shuffled labels, test in ERCB) → median 0.516, 95th pct 0.871, max 0.958, **p = 0.055**. The second null is the stricter test and it is **borderline** — with only 21 ERCB MN cases a random 5-gene model occasionally transfers well. Report both, and do not claim the strict refit null is significant at 0.05.

> **Honest caveats to carry into the manuscript**
> 1. **The panel discriminates but is not calibrated.** Calibration-in-the-large is +0.63 in ERCB and −0.22 in RNA-seq, with calibration slopes 6–45; at a fixed 0.5 cut, ERCB sens/spec is 1.00/0.22 and RNA-seq is 0.00/1.00. Cohort-specific recalibration is required before any clinical probability use. Present it as a *discrimination* biomarker.
> 2. **ICAM1 flips sign**: it is MN-upregulated in the meta-analysis (logFC +1.10) but has a negative model coefficient (OR 0.39) — a suppression/collinearity effect among the five genes. Say so rather than describing it as protective.
> 3. **TRPC6 is not individually significant** (p = 0.30) yet is retained because all three selectors and the podocyte identity support it; the panel is a joint model.
> 4. The 5-gene panel is 2,676× smaller than the gate model for 0.034 less AUC (0.881 vs 0.915) — that trade is the paper's Phase E claim, not superiority.

### 2.12 Panel robustness, transfer and recalibration ✅ (2026-10-19)

> ⚠️ **SUPERSEDED — see §2.19** (TI 0.449/0.539; LCM 0.929; other-CDF 0.848; Brier 0.845 → 0.096).

`R/09c_panel_robustness.R` runs the five checks a reviewer will ask for, each on data already in hand.

| Test | n (MN) | AUC | 95 % CI | Reading |
|---|---|---|---|---|
| ERCB glomerular (reference) | 170 (21) | 0.881 | 0.812–0.951 | as in §2.11 |
| **ERCB tubulointerstitial** | 169 (18) | **0.509** | 0.367–0.651 | **chance — the panel is compartment-specific** |
| **RNA-seq tubulointerstitial** | 199 (59) | **0.556** | 0.469–0.642 | chance again (SLAMF8 absent, scored at the mean) |
| **LCM glomeruli, GSE216841 (third cohort)** | 26 (12) | **0.941** | 0.823–1.000 | independent sample source, small |
| **GSE99339 (ERCB patients, other Brainarray CDF)** | 187 (21) | **0.903** | 0.838–0.968 | **probe-set independence**: same cohort, different annotation |

**The tubulointerstitial null is the strongest specificity result in the paper.** A glomerular signature that scored 0.88 in glomeruli drops to chance in the interstitium of the same patients — consistent with §2.9 (no TI module tracks MN) and with the Phase B caveat that TI separates disease from donors, not MN from its mimics.

**Recalibration (fitted on one half of ERCB, scored on the held-out half):**

| Version | Calibration-in-the-large | Slope | Brier | AUC |
|---|---|---|---|---|
| before | +0.644 | 44.8 | 0.533 | 0.900 |
| **after** | **−0.004** | 5.45 | **0.084** | 0.900 |

A two-parameter (intercept + slope) recalibration removes the systematic offset and cuts the Brier score 6.3-fold **without touching discrimination**. So the miscalibration flagged in §2.11 is fixable per platform; the slope stays >1, i.e. predictions remain over-dispersed, so report recalibrated probabilities rather than raw ones.

**Cross-platform training:** fitting the five genes on the RNA-seq cohort and testing on ERCB gives AUC **0.870** with **5/5 coefficient signs** concordant. The coefficients themselves are platform-portable, not just the rankings.

### 2.13 Where the panel genes live — and one red flag ✅ (2026-10-19)

> ⚠️ **Partly superseded — see §2.19.** ICAM1 is no longer in the panel; the MYOZ2 caveat stands.

`R/09d_panel_mechanism.R` maps the five genes onto the kidney atlas, tests them in MN-vs-control pseudobulk, and asks whether a module-level model would do the same job.

| Gene | Atlas max CPM | Top atlas cell type | Bulk detection | Note |
|---|---|---|---|---|
| **TRPC6** | 44.5 | **Podocyte** | 100 % | the mechanistic anchor (podocyte calcium channel) |
| **ICAM1** | 108.9 | Myeloid (also endothelium) | 100 % | endothelial/immune adhesion |
| **PPP1R3C** | 10.8 | Mesangium | 100 % | mesangial |
| SLAMF8 | 1.0 | B | **9.5 % (RNA-seq), 27 % (LCM)** | weakest link — barely detected on two platforms |
| **MYOZ2** | 16.4 | pDC (unreliable column) | 56–100 % | **93.8 % of atlas cells are zero — no resolvable kidney source** |

**MYOZ2 deserves explicit scrutiny** because it is the single strongest gene in ERCB (alone: AUC 0.903) yet absent from kidney parenchyma in the atlas. It is **not** a probe artefact: it is detected in 56–100 % of bulk samples and its MN-vs-other median difference is consistent in all four cohorts (disc 7.12 vs 5.75; ERCB 9.38 vs 7.67; RNA-seq 6.56 vs 4.96; LCM 8.61 vs 4.36). The honest statement is **tissue-level reproducibility with an unresolved cellular source** (vascular smooth muscle and infiltrating cells are not represented in the 16-cell-type atlas), and it is the first gene to check by immunofluorescence.

**Leave-one-gene-out** (refit on discovery, evaluated externally):

| Model | ERCB | RNA-seq |
|---|---|---|
| all 5 genes | 0.881 | 0.894 |
| without SLAMF8 | **0.900** | 0.895 |
| without TRPC6 | 0.895 | 0.881 |
| without MYOZ2 | 0.804 | 0.877 |
| without ICAM1 | 0.816 | 0.900 |
| without PPP1R3C | 0.873 | 0.889 |
| only MYOZ2 | **0.903** | 0.876 |
| only ICAM1 | 0.847 | 0.564 |
| only TRPC6 | 0.735 | 0.828 |
| only PPP1R3C | 0.780 | 0.738 |
| only SLAMF8 | 0.563 | 0.760 |

No single gene is indispensable (every 4-gene model stays above 0.75 in ERCB; the weakest is without MYOZ2). SLAMF8 adds nothing — dropping it slightly improves both cohorts — so it is defensible to report a **4-gene version (TRPC6, MYOZ2, ICAM1, PPP1R3C)** as a footnoted sensitivity model. Do **not** silently redefine the locked panel on validation performance; the 5-gene result stays primary and the 4-gene model is reported as sensitivity.

**Module scores vs the gene panel** (pink + black + red eigengenes, 3 features):

| Cohort | 5-gene panel | 3 module scores |
|---|---|---|
| Discovery | **0.935** | 0.834 |
| ERCB | **0.881** | 0.750 |
| RNA-seq | 0.894 | **0.904** |

The panel is *not* merely a re-expression of the module summary: it beats the module model in both array cohorts and loses slightly on RNA-seq.

**Pseudobulk MN vs control for the panel genes:** ICAM1 is nominally up in PT (+1.50), CD_PC (+1.76), TAL (+1.52) and endothelium (+2.21) at raw p ≈ 0.02 — directionally consistent with an inflammatory-adhesion mechanism — but **no comparison survives BH** with 15 vs 10 donors. Report as hypothesis-grade.

### 2.14 Panel sensitivity: comparators, pipeline independence, technical confounds, four-gene model ✅ (2026-10-19)

> ⚠️ **SUPERSEDED — see §2.19.** FSGS is now 0.632; the "four-gene" sensitivity model is now a
> three-gene model without SLAMF8 (ERCB 0.845) and is comparable, not superior.

`R/09e_panel_sensitivity.R` answers the four questions a reviewer reaches for first.

**A. Per-comparator AUC** (MN versus one comparator at a time, locked panel):

| ERCB comparator | n | AUC | | RNA-seq comparator | n | AUC |
|---|---|---|---|---|---|---|
| FSGS_MCD | 6 | 0.556 | | Other | 29 | 0.857 |
| **FSGS** | 18 | **0.717** | | **FSGS** | 93 | **0.891** |
| DN | 12 | 0.861 | | MCD | 89 | 0.910 |
| RPGN | 22 | 0.872 | | | | |
| SLE | 32 | 0.893 | | | | |
| MCD | 14 | 0.908 | | | | |
| IgAN | 27 | 0.979 | | | | |
| HT / TMD | 15 / 3 | 0.987 / 1.000 | | | | |

**FSGS is the weak comparator (0.72 in ERCB)** — reported, not hidden; it mirrors the gate's per-diagnosis pattern (FSGS 0.80) and the clinical reality that FSGS is the hardest mimic.

**B. Pipeline independence** — the five genes were re-scored in the submitters' own processed matrix (`GSE200828_processed_data.txt.gz`, Entrez-keyed, disease in the column names), which never touches our RMA/Brainarray path:

| Gene | submitters' matrix | our RMA | Δ |
|---|---|---|---|
| MYOZ2 | 0.865 | 0.866 | −0.001 |
| TRPC6 | 0.857 | 0.863 | −0.006 |
| SLAMF8 | 0.800 | 0.808 | −0.008 |
| PPP1R3C | 0.796 | 0.786 | +0.010 |
| ICAM1 | 0.319 | 0.277 | +0.041 |

The single-gene separation reproduces to within 0.04, so **the panel is not an artefact of our preprocessing** — the decisive check for MYOZ2 (see §2.13).

**C. Technical confound** — panel score versus expression PCs:

| Cohort | R²(score ~ PC1–5) | AUC unadjusted | AUC adjusted for PC1–5 |
|---|---|---|---|
| Discovery | 0.343 | 0.935 | 0.944 |
| ERCB | 0.315 | 0.881 | 0.978 |

The score is not driven by global expression structure, and adjusting does not remove the signal. **Caveat:** the adjusted AUC is fitted on the same cohort it is evaluated in, so treat 0.978 as a diagnostic of confounding (absence of it), not as a performance claim.

**D. Four-gene sensitivity model** (drop SLAMF8 — detected in only 9.5 % of RNA-seq samples):

| Cohort | AUC | 95 % CI | calibration-in-large |
|---|---|---|---|
| Discovery | 0.927 | 0.879–0.974 | 0.00 |
| **ERCB external** | **0.900** | 0.830–0.970 | +0.339 |
| RNA-seq | 0.895 | 0.854–0.935 | −0.213 |
| LCM GSE216841 | 0.982 | 0.943–1.000 | +0.174 |

Permutation *p* = 0.005 in ERCB. The four-gene model is **non-inferior** to the locked five-gene panel in every validation cohort and needs no gene that two platforms fail to measure. The five-gene panel remains the primary (pre-specified) result; the four-gene model is the deployable fallback and is reported as such.

### 2.15 PLA2R strata — B12 resolved by re-parsing ✅ (2026-10-19)

> ⚠️ **SUPERSEDED — see §2.19.** Now 0.976 / 0.881 vs MCD and 0.833 (0.554–1.000) within MN — the
> within-MN contrast is uninterpretable at 6 vs 6.

B12 recorded PLA2R status as unavailable. It is in fact present in the `GSE216841` series matrix as a per-sample characteristic (`phenotype: PLA2R_Positive` / `PLA2R_Negative`), split across characteristic lines; `01a` did not harvest it. `R/09f_pla2r_strata.R` parses it (12/12 MN samples: **6 positive, 6 negative**), writes `metadata/GSE216841_pla2r.tsv` for reuse, and tests the panel by stratum:

| Comparison | n (MN/comparator) | AUC | 95 % CI |
|---|---|---|---|
| MN vs MCD (all) | 26 (12/14) | 0.941 | 0.823–1.000 |
| **MN PLA2R+ vs MCD** | 20 (6/14) | **0.941** | 0.822–1.000 |
| **MN PLA2R− vs MCD** | 20 (6/14) | **0.941** | 0.822–1.000 |
| within MN: PLA2R+ vs PLA2R− | 12 (6/6) | 0.528 | 0.152–0.903 |

Reading: **the panel detects MN irrespective of antigen status** (both strata 0.94, each with 79/84 concordant pairs — verified as a genuine coincidence of rank structure, not a coding slip) and **does not track PLA2R subtype** (0.53 within MN). That is the statement the manuscript needs, and it no longer depends on an external table.

### 2.16 Deliverables produced ✅ (2026-10-19)

| Deliverable | Where | Notes |
|---|---|---|
| **Figures 1–6 + S1–S3** | `output/figures/` (main), `output/figures/supplementary/` | Each plate is a **vector PDF** composited from the source panel PDFs via LaTeX (`tools/make_figures.py`), then Ghostscript-rewritten so **every font is embedded**, plus a 300 dpi JPEG preview. 9/9 built, 0 unembedded fonts |
| **Manuscript draft (EBioMedicine format)** | `Manuscript/Paper1_manuscript.md` + `.docx` | Restructured 2026-10-20: **Summary** (Background/Methods/Findings/Interpretation/Funding, 370 words) → Keywords → **Research in context** (230 words) → Introduction → Methods (incl. Ethics, Role of funders) → Results → Discussion → Limitations → Conclusion → Contributors (CRediT) → Data sharing statement → Declaration of interests → Acknowledgements → References → Figure legends. Main text **3,752 words**; 5,780 words including references, summary and legends. `.docx` uses `reference/tnr_fangsong_redcode.docx` (TNR + 仿宋 + red code) |
| **Reference list** | `Manuscript/references.tsv` (44 entries) | Every citation verified against Europe PMC `resultType=core` or Crossref — authors, title, journal abbreviation, year, volume(issue), pages and DOI; no entry typed from memory. Manuscript carries **44 numbered citations in strict order of first appearance** (checked by `check_cites` logic reproduced in the handoff) |
| **Cover letter** | `Manuscript/Cover_letter.md` + `.docx` | 519 words, addressed to *EBioMedicine*: title, design, the four external AUCs, compartment specificity, PLA2R strata, the two negative results, fit to the journal, declarations. Author list and reviewer suggestions are marked placeholders |
| **IF validation plan** | `Manuscript/IF_validation_plan.md` + `.docx` | Pre-specified success criteria per protein, sample sizes, multiplex protocol, lineage anchors, blinded quantification, five named pitfalls (autofluorescence, podocyte loss with chronic damage, MYOZ2 muscle contamination) |
| **Data-request instructions** | `Manuscript/DATA_REQUESTS.md` + `.docx` | Exact route for B2 (NIDDK Central Repository R4R, DOI 10.58020/f2ae-n094; DACC NEPTUNE-STUDY@umich.edu), B7 (subject→sample mapping via DACC), a paste-ready project description, and the variable list to request |
| **Reproducibility package** | `environment/renv.lock` (1,669 packages), `environment/packages_*.csv`, `sessionInfo.txt`, `README.md` | `renv.lock` is a record, not an active renv project (the analysis keeps resolving its own library in `R/_lib.R`); locally built Brainarray CDFs are marked `Source = Local` and rebuilt by `00d`/`02a` |

### 2.17 P1 manuscript assembly — journal chosen, references verified, EBioMedicine format ✅ (2026-10-20)

Session started by re-reading this handoff and re-checking the recorded state: only
`logs/phase_{a,b,downstream,scrna,d,e}.out` were newer than the document, and they report
**51 / 86 / 27 / 42 / 117 / 102, all exit 0** — i.e. no work had happened after the last handoff.
No 0-byte outputs outside `.reasonix/` task locks.

**One environment drift was found and fixed.** The volume is now mounted at
`/run/media/**yang**/MY_Data` (it was `/run/media/meng`). `R/_lib.R` is path-independent so nothing
broke, but 8 paths in this document and a hard-coded `ROOT=` in three `R/00*_download*.sh` helpers
were stale. All corrected; the helpers now resolve
`ROOT="${PAPER1_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"` and pass `bash -n`.

**Decisions taken (user, 2026-10-20).** Target journal = ***EBioMedicine***; the author list is
**"to be confirmed later"**, so the manuscript carries `[further authors to be confirmed before
submission]` and CRediT for Meng only.

**Format was derived from a live article, not from the publisher PDF** (P-50). The Lancet/Elsevier
author-information PDFs return **403 + Cloudflare challenge**; the authoritative substitute is an
open-access article in the journal — `PMC13657036` `fullTextXML` from Europe PMC gave the real section
order (Summary → Keywords → Research in context → … → Contributors → Data sharing statement →
Declaration of interests → Acknowledgements) and the real reference style. Its abstract is **329
words**, which set the Summary target.

**References: 44 entries, machine-verified.** `Manuscript/references.tsv` holds authors, title, journal
abbreviation, year, volume, issue, pages, DOI and PMID for every entry; nothing was typed from memory.
Each was fetched from Europe PMC `resultType=core` (or Crossref for non-PubMed items: Benjamini &
Hochberg 1995, Boruta JSS, glmnet JSS). **Four title searches returned the wrong paper** and had to be
re-fetched by DOI — NELL-1 (returned case reports), Lake kidney atlas (title search empty; correct DOI
`10.1038/s41586-023-05769-3`), Liberzon Hallmark, Menon JCI Insight — see P-51.

**Citations inserted and numbered.** 23 anchored insertions (`[@key]` placeholders) covering antigen
discovery, guideline diagnosis, treatment response, NEPTUNE/ERCB redundancy, every software method, the
FSGS comparator and the atlas; a generator then assigned numbers in order of first appearance and
rewrote the reference section. The insertion script **asserted every anchor occurred exactly once**
before writing (the P-27 lesson). Final check: 0 unresolved placeholders, 24 citation groups, 44
distinct references, first-appearance order exactly 1..44, list contiguous.

**Structure rewritten for *EBioMedicine*.** Abstract → **Summary** (Background/Methods/Findings/
Interpretation/**Funding**, 370 words); visible **Keywords** line; **Research in context** panel
(230 words); `Ethics` and `Role of funders` moved into Methods; `Data and code availability` replaced by
**Contributors** (CRediT), **Data sharing statement**, **Declaration of interests**, **Acknowledgements**;
**Figure legends moved after References**. Main text **3,752 words**; 5,780 words rendered including
references, summary and legends.

**One honesty fix worth recording.** The first draft of "Evidence before this study" claimed a PubMed/GEO
search with named terms. No such systematic search was run in this session, so the paragraph was rewritten
to state only what is true (the GEO cohorts plus the literature retrieved for this paper, October 2026).
Do not reintroduce a search-methodology claim that was not executed.

**Cover letter.** `Manuscript/Cover_letter.md` + `.docx` (519 words, addressed to *EBioMedicine*): title,
design, the four external AUCs, compartment specificity, PLA2R strata, the two reported negative results,
fit to the journal, and the declarations. Author list and suggested reviewers are marked placeholders.

**Submission assets checked, not assumed.** 9/9 figure PDFs report **0 non-embedded fonts** (`pdffonts`);
all 9 JPEG previews are **2245 × 2953 px = 300 dpi at 190 mm** (634 dpi at 90 mm). All three
`Manuscript/*.docx` (manuscript 32,085 B, cover letter 12,656 B, handoff) plus `handoff.docx` pass
`unzip -t`; all were built with
`pandoc <in>.md --reference-doc=reference/tnr_fangsong_redcode.docx -o <out>.docx`.

**Housekeeping.** Verifier counts made consistent (`99_verify_phase_e.R` = **102**, not the 55/80 that
§6.5 and §7 still claimed). Scratch scripts (`build_refs.py`, `apply_cites.py`, `check_cites.py`,
`backmatter.py`) were run from `/tmp` and are **deliberately not in the repository**; the authoritative
record of the bibliography is `Manuscript/references.tsv`.

**Not done (deliberately), and why:** the Summary was not trimmed below 370 words (each further cut costs
a reported number); the exact EBioMedicine word/figure caps were **not** asserted because the publisher
PDF is unreachable — mark them "to confirm". The 4-gene panel remains a sensitivity model; nothing in this
session touched any analysis result, and no verifier was re-run because no analysis code changed.

### 2.18 Full audit of scripts, parameters, tables and figures — 37 findings, fixes applied, re-run IN PROGRESS ⏳ (2026-10-20)

The user asked for a strict, expert revision of **every** script, parameter, table and figure. The
work was done as an audit-then-fix cycle, and every suspicion was verified against the real data
before it was recorded (two suspicions were disproved). The full finding log, with evidence, is in
**`AUDIT_2026-10-20.md`** (196 lines; section 7, the post-run numbers, is still to be filled).

**Scope actually covered:** 47 R scripts read end-to-end, 106 result CSVs checked by script, 223
figure files checked for pairing/naming/fonts/grid lines, plus visual inspection of the plates.

**Findings:** 37 total — 5 high, 14 medium, 18 low/clarifying. The five high-severity items:

| ID | Where | What was wrong | Evidence |
|---|---|---|---|
| F01 | `09a` | Nested CV was not nested: a fresh `caret_folds()` was drawn inside every (repeat, fold) task, so the five "folds" of a repeat overlapped (mean 19.3 %, max 31 %) and 36/147 samples were scored twice | `/tmp/audit_folds.R` |
| F02 | `04c` | The DerSimonian–Laird fallback was wired to `tryCatch(error=)`, but metafor reports non-convergence as a **warning**; it never fired and no counter existed. The manuscript's "724 genes recovered by DL" has no supporting output | `/tmp/test_reml.R` |
| F03 | `04c` | I² = 0 at k = 2 is a boundary solution, not evidence of cross-cohort agreement | `/tmp/audit_meta.R` |
| F04 | `07f` | Signature and mixtures were quantile-normalised **together**, flattening the cell-type profiles the τ filter exists to preserve | code review |
| F05 | `07b` | Resume guard still checked only `file.exists()`, i.e. P-28's 0-byte trap was documented but never fixed in code | code review |
| F06 | `00d` | Script did not parse (unbalanced parenthesis) — found by the new parse gate | `Rscript /tmp/parse_all.R` |

**Backup of the pre-audit state:** `R/_pre_audit_backup/` (47 scripts, byte-identical copies).

**Scripts changed (17):** `00d`, `01b`, `02b`, `03d`, `04b`, `04c`, `05a`, `07b`, `07f`, `08a`,
`08b`, `08c`, `08d`, `08e`, `08f`, `09a`, `09b`. All 47 scripts parse (`/tmp/parse_all.R` → 0
failures). The changes are listed per script in `AUDIT_2026-10-20.md` §6.

**Decisions taken with the user (2026-10-20):**
1. **Full revision and recalculation** (option A), i.e. also fix the statistical defects, accept that
   the reported numbers change, and re-run the chain.
2. **Single-cell arm:** re-run `07f` only; do **not** re-run `07b`/`07c` (QC threshold and 6,000-cell
   cap stay as they are and are documented as limitations).
3. **The re-run was stopped by the user at 12:18 and is to STAY stopped for now** (decision taken
   immediately after the stop). Nothing restarts automatically; the resume command is
   `bash tools/rerun_from_08d.sh`, to be run only when the user asks.
4. **Verifier assertions may be updated to the new numbers once the chain is re-run** — the user
   authorised this explicitly, on condition that each changed expectation is documented with the
   reason it changed (e.g. ERCB AUC > 0.8 → the corrected value, citing F07). This authorisation does
   **not** cover relaxing an assertion merely because it failed.

**Re-run status (this is the part a new session must not assume is finished).**
Driver scripts: `tools/rerun_after_audit.sh` (full), `tools/rerun_from_05a.sh`,
`tools/rerun_from_07f.sh` (resume points). Log: `logs/audit_rerun.log`.

| Step | Status | Notes |
|---|---|---|
| 01b, 03d | ✅ exit 0 | new `overlap_all_diseases.csv`; flag reconciliation assertion passes (200 = 200) |
| 04b, 04c, 04d | ✅ exit 0 | meta-DEG counts unchanged (584 / 18 / 60 / 16 / 18 glom; 323 / 2 / 4 / 7 / 1 TI) |
| 05a, 05b, 06a | ✅ exit 0 | **external gate AUC changed** — see below |
| 07f | ✅ exit 0 | deconvolution recomputed |
| 08a, 08b, 08c | ✅ exit 0 | stability unchanged (glom Jaccard 0.292 / ME 0.961; TI 0.316 / 0.919) |
| **08d** | ⛔ **STOPPED by user at 12:18** | the glom primary test completed and wrote its table at 12:15; the ti test, the secondary cross-platform test and the summary table never ran |

> **⚠️ MIXED OUTPUT STATE — do not verify, quote or submit from `output/tables` right now.**
> Regenerated on 2026-10-20 (new pipeline): `deg/`, `meta/`, `meta_summary.csv`,
> `separability_*.csv`, `gate_robustness.csv`, `deconvolution_*.csv`, `wgcna_modules_*`,
> `wgcna_module_trait_*`, `wgcna_stability_*`, `wgcna_network_summary.csv`,
> **`wgcna_preservation_glom.csv` (12:15)**.
> Still from 2026-10-19 (old pipeline): **`wgcna_preservation_ti.csv`,
> `wgcna_preservation_secondary.csv`, `wgcna_preservation_summary.csv`**, `wgcna_hub_genes.csv`,
> `phaseE_candidate_space.csv`, `phaseE_candidate_core.csv`, every `09a_*` … `09f_*` table, and all
> `output/models_panel/` figures.
> Consequences: 08e would mix a new glom preservation verdict with an old ti verdict; 09a reads the old
> candidate space; the six verifiers assert the pre-audit numbers. **The chain is internally
> inconsistent until the resume script finishes.**
| 08e, 08f | pending | 08f must run with `--force` |
| 09a → 09f | pending | panel, validation, robustness, mechanism, sensitivity, PLA2R |
| six verifiers | pending | expectations still hold the OLD numbers; they will need updating |
| figures + manuscript | pending | plates must be rebuilt; the gate panel and Methods text need the new values |

**Failures hit during the re-run, and their fixes (all three are real lessons, see §5):**
1. `04c` died with `replacement has 0 rows` — the new `used_dl` column was added to the k ≥ 2 branch
   but not to the k = 1 early return. Fixed; `used_dl = FALSE` added there.
2. `05a` died with *"No case observation"* — the new training-statistics `std()` dropped the old
   `z[!is.finite(z)] <- 0` guard, so all-NA genes reached glmnet/pROC. Fixed for both `std()` and
   `astd()`.
3. `07f` died in `reshape()` with *"'idvar' must uniquely identify records"* — the new `zero_sum`
   column was picked up by `ct_cols`, and NA groups confused `reshape`. Fixed by excluding
   `zero_sum` and mapping NA group to `"Unknown"` for plotting only.

**Key intermediate results that ALREADY differ from the manuscript (must be propagated):**

1. **The separability gate's external performance dropped substantially** after cohort-internal
   z-scoring was replaced by training-set statistics (F07):
   ERCB AUC is now **0.674–0.837** across feature sizes (K = 200 lowest, K = 1000 highest;
   K = all 0.765), versus 0.866–0.915 before. The independent RNA-seq cohort stays at 0.904–0.933.
   *(Historical entry: superseded on 2026-11-18 by F62/§2.20 — training statistics alone do not
   transfer across array platforms, so 05a/05b now also harmonise each validation cohort; the current
   range is 0.861–0.922. Do not quote these numbers.)*
   Nested-CV (10 × 5, correct partition) is 0.911–0.926. **This is a correction, not a regression to
   be reverted**: the old number used the validation cohort's own distribution. The manuscript's
   gate and Discussion paragraphs, and the cover letter, still quote the old range.
2. **`n_dl_fallback = 0` for every meta** — direct evidence that F02's "724 genes" claim was wrong.
   Every meta now also carries `pooling` (the TI MN-vs-control contrast is explicitly
   "single cohort (no pooling)", k = 1) and `median_tau2_meta_deg`.
3. **Deconvolution after separate normalisation:** 2,450 of 14,751 genes retained at τ > 0.5; **0**
   degenerate fits; **67 of 316** testable comparisons significant (was 55 of 321); the
   filtered-vs-all-gene signatures still disagree on the dominant cell type in **1,314 of 1,328**
   samples, so that fragility headline does not change.
4. **WGCNA is stable against these changes:** 08a/08b/08c reproduce the module sizes, the MN modules
   and the split-half stability almost exactly (the deconvolution traits changed, the network did not).

### 2.19 The audit re-run COMPLETED and P0 closed — current numbers ✅ (2026-11-17)

**Session start.** `handoff.md` was found at the project root and read in full; nothing had run after it
except `AUDIT_2026-10-20.md` §7 being filled one minute later (so the old "§7 pending" note was stale).
No Rscript was running, no 0-byte outputs existed.

**Environment drift, found before anything else.** The volume is mounted at
`/run/media/**meng**/MY_Data` again (it was `yang` on 2026-10-20), and **all four `tools/rerun_*.sh`
drivers hard-coded `cd /run/media/yang/...`** — the documented resume command would have `cd`, failed
and `exit 1` before starting 08d. Each driver now resolves its root:

```bash
ROOT="${PAPER1_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT" || exit 1
```

**Re-run.** `bash tools/rerun_from_08d.sh` launched from `cwd=/tmp` (proving the path fix): 08d 21 min,
08e, 08f `--force` 6 min, 09a 23 s, then 09b–09f — **every step exit 0**, `ALL STEPS DONE` 18:15:19.
The mixed output state in §2.18 is gone; every 08d-onward table is regenerated.

**The panel changed — this is the headline consequence of the audit.** Correcting the nested CV to
draw one fold partition per repeat (F01/P-53) changed the consensus panel from five genes to **four:
SLAMF8, TRPC6, MYOZ2, PPP1R3C**. ICAM1 no longer reached consensus (LASSO 0.70/Boruta 0.64 but
frequency below the 0.60 rule once the folds were disjoint) and its fitted coefficient had the
opposite sign to its marginal effect, so the sign-flip caveat in §2.11 is retired.

| Quantity | Pre-audit | **Current** |
|---|---|---|
| Consensus panel | SLAMF8, TRPC6, MYOZ2, ICAM1, PPP1R3C | **SLAMF8, TRPC6, MYOZ2, PPP1R3C** |
| Nested-CV mean AUC (LASSO / SVM-RFE / Boruta) | 0.847 / 0.859 / 0.807 | **0.846 / 0.843 / 0.814** |
| Discovery apparent AUC | 0.935 | **0.922 (0.878–0.967)** |
| **ERCB external AUC** | **0.881 (0.812–0.951)** | **0.816 (0.721–0.910)** |
| RNA-seq independent | 0.894 (0.856–0.932) | **0.901 (0.862–0.939)** |
| LCM GSE216841 (third cohort) | 0.941 | **0.929 (0.824–1.000)** |
| Other Brainarray CDF (GSE99339) | 0.903 (0.838–0.968) | **0.848 (0.759–0.938)** |
| ERCB incl. donors (secondary) | 0.897 | **0.837 (0.749–0.925)** |
| ERCB tubulointerstitium | 0.509 | **0.449 (0.310–0.589)** |
| RNA-seq tubulointerstitium | 0.556 | **0.539 (0.447–0.630)** |
| Label-permutation *p* | 0.005 | **0.005** |
| Refit-permutation *p* | 0.055 (borderline) | **0.11 — NOT significant, report as such** |
| Recalibration CIL before → after | +0.644 → −0.004 | **+0.862 → −0.003** |
| Brier before → after | 0.533 → 0.084 | **0.845 → 0.096** |
| Cross-platform coefficient signs | 5/5 | **4/4** (RNA-seq → ERCB AUC 0.850) |
| Sensitivity model without SLAMF8 | 4 genes, ERCB 0.900 | **3 genes, ERCB 0.845 (0.753–0.937), RNA-seq 0.891, LCM 0.941** |
| Leave-one-gene-out, worst (ERCB) | 0.804 without MYOZ2 | **0.738 without MYOZ2**; MYOZ2 alone 0.903 |
| Per-comparator ERCB: FSGS / FSGS_MCD | 0.717 / 0.556 | **0.632 / 0.484** |
| Per-comparator ERCB: TMD / HT / IgAN / SLE / MCD / RPGN / DN | — | **0.968 / 0.959 / 0.940 / 0.872 / 0.878 / 0.714 / 0.722** |
| PLA2R+ vs MCD / PLA2R− vs MCD | 0.941 / 0.941 | **0.976 (0.921–1.000) / 0.881 (0.712–1.000)** |
| Within-MN PLA2R+ vs PLA2R− | 0.528 | **0.833 (0.554–1.000), n = 6 vs 6 → uninterpretable** |
| Gate external ERCB AUC | 0.866–0.915 | ~~0.674–0.837~~ → **0.861–0.922** (§2.20, F62) |
| Gate secondary endpoint (incl. donors) | 0.911 | ~~0.6405~~ → **0.9007** (§2.20, F62) |
| Deconvolution significant comparisons | 55 / 321 | **67 / 316** |
| Signature instability | 1,325 / 1,328 | **1,314 / 1,328** |
| Glomerular MN module adjusted *p* | pink 3.6e-6, red 4.1e-4 | **pink 1.0e-5, red 6.6e-4** (both still significant) |
| ERCB preservation | 9/9 glom, 7/7 TI | **unchanged** (median Z 29.0 / 26.6; grey and gold still fail) |
| Meta-DEGs, all contrasts | §2.4 | **unchanged** (584 / 18 / 60 / 16 / 18 glom; 323 / 2 / 4 / 7 / 1 TI), `n_dl_fallback = 0` everywhere |

**A second deconvolution correction worth carrying:** the podocyte contrast is significant against
the *other glomerulopathies* (discovery 0.637 vs 0.583, adjusted *p* = 3.3 × 10⁻⁴; RNA-seq 0.847 vs
0.730–0.753, *p* = 1.0 × 10⁻⁷…2.9 × 10⁻⁵) but **not** against healthy glomeruli in discovery
(0.637 vs 0.615, adjusted *p* = 0.29). The claim is "podocyte-rich relative to other nephropathies",
not "podocyte loss".

#### Code changes made in this session (all documented in-file)

| File | Change |
|---|---|
| `tools/rerun_{after_audit,from_05a,from_07f,from_08d}.sh` | path-independent `ROOT=…` instead of the hard-coded `yang` mount |
| `R/09d_panel_mechanism.R` | panel-vs-module label and figure title derived from `length(panel)` (was hard-coded "5-gene") |
| `R/09e_panel_sensitivity.R` | "without SLAMF8" model labelled/sized from `setdiff(panel, "SLAMF8")`; outputs renamed `09e_no_slamf8_model.csv` / `09e_no_slamf8_roc.{pdf,jpeg}`; comparator and authors'-matrix axis windows computed from the data (no cropped whiskers) |
| `R/09b_panel_validation.R` | ROC subtitle split onto two lines (it overflowed the panel); Methods comment de-hard-coded |
| `R/09c_panel_robustness.R` | transfer-plot subtitle uses `length(panel)` |
| `R/04c_meta_analysis.R` | forest and volcano titles use `gsub("_", " ")` so contrasts read "MN vs Control" |
| `tools/make_figures.py` | content-driven paper height (Ghostscript bbox, staged via `/tmp`) + fixed 82 mm panels; Fig 5 title built from the locked panel size |
| `R/99_verify_downstream.R`, `R/99_verify_phase_e.R` | expectations re-expressed with a documented reason each (see `AUDIT_2026-10-20.md` §7.8) |

#### Verification state (final sweep 2026-11-17 18:30, repeated after the last edit)

- `99_verify_phase_a` 51/51, `phase_b` 86/86, `downstream` **27/27**, `scrna` 42/42,
  `phase_d` 117/117, `phase_e` **103/103** — every script exits 0 (`TOTAL_NONZERO_EXITS=0`, run
  twice with identical results, before and after the manuscript/handoff edits).
- Syntax gate: `parse()` over every `R/*.R` → **0 failures**; `py_compile tools/make_figures.py` → OK;
  `bash -n` over all four `tools/*.sh` → OK.
- Nine plates rebuilt (`output/figures`, `output/figures/supplementary`): **0 non-embedded fonts**,
  previews 2245 px wide = 300 dpi at 190 mm, trailing empty band **4–10 %** (was 30–35 %).
- `Manuscript/Paper1_manuscript.md` (+`.docx`, 33.6 kB), `Cover_letter.md` (+`.docx`, 12.7 kB) and
  `AUDIT_2026-10-20.md` (+`.docx`, 24.3 kB) all carry the current numbers; each `.docx` is newer than
  its `.md` and its rendered text contains the new markers with no stale ones (the only pre-audit
  values inside `AUDIT_2026-10-20.docx` are the "before" column of its own before/after table); fonts
  TNR + 仿宋 with red `VerbatimChar` code; all archives pass `unzip -t`.
- `handoff.docx` in the project root is **stale (2026-10-20)** and is deliberately not updated.
- Change inspection without a VCS: `handoff.md` was diffed against its pre-edit copy in `/tmp`
  (241 changed lines, all intended), and every edited block was read back.

#### If the panel must be revisited

- The consensus rule lives in `R/09a_panel_selection.R` (`n_methods_freq` = per-method frequency
  ≥ 0.60, keep rows with ≥ 2 methods) and the fold partition in the same file — the *only* thing that
  changed the panel from 5 to 4 genes was drawing the partition once per repeat (F01/P-53).
- Everything downstream of the panel reads `output/tables/09a_panel_consensus.csv`, and the labels /
  figure titles are now derived from its row count, so a different panel size propagates without
  editing wording. Re-run `09a → 09b → 09c → 09d → 09e → 09f`, then `python3 tools/make_figures.py`.

---

### 2.20 F62 — the gate's cross-platform transfer, corrected (2026-11-18) ✅

**How it started.** A P1 manuscript-polish pass, auditing every quoted number against the current
tables. Everything in Phase E matched (panel AUCs/CIs, ORs, both permutation nulls, recalibration,
three-gene model, per-comparator AUCs, PC confound, PLA2R strata, authors'-matrix, leave-one-out,
pseudobulk, module statistics, preservation). Four defects did not:

| # | Where | Was | Now |
|---|---|---|---|
| 1 | Methods, *Software and reproducibility* | "80 assertions" | **106** |
| 2 | Methods, *Preprocessing and QC* | "ERCB AUC improved from 0.903 to 0.941" | **0.914 → 0.939** |
| 3 | Results, gate paragraph | permutation q95 0.703, p = 0.048, ERCB 0.674–0.837, FSGS 0.80 | see §2.20 numbers below |
| 4 | Discussion ¶4 | "external AUC 0.88 — 0.90" | **0.86–0.92** |

**Root cause: `R/05b_gate_robustness.R` still z-scored each cohort against its own mean/sd.** The
2026-10-20 audit (F07) fixed that in `05a` only; `05b` kept a private `zscore()`, so the permutation
null, the QC-removal AUC, the reverse direction and every per-diagnosis AUC were never re-derived.

**Fixing 05b alone made things worse, which is the finding.** With training statistics and no other
change, the reference ERCB AUC fell to 0.6504, the permutation *p* rose to 0.095 (not significant),
and the per-diagnosis AUCs inverted: RPGN **0.05**, FSGS_MCD 0.01–0.14, FSGS 0.40, DN 0.48 while
HT/IgA/SLE/TMD stayed 0.95–0.99. Scratch diagnostics (`/tmp/audit_05b_*.R`, not in the repo) showed
why:

- 12 RNG streams × the K = all model: ERCB AUC **0.638–0.679** (median 0.645) — so 05a's 0.765 for
  the same configuration was a lucky draw, and the *inversions are stable*, not noise;
- the model keeps only 35–198 non-zero coefficients out of 19,168;
- **|z| after standardising ERCB with discovery μ/σ: median 3.98, 99.9 % 22.1, max 35.4** — versus
  median 0.64 in discovery. Discovery is HuGene-2.1-ST, ERCB is HG-U133: training statistics do not
  transfer across platforms, the linear predictor saturates, and the "transfer" measures platform
  shift.

**Fix (user's decision, 2026-11-18): harmonise, then re-derive.** Both gate scripts now map each
validation cohort onto the discovery reference distribution by rank-preserving quantile mapping
(`quantile_map()`, gene by gene, no outcome information) *before* applying training μ/σ. In `05b` the
scaling now happens inside `run_fit()`, so the reverse-direction test is fitted on its own training
set, and a reconciliation line compares 05b's reference AUC against 05a's K = all row.

**Current gate numbers (all script-generated):**

| Quantity | F07-only (2026-11-17) | **Harmonised (current)** |
|---|---|---|
| Nested-CV AUC (10 × 5, disjoint folds) | 0.911–0.926 | **0.911–0.926** (unchanged, discovery only) |
| External ERCB AUC across K | 0.674–0.837 | **0.861–0.922** (K = all 0.9217, CI 0.862–0.982) |
| RNA-seq independent AUC | 0.904–0.933 | **0.928–0.939** |
| Secondary endpoint (incl. donors) | 0.6405 | **0.9007** |
| Permutation null median / 95th | 0.500 / 0.703 | **0.500 / 0.674** |
| Empirical *p* | 0.048 | **0.048** |
| QC-flagged removed | 0.941 (transductive) | **0.939** |
| Reverse direction (ERCB → discovery) | 0.849 (transductive) | **0.772** |
| Per-diagnosis AUC range | 0.67–0.98 | **0.71–0.98** (FSGS 0.83, FSGS_MCD 0.71, RPGN 0.95, DN 0.89) |
| 05b vs 05a agreement | — | **Δ = 0.008** (reconciliation line in `logs/gate_robustness.out`) |

**Consequence for the paper's claim, and it is a downgrade, not an upgrade.** The 13,379-gene model
now beats the four-gene panel in ERCB (0.922 versus 0.816). The manuscript therefore says the panel
*approaches without matching* the whole-transcriptome model, and that its value is interpretability
and transportability — the pre-audit "same range as a 13,379-gene model" framing is gone from the
Summary, *Research in context*, Results and the cover letter.

**Panel checked against the same convention (09e section E, new table).** The panel keeps its
original simple path, but harmonisation was run as a documented sensitivity: ERCB AUC 0.8156 →
**0.8169**, RNA-seq 0.9006 → **0.8986**, calibration-in-the-large +0.8588 → **+0.1973** (ERCB) and
−0.2148 → **+0.1055** (RNA-seq), score range 0.75–1.00 → 0.00–0.99. So the panel answer does not
depend on the standardisation choice, while the saturated score range that produced a calibration
slope of ~100 does. `output/tables/09e_harmonisation_check.csv`.

**Code changed this session:** `R/05a_separability_gate.R` (harmonisation for X_ex, X_xp and the
secondary endpoint), `R/05b_gate_robustness.R` (scaling inside `run_fit`, quantile mapping,
reconciliation line), `R/09e_panel_sensitivity.R` (new section E + table), `R/09d_panel_mechanism.R`
(a hard-coded "all 5 genes" label → `length(panel)`), `R/99_verify_downstream.R` (documented reasons
refreshed; bounds *tightened* 0.65 → 0.8 and 0.55 → 0.75 because the transfer is now valid rather
than weakened), `R/99_verify_phase_e.R` (**+3 assertions** for the harmonisation check → 106).

**Verification state:** all six verifiers exit 0 — 51 / 86 / 27 / 42 / 117 / **106**. `parse()` OK on
every edited script. Manuscript `.docx`, cover letter `.docx` and the new `SUBMISSION_CHECKLIST.docx`
rebuilt with `reference/tnr_fangsong_redcode.docx`; the rendered text contains the new values and no
stale ones (a `pdftotext`-equivalent check on the .docx text). No composite figure plate contains a
gate panel, so no plate was rebuilt.

**Two traps recorded as P-62/P-63 (§5).** And one caution for the next session: the fallback of
"just use training statistics" is not a safe default across platforms — check the standardised
|z| distribution (median ≈ 1 expected) before believing any cross-cohort transfer.

> Housekeeping slip, recorded honestly: `handoff.docx` in the project root was regenerated from
> `handoff.md` by mistake during this session (the standing rule is `handoff.md` only). `handoff.md`
> is authoritative; `handoff.docx` may be deleted at any time and must never be treated as a source.

---

### 2.21 BAP audit of scripts, tables and figures (`bioanalyst` v1.1.0) — 2026-11-18 ✅

Ran the governed bioinformatics protocol against the finished project (existing tree, so §11
script-revision mode rather than scaffolding a parallel layout). Deliverables:
**`Docs/BAP_AUDIT_2026-11-18.md`**, **`Docs/Findings_Ledger.csv`** (15 graded claims),
**`Docs/CONSOLIDATED_RESULTS.md`**.

**Verdict:** *reliable with important limitations* — the primary claim is E2 across four independent
sample sources; no Critical defect; one Major metadata defect found and fixed.

| Gate | Verdict | Key evidence |
|---|---|---|
| G1 annotation | PASS after fix | 2,165 samples / 453 MN / **178 unique patients**; no accession is single-disease (no batch×group confounding); no analysis set repeats a sample |
| G2 design | PASS | model, thresholds, multiplicity and discovery/validation split written down in the manuscript + `AUDIT_2026-10-20.md` |
| S4 QC | PASS | 7/7 sets `retain`, 200 flagged retained, register reconciles 200 = 200 |
| Counts integrity | PASS | DGEList objects integer; arrays use limma on log-intensity |

**BAP-01 (Major, fixed).** `sample_metadata.csv`'s `gsm` column doubles as the matrix key and mixes
provenance: **1,623 real `GSM####`**, **473 `Sample_N` local indices** (GSE197307 282, GSE182380
199), **69 submitter ids**. The two RNA-seq series both number 1..n, so **134 labels collide** — the
P-2 trap. Not a science error (participant linkage is already refused for those cohorts,
`MN_overlapping_discovery = NA`), but a reader or a future join could read `Sample_1` as one GEO
sample in two cohorts. Fixed **additively** in `R/01a_build_sample_metadata.R` (`id_type`,
`id_is_globally_unique`); pre-existing columns verified `identical()` after the re-run, and
`99_verify_phase_a.R` gained **3 assertions** (→ **54**).

**Minor, closed:** 3 orphan rasters (`Fig2/5/6 …-2.jpg`) + `tools/__pycache__/` removed → **177 PDF
/ 177 raster**, exact pairing. **Minor, open (needs a decision, not an audit):** `R/07b` calls
`scDblFinder()` without `set.seed()` (doublet calls not byte-reproducible; atlas numbers frozen, so
only affects a re-run), and `07b`/`08f` lack a written `# Parallelism: not needed — …` waiver
(BAP R4).

**No claim-verb violation:** a scan for `drives/regulates/proves/mechanism/causes/demonstrates`
returned 4 hits — two background epidemiology, one "the reimplementation was *validated* against the
submitters' matrix", one inside a reference title.

**No new data downloaded, deliberately.** Nothing a public dataset could add would change a claim:
the two real gaps (clinical variables B2, participant map B7) are **access-controlled**, and
`Manuscript/DATA_REQUESTS.md` already holds the route. More GEO series would add samples, not
variables, and would re-import the redundancy this project spent Phase 1 removing.

**Verification after this audit:** six verifiers exit 0 — **54** / 86 / 42 / 117 / 106 / 27.

### 2.22 Figure set re-laid-out for submission — 2026-11-18 ✅

The manuscript's figure *objects* were already submission-grade, but the figure *set* was not: four of
the nine plates were never cited in the Results (Fig 3, Fig 6, S1, S2), the separability gate — the
paper's central quantitative claim — had **five finished panels and no plate**, and Methods promised
"decision-curve analysis and calibration assessment" whose panels were never shown. Fixed as follows.

| Plate | Change |
|---|---|
| **Fig 3** | now the **separability gate**: (A) ERCB ROC by *K*, (B) AUC vs *K*, (C) permutation null, (D) per-comparator AUC. Functional profiling moved out. |
| **Fig S2** | extended to 6 panels: adds **calibration** (C) and **decision curve** (D), so the Methods claim is now shown |
| **Fig S4** (new) | functional profiling (the demoted panels) |
| **Fig S5** (new) | panel-gene atlas localisation (A) + PLA2R strata (B) — two Results subsections that previously had no figure |
| text | all 11 plates now cited at least once; a new Results sentence reports the pathway result (EMT NES +2.38, adj *p* = 4.1 × 10⁻¹⁰; apical junction +2.20; fatty-acid −2.08) and cites Fig S4 — that analysis had been invisible in Results entirely |

**Two real defects found while verifying, both fixed in `tools/make_figures.py`:**
1. **Four plates were two pages** (Fig 2, 5, 6, S2) — the bbox-driven page-height loop could shrink a
   plate below its row budget, pushing the last row to a second page (and leaving the tell-tale
   `-2.jpg` preview). Added a `_page_count()` guard that grows the page from the row budget until
   LaTeX stops breaking. All 11 plates are now **1 page**.
2. The `pdftoppm` page previews (`<plate>-N.jpg`) were never cleaned up; the build now deletes them.

**Verified:** 11/11 plates built, `pdfinfo` pages = 1 for each, `pdffonts` reports **0 non-embedded
fonts** per plate, JPEG previews 2245 px = 300 dpi at 190 mm, **179 figure PDFs / 179 rasters** with no
strays. Manuscript `.docx` rebuilt (35,383 B, newer than the `.md`). Main text is now **4,254 words**
(Summary 377) — still under any plausible cap but worth watching.

> Also this step: `R/handoff.md` turned out to be a **stale 2026-10-18 copy** of this document sitting
> inside the code directory — a trap for a future session. It was moved to
> `Docs/handoff_2026-10-18_superseded.md` (kept for provenance, name says what it is).

> Counts in this entry are **as of that step**. §2.23 then added `05c` and Fig S6, so the current
> inventory is **12 plates / 181 figure PDFs / 221 tables / 48 scripts**, main text **4,336 words**.

### 2.23 P0 consolidation — gate stability across fold partitions (2026-11-18) ✅

The F62 fix was validated once. This closes that loop: **`R/05c_gate_seed_stability.R`** refits the
locked elastic-net gate under **10 independent RNG streams** at two feature counts and reports the
spread, so "ERCB 0.861–0.922 across K" is a range over *both* feature counts and fold partitions.

| Feature set | ERCB AUC across 10 partitions | Range | SD |
|---|---|---|---|
| *K* = 1000 | 0.9003 – 0.9032 | **0.0029** | 0.0008 |
| *K* = all | 0.9118 – 0.9281 | **0.0163** | 0.005 |

**Interpretation.** The harmonised gate is partition-stable; the broken (unharmonised) path was not
(0.638–0.679 across twelve draws, with RPGN/FSGS_MCD AUCs inverted). One nuance worth carrying: the
number of non-zero coefficients still varies 35–198 across draws, so **the prediction is stable where
the selected gene set is not** — which is exactly why Phase E locked the panel by selection frequency
across folds instead of by one fit.

**Outputs:** `output/tables/05c_gate_seed_stability.csv`, `output/models_gate/05c_gate_seed_stability.{pdf,jpeg}`,
new **Fig S6** plate (single panel), a Results sentence after the gate paragraph, and **5 new verifier
assertions** in `99_verify_downstream.R` (27 → **32**). Word count now 4,336 main text / 377 Summary.

**Note for a re-run:** `05c` depends only on `work/sets/`, so it can be re-run in ~6 minutes without
touching the rest of the pipeline.
### 2.24 Analysis strategy v2 — repositioned around a methodological gap (2026-11-18) ✅

`Paper1_strategy_v2_2026-11-18.md` (+ `.docx`) supersedes the 2026-10-15 plan. **No analysis was
re-run**; this is a strategy document, committed to git (`5a77f69`).

**Field scan (verified online this session):** 2024–2025 MN bioinformatics papers keep publishing
`limma → WGCNA → LASSO → AUC > 0.9` on the same accessions — Ren Fail 2025 (PMID 39806779, discovery
GSE115857+GSE200828, "validation" GSE133288) and Front Immunol 2025 (PMID 41132671, GSE108109) —
with no patient de-duplication, no stated standardisation convention, and in the first case a
"validation" cohort from the same NEPTUNE consortium as the discovery set.

**Repositioning:** four tracks — A = Paper 1 frozen for submission; **B = the flagship (new)**: a
leakage-audited benchmark of *published* MN panels with a 7-rung attrition ladder (R0 as published →
R6 permutation null) and a pre-registered falsifiable hypothesis (median R5 AUC < 0.70); C = MR/pQTL +
IF causal anchor; D = outcome modelling (blocked on B2). Plus gap register **G1–G8**.

**Tool upgrades proposed with adoption tests**, verified this session unless marked 【verify】:
TabPFN v2 (Nature 2025, PMID 39780007) — *with* the TabBench-Bio caveat that penalised logistic
regression is statistically indistinguishable on biomedical tables, so a tie must be reported;
BayesPrism/**InstaPrism** (Bioinformatics 2024, PMID 38970377; BayesPrism top-ranked in the
heterogeneous-pseudobulk benchmark, Genome Biol 2024, PMID 38956606); integrated sc+sn reference for
deconvolution (2026, PMC13106970); **decoupleR** for regulator activity. Each requires a pre-specified
ΔAUC margin before adoption.

**Deliberate boundary:** no new methods enter Paper 1 — it is frozen at its current claims, so the
submission and the repository cannot diverge.

### 2.25 Track B started — the panel attrition ladder, and a falsified hypothesis (2026-11-18) ✅

New code: `R/10b_harmonised_matrices.R` (caches the shared harmonisation so the benchmark and the gate
use ONE implementation — P-63), `R/10c_panel_reproduction.R` (the 7-rung ladder),
`R/99_verify_trackb.R` (**27 assertions, green**: identifier = PMID or DOI, ≥2 genes per panel,
exclusions logged with reasons). Register: `metadata/trackB_panel_register.csv`
(6 panels with PMID/DOI/cohorts/gene lists/reported performance, each flagged for whether it stated
de-duplication and standardisation — all 6: no). Two further candidates are logged as exclusions with
reasons in `metadata/trackB_exclusions.csv`. Outputs: `output/tables/10c_attrition.csv` (42 rows),
`10c_panel_summary.csv`, `output/models_trackb/10c_attrition.{pdf,jpeg}`.

**The pre-registered hypothesis was FALSIFIED, and that is the result.** I predicted median R5 AUC
< 0.70; the 6-panel pilot gave **0.821** (median R1→R5 change −0.041). Published panels DO transfer.
P06 is the sharpest illustration: apparent AUC 0.590 but external 0.813, with RNA-seq at 0.435 — a
weak-looking panel that "succeeds" externally and fails on the other platform, exactly what a random
small gene set does.

| Panel | R1 | R3 | R4 transductive | R5 harmonised | R5 RNA-seq | R6 refit-null p |
|---|---|---|---|---|---|---|
| P01 Huang 2025 (PMID 39806779) | 0.835 | 0.811 | 0.871 | 0.893 | 0.836 | 0.124 |
| P02 Xu 2025 (41132671) | 0.806 | 0.788 | 0.834 | 0.830 | 0.799 | 0.095 |
| P03 Yang 2024 (38846934) | 0.825 | 0.795 | 0.571 | 0.675 | 0.866 | 0.139 |
| P04 Han 2023 (37974210) | 0.712 | 0.672 | 0.713 | 0.672 | 0.790 | 0.269 |
| P05 Han 2024 (38029552) | 0.764 | 0.706 | 0.836 | 0.831 | 0.686 | 0.055 |
| P06 Zhang 2023 (DOI 10.3389/fimmu.2023.1252347) | 0.590 | 0.567 | 0.810 | 0.813 | 0.435 | 0.075 |

**Replacement claim (stronger): 0/6 panels beats its own refit permutation null (p = 0.055–0.269) while
all reach R5 0.67–0.89 — so MN-versus-other is an axis any small gene set captures, and no MN paper
reports that null. Our own four-gene panel (p = 0.11) sits in exactly the same place.** Secondary:
the standardisation convention moves R4→R5 by up to 0.22 in *both* directions (N2 now quantified
across independent panels); within-compartment de-duplication is a no-op (147 samples = 147 patients),
so the field's independence failure is cross-compartment (P01 validates glomerular discovery in the
NEPTUNE TI series). 3/5 panels lose ≥1 published gene in a shared 3-cohort gene space (P04 loses two
non-coding RNAs).

`Paper1_strategy_v2_2026-11-18.md` §3.3 now records the falsification verbatim and §3.5 carries the
pilot table and the revised N1/N2 rows. Next: expand the register to ≥10 panels before submission.

### 2.26 Track B — panel-size frontier: "why four genes?" answered (2026-11-18) ✅

`R/10d_panel_size_frontier.R` → `output/tables/10d_panel_size_frontier.csv` +
`output/models_trackb/10d_panel_size_frontier.{pdf,jpeg}`. Genes ranked by moderated *t* **inside every
fold** (no leakage), K = 2…50, external bars lock each K on full discovery. The locked panel's CV and
ERCB values are **read from `09a_panel_cv_summary.csv` and `09b_panel_external.csv`** with `stopifnot`
guards — not transcribed (a first draft hard-coded them from memory; that was caught and replaced).

| K | 2 | 3 | 4 | 6 | 8 | 12 | 20 | 50 |
|---|---|---|---|---|---|---|---|---|
| internal CV | 0.791 | 0.829 | 0.838 | 0.856 | 0.873 | **0.875** | 0.825 | 0.785 |
| external ERCB | **0.858** | 0.826 | 0.821 | 0.802 | 0.767 | 0.827 | 0.843 | 0.690 |

**Findings.** External discrimination is flat from K = 2 to K = 20 and **collapses at K = 50** (0.690),
while internal CV peaks at K ≈ 12 and then overfits — so four genes sits inside the external plateau and
adding genes buys nothing externally. The three-method locked panel matches a purely statistical top-4
at the same K (CV 0.846 vs 0.838; ERCB 0.816 vs 0.821), i.e. the biological selection did not degrade
the model. This is the figure that answers a reviewer's "why four when your 19k model scores 0.922?".

**Verifier:** `99_verify_trackb.R` extended to **37 assertions**, green — including that the locked-panel
row reproduces 09a/09b exactly, that CV peaks before K = 50, that ERCB degrades from K = 4 to K = 50, and
that both Track B figures exist as PDF+JPEG.

**Environment note for the deferred adoption tests:** `BayesPrism` and `MuSiC` are **installed**;
`decoupleR`, `InstaPrism` and Python `tabpfn` are **not** (checked across all `~/envs/*`). So the
two-method deconvolution (G6) can proceed without provisioning; the TabPFN and decoupleR tests need the
user to install packages first.


---

## 3. Current Blockers (Where we are stuck)

### B1 — `.docx` generation — RESOLVED (2026-10-18)
- `pandoc` is installed at `/usr/bin/pandoc`. A project reference template is built by `tools/make_reference_doc.py` → `reference/tnr_fangsong_redcode.docx` (CJK = 仿宋, Latin = Times New Roman, code = red). Convert with:
  `pandoc <in>.md --reference-doc=reference/tnr_fangsong_redcode.docx -o <out>.docx`

### B3 — Deconvolution — RESOLVED (2026-10-19)
- `R/07e` + `R/07f` build a kidney signature matrix from the scRNA atlas and deconvolve all 7 bulk sets with NNLS — no token, fully reproducible. Primary = specificity-filtered signature; all-gene fit retained as a sensitivity table.
- CIBERSORTx is now an **optional** cross-check only; `xCell`/`MCP-counter`/`quanTIseq` remain absent.
- **Residual caveat (not a blocker):** the pDC and Neutrophil signature columns rest on 97 and 535 atlas cells and carry shared/ambient signal, so those two fractions are unreliable. Fixing that means revisiting the atlas annotation (see B8).

### B2 — Clinical variables not acquired
- **Symptom:** NEPTUNE baseline eGFR, UPCR, and longitudinal remission/relapse status are not in any downloaded file.
- **Impact:** Phase C's eGFR/UPCR correlation half and Phase F step 21 (subtype–outcome Cox models) cannot run.
- **Needs:** NEPTUNE data-sharing application, or per-sample annotations from published supplements. 【To be confirmed with the user which route he prefers】

### B4 — PCA outlier decisions — RESOLVED (2026-10-18)
200 flagged samples across the 7 sets were decided to be **retained, documented only**. See `R/03d_qc_policy.R` and `qc_outlier_policy.csv`. No downstream result excludes them.

### B12 — PLA2R status — RESOLVED (2026-10-19; interpretation updated 2026-11-17)

> Current reading (§2.19): **0.976** PLA2R+ vs MCD, **0.881** PLA2R− vs MCD, **0.833 (0.554–1.000)** within MN on 6 vs 6. The within-MN contrast is *not* interpretable, and the manuscript now makes no claim about PLA2R tracking in either direction.
- **Original symptom:** the harmonised `pla2r_lcm_samples.csv` carried no PLA2R+/− field, so the strata analysis looked blocked.
- **Cause:** the annotation *is* in `GSE216841_series_matrix.txt.gz` (`phenotype: PLA2R_Positive` / `PLA2R_Negative`), but it is spread over several `!Sample_characteristics_ch1` lines and `01a` only harvested disease. The samples are keyed by the submitters' ids (`MEMBR1`, `MCD3`, `N_CTRL1`), not by GSM, so the join runs through `!Sample_description`.
- **Fix:** `R/09f_pla2r_strata.R` parses it into `metadata/GSE216841_pla2r.tsv` (6 positive / 6 negative) and reports the strata. **Current values are 0.976 / 0.881 versus MCD and 0.833 within MN (see §2.19); the 0.941 / 0.528 pair in §2.15 is the pre-audit record.**
- **Lesson:** before declaring an annotation unavailable, print every `!Sample_characteristics_ch1` line — they are per-sample vectors, and a given field may sit on any line.

### B5 — GSE73953 (PBMC arm) not normalised
- **Symptom:** it is an Agilent one-channel array (`GPL4133`, plain `.txt`, no CEL) and was never added to `02b_rma_arrays.R`.
- **Impact:** the optional peripheral cross-check is unavailable. Low priority (8 MN samples).

### B6 — Two series matrices are incomplete
- **Symptom:** `GSE108112_series_matrix.txt.gz.part` (18 MB) and `GSE133288_series_matrix.txt.gz.part` (22 MB) fail `gzip -t`.
- **Impact:** **none** — both cohorts' phenotypes were recovered from CEL filenames (`B1_TI_MN_94`), and GSE108112 is dropped by design anyway.

### B7 — RNA-seq cohorts cannot be linked to participants (data limitation, not a bug)
- **Symptom:** GSE197307/GSE182380 carry no participant id; titles are local sample indices. No age/sex/eGFR/UPCR either.
- **Impact:** the subtype analysis cannot use shared-patient replication, and "the same patients measured twice" is withdrawn. It does **not** block Phase D or E.
- **Fix if needed:** ask the NEPTUNE authors (or the GEO submitter) for a sample→participant map.

### B8 — Atlas annotation is only cell-type level, and small clusters are shaky
- **Symptom:** the 16 cell types come from max module score, not from a reference-based transfer. `pDC` (97 cells) and `Neutrophil` (535 cells) are small, and their profiles are contaminated (the Neutrophil column carries the atlas's highest S100A9 and LYZ).
- **Impact:** those two deconvolution fractions and the corresponding pseudobulk DE rows are unreliable; both are flagged in §2.7 and are now **excluded from the Phase D trait set and from the cell-type-attribution evidence line**, so no candidate gene is credited with a Neutrophil/pDC identity.
- **User decision 2026-10-19:** **do not** re-annotate now — caveat in the manuscript instead. Revisit only if a reviewer objects, or if a labelled kidney reference is added.

### B9 — The gate's "all" feature set is not a compact panel — HANDLED in Phase D (still true for the panel)
- **Symptom:** `separability_selected_genes.csv` stores K = "all" as **13,379 genes** (the full elastic-net model), not 4-10 markers. `panel_celltype_specificity.csv` therefore has 12,518 rows.
- **Impact:** the "panel" attribution table is really a whole-transcriptome attribution. Do not present it as a marker panel; the panel is Phase E's job.
- **Phase D handling:** `08e` counts a gate feature as evidence **only** at a reduced K (≤ 1,000), so K = "all" cannot inflate the candidate space; `99_verify_phase_d.R` asserts the space stays under 5,000 genes. Without this the candidate table was 26,864 rows, i.e. the transcriptome.

### B10 — The 2026-10-19 session was interrupted mid-`07b` (RESOLVED)
- **Symptom:** `work/scrna/qc/GSE241302__GSM7720823_R2_MN.rds` and its `_qc.rds` were **0 bytes**, and `output/tables/scrna_qc.csv` was missing, so `99_verify_scrna.R` exited 1.
- **Cause:** the run died on the largest sample (23,503 cells) — consistent with the 14 GB box, not with the current 61 GB machine.
- **Fix applied:** deleted the two empty files (the resume guard is `file.exists()`, so it would otherwise have treated them as done — see P-28), re-ran `07b`; it then also processed R3_MN and R1_NC, which had never run.

### B13 — `07b` doublet calls are unseeded (open, minor)
- **Symptom:** `R/07b_scrna_qc.R` calls `scDblFinder()` without `set.seed()`, so the doublet assignment is not byte-reproducible.
- **Impact:** none on the current results — the atlas is frozen and every object is saved — but a re-run of `07b` could shift a handful of cells.
- **Not tried/fixed on purpose:** fixing it means re-running `07b` (≈25 min for 49 samples), which is a pipeline re-run, not an audit; needs an explicit request.

### B14 — Two heavy serial steps have no written parallelism waiver (open, minor)
- **Symptom:** `07b` (49 samples, ≈25 min) and `08f` (module ORA, ≈12 min) run serially, and 38 of 48 scripts carry no `# Parallelism: not needed — <reason>` line.
- **Impact:** style/protocol compliance only (BAP R4); no numerical consequence.
- **Not fixed on purpose:** adding waivers is trivial but rewriting settled, verified pipeline steps to parallelise them risks changing summation order for no scientific gain.

---

## 4. Next Steps (What to do next)

**Reading key:** the discovery pipeline is closed — those items are labelled **D1–D4** (delivered). The
**open** items are labelled **P0 → P7** in priority order; P0/P0b are outside computation (wet lab and
an external data request), P1 is writing, and everything below P2 is optional.

### ✅ D1 — Phase D: WGCNA, module identity, preservation, hubs (closed 2026-10-19)
- Delivered: 9 glomerular + 7 TI modules, all preserved in ERCB; pink (84 genes) and red (144 genes) carry the MN program and survive composition adjustment; 708 hubs; 515-gene core candidate set. Details in §2.9.
- Re-run order if a parameter changes: `Rscript R/08a_wgcna_network.R` → `08b` → `08c` → `08d` → `08e`; 08c/08d read the 08a objects, 08e reads all of them. Total ≈ 60 min wall clock, dominated by 08d (≈ 15 min per compartment at 100 permutations).

### ✅ D2 — Phase E: biomarker panel, robustness and mechanism (closed 2026-10-19; re-derived 2026-11-17)

> The panel is now **4 genes** — read §2.19, not the numbers in this entry.
- Delivered: `R/09a_panel_selection.R` + `R/09b_panel_validation.R` → a locked **5-gene panel** (SLAMF8, TRPC6, MYOZ2, ICAM1, PPP1R3C) with **ERCB AUC 0.881 (0.812–0.951)** and RNA-seq AUC 0.894 (0.856–0.932). Details, caveats and the two permutation nulls in §2.11.
- Robustness, transfer, recalibration and mechanism follow in `R/09c_panel_robustness.R` and `R/09d_panel_mechanism.R` (§2.12, §2.13): the panel is compartment-specific (TI AUC 0.51/0.56), replicates in a third cohort (LCM 0.941) and across a different CDF revision (0.903), recalibrates cleanly (Brier 0.533 → 0.084), and transfers coefficients across platforms (5/5 signs).
- Re-run order: `Rscript R/09a_panel_selection.R` → `Rscript R/09b_panel_validation.R` → `Rscript R/09c_panel_robustness.R` → `Rscript R/09d_panel_mechanism.R` → `Rscript R/99_verify_phase_e.R`. ~1 min each; all seeded (20261019).
- If the panel is revisited, keep the **ERCB measurability constraint** in 09a (see P-41) or the model cannot be applied to the validation cohort.
- `caret`, `glmnet`, `e1071`, `randomForest`, `xgboost`, `Boruta` are all installed.

### ✅ D5 — Audit re-run, verifier update, figures and number propagation (closed 2026-11-17)
- The 37-finding audit's pipeline is now the pipeline: 08d–09f re-ran to completion, the six verifiers
  pass (51 / 86 / 27 / 42 / 117 / 103) and every reported number in §2.19 comes from it.
- Panel re-derived: **4 genes** (SLAMF8, TRPC6, MYOZ2, PPP1R3C), ERCB 0.816, RNA-seq 0.901, LCM 0.929,
  other-CDF 0.848; refit-permutation *p* = 0.11 (not significant — state it).
- Re-run order if anything upstream changes: `Rscript R/08d_module_preservation.R` → `08e` → `08f
  --force` → `09a` → `09b` → `09c` → `09d` → `09e` → `09f`, or just `bash tools/rerun_from_08d.sh`.

### P0 — Protein-level validation of the four panel proteins 【owner: user, wet lab — plan written】
- **Resume from:** `Manuscript/IF_validation_plan.md`. Pre-specified success criteria per target, sample sizes (MN 10–15 with a PLA2R+/− split, FSGS 6–8, MCD 4–6, normal 4–5), multiplex protocol with lineage anchors (nephrin/podocin, CD31, PDGFRB/desmin, CD20), blinded per-glomerulus quantification.
- **Order:** MYOZ2 first (its kidney cell source is unresolved and it is the strongest single gene, §2.13), then TRPC6 (podocyte), PPP1R3C (mesangial), SLAMF8 last (barely detected on two platforms).
- **Feeds back into:** whether the paper keeps the four-gene panel or drops to the three-gene version without SLAMF8 (§2.19). ICAM1 is no longer in the panel, so it does **not** need validation for this paper.

### P0b — Submit the B2 / B7 data request 【owner: user, long lead time】
- **Resume from:** `Manuscript/DATA_REQUESTS.md`. Route: NIDDK Central Repository R4R (study DOI `10.58020/f2ae-n094`), DACC contact `NEPTUNE-STUDY@umich.edu`; ask for the subject→sample mapping table in the same request (B7).
- This is the only thing standing between the current biomarker manuscript and a prognostic one — start it early.

### P1 — Manuscript revision and submission 【owner: agent + user — journal chosen, draft formatted】
- **Target journal: *EBioMedicine*** (user decision 2026-10-20). Format was derived from a live open-access
  article's JATS XML, not from the publisher PDF, because Elsevier/Lancet PDFs return **403 behind Cloudflare**
  (see P-50).
- **Done 2026-11-17 (P0 chain):** every quoted number regenerated from the audited pipeline and propagated
  into the manuscript, the cover letter and `AUDIT_2026-10-20.md`; plates rebuilt with the F37 fixes.
- **Done 2026-10-20:** 44 verified references (`Manuscript/references.tsv`, every row with a DOI) and numbered
  inline citations; EBioMedicine structure (Summary with Funding, Keywords, Research in context, Ethics and Role
  of funders in Methods, Contributors/CRediT, Data sharing statement, Declaration of interests, Acknowledgements,
  References, then Figure legends); 519-word cover letter; both `.docx` rebuilt with
  `pandoc --reference-doc=reference/tnr_fangsong_redcode.docx` and archive-checked; figures confirmed
  submission-ready (**9/9 vector PDFs with 0 non-embedded fonts**; previews are 2245 px wide = **300 dpi at
  190 mm**, page heights now content-sized rather than a fixed 250 mm).
- **Still to do (user-owned):** confirm the **author list, affiliations, ORCIDs and CRediT** roles (currently
  "further authors to be confirmed"), confirm the funding statement (now "no specific funding"), supply three
  suggested reviewers, and sign off the declarations.
- **Still to do (agent, on request):** trim the Summary from 370 → ~330 words if the submission system enforces a
  hard cap; confirm the exact EBioMedicine word/figure caps from the author information once it is reachable;
  renumber citations if the author list changes the reference set.
- **Numbers are frozen:** every value is regenerated by scripts and asserted by the verifiers, so a drift fails
  verification first.

### P1b — Optional consolidation analyses — OFFERED AND DEFERRED (user decision 2026-11-18)
Scoped, costed, and explicitly **not run** — the user chose only the P0 item (§2.23). Recorded so a
future session does not re-derive the reasoning:

| Priority | Analysis | Why it would consolidate | Cost |
|---|---|---|---|
| P1 | **Panel-size sweep** (2, 3, 4, 6, 8, 12, 20 genes, nested CV → ERCB AUC) | answers "why four genes when your 19k model scores 0.922?"; turns the panel into a documented parsimony frontier | moderate — `09a_panel_selection_frequency.csv` already holds per-gene frequencies, so it is a re-fit per size, not a re-selection |
| P1 | **MYOZ2 provenance test** (correlate with VSMC/mesangial markers ACTA2/MYH11/TAGLN/PDGFRB across the four cohorts; test in the snRNA subset) | converts "cellular source unresolved" from an admission into a *tested* statement — the obvious reviewer objection | cheap, all data in hand |
| P1 | **Third deconvolution method** (nu-SVR/CIBERSORT-style, or a top-500 specific-gene signature) | the BAP protocol requires two methods for deconvolution claims; the paper currently has two signature variants of one method, and the two disagree on the dominant cell type in 1,314/1,328 samples | moderate |
| P2 | Leave-one-cohort-out meta on the 584 meta-DEGs | shows no single cohort drives the signature | cheap |
| P2 | Subtyping power / ARI bootstrap | a negative result is stronger with "we could have detected ARI ≥ X" | cheap |
| P2 | RNA-seq cohort relatedness check (per-disease composition across GSE197307/GSE182380) | the `Sample_N` collision is settled as metadata, but a reviewer may still ask whether the two series share patients; only B7 can settle it definitively | cheap |

**Not worth doing (stated for the record):** adding more GEO series (adds samples, not the missing
variables, and re-imports the redundancy Phase 1 removed); raising the refit-permutation count from 200
(it tightens the estimate of *p* = 0.11 but cannot fix a null whose median is 0.537 with 21 ERCB cases —
a power limit already reported).

### P2 — Phase F: subtypes 【DEFERRED — analysis done, result negative; blocked on B2 only】
- **Resume only with new data:** stable in the array cohort (k = 2, PAC 0.0015) but absent in RNA-seq and technically confounded; cross-platform Jaccard 0.007.
- **Needs:** eGFR/UPCR/remission (B2). PLA2R no longer blocks anything — it is parsed into `metadata/GSE216841_pla2r.tsv` and showed the panel is antigen-independent (§2.15). A participant map (B7) would still help replication.

### ✅ D3 — Phase C: deconvolution (closed 2026-10-19)
- Only the eGFR/UPCR correlation half remains, and it waits on B2.
- Optional follow-ups: `immunedeconv` (xCell/MCP-counter/quanTIseq) as an orthogonal check; CIBERSORTx is no longer required. Revisit the atlas annotation only if a reviewer objects (B8).

### P4 — Acquire clinical variables 【owner: user — step-by-step route written】
- See B2/B7 and **`Manuscript/DATA_REQUESTS.md`** for the exact request path (NIDDK Central Repository R4R, DOI 10.58020/f2ae-n094; DACC contact NEPTUNE-STUDY@umich.edu), a paste-ready project description and the variable list to request. Blocks the outcome analyses in Phases C and F.
- When the export arrives: `metadata/neptune_clinical.csv` -> append to `wgcna_traits.csv` -> re-run `08b`; update the trait-count assertion in `99_verify_phase_d.R`.

### P5 — Optional: normalise GSE73953 【owner: agent】
- Add an Agilent branch to `R/02b_rma_arrays.R` or a small sibling script reading the `.txt` files in `Datasets/GSE73953/raw/`.

### ✅ D4 — Reproducibility package (closed 2026-10-19)
- `environment/renv.lock` (1,669 packages, R 4.6.1, CRAN + Bioconductor), per-library inventories, `sessionInfo.txt`, and `README.md` with the rebuild recipe. Regenerate with `Rscript tools/make_environment_record.R`.
- The seven locally built CDF packages are `Source = Local` (no repository carries them); they are rebuilt from the `.cdf.gz` files inside the GEO archives by `00d` + `02a`.

### P7 — Wire Paper_2 / Paper_3 to the shared dataset collection 【owner: user, when those projects start】
- Paper_2 (single-cell/spatial) already shares `GSE271188`, `GSE171458`, `GSE241302` with Paper 1 and still needs its 5 tier-2 accessions (`GSE286184`, `GSE283826`, `GSE341102`, `GSE199850`, `GSE193512`).
- Pattern: `DATA <- Sys.getenv("PAPER2_DATA", unset = file.path(dirname(PROJ), "Datasets"))` — the same shape as `R/_lib.R`. Paper_3 (MR/drug targets) needs none of these accessions; it needs eQTL/GWAS summary statistics instead.
- Keep `MN/Datasets/manifest.tsv` updated with the `projects` column whenever an accession gains a second user.

---

## 5. Pitfalls & Anti‑patterns (Never repeat these)

### P-1 — Do not trust the plan's cohort arithmetic
The plan asserted ">500 MN samples" and instructed merging releases that are subsets of one another. **Always compute patient-level overlap across *all* releases before building any merged matrix.** Concretely: 13 releases → 178 patients, and 4 releases are 100 % contained in another.

### P-2 — Never infer cross-cohort patient overlap from sample titles
This mistake was made **twice**, so read both halves.

1. The first estimate ("13 shared patients") compared MN patients within one compartment only. Overlap must be computed over **every disease class, compartment and accession**.
2. The "corrected" figure (50/51 MN patients, 208/282 samples) was *also* wrong: GSE197307/GSE182380 sample titles are **local sample indices** (`Sample_1`…`Sample_199` in GSE182380 — exactly 1..199) and their series matrices carry **only disease + tissue, no participant id**. Intersecting those numbers with NEPTUNE participant ids produces **identifier collisions**, not shared patients.

**Rule:** before claiming any cross-release overlap, prove the identifier exists. Check the series matrix for a participant/id field, check whether the ids are contiguous 1..n (local) or sparse (>n, with gaps), and inspect the per-disease overlap table — a large total overlap with near-zero same-disease overlap is the signature of collision. If no participant field exists, state "patient linkage unavailable" and treat the cohorts as independent.

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
They share samples but not genes. Align on the **union of genes** and fill `NA` where absent (implemented in `R/02b_rma_arrays.R`), then drop incomplete genes before PCA/correlation (implemented in `R/03c_qc.R`). **Phase B consequence:** those `NA`s must also be dropped before limma — a `lmFit` on a matrix with `NA` returns `sigma = NA` for every affected gene and `eBayes(robust = TRUE)` then dies with *"covariate contains NA or infinite values"*. `R/04b_deg_limma.R` therefore filters `rowSums(is.na(ms)) == 0`.

### P-11 — Do not redirect a script's stdout to the same file its `writeLines()` writes
The script's own log write truncates the captured stdout. Use `.out` for the shell redirect and let the script write its `.log`.

### P-12 — Do not assume count-style filenames contain counts
`GSE216841_MNd_ncounts_annot.txt.gz` holds **already-normalised decimals**; it also carries `hgnc_GRCh38p12` and `description` columns beside the gene id, so a naive "first column = id, rest = samples" reader silently treats text columns as samples (this produced `NA`s and 36 "samples" instead of 34).

### P-13 — Do not use raw integer donor ids as patient ids
Living donors are numbered 1..n independently of patients, so `1_Living_donor` collides with patient `1`. Namespace them (`LD-<n>`). Symptom before the fix: controls silently deduplicated away and comparator counts were wrong.

### P-14 — Data acquisition is manual in this project
Per a standing user rule, **never download datasets**. Present a concrete request list and let the user download. (Package installation is *not* covered by that rule and was performed directly; MSigDB GMT files are reference annotation and were fetched as such.)

### P-15 — `topTable(fit, coef = "MN - Control")` does not work on a `~0 + group` fit
The coefficients are group means, not contrasts. Build the contrast with `makeContrasts()` and run `contrasts.fit()` first; then `topTable(coef = 1)`.

### P-16 — An `mclapply` worker that throws returns a `try-error`, not `NULL`
Filtering with `!vapply(res, is.null, ...)` lets it through and `do.call(rbind, …)` then recycles the error message into the result table — this silently produced 724 rows whose "gene" was `Error in rma(...): Fisher scoring algorithm did not converge`. Filter on `is.data.frame(z) && nrow(z) == 1` instead, and wrap the model call in `tryCatch`.

### P-17 — `metafor::rma(method = "REML")` fails for some genes
Fisher scoring does not converge for genes with near-identical SEs and zero heterogeneity. Always wrap it and fall back to `method = "DL"`; report how many genes needed the fallback.

### P-18 — `enrichplot::gseaplot2()` does not accept an `fgsea` result
Neither the `data.table` nor a collapsed `data.frame` works (`object@geneList` dispatch error). Use `fgsea::plotEnrichment(pathways[[set]], stats)` for running-score plots and keep `fgsea` output as a table.

### P-19 — `msigdbr` ≥ 26 downloads its gene sets on first use
On this machine that download stalled. Pre-fetch the MSigDB v2025.1.Hs symbol GMT files into `work/msigdb/` and parse them directly (long format: one row per term–gene pair) — that single object then feeds both `enricher()` and `fgsea()`.

### P-20 — Verification scripts must not count "every PDF in a directory"
`99_verify_phase_a.R` originally asserted `length(list.files(output/qc, "\\.pdf$")) == 7`, which broke the moment Phase B wrote its PCA plates into the same folder. Scope figure assertions to the generating step's filename prefix (`^03c_.*_pca`).

### P-21 — The foreground shell has a 2-minute cap
`Rscript R/04c_meta_analysis.R` and `04d` need 4–6 minutes. Run them with `run_in_background` and read the `.out` log; a foreground call silently times out at 2 min regardless of the timeout requested.

### P-22 — `pROC::roc()` needs a response with the two levels it is told to expect
Passing a multi-class factor plus `levels = c("rest", "MN")` fails with *"No control observation"* / *"No case observation"* — the factor has no `rest` level. Build the binary response first (`factor(ifelse(y == "MN", "MN", "rest"), levels = c("rest", "MN"))`) and then call `roc(direction = "<")`. Note that ERCB labels MN as **`MGN`** in `disease_raw`, so case status must come from the harmonised column.

### P-23 — `apply()` on an `sapply()`-built centroid matrix loses column names
`sapply(1:2, ...)` over an unnamed numeric vector returns an unnamed matrix, so `colnames(cc)` is `NULL` and `colnames(cc)[apply(cc, 1, which.max)]` yields a zero-length vector. Set `colnames(cents) <- paste0("C", cls)` explicitly.

### P-24 — PAC is the stability gate, not silhouette
ConsensusClusterPlus "picked" k = 5 for the RNA-seq cohort by the smallest PAC (0.215) — still far above the 0.1 threshold, i.e. there is **no stable partition at all**. Always report PAC per k and refuse to name a k whose PAC exceeds 0.1; a small silhouette is not evidence of a subtype.

### P-25 — A stable subtype split can be a quality split
The array k = 2 partition (PAC 0.0015) tracks sample–sample correlation (0.686 vs 0.638, p = 4e-4) and mean expression (p ≈ 0). With no clinical covariate available that confound cannot be removed — check technical covariates (median correlation, mean expression, detection rate) against every cluster before interpreting it.

### P-26 — `geom_text(stat = "count")` does not accept `fill`
Building a confusion heatmap straight from `aes(fill = after_stat(count))` fails with *"must only have an x or y aesthetic"*. Precompute `as.data.frame(table(...))` and map `fill`/`label` to `Freq`.

### P-27 — Python `str.replace()` with a missing anchor fails silently
Two insertions into R scripts did nothing because the anchor text had been rephrased earlier in the same session, and the following R run died on a missing object. Always `assert anchor in text` (or check the printed diff) before writing the file back.

### P-28 — A 0-byte output file passes a `file.exists()` resume guard
`07b` skips a sample when `file.exists(obj) && file.exists(qc_row)`. The interrupted run left **both files present and 0 bytes**, so a naive re-run would have printed `[resume]` for that sample and then failed much later, in the merge. **Always check `file.size() > 0` (or `-size 0`) before trusting a resume cache**, and `find <dir> -size 0` after any crash.

### P-29 — Seurat 5 `merge()` keeps one layer per sample, and `AddModuleScore` iterates all of them
`merge()` produces `counts.<sample>` / `data.<sample>` layers, each carrying **only that sample's genes**. `AddModuleScore.StdAssay` loops over `Layers(object, search = "data")`, so a marker missing from a single sample triggers *"The following feature lists do not have enough features present in the object"* and aborts. Call **`JoinLayers()` straight after `merge()`** (it fills 0 for genes a sample never detected). Related trap: `ctrl` is drawn **from every expression bin**, so it must be below the **smallest** bin, not the average one — with a 2,019-gene atlas the default `ctrl = 100 × 24 bins` is impossible (we now compute the bins exactly and set `ctrl = 0.6 × min(bin)`).

### P-30 — `merge()` renames colliding cell barcodes, breaking the link back to the samples
10x barcodes repeat across samples, so Seurat appended `_1`, `_2` … to **every** cell. `atlas_cells.csv` then held `CCGAAGTTCCGT_1` while the per-sample object held `CCGAAGTTCCGT`, and `intersect()` returned **zero** cells — `07d` silently aggregated nothing and reported "no cell type had >= 2 donors per group". Fix: rename before merging with the sample tag as prefix (`RenameCells(obj, paste0(tag, "__", barcode))`), and rebuild the same key in the consumer.

### P-31 — `do.call(cbind, ...)` on named vectors of unequal length aligns by POSITION
The per-cell-type pseudobulk matrices do not share a gene universe (13,868-28,342 genes here, because a cell type is only present in some samples), so `cbind` recycled the shorter vectors and **silently scrambled the signature matrix** — R only warned *"number of rows of result is not a multiple of vector length"*. Always align on `union(names(...))` and fill 0. Symptom to watch for: a downstream table whose numbers shift when an unrelated step changes.

### P-32 — NNLS deconvolution returns corner solutions unless the two sides share a scale
Bulk matrices (`2^` of RMA log2, or `2^` of voom log2 CPM) and a CPM-derived signature are on different scales, so NNLS dumped **weight 1.000 on one cell type for every discovery sample**. Three separate traps here:
1. **Fix = quantile-normalise the mixtures together with the signature columns** (`limma::normalizeQuantiles(cbind(sig, mix))`), then unit-sum both sides. Do **not** normalise to the signature's *mean* profile — the mixture then becomes a permutation of that mean and a single column fits it almost perfectly (RMS 2e-4, still degenerate).
2. **Filter the signature to cell-type-specific genes.** On the full shared gene set every column correlates 0.99 with the grand mean, so the fit is driven by housekeeping genes and every sample reads 66-88 % "Neutrophil". `tau > 0.5` restores podocyte dominance in glomerular samples. Top-N-by-abundance marker selection does **not** work — it selects the same ubiquitous genes.
3. **Compute tau on the CPM scale, before column normalisation.** After `sweep(sig, 2, colSums(sig), "/")` the values are ~1e-5 and `log2(x + 1)` collapses to ~0, which marks 6,342 genes as specific instead of 2,450 and quietly undoes fix 2.

### P-33 — Small API differences that cost a full re-run each
- `[.DGEList` has **no `drop` argument**: `y[keep_g, , drop = FALSE]` errors with *"unused argument"*. Use `y[keep_g, ]`.
- `ncol()` on a Seurat 5 object returns a **double**, so `vapply(sub, ncol, integer(1))` fails; wrap with `as.integer()`.
- `wilcox.test()` on two all-zero groups returns `NaN`, not p = 1. Flag those comparisons as untestable (`testable = FALSE`) instead of letting `NA` propagate into `p.adjust` and the summary line.

### P-34 — This volume is **exfat**: no hardlinks, no symlinks, and `ln -f` DELETES the target first
`/run/media/meng/MY_Data` is formatted exfat (`mount` shows `type exfat`, `fmask=0022`). Consequences that bit us on 2026-10-19:
1. **Hardlinks are impossible** (`ln: Operation not permitted` / *Permission denied*), and so are **symlinks**. So two projects cannot share one physical file — a shared dataset must either be pointed at by path or genuinely copied. Case-insensitive too.
2. **`ln -f existing dest` is destructive here:** `-f` unlinks the destination and only then calls `link()`, which fails on exfat. It silently removed `Paper_1/scrna/raw/GSE171458/GSE171458_RAW.tar`; recovered with `cp -p` from the master. **Never use `ln -f` on this volume.**
3. Tooling that writes atomically via *temp file + link/rename* also fails (Reasonix's `write_file` returns `publish new file … link … operation not permitted`). Write new files with a plain shell heredoc on this volume.
4. `mv` **within** the same exfat filesystem is a cheap rename (the 7.9 GB `Datasets` move was instant), but moving across filesystems would copy.
5. **A copy onto exfat silently dereferences symlinks.** `scrna/extracted/<GSE>/samples/<s>/{barcodes,features,matrix}` were symlinks created by 07a on the previous filesystem; after the project was copied here they were real files with **different inodes**, duplicating 3.5 GB. `du` vs `ls -l` is how it was caught (no `->`, and `extracted/` was 2× the archive size). Fixed on 2026-10-19 by having 07a record dir + prefix and 07b call `ReadMtx`; the `samples/` folders were deleted.

### P-35 — WGCNA 1.74's `modulePreservation` will not take your network object
Six failed calls before this worked. The actual API in 1.74:
1. `referenceNetworks` must be a **numeric index** into `multiData` (it is coerced with `as.numeric`); passing `list(net)` dies with *"'list' object cannot be coerced to type 'double'"*.
2. Your modules go in through **`multiColor = list(<colour vector>, NA)`** — the function re-derives connectivity from `multiData[[1]]$data` and uses the colours you give it.
3. `multiData` is a **2-level** list `list(set = list(data = samples x genes))`. Three-level nesting (set → tissue → data) errors with *"'data' must be of a vector type, was 'NULL'"*.
4. `includeModules` and `verbose` are **not** arguments.
5. The output is the **`.qual` family** (`Zsummary.qual`, `Z.meanAdj.qual`, `Z.meanSignAwareKME.qual`, …), not the legacy `Zsummary.pres`; q-values need `calculateQvalue = TRUE`.
6. **The statistics carry no module labels.** Map them back by matching the reported `moduleSize` against the module sizes recomputed on the shared gene set; the gold module is inserted into the array (here at position 4). Always assert the mapping is unambiguous.

### P-36 — `bicorAndPvalue` returns `$bicor`, not `$cor`
`corAndPvalue(..., method = "bicor")` does not exist (it accepts only pearson/kendall/spearman and dies in `match.arg`). Use `bicorAndPvalue(x, y, use = "p")` and read **`$bicor`** — `$cor` is `NULL`, which then surfaces as *"non-numeric argument to mathematical function"* one line later. Related: **`mad()` is 0 whenever more than half the values are identical**, so a binary trait at 33 % prevalence and every tubule-specific fraction in glomeruli score MAD = 0 and silently fall back to Pearson; record the method per pair instead of claiming bicor throughout.

### P-37 — A module colour vector must be named by GENE, never by sample
`setNames(colors, rownames(dat))` (samples) instead of `colnames(dat)` (genes) makes `moduleEigengenes()` mis-align the colours against the expression columns, and **both** the Jaccard and the eigengene correlations turn to noise (they read 0.00–0.25 instead of 0.55–0.99). Symptom: a stability panel in which nothing is reproducible while module–trait correlations are strong. Check the direction of the matrix before naming anything: `dat` is samples × genes.

### P-38 — A split-half stability test must re-select the soft power inside each half
Re-using the full-cohort power on half the samples over-sparsifies the network (the same threshold cuts weaker correlations), which manufactures instability. `08c` now runs `pickSoftThreshold` per subsample with the identical rule and records `power_ref` / `power_alt`. Report membership (Jaccard) **and** eigengene reproducibility separately: axes can be near-perfect (bicor 0.96) while membership is only partly stable (Jaccard 0.29).

### P-39 — `numericLabels = TRUE` makes `signedKME` return `kME_ME1`, not `kME<colour>`
With numeric module labels the eigengenes are `ME1…MEk` and `signedKME` names its columns `kME_ME1…`. Looking them up as `paste0("kME", labels2colors(colors[i]))` finds nothing, the `if (!col %in% colnames(kme)) next` guard hides it, and **every kME is silently NA** — the module tables look fine until hub selection returns 0 genes. Normalise with `colnames(kme) <- paste0("kME", sub("^ME", "", colnames(mes$eigengenes)))` and `stopifnot(sum(is.na(own_kme[colors != 0])) == 0)`.

### P-40 — `bulk_deg_celltype_attribution.csv` is keyed `meta_<Glom|TI>_<contrast>`
The meta tables carry the bare contrast (`MN_vs_Control`), the attribution table prefixes it (`meta_Glom_MN_vs_Control`), so a naive `%in%` join returns **zero rows silently** and every candidate loses its cell-type evidence. Build the full key (`paste0("meta_", compartment_label, "_", contrast)`) and print the join counts.

### P-41 — A panel gene must be measurable in the cohort you validate in
`GPC6` survived selection (LASSO 0.74, Boruta 0.98) and then turned out to be **NA in 125 of 196** ERCB glomerular samples — the GSE104948 release is a two-platform union, so some genes exist on one platform only. The first validation run silently produced 120 NA scores and pROC reported the unhelpful *"No case observation"*. **Make coverage a selection constraint** (≥ 95 % non-NA in every cohort the model will be applied to) rather than dropping a gene after seeing the validation data, and assert the coverage in the verifier.

### P-42 — Never report sensitivity/specificity at the Youden point of a low-prevalence cohort
With 21 MN among 170 ERCB samples the Youden-optimal threshold came out at 0.96, which looks absurd and hides the real problem: the model is **over-confident** (calibration-in-the-large +0.63 in ERCB, −0.22 in RNA-seq, slopes 6–45). Report sens/spec at a **fixed 0.5** as well, report calibration-in-the-large and the calibration slope, and describe such a model as a discrimination biomarker that needs recalibration — never as a ready risk calculator.

### P-43 — Compute null summaries per null, not once
The first version of the permutation table reused `median(null_ext)` and `quantile(null_ext, 0.95)` for **both** null rows, so the refit null appeared to have a 0.62 95th percentile while actually reaching 0.87. Build each row from its own distribution, and include `null_max` so a borderline p value is visible rather than implied.

### P-44 — `pROC` reports "No case observation" when the predictor is NA, not when the response is empty
A cohort with both classes and a predictor containing NAs gives *"No case observation"* from `roc()`, which sends you looking at the labels. Check `sum(is.na(score))` first: impute deliberately (or fail loudly), then call `roc()`.

### P-45 — Gene-wise summaries on a genes × samples matrix are ROW-wise
`colMeans(ms[genes, ])` and `apply(ms[genes, ], 2, sd)` on a genes × samples matrix return **one value per sample**, not per gene. Feeding those into `sweep(t(ms), 2, ...)` recycles by position, the z-scores become NaN, a `zz[!is.finite(zz)] <- 0` guard turns the whole matrix into zeros — and every downstream AUC comes out as exactly **0.5** with no error. Use `rowMeans()` / `apply(..., 1, sd)`, and treat a suspiciously exact 0.500 as a signal that the design matrix is constant.

### P-46 — GSE99339 is split by platform and codes MN as **MGN**
Its two series matrices (`GPL19109` = Plus2, `GPL19184` = 133A) together hold the cohort; the MN samples are titled `…-Glom-MGN…`, not `-MN…`. A file glob of `GPL1910[94]` matches only the first platform and a regex looking for `-MN[0-9]+` matches nothing, which silently produced "187 samples, **0 MN**" and made the annotation-robustness test impossible. Check the token census of the CEL filenames when a parsed cohort comes back with zero cases.

### P-47 — `merge()` renames colliding columns instead of warning
`merge(diag, atlas, by = "gene")` where both frames have a `note` column returns `note.x` / `note.y`; any later `d$note` is then `NULL` and an assertion like "low-detection genes are flagged" fails for the wrong reason. Give the derived column a distinct name (`atlas_note`) rather than relying on the merge suffix.

### P-48 — GEO characteristics are per-sample vectors spread over several lines
`!Sample_characteristics_ch1` appears once per *field group*, each line carrying one entry per
sample. In `GSE216841` the PLA2R call sits on a different line from tissue and cell type, so
`unlist()`-ing the lines and testing each scalar finds almost nothing (1 of 34 samples). Iterate over
**lines**, compare `length(line)` to the sample count, and assign by position. Related trap: that
cohort is keyed by submitter ids (`MEMBR1`, `MCD3`, `N_CTRL1`) rather than GSM, so the join key has
to come from `!Sample_description` (`MEMBR1.MN.PLA2Rpos` → `MEMBR1`).

### P-49 — Ghostscript cannot read this exfat workspace, and R's `pdf()` leaves Helvetica un-embedded
Two separate traps behind figure submission. `gs` fails with *"undefinedfilename … Permission
denied"* for any path on `/run/media` even as the owning user, while the identical command works on
`/tmp` — so stage inputs and outputs through `/tmp`. Separately, `ggsave(device = "pdf")` writes
Helvetica **without embedding** it, which journals reject; the fix used here is
`gs -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dSubsetFonts=true` over the composited plate,
which embeds every font (verified with `pdffonts`).

### P-50 — Publisher author-guideline PDFs are Cloudflare-blocked; read a live article's JATS XML instead
`https://…thelancet.com/…/ebm-info-for-authors.pdf` and the Elsevier guide-for-authors page both return **403 +
"Just a moment…"** to curl (and a 1.2 MB HTML challenge page when it looks like it worked — always `file` the
download). The reliable source of truth is an open-access article in the target journal: fetch its
`fullTextXML` from **Europe PMC** (`/PMC<id>/fullTextXML`) and read the actual section titles, back-matter order
and reference style. That is how the EBioMedicine structure used here was established, and it is current by
construction. Corollary: do not quote word/figure caps you could not read — mark them "to confirm".

### P-51 — Never hand-type a citation; fetch it by DOI and check the title
A title-string search is not an identifier. In this session four queries returned the **wrong paper**:
`NELL-1 + membranous nephropathy` returned three case reports instead of Sethi 2020; the Lake kidney-atlas
title search returned nothing at all (the DOI is `10.1038/s41586-023-05769-3`, and its title contains
"states and **niches**" — a missing word); Liberzon and Menon likewise needed a second pass. **Practice:**
fetch by identifier (`Europe PMC DOI:"…"` or `https://api.crossref.org/works/<DOI>`) and assert the
returned title matches what you expect. Store authors/journal/year/volume/issue/pages/DOI in a TSV
(`Manuscript/references.tsv`) and generate the list from it — never retype publisher metadata, and never
trust a first search hit.

### P-52 — String-built reference lists accumulate double periods; assert on the output
`authors.rstrip('.')` followed by rejoining with `'et al.'` and then a `'%s. '` template produced
`… , et al.. Membranous nephropathy.` for every multi-author entry. The fix is to strip punctuation
**once at the end** (or use `'et al'` without the period), but the durable lesson is to **assert on the
generated text** — `assert '..' not in generated` caught it immediately. Same family as P-27: when a
script assembles prose, verify the assembled result, not just that the script ran.

### P-53 — A nested CV must draw its folds once per repeat, outside the fold loop
`09a` reseeded inside `one_fold()` and called `caret_folds()` there, so every (repeat, fold) task had
its own random partition. The five "folds" of a repeat overlapped by 19 % on average and 36/147 samples
were scored twice, while the output table still said "10 repeats × 5 folds". **Practice:** generate
`folds[[r]]` once per repeat, pass it into the fold function, and record the structure in the run
info. Verify with a partition check (pairwise test-set overlap must be 0), not by eyeballing the AUC.

### P-54 — `tryCatch(error=)` does not catch a non-converging mixed model
metafor 5.2.1 reports Fisher-scoring non-convergence with a **warning** and still returns a fit; four
boundary cases produced no error at all. A fallback wired only to `error=` is dead code, and any
manuscript sentence that relies on it is unsupported. **Practice:** wrap in `withCallingHandlers`,
match the warning text, record a per-gene flag, and report the count. If a previous run claimed a
fallback count, treat the claim as unverified until a counter exists.

### P-55 — Never quantile-normalise a signature matrix together with its mixtures
`normalizeQuantiles(cbind(sig, mix))` replaces every signature value by a cross-column quantile mean;
with hundreds of samples against 16 signature columns the sample distribution dominates and the
cell-type specificity is flattened — exactly what the τ filter is there to preserve. **Practice:**
normalise the two sides separately, then unit-sum each; keep the all-gene signature as the explicit
sensitivity arm.

### P-56 — Adding a column to a per-gene table means adding it to EVERY branch
`04c`'s new `used_dl` was added to the k ≥ 2 data.frame but not to the k = 1 early `return()`, and the
numeric-coercion loop then failed with `replacement has 0 rows, data has 7616`. **Practice:** after
adding a column, grep every `return(data.frame(` / `data.frame(` in the function, and assert the
column set is identical across branches before the loop that coerces types.

### P-57 — `reshape()`'s `varying` must exclude every new non-measurement column
Adding `zero_sum` to the deconvolution table silently put a logical column into `varying`, and the NA
group from the new UNKNOWN exclusion made `reshape` report *"'idvar' must uniquely identify
records"*. **Practice:** define `ct_cols` from an explicit allow-list (the signature column names),
never as "everything except these five"; map NA groups to a plotting label.

### P-58 — Fixing a transductive standardisation WILL lower external performance
Replacing per-cohort z-scores with training-set statistics dropped the gate's ERCB AUC from
0.87–0.92 to 0.67–0.84 while the RNA-seq cohort stayed at 0.90–0.93. That is the honest number: the
old figure used the validation cohort's own mean/sd. **Never revert to cohort-internal scaling to
recover the higher AUC**, and never present the old range alongside the new pipeline. Say in the
Discussion that array-to-array transfer is the hard direction; cross-platform transfer to RNA-seq is
easier.

### P-59 — A driver script that hard-codes the absolute project path dies on the next mount
The volume appears as `/run/media/{meng,yang}/MY_Data` depending on the session, and all four
`tools/rerun_*.sh` drivers hard-coded one of them, so the documented resume command `exit 1`-ed before
running anything. **Practice:** in every shell driver resolve the root the way `R/_lib.R` does —
`ROOT="${PAPER1_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"; cd "$ROOT"` — and test it by
launching from an unrelated cwd (e.g. `/tmp`); a driver that only works from the project directory is
untested.

### P-60 — A result that is data-driven must not be hard-coded into labels or filenames
The panel size is decided by 09a's consensus rule, but "5-gene" was baked into 09c/09d/09e titles, the
Fig 5 plate title and the filename `09e_four_gene_model.csv`. When the corrected CV produced a
four-gene panel, the tables and figures silently mis-described the model. **Practice:** derive the
label from `length(panel)`, name outputs by role (`09e_no_slamf8_model.csv`) rather than by a count,
and let the compositor read the panel size from `09a_panel_consensus.csv`.

### P-61 — Do not verify a composed figure by eye; measure it
Layout defects (30–35 % of every page blank, a clipped subtitle, cropped CI whiskers) were invisible in
a thumbnail and obvious in numbers. **Practice:** (a) get the content height from
`gs -sDEVICE=bbox` (staged through `/tmp`, because `gs` cannot read this exfat volume — P-49) and size
the page from it; (b) measure the trailing blank band by rasterising and scanning rows for ink;
(c) compare the ink bounding box with the MediaBox to catch text that overflows its panel;
(d) grep the extracted text for the titles you expect (`pdftotext`) to catch `MN_vs_Control`-style
labels. A fixed `ylim` is a bug everywhere — compute axis windows from the data.

### P-62 — "Use training statistics" is not a cross-platform transfer
Removing a transductive z-score is necessary but **not sufficient** when the training and validation
cohorts sit on different measurement platforms. Discovery (HuGene-2.1-ST) against ERCB (HG-U133):
standardising ERCB with discovery μ/σ puts the median gene at **|z| ≈ 3.98** (99.9 % 22.1, max 35.4)
instead of ≈ 0.64, the linear predictor saturates, and with 19k collinear features the model picks
only 35–198 non-zero coefficients and ranks four comparator diagnoses *above* MN (RPGN 0.05,
FSGS/MCD 0.01–0.14). The AUC also becomes K-dependent (0.638–0.679 at K = all across 12 RNG streams).
**Practice:** after any cross-cohort standardisation, check the standardised |z| distribution — a
median near 1 is the sanity condition; 4 means the transform is not a scale transfer. Then harmonise
(rank-preserving quantile mapping onto the reference distribution, label-blind) *before* applying
training statistics. Do not present a lower AUC as "the honest number" before this check: the
transductive path can be optimistic, but the naive correction can be meaningless.

### P-63 — A re-run is not a re-derivation
The 2026-10-20 audit changed the gate's standardisation rule in `05a`, re-ran the whole chain, and
reported every step green — but `05b` carried its own private copy of the z-scoring helper, so its
permutation null, QC-removal AUC, reverse direction and per-diagnosis AUCs silently described the
*old* pipeline while sitting in tables dated the same day as the rest. **Practice:** when a shared
statistical rule changes, grep every script for a private implementation of it (`zscore(`, `std(`,
`scale(`) rather than assuming the rule lives in one place, and prefer a reconciliation print that
compares two scripts' versions of the same quantity (05b now prints its reference AUC against 05a's
K = all row, with a warning above 0.05 AUC).

### P-64 — A height-fitting loop that only measures can shrink a figure onto a second page
The plate compositor sets the PDF page height from a `gs -sDEVICE=bbox` measurement of the rendered
content. For four plates (Fig 2, 5, 6, S2) the measured "ink height" came out far below the row
budget, the page was sized to ~92 mm for two rows of 82 mm panels, and the last row was pushed onto a
**second page** — which journals reject, and which is what the stray `pdftoppm` preview files
(`<plate>-2.jpg`) were actually telling us. **Practice:** after any content-driven page sizing, assert
the result (`pdfinfo` → `Pages: 1`) and regrow the page from the *row budget*, not from the previous
(shrunk) value; delete the `<plate>-N.jpg` preview leftovers so a two-page plate cannot hide.

### P-65 — A column that doubles as an analysis key must record its provenance
`sample_metadata.csv` uses one `gsm` column as the matrix column key and the GEO accession. Only the
array series carry a real accession: 1,623 `GSM####`, but 473 `Sample_N` local indices (GSE197307,
GSE182380) and 69 submitter ids. Because both RNA-seq series number 1..n, **134 labels collide** — a
naive join or a reader concludes the two cohorts share 134 samples, when the analysis already refuses
participant linkage for them. **Practice:** whenever a key column is shared between file-level
provenance and matrix indexing, add `id_type` + `id_is_globally_unique` columns and assert them in the
verifier; never let a local index masquerade as an accession.

### P-66 — Count-bearing prose in `Docs/` goes stale the moment an artefact is added
Three times this session a documentation file kept a true-at-the-time count that the next addition
falsified (177 → 181 figure PDFs, 9 → 12 plates, 47 → 48 scripts, 220 → 221 tables). **Practice:**
re-derive every count in prose from the filesystem (`ls R/*.R | wc -l`, `find output -name '*.pdf' |
wc -l`) when you touch the file, and prefer writing counts as "as of <date>" where they are quoted.

### P-67 — Two "checks" that lied, and the habit that caught both
`pdffonts` output puts the `emb` column at `$(NF-4)`, not `$4` (field 4 is the *encoding*); reading the
wrong field reported "2 non-embedded fonts" on a plate where every font was embedded. Separately, an
`awk -F','` census of a quoted CSV with embedded commas picked the wrong columns and reported a
2/8/5 grade split where the true split was 1/10/4. **Practice:** calibrate any checker against a known
value before reporting it (BAP §11.2), and parse CSVs with a real reader — a checker you have not
calibrated is a random number generator.


---

## 6. Key Resources & Environment Info

### 6.1 Compute
| Item | Value |
|---|---|
| OS / shell | Linux amd64, bash |
| CPU | **20** threads |
| RAM | **61 GB** (~57 GB available) — the old 14 GB constraint is gone |
| Disk | 1.8 TB free on `/run/media/meng/MY_Data` |
| Project footprint | **~33 GB** (GEO data no longer counted here): `work/` 23 GB (of which `work/wgcna/` 57 MB), `scrna/` 3.3 GB, `.Rlib` 4.2 GB, `output/` ~650 MB (of which `output/wgcna/` 47 MB), `logs/` ~55 MB, `reference`+`metadata`+`tools` 6 MB |
| Shared data root | `/run/media/meng/MY_Data/MN/Datasets` — **18 GB** (9.8 GB archives + 7.9 GB extracted), shared with Paper_2/Paper_3 |
| Volume filesystem | **exfat** — no hardlinks, no symlinks, case-insensitive (see P-34) |
| Audit artefacts | `AUDIT_2026-10-20.md` (37 findings), `Docs/BAP_AUDIT_2026-11-18.md` (bioanalyst protocol: gate verdicts + Critical/Major/Minor register), `Docs/Findings_Ledger.csv` (15 graded claims: 10 E2 / 4 E1 / 0 E4), `Docs/CONSOLIDATED_RESULTS.md`, `R/_pre_audit_backup/` (pre-audit scripts), `tools/rerun_*.sh` (4 drivers), `logs/audit_rerun.log` |
| Mount name | the volume is mounted at `/run/media/**meng**/MY_Data` (it has alternated `meng` ↔ `yang` between sessions). Every path in this document, the three `R/00*_download*.sh` helpers **and the four `tools/rerun_*.sh` drivers** are path-independent since 2026-11-17 — see P-59 |
| Foreground shell cap | **2 minutes** — anything longer must run with `run_in_background` (P-21) |

### 6.2 Languages and libraries
- **R 4.6.1 (2026-06-24) "Happy Hop"**; `pandoc` at `/usr/bin/pandoc`.
- **Two libraries are in play** — this is essential to reproduce anything:
  1. Workspace library `<PROJ>/.Rlib` — prepended by `R/_lib.R` (`.libPaths(c(WS_LIB, .libPaths()))`). Holds `makecdfenv`, the 5 Brainarray CDF packages, `oligo`, `oligoClasses`, `affyPLM`, `arrayQualityMetrics`, `WGCNA` 1.74, `ConsensusClusterPlus` 1.76, `caret`, `glmnet` 5.1, `randomForest`, `xgboost`, `Boruta`, `metafor` 5.2.1, `msigdbr`, `rms`, `babelgene`, `hugene21sttranscriptcluster.db`, `hgu133b.db`, `pd.hugene.2.1.st`.
  2. User library `/home/meng/R/x86_64-pc-linux-gnu-library/4.6` — **1,701 packages**. Confirmed present for the scRNA arm: `Seurat` 5.6.0, `scDblFinder` 1.26.9, `harmony` 2.0.5, `muscat` 1.26.0, `miloR` 2.8.1, `nnls` 1.6; plus `limma` 3.68.5, `edgeR`, `DESeq2`, `sva`, `GEOquery`, `preprocessCore`, `affy`, `clusterProfiler`, `org.Hs.eg.db`, `fgsea` 1.38.0, `ComplexHeatmap` 2.28.0, `survival` 3.8.12.
  3. `/usr/local/lib/R/site-library` 135 packages; `/usr/lib/R/site-library` empty; `/usr/lib/R/library` base only.
- `R/_lib.R` resolves the project root from its own path, so the project can be moved between machines/drives without editing scripts (override with `PAPER1_PROJ`). **Scratch scripts outside `R/` must set `PAPER1_PROJ`, or `_lib.R` resolves `PROJ` to `/`.**
- **Not installed:** `ReactomePA` (covered via the MSigDB C2:CP:REACTOME GMT), `immunedeconv`/`xCell`/`quanTIseq`, `msigdbdf`, and `liana` (not published for Bioc 3.23 — cell–cell communication belongs to Paper 2).
- **R package mirrors used:** `https://mirrors.tuna.tsinghua.edu.cn/CRAN/` and `https://mirrors.tuna.tsinghua.edu.cn/bioconductor`.
- **Reproducibility record: `environment/renv.lock` (1,669 packages, R 4.6.1, CRAN + Bioconductor),**
  `environment/sessionInfo.txt`, per-library package inventories and `environment/README.md`
  (rebuild recipe). Generated by `Rscript tools/make_environment_record.R`; see D4 in §4.
- **Version control: git, since 2026-11-18.** `git log` is authoritative; the substantive commit is
  `6f8cebe` (initial import: code, `Docs/`, `Manuscript/`, figure PDFs, cited tables), followed by
  `c7ba326` (ignore `Rplots.pdf`) and small handoff corrections. **326 files**,
  24.5 MB of content; `work/` (23 GB), `.Rlib/` (4.2 GB), `scrna/extracted/` (3.3 GB), `logs/` and five
  bulk CSVs >3 MB are ignored as regenerable — **a checkout alone is not a reproducible run**, it needs
  `Datasets/` + a re-run of the numbered scripts. `core.fileMode = false` and `core.autocrlf = false`
  are set locally because the volume is exfat. `handoff.md` remains the progress carrier: commit it at
  every checkpoint (and still copy `environment/` + `Manuscript/` off-volume for a disk-failure hedge).

### 6.3 Data locations
| Path | Contents |
|---|---|
| **`$DATA` = `/run/media/meng/MY_Data/MN/Datasets`** | **Shared dataset collection** for Paper_1/2/3, 18 GB. Flat, accession-keyed: `<GSE>_RAW.tar` (archive) next to `<GSE>/{raw,processed}/` (extracted). Indexed by `manifest.tsv` (accession → layer, projects, status) and explained in `README.md`/`.docx`. Paper 1 reaches it through `DATA` in `R/_lib.R`; override with `PAPER1_DATA`. |
| `$DATA/<GSE>/raw/` | Extracted CEL files + the Brainarray `.cdf.gz` / mapping files shipped in the tars |
| `$DATA/<GSE>/processed/` | Repaired count matrices (GSE197307/GSE182380/GSE216841) and the authors' `GSE200818/200828_processed_data.txt.gz` |
| `$DATA/<GSE>/*_series_matrix.txt.gz` | GEO metadata (11 of 13 intact; see B6) |
| `metadata/ensg2symbol.tsv` | 19,847-entry Ensembl→symbol map harvested from GSE216841 |
| `work/rma/` | 7 normalised array cohorts (`.rds` + `.txt.gz`) |
| `work/rnaseq/` | 3 normalised RNA-seq cohorts (+ `DGEList` objects for the two count cohorts) |
| `work/sets/` | The 7 cohort sets: `*_expr.rds`, `*_expr_symbol.rds`, `*_samples.csv` |
| `work/combat/` | ComBat-adjusted discovery TI matrices (Entrez + symbol) |
| `work/msigdb/` | MSigDB v2025.1.Hs symbol GMTs (Hallmark, C7, Reactome, KEGG-legacy) |
| `work/cel/` | **16 GB** of decompressed CELs (the single largest reclaimable item) — deletable once `work/rma/` is final |
| `work/rma_cache/` | 1.1 GB chunk caches — deletable once `work/rma/` is final |
| `work/cdf/` | Built CDF source packages (provenance for the `.Rlib` installs) |
| `output/tables/` | 22 result tables at top level |
| `output/tables/deg/` | 32 per-cohort × per-contrast limma tables |
| `output/tables/meta/` | 30 meta tables (10 primary + 10 meta-DEG subsets + 10 all-cohorts sensitivity) |
| `output/tables/enrichment/` | 50 ORA/GSEA result tables |
| `output/qc/` | 7 QC PCA (03c) + 3 ComBat PCA (04a) plates, PDF + JPEG, plus per-sample QC tables |
| `output/deg/` | 53 DEG/meta figures (volcano, forest, heatmap) |
| `output/enrichment/` | 37 enrichment figures (dotplots, GSEA bars, running scores) |
| `output/models_gate/` | 5 gate figures (ROC, AUC vs K, confusion, permutation null, per-diagnosis AUC), PDF + JPEG |
| `output/subtypes/` | 11 subtype figures (consensus matrices, PAC curve, heatmaps, centroid concordance, ERCB assignment), PDF + JPEG |
| `output/tables/separability_*.csv` | gate results: nested CV, external transfer, confusion, selected genes, secondary endpoint |
| `output/tables/gate_robustness.csv`, `gate_auc_by_comparator.csv` | permutation null / QC sensitivity / reverse direction, per-diagnosis AUCs |
| `output/tables/subtype_*.csv` | stability (PAC), clusters, defining genes, signature concordance, QC check, ERCB assignment |
| `scrna/` | scRNA arm: `INTEGRATION_PLAN`, `DATA_REQUEST`, `manifest.tsv` (archives resolve to `../Datasets/`), `readiness.csv` (**READY 3/3**), `samples.csv` (49: `path` + `prefix`), `sample_metadata.tsv`, `metadata/<GSE>_series_matrix.txt.gz`, `raw/` (empty — archives live in the shared folder), `extracted/<GSE>/` (3.3 GB untarred, read in place) |
| `work/scrna/qc/` | 49 per-sample QC'd Seurat objects, all genes (needed for pseudobulk) + 49 `_qc.rds` metric rows + `cells_qc.csv` |
| `work/scrna/atlas.rds` (978 MB), `atlas_cells.csv` (24 MB) | Harmony-integrated atlas: 236,107 cells × 2,019 genes, 35 clusters, 16 cell types |
| `work/scrna/pseudobulk/<celltype>.rds` | 16 donor × gene pseudobulk count matrices + donor tables |
| `output/scrna/` | scRNA figures: 07b sample QC, 07c UMAP (cell type + disease) + marker dotplot, 07d composition, 07e bulk-vs-pseudobulk, 07f deconvolution stacked + boxplots — all PDF + JPEG |
| `output/tables/scrna_qc.csv`, `scrna_celltype_counts.csv`, `scrna_pseudobulk_de.csv` (534,128 rows), `scrna_detection_fraction.csv`, `scrna_composition.csv`, `scrna_composition_tests.csv` | 07b–07d result tables |
| `output/tables/signature_matrix.csv`, `bulk_deg_celltype_attribution.csv`, `panel_celltype_specificity.csv`, `bulk_vs_pseudobulk_concordance.csv` | 07e result tables |
| `output/tables/deconvolution_fractions.csv` (primary), `deconvolution_tests.csv`, **`deconvolution_fractions_allgenes.csv` (sensitivity)**, `deconvolution_signature_sensitivity.csv` | 07f Phase C result tables |
| `work/wgcna/` | 08a–08e objects: `datExpr_<comp>.rds` (network input), `net_<comp>.rds` (network + eigengenes + traits + module table), `me_<comp>.csv`, `traits_<comp>.csv`, `preservation_<comp>.rds`, and `permutedStats-actualModules_partial.RData` (modulePreservation's permutation cache — it writes this into the working directory, so run 08d from the project root and move it here afterwards) |
| `output/tables/wgcna_network_summary.csv`, `wgcna_traits.csv`, `wgcna_softpower_<comp>.csv` | 08a: network parameters, the 422-row trait matrix, scale-free fit grid |
| `output/tables/wgcna_modules_<comp>.csv`, `wgcna_module_sizes_<comp>.csv`, `wgcna_module_stability_<comp>.csv`, `wgcna_outlier_sensitivity.csv` | 08a: gene→module→kME, module sizes, per-module Jaccard vs the outlier-removed refit |
| `output/tables/wgcna_module_trait_<comp>.csv`, `..._adjusted_<comp>.csv`, `wgcna_module_significance_<comp>.csv`, `wgcna_trait_diagnostics_<comp>.csv` | 08b: bicor + BH, composition-adjusted limma fit, gene-level module significance, per-trait MAD/method |
| `output/tables/wgcna_stability_<comp>.csv`, `wgcna_stability_modules_<comp>.csv`, `wgcna_stability_summary.csv` | 08c: split-half and flagged-removal reproducibility |
| `output/tables/wgcna_preservation_<comp>.csv`, `wgcna_preservation_secondary.csv`, `wgcna_preservation_summary.csv` | 08d: ERCB preservation (`Zsummary.qual`, `qsummary`, medianRank, verdict) |
| `output/tables/wgcna_hub_genes.csv` (708), `phaseE_candidate_space.csv` (3,279), **`phaseE_candidate_core.csv` (515)**, `wgcna_hub_gene_overlap.csv` | 08e: hubs of preserved modules and the Phase E candidate space |
| `output/wgcna/` | 08a–08f figures (PDF + JPEG): sample/module dendrograms, soft power, eigengene correlation, module–trait heatmaps, module significance, ME boxplots, stability bars, preservation Z, hub networks, kME distributions, module-identity dotplots |
| `output/tables/08f_module_enrichment_all.csv` (70,117), `08f_module_enrichment_top.csv`, `08f_module_programme_labels.csv` (16) | 08f: per-module ORA over GO/Hallmark/Reactome/KEGG and the programme label assigned to each module |
| `output/tables/09a_panel_nestedcv.csv` (150), `09a_panel_selection_frequency.csv`, `09a_panel_consensus.csv` (5), `09a_panel_cv_summary.csv`, `09a_panel_runinfo.csv` | 09a: nested-CV results, per-gene selection frequency, the locked panel, design record |
| `output/tables/09b_panel_coefficients.csv`, `09b_panel_external.csv`, `09b_panel_permutation.csv`, `09b_panel_decision_curve.csv` | 09b: ORs with CIs, per-cohort AUC/CI/sens/spec/calibration, the two nulls, net benefit |
| `output/models_panel/` | 09a–09d figures: selection-frequency tile, consensus counts, nested-CV AUC, ROC, score-by-diagnosis, coefficient forest, decision curve, calibration, robustness ROC, recalibration, gene-by-cell-type atlas map, pseudobulk heatmap, module-vs-panel, leave-one-gene-out |
| `output/tables/09c_panel_robustness.csv`, `09c_panel_recalibration.csv`, `09c_crossplatform_coefficients.csv` | 09c: transfer/robustness AUCs, before/after recalibration, platform-to-platform coefficient comparison |
| `output/tables/09d_panel_expression_diagnostics.csv`, `09d_panel_pseudobulk.csv`, `09d_module_vs_panel.csv`, `09d_leave_one_out.csv` | 09d: detection + atlas localisation, pseudobulk MN-vs-control, module-vs-gene comparison, leave-one-gene-out |
| `output/tables/09e_panel_by_comparator.csv`, `09e_authors_matrix_check.csv`, `09e_technical_confound.csv`, **`09e_no_slamf8_model.csv`** (renamed 2026-11-17 from `09e_four_gene_model.csv`) | 09e: per-diagnosis AUC, submitters'-matrix reproduction, PC confound, the model without SLAMF8 (3 genes with the current panel) |
| `output/tables/09f_pla2r_strata.csv`, `09f_pla2r_expression.csv`; `metadata/GSE216841_pla2r.tsv` | 09f: PLA2R-stratified panel performance and the reusable PLA2R annotation |
| `output/figures/` (+ `supplementary/`) | **12 plates: Figures 1–6 and S1–S6** as vector PDFs with embedded fonts plus 300 dpi JPEG previews (`work/figtex/` holds the LaTeX sources). Fig 3 = the separability gate; S4 = functional profiling; S5 = panel-gene localisation + PLA2R; S6 = gate partition stability |
| `Manuscript/` | `Paper1_manuscript.md/.docx` (**4,520 words**), `IF_validation_plan.md/.docx` (1,051 words), `DATA_REQUESTS.md/.docx` (957 words) |
| `environment/` | `renv.lock` (1,669 packages), `packages_workspace|user|site_library|library.csv`, `sessionInfo.txt`, `README.md` (rebuild recipe) |
| `Manuscript/` | `Paper1_manuscript.md/.docx`, `IF_validation_plan.md/.docx`, `DATA_REQUESTS.md/.docx` — all rebuilt from the `.md` with `pandoc --reference-doc=reference/tnr_fangsong_redcode.docx` |
| `output/figures/` (+ `supplementary/`) | Figures 1–6 and S1–S6: vector PDFs assembled by `tools/make_figures.py` from the source panel PDFs, fonts embedded, one page each, plus 300 dpi JPEG previews. LaTeX sources in `work/figtex/` |
| `tools/make_figures.py`, `tools/make_environment_record.R` | Figure compositor (LaTeX + Ghostscript, staged through `/tmp`) and the environment-record generator |
| `metadata/GSE216841_pla2r.tsv` | PLA2R status per sample (`gsm`, `sample_key`, `PLA2R`) recovered from the series matrix — the annotation that resolves B12 |
| `logs/` | Run logs: `rma_arrays`, `rnaseq_normalize`, `qc`, `build_cohort_sets`, `repair_counts`, `qc_policy`, `combat`, `deg_limma`, `meta_analysis`, `enrichment`, `separability_gate`, `gate_robustness`, `subtypes`, `scrna_install`, `scrna_manifest`, `scrna_qc`, `scrna_atlas`, `scrna_pseudobulk`, `scrna_mapping`, `scrna_deconvolution`, plus `verify_phase_a/b`, `verify_downstream`, `verify_scrna` |

**Reference/annotation:** Brainarray Human ENTREZG custom CDFs, version 19 (shipped inside the GEO tars); HGNC symbols via `org.Hs.eg.db`; gene sets MSigDB v2025.1.Hs. No genome build is used — all data are pre-quantified expression matrices.

### 6.4 Script map
| Script | Role | Runtime |
|---|---|---|
| `R/03d_qc_policy.R` | freeze the QC outlier policy | seconds |
| `R/04a_combat_integrate.R` | ComBat + PCA/silhouette verification | ~30 s |
| `R/04b_deg_limma.R` | per-cohort limma DEGs + ComBat sensitivity | ~2 min |
| `R/04c_meta_analysis.R` | random-effects meta + ERCB validation + figures | ~4.5 min (**run in background**) |
| `R/04d_enrichment.R` | ORA + GSEA + expected-signal audit | ~5 min (**run in background**) |
| `R/05a_separability_gate.R` | nested-CV + ERCB transfer of the MN-vs-other-GN classifier | ~30 s |
| `R/05b_gate_robustness.R` | permutation null, QC sensitivity, reverse direction, per-diagnosis AUC | ~2 min |
| `R/06a_subtypes.R` | exploratory ConsensusClusterPlus per platform + concordance + ERCB mapping | ~30 s |
| `R/99_verify_phase_a.R` | 51 Phase A assertions | ~10 s |
| `R/99_verify_phase_b.R` | 86 Phase B assertions | ~5 s |
| `R/99_verify_downstream.R` | 32 gate + subtype + gate-stability assertions | ~5 s |
| `R/07a_scrna_manifest.R` | archive validation + 10x inventory as dir+prefix pairs (never downloads, never links) | seconds |
| `R/07a0_scrna_metadata.R` | series matrix -> `scrna/sample_metadata.tsv` (disease labels per sample) | ~10 s |
| `R/07b_scrna_qc.R` | per-sample import (`ReadMtx` by prefix), MAD QC, doublets, normalisation (resumable) | **~25 min for 49 samples** (`run_in_background`; ~1–3 min/sample, the 23.5 k-cell sample is the slow one) |
| `R/07c_scrna_integrate_annotate.R` | HVG merge + JoinLayers + ScaleData/PCA + Harmony + clustering + marker annotation | **~25 min** for 236 k cells at 2,019 genes |
| `R/07d_scrna_pseudobulk.R` | donor × cell-type pseudobulk, limma-voom DE, composition | **~7 min** |
| `R/07e_bulk_scrna_mapping.R` | signature matrix (name-aligned), tau attribution, bulk↔pseudobulk concordance | **~1 min** |
| `R/07f_deconvolution.R` | NNLS deconvolution of the 7 bulk sets, primary + all-gene sensitivity | **~2 min** |
| `R/07z_install_scrna_pkgs.R` | single-cell package install (7 core packages; liana unavailable in Bioc 3.23) | done 2026-10-18 |
| `R/99_verify_scrna.R` | scRNA infrastructure (16) + data assertions | ~5 s |
| `R/08a_wgcna_network.R` | signed WGCNA per compartment + soft power + outlier sensitivity | **~4 min** (both compartments) |
| `R/08b_module_trait.R` | eigengene ↔ trait correlations, composition-adjusted fit, module significance | **~1.5 min** |
| `R/08c_wgcna_stability.R` | 3 stratified split halves + flagged removal, per-subsample power | **~25 min** (`run_in_background`) |
| `R/08d_module_preservation.R` | ERCB preservation, 100 permutations (+50 cross-platform) | **~15 min per primary test** (`run_in_background`) |
| `R/08e_hub_genes.R` | hub genes of preserved modules + Phase E candidate space + hub networks | **~1.5 min** |
| `R/99_verify_phase_d.R` | 117 Phase D assertions | ~10 s |
| `R/08f_module_enrichment.R` | ORA per module + programme labels (resumable, `--force` to recompute) | **~12 min** first run, seconds on resume |
| `R/09a_panel_selection.R` | nested 10×5 CV, three selectors, consensus panel | **~1 min** (10 workers) |
| `R/09b_panel_validation.R` | locked model, ERCB/RNA-seq evaluation, 2 permutation nulls, DCA, calibration | **~1 min** |
| `R/09c_panel_robustness.R` | compartment transfer, third cohort, CDF-revision robustness, recalibration, cross-platform training | **~1 min** |
| `R/09d_panel_mechanism.R` | panel-gene diagnostics, atlas localisation, pseudobulk, module-vs-panel, leave-one-gene-out | **~1 min** |
| `R/09e_panel_sensitivity.R` | per-comparator AUC, pipeline-independence check, PC/batch confound, four-gene model | **~1 min** |
| `R/09f_pla2r_strata.R` | PLA2R characteristic parsing + stratified panel performance | **~30 s** |
| `R/05c_gate_seed_stability.R` | refits the locked gate under **10 independent RNG streams** at *K* = 1000 and *K* = all; reports the AUC spread (the F62 "is it one lucky draw?" check) | **~6 min** (`run_in_background`) |
| `tools/make_figures.py` | composite Figures 1–6 + S1–S6 as vector PDF (LaTeX) with embedded fonts + JPEG previews; page height measured per plate with `gs -sDEVICE=bbox`, then a **one-page guard** regrows from the row budget if LaTeX breaks the plate | **~2–3 min** |
| `tools/make_environment_record.R` | `environment/renv.lock`, per-library package CSVs, sessionInfo, rebuild README | ~10 s |
| `R/99_verify_phase_e.R` | **106** Phase E + module-identity + robustness/mechanism/sensitivity/PLA2R assertions | ~10 s |

### 6.5 Cross-machine resume
1. Copy the whole project directory (at minimum `R/`, `Datasets/`, `metadata/`, `output/`, `reference/`, and the `*.rds` under `work/`; `work/cel` and `work/rma_cache` are regenerable).
2. Install R 4.6.x, then restore the user library packages listed in §6.2(2) from CRAN/Bioconductor.
3. Rebuild the workspace library: `Rscript R/00d_install_cdf.R` then `Rscript R/02a_build_cdf.R` followed by `R CMD INSTALL -l .Rlib work/cdf/<pkg>` for each of the 5 CDF packages (they are built from the `.cdf.gz` files inside `Datasets/*/raw/`, so no network needed).
4. Sanity check — **all six verifiers must exit 0**: `Rscript R/99_verify_phase_a.R` (51 checks), `Rscript R/99_verify_phase_b.R` (86), `Rscript R/99_verify_downstream.R` (27), `Rscript R/99_verify_scrna.R` (42), `Rscript R/99_verify_phase_d.R` (117), `Rscript R/99_verify_phase_e.R` (**103**).
5. The scRNA arm needs only `work/scrna/qc/*.rds` plus `atlas_cells.csv` to re-run 07c→07f; `scrna/extracted/` (3.3 GB, read in place) is regenerated by `Rscript R/07a_scrna_manifest.R` from the archives in `../Datasets/` — `scrna/raw/` is intentionally empty.
6. **Figure rebuild needs TeX Live and Ghostscript**: `python3 tools/make_figures.py` composites the plates with `pdflatex` and rasterises previews with `pdftoppm`; the Ghostscript embedding pass is staged through `/tmp` because `gs` cannot read this exfat volume (P-49).
7. **`99_verify_phase_e.R` also covers 09c/09d/09e/09f** (compartment specificity, third cohort, CDF-revision robustness, recalibration, cross-platform sign concordance, panel-gene localisation, leave-one-gene-out, module-vs-panel).
7. **Data is no longer inside the project.** On another machine, copy `<MN>/Datasets` too and either keep it as a sibling of `Paper_1/` (auto-detected) or set `PAPER1_DATA=/path/to/Datasets`.

### 6.6 API keys / credentials
None used. CIBERSORTx (now only an optional cross-check, see B3) would need a **licence token** from the CIBERSORTx site — it does not exist on this machine.

---

## 7. Startup Instruction for the Next Session

> **When starting a new session, the very first message must be: "Read handoff.md first."**

**Checkpoint of this refresh (2026-11-18 12:15):** the analysis is **consistent, verified and
submission-ready from the computational side.** Three things a new session must not do:
(1) quote the pre-audit figures in §2.5 and §2.11–§2.15 — **the panel is four genes**;
(2) quote the F07-only gate range — the gate's ERCB AUC is **0.861–0.922** under the harmonised
transfer and **partition-stable** (0.912–0.928 over 10 draws, Fig S6), **not** 0.674–0.837 and **not**
the old transductive 0.866–0.915;
(3) read the `gsm` column of `sample_metadata.csv` as a GEO accession for the RNA-seq cohorts — it holds
`Sample_N` local indices there (see `id_type`), and 134 of them collide between the two RNA-seq series.

Current state: six verifiers exit 0 — **54 / 86 / 42 / 117 / 106 / 32**; 12 plates, every one **one
page** with all fonts embedded and every one cited in the Results; 221 result tables; 48 R scripts;
`Docs/` holds the BAP audit, the graded Findings Ledger (10 E2 / 4 E1 / 0 E4) and the consolidated
results. The refit-based permutation null for the panel is *p* = 0.11, i.e. **not significant**, and the
manuscript says so explicitly.

The new session should restate the current state from this document to confirm understanding —
specifically: 13 accessions → 2,165 samples → 453 MN samples → **178 unique MN patients** → 7
non-redundant cohort sets → Phase A 54/54 → Phase B 86/86 (584 glomerular MN-vs-control meta-DEGs,
`n_dl_fallback = 0`) → separability gate passes with **nested-CV AUC 0.911–0.926 and external ERCB
0.861–0.922** (0.922 at K = all; permutation *p* = 0.048) → exploratory subtypes **negative** →
scRNA arm: 432,378 cells → 236,107-cell atlas, 16 cell types → deconvolution 1,328 samples, **67/316**
comparisons significant, podocyte-rich versus other nephropathies (discovery adjusted *p* = 3.3 × 10⁻⁴)
but **not** versus healthy glomeruli (adjusted *p* = 0.29) → Phase D: 9 glomerular (power 18) + 7 TI
(power 10) modules, **all preserved in ERCB** (median Z 29.0 / 26.6), MN signal in **pink** (84 genes,
adjusted *p* = 1.0 × 10⁻⁵) and **red** (144 genes, adjusted *p* = 6.6 × 10⁻⁴), 708 hubs → 515-gene core
→ Phase E: nested 10 × 5 CV with disjoint per-repeat folds over 287 ERCB-measurable core genes gives a
locked **four-gene panel (SLAMF8, TRPC6, MYOZ2, PPP1R3C)** with ERCB AUC **0.816 (0.721–0.910)**,
RNA-seq 0.901, LCM 0.929, other-CDF 0.848, tubulointerstitium 0.449/0.539 (compartment-specific),
label-permutation *p* = 0.005 and **refit null *p* = 0.11** → recalibration 0.862 → −0.003 with Brier
0.845 → 0.096 → a three-gene model without SLAMF8 is comparable (ERCB 0.845), MYOZ2 alone reaches 0.903
→ PLA2R: 0.976 / 0.881 versus MCD, 0.833 within MN on 6 vs 6 (uninterpretable) → deliverables: **twelve**
vector plates (Fig 1–6, S1–S6) each one page with 0 non-embedded fonts, EBioMedicine-format manuscript
(377-word Summary, 44 verified references), cover letter, IF plan, data-request pack, `renv.lock` → GEO inputs live in the shared
`MN/Datasets` → **the remaining work is external: B2 clinical variables, B7 participant map, and IF
validation of MYOZ2 first.**

If the previous session was interrupted, first record what actually completed before the interruption. The last verified state is: `Rscript R/99_verify_phase_a.R`, `99_verify_phase_b.R`, `99_verify_downstream.R`, `99_verify_scrna.R`, `99_verify_phase_d.R` **and** `99_verify_phase_e.R` all exit 0 (**54 / 86 / 32 / 42 / 117 / 106** checks). If any now fails, diagnose that before any new analysis — and check for 0-byte outputs first (P-28).
