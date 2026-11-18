# Handoff — Paper 1 (frozen) + Track B (active)

**Project:** machine-learning biomarker discovery for membranous nephropathy (MN) from multi-cohort bulk transcriptomics, plus Track B — a leakage-audited benchmark of *published* MN panels.
**Project root:** `/run/media/meng/MY_Data/MN/Paper_1` — scripts resolve the root from `R/_lib.R` (override `PAPER1_PROJ`).
**Data root:** `/run/media/meng/MY_Data/MN/Datasets` (shared with Paper_2/Paper_3) — `R/_lib.R` sets `DATA`; override `PAPER1_DATA`. The GEO inputs live outside the project.
**Last updated:** 2026-11-18. **This file keeps only what a future session needs** (current numbers, open work, claim limits, pitfalls, environment). The full session-by-session record is in git: `git show 9230c6a:handoff.md`.

---

## 1. Where things stand

**Paper 1** (Phases 0, 1, A–E) is complete, re-derived under the audited pipeline and **frozen — add no methods to it**. Claims, figures, tables and seven verifiers are final.

**Track B** is the live workstream: **12** published MN panels carried down a 7-rung attrition ladder (R1 apparent → R2 de-duplication → R3 nested CV → R4 transductive external → R5 harmonised external → R5r RNA-seq → R6 refit null). The pre-registered hypothesis (median R5 AUC < 0.70) was tested on the original 6-panel registry and **falsified** (median 0.821). The registry was expanded to 12 panels on 2026-11-18 (B15): median R5 is now **0.674** (8/12 below 0.70), so transfer is real for some panels but far from uniform. The robust expanded claim is **0/12 beats its own refit permutation null** (*p* = 0.055–1.000), and external discrimination still does not improve with panel size.

**Chain to restate on resume:** 13 accessions → 2,165 samples → 453 MN samples → **178 unique MN patients** → 7 non-redundant cohort sets → Phase A 54/54 → Phase B 86/86 → gate passes (nested CV 0.911–0.926, ERCB 0.861–0.922, *p* = 0.048) → subtypes **negative** → scRNA 432,378 cells → 236,107-cell atlas, 16 cell types → deconvolution 1,328 samples, 67/316 significant → Phase D modules **all preserved**, MN signal in pink + red → Phase E **four-gene panel**, ERCB 0.816 → recalibrated (Brier 0.845 → 0.096).

**Verified state (re-derived 2026-11-19, not copied):** seven verifiers exit 0 — **54 / 86 / 32 / 42 / 117 / 106 / 37**; **12 plates** (main Fig 1–6: 5, 6, 5, 5, 5, 5 panels = 31; supplementary S1–S6: 5 each = 30; **61 panels**), each one page, 0 non-embedded fonts; **183 figure PDFs / 183 rasters**; **225 result CSVs** (113 at `output/tables/` top level); **52 files in `R/`** (44 pipeline + 7 verifiers + `_lib.R`; Track B = `10b`, `10c`, `10d`, `99_verify_trackb`); project ~32 GB.

**Version control:** git since 2026-11-18; `git log` is authoritative. `core.fileMode`/`core.autocrlf` = `false` (exfat). `work/`, `.Rlib/`, `scrna/extracted/`, `logs/` and five bulk CSVs > 3 MB are gitignored as regenerable — a checkout alone is not a reproducible run. Commit `handoff.md` at every checkpoint.

---

## 2. Numbers that must be quoted (current only)

Every value below is script-generated and asserted by a verifier; if prose and table disagree, the table wins. Do **not** resurrect the pre-audit numbers (five-gene panel, ERCB 0.88, gate 0.866–0.915, "724 genes by DL fallback") — they are wrong, not alternative.

**Paper 1 — locked results**

| Block | Quantity | Value |
|---|---|---|
| Design | accessions / samples / MN samples / unique MN patients | 13 / 2,165 / 453 / **178** |
| Design | cohort sets | 7 (2 discovery array, 2 ERCB validation, 2 RNA-seq cross-platform, 1 LCM) |
| QC | flagged samples | 200, **retained** (annotation only, policy frozen) |
| Phase B | glom meta-DEGs (Ctrl / FSGS / MCD / Other / nonMN) | **584 / 18 / 60 / 16 / 18** |
| Phase B | TI meta-DEGs (same contrasts) | **323 / 2 / 4 / 7 / 1**; MN vs Control is *k* = 1 (no TI donors in RNA-seq) |
| Phase B | DL fallback count | **0** everywhere (the "724 genes" claim was wrong) |
| Phase B | ERCB replication of meta-DEGs | 55.3 % glom / 32.4 % TI (MN vs Control) |
| Gate | nested CV (10 × 5, disjoint) | **0.911–0.926** |
| Gate | external ERCB across *K* | **0.861–0.922** (*K* = all 0.9217, CI 0.862–0.982) |
| Gate | RNA-seq independent / secondary endpoint (incl. donors) | **0.928–0.939** / **0.9007** |
| Gate | permutation null median / 95th / *p* | 0.500 / 0.674 / **0.048** |
| Gate | QC-flagged removed / reverse direction | 0.939 / 0.772 |
| Gate | partition stability (10 RNG streams) | *K* = 1000 **0.9003–0.9032**; *K* = all **0.9118–0.9281** |
| Gate | per-diagnosis ERCB AUC | 0.71–0.98 (FSGS 0.83, FSGS_MCD 0.71, RPGN 0.95, DN 0.89) |
| Subtypes | verdict | **negative** (array *k* = 2, PAC 0.0015 but technically confounded; RNA-seq none; Jaccard 0.007) |
| scRNA | cells / atlas / cell types | 432,378 / 236,107 × 2,019 genes / **16** |
| Deconvolution | samples / degenerate fits / significant comparisons | 1,328 / **0** / **67 of 316** |
| Deconvolution | signature fragility | filtered vs all-gene signature disagree on dominant cell type in **1,314 / 1,328** |
| Deconvolution | podocyte | rich vs **other nephropathies** *p* = 3.3e-4; **not** vs healthy glomeruli (*p* = 0.29) |
| Phase D | modules / preservation | 9 glom (power 18) + 7 TI (power 10); **all real modules preserved** (median Z 29.0 / 26.6); grey + gold fail |
| Phase D | MN modules | **pink** 84 genes adj *p* = 1.0e-5; **red** 144 genes adj *p* = 6.6e-4 (both survive composition adjustment) |
| Phase D | hubs → core candidate space | 708 hubs → **515-gene core** (3,279 with provenance) |
| Phase E | panel | **SLAMF8, TRPC6, MYOZ2, PPP1R3C** (287 ERCB-measurable core genes, nested 10 × 5) |
| Phase E | nested CV (LASSO / SVM-RFE / Boruta) | 0.846 / 0.843 / 0.814 |
| Phase E | discovery apparent | 0.922 (0.878–0.967) |
| Phase E | **ERCB external** | **0.816 (0.721–0.910)** |
| Phase E | RNA-seq / LCM / other-CDF / incl. donors | 0.901 / 0.929 / 0.848 / 0.837 |
| Phase E | tubulointerstitium (ERCB / RNA-seq) | **0.449 / 0.539** — chance; compartment specificity |
| Phase E | label-permutation *p* / refit-permutation *p* | **0.005** / **0.11 — NOT significant; say so** |
| Phase E | three-gene model (no SLAMF8) | ERCB **0.845 (0.753–0.937)**, RNA-seq 0.891, LCM 0.941 |
| Phase E | single genes | MYOZ2 alone 0.903; worst leave-one-out 0.738 (without MYOZ2) |
| Phase E | calibration-in-the-large / Brier | +0.862 → **−0.003** / 0.845 → **0.096** |
| Phase E | cross-platform transfer (RNA-seq → ERCB) | AUC 0.850, **4/4 coefficient signs** concordant |
| Phase E | per-comparator ERCB | FSGS **0.632**, FSGS_MCD 0.484, RPGN 0.714, DN 0.722, MCD 0.878, SLE 0.872, IgAN 0.940, HT 0.959, TMD 0.968 |
| Phase E | PLA2R | 0.976 (PLA2R+) / 0.881 (PLA2R−) vs MCD; **0.833 within MN, 6 vs 6 → uninterpretable** |
| Manuscript | Summary / main text (Intro+Results+Discussion) / references | 412 words / 3,929 words / 44 verified |

**Track B — current**

| Quantity | Value |
|---|---|
| Panels / ladder | **12** registered (P01–P12) × 7 rungs = **84 rows** (`output/tables/10c_attrition.csv`) |
| Median AUC by rung | R1 **0.675** / R3 **0.625** / **R5 0.674** (range 0.359–0.893; **8/12 below 0.70**) |
| Registry split (descriptive) | P01–P06 median R5 **0.821** vs P07–P12 **0.500** — the six added panels transfer worse |
| Refit null | **0/12 significant** (*p* = 0.055–1.000); our own panel *p* = 0.11 sits in the same place |
| Standardisation effect | R4 → R5 moves AUC by up to **0.105** **in both directions** |
| Panel-size frontier | internal CV peaks at *K* ≈ 12 then overfits; external AUC **flat 0.767–0.858 from *K* = 2 to 20 and collapses at *K* = 50 (0.690)** |
| Frontier identities | top-4 by moderated *t* = SULF2, MYOZ2, TRPC6, ARHGAP20 (shares only 2 with the locked panel); best external = 2-gene SULF2 + MYOZ2, 0.858 |

Artefacts for the numbers: `output/tables/09a_*…09f_*`, `separability_*`, `gate_robustness.csv`, `05c_gate_seed_stability.csv`, `10c_*`, `10d_*`, plus `Docs/CONSOLIDATED_RESULTS.md`, `Docs/BAP_AUDIT_2026-11-18.md`, `Docs/Findings_Ledger.csv` (15 graded claims), `AUDIT_2026-10-20.md` (37 findings + per-change log).

---

## 3. Open work

| # | Item | Owner | What to do / resume from |
|---|---|---|---|
| B15 | **Track B panel registry** | **done 2026-11-18** | Expanded 6 → **12** panels (P07–P12: Luo 2026, Fan 2025, Gui 2024, Shao 2025, Hui 2025, Zhao 2026) from PubMed/PMC full text; **15 exclusions** logged with reasons in `metadata/trackB_exclusions.csv`; method in `metadata/trackB_retrieval_log.md`. The ≥ 10 bar is met. Remaining user decision: submit the 12-panel benchmark as is, or fund a further retrieval pass (no qualifying paywalled candidate is currently identified). |
| P0 | **IF validation of the panel proteins** | user (wet lab) | `Manuscript/IF_validation_plan.md`: pre-specified criteria per target, MN 10–15 (PLA2R+/− split), FSGS 6–8, MCD 4–6, normal 4–5, lineage anchors, blinded per-glomerulus quantification. Order: **MYOZ2** (unresolved cellular source, strongest single gene) → TRPC6 → PPP1R3C → SLAMF8 (barely detected on 2 platforms). ICAM1 is out of the panel and needs nothing. |
| P0b | **Submit the B2 / B7 data request** | user (long lead time) | `Manuscript/DATA_REQUESTS.md`: NIDDK Central Repository R4R, study DOI `10.58020/f2ae-n094`; DACC `NEPTUNE-STUDY@umich.edu`; ask for the subject→sample map in the same request. This is the only thing between the biomarker paper and a prognostic one — start early. |
| P1 | **Submission to *EBioMedicine*** | agent + user | User: authors/affiliations/ORCIDs/CRediT (currently "to be confirmed"), funding statement, 3 suggested reviewers, sign-off on declarations. Agent on request: trim Summary 377 → ~330 words if the system enforces a cap; confirm exact caps (publisher PDF is Cloudflare-blocked, P-50); renumber citations if the author list changes. Checklist: `Manuscript/SUBMISSION_CHECKLIST.md`. |
| P1b | **Optional consolidation analyses — DEFERRED (user, 2026-11-18)** | agent on request | (a) **MYOZ2 provenance test** — correlate with ACTA2/MYH11/TAGLN/PDGFRB across the 4 cohorts + snRNA subset; converts "source unresolved" into a tested statement (cheap). (b) **third deconvolution method** (nu-SVR/CIBERSORT-style, or top-500 specific-gene signature) — BAP wants two methods (moderate). (c) leave-one-cohort-out meta on the 584 meta-DEGs, subtyping power/ARI bootstrap, RNA-seq cohort relatedness check (all cheap). **Not worth doing:** more GEO series (adds samples, not variables, re-imports the redundancy Phase 1 removed); raising the permutation count (cannot fix a null with median 0.537 at 21 ERCB cases). |
| P2 | **Phase F subtypes** | deferred | Analysis done and **negative**; resume only with new data. Blocked on B2 (eGFR/UPCR/remission). PLA2R no longer blocks anything. |
| P4 | **Clinical variables** | user | Route in `Manuscript/DATA_REQUESTS.md`. On arrival: `metadata/neptune_clinical.csv` → append to `wgcna_traits.csv` → re-run `08b`, update the trait-count assertion in `99_verify_phase_d.R`. |
| B13 | `07b` doublets unseeded | agent on request | `scDblFinder()` has no `set.seed()`; only matters if `07b` is re-run (~25 min). Atlas numbers are frozen. |
| B14 | Missing parallelism waivers | agent, low priority | Most pipeline scripts lack a `# Parallelism: not needed — <reason>` line; style/protocol only. Do not parallelise settled steps (summation-order risk). |
| B16 | Adoption-test packages absent | user (install) | `decoupleR`, `InstaPrism`, Python `tabpfn` are **not** installed; `BayesPrism`, `MuSiC`, `hdWGCNA` **are**. So the two-method deconvolution (G6) can start; TabPFN and decoupleR cannot. |
| P5 | GSE73953 (PBMC, Agilent) | agent, optional | Never normalised; add an Agilent branch to `02b` (8 MN samples, low priority). |
| P7 | Wire Paper_2 / Paper_3 to `MN/Datasets` | user, when they start | Paper_2 needs 5 tier-2 accessions (`GSE286184`, `GSE283826`, `GSE341102`, `GSE199850`, `GSE193512`); Paper_3 needs eQTL/GWAS summary stats instead. Keep `Datasets/manifest.tsv`'s `projects` column current. |

---

## 4. Claim limits (do not overstate)

- **The gate model beats the panel** in ERCB (0.922 vs 0.816). The panel's claim is interpretability and transportability, **not** superiority.
- The **refit-permutation null is not significant** (*p* = 0.11, 21 ERCB cases). Report both nulls; never claim 0.05 significance for the strict one.
- **Array-to-array transfer is the hard direction**; cross-platform transfer to RNA-seq is easier. Never revert to cohort-internal z-scoring to recover the higher old AUC (P-58, P-62).
- **Tubulointerstitial = chance** (0.449 / 0.539). TI separates disease from donors, not MN from its mimics; there is **no TI MN module**. The TI null is the strongest specificity result — present it as such.
- **FSGS is the hard comparator** (0.632 ERCB). Report per-diagnosis AUCs; do not claim a diagnostic multi-classifier (the multiclass model's balanced accuracy is 0.44) — it is MN **detection** with one decision boundary.
- **MYOZ2** is the strongest single gene (0.903) but has **no resolvable kidney cellular source** (93.8 % of atlas cells zero). Not a probe artefact: bulk-detected, consistent in all four cohorts. State "tissue-level reproducibility with unresolved cellular source", flag IF first.
- **SLAMF8** is detected in only 9.5 % (RNA-seq) / 27 % (LCM); the three-gene model without it is comparable (ERCB 0.845).
- **Deconvolution** is fragile: two signature variants disagree on the dominant cell type in 1,314/1,328 samples. Report the filtered signature as primary, all-gene as sensitivity. pDC (97 cells) and Neutrophil (535) fractions are unreliable and excluded from traits/attribution.
- **Subtypes are a negative result**; the array *k* = 2 split tracks technical quality (correlation, mean expression).
- **Bulk ↔ pseudobulk concordance is weak** — best Pearson 0.29 (endothelium), 0.27 TAL, 0.15 podocyte, negative for B and T_NK; ~55–60 % direction agreement. A limit of comparing whole-tissue bulk with one cell type's pseudobulk, not a bug. Pseudobulk DE is hypothesis-grade (MN vs control 15 vs 10 donors; 6 MN donors in the reference cohort).
- **No participant linkage** for GSE197307/GSE182380 (no id field). Never infer overlap from sample titles (P-2). `sample_metadata.csv`'s `gsm` column mixes 1,623 real `GSM####`, 473 `Sample_N` indices and 69 submitter ids — read `id_type` / `id_is_globally_unique`.
- **PLA2R within-MN** (6 vs 6) is uninterpretable; the panel is antigen-independent (0.976 / 0.881 vs MCD).
- **Track B now has 12 panels**, not 6. The pre-registration was falsified on the original 6; the expanded registry's median R5 is 0.674, but do **not** present that as confirmation of the pre-registered hypothesis — the registry was expanded after the fact. The defensible claims are: **0/12** panels beat their refit null; transfer is heterogeneous (P01–P06 median 0.821 vs P07–P12 0.500); panel size does not predict external AUC. Panels are vindicated on performance, not gene identity.
- **Not claimed anywhere:** a PubMed/GEO systematic search ("Evidence before this study" was rewritten to state only the retrieved literature); EBioMedicine word/figure caps (unverified — mark "to confirm").

---

## 5. Running, verifying, committing

- **Verifiers (all must exit 0):** `99_verify_phase_a` 54, `phase_b` 86, `downstream` 32, `scrna` 42, `phase_d` 117, `phase_e` 106, `trackb` 37. Each ~5–20 s.
- **Re-run orders:** `08a → 08b → 08c → 08d → 08e` (08f with `--force`), or the whole 08d-onward chain via `bash tools/rerun_from_08d.sh`; panel revisit = `09a → 09b → 09c → 09d → 09e → 09f` then `python3 tools/make_figures.py`; Track B = `10b → 10c → 10d → 99_verify_trackb.R`.
- **Foreground shell caps at 2 min** — run `04c`, `04d`, `07b`, `07c`, `07d`, `08c`, `08d`, `08f`, `05c`, `10d` in the background and read the `.out` log (P-21).
- **Figure rebuild** needs TeX Live + Ghostscript: `python3 tools/make_figures.py` (composites via `pdflatex`, embeds fonts via `gs`, staged through `/tmp` because `gs` cannot read this exfat volume). Assert `pdfinfo → Pages: 1` per plate. Main plates carry **5–6 decisive panels**; **S1–S6 carry 5 support panels each** (`panel_h` = 62 mm): S1 supports cohort/gate/atlas QC, S2–S3 the module figure, S4–S5 the panel figure, and S6 is the stand-alone negative subtyping result. Every panel is cited in order in Results (checked by a citation-order script).
- **`handoff.docx` in the project root is stale** — `handoff.md` is the only source; do not regenerate or read the `.docx`.
- **Mount drift:** the volume alternates `/run/media/{meng,yang}/MY_Data`. `R/_lib.R`, the `00*_download*.sh` helpers and `tools/rerun_*.sh` are path-independent — keep it that way (P-59).

---

## 6. Pitfalls (never repeat)

- **P-1** Compute patient-level overlap across *all* releases before merging any matrix (13 releases → 178 patients; 4 releases 100 % contained in another).
- **P-2** Never infer cross-cohort overlap from sample titles. Check for a participant field, check whether ids are contiguous 1..n (local), inspect the per-disease overlap table. No field → state "patient linkage unavailable", treat cohorts as independent.
- **P-3** `affy` cannot process HuGene-2.1-ST. Use `oligo` + Brainarray probe indices (`02b`).
- **P-4** `preprocessCore::rcModelPLM()$Estimates` has length `nrow + ncol`; it is not the sample-effect vector. Use `stats::medpolish()` `overall + col` (gave r = 0.986).
- **P-5** Never `match()` per gene against the 1.02 M-row probe index (436 s). Build one lookup vector.
- **P-6** Never process a 244-array ST cohort in one gulp (12.1 GB RSS). Chunked two-pass quantile, `CHUNK = 40`.
- **P-7** CDF package names must end in `cdf`.
- **P-8** `make.cdf.package()` only *builds* a source package (then `R CMD INSTALL -l .Rlib`) and resolves input through `getwd()` — pass a bare filename after `setwd()`.
- **P-9** Do not derive array platform from the series matrix (two-platform releases share GSM ids). Read each CEL header (`affyio::read.celfile.header`).
- **P-10** Do not `rbind()` two platforms: align on the gene union, `NA`-fill, and drop incomplete genes before PCA *and* before limma (`lmFit` with `NA` → `sigma = NA` → `eBayes` dies).
- **P-11** Do not redirect stdout to the file the script's `writeLines()` writes.
- **P-12** Do not assume count-style filenames hold counts (`GSE216841_MNd_ncounts_annot.txt.gz` is normalised decimals plus text columns).
- **P-13** Do not use raw donor integers as patient ids — namespace them (`LD-<n>`).
- **P-14** Data acquisition is manual: never download datasets; present a request list. (Package installs and MSigDB GMT annotation were fetched directly and are not covered.)
- **P-15** `topTable(coef = "MN - Control")` does not work on `~0 + group`; use `makeContrasts()` + `contrasts.fit()`.
- **P-16** An `mclapply` worker that throws returns a `try-error`, not `NULL`; `do.call(rbind, …)` then recycles the error message into the table. Filter on `is.data.frame(z) && nrow(z) == 1` and `tryCatch` the model call.
- **P-17** `metafor::rma(method = "REML")` may not converge → wrap, fall back to `"DL"`, count and report the fallbacks.
- **P-18** `enrichplot::gseaplot2()` rejects `fgsea` results; use `fgsea::plotEnrichment`.
- **P-19** `msigdbr` ≥ 26 downloads on first use; parse the pre-fetched GMTs in `work/msigdb/` instead.
- **P-20** Verifier figure assertions must scope to the generating step's filename prefix, never "every PDF in a directory".
- **P-21** The foreground shell caps at 2 minutes regardless of the requested timeout — use `run_in_background`.
- **P-22** `pROC::roc()` needs a 2-level response; ERCB codes MN as `MGN` in `disease_raw`.
- **P-23** `apply()` over an `sapply()`-built centroid matrix has no dimnames; set `colnames()` explicitly.
- **P-24** PAC (≤ 0.1) is the stability gate, not silhouette.
- **P-25** A stable split can be a technical split — test median sample correlation, mean expression and detection rate against every cluster.
- **P-26** `geom_text(stat = "count")` rejects `fill`; precompute `as.data.frame(table(...))`.
- **P-27** Python `str.replace()` with a missing anchor fails silently — assert the anchor exists, and assert on the result (P-52).
- **P-28** A 0-byte file passes `file.exists()`. Check `file.size() > 0`; `find <dir> -size 0` after any crash.
- **P-29** Seurat 5 `merge()` keeps per-sample layers → call `JoinLayers()` immediately. `AddModuleScore`'s `ctrl` is drawn from every bin, so it must sit below the **smallest** bin (2,019-gene atlas: compute bins exactly).
- **P-30** `merge()` renames colliding barcodes (`_1`, `_2`) → rename with a sample prefix before merging and rebuild the identical key downstream.
- **P-31** `do.call(cbind, …)` on named vectors of unequal length aligns by **position** → align on the union, fill 0.
- **P-32** NNLS returns corner solutions unless scaled: (a) quantile-normalise mixtures **together with** the signature columns (never to the signature *mean*), then unit-sum both; (b) filter to τ > 0.5 cell-type-specific genes; (c) compute τ on the **CPM scale before** column normalisation.
- **P-33** API traps: `[.DGEList` has no `drop`; `ncol()` on Seurat 5 is a double; `wilcox.test()` of two all-zero groups returns `NaN` (flag `testable = FALSE`).
- **P-34** This volume is **exfat**: no hardlinks/symlinks, case-insensitive; **`ln -f` deletes the destination first — never use it**; atomic temp-file+link writers (including `write_file`) fail — write with a plain shell heredoc; `mv` within the volume is instant; copying onto exfat **dereferences symlinks** (duplicates, new inodes — catch with `du` vs `ls -l`).
- **P-35** WGCNA 1.74 `modulePreservation`: `referenceNetworks` = numeric index; own modules via `multiColor = list(colours, NA)`; `multiData` = two-level `list(set = list(data = samples × genes))`; no `includeModules`/`verbose`; output is the `.qual` family (`calculateQvalue = TRUE`); statistics carry no labels — map back by `moduleSize` and assert the mapping is unambiguous.
- **P-36** `bicorAndPvalue()` returns `$bicor` (`corAndPvalue(method = "bicor")` does not exist). `mad()` = 0 when > half the values are identical → record the correlation method per pair.
- **P-37** A module colour vector must be named by **gene** (`colnames(dat)`), never by sample.
- **P-38** A split-half stability test must re-select the soft power inside each half; report Jaccard **and** eigengene bicor separately (axes can be near-perfect while membership is 0.3).
- **P-39** `numericLabels = TRUE` makes `signedKME` return `kME_ME1` — renormalise names and `stopifnot` that no own-gene kME is `NA` (otherwise hub selection silently returns 0 genes).
- **P-40** `bulk_deg_celltype_attribution.csv` is keyed `meta_<Glom|TI>_<contrast>`; build the full key and print join counts.
- **P-41** A panel gene must be measurable in every validation cohort (GPC6 was NA in 125/196 ERCB samples) — make ≥ 95 % coverage a **selection** constraint (and assert it in the verifier), never a post-hoc drop.
- **P-42** Never report sens/spec at the Youden point of a low-prevalence cohort. Report a fixed 0.5 threshold, calibration-in-the-large and slope; describe the model as a *discrimination* biomarker.
- **P-43** Build each permutation-null summary from its own distribution (and include `null_max`), never reuse one row's quantiles.
- **P-44** `pROC`'s "No case observation" can mean the **predictor** is `NA` — check `sum(is.na(score))` first.
- **P-45** On a genes × samples matrix use row-wise summaries (`rowMeans`, `apply(…, 1, sd)`). A suspiciously exact AUC of 0.500 means the design matrix went constant.
- **P-46** GSE99339 is split by platform and codes MN as **`MGN`**. When a cohort parses to zero cases, print the filename token census.
- **P-47** `merge()` renames colliding columns (`note.x`/`note.y`) → give derived columns distinct names.
- **P-48** GEO `!Sample_characteristics_ch1` is one line per field group, each a per-sample vector: iterate over lines, match `length(line)` to the sample count, assign by position. The join key may be `!Sample_description` (submitter ids, not GSM).
- **P-49** `gs` cannot read `/run/media` (stage through `/tmp`), and R's `pdf()` leaves Helvetica un-embedded → run `gs -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dSubsetFonts=true` over each plate and verify with `pdffonts`.
- **P-50** Publisher guideline PDFs are Cloudflare-403; read a live open-access article's JATS XML via Europe PMC `fullTextXML` instead, and never quote a cap you could not read.
- **P-51** Never hand-type a citation: fetch by DOI/PMID and assert the title; keep `Manuscript/references.tsv` as the only source.
- **P-52** When a script assembles prose, assert on the assembled text (e.g. `'..' not in text`), not just that it ran.
- **P-53** A nested CV must draw its fold partition once per repeat, outside the fold function; verify pairwise test-set overlap = 0, not the AUC.
- **P-54** `tryCatch(error=)` misses warnings — metafor reports non-convergence as a **warning** and still returns a fit. Use `withCallingHandlers`, record a per-gene flag, report the count; treat any uncounted fallback claim as unverified.
- **P-55** Never quantile-normalise a signature matrix together with its mixtures (it flattens the specificity the τ filter preserves) — see P-32 for the correct scaling.
- **P-56** A new column must be added to **every** branch of a per-gene builder (grep all `return(data.frame(`) — then assert identical column sets.
- **P-57** `reshape()`'s `varying` must come from an explicit allow-list, never "everything except these".
- **P-58** Fixing a transductive standardisation **will** lower external AUC. Never revert to cohort-internal scaling.
- **P-59** Drivers must resolve `ROOT="${PAPER1_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"` and be tested from an unrelated cwd.
- **P-60** Data-driven labels and filenames must be derived (`length(panel)`), never hard-coded; name outputs by role.
- **P-61** Verify a composed figure by measurement: `gs -sDEVICE=bbox` page sizing, trailing blank band scan, ink box vs MediaBox, `pdftotext` title check. Compute axis windows from the data.
- **P-62** "Use training statistics" is **not** a cross-platform transfer: check the standardised |z| median ≈ 1 first, then quantile-map (label-blind) onto the reference distribution before applying training μ/σ.
- **P-63** A re-run is not a re-derivation: grep every script for a private copy of a changed rule (`zscore(`, `std(`, `scale(`) and print a reconciliation of two scripts' version of the same quantity.
- **P-64** Content-driven page sizing can shrink a plate onto a second page: assert `pdfinfo → Pages: 1`, regrow from the row budget, delete `<plate>-N.jpg` leftovers.
- **P-65** A key column that doubles as provenance must carry `id_type` / `id_is_globally_unique`, asserted in the verifier.
- **P-66** Counts in prose go stale the moment an artefact is added — re-derive from the filesystem when you touch the file, and write "as of <date>".
- **P-67** Calibrate any checker against a known value before reporting it; parse CSVs with a real reader (`pdffonts`' `emb` column is `$(NF-4)`, not `$4`).
- **P-68** Never type a number into a script: read it from the artefact with a `stopifnot()` and let the verifier assert reproduction (`abs(x − source) < 1e-6`).
- **P-69** `check(isTRUE(x))` silently fails on vector predicates — wrap in `all()` / `any()`.
- **P-70** A literature-derived registry needs a code-enforced inclusion rule (≥ 2 named genes, PMID or DOI, reported performance) plus an exclusion log with a reason per rejected candidate.

---

## 7. Environment

**Compute:** Linux amd64; 20 threads; 61 GB RAM; 1.8 TB free on the volume; project ~32 GB (`work/` 23 GB, `.Rlib` 4.2 GB, `scrna/` 3.3 GB, `output/` ~650 MB); shared `Datasets` 18 GB; **exfat** (P-34). `work/cel/` (16 GB) and `work/rma_cache/` (1.1 GB) are deletable once `work/rma/` is final.

**R 4.6.1** with two libraries: `<PROJ>/.Rlib` (prepended by `_lib.R`; makecdfenv, the 5 Brainarray CDFs, oligo, WGCNA 1.74, ConsensusClusterPlus, caret, glmnet, randomForest, xgboost, Boruta, metafor 5.2.1, msigdbr, rms, babelgene) and the user library `/home/meng/R/x86_64-pc-linux-gnu-library/4.6` (~1,701 pkgs; Seurat 5.6.0, scDblFinder 1.26.9, harmony 2.0.5, muscat, nnls, limma, edgeR, DESeq2, sva, clusterProfiler, fgsea, ComplexHeatmap, survival). `pandoc` at `/usr/bin/pandoc`. **Absent:** ReactomePA, immunedeconv/xCell/quanTIseq, msigdbdf, liana, decoupleR, InstaPrism, tabpfn. Mirrors: TUNA CRAN/Bioconductor. Reproducibility record: `environment/renv.lock` (1,669 pkgs) + `environment/README.md`; regenerate with `Rscript tools/make_environment_record.R`.

**Key paths**

| Path | Contents |
|---|---|
| `$DATA` = `MN/Datasets` | shared GEO collection: `<GSE>_RAW.tar` + `<GSE>/{raw,processed}/`, indexed by `manifest.tsv` |
| `work/sets/` | the 7 cohort sets (`*_expr.rds`, `*_expr_symbol.rds`, `*_samples.csv`) |
| `work/rma/`, `work/rnaseq/`, `work/combat/`, `work/msigdb/` | normalised cohorts, ComBat TI matrices, MSigDB v2025.1.Hs GMTs |
| `work/scrna/{qc,atlas.rds,pseudobulk/}` | per-sample QC objects, 236k-cell atlas, 16 donor × gene matrices |
| `work/wgcna/` | networks, eigengenes, traits, preservation objects (`modulePreservation` writes its cache into the cwd — run 08d from the project root) |
| `work/bench/10b_harmonised.rds` | Track B matrices: `z_harmonised` (primary) and `z_transductive` (kept only to quantify the difference) |
| `output/tables/` (+ `deg/`, `meta/`, `enrichment/`) | 225 CSVs — the numbers above live here |
| `output/figures/`, `output/figures/supplementary/` | 12 plates (PDF + 300 dpi JPEG); LaTeX sources in `work/figtex/` |
| `output/models_{gate,panel,trackb}/`, `output/{qc,deg,enrichment,wgcna,scrna,subtypes}/` | per-step figures (PDF + JPEG) |
| `metadata/` | `trackB_panel_register.csv`, `trackB_exclusions.csv`, `GSE216841_pla2r.tsv`, `ensg2symbol.tsv` |
| `Manuscript/` | manuscript, cover letter, IF plan, data requests, submission checklist, `references.tsv` (all `.md` + `.docx`) |
| `Docs/` | `BAP_AUDIT_2026-11-18.md`, `Findings_Ledger.csv`, `CONSOLIDATED_RESULTS.md`, superseded 2026-10-18 handoff |
| `R/_pre_audit_backup/` | pre-audit scripts (provenance) |

**Script map** (`*` = run in background, > 2 min)

| Scripts | Role |
|---|---|
| `03d`, `04a`, `04b`, `04c*`, `04d*` | QC policy; ComBat; limma DEGs; random-effects meta + ERCB validation; ORA/GSEA |
| `05a`, `05b`, `05c*` | separability gate (harmonised transfer); permutation/QC/reverse/per-diagnosis robustness; gate partition stability |
| `06a` | exploratory subtypes (negative; kept as a negative result) |
| `07a`, `07a0`, `07b*`, `07c*`, `07d*`, `07e`, `07f`, `07z` | scRNA manifest, metadata, per-sample QC (~25 min), atlas (~25 min), pseudobulk (~7 min), bulk↔scRNA bridge, NNLS deconvolution, package install |
| `08a`, `08b`, `08c*`, `08d*`, `08e`, `08f*` | WGCNA networks, module–trait, stability (~25 min), preservation (~15 min/test), hubs/candidate space, module ORA (resumable, `--force`) |
| `09a`–`09f` | panel selection, validation, robustness, mechanism, sensitivity (+harmonisation check), PLA2R strata |
| `10b`, `10c`, `10d` | Track B: harmonised matrix cache, 7-rung attrition ladder, panel-size frontier |
| `99_verify_*` (7 scripts) | 54 / 86 / 32 / 42 / 117 / 106 / 37 assertions |
| `tools/make_figures.py`, `tools/make_environment_record.R`, `tools/rerun_*.sh` | plate compositor (~2–3 min), environment record, resume drivers |

**Cross-machine resume:** copy the project (at least `R/`, `metadata/`, `output/`, `reference/`, `Manuscript/`, `Docs/` and the `.rds` under `work/`) plus `MN/Datasets`; install R 4.6.x; rebuild `.Rlib` with `00d` → `02a` → `R CMD INSTALL -l .Rlib work/cdf/<pkg>` (×5, no network needed); run all seven verifiers; `scrna/extracted/` regenerates from `Rscript R/07a_scrna_manifest.R`; figure rebuild needs TeX Live + Ghostscript. **No API keys or credentials are used** (CIBERSORTx is an optional cross-check only and would need a licence token).

---

## 8. Startup instruction

1. **Read `handoff.md` first**, then confirm the state by running the seven verifiers (expect **54 / 86 / 32 / 42 / 117 / 106 / 37**, all exit 0). Check for 0-byte outputs before trusting any resume (P-28), and re-derive any count you are about to quote (P-66).
2. Do not quote pre-audit figures, do not add methods to Paper 1, do not hard-code a number a table already stores (P-68), and do not read `gsm` as a GEO accession for the RNA-seq cohorts.
3. **B15 is resolved** (12-panel registry, 2026-11-18); the open decision is whether the 12-panel benchmark ships as is. The remaining user decision is **P1** (author list, funding, reviewers, declarations). Ask which to proceed with; everything else is either blocked externally (B2, B7) or wet-lab work (P0).
