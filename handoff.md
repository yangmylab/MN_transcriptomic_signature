# Handoff — Project Handoff Document

**Project:** machine-learning biomarker discovery for membranous nephropathy (MN) from multi-cohort bulk transcriptomics, plus Track B, a leakage-audited benchmark of published MN panels.

- **Project root:** `/run/media/meng/MY_Data/MN/Paper_1`
- **Shared data root:** `/run/media/meng/MY_Data/MN/Datasets`
- **Branch:** `master`
- **Last updated:** 2026-11-19
- **Current stage:** Paper 1 analyses frozen; manuscript package recomposed and verified; Track B active.
- **Checkpoint review:** no evidence of an interrupted, unrecorded analysis session.

---

## 1. Objective (What we are doing)

- Deliver a submission-ready manuscript reporting a compartment-resolved glomerular transcriptomic signature and a four-gene classifier for MN, with every claim externally validated.
- Preserve the audited Paper 1 claim set: interpretability and transportability, not superiority over the whole-transcriptome gate, no prognostic claim, and no unsupported subtype or cell-of-origin claim.
- Complete Track B, a leakage-audited benchmark of published MN panels, with explicit attrition, external transfer and refit-permutation analyses.
- Success criteria are seven passing verifiers, the frozen Paper 1 results, a logically organized manuscript and figure package, and submission declarations completed by the user.

## 2. Completed Items (What has been done)

### 2.1 Paper 1 analysis chain — complete and frozen

The verified chain is:

**13 accessions → 2,165 samples → 453 MN samples → 178 unique MN patients → 7 non-redundant cohort sets → Phase A 54/54 → Phase B 86/86 → separability gate passed → subtypes negative → 432,378 cells → 236,107-cell atlas, 16 cell types → deconvolution 1,328 samples, 67/316 significant → preserved WGCNA modules → four-gene panel, ERCB AUC 0.816 → recalibrated Brier 0.096.**

Paper 1 must not gain new methods. Current script-generated values are:

| Block | Current verified value |
|---|---|
| Design | 13 accessions; 2,165 samples; 453 MN samples; **178 unique patients**; 7 cohort sets |
| QC | 200 flagged samples retained as annotation only |
| Glomerular meta-DEGs | Control 584; FSGS 18; MCD 60; Other 16; non-MN 18 |
| TI meta-DEGs | Control 323; FSGS 2; MCD 4; Other 7; non-MN 1; TI control has *k* = 1 |
| DL fallback | **0**; never resurrect the incorrect “724 genes” statement |
| ERCB replication | 55.3% glomerular / 32.4% TI for MN versus control |
| Gate | Nested-CV AUC 0.911–0.926; external ERCB 0.861–0.922; RNA-seq 0.928–0.939 |
| Gate robustness | Permutation null median 0.500, 95th 0.674, *p* = 0.048; QC-removed 0.939; reverse 0.772 |
| Gate partition stability | *K* = 1000: 0.9003–0.9032; all: 0.9118–0.9281 |
| Subtypes | Negative: array *k* = 2, PAC 0.00077 but technically confounded; RNA-seq minimum PAC 0.211 and no stable partition; 1 shared defining gene, Jaccard 0.0011 |
| scRNA atlas | 432,378 retained cells; 236,107 × 2,019 genes; 16 cell types |
| Deconvolution | 1,328 samples; 0 degenerate fits; 67/316 significant comparisons |
| Signature fragility | Filtered and all-gene signatures disagree on dominant cell type in 1,314/1,328 samples |
| Podocyte result | Rich versus other nephropathies, adjusted *p* = 3.3e-4; not versus healthy glomeruli, *p* = 0.29 |
| Modules | 9 glomerular (power 18) + 7 TI (power 10); all real modules preserved; grey and gold fail |
| MN modules | Pink 84 genes, adjusted *p* = 1.0e-5; red 144 genes, adjusted *p* = 6.6e-4 |
| Candidate space | 708 hubs → 515-gene core; 3,279 with provenance |
| Locked panel | **SLAMF8, TRPC6, MYOZ2, PPP1R3C** |
| Candidate accounting | 515 cross-compartment core → 353 glomerular → 66 excluded for <95% ERCB coverage → **287 tested** |
| Nested CV | LASSO 0.846; SVM-RFE 0.843; Boruta 0.814 |
| Panel external performance | ERCB **0.816 (0.721–0.910)**; RNA-seq 0.901; LCM 0.929; alternative CDF 0.848; including donors 0.837 |
| Compartment specificity | ERCB TI 0.449; RNA-seq TI 0.539 — chance in both |
| Permutation tests | ERCB-label *p* = 0.005; strict refit *p* = **0.11, not significant** |
| Three-gene sensitivity | ERCB 0.845 (0.753–0.937); RNA-seq 0.891; LCM 0.941 |
| Single genes | MYOZ2 0.903; worst leave-one-out 0.738 without MYOZ2 |
| Calibration | Intercept +0.862 → −0.003; Brier 0.845 → 0.096 |
| Cross-platform transfer | AUC 0.850; all four coefficient signs concordant |
| Comparator limitation | FSGS 0.632; FSGS_MCD 0.484; do not claim a diagnostic multiclassifier |
| PLA2R | PLA2R+ 0.976 / PLA2R− 0.881 versus MCD; within-MN 0.833 on 6 versus 6 is uninterpretable |

The gate outperforms the panel in ERCB (0.922 versus 0.816); the panel’s role is compact interpretability and transportability, not superiority.

### 2.2 Manuscript and figures — recomposed and verified

- The manuscript is at `Manuscript/Paper1_manuscript.md` and `Manuscript/Paper1_manuscript.docx`.
- Current architecture is **13 plates / 65 panels**:
  - Main Fig 1–6 contain 5, 6, 5, 5, 5 and 5 panels (31 total).
  - Fig S1 contains four routine bulk-cohort PCA panels; S2–S7 contain five each (34 total).
- Main Fig 1 is now **cohort design and analytical independence**, not a PCA collection:
  - A: study flow and patient-level de-duplication.
  - B: patient redundancy and release containment.
  - C: glomerular allocation and independence.
  - D: ComBat before/after PCA.
  - E: quantitative batch R² and silhouette with disease structure preserved.
- Routine discovery/ERCB/RNA-seq/LCM PCA moved to **Fig S1A–D**; previous S1–S6 became S2–S7.
- Fig 2 order is control meta-analysis → representative forests → disease mimics → cross-cohort concordance.
- Panel-specific supplementary results are grouped with Fig 5 before the atlas.
- All **65/65** panels are cited and interpreted in first-appearance order in Results: main 1→6, supplementary S1→S7, letters sequential within each plate.
- `Manuscript/Figures/` contains **26 files** (13 PDF + 13 JPEG), byte-identical to `output/figures/`.
- Every plate is one page, has embedded fonts and passes the 5-point glyph audit.
- All 44 references are cited; the reference sequence is 1–44.
- The DOCX has a live PAGE footer, passes OfficeCLI schema validation, and renders as **18 pages without a blank or orphaned final page**.
- Current text size: Summary 390 whitespace tokens; Introduction + Results + Discussion 3,674 whitespace tokens.

### 2.3 Track B — active but fully computed

- The registry contains **12 published panels**, expanded from six on 2026-11-18; 15 exclusions have reasons in `metadata/trackB_exclusions.csv` and methods in `metadata/trackB_retrieval_log.md`.
- The seven-rung ladder (R1 apparent → R2 de-duplication → R3 nested CV → R4 transductive external → R5 harmonised external → R5r RNA-seq → R6 refit null) is complete: **12 × 7 = 84 rows** in `output/tables/10c_attrition.csv`.
- The original six-panel preregistered hypothesis was **falsified** (median R5 0.821). The expanded, post-hoc registry has median R5 **0.674** (8/12 below 0.70), which must not be presented as confirmation of the original preregistration.
- Registry split is descriptive: P01–P06 median R5 0.821 versus P07–P12 0.500; transfer is heterogeneous.
- The robust expanded claim is **0/12 panels beating their refit null** (*p* = 0.055–1.000). Our own four-gene panel has *p* = 0.11 in the same context.
- R4 → R5 changes AUC by up to 0.105 in either direction.
- Internal CV peaks near *K* ≈ 12; external AUC is approximately flat at 0.767–0.858 from *K* = 2 to 20 and falls at *K* = 50 (0.690).
- Best external two-gene frontier is SULF2 + MYOZ2, AUC 0.858.

### 2.4 MYOZ2 provenance — completed as a separate P1b test

`Rscript tools/moyz2_provenance.R` and `Docs/MYOZ2_PROVENANCE.md` established:

- Bulk MYOZ2 does **not** track the SMC/perivascular programme in all four cohorts; SMC composite rho 0.01–0.45.
- NPHS2 association is partly disease-driven; within-MN associations are weak and not consistently significant.
- snRNA detects MYOZ2 in 1,686/210,879 cells (0.80%); scRNA detects 57/21,579 (0.26%). Depth OR is 8.69 per log10 counts (7.47–10.10).
- There is no cell-type enrichment after depth and cell-type adjustment.
- Verdict remains **tissue-level reproducibility with unresolved cellular source**. It is not a vascular/SMC readout and must not be upgraded to “podocyte-derived”.

### 2.5 Verification completed at this checkpoint

| Verifier | Result |
|---|---:|
| `99_verify_phase_a.R` | 54 passed, 0 failed |
| `99_verify_phase_b.R` | 86 passed, 0 failed |
| `99_verify_downstream.R` | 32 passed, 0 failed |
| `99_verify_scrna.R` | 42 passed, 0 failed |
| `99_verify_phase_d.R` | 117 passed, 0 failed |
| `99_verify_phase_e.R` | 106 passed, 0 failed |
| `99_verify_trackb.R` | 37 passed, 0 failed |

Additional verified state: 191 figure PDFs / 191 rasters; 240 result CSVs (121 at `output/tables/` top level); 52 files in `R/`; no zero-byte outputs.

## 3. Current Blockers (Where we are stuck)

### 3.1 P0 — Protein-level validation is not yet performed

- **Symptom:** no IF/protein validation of the four-gene panel; MYOZ2 cellular origin remains unresolved.
- **Impact:** limits translational claims and biological interpretation, especially for MYOZ2 and low-detection SLAMF8.
- **Already tried:** bulk reproducibility checks, atlas localisation, donor-level pseudobulk and the dedicated MYOZ2 provenance analysis. These did not identify a defensible cell of origin.
- **Owner / dependency:** user wet lab. Resume from `Manuscript/IF_validation_plan.md`; priority MYOZ2 → TRPC6 → PPP1R3C → SLAMF8. ICAM1 is not in the locked panel.

### 3.2 P0b — Clinical outcomes and subject-to-sample mapping are missing

- **Symptom:** eGFR, UPCR, remission/relapse and a subject→sample map are unavailable in downloaded public files.
- **Impact:** prevents prognostic modelling, subtype–outcome analysis and definitive participant linkage.
- **Already tried:** local metadata review; PLA2R was recovered only for GSE216841 (6 positive, 6 negative), which is insufficient for stratum-level inference.
- **Owner / dependency:** user must submit the long-lead request described in `Manuscript/DATA_REQUESTS.md`: NIDDK Central Repository R4R, study DOI `10.58020/f2ae-n094`, DACC `NEPTUNE-STUDY@umich.edu`, including the subject→sample map.

### 3.3 P1 — Submission metadata and exact journal limits are incomplete

- **Symptom:** authors, affiliations, ORCIDs, CRediT roles, funding wording, three suggested reviewers and declaration sign-off are still 【To be confirmed — user】. Exact EBioMedicine caps are also 【To be confirmed — user/editor】.
- **Impact:** the package cannot be submitted without inventing metadata or risking a format violation.
- **Already tried:** manuscript and figures are fully formatted and verified. Publisher guideline PDF retrieval was Cloudflare-blocked; no cap has been guessed.
- **Owner / dependency:** user supplies metadata; agent can then update `Manuscript/Paper1_manuscript.md`, references, cover letter and `Manuscript/SUBMISSION_CHECKLIST.md`.

### 3.4 Track B registry decision is unresolved

- **Symptom:** the user must decide whether to ship the 12-panel benchmark or fund another literature-retrieval pass.
- **Impact:** changes the scope and timing of the Track B manuscript, not the computed results.
- **Already tried:** 12 qualifying panels and 15 reasoned exclusions are registered; no qualifying paywalled candidate is currently identified.
- **Owner / dependency:** user decision. Resume from `metadata/trackB_panel_register.csv` and `metadata/trackB_retrieval_log.md`.

### 3.5 Optional independent deconvolution and adoption tests are limited by packages

- **Symptom:** `decoupleR`, `InstaPrism` and Python `tabpfn` are absent. `BayesPrism`, `MuSiC` and `hdWGCNA` are installed.
- **Impact:** a third deconvolution method is possible but TabPFN/decoupleR adoption tests cannot start until installation.
- **Already tried:** NNLS with filtered and all-gene signatures; fragility is documented rather than hidden.
- **Owner / dependency:** user installation if those tests are desired. Any third deconvolution is optional P1b sensitivity work, not an addition to frozen Paper 1.

### 3.6 Deferred / low-priority issues

- **P2 subtypes:** already negative and must not be reopened without new clinical data; no technical bug.
- **B13:** `scDblFinder()` lacks `set.seed()` in `07b`; matters only if the atlas is rerun (~25 min). Current atlas is frozen.
- **B14:** many scripts lack a `# Parallelism: not needed — <reason>` waiver; documentation/style only. Do not parallelise settled steps and risk summation-order changes.
- **P5:** GSE73953 PBMC/Agilent remains unnormalised; optional eight-MN-sample branch in `02b`.
- **P7:** Paper_2 and Paper_3 are not yet wired to the shared dataset manifest; handle only when those projects start.

## 4. Next Steps (What to do next)

1. **P0 — Start IF validation**
   - **Owner:** user wet lab.
   - **Trigger:** laboratory scheduling; no further computational dependency.
   - **Resume:** `Manuscript/IF_validation_plan.md`, target order MYOZ2 → TRPC6 → PPP1R3C → SLAMF8, with blinded per-glomerulus quantification and prespecified lineage anchors.

2. **P0b — Submit the B2/B7 request now**
   - **Owner:** user.
   - **Trigger:** immediate, because of long lead time.
   - **Resume:** `Manuscript/DATA_REQUESTS.md`. On arrival, write `metadata/neptune_clinical.csv`, append to `wgcna_traits.csv`, rerun `08b`, and update the trait-count assertion in `R/99_verify_phase_d.R`.

3. **P1 — Finish the EBioMedicine submission package**
   - **Owner:** user for metadata; agent for document updates on request.
   - **Trigger:** before journal submission.
   - **Dependencies:** blocker 3.3.
   - **Resume:** `Manuscript/SUBMISSION_CHECKLIST.md`, `Manuscript/Paper1_manuscript.md`, `Manuscript/Cover_letter.md`, and `Manuscript/references.tsv`. Confirm caps; trim the Summary only if the live submission system requires it. Renumber citations if the author list changes.

4. **P1 — Decide whether Track B ships at 12 panels**
   - **Owner:** user.
   - **Trigger:** before Track B manuscript finalisation.
   - **Resume:** `metadata/trackB_panel_register.csv` and `metadata/trackB_exclusions.csv`. If another retrieval is funded, resume from `metadata/trackB_retrieval_log.md`, then rerun `10b → 10c → 10d → 99_verify_trackb.R`.

5. **P1b — Optional consolidation only if explicitly requested**
   - **Owner:** agent.
   - **Trigger:** user elects a third deconvolution method or sensitivity analyses.
   - **Resume:** `R/07f_deconvolution.R`, `Docs/BAP_AUDIT_2026-11-18.md`, and installed BayesPrism/MuSiC options. Other cheap options are leave-one-cohort-out meta-analysis, subtype power/ARI bootstrap and RNA-seq relatedness checks. Do not add GEO series or increase permutation counts merely to seek a more favourable result.

6. **P2 — Reopen Phase F only with new clinical data**
   - **Owner:** user supplies data; agent reruns on request.
   - **Trigger:** arrival of eGFR/UPCR/remission or other genuinely new molecular data.
   - **Resume:** current negative results in `output/tables/subtype_*` and the updated negative-subtype legend in `Manuscript/Paper1_manuscript.md`.

7. **P2 — Housekeeping only when relevant**
   - Seed `scDblFinder` before any `07b` rerun; add parallelism waivers without changing settled computation; install adoption-test packages only if requested; add the Agilent PBMC branch only if prioritised. Paper_2 needs GSE286184, GSE283826, GSE341102, GSE199850 and GSE193512; Paper_3 needs eQTL/GWAS summary statistics. Keep `Datasets/manifest.tsv` current.

## 5. Pitfalls & Anti-patterns (Never repeat these)

### 5.1 Claim limits that must not be overstated

- The gate beats the panel in ERCB; never claim panel superiority.
- The strict refit-permutation null is not significant (*p* = 0.11). Report both nulls.
- Array-to-array transfer is the hard direction. Never restore cohort-internal z-scoring to recover old AUCs.
- Tubulointerstitial discrimination is chance (0.449/0.539); there is no TI MN module.
- FSGS is the hard comparator (0.632). The multiclass balanced accuracy is 0.44, so this is MN detection, not a diagnostic multiclassifier.
- MYOZ2 has tissue-level reproducibility with unresolved cellular source, not a proven podocyte or SMC origin.
- SLAMF8 detection is only 9.5% in RNA-seq and 27% in LCM; report the three-gene sensitivity model.
- Deconvolution is fragile; use the filtered signature as primary, all-gene as sensitivity, and exclude pDC/neutrophil from traits and attribution.
- Subtypes are a negative result; the array split tracks technical quality.
- Bulk↔pseudobulk concordance is weak (best Pearson 0.29); pseudobulk DE is hypothesis-grade.
- GSE197307/GSE182380 have no participant linkage. Never infer overlap from titles; use `id_type` and `id_is_globally_unique`.
- Within-MN PLA2R (6 versus 6) is uninterpretable; do not claim PLA2R tracking.
- Expanded Track B is post hoc. Do not present its median R5 0.674 as preregistered confirmation.
- Do not claim a systematic PubMed/GEO search or an unverified EBioMedicine cap.

### 5.2 Technical pitfalls to carry forward

- **P-1** Compute patient overlap across all releases before merging any matrix.
- **P-2** Never infer cross-cohort overlap from sample titles. Check participant fields, local 1..n indices and per-disease overlap tables; absent linkage means independent cohorts.
- **P-3** Use `oligo` + Brainarray indices for HuGene-2.1-ST, not `affy`.
- **P-4** `preprocessCore::rcModelPLM()$Estimates` is not the sample-effect vector; use `stats::medpolish()` overall + column effects.
- **P-5** Never `match()` each gene against the 1.02M-row probe index; build one lookup vector.
- **P-6** Chunk large ST cohorts (`CHUNK = 40`) instead of processing 244 arrays at once.
- **P-7** CDF package names must end in `cdf`.
- **P-8** `make.cdf.package()` builds a source package; install it and pass a bare filename after `setwd()`.
- **P-9** Read CEL headers for platform identity, not series-matrix metadata shared across platforms.
- **P-10** Align platform gene unions with NA fill and drop incomplete genes before PCA and limma; never naively `rbind()`.
- **P-11** Never redirect stdout to the file a script writes with `writeLines()`.
- **P-12** Count-style filenames may contain normalised decimals and text columns, not raw counts.
- **P-13** Namespace living donors as `LD-<n>`; never use raw donor integers as globally unique IDs.
- **P-14** Data acquisition is manual: never download datasets; present a request list. Package/MSigDB annotation fetching is separate.
- **P-15** On `~0 + group`, use `makeContrasts()` + `contrasts.fit()`, not `topTable(coef = "MN - Control")`.
- **P-16** An `mclapply` error returns `try-error`, not NULL; filter valid one-row data frames and wrap model calls in `tryCatch`.
- **P-17** `metafor::rma(method = "REML")` may not converge; wrap and count DL fallbacks.
- **P-18** Use `fgsea::plotEnrichment`, not `enrichplot::gseaplot2()` for fgsea results.
- **P-19** Parse prefetched GMTs in `work/msigdb/`; avoid first-use `msigdbr` downloads.
- **P-20** Scope verifier figure assertions to the generating step’s filename prefix.
- **P-21** Foreground shells cap at two minutes; run long jobs in the background and inspect `.out` logs.
- **P-22** ERCB codes MN as `MGN`; `pROC::roc()` needs a two-level response.
- **P-23** Set centroid-matrix `colnames()` explicitly after `sapply()`.
- **P-24** PAC ≤ 0.1 is the subtype stability gate, not silhouette.
- **P-25** Stable splits can be technical; test correlation, mean expression and detection rate against clusters.
- **P-26** Precompute counts for text labels; `geom_text(stat = "count")` rejects `fill`.
- **P-27** Missing string-replacement anchors fail silently; assert both anchor and assembled result.
- **P-28** `file.exists()` accepts zero-byte files; check `file.size() > 0` and run `find <dir> -size 0` after crashes.
- **P-29** Call `JoinLayers()` immediately after Seurat 5 `merge()`; compute AddModuleScore controls below the smallest atlas bin.
- **P-30** Prefix barcodes before merging and rebuild the same key downstream to avoid `_1`, `_2` collisions.
- **P-31** Align named vectors on the union and fill zero before `cbind()`; never align unequal vectors by position.
- **P-32** NNLS needs correct scaling: joint quantile normalisation, unit sums, τ > 0.5 filtering, and τ computed on CPM before column normalisation.
- **P-33** `[.DGEList` has no `drop`; Seurat `ncol()` can be double; all-zero Wilcoxon groups return NaN and must be marked untestable.
- **P-34** exfat has no hardlinks/symlinks and is case-insensitive. Never use `ln -f`; use `mv` within the volume and stage link-atomic writers through safe paths.
- **P-35** WGCNA 1.74 `modulePreservation` requires numeric reference indices, two-level `multiData`, `multiColor`, `.qual` output and unambiguous module-size mapping.
- **P-36** `bicorAndPvalue()` returns `$bicor`; `mad()` can be zero with >50% identical values, so record the correlation method per pair.
- **P-37** Name module colour vectors by gene, never by sample.
- **P-38** Split-half tests must reselect soft power in each half and report Jaccard separately from eigengene bicor.
- **P-39** With numeric labels, normalise `kME_ME1` names and assert no own-gene kME is NA.
- **P-40** Build full `meta_<Glom|TI>_<contrast>` attribution keys and print join counts.
- **P-41** Require ≥95% ERCB coverage at panel selection and assert it; never drop nonmeasurable genes post hoc.
- **P-42** Do not report Youden sens/spec for a low-prevalence cohort; use a fixed 0.5 threshold, calibration and slope, and describe discrimination rather than a risk instrument.
- **P-43** Summarise each permutation null from its own distribution, including `null_max`.
- **P-44** “No case observation” can mean NA predictors; check `sum(is.na(score))` first.
- **P-45** Use row-wise summaries on genes × samples matrices; an exact AUC 0.500 can indicate a constant design.
- **P-46** GSE99339 is split by platform and codes MN as `MGN`; inspect filename tokens when parsing yields zero cases.
- **P-47** Avoid `merge()` collisions such as `note.x/note.y` by naming derived columns distinctly.
- **P-48** GEO characteristic lines are per-field vectors; match line length to sample count and assign by position. `!Sample_description` can contain submitter IDs.
- **P-49** Ghostscript cannot read `/run/media`; stage through `/tmp`, embed all fonts and verify with `pdffonts`.
- **P-50** Publisher PDFs can be Cloudflare-403; use live open-access JATS XML and never quote an unread cap.
- **P-51** Never hand-type citations; resolve DOI/PMID and assert title. `Manuscript/references.tsv` is the only source.
- **P-52** Assert on assembled prose, not just on successful script completion.
- **P-53** Draw each nested-CV partition once per repeat and prove test folds are disjoint, rather than checking only AUC.
- **P-54** `tryCatch(error=)` misses metafor warnings; use `withCallingHandlers`, record per-gene flags and report counts.
- **P-55** Do not quantile-normalise a signature with its mixtures; it destroys the specificity retained by τ filtering.
- **P-56** Add a new column to every branch of a per-gene builder and assert identical column sets.
- **P-57** `reshape()` varying columns must come from an explicit allow-list.
- **P-58** Fixing transductive standardisation can lower external AUC; never revert to cohort-internal scaling.
- **P-59** Drivers must resolve `ROOT="${PAPER1_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"` and be tested from an unrelated cwd.
- **P-60** Derive labels and filenames from data (for example `length(panel)`), never hard-code panel size.
- **P-61** Verify composed figures by measurement: bbox sizing, trailing blank bands, ink box/MediaBox, title text and data-derived windows.
- **P-62** Training statistics alone are not cross-platform transfer; first check |z| median, then label-blind quantile-map before training μ/σ.
- **P-63** A rerun is not re-derivation: grep for private copies of changed scaling rules and reconcile duplicate quantities.
- **P-64** Content-driven page sizing can create a second page; assert one page and delete `<plate>-N.jpg` leftovers.
- **P-65** Provenance keys require `id_type` and `id_is_globally_unique`, asserted in verifiers.
- **P-66** Re-derive counts from the filesystem whenever artefacts change and state “as of <date>”.
- **P-67** Calibrate checkers against known values and parse CSVs with a real reader; `pdffonts` embedding is `$(NF-4)`.
- **P-68** Never type a number into a script when an artefact supplies it; read it and assert reproduction.
- **P-69** `check(isTRUE(x))` fails on vector predicates; wrap with `all()` or `any()`.
- **P-70** A literature registry needs code-enforced inclusion (≥2 named genes, PMID/DOI, reported performance) and a reason for every exclusion.
- **P-71** `ggsave(device = "jpeg")` fails on R 4.6.1; use `device = grDevices::jpeg`. `tools/make_figures.py` is unaffected.
- **P-72** Never infer a panel’s content from its filename or an older legend. Inspect title, axes and source table before citing: meta versus discovery, top-40 versus full set, discovery versus all sets, centroid versus overlap, and blocks versus stable subtypes.
- **P-73** Correct code does not prove current output. After label changes, rerun and inspect regenerated CSVs and figures; the stale “all 5 genes” artefact was corrected by rerunning `09d` and its verifier.
- **P-74** Main figures must defend claims, not merely show routine QC. Keep design, independence and harmonisation in main Fig 1 and routine PCA in the supplement. When inserting a supplementary plate, renumber downstream references in descending order, update Results and legends together, and delete superseded files.

## 6. Key Resources & Environment Info

### 6.1 Repository and path configuration

- Branch: `master`; git log is authoritative. `core.fileMode` and `core.autocrlf` are `false` for exfat.
- Project root: `/run/media/meng/MY_Data/MN/Paper_1`.
- Data root: `/run/media/meng/MY_Data/MN/Datasets`, shared with Paper_2/Paper_3.
- `R/_lib.R` resolves the root and data paths. Overrides: `PAPER1_PROJ` and `PAPER1_DATA`.
- The mount can alternate `/run/media/{meng,yang}/MY_Data`; keep scripts path-independent.
- `work/`, `.Rlib/`, `scrna/extracted/`, `logs/` and five bulk CSVs >3 MB are gitignored as regenerable. A checkout alone is not a reproducible run.
- No API keys or credentials are used. CIBERSORTx would require a licence token only if requested as an optional cross-check.

### 6.2 Languages and packages

- **R 4.6.1**, Bioconductor 3.23.
- **Python 3.14.4**; plotting/composition uses the project’s Python + LaTeX/TikZ workflow rather than matplotlib.
- **Pandoc 3.8.3**, **OfficeCLI 1.0.144**, TeX Live, Ghostscript and LibreOffice are used for manuscript/figure production and QA.
- Two R libraries:
  - `<PROJ>/.Rlib`: makecdfenv, five Brainarray CDFs, oligo, WGCNA 1.74, ConsensusClusterPlus, caret, glmnet, randomForest, xgboost, Boruta, metafor 5.2.1, msigdbr, rms, babelgene.
  - `/home/meng/R/x86_64-pc-linux-gnu-library/4.6`: ~1,701 packages, including Seurat 5.6.0, scDblFinder 1.26.9, harmony 2.0.5, muscat, nnls, limma, edgeR, DESeq2, sva, clusterProfiler, fgsea, ComplexHeatmap and survival.
- Installed optional packages: BayesPrism, MuSiC, hdWGCNA.
- Absent: ReactomePA, immunedeconv/xCell/quanTIseq, msigdbdf, liana, decoupleR, InstaPrism, tabpfn.
- Mirrors: TUNA CRAN/Bioconductor; China-network preference is ghproxy/Aliyun where applicable.

### 6.3 Reproducibility record and restoration

- Lockfile: `environment/renv.lock` (1,669 packages).
- Environment notes: `environment/README.md`, `environment/sessionInfo.txt`, `environment/packages_*.csv`.
- Regenerate the record with `Rscript tools/make_environment_record.R`.
- `renv` is a reproducibility record and is not activated at runtime.
- Restore on another machine:

```r
install.packages("renv")
renv::restore(lockfile = "environment/renv.lock")
```

- Rebuild local Brainarray CDFs from GEO archives (no network needed):

```bash
Rscript R/00d_install_cdf.R
Rscript R/02a_build_cdf.R
for p in work/cdf/*; do R CMD INSTALL -l .Rlib "$p"; done
export PAPER1_DATA=/path/to/MN/Datasets
```

- Cross-machine copy: include at least `R/`, `metadata/`, `output/`, `reference/`, `Manuscript/`, `Docs/` and the `.rds` files under `work/`, plus `MN/Datasets`. Install R 4.6.x and Bioconductor 3.23, restore packages and rebuild `.Rlib` as above. Regenerate `scrna/extracted/` with `Rscript R/07a_scrna_manifest.R`. Install TeX Live + Ghostscript for figures, then run all seven verifiers.

### 6.4 Compute resources and key paths

- Linux amd64, 20 threads, 61 GB RAM, 1.8 TB free; local workstation only, no HPC job files.
- Project ~32 GB: `work/` 23 GB, `.Rlib/` 4.2 GB, `scrna/` 3.3 GB, `output/` ~650 MB; shared `Datasets` 18 GB; exfat.
- `work/cel/` (16 GB) and `work/rma_cache/` (1.1 GB) may be deleted once `work/rma/` is final.

| Path | Contents |
|---|---|
| `$DATA` | Shared GEO collection: `<GSE>_RAW.tar`, `<GSE>/{raw,processed}/`, `manifest.tsv` |
| `work/sets/` | Seven cohort sets: `*_expr.rds`, `*_expr_symbol.rds`, `*_samples.csv` |
| `work/rma/`, `work/rnaseq/`, `work/combat/`, `work/msigdb/` | Normalised cohorts, ComBat TI matrices, MSigDB v2025.1.Hs GMTs |
| `work/scrna/` | QC objects, 236k-cell atlas, donor × gene pseudobulk matrices |
| `work/wgcna/` | Networks, eigengenes, traits, preservation objects; run `08d` from project root |
| `work/bench/10b_harmonised.rds` | Track B `z_harmonised` primary and `z_transductive` comparison matrices |
| `output/tables/` and subdirectories | 240 result CSVs, 121 at top level; current numbers live here |
| `output/figures/`, `output/figures/supplementary/` | 13 plates, PDF + 300-dpi JPEG; LaTeX in `work/figtex/` |
| `output/models_*`, `output/{qc,deg,enrichment,wgcna,scrna,subtypes}/` | Per-step PDF + JPEG figures |
| `metadata/` | Track B register/exclusions, PLA2R table, Ensembl-symbol map |
| `Manuscript/` | Manuscript, cover letter, IF plan, data requests, checklist, `references.tsv` |
| `Docs/` | Consolidated results, BAP audit, findings ledger, MYOZ2 provenance, superseded handoff |
| `R/_pre_audit_backup/` | Pre-audit scripts for provenance only |

### 6.5 Execution map

| Scripts | Role |
|---|---|
| `03d`, `04a`, `04b`, `04c*`, `04d*` | QC policy; ComBat; limma; random-effects meta + ERCB; ORA/GSEA |
| `05a`, `05b`, `05c*` | Separability gate, robustness, partition stability |
| `06a` | Exploratory negative subtyping |
| `07a`, `07a0`, `07b*`, `07c*`, `07d*`, `07e`, `07f`, `07z` | scRNA manifest, QC, atlas, pseudobulk, bridge, deconvolution, packages |
| `08a`, `08b`, `08c*`, `08d*`, `08e`, `08f*` | WGCNA, traits, stability, preservation, hubs, enrichment |
| `09a`–`09f` | Panel selection, validation, robustness, mechanism, sensitivity, PLA2R |
| `10b`, `10c`, `10d` | Track B matrices, attrition ladder, panel-size frontier |
| `99_verify_*` | 54 / 86 / 32 / 42 / 117 / 106 / 37 assertions |
| `tools/make_cohort_design.py` | Data-derived main-Fig-1 panels from existing metadata/result tables |
| `tools/make_figures.py` | Plate compositor; automatically calls the cohort-design builder |
| `tools/moyz2_provenance.R` | Separate MYOZ2 provenance test, ~4 min cold, cached thereafter |
| `tools/make_environment_record.R`, `tools/rerun_*.sh` | Environment record and resume drivers |

`*` means run in background (>2 min), especially `04c`, `04d`, `07b`, `07c`, `07d`, `08c`, `08d`, `08f`, `05c` and `10d`.

### 6.6 Rebuild and verification commands

- Panel revisit: `09a → 09b → 09c → 09d → 09e → 09f`, then `python3 tools/make_figures.py`.
- WGCNA revisit: `08a → 08b → 08c → 08d → 08e` (`08f --force`), or `bash tools/rerun_from_08d.sh`.
- Track B: `10b → 10c → 10d → 99_verify_trackb.R`.
- MYOZ2: `Rscript tools/moyz2_provenance.R` (`--refresh` rebuilds the cell cache).
- Figure rebuild: `python3 tools/make_figures.py`; it invokes `make_cohort_design.py`, composites with `pdflatex`, embeds fonts with Ghostscript through `/tmp`, and writes PDF + JPEG. Assert one page per plate and run the 5-point glyph audit.
- Submission copies: remove old `Manuscript/Figures/*.pdf` and `*.jpeg`, then copy `output/figures/Fig*` and `output/figures/supplementary/FigS*` into `Manuscript/Figures/` (26 files).
- Manuscript rebuild:

```bash
pandoc Manuscript/Paper1_manuscript.md \
  --reference-doc=reference/tnr_fangsong_redcode.docx \
  -o Manuscript/Paper1_manuscript.docx
```

- After every Pandoc rebuild, reapply the OfficeCLI finishing steps: normalise `VerbatimChar` (Consolas/FangSong, 11 pt, red) so `rFonts` precedes colour/size; add `TableNormal`; add a centred live PAGE footer; set `updateFields=true`; format paragraphs after the Figure legends heading (currently `/body/p[150]`–`[162]`) at 10 pt, 3 pt before/after, 1× spacing; save and close.
- Final gates: `officecli validate` = no errors; `officecli view … issues --type structure` = 0; LibreOffice render = 18 pages with no blank/orphan page; plain-text check = 13 legends, all 44 references cited, all 65 panel citations sequential.
- `handoff.docx` is stale: `handoff.md` is the only checkpoint source; do not read or regenerate the DOCX.

## 7. Startup Instruction for the Next Session

> **When starting a new session, the very first message must be: "Read handoff.md first."**

1. Restate the current state from this document before making changes. If a newer project file postdates this checkpoint, first determine whether it represents unrecorded work.
2. Confirm the frozen state with all seven verifiers, expecting **54 / 86 / 32 / 42 / 117 / 106 / 37**, all exit 0. Check for zero-byte outputs and re-derive any count before quoting it.
3. Do not quote pre-audit results, do not add methods to Paper 1, do not hard-code values already stored in result tables, and do not treat `gsm` as a GEO accession for RNA-seq cohorts.
4. B15 is resolved at 12 panels. The live decisions are whether Track B ships at 12 and the P1 submission metadata (authors, funding, reviewers and declarations). P0 IF work and P0b clinical-data requests are externally blocked. Ask which item the user wants to proceed with rather than opening new analyses automatically.
5. For persistence, commit this `handoff.md` with the completed manuscript/figure changes or copy it to a backed-up location. Do not commit generated work without reviewing the diff first.
