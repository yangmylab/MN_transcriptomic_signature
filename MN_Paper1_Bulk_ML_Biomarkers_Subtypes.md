# Paper 1 — Machine-Learning Biomarker Discovery & Molecular Subtyping of Membranous Nephropathy from Multi-Cohort Bulk Transcriptomics

**Working title:** *Machine-learning–guided deconvolution of glomerular and tubulointerstitial transcriptomes reveals immune–complement biomarkers and molecular subtypes of membranous nephropathy*

**Document purpose:** Independent analysis strategy + detailed pipeline for high-impact SCI Paper 1. All dataset accessions below were verified against NCBI GEO (2026-10-15).

---

## 1. Research Question & Rationale

Membranous nephropathy (MN) is an antibody-mediated glomerular disease (anti-PLA2R in ~70–80%, also THSD7A/NELL1/SEMA3B/EXT) with heterogeneous outcomes: one-third spontaneous remission, one-third persistent proteinuria, one-third progressive kidney failure. Current prognostication relies on proteinuria trajectory and anti-PLA2R titers — no tissue-molecular biomarker is clinically used.

**Central hypothesis:** The glomerular and tubulointerstitial (TI) transcriptomes of MN contain (a) a compact, compartment-resolved gene signature that robustly separates MN from other proteinuric glomerulopathies (FSGS, MCD, IgAN, LN, ANCA-GN) and healthy tissue, and (b) reproducible molecular subtypes (e.g., immune/complement-high vs. fibrosis/metabolic-high) with distinct prognostic trajectories and druggable pathway profiles.

**Why this is publishable at high impact:**
- Largest assembled MN bulk-transcriptomic resource to date: **>500 MN samples** across 8+ independent cohorts, two micro-dissected compartments, both microarray and RNA-seq — no prior study has harmonized all of them.
- Clinical anchoring: NEPTUNE cohorts carry longitudinal proteinuria/eGFR/remission data → biomarkers tied to *outcomes*, not just diagnosis.
- Cross-platform, cross-compartment design is intrinsically robust to the "single-cohort overfitting" criticism that kills most biomarker papers.

---

## 2. Datasets (all verified public, GEO)

| Accession | Source cohort | Compartment | Platform | MN n | Comparator groups | Role |
|---|---|---|---|---|---|---|
| GSE133288 | NEPTUNE | Tubulointerstitium | Affy HuGene 2.1 ST (GPL19983) | 48 | 69 FSGS, 54 MCD, 68 other, 5 donors | **Discovery (TI)** |
| GSE200828 | NEPTUNE screening | Glomeruli | GPL19983 | 51 | 46 other, 31 FSGS, 19 MCD, 6 donors | **Discovery (glom)** |
| GSE200818 | NEPTUNE screening | Tubulointerstitium | GPL19983 | 48 | 56 other, 54 FSGS, 30 MCD, 5 donors | Discovery support |
| GSE108109 | ERCB/NEPTUNE | Glomeruli | GPL19983 | 44 | 30 FSGS, 16 MCD, 15 ANCA, 6 donors | Validation |
| GSE108112 | ERCB/NEPTUNE | Tubulointerstitium | GPL19983 | 43 | 57 ANCA, 46 FSGS, 18 MCD, 5 donors | Validation |
| GSE104948 | ERCB | Glomeruli | Affy U133 Plus2/A (GPL22945/24120) | 21 | SLE 32, IgAN 27, FSGS 24, ANCA 22, MCD 20, DN 15, controls | External validation |
| GSE104954 | ERCB | Tubulointerstitium | same | 17 | similar mix | External validation |
| GSE99339 | ERCB (hypoxia) | Glomeruli | Brainarray (GPL19109/19184) | 21 | SLE 30, IgAN 26, RPGN 23, FSGS 22, MCD 19, DN 14 | External validation |
| GSE99325 | ERCB (hypoxia) | Tubulointerstitium | same | 18 | similar mix | External validation |
| GSE197307 | NEPTUNE + donors | Glomeruli | NextSeq 500 **RNA-seq** (GPL18573) | 62 | FSGS 93, MCD 89, other 29, **8 donors** | **RNA-seq validation (glom)** — has JASN 2022 MN paper (PMID 35477557) |
| GSE182380 | NEPTUNE | Tubulointerstitium | NextSeq 500 **RNA-seq** | 59 | FSGS 70, MCD 32, other 37 | RNA-seq validation (TI) |
| GSE216841 | LCM-FFPE | Glomeruli | HiSeq 4000 RNA-seq | 12 (6 PLA2R+, 6 PLA2R−) | MCD 14, normal 8 | PLA2R-status sub-analysis |
| GSE73953 | — | PBMC | Agilent 4x44K | 8 | IgAN 15, healthy 2 | Peripheral cross-check (optional) |

⚠️ Do **not** use: GSE32591 (lupus nephritis, not MN), GSE104066 (no MN samples), GSE99340/GSE47185 (SuperSeries wrappers only).

**Clinical variables (NEPTUNE):** baseline eGFR, UPCR, and longitudinal remission/relapse status — request/download via NEPTUNE data sharing or use published per-sample annotations.

---

## 3. Analysis Pipeline (step by step)

### Phase A — Data assembly & preprocessing (Week 1–2)

1. Download raw CEL files (microarray: GSE133288/200828/200818/108109/108112/104948/104954/99339/99325) and raw count matrices (RNA-seq: GSE197307/182380/216841) via GEOquery/`wget` from `ftp.ncbi.nlm.nih.gov/geo/series/`.
2. **Microarray (Affymetrix):** RMA normalization per platform family using `oligo::rma()` (HuGene 2.1 ST) and `affy::rma()` (U133), with Brainarray ENTREZG CDFs for ERCB sets so probes collapse to Entrez IDs.
3. **RNA-seq counts:** filter low counts (edgeR `filterByExpr`), TMM normalization, `voom` for linear modeling; keep the RNA-seq cohorts as an *independent validation layer* (never merged with microarrays).
4. Build two master discovery matrices (glom: GSE200828; TI: GSE133288) and four validation matrices; map all to HGNC symbols; log2-transform where needed.
5. QC: hierarchical clustering + PCA outlier screen, arrayQualityMetrics reports; document every excluded sample.

### Phase B — Batch integration & differential expression (Week 2–4)

6. Within-platform meta-integration: merge GPL19983 cohorts (GSE133288+200818 for TI; GSE200828+108109 for glom) with `sva::ComBat_seq`-equivalent for microarray (ComBat on log-intensities), study as batch, disease as biological covariate. Verify batch removal by PCA silhouette.
7. **DEG analysis** with limma (~trend/eBayes), MN vs. (a) healthy donors and (b) each other glomerulopathy separately, per compartment. Thresholds: |log2FC| > 1, adj. P (BH) < 0.05.
8. **Meta-analysis across cohorts:** combine effect sizes with random-effects model (`metafor`) on the per-cohort limma results → "meta-DEG" list with heterogeneity statistics (I²); this is the reviewer-proof core signature.
9. Functional profiling: clusterProfiler (GO/KEGG/Reactome), GSEA with Hallmark + immunologic (C7) collections, DOSE; expected signals: complement cascade, FcγR-mediated phagocytosis, antigen presentation, ECM organization, hypoxia (validated by GSE99339/99325 design).

### Phase C — Immune/stromal deconvolution (Week 4–5)

10. Estimate cell-type fractions per sample: **CIBERSORTx** (LM22 + a kidney-reference signature built later from Paper 2's single-cell atlas — cite as cross-paper synergy), **xCell**, **MCP-counter**, **quanTIseq**. Compare MN vs. controls/other GNs (Wilcoxon, BH-adjusted).
11. Correlate infiltrates with eGFR/UPCR (Spearman) in NEPTUNE samples; test association with remission status (logistic regression adjusted for age/sex).

### Phase D — WGCNA & module–trait mapping (Week 5–7)

12. Signed WGCNA per compartment (soft-power by scale-free fit R² > 0.85; dynamic tree cut, mergeCutHeight 0.25).
13. Module eigengenes ↔ traits: MN status, eGFR, UPCR, remission, deconvolution fractions. Prioritize modules significant in *both* discovery cohorts.
14. Module preservation analysis (Zsummary > 10) across validation cohorts → only preserved modules proceed. Hub genes = |kME| > 0.8.

### Phase E — Machine-learning biomarker panel (Week 7–10)

15. Candidate space: intersection of meta-DEGs ∩ preserved WGCNA hub genes (expect ~50–150 genes).
16. Three orthogonal selectors on the discovery glom+TI sets (nested 5-fold CV, repeated 10×):
    - **LASSO** (glmnet, λ.1se),
    - **SVM-RFE** (caret/e1071),
    - **Boruta** (all-relevant, random forest),
    - keep genes selected by ≥2 methods → core panel (target 4–10 genes).
17. Final classifier: logistic regression + random forest + XGBoost compared; pick by cross-validated AUC on discovery; **locked model** then evaluated once on each validation cohort (microarray *and* RNA-seq layers).
18. Metrics: AUC with DeLong CIs, calibration curves, decision-curve analysis; also one-vs-rest performance for MN vs. each other GN (clinically meaningful differential diagnosis).
19. Build a nomogram (rms) combining the gene score + UPCR + anti-PLA2R status (where available) vs. gene score alone — quantify incremental value (IDI/NRI).

### Phase F — Molecular subtyping (Week 10–12)

20. Consensus clustering (ConsensusClusterPlus, k = 2–6, 1,000 resamples) of MN samples only, on the meta-DEG space, separately per compartment; choose k by consensus CDF/Δ-area and silhouette.
21. Characterize subtypes: pathway scores (GSVA Hallmark/Reactome), deconvolution profiles, clinical outcomes (remission-free survival, Cox PH; eGFR slope, linear mixed models), and PLA2R+/− distribution (GSE216841).
22. Reproduce subtypes in GSE197307/182380 (RNA-seq) using a nearest-centroid classifier trained on discovery centroids.
23. **Subtype-directed drug repositioning:** query each subtype's up/down signature against LINCS L1000 (signatureSearch / cmapR) and DSigDB/Enrichr; rank perturbagens that reverse the subtype signature; flag known MN-relevant agents (rituximab targets, calcineurin inhibitors, complement inhibitors, anti-CD38) as positive controls.

### Phase G — Validation & translational wrap-up (Week 12–14)

24. Literature/omics cross-validation: overlap panel genes with Xie 2020 GWAS loci (Paper 3 synergy), with the JASN 2022 NEPTUNE glomerular MN signature (PMID 35477557), and with urinary biomarker literature.
25. Orthogonal validation plan (for the discussion + future wet-lab): qPCR/IHC of top panel genes on an in-house MN biopsy cohort; staining quantification vs. transcript score.
26. Reproducibility package: renv lockfile, scripted end-to-end pipeline, all random seeds fixed, figures generated by code only.

---

## 4. Statistical & Reporting Standards

- Multiplicity: BH-FDR < 0.05 throughout; effect sizes always reported with 95% CI.
- Class imbalance: report PR-AUC alongside ROC-AUC for MN-vs-each-GN tasks.
- No test-set peeking: validation cohorts touched exactly once with the locked model.
- Reporting: follow TRIPOD for the diagnostic model; STROBE for subtype–outcome analyses.

## 5. Planned Figures (6 main + supplement)

1. **Fig 1** — Study design flowchart (cohorts, compartments, phases) + PCA before/after batch correction.
2. **Fig 2** — Meta-DEG landscape: volcano (glom & TI), cross-compartment UpSet, GSEA hallmark ridge plots.
3. **Fig 3** — Immune/stromal deconvolution heatmaps + correlations with eGFR/UPCR.
4. **Fig 4** — WGCNA: dendrogram, module–trait heatmap, preservation Zsummary, hub-gene networks.
5. **Fig 5** — ML panel: feature-selection Venn, cross-validated ROCs (discovery + 4 validations), nomogram, decision curves.
6. **Fig 6** — Subtypes: consensus matrices, Kaplan–Meier/eGFR-slope by subtype, LINCS drug-reversal bubble plot.

## 6. Risks & Mitigations

| Risk | Mitigation |
|---|---|
| Small healthy-donor n in some cohorts | Use ERCB tumor-nephrectomy controls; sensitivity analysis excluding tiny cohorts |
| Platform heterogeneity | Two-layer design: microarray discovery → RNA-seq validation; never pool raw values across platforms |
| Biomarker overfitting | Nested CV, locked model, single-touch validation, TRIPOD |
| Subtypes driven by batch | Subtype discovery only after ComBat; require RNA-seq reproduction |

## 7. Target Journals & Timeline

- **Primary targets:** *Kidney International*, *JASN*, *Clinical Journal of the ASN*, *EBioMedicine*, *Briefings in Bioinformatics*.
- **Timeline:** 14 weeks analysis + 4–6 weeks writing → submission ready in ~4.5–5 months.

---

*Author: analysis strategy prepared for Prof. Lingzhang Meng, 2026-10-15. Synergies: deconvolution signature from Paper 2 (single-cell atlas); GWAS overlap from Paper 3 (MR).*
