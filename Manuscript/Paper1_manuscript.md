---
title: "A compartment-resolved glomerular transcriptomic signature and a four-gene classifier for membranous nephropathy across three independent cohorts"
author:
  - "Lingzhang Meng^1,2,\\*"
  - "[further authors to be confirmed before submission]"
affiliation:
  - "^1 Research Center for Medical Big Data and Artificial Intelligence, The People's Hospital of Guangxi Zhuang Autonomous Region, Nanning, China"
  - "^2 Guangxi Academy of Medical Sciences, Nanning, China"
  - "^\\* Correspondence: Lingzhang Meng, Research Center for Medical Big Data and Artificial Intelligence, The People's Hospital of Guangxi Zhuang Autonomous Region, No. 6 Taoyuan Road, Nanning 530021, China (m.lingzhang@hotmail.com; lzmeng@gxams.org.cn)"
date: "Draft — 2026-11-19"
---


# Summary

**Background.** Membranous nephropathy (MN) is an antibody-mediated glomerular disease diagnosed by
biopsy immunofluorescence and phospholipase A2 receptor (PLA2R) serology. No tissue transcriptomic
biomarker is in clinical use, and it is unclear whether the glomerular transcriptome holds a compact
signature separating MN from the proteinuric glomerulopathies it mimics.

**Methods.** We harmonised 13 public GEO accessions (2,165 samples; 453 MN samples from 178 unique
patients) into seven non-redundant cohort sets: NEPTUNE microarrays for discovery, the European Renal
cDNA Bank (ERCB) for independent external validation (no shared participants) and two RNA-seq cohorts
for cross-platform replication. Compartments were analysed separately, using limma and
random-effects meta-analysis, an elastic-net separability gate, signed WGCNA with ERCB preservation
testing, and a 236,107-cell kidney atlas for cell-type attribution and NNLS deconvolution. A panel was selected by LASSO, SVM-RFE and Boruta inside nested 10 × 5-fold
cross-validation and locked before any validation data were used.

**Findings.** Discovery separated MN from healthy glomeruli (584 meta-differentially expressed genes)
but barely from other glomerulopathies at threshold (18 genes versus FSGS), yet discrimination was
strong (nested-CV AUC 0.91–0.93; external ERCB AUC 0.86–0.92 once the validation cohort is
harmonised to the discovery reference distribution and scored with training coefficients only;
permutation *p* = 0.048). Two small glomerular modules
carried the signal — podocyte surface/calcium-binding (84 genes; adjusted *p* = 1.0 × 10⁻⁵) and
lysosomal/glycosylation (144 genes; adjusted *p* = 6.6 × 10⁻⁴) — and both were preserved in ERCB
(Z = 21–85). A four-gene panel (*SLAMF8*, *TRPC6*, *MYOZ2*, *PPP1R3C*) reached AUC 0.816
(95 % CI 0.721–0.910) in ERCB, 0.901 in RNA-seq, 0.929 in a laser-capture cohort and 0.848 across
array annotation (label-permutation *p* = 0.005). It was compartment-specific (tubulointerstitial
0.45/0.54), both PLA2R strata were detected against MCD (0.98/0.88), and a two-parameter
recalibration cut the Brier score from 0.845 to 0.096 at unchanged AUC. Exploratory subtyping was
negative.

**Interpretation.** The MN glomerular transcriptome contains a reproducible, ERCB-preserved podocyte
surface/calcium and lysosomal/glycosylation signature, from which a four-gene classifier approaches —
without matching — the external discrimination of a 13,379-gene model (0.82 versus 0.92 in ERCB). It
is glomerular-compartment-specific and a discrimination rather than a calibrated risk
instrument; no tubulointerstitial or subtype-based claim is supported.

**Funding.** No specific funding; the author had full access to all data and final responsibility for
the decision to submit.

**Keywords:** membranous nephropathy; transcriptomics; biomarker panel; WGCNA; machine learning;
external validation; single-cell deconvolution

# Research in context

**Evidence before this study.** The evidence base for this study is the set of public MN expression
cohorts deposited in GEO together with the literature cited in this paper (references retrieved from
PubMed/Europe PMC, October 2026). Published MN cohorts are small at the glomerular level, several are
re-releases of the same NEPTUNE patients, and proposed marker sets have almost always been fitted and
evaluated in the same cohort. Thresholded lists of differentially
expressed genes against FSGS are short, which has been read as evidence that MN is transcriptomically
close to its mimics.

**Added value of this study.** By resolving patient-level redundancy across every public cohort
before any modelling, then separating discovery, external microarray validation and independent
RNA-seq replication, we show that MN is separable from its mimics far more strongly than differential
expression suggests, that the signal is confined to the glomerular compartment and to two preserved,
composition-adjusted co-expression modules, and that four genes retain external discrimination
approaching that of a 13,379-gene classifier. We also report, rather than suppress, a negative subtyping
result and a miscalibration that requires per-platform correction.

**Implications of all the available evidence.** A compact glomerular classifier for MN is achievable
from public data, but its translation requires protein-level confirmation, per-platform
recalibration, and clinical annotation that current public cohorts do not carry. Negative
"MN-versus-FSGS" differential-expression results should not be read as evidence that the diseases are
transcriptomically indistinguishable.

# Introduction

Membranous nephropathy is the commonest cause of nephrotic syndrome in white adults and a leading
cause of end-stage kidney disease among the glomerulonephritides [1,2]. The identification of PLA2R1 as
the target antigen in roughly 70–80 % of primary cases transformed its serological diagnosis [3,4], and
thrombospondin type-1 domain-containing 7A (THSD7A) and several minor antigens account for part of
the remainder [5,6,7,8,9]. Diagnosis nevertheless still depends on a kidney biopsy read by
immunofluorescence and electron microscopy [10,11], and the biopsy report is morphological rather than
molecular: it does not say which pathways are active, and it does not predict which patients will
remit spontaneously, respond to rituximab, or progress [12,13].

Transcriptomic profiling of glomerular tissue has been used for over a decade to describe MN
biology, and several groups have proposed marker sets [14,15,16]. Three problems recur. First, public MN
cohorts overlap heavily — the same NEPTUNE patients have been re-released in several accessions [17,18],
and treating them as independent inflates apparent replication. Second, cohort sizes at the
glomerular level are small, so thresholded differential-expression lists are unstable. Third, a
list of differentially expressed genes is not a classifier, and a classifier fitted and evaluated
on the same cohort is not evidence.

We therefore asked a deliberately narrow question: does the glomerular transcriptome of MN contain
(i) a compact gene signature that separates MN from the other proteinuric glomerulopathies and
healthy tissue and (ii) reproducible molecular subtypes with distinct biology? To answer it we
first established non-redundancy at patient level across every public cohort we could assemble,
then separated discovery from validation before any modelling, and finally required every claim to
survive an independent cohort and, for modules, an explicit preservation test.

Two negative results shaped the paper. Thresholded differential expression greatly understates the
separability of MN from its mimics, so a machine-learning gate rather than a DEG count was used to
decide whether a biomarker programme was worth pursuing. And exploratory subtyping did not
reproduce, so it is reported as a negative finding rather than omitted.

# Methods

## Cohort assembly, harmonisation and de-duplication

Thirteen GEO accessions were inventoried (2,165 samples). Sample-level metadata were harmonised to
a canonical disease vocabulary (ERCB labels MN as `MGN`), and living donors were namespaced
(`LD-<n>`) because they are numbered independently of patients and otherwise collide with patient
identifiers. Patient-level overlap was computed across every disease class, compartment and
accession before any matrix was merged. Four releases were 100 % contained in another
(`GSE108112` ⊂ `GSE133288`; `GSE108109` ⊂ `GSE200828`; `GSE104948` ≡ `GSE99339`;
`GSE104954` ≡ `GSE99325`) and were dropped or retained only as a normalisation-robustness layer.
The final resource comprised 453 MN samples from **178 unique patients**, with redundancy factors of
2.98× (NEPTUNE-numbered) and 2.00× (ERCB-numbered).

Seven non-redundant cohort sets were analysed separately by compartment: glomerular discovery
(`GSE200828`; 153 samples, 51 MN), tubulointerstitial discovery (`GSE133288` + `GSE200818`; 269
samples, 51 MN), glomerular and tubulointerstitial ERCB validation (`GSE104948`, `GSE104954`; 196
and 195 samples, 21 and 18 MN; **0 % participant overlap with discovery**), two RNA-seq cohorts
(`GSE197307`, `GSE182380`; 282 and 199 samples) and a laser-capture cohort (`GSE216841`; 34
samples, 12 MN). The RNA-seq series carry no participant identifier — their sample titles are local
indices — so they are treated as independent cohorts, never as the same patients measured twice.

## Preprocessing and quality control

Affymetrix data were normalised from CEL files with Brainarray Entrez custom CDFs (v19) built from
the `.cdf.gz` files shipped inside the GEO archives: `oligo` background correction and quantile
normalisation followed by gene-level median polish for HuGene-2.1-ST arrays, and `affy::rma()` for
HG-U133 arrays [19,20,21]. The HuGene-2.1-ST reimplementation was validated against the submitters' own
processed matrix across all 25,582 genes (per-sample Pearson *r* = 0.985–0.987; overall 0.9859).
RNA-seq counts were filtered with `edgeR::filterByExpr`, TMM-normalised and transformed with
`limma::voom` [22,23]. One count matrix required a deterministic repair (every Ensembl identifier had lost
its `G`); the correction is guarded by an assertion.

Batch structure was removed with `sva::ComBat(batch = accession, mod = ~group)` [24] on the only
multi-series discovery set, reducing the batch silhouette from 0.097 to −0.033 and batch R² on
PC1–5 from 0.168 to 0.0003 while leaving the disease silhouette unchanged.

Quality control used PCA on the 2,000 most variable genes plus a sample–sample correlation screen
(flag below median − 3 × MAD). **All 200 flagged samples were retained**, with flags recorded as
annotation only; sensitivity analyses confirmed that removing them does not weaken any result
(the gate's ERCB AUC rose from 0.914 to 0.939 when they were excluded).

## Differential expression, meta-analysis and functional profiling

Per-cohort models used `limma` [25] with `~0 + group`, contrast fits, and `eBayes(trend = TRUE,
robust = TRUE)`; genes were required to be detected in ≥ 25 % of samples and complete. Thresholds
were |log2FC| > 1 and BH < 0.05 [26]. Random-effects meta-analysis used `metafor::rma(method = "REML")`
with a DerSimonian–Laird fallback where Fisher scoring did not converge (the fallback is triggered on
metafor's convergence *warning*, and no gene in this dataset required it); a meta-DEG required
|pooled log2FC| > 1, pooled BH < 0.05 and concordant direction in every contributing cohort. A
contrast supported by a single cohort is labelled as such rather than pooled (the tubulointerstitial
MN-versus-control comparison, k = 1). ERCB
was held out of the primary meta and used to test every meta-DEG.

Enrichment used MSigDB v2025.1.Hs (Hallmark, C7, Reactome, KEGG-legacy) [27] and `org.Hs.eg.db` GO
terms, with the genes actually tested as the universe; GSEA used `fgsea` [28] on pooled z-scores.

## Separability gate

Regularised elastic-net logistic regression (`glmnet` [29], α = 0.5, λ by inner 5-fold CV) was trained on
discovery glomerular samples with features reselected inside every outer fold (K = 50, 100, 200,
500, 1,000, all) and tested blind in ERCB. Before scoring, each validation cohort was harmonised to
the discovery reference distribution by rank-preserving quantile mapping applied gene by gene, and
then standardised with the **training-set** means and standard deviations; the mapping uses no
outcome information, whereas cohort-internal z-scores would let the validation cohort's own
distribution set the model's feature weights. This two-step transfer is necessary because the
discovery and validation cohorts were assayed on different array platforms (HuGene-2.1-ST versus
HG-U133): without it the median gene in ERCB sits at \|z\| ≈ 4 against the discovery standard
deviation instead of ≈ 0.6, the linear predictor saturates, and the comparison measures platform
shift rather than biology (0.64–0.77 across feature sizes, with per-diagnosis AUCs as low as 0.05).
Robustness included a
20-shuffle permutation null, removal of QC-flagged samples, the reverse direction and per-diagnosis
AUCs.

## Co-expression network analysis

Signed networks were built per compartment with `WGCNA` [30] (`networkType = "signed"`,
`TOMType = "signed"`, `corType = "bicor"`, `maxPOutliers = 0.05`, `minModuleSize = 30`,
`deepSplit = 2`, `mergeCutHeight = 0.25`, seed 20261019) on genes detected in ≥ 25 % of samples and
ranked by median absolute deviation (top 8,000; single block). Soft power was the smallest value
attaining scale-free R² ≥ 0.85. Module eigengenes were related to MN status, group and the
deconvolution cell-type fractions by bicor with BH correction, and, separately, by a `limma` fit
adjusted for podocyte, fibroblast and myeloid fractions. Reproducibility was assessed with three
diagnosis-stratified split halves plus a flagged-sample-removal refit, each re-selecting its own
soft power, and scored by gene-membership Jaccard and eigengene correlation. Module preservation
was tested with `modulePreservation` (100 permutations; 50 for the cross-platform test) with the
discovery modules as reference; a module was called preserved at Zsummary > 10 and *q* < 0.05.
Modules were annotated by ORA against GO, Hallmark, Reactome and KEGG.

## Single-cell arm

A kidney single-cell/nucleus atlas (37 snRNA + 12 scRNA samples) was quality-controlled per sample
(MAD-based filtering, percent.mt < 20 %, `scDblFinder` for scRNA), integrated with Harmony [31,32,33], and
annotated into 16 cell types from marker modules (236,107 cells × 2,019 genes). Donor-level
pseudobulk differential expression used `limma-voom` with donor as the unit of replication.
Bulk↔single-cell mapping produced a kidney-derived signature matrix, and NNLS deconvolution was
applied to all seven bulk sets with a cell-type-specificity-filtered signature (primary) and an
all-gene signature (sensitivity) [34,35].

## Biomarker panel selection and validation

Candidate genes were the Phase D core set — genes supported by at least two of four independent
evidence lines (preserved-module hub with |kME| > 0.8, meta-DEG, elastic-net gate feature at a
reduced *K*, cell-type-specific attribution with τ > 0.5) — intersected with genes measurable in
≥ 95 % of ERCB samples. Three selectors (LASSO with λ.1se; linear SVM-RFE with panel size chosen by
inner CV; Boruta [36] with 100 runs) were run **inside** each fold of a 10 × 5 nested cross-validation
(10 repeats × 5 folds, with the fold partition drawn **once per repeat** so that the five test folds of
a repeat are disjoint), with scaling fitted on training folds only. A gene entered the panel if selected by ≥ 2 methods at
frequency ≥ 0.60 across all 50 folds. The resulting model was fitted once on all discovery samples
and evaluated without further tuning in ERCB, the RNA-seq cohort, the laser-capture cohort and the
alternative-CDF reanalysis, with DeLong confidence intervals, permutation nulls (200 shuffles
each), decision-curve analysis and calibration assessment [37,38,39].

## Antigen status

PLA2R status for the laser-capture cohort was recovered from the per-sample `phenotype:`
characteristic in the `GSE216841` series matrix (6 PLA2R-positive and 6 PLA2R-negative of 12 MN
samples), joined to the processed matrix through the submitter sample identifiers recorded in the
sample description field. The panel was evaluated against MCD within each stratum, and
PLA2R-positive versus PLA2R-negative was compared within MN.

## Ethics

This study reanalyses publicly available, de-identified data deposited in GEO; no ethics approval or
informed consent was required.

## Role of funders

No funder had any role in study design, data collection, analysis, interpretation or the decision to
submit.

## Software and reproducibility

Analyses used R 4.6.1 with Bioconductor 3.23 (limma, edgeR, sva, WGCNA, fgsea, clusterProfiler,
Seurat, scDblFinder, harmony, nnls) and the workspace library holding the Brainarray CDF packages.
Every result is regenerated by numbered scripts (`R/00…` → `R/09…`) with six verification scripts
(51 + 86 + 27 + 42 + 117 + 106 assertions) that fail loudly if any reported number changes.

# Results

## Thirteen releases collapse to seven non-redundant cohort sets

Patient-level de-duplication was the first analytical step and changed the design. Of 2,165 samples,
453 were MN, but they came from 178 patients; NEPTUNE-numbered samples were redundant 2.98-fold and
ERCB-numbered 2.00-fold, and four accessions were entirely contained in another. Merging by
accession, as a naïve meta-analysis would, would have counted a substantial fraction of patients
several times.

The retained sets were checked for structure before any modelling (Fig 1). Discovery glomerular
samples separated by diagnosis on the first two principal components (Fig 1A), as did the independent
ERCB glomerular cohort (Fig 1B) and the RNA-seq glomerular cohort (Fig 1C); no single sample
dominated an axis. In the tubulointerstitial compartment, ComBat removed the series effect between
the two discovery series while preserving disease separation (Fig 1D); the same diagnosis-led
structure was present in the ERCB tubulointerstitial set (Fig S1A) and the RNA-seq tubulointerstitial
set (Fig S1B). The laser-capture cohort used for the PLA2R strata contained both PLA2R-positive and
PLA2R-negative MN and separated from MCD on the first component (Fig 1E). After de-duplication,
discovery (array), external validation (ERCB
array, 0 % shared participants) and the RNA-seq cohorts are three genuinely independent sample
sources.

## Thresholded differential expression understates separability

Discovery glomerular MN differed from healthy glomeruli by 904 genes (|log2FC| > 1, BH < 0.05; Fig
2A), and the random-effects meta-analysis of the discovery and RNA-seq cohorts produced **584
glomerular meta-DEGs** (353 up, 231 down) with median I² = 0, i.e. the two platforms agreed (Fig 2B).
Representative forest plots show that agreement gene by gene (Fig 2C, *CCNL1*; Fig 2D, *HBA2*).
Against other glomerulopathies, however, thresholded signal was thin: 18 meta-DEGs versus FSGS (Fig
2E) and 18 versus all non-MN disease combined (Fig 2F), with 60 versus MCD and 16 versus other
glomerulonephritis in the underlying contrasts.

ERCB replicated 55.3 % of glomerular MN-versus-control meta-DEGs (252/456; BH < 0.05, same sign) but
only 32.4 % in the tubulointerstitial compartment, and tubulointerstitial MN-versus-FSGS was null in
every cohort. The glomerular compartment, not the interstitium, carries the discriminating signal.

Given 584 versus 18 genes, we tested the claim the paper actually needed — that MN is
distinguishable from its mimics — rather than inferring it from DEG counts. The elastic-net gate
achieved nested-CV AUC 0.911–0.926, and its ERCB ROC was high across feature counts (Fig 3A).
**External ERCB AUC was 0.861–0.922 across feature sizes (0.922 at the full feature set; 95 % CI
0.862–0.982)**, with 0.928–0.939 in the independent RNA-seq cohort, once each validation cohort was
harmonised to the discovery reference distribution and scored with training coefficients only (Fig
3B). That harmonisation is not cosmetic. Removing the transductive z-scoring alone — the intuitive
"training statistics only" fix — left the model measuring platform shift rather than disease: with
discovery standard deviations applied to raw HG-U133 values the median gene falls at \|z\| ≈ 4
instead of ≈ 0.6, the linear predictor saturates, the ERCB AUC becomes K-dependent (0.638–0.679
across twelve refits at the full feature set) and four comparator diagnoses rank *above* MN (RPGN
0.05, FSGS/MCD 0.01–0.14, FSGS 0.40, DN 0.48). Rank-preserving quantile mapping onto the discovery
reference, which is label-blind because it transfers each gene's marginal distribution, removes that
artefact; an independent re-fit of the same model agrees with the sweep to within 0.01 AUC.
Permutation nulls (median AUC 0.500, 95th percentile 0.674) never reached the observed value
(empirical *p* = 0.048; Fig 3C), and removing QC-flagged samples raised ERCB AUC to 0.939. FSGS and
the mixed FSGS/MCD group remained the hardest comparators (0.83 and 0.71), as expected clinically,
and every single comparator was separated (0.71–0.98; Fig 3D). The multiclass extension did not
resolve diagnoses: at *K* = all it recovered 20 of 21 ERCB MN samples but assigned MN to 50 non-MN
samples, giving a balanced accuracy of 0.44 (Fig 3E), which is the behaviour of a single
MN-detection boundary rather than a diagnostic classifier. Refitting the same model under ten
independent fold partitions gave ERCB AUC 0.912–0.928 at the full feature set and 0.900–0.903 at
*K* = 1000 (Fig S1C), so the reported value is not the product of a fortunate partition. The number of
non-zero coefficients did vary across those draws (35–198), i.e. the prediction is stable even where
the selected gene set is not — which is why the panel in Phase E was locked by selection frequency
across folds rather than by any single fit.

## Two glomerular modules carry the MN programme

Signed networks were built separately per compartment (glomerular: 7,882 genes × 153 samples, soft
power 18, scale-free R² = 0.872, 9 modules; tubulointerstitial: 7,616 × 269, power 10, 7 modules).
The glomerular dendrogram shows the nine modules (Fig 4A); the soft-threshold fit and the split-half
membership stability that justified this construction are in Fig S2A and S2B. Four glomerular
modules associated with MN at BH < 0.05 (Fig 4B). Adjusting for podocyte, fibroblast and myeloid
fractions — the guard against a pure cellular-composition artefact — left two significant: a **pink**
module of 84 genes (bicor *r* = 0.55, BH *p* = 6.5 × 10⁻¹², adjusted *p* = 1.0 × 10⁻⁵) and a **red**
module of 144 genes (*r* = 0.46, adjusted *p* = 6.6 × 10⁻⁴). The 2,304-gene turquoise module was the
largest and showed no MN association (*r* = 0.034), a useful internal negative control. The
corresponding eigengenes show pink and red highest in MN — and also elevated in FSGS, the hard
comparator — whereas the black module tracks FSGS and other nephropathies rather than MN and the
green module is flat across diagnoses (Fig 4C), so the MN-associated modules are not simply
mirroring a global expression axis.

**All nine glomerular modules were preserved in ERCB** (Z = 21.4–84.6, *q* < 10⁻¹⁰⁰; Fig 4D) and all
seven tubulointerstitial modules likewise; the unassigned (grey) and random (gold) modules failed in
every test (Z = −27 to +1.2), confirming the test discriminates rather than rubber-stamping.
Functional annotation then gave the two MN-associated modules identities rather than colours (Fig
4E): the pink module is a **cell-surface / calcium-ion-binding /
external-encapsulating-structure** programme — a podocyte surface and matrix-interface signature —
and the red module a **metal-ion transport / lysosomal / glycosylation** programme. A third module
(black, 113 genes) showed ECM-organisation and epithelial–mesenchymal-transition terms but lost
significance after composition adjustment, marking it as an injury-response rather than an
MN-intrinsic programme.

At pathway level the glomerular contrast was dominated by epithelial–mesenchymal transition (NES
+2.38, adjusted *p* = 4.1 × 10⁻¹⁰), apical junction (+2.20) and myogenesis (+2.23; both adjusted
*p* = 3.4 × 10⁻⁷), with fatty-acid metabolism (−2.08, adjusted *p* = 4.3 × 10⁻⁶) and oxidative
phosphorylation (−1.92) suppressed. The underlying enrichments are shown in Fig S3: GO
biological-process ORA of the meta-DEGs (S3A), C7 immunologic-signature GSEA (S3B), the
apical-junction running score (S3C), the epithelial–mesenchymal-transition running score (S3D) and
the MN-versus-non-MN ORA (S3E). These are membership-based enrichments of the ranked gene list, not
measurements of pathway activity.

Split-half reproducibility showed that co-expression **axes** are highly reproducible (eigengene
bicor 0.92–0.96) while gene **membership** is only partly so (Jaccard ≈ 0.3; Fig S2B) — which is why
hub genes were admitted only after independent preservation, their network and module-membership
distributions being shown in Fig S2C and Fig S2D, and why the panel was not selected on membership
alone; per-module annotation of the preserved modules is in Fig S2E. No tubulointerstitial module
tracked MN.

## A four-gene panel retains external discrimination

Selection frequency across nested CV folds for the three selectors is shown in Fig 5A. From 287
candidate genes (of 353 core genes; 68 dropped because they were not measurable in ≥ 95 % of ERCB
samples), nested CV gave mean AUC 0.846 (LASSO), 0.843 (SVM-RFE) and 0.814 (Boruta). Genes
selected by at least two of the three methods at a per-method frequency of 0.60 or more formed the
consensus panel: **SLAMF8, TRPC6, MYOZ2 and PPP1R3C**. *ICAM1*, which appeared in an earlier version
of this panel, was dropped once the cross-validation folds were drawn correctly (one partition per
repeat instead of a fresh random partition per fold); it no longer met the consensus frequency and
its fitted coefficient had the opposite sign to its marginal effect.

| Cohort | n (MN/other) | AUC (95 % CI) |
|---|---|---|
| Discovery (apparent) | 147 (51/96) | 0.922 (0.878–0.967) |
| **ERCB external** | 170 (21/149) | **0.816 (0.721–0.910)** |
| **RNA-seq independent** | 273 (62/211) | **0.901 (0.862–0.939)** |
| Laser-capture GSE216841 | 26 (12/14) | 0.929 (0.824–1.000) |
| ERCB reprocessed with another CDF | 187 (21/166) | 0.848 (0.759–0.938) |
| ERCB incl. healthy donors (secondary) | 196 (21/175) | 0.837 (0.749–0.925) |
| Reference: 13,379-gene elastic net | 170 (21/149) | 0.922 (0.862–0.982) |

Its external ROC curves are in Fig 5B. All four genes contribute positively in the locked model (Fig
5C): odds ratios 2.83 (95 % CI 1.28–6.24) for *SLAMF8*, 3.40 (1.06–10.93) for *TRPC6*, 3.32
(1.20–9.22) for *MYOZ2* and 2.59 (1.31–5.12) for *PPP1R3C* per SD of expression. The four-gene model
uses 1/3,345 of the whole-transcriptome classifier's genes and is below it in ERCB (0.816 versus
0.922) while being close in RNA-seq (0.901 versus 0.939): the compact panel buys interpretability and
transportability, not superiority. The predicted probabilities make the calibration problem visible —
scores saturate near 1 in ERCB and near 0 in RNA-seq while preserving ranking (Fig 5D). ERCB-label
permutations gave *p* = 0.005; a stricter null that refits the panel on shuffled discovery labels
reached *p* = 0.11, which is **not** significant, and we report it as such — with 21 ERCB cases a
randomly selected four-gene model transfers well often enough that the strict null cannot exclude it.

Against individual comparator diagnoses the panel was strongest for the immune-complex and
hypertensive groups (TMD 0.97, HT 0.96, IgAN 0.94, SLE 0.87, MCD 0.88) and weakest for **FSGS
(0.63; 0.44–0.83)** and mixed FSGS/MCD (0.48; Fig 5E), mirroring the gate and the clinical
difficulty [42].

## The atlas assigns the modules to cells, and the deconvolution travels with them

A 236,107-cell kidney atlas (16 cell types) [33,40,41] integrated 49 samples across MN, IgAN, MCD, DKD and
control (Fig 6A); canonical markers supporting every annotation are in Fig S1D and per-sample quality
after filtering in Fig S1E. MN-versus-control
cell composition did not differ after BH correction (smallest adjusted *p* = 0.218; Fig 6B), so the
transcriptomic difference is not a simple abundance shift. NNLS deconvolution of all seven bulk sets
(1,328 samples) assigned glomerular MN samples a podocyte-dominant profile (Fig 6C), and 67 of 316
testable cell-type comparisons survived BH correction (Fig 6D). The podocyte contrast was significant
against the other glomerulopathies — discovery MN 0.637 versus other glomerulonephritis 0.583
(adjusted *p* = 3.3 × 10⁻⁴) and versus FSGS 0.587 (*p* = 0.029), and in RNA-seq MN 0.847 versus MCD
0.745, FSGS 0.753 and other disease 0.730 (*p* = 1.0 × 10⁻⁷ to 2.9 × 10⁻⁵) — but **not** against
healthy glomeruli in discovery (0.637 versus 0.615, adjusted *p* = 0.29), so the finding is
"podocyte-rich relative to other nephropathies", not "podocyte loss". The deconvolution is fragile in
detail — a specificity-filtered and an all-gene signature disagreed on the dominant cell type in
1,314 of 1,328 samples — and bulk-to-pseudobulk concordance was weak (best Pearson 0.29; Fig 6E),
a limit of comparing whole tissue with one cell type's pseudobulk, so only the podocyte result is
interpreted; pDC and neutrophil columns rest on 97 and 535 atlas cells and were excluded from all
downstream trait and attribution analyses.

## Specificity, portability and recalibration

Three further tests constrained interpretation. First, the panel was **compartment-specific**: AUC
fell to 0.449 (0.310–0.589) in ERCB tubulointerstitium and 0.539 (0.447–0.630) in the RNA-seq
tubulointerstitium — chance in both, despite the same patients and platform (Fig S4A). Reprocessing
ERCB through an alternative CDF left discrimination at 0.848 (0.759–0.938), so the result is not a
single-annotation artefact (Fig S4A). Second, the coefficients were platform-portable: refitting the
four genes on RNA-seq and testing in ERCB gave AUC 0.850 with all four signs concordant. Third,
discrimination was portable but **absolute risk was not**: calibration-in-the-large was +0.86 in ERCB
and −0.21 in RNA-seq (Fig S4B). A two-parameter recalibration (intercept and slope) fitted on one
half of ERCB and applied to the other removed the offset (+0.862 → −0.003) and reduced the Brier
score from 0.845 to 0.096 without changing AUC (Fig S4C), and decision-curve analysis showed net
benefit across clinically relevant thresholds (Fig S4D). No single gene was indispensable: every
leave-one-gene-out model exceeded 0.70 in ERCB (Fig S4E), and *MYOZ2* alone matched the full panel
(0.903). Because *SLAMF8* was detected in only 9.5 % of RNA-seq samples, a prespecified three-gene
sensitivity analysis is reported below.

## Where the panel genes live — and one caveat

Atlas localisation placed *TRPC6* in podocytes (44.5 CPM) [43], *PPP1R3C* in mesangium (10.8) and
*SLAMF8* in B cells (1.0), whereas *MYOZ2* was the exception: 93.8 % of atlas cells were zero and its
top column was the unreliable pDC (Fig S5A). MN-versus-control pseudobulk for the panel genes showed
no comparison surviving BH correction with 15 versus 10 donors (smallest raw *p* = 0.00095, *TRPC6*
in mast cells; adjusted *p* = 0.29; Fig S5B), so cell-type-level expression of the panel remains
hypothesis-generating. Because *MYOZ2* is the strongest single gene in ERCB (AUC 0.903 alone), its
signal was tested independently in the submitters' own processed matrix in the sensitivity analyses
below. Its cellular source in kidney is unresolved, and it is the first candidate for protein-level
confirmation.

## Sensitivity: comparators, preprocessing and a three-gene alternative

Because *SLAMF8* is detected in only 9.5 % of RNA-seq samples, a prespecified model without it was
tested: as a three-gene assay it is comparable to the locked panel (ERCB 0.845, 95 % CI 0.753–0.937;
RNA-seq 0.891; laser-capture 0.941; permutation *p* = 0.005; Fig S5C), so a three-gene format is
defensible when a platform cannot measure *SLAMF8* — but it is no longer a strict improvement, so the
locked four-gene panel remains the primary result.

To exclude a preprocessing artefact we re-derived each gene's single-marker separation in the
submitters' own processed matrix, which shares no step with our pipeline; the two estimates agreed to
within 0.01 for all four genes (*MYOZ2* 0.865 versus 0.866; *TRPC6* 0.857 versus 0.863; *SLAMF8*
0.800 versus 0.808; *PPP1R3C* 0.796 versus 0.786; Fig S5D). The panel score also carried little
global expression structure (R² = 0.37 in discovery and 0.21 in ERCB against PC1–5), and adjusting
for those components did not reduce discrimination.

Finally, the panel was re-scored under the cross-platform harmonisation used for the gate. The
answer does not depend on that convention: ERCB discrimination is 0.817 versus 0.816 unadjusted and
RNA-seq 0.899 versus 0.901, while calibration-in-the-large improves from +0.86 to +0.20 in ERCB and
from −0.21 to +0.11 in RNA-seq, and the score range stops saturating at 0.75–1.00. The panel is
therefore reported on its original, simpler path, with harmonisation as a documented robustness check.

## PLA2R stratum

Within the laser-capture cohort the panel detected MN against MCD in PLA2R-positive (0.976;
0.921–1.000, *n* = 6) and PLA2R-negative disease (0.881; 0.712–1.000, *n* = 6), and against all MCD
at 0.929 (0.824–1.000). The within-MN contrast (PLA2R-positive versus PLA2R-negative) was 0.833
(0.554–1.000) on six cases per group: it is numerically higher than the pre-audit panel's 0.528 but
the confidence interval spans chance and the sample is far too small to interpret. We therefore make
no claim about PLA2R tracking in either direction (Fig S5E).

## Exploratory subtyping does not reproduce

Consensus clustering [44] of the glomerular discovery cohort gave a stable *k* = 2 partition (PAC
0.0015; Fig S6A), but the RNA-seq cohort had no stable partition at all (PAC 0.215 at *k* = 5; Fig
S6B), cross-platform subtype-defining genes shared only 6 of 889 (Jaccard 0.007; Fig S6C), maximal
centroid correlation was 0.50, and the array partition tracked technical covariates (median
sample–sample correlation 0.686 versus 0.638, *p* = 4 × 10⁻⁴; Fig S6D). The RNA-seq consensus matrix
at *k* = 3 shows no block structure (Fig S6E). We therefore report molecular subtyping of MN as not
reproducible in these data and do not present subtype–outcome claims.

# Discussion

Three findings organise the paper.

First, the glomerular transcriptome of MN is separable from other proteinuric glomerulopathies to a
degree that thresholded differential expression hides. Only 18 genes passed a conventional threshold
versus FSGS while a multivariate classifier reached external AUC 0.86–0.92 (0.92 at the full feature
set). Discrimination lives
in many genes of moderate effect, which is exactly what a machine-learning panel exploits and what a
DEG list discards. This has a practical consequence for the field: negative "MN versus FSGS DEG"
results should not be read as evidence that the diseases are transcriptomically indistinguishable.

Second, the signal is compartment-resolved and modular. It sits in two preserved, composition-adjusted
glomerular modules — a podocyte surface/calcium programme and a lysosomal/glycosylation programme —
and collapses to chance in the tubulointerstitium of the same patients. That combination (module
preservation in an independent cohort plus compartment specificity plus insensitivity to cell
composition adjustment) is stronger than any of the three alone, and it explains why
tubulointerstitial MN studies have produced inconsistent marker lists.

Third, a four-gene classifier is achievable but its honest description is narrow. It discriminates
(AUC 0.82 external, 0.90 cross-platform, 0.93 in laser-capture tissue, 0.85 across array
annotation) and it does not calibrate across cohorts without refitting; recalibration works but must
be done per platform. Two of its four genes need protein-level confirmation — *MYOZ2* because its
kidney cell of origin is unresolved despite reproducible tissue-level signal, and *SLAMF8* because it
is barely detected on two platforms and adds little to discrimination, so a three-gene version
without it is a defensible alternative rather than an improvement.

Our negative results are part of the contribution. Exploratory subtyping did not reproduce beyond the
discovery platform and tracked technical covariates; we report it rather than publishing a
platform-specific partition. Tubulointerstitial differential expression replicated poorly in ERCB
(32 %) and was null against FSGS, so no interstitial signature is claimed.

A fourth point concerns how easily transcriptomic discrimination is overstated, and here it cuts both
ways. Standardising each cohort against its own mean and standard deviation is transductive: it lets
the validation cohort's distribution set the model's feature weights. Removing it is necessary but
not sufficient, because the discovery and validation cohorts were assayed on different array
platforms. Our own first correction removed the transductive step and left training statistics in
place; the gate's ERCB AUC promptly fell to 0.64–0.77 and four comparator diagnoses came out ranking
above MN — a transfer that was measuring platform shift, and which a less careful analysis would have
published either as a weak result or, with the transductive step retained, as a strong one. What
makes the transfer valid is the pair of steps: harmonise the validation cohort to the discovery
reference distribution without using outcome information, then apply training-set coefficients. Any
MN marker paper claiming external validation should state both. The same check applied to the
four-gene panel is reassuring: harmonisation leaves its ERCB discrimination unchanged (0.82) and
improves its calibration, so the panel result does not depend on the standardisation choice.

## Limitations

Clinical outcomes were unavailable: baseline eGFR, proteinuria, remission and relapse status are not
in any downloaded file, so no prognostic model, subtype–outcome analysis or treatment-response claim
is made. PLA2R status was recoverable only for the 12 MN samples of the laser-capture cohort (6
positive, 6 negative), which is far too few to estimate stratum-specific performance, and the
within-MN contrast is uninterpretable. The RNA-seq cohorts carry no participant
identifiers, precluding
patient-level replication. The single-cell reference contains six MN donors, so pseudobulk
differential expression is hypothesis-grade, and single-nucleus data under-represent cytoplasmic
transcripts. Two atlas cell types (pDC, neutrophil) rest on fewer than 600 cells and were excluded
from trait and attribution analyses. The harmonisation that makes the external gate valid is a
batch-level operation — it needs the whole validation cohort to estimate each gene's marginal
distribution — so deploying the classifier on a single new sample would require either a reference
batch or a per-platform transform fitted in advance. Finally, ERCB validation rests on 21 MN samples, so confidence
intervals are wide and the stricter refit-based permutation null does not reach significance
(*p* = 0.11).

## Conclusion

MN carries a reproducible, ERCB-preserved, compartment-specific glomerular transcriptomic signature
concentrated in podocyte surface/calcium and lysosomal/glycosylation programmes, from which a
four-gene classifier retains external discrimination below, but in the same range as, a
whole-transcriptome model. The classifier is a discrimination instrument that requires per-platform
recalibration, and its translation now depends on protein-level validation and on clinical annotation
that the current public data do not contain.

# Contributors

Lingzhang Meng conceived the study, assembled and harmonised the public cohorts, wrote and ran the
analysis code, performed the statistical analyses, generated the figures and drafted the manuscript.
[Further authors to be added: contributions to be recorded here in CRediT format before submission.]

# Data sharing statement

All primary data are public GEO accessions (listed in the Methods); no new data were generated.
Analysis code, harmonised metadata, all result tables and figure-generating scripts are retained in
the project repository, together with six verification scripts that reproduce every reported number.
The harmonised sample metadata, the seven cohort definitions, the derived network objects and the
result tables are available from the corresponding author on reasonable request and are intended for
deposition in a public repository at acceptance.

# Declaration of interests

The author declares no competing interests.

# Acknowledgements

The author thanks the investigators of the Nephrotic Syndrome Study Network (NEPTUNE), the European
Renal cDNA Bank and the submitters of every GEO series reanalysed here for making their data public.

# References

*Vancouver style, numbered in order of first citation. Every entry is held with its identifier in `Manuscript/references.tsv`, from which this list was generated; verify against the DOI at proof stage.*

1. Couser WG. Primary Membranous Nephropathy. Clin J Am Soc Nephrol. 2017;12(6):983–997. doi:10.2215/cjn.11761116
2. Ronco P, Beck L, Debiec H, Fervenza FC, Hou FF, Jha V, et al. Membranous nephropathy. Nat Rev Dis Primers. 2021;7(1):69. doi:10.1038/s41572-021-00303-z
3. Beck LH, Bonegio RG, Lambeau G, Beck DM, Powell DW, Cummins TD, et al. M-type phospholipase A2 receptor as target antigen in idiopathic membranous nephropathy. N Engl J Med. 2009;361(1):11–21. doi:10.1056/nejmoa0810457
4. Debiec H, Ronco P. PLA2R autoantibodies and PLA2R glomerular deposits in membranous nephropathy. N Engl J Med. 2011;364(7):689–690. doi:10.1056/nejmc1011678
5. Tomas NM, Beck LH, Meyer-Schwesinger C, Seitz-Polski B, Ma H, Zahner G, et al. Thrombospondin type-1 domain-containing 7A in idiopathic membranous nephropathy. N Engl J Med. 2014;371(24):2277–2287. doi:10.1056/nejmoa1409354
6. Sethi S, Debiec H, Madden B, Charlesworth MC, Morelle J, Gross L, et al. Neural epidermal growth factor-like 1 protein (NELL-1) associated membranous nephropathy. Kidney Int. 2020;97(1):163–174. doi:10.1016/j.kint.2019.09.014
7. Sethi S, Madden BJ, Debiec H, Charlesworth MC, Gross L, Ravindran A, et al. Exostosin 1/Exostosin 2-Associated Membranous Nephropathy. J Am Soc Nephrol. 2019;30(6):1123–1136. doi:10.1681/asn.2018080852
8. Al-Rabadi LF, Caza T, Trivin-Avillach C, Rodan AR, Andeen N, Hayashi N, et al. Serine Protease HTRA1 as a Novel Target Antigen in Primary Membranous Nephropathy. J Am Soc Nephrol. 2021;32(7):1666–1681. doi:10.1681/asn.2020101395
9. Caza TN, Hassen SI, Kuperman M, Sharma SG, Dvanajscak Z, Arthur J, et al. Neural cell adhesion molecule 1 is a novel autoantigen in membranous lupus nephritis. Kidney Int. 2021;100(1):171–181. doi:10.1016/j.kint.2020.09.016
10. Kidney Disease: Improving Global Outcomes (KDIGO) Glomerular Diseases Work Group. KDIGO 2021 Clinical Practice Guideline for the Management of Glomerular Diseases. Kidney Int. 2021;100(4S):S1–S276. doi:10.1016/j.kint.2021.05.021
11. Sethi S, Haas M, Markowitz GS, D'Agati VD, Rennke HG, Jennette JC, et al. Mayo Clinic/Renal Pathology Society Consensus Report on Pathologic Classification, Diagnosis, and Reporting of GN. J Am Soc Nephrol. 2016;27(5):1278–1287. doi:10.1681/asn.2015060612
12. Fervenza FC, Appel GB, Barbour SJ, Rovin BH, Lafayette RA, Aslam N, et al. Rituximab or Cyclosporine in the Treatment of Membranous Nephropathy. N Engl J Med. 2019;381(1):36–46. doi:10.1056/nejmoa1814427
13. Beck LH, Fervenza FC, Beck DM, Bonegio RG, Malik FA, Erickson SB, et al. Rituximab-induced depletion of anti-PLA2R autoantibodies predicts response in membranous nephropathy. J Am Soc Nephrol. 2011;22(8):1543–1550. doi:10.1681/asn.2010111125
14. Ju W, Greene CS, Eichinger F, Nair V, Hodgin JB, Bitzer M, et al. Defining cell-type specificity at the transcriptional level in human disease. Genome Res. 2013;23(11):1862–1873. doi:10.1101/gr.155697.113
15. Hoxha E, Kneißler U, Stege G, Zahner G, Thiele I, Panzer U, et al. Enhanced expression of the M-type phospholipase A2 receptor in glomeruli correlates with serum receptor antibodies in primary membranous nephropathy. Kidney Int. 2012;82(7):797–804. doi:10.1038/ki.2012.209
16. Ronco P, Debiec H. Pathophysiological advances in membranous nephropathy: time for a shift in patient's care. Lancet. 2015;385(9981):1983–1992. doi:10.1016/s0140-6736(15)60731-0
17. Gadegbeku CA, Gipson DS, Holzman LB, Ojo AO, Song PX, Barisoni L, et al. Design of the Nephrotic Syndrome Study Network (NEPTUNE) to evaluate primary glomerular nephropathy by a multidisciplinary approach. Kidney Int. 2013;83(4):749–756. doi:10.1038/ki.2012.428
18. Cohen CD, Frach K, Schlöndorff D, Kretzler M. Quantitative gene expression analysis in renal biopsies: a novel protocol for a high-throughput multicenter application. Kidney Int. 2002;61(1):133–140. doi:10.1046/j.1523-1755.2002.00113.x
19. Carvalho BS, Irizarry RA. A framework for oligonucleotide microarray preprocessing. Bioinformatics. 2010;26(19):2363–2367. doi:10.1093/bioinformatics/btq431
20. Gautier L, Cope L, Bolstad BM, Irizarry RA. affy--analysis of Affymetrix GeneChip data at the probe level. Bioinformatics. 2004;20(3):307–315. doi:10.1093/bioinformatics/btg405
21. Irizarry RA, Hobbs B, Collin F, Beazer-Barclay YD, Antonellis KJ, Scherf U, et al. Exploration, normalization, and summaries of high density oligonucleotide array probe level data. Biostatistics. 2003;4(2):249–264. doi:10.1093/biostatistics/4.2.249
22. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26(1):139–140. doi:10.1093/bioinformatics/btp616
23. Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 2014;15(2):R29. doi:10.1186/gb-2014-15-2-r29
24. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007;8(1):118–127. doi:10.1093/biostatistics/kxj037
25. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47. doi:10.1093/nar/gkv007
26. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
27. Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1(6):417–425. doi:10.1016/j.cels.2015.12.004
28. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102(43):15545–15550. doi:10.1073/pnas.0506580102
29. Friedman J, Hastie T, Tibshirani R. Regularization Paths for Generalized Linear Models via Coordinate Descent. J Stat Softw. 2010;33(1):1–22. doi:10.18637/jss.v033.i01
30. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559. doi:10.1186/1471-2105-9-559
31. Hao Y, Stuart T, Kowalski MH, Choudhary S, Hoffman P, Hartman A, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol. 2024;42(2):293–304. doi:10.1038/s41587-023-01767-y
32. Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Zhang F, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods. 2019;16(12):1289–1296. doi:10.1038/s41592-019-0619-0
33. Lake BB, Menon R, Winfree S, Hu Q, Melo Ferreira R, Kalhor K, et al. An atlas of healthy and injured cell states and niches in the human kidney. Nature. 2023;619(7970):585–594. doi:10.1038/s41586-023-05769-3
34. Newman AM, Liu CL, Green MR, Gentles AJ, Feng W, Xu Y, et al. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods. 2015;12(5):453–457. doi:10.1038/nmeth.3337
35. Newman AM, Steen CB, Liu CL, Gentles AJ, Chaudhuri AA, Scherer F, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol. 2019;37(7):773–782. doi:10.1038/s41587-019-0114-2
36. Kursa MB, Rudnicki WR. Feature selection with the Boruta package. J Stat Softw. 2010;36(11):1–13. doi:10.18637/jss.v036.i11
37. Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez JC, et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics. 2011;12:77. doi:10.1186/1471-2105-12-77
38. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics. 1988;44(3):837–845. doi:10.2307/2531595
39. Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating prediction models. Med Decis Making. 2006;26(6):565–574. doi:10.1177/0272989x06295361
40. Wilson PC, Wu H, Kirita Y, Uchimura K, Ledru N, Rennke HG, et al. The single-cell transcriptomic landscape of early human diabetic nephropathy. Proc Natl Acad Sci U S A. 2019;116(39):19619–19625. doi:10.1073/pnas.1908706116
41. Menon R, Otto EA, Hoover P, Eddy S, Mariani L, Godfrey B, et al. Single cell transcriptomics identifies focal segmental glomerulosclerosis remission endothelial biomarker. JCI Insight. 2020;5(6):133267. doi:10.1172/jci.insight.133267
42. D'Agati VD, Kaskel FJ, Falk RJ. Focal segmental glomerulosclerosis. N Engl J Med. 2011;365(25):2398–2411. doi:10.1056/nejmra1106556
43. Winn MP, Conlon PJ, Lynn KL, Farrington MK, Creazzo T, Hawkins AF, et al. A mutation in the TRPC6 cation channel causes familial focal segmental glomerulosclerosis. Science. 2005;308(5729):1801–1804. doi:10.1126/science.1106215
44. Wilkerson MD, Hayes DN. ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. Bioinformatics. 2010;26(12):1572–1573. doi:10.1093/bioinformatics/btq170

# Figure legends

**Figure 1. Cohort harmonisation and expression quality control.** (A) PCA of the glomerular
discovery cohort (top 2,000 variable genes), coloured by diagnosis. (B) PCA of the independent ERCB
glomerular validation cohort. (C) PCA of the independent RNA-seq glomerular cohort. (D) PCA of the
two-series tubulointerstitial discovery set before and after ComBat correction. (E) PCA of the
laser-capture cohort used for the PLA2R-stratified analysis.

**Figure 2. Differential expression and random-effects meta-analysis.** (A) Volcano plot, glomerular
MN versus control in discovery. (B) Heatmap of the 584 glomerular meta-DEGs across the discovery and
RNA-seq cohorts. (C, D) Representative random-effects forest plots for *CCNL1* and *HBA2*. (E)
Volcano plot, meta-analysis of MN versus FSGS. (F) Volcano plot, meta-analysis of MN versus all
non-MN disease.

**Figure 3. Separability of MN from its mimics, and why a DEG count hides it.** (A) ROC of the
discovery-trained elastic-net gate in ERCB across feature counts *K*. (B) External ERCB AUC as a
function of *K* (training-statistic standardisation with platform harmonisation), with the
independent RNA-seq cohort shown for comparison. (C) Permutation null (20 discovery-label shuffles)
against the observed ERCB AUC. (D) AUC against each individual comparator diagnosis in ERCB. (E)
Multiclass confusion matrix in ERCB at *K* = all (rows, predicted class; columns, true diagnosis).

**Figure 4. Glomerular co-expression modules, trait association and ERCB preservation.** (A) Gene
dendrogram with module colours. (B) Module–trait relationships (bicor, BH-adjusted). (C) Module
eigengenes by diagnosis for the black, green, pink and red modules. (D) Module preservation in ERCB
(Zsummary versus module size). (E) Functional identity of the two MN-associated modules.

**Figure 5. A four-gene glomerular panel.** (A) Selection frequency across nested CV folds for the
three selectors. (B) ROC curves of the locked panel in the external cohorts. (C) Locked-model
coefficients with 95 % confidence intervals. (D) Predicted probability of MN by diagnosis in
discovery, ERCB and RNA-seq. (E) AUC against each comparator diagnosis.

**Figure 6. Kidney single-cell atlas and bulk deconvolution.** (A) Atlas cell types (UMAP). (B)
Cell-type composition by disease. (C) NNLS deconvolution of all seven bulk cohorts (stacked
proportions). (D) Deconvolved cell-type proportions by disease group. (E) Bulk-to-pseudobulk
concordance.

**Figure S1. Cohort QC, gate stability and single-cell annotation.** (A) PCA of the ERCB
tubulointerstitial cohort. (B) PCA of the RNA-seq tubulointerstitial cohort. (C) ERCB AUC of the
locked gate refitted under ten independent random-number streams at *K* = 1000 and at the full
feature set (horizontal bar = median). (D) Canonical marker expression by atlas cell type. (E)
Per-sample single-cell quality control after filtering.

**Figure S2. Module reproducibility, hub genes and annotation detail.** (A) Soft-threshold
diagnostics. (B) Split-half stability of module eigengenes and gene membership. (C) Hub-gene network
of the preserved modules. (D) Module membership (|kME|) of hub genes by module, coloured by
preservation status. (E) Per-module enrichment for the preserved modules.

**Figure S3. Functional profiling of the glomerular signature.** (A) GO biological-process ORA of
the glomerular meta-DEGs. (B) C7 immunologic-signature GSEA. (C) Running enrichment score for
apical-junction genes. (D) Running enrichment score for epithelial–mesenchymal-transition genes. (E)
ORA for the MN-versus-non-MN contrast.

**Figure S4. Panel robustness, calibration and transfer.** (A) ROC of the panel in the
tubulointerstitial compartment and after alternative-CDF reprocessing. (B) Calibration of the locked
panel before recalibration. (C) Recalibration fitted on one half of ERCB and scored on the other.
(D) Decision-curve analysis. (E) Leave-one-gene-out sensitivity.

**Figure S5. Where the panel genes sit, panel sensitivity and PLA2R stratum.** (A) Panel-gene
expression across the 16 atlas cell types (max CPM per cell type). (B) Panel genes in
MN-versus-control pseudobulk per cell type. (C) Locked four-gene versus three-gene (without *SLAMF8*)
ROC in ERCB. (D) Single-gene AUC computed in our pipeline versus the submitters' processed matrix.
(E) Panel performance within the laser-capture cohort: MN versus MCD overall, by PLA2R stratum, and
within MN.

**Figure S6. Exploratory subtyping is not reproducible.** (A) Array consensus matrix at *k* = 2.
(B) PAC stability curve. (C) Centroid concordance between platforms. (D) Array subtype heatmap with
its technical confounder. (E) RNA-seq consensus matrix at *k* = 3.
