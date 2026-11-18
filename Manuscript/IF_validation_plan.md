---
title: "Protein-level validation plan for the five-gene MN panel"
subtitle: "Immunofluorescence / immunohistochemistry on archival glomerular tissue"
author: "Lingzhang Meng"
date: "2026-10-19"
---

# Purpose and decision rules

The bulk and single-cell analyses establish five candidate transcripts. This plan tests whether the
corresponding **proteins** are present in the expected glomerular compartment and whether their
abundance differs between MN and its mimics. It is designed so that the result can be reported
whatever it shows.

**Pre-specified success criteria**

| Target | Expected localisation (from the atlas) | Success criterion |
|---|---|---|
| TRPC6 | podocyte (44.5 CPM) | podocyte signal in ≥ 70 % of MN biopsies and higher than FSGS/MCD (pre-specified ordering test) |
| ICAM1 | endothelium, myeloid (108.9) | endothelial/glomerular capillary signal, MN > FSGS/MCD |
| PPP1R3C | mesangium (10.8) | mesangial signal, MN ≥ FSGS/MCD |
| **MYOZ2** | **unresolved** (93.8 % of cells zero) | **exploratory**: any reproducible glomerular pattern; absence is an equally reportable outcome |
| SLAMF8 | B cells (1.0) | exploratory; if absent, it supports dropping it and reporting the four-gene model |

A target that fails is not a failed experiment: it converts a transcript-level panel into an
honest statement about which proteins are measurable, and it justifies the four-gene sensitivity
model already reported.

# Sample selection

- **MN**: 10–15 biopsies with residual tissue, ideally 6–8 PLA2R-positive and 4–6 PLA2R-negative
  (the biobank's serology and glomerular PLA2R staining define the strata).
- **Disease controls**: FSGS 6–8, minimal change disease 4–6, IgA nephropathy 4 (optional
  specificity check).
- **Non-disease control**: 4–5 tumour-nephrectomy or donor biopsies with histologically normal
  glomeruli.
- Matching on age, sex and estimated glomerular filtration rate at biopsy where the clinical record
  allows, because podocyte protein abundance falls with chronic damage.
- Block selection: use cases with ≥ 3 glomeruli per section and without global sclerosis in > 50 %
  of glomeruli, so that absence of signal can be interpreted.

# Antibodies and controls

Use only antibodies with published validation in kidney or an explicit provider RRID; record lot
numbers. Verify each candidate first against a positive-control tissue before touching study
sections:

| Target | Positive control tissue | Notes |
|---|---|---|
| TRPC6 | kidney (podocytes), lung | membrane pattern expected; beware podocyte foot-process staining in FSGS |
| ICAM1 | tonsil (endothelium), inflamed kidney | strong endothelial signal; run isotype control on same section |
| PPP1R3C | kidney (mesangium), skeletal muscle | cytoplasmic |
| MYOZ2 | **skeletal muscle** | establishes that the antibody works at all; a negative kidney result is then interpretable |
| SLAMF8 | tonsil (B-cell zones) | low abundance; consider tyramide signal amplification |

Required controls on every run: isotype-matched primary, omission of primary, and an internal
positive structure (vascular endothelium for ICAM1; tubular basement membrane/epithelium as a
general staining control). If the antibody datasheet documents peptide-competition validation, run
one absorption control per target.

# Protocol outline (multiplex immunofluorescence)

1. FFPE sections, 3–4 µm, on charged slides; bake 60 °C 60 min; dewax and rehydrate.
2. Antigen retrieval calibrated per antibody (citrate pH 6, 95 °C, 20–40 min for cytoplasmic
   targets; Tris-EDTA pH 9 for membrane targets). Do not assume one condition fits all five.
3. Block: 5 % normal serum + 1 % BSA + 0.1 % Triton X-100, 30 min. Add a commercial
   Fc-receptor/endogenous-fluorescence block if the tissue shows autofluorescence.
4. Primary antibodies overnight at 4 °C in blocking buffer. Include a lineage anchor in the same
   panel so the compartment is unambiguous:
   - podocyte: nephrin or podocin (or gp135); endothelial: CD31/ERG; mesangial: PDGFRB or desmin;
     B cells: CD20.
5. Secondary detection: species-appropriate Alexa 488 / 555 / 647; for low-abundance targets use
   tyramide amplification with a stripping step between cycles.
6. Counterstain DAPI, mount in anti-fade medium, image on a confocal or high-resolution widefield
   system with identical settings across the cohort (record gain/exposure per channel).
7. **Acquisition discipline:** capture ≥ 5 glomeruli per case at 40× and at least one full-glomerulus
   stack; save raw channels so that segmentation is not baked into the image.

# Quantification and statistics

- Blinded, two-observer scoring; report inter-observer agreement (Cohen's κ) before any group
  comparison.
- Per glomerulus, measure the fraction of area positive inside the lineage-marker mask
  (podocyte/endothelial/mesangial), not whole-glomerulus mean intensity, so that the readout is
  cell-type specific.
- Compare groups by linear mixed-effects model with case as a random effect (multiple glomeruli per
  case) and section batch as a covariate, following the same composition-aware logic used for the
  transcript modules.
- Pre-specified primary comparison: MN versus pooled FSGS/MCD. Secondary: MN versus normal
  glomeruli; PLA2R-positive versus PLA2R-negative.
- Report effect sizes with confidence intervals, not only *p* values, and state the number of
  glomeruli per group (not just cases).

# Practical pitfalls

1. **Podocyte protein loss with chronic damage** can create a spurious "lower in MN" result; that is
   why groups must be balanced on eGFR and global sclerosis.
2. **Autofluorescence** is strong in fibrotic glomeruli — use spectral unmixing or a
   autofluorescence-quenching step (Sudan Black, 0.1 %), and always run a no-primary section per case.
3. **Antibody-dependent signal on capillaries** (ICAM1) can be mistaken for glomerular expression;
   require endothelial colocalisation.
4. **MYOZ2** is a striated-muscle gene: contamination by capsule, arteriolar smooth muscle or muscle
   fibres will produce false glomerular signal. Score only well-separated glomeruli and record
   capsule distance.
5. Do not attempt all five targets on one slide as a first pass — validate singly, then multiplex
   three at a time.

# Orthogonal (non-tissue) options

- **Plasma protein measurement** for the secreted/membrane candidates (sILAMF8, sICAM1) by ELISA in
  the same biobank's serum, which would separate a tissue-expression result from a systemic one.
- **RNAscope in situ hybridisation** for MYOZ2, which is more sensitive than IF at low abundance and
  localises transcripts to a specific cell in the same section.
- Single-cell or spatial profiling of 2–3 MN biopsies, if the biobank allows, which would resolve the
  MYOZ2 source directly.

# Expected deliverable

A single supplementary figure with five rows (targets), four columns (MN, FSGS, MCD, normal) and
per-glomerulus quantification, plus one table of effect sizes with confidence intervals and the
κ for observer agreement. Whatever the outcome, the result should be reported — including a
negative MYOZ2 result, which would justify presenting the four-gene model as the deployable assay.
