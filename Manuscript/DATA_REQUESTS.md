---
title: "How to obtain the two remaining external inputs"
subtitle: "NEPTUNE clinical variables (B2) and the participant map (B7)"
author: "Lingzhang Meng"
date: "2026-10-19"
---

# Status

| Blocker | Status | Effort |
|---|---|---|
| **B12 — PLA2R status** | **RESOLVED locally (2026-10-19)** — the `phenotype: PLA2R_Positive/Negative` characteristic is present in `GSE216841_series_matrix.txt.gz`; `R/09f_pla2r_strata.R` parses it into `metadata/GSE216841_pla2r.tsv` (6 positive / 6 negative). No external request needed | none |
| **B2 — clinical variables** (eGFR, UPCR, remission) | Blocked: not present in any public file | one application, 2–8 weeks |
| **B7 — participant map for the RNA-seq cohorts** | Blocked: NEPTUNE-linked, held by the DACC | same application as B2 |

# B2 — NEPTUNE clinical variables

The NEPTUNE resource is distributed through the **NIDDK Central Repository** under a "Resources for
Research (R4R)" request, not from GEO.

1. Register at **https://repository.niddk.nih.gov/** and locate the study
   *Nephrotic Syndrome Study Network (NEPTUNE)* (Version 1; DOI `10.58020/f2ae-n094`).
2. Submit an R4R request selecting the **phenotype/clinical** data package (the repository holds the
   subject-level phenotype tables; the WGS is separate, under dbGaP `phs003210.v1.p1`, and is not
   needed for this project).
3. Expect a data-use agreement. Requests typically ask for a short analysis plan — the project
   description below can be pasted in.
4. If you need the **glomerular expression data linked to phenotype**, that is the dbGaP/GEO path
   (`GSE197307`), and the link to clinical variables must come through the DACC (next section).
5. For an **ancillary study** (e.g. adding a new assay rather than reusing data), submit the
   NEPTUNE Ancillary Studies application; the DACC handles review and data transfer. Contact:
   **NEPTUNE-STUDY@umich.edu**; the ancillary application form is distributed by the DACC
   (historically hosted via the University of Michigan / nephcure-NEPTUNE pages).

**Project description to paste into the request**

> We are completing a multi-cohort transcriptomic study of membranous nephropathy (MN) that uses the
> public NEPTUNE glomerular and tubulointerstitial expression data (GEO GSE200828, GSE133288,
> GSE200818, GSE197307, GSE182380). We have defined a five-gene glomerular classifier and two
> ERCB-preserved co-expression modules, and now need the clinical annotation to test whether module
> activity and the classifier score associate with baseline kidney function, proteinuria and
> remission status.

**Variables requested (subject level, one row per participant)**

| Variable | Format | Purpose |
|---|---|---|
| `subject_id` / `sample_id` | string | join key to the expression sample |
| `visit` (baseline, month 4, 12, 24…) | string | longitudinal structure |
| `age_at_biopsy`, `sex`, `race`, `ethnicity` | numeric/string | covariates |
| **eGFR** (with the equation used, e.g. CKD-EPI 2009/2021, bedside Schwartz for children) | numeric | primary outcome axis |
| **UPCR** or 24 h proteinuria | numeric + unit | activity axis |
| serum albumin, serum creatinine | numeric | context |
| **remission / relapse status and dates** (complete, partial, none) | string + date | outcome endpoint |
| treatment at each visit (RAAS blockade, rituximab, cyclophosphamide, calcineurin inhibitor, steroids) | string | confounder |
| PLA2R antibody status and date (if held) | string + date | antigen strata (cf. our `metadata/GSE216841_pla2r.tsv`) |
| biopsy diagnosis and cohort (MN / FSGS / MCD / other) | string | grouping |

Ask for **de-identified** data with dates as study-relative offsets, which is the repository's
standard format.

# B7 — participant map for the RNA-seq cohorts

The two RNA-seq series we use (`GSE197307`, `GSE182380`) carry **no participant identifier**: their
sample titles are local indices (`Sample_1` … `Sample_199`), so patient-level replication cannot be
tested. Two routes:

1. **Preferred** — ask the NEPTUNE DACC for a **GSM → NEPTUNE sample ID → subject ID** crosswalk for
   the glomerular RNA-seq series. The NIDDK Central Repository phenotype release includes a
   subject→sample mapping table, so this usually comes with the B2 request; state explicitly that
   you need the *mapping table* alongside the clinical variables.
2. **Fallback** — contact the GEO submitters of `GSE182380` directly (the corresponding author is
   listed on the GEO record) and ask for the sample-to-patient key for their series.

A related trap to record: a further NEPTUNE release exists — **`GSE219185`** (tubulointerstitial
RNA-seq, voom-transformed). It is from the *same* NEPTUNE participants as `GSE197307`, so it must
**never** be used as independent validation. It is only useful once the participant map exists, and
then only for compartment contrasts within a patient.

# What to do when the data arrive

1. Drop the clinical export into `Datasets/NEPTUNE_clinical/` (or any path; it stays outside the
   repository because it is access-controlled).
2. Convert to a one-row-per-sample table with the columns above and save as
   `metadata/neptune_clinical.csv`.
3. The analyses that are already wired for it:
   - module–trait: append eGFR and UPCR as columns to `wgcna_traits.csv` and re-run
     `R/08b_module_trait.R` (no recomputation of the networks needed);
   - deconvolution–function: correlate `deconvolution_fractions.csv` with eGFR/UPCR;
   - outcome models: Cox/logistic models of remission by panel score and module eigengene, with
     age, sex, eGFR and treatment as covariates;
   - Phase F (subtypes) only becomes interpretable with these variables, which is why it is
     deferred rather than abandoned.
4. Re-run the verifiers: the trait-count assertions in `99_verify_phase_d.R` will need the new
   expected column count, so update that assertion in the same commit.

# Practical notes

- The repository route is free but requires a signed data-use agreement and institutional
  signature; start it early, because it is the long pole for the paper's prognostic claims.
- Do **not** put access-controlled clinical data in a public repository when the manuscript is
  submitted; report it as "available from the NIDDK Central Repository under a data-use agreement"
  in the Data Availability statement.
- If the application stalls, the paper is still publishable as it stands: every claim in the current
  draft is supported by data already in hand (see the Limitations section of the manuscript).
