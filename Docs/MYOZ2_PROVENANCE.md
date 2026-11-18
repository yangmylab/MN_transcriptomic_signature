# MYOZ2 provenance test (P1b, 2026-11-19)

**Status:** done. Deferred consolidation analysis, **not** a Paper 1 method and **not** a change to any
frozen claim. Script: `tools/moyz2_provenance.R`; artefacts listed at the end. The manuscript wording
"tissue-level reproducibility with unresolved cellular source" remains correct — this test adds a
*tested* statement about what the bulk signal is **not**, and quantifies how little is resolvable at
single-cell level.

**Question.** MYOZ2 is the strongest single panel gene in bulk (ERCB AUC 0.903; leave-one-out without
it 0.738) but is essentially absent from kidney single-cell/single-nucleus profiles. Is the bulk
signal a perivascular / smooth-muscle / mesangial readout (ACTA2, MYH11, TAGLN, PDGFRB), a podocyte
readout (NPHS1, NPHS2), a fibroblast readout (COL1A1), or none of them?

**Pre-specified decision rule** (in the script header, before the run): "tracks the SMC programme" =
≥ 3 of the 4 SMC markers with Spearman rho ≥ 0.5 and BH q < 0.05 in ≥ 3 of the 4 cohorts; "does not
track" = fewer than 2 such markers in ≥ 3 of the 4 cohorts; otherwise "cohort-dependent". The rule
is applied separately to all non-UNKNOWN samples and to MN-only samples. The empirical null is the
distribution of rho(MYOZ2, every gene) in the same cohort, so each marker gets a percentile and a
BH q over all genes; a partial Spearman controlling for each sample's median expression is reported
as a sensitivity analysis.

## 1. Bulk — four cohorts

Cohorts: discovery array (n = 153), ERCB array (n = 196), RNA-seq cross-platform (n = 281, one
UNKNOWN dropped), LCM (n = 34).

Spearman rho, MYOZ2 vs marker (all samples; MN-only in parentheses):

| Marker | Discovery array | ERCB array | RNA-seq | LCM |
|---|---|---|---|---|
| ACTA2 | 0.14 (0.01) | 0.15 (−0.07) | 0.32 (0.06) | −0.16 (0.33) |
| MYH11 | 0.08 (0.05) | 0.17 (−0.22) | 0.23 (−0.04) | 0.17 (−0.02) |
| TAGLN | 0.24 (0.44) | 0.23 (0.11) | **0.52** (0.28) | −0.05 (0.20) |
| PDGFRB | 0.27 (0.22) | 0.21 (0.15) | 0.23 (−0.01) | 0.10 (−0.64) |
| SMC composite (mean z) | 0.27 (0.32) | 0.27 (0.01) | 0.45 (0.11) | 0.03 (0.09) |
| Podocyte composite | 0.35 (0.43) | 0.18 (0.37) | 0.37 (0.31) | 0.30 (−0.12) |
| NPHS2 | 0.35 (0.31) | 0.28 (0.52) | 0.53 (0.37) | 0.81 (0.40) |
| COL1A1 | 0.03 (0.10) | −0.21 (−0.23) | 0.08 (−0.10) | −0.39 (−0.40) |

- **Verdict: does not track the SMC/perivascular programme — 4/4 cohorts, both sample sets.**
  Only one individual marker in one cohort reaches the rule's threshold (TAGLN, RNA-seq 0.52).
  The SMC composite is 0.01–0.45 across cohorts and 0.03 in LCM. Partial correlations are similar
  (e.g. discovery SMC composite 0.28; RNA-seq 0.48), so this is not a global-expression artefact.
- **MYOZ2 is not even the panel's SMC proxy:** PPP1R3C has higher SMC-composite correlations
  (ERCB 0.51, RNA-seq 0.47) and TRPC6 has higher podocyte-composite correlations (ERCB 0.59,
  RNA-seq 0.51, LCM 0.62) — see `09g_moyz2_panel_reference.csv`.
- **A modest podocyte-content association exists but is largely disease-driven.** NPHS2 rho rises to
  0.53 (RNA-seq) and 0.81 (LCM); the LCM value is robust to leave-one-out (0.79–0.85) but comes from
  the MN-vs-control contrast: within LCM MN only it is 0.40 (95% CI −0.23 to 0.79, n = 12). RNA-seq
  shows the same attenuation (0.53 all → 0.37 MN-only, q = 0.065). Do not present this as a
  within-disease coupling.

## 2. Cell level — snRNA subset vs scRNA

Read from the per-sample QC objects (the atlas HVG assay does not carry MYOZ2), keyed to
`work/scrna/atlas_cells.csv`: 232,458 matched cells (snRNA 210,879 from GSE271188; scRNA 21,579 from
GSE171458 + GSE241302).

- **Detection is rare:** snRNA 1,686 / 210,879 cells (0.80 %); scRNA 57 / 21,579 (0.26 %).
- **Detection is dominated by sequencing depth, not cell identity:** median total counts 4,874 in
  MYOZ2+ vs 2,971 in MYOZ2− cells; logistic OR 8.69 per log10 counts (95% CI 7.47–10.10,
  p = 2.4e-172).
- **No cell-type enrichment after depth adjustment.** Raw snRNA detection is 0.65–1.08 % across all
  14 represented cell types; the highest raw enrichment is Podocyte (4.3×, exact p = 0.006,
  BH q = 0.097, 64 positive cells). Adjusted for log counts + cell type (ref PT), cell-type ORs are
  0.76–1.24; only DCT (0.82, p = 0.03) and PEC (0.76, p = 0.01) are lower, none higher.
- **Co-detection with podocyte markers survives adjustment, weakly.** Adjusted for depth and cell
  type, MYOZ2+ cells have OR 2.11 (NPHS1) and 2.58 (NPHS2); SMC markers are 1.12 (ACTA2, ns),
  1.28 (MYH11), 1.30 (TAGLN), 1.77 (PDGFRB); fibroblast/endothelial markers are ns. Absolute numbers
  are small (1,686 positive cells, 64 podocyte-positive), so this is a hypothesis, not a localisation.
- **Donor-level snRNA pseudobulk (37 donors):** MYOZ2 CPM correlates with TAGLN 0.54, RGS5 0.50,
  MYH11 0.43, ACTA2 0.36, PDGFRB −0.02, NPHS1 −0.01, NPHS2 −0.06, COL1A1 −0.09. The six MN donors
  alone give TAGLN 0.94 and RGS5 0.94 (p = 0.005 each) but n = 6; descriptive only.

## 3. Statement this supports

> MYOZ2's bulk signal is **not** explained by perivascular / smooth-muscle / mesangial tissue content
> (0/4 cohorts; rho 0.01–0.45 for the SMC composite), and MYOZ2 is **not assignable to a kidney cell
> type**: it is detected in ~0.8 % of single-nucleus profiles, in every cell type at a similar rate,
> with detection dominated by sequencing depth and no enrichment after adjustment. A weak,
> partly disease-driven association with glomerular podocyte content (NPHS2 / podocyte composite) is
> present, but it does not establish a podocyte source. The cellular source remains unresolved.

What would change the statement: the pre-specified IF plan (`Manuscript/IF_validation_plan.md`,
MYOZ2 first), or deeper/probe-based single-cell measurement. Neither is available now.

## Artefacts

`output/tables/09g_moyz2_bulk_correlations.csv`, `09g_moyz2_bulk_summary.csv`,
`09g_moyz2_panel_reference.csv`, `09g_moyz2_cell_detection.csv`,
`09g_moyz2_cell_codetection.csv`, `09g_moyz2_cell_adjusted.csv`, `09g_moyz2_depth_check.csv`,
`09g_moyz2_donor_correlation.csv`;
`output/models_panel/09g_moyz2_{bulk_correlations,bulk_scatter,cell_detection}.{pdf,jpeg}`.
Resume cache: `work/scrna/09g_moyz2_cell_matrix.rds` (rebuild with `--refresh`).
