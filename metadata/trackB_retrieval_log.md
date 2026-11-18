# Track B retrieval log — B15 registry expansion (2026-11-18)

**Goal:** raise the published-panel registry from 6 to ≥ 10 so the attrition ladder is not a
six-panel pilot (B15).

**Scope rule (code-enforced, P-70):** a panel must name **≥ 2 genes**, carry a **DOI** (PMID
preferred), report **performance**, and be derived in **human kidney tissue** (bulk microarray or
RNA-seq). Rejections are recorded with a reason in `metadata/trackB_exclusions.csv`.

**Sources:** PubMed E-utilities via the PubMed MCP server; full text from PMC (structured JATS).
Europe PMC search returned HTTP 503 on 2026-11-18 and contributed nothing.

## Queries run (PubMed, publication dates 2018/2019–2026)

| # | Query (abridged) | Hits |
|---|---|---|
| A | `"membranous nephropathy"[TIAB]` AND (diagnostic OR diagnosis OR biomarker) AND (LASSO OR WGCNA OR "machine learning" OR "random forest" OR "SVM-RFE" OR nomogram OR "gene signature" OR "diagnostic model") AND (microarray OR transcriptome OR "RNA-seq" OR GEO) | 25 |
| B | `"membranous nephropathy"[TIAB]` AND ("hub genes" OR "immune infiltration" OR ferroptosis OR cuproptosis OR pyroptosis OR m6A OR ceRNA) AND (ROC OR AUC OR diagnostic OR diagnosis) | 27 |
| C | `"membranous nephropathy"[TIAB]` AND ("diagnostic model" OR nomogram OR "risk score" OR "gene signature" OR "mRNA signature" OR "biomarker panel" OR "diagnostic signature" OR "diagnostic biomarker") | 62 (60 reviewed) |
| D | `"membranous nephropathy"[TIAB]` AND ("differentially expressed genes" OR DEGs OR "expression profile" OR transcriptomic*) AND ("logistic regression" OR "least absolute shrinkage" OR LASSO OR "support vector machine" OR "random forest" OR "classification model") | 10 |

Overlap was heavy; the union was screened by title/abstract and full text where a panel was plausible.

## Decisions

**Included (P07–P12):** Luo 2026 (PMID 42021359), Fan 2025 (40539070), Gui 2024 (39524772),
Shao 2025 (40564208), Hui 2025 (40489038), Zhao 2026 (42801211). Gene lists, discovery/validation
cohorts and the performance text as reported are in `metadata/trackB_panel_register.csv`.

**Excluded (full reasons in `metadata/trackB_exclusions.csv`):**
- Wu 2026 (41948325) — 8-gene signature is pan-immunonephropathy, not MN-specific.
- Dong 2022 (36105034) — 3 genes but no DOI and no retrievable full-text body.
- Wu 2019 (31160422) — lipid-metabolism hub genes, no diagnostic model/AUC.
- Wang 2026 (41638545), Zee 2025 (40019793) — clinical/clinicopathological models, not expression panels.
- Sealfon 2022 (35477557) — molecular characterisation, not a diagnostic panel.
- Yang 2025 (42147989) — ITP blood B cells; MN only via GeneCards overlap.
- SPP1 (41080598), CALR (40568201), HOPX (42309401), RACGAP1 (39068752), Ets2 (35372448),
  ITGB2 (40697036), SLC1A5 (41189572) — single-gene studies, below the 2-gene threshold.

## Measurability check (before registration)

Against `work/bench/10b_harmonised.rds` (13,379 shared genes): Luo 3/3, Fan 6/7, Gui 2/2,
Shao 4/4, Hui 4/4, Zhao 1/2 (TLR10 absent). Absent genes are recorded in each row's `notes`;
10c uses the measurable subset and reports `n_genes_used`.

## Limits

This is a PubMed/PMC retrieval, not a systematic review. Paywalled supplements and publisher-only
tables were not accessible; no candidate that clearly meets the rule is known to remain behind one.
