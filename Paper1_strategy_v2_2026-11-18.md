# Paper 1 — Analysis Strategy v2 (2026-11-18)

**Supersedes:** `MN_Paper1_Bulk_ML_Biomarkers_Subtypes.docx` (dated 2026-10-15, accessions verified
2026-10-15). That document remains the *historical* plan; this one records what the completed analysis
taught us, what the field has published since, and how to make the programme genuinely novel.
**Companion state document:** `handoff.md` (§2.19–§2.23, §3, §4). **Evidence ledger:**
`Docs/Findings_Ledger.csv`.

---

## 0. Why the strategy needs revising (not because the analysis failed)

The 2026-10-15 plan assumed novelty from **scale** ("largest assembled MN resource") and **design**
("cross-platform, cross-compartment"). Both are now delivered — 13 accessions, 2,165 samples, 178
unique patients, four independent sample sources, 12 figure plates. What the completed work has
actually produced, and what the plan did not anticipate, is a *methodological* result:

> **The external AUC of an MN classifier depends more on how the validation cohort is standardised
> than on which genes are chosen.** Removing a transductive z-score and then *not* replacing it with a
> legitimate cross-platform harmonisation dropped the gate's ERCB AUC to 0.64–0.77 and inverted four
> comparator diagnoses; adding rank-preserving quantile mapping restored 0.861–0.922
> (`handoff.md` §2.20, F62).

Meanwhile the 2024–2025 MN literature continues to publish `limma → WGCNA → LASSO → AUC > 0.9` on the
same accessions, almost always without patient de-duplication, without a stated standardisation
convention, and sometimes "validating" in `GSE133288` — a NEPTUNE tubulointerstitial series — against a
discovery set assembled from NEPTUNE glomerular series, i.e. the same patient population. Two recent
examples illustrate the pattern precisely (Ren Fail 2025, PMID 39806779, discovery GSE115857+GSE200828,
"validation" GSE133288; Front Immunol 2025, PMID 41132671, discovery GSE108109, panel genes CCND1 /
CEBPD / COL10A1 / BMP2).

**That is the gap worth occupying.** The programme's differentiator is not another panel — it is the
demonstration, on the field's own data, of how much reported MN biomarker performance survives
rigorous validation, plus the design rules that make it survive.

---

## 1. Field audit and gap register

### 1.1 What the current literature does (2023–2026)

| Element | Typical practice | Consequence |
|---|---|---|
| Cohort | 1–2 GEO accessions, often NEPTUNE arrays | Small, redundant, single-population |
| De-duplication | never done | 2–3× patient reuse inflates n and narrows CIs |
| Platform handling | not stated | AUC depends on an unreported choice (see §0) |
| Feature selection | LASSO on the full cohort, then ROC on the same cohort | optimistic bias; no outer loop |
| "Validation" | a second accession from the same consortium | not independent |
| Panel size | 4–12 genes, chosen by CV curve | unstable: our own 19k model used 35–198 coefficients across partitions |
| Cell composition | CIBERSORT LM22 (a **blood** signature) on **kidney** tissue | mis-specified reference |
| Endpoint | diagnosis, occasionally PLA2R correlation | no outcome, no protein, no genetics |

### 1.2 Gap register (each gap is addressable with existing data unless marked)

| ID | Gap | Why it is a gap (evidence) | Fix |
|---|---|---|---|
| **G1** | **No leakage-audited benchmark of published MN panels.** Nobody has measured how much of the repeatedly reported AUC > 0.9 survives de-duplication + honest external validation | 2024–25 papers cite AUC > 0.9 on the same accessions we show is platform-sensitive | **Track B, the flagship** (§3) |
| **G2** | Patient-level redundancy ignored | our resource: 453 MN samples = **178 patients** (2.7 samples/patient) | already solved here; becomes a benchmark axis |
| **G3** | Cross-platform standardisation never stated | F62: 0.64–0.77 (invalid) vs 0.861–0.922 (harmonised) | harmonisation as a stated, tested step |
| **G4** | Compartment conflation | our panel: glomerular 0.816, **tubulointerstitial 0.449** | compartment-specific claims only |
| **G5** | **Antigen-defined MN is not stratified.** PLA2R-negative (~20–30 %), THSD7A, NELL1, HTRA1, SEMA3B, NCAM1 define clinically distinct disease | our PLA2R strata: 0.976 / 0.881, but n = 6 + 6 — underpowered | antigen-aware design; request antigen status |
| **G6** | Deconvolution uses blood signatures on kidney; no kidney-specific Bayesian benchmark in MN | our two signatures disagree on the dominant cell type in **1,314/1,328** samples | BayesPrism/InstaPrism with an integrated sc+sn kidney reference (§4, verified tools) |
| **G7** | **No protein or genetic anchor.** Every MN panel stops at tissue mRNA | our own E4 = 0 in the ledger | pQTL/UKB-PPP anchoring + MR (§4); IF in-house |
| **G8** | No outcome-linked tissue biomarker (eGFR/UPCR/remission) | NEPTUNE clinical variables access-controlled (B2) | data request; then Cox/ΔAUC over clinical model |

**Orthogonal gap (opportunistic):** no **spatial** validation of glomerular module localisation in MN.
Spatial kidney atlases now exist; a 3–5 sample spatial run would convert "podocyte-rich by
deconvolution" into "signal localised in situ". 【Spatial platform access to be confirmed】

---

## 2. Repositioned programme

| Track | Question | Deliverable | State |
|---|---|---|---|
| **A** — current Paper 1 | Is there a reproducible, compartment-resolved, externally validated MN signature? | the manuscript as it stands (EBioMedicine format) | **submission-ready** (computational side) |
| **B** — Paper 1b (new, flagship) | *How much of the published MN biomarker performance survives rigorous validation — and what design rules make a panel survive?* | benchmark + attrition table + design checklist | starts from existing data + published gene lists |
| **C** — Paper 2 | Do the panel genes *cause* MN, and can the modules be drugged? | MR/colocalisation + LINCS repurposing + IF | needs GWAS/pQTL access + wet lab |
| **D** — Paper 3 | Does the tissue signature add to clinical prognostication? | outcome model (Cox, ΔC-index, DCA) | blocked on B2 |

**Do not** fold Track B into Track A. Paper 1 is a discovery paper; Track B is a meta-research paper
that *uses* Paper 1's pipeline as its instrument. Submitting them separately doubles the contribution;
merging them makes both halves thinner and invites "the authors benchmarked themselves".

---

## 3. Track B — the flagship design (this is the novelty)

### 3.1 Question
Of the diagnostic performance reported for published MN transcriptomic panels, what fraction survives
(a) patient de-duplication, (b) a stated cross-platform harmonisation, (c) nested cross-validation,
(d) a truly independent cohort, and (e) a refit-based permutation null?

### 3.2 Material
- **Published panels** to reproduce: minimum 8, target 12, from 2020–2026, each with a gene list and a
  reported AUC. Candidate sources: the two 2025 papers above, plus prior PLA2R-association, immune-
  infiltration and WGCNA–LASSO MN papers. 【Selection rule: paper reports ≥3 named genes AND an AUC;
  every included panel gets an ID and a row in the register — no cherry-picking】
- **Cohorts**: our existing 7 non-redundant sets (discovery NEPTUNE glomerular; ERCB glomerular and TI;
  two RNA-seq; LCM; alternative-CDF reprocessing).
- **Nothing new needs downloading** for the core benchmark.

### 3.3 The attrition ladder (the paper's central figure)
For every panel, compute AUC at each rung and report the **drop**:

| Rung | What changes | Expected direction |
|---|---|---|
| R0 | as published (their cohorts, their statement) | reported AUC |
| R1 | our cohorts, no de-duplication, cohort-internal z | reproduces ≈ R0 |
| R2 | + patient de-duplication | ↓ (redundancy removed) |
| R3 | + nested CV / out-of-fold selection | ↓↓ (selection optimism removed) |
| R4 | + training-statistics + quantile-mapped harmonisation | ↓ or ↑ (platform artefact removed) |
| R5 | + independent ERCB / RNA-seq cohort | ↓↓↓ (the real number) |
| R6 | + refit permutation null | significance floor |

**Pre-registered hypothesis:** median R5 AUC across published panels < 0.70, i.e. most reported MN
panel performance does not survive; and the R3→R4 step is the single largest unexplained variance
source in the literature. Both are falsifiable and both are publishable either way.

### 3.4 Secondary analyses (cheap, high value)
- **Panel-size frontier**: nested CV at K = 2, 3, 4, 6, 8, 12, 20 → external AUC. Answers "why four
  genes when 19k scores higher?" for our own panel *and* for the field.
- **Design checklist**: derive 8–10 rules from the ladder (de-duplicate; state the standardisation;
  outer loop for selection; independent cohort; permutation null; compartment declared; report CI;
  report n patients not samples). Publish as a boxed table — that is what reviewers and editors cite.

---

## 4. Tool and method upgrades (2023–2026), with evidence

Verified in this session unless marked 【verify】. **Rule: adopt a tool only after it beats the current
method on our own data; a newer tool is not a better tool.**

| Layer | Current in Paper 1 | Candidate (year) | Evidence / why | Adoption test |
|---|---|---|---|---|
| Tabular classifier | elastic net / logistic | **TabPFN v2** (Nature 2025, PMID 39780007) — tabular foundation model, no tuning, handles missing/outliers | strong on small-n tables; **but** the TabBench-Bio biomedical benchmark (43 datasets) finds **logistic regression statistically indistinguishable** from TabPFN v2.5 | 5-fold CV on discovery with identical folds; adopt only if ΔAUC > 0.02 **and** calibration no worse. Report the null result if it ties — that is itself a finding |
| Deconvolution | NNLS on a filtered snRNA signature | **BayesPrism / InstaPrism** (Bioinformatics 2024, PMID 38970377; BayesPrism top-ranked in the heterogeneous-pseudobulk benchmark, Genome Biol 2024, PMID 38956606) | Bayesian, returns uncertainty; InstaPrism = fast deterministic equivalent | re-deconvolve all 7 sets; require the podocyte result to hold and the two-method disagreement (currently 1,314/1,328) to shrink |
| Reference construction | single snRNA atlas | **integrated sc + sn reference** (2026, PMC13106970) | snRNA and scRNA differ by cell type; integrating improves bulk deconvolution | build both references, compare on held-out pseudobulk |
| Pathway / TF activity | GSVA, GSEA only | **decoupleR** (R + Python; benchmark: consensus / ULM / multivariate linear model best) | turns modules into *regulator activity* rather than gene-set overlap | run on our 584 meta-DEGs; require top regulators to be concordant across cohorts |
| Drug repurposing | not done | **LINCS L1000 connectivity via `signatureSearch`** 【verify version】 | maps module signature to perturbagens with reversal score | report top 20 with MoA; must be contrasted against the drugs MN actually uses (rituximab, cyclophosphamide, CNIs, rituximab-alternatives) |
| Causal anchor | not done | **MR + colocalisation** (TwoSampleMR, MR-PRESSO, coloc 【verify versions】) with eQTL (eQTLGen/GTEx) and **pQTL (UKB-PPP, deCODE)** | converts "associated" into "genetically supported"; directly addresses G7 | 4 genes × cis-eQTL/pQTL → MN GWAS; colocalisation PP.H4 > 0.8 to claim |
| Protein bridge | wet lab only | **pQTL-anchored protein levels** in kidney-disease GWAS 【verify availability】 | a computational E3 step before IF | report as E3 only if colocalised |
| Annotation (scRNA) | marker score | foundation-model / LLM annotators (scGPT, CellPLM, GPTCelltype) 【verify】 | faster, but no benchmark beating marker scores on kidney | only if re-annotation is ever needed (B8 currently accepted) |
| Spatial | none | cell2location / Tangram 【verify】 | localises modules in situ | 3–5 MN biopsies 【access to confirm】 |

**Anti-pattern to avoid:** stacking new tools as decoration. Each row above must either change a
number in the ledger or be reported as a null.

---

## 5. Revised pipeline (what actually changes)

### Phase 0 — Freeze Track A
No further analysis on Paper 1 except reviewer response. Its claims are E2-graded, its figures are
submission-ready, and its verifier suite is green (54 / 86 / 42 / 117 / 106 / 32). New methods go into
Track B, not into Paper 1 — otherwise the submitted paper and the repository diverge.

### Phase B1 — Benchmark scaffolding (≈1 week, no new data)
1. `10a_panel_register.csv`: every published panel — paper ID, DOI, reported AUC, cohort, panel genes,
   selection method, whether de-duplicated, whether standardisation stated.
2. `10b_harmonised_matrices.R`: reuse `05a`'s `quantile_map()` + `std()`; export harmonised ERCB/RNA-seq
   matrices so the benchmark and the gate share one transform (P-63: one implementation, not two).
3. `10c_panel_reproduction.R`: implement the attrition ladder (§3.3) for every panel; write
   `10c_attrition.csv` with one row per panel × rung.
4. Verifier `99_verify_trackb.R`: every panel has 7 rungs; no rung missing; the ladder is monotone in
   the expected direction or the deviation is explained.

### Phase B2 — Method upgrades with adoption tests (§4)
`10d_tabpfn_benchmark.R`, `10e_bayesprism_deconvolution.R`, `10f_decoupleR_activity.R`,
`10g_lincs_repurposing.R`. Each writes a comparison table against the incumbent and a one-line verdict
(adopt / reject / tie). **Ties are reported, not hidden** — that is the honesty the field lacks.

### Phase B3 — Antigen-aware extension
Extend the panel evaluation to antigen strata (PLA2R+ / PLA2R− now; THSD7A/NELL1/HTRA1 when status is
obtained). Pre-specify that within-MN contrasts need n ≥ 20 per group to be interpretable — our 6 vs 6
is reported as uninterpretable and stays that way.

### Phase C — Causal and protein anchors (needs access)
MR + coloc for the four genes; pQTL anchoring; then IF (MYOZ2 first) — the existing
`Manuscript/IF_validation_plan.md` is unchanged and remains the E4 gate.

---

## 6. Novelty claims, and the experiment that could kill each

| # | Claim | Killer experiment | Current status |
|---|---|---|---|
| N1 | Published MN panel AUCs do not survive leakage-free validation | the attrition ladder itself — if R5 ≈ R0, the claim dies | untested, cheap, pre-register the hypothesis |
| N2 | The standardisation convention, not gene choice, dominates external AUC | R3 vs R4 on 12 panels | supported by F62 (single model, 0.64 → 0.92) |
| N3 | MN signal is glomerular-specific | TI transfer (0.449) | **done, E2** |
| N4 | A four-gene panel transports across four sources | ERCB/RNA-seq/LCM/alt-CDF | **done, E2**; refit null p = 0.11 is the honest limit |
| N5 | Panel genes are causally implicated | MR/coloc | untested |
| N6 | Foundation models do not beat penalised logistic regression on n ≈ 150 transcriptomes | TabPFN benchmark | plausibly a *null* — publish as such |

**Do not claim:** that our panel is better than published panels (it may not be — the benchmark will
tell); that MYOZ2 is podocyte-derived (unresolved); that any subtype exists (negative, E1).

---

## 7. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| Benchmark is perceived as attacking other groups | medium | frame as *field-level* method audit; cite panels neutrally by DOI; invite authors to contribute panels |
| Panels cannot be reproduced (gene lists incomplete) | medium | inclusion rule requires ≥3 named genes + an AUC; report exclusions with reason |
| A published panel *does* survive → N1 weakens | medium | the paper still delivers the design rules and the R3/R4 variance result; pre-register both outcomes |
| MR yields null | medium-high | null MR is reportable and honest; it reframes the panel as a diagnostic, not causal, marker |
| Reviewer demands clinical outcome data | high | B2 request must be submitted early (see §8) |
| Tool churn consumes time with no gain | medium | every adoption requires a pre-specified ΔAUC margin (§4) |

---

## 8. Order of work

| Priority | Action | Needs | Effort |
|---|---|---|---|
| **P0** | Submit Paper 1 as it stands; freeze its analysis | user sign-offs in `Manuscript/SUBMISSION_CHECKLIST.md` | — |
| **P0b** | Submit the B2/B7 data request (NIDDK R4R, DOI `10.58020/f2ae-n094`) | user | 1 afternoon |
| **P1** | Build the panel register + attrition ladder (Track B core) | nothing new | 1–2 weeks |
| **P1b** | Adoption tests: TabPFN, BayesPrism/InstaPrism, decoupleR | nothing new | 1–2 weeks |
| **P2** | LINCS repurposing on the module signature | nothing new | 3–5 days |
| **P2b** | MR + colocalisation (needs GWAS/eQTL/pQTL access) | access 【confirm】 | 2–3 weeks |
| **P3** | Spatial validation (3–5 biopsies) | tissue + platform 【confirm】 | collaboration |
| **E4** | IF of MYOZ2 → TRPC6 → PPP1R3C → SLAMF8 | in-house imaging | ongoing |

---

## 9. What this strategy does **not** change

- The data resource, harmonisation and QC-first discipline (`bioanalyst` gates, `Docs/BAP_AUDIT_2026-11-18.md`).
- The four-gene panel as Paper 1's result, with its honest limits (refit null p = 0.11; below the
  19k-gene model at 0.816 vs 0.922; needs per-platform recalibration).
- The evidence grading: nothing becomes E3/E4 without protein or genetic support.

---

*Prepared 2026-11-18. Literature and tool claims marked 【verify】 were named from domain knowledge and
were not confirmed online in this session; confirm version and licence before implementation.*
