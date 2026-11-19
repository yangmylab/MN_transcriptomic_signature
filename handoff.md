# Handoff — Project Handoff Document

**Project:** machine-learning biomarker discovery for membranous nephropathy (MN) from multi-cohort bulk transcriptomics, plus Track B, a leakage-audited benchmark of published MN panels, and a public code archive with a Zenodo DOI.

- **Project root:** `/run/media/meng/MY_Data/MN/Paper_1`
- **Shared data root:** `/run/media/meng/MY_Data/MN/Datasets`
- **Git branch:** `main` (renamed from `master` this session), tracks `origin/main`
- **Remote:** `git@github-mn:yangmylab/MN_transcriptomic_signature.git`
- **Archive DOI:** `10.5281/zenodo.24358387` (version); concept DOI `10.5281/zenodo.24358386`
- **Last updated:** 2026-11-19 (evening)
- **Current stage:** Paper 1 analyses frozen and re-verified; manuscript recomposed, humanised, references rebuilt with PMID/DOI; code repository public and archived; Track B computed but scope decision outstanding.

---

## 1. Objective (What we are doing)

- Deliver a submission-ready manuscript reporting a compartment-resolved glomerular transcriptomic signature and a four-gene classifier for MN, with every claim externally validated.
- Preserve the audited Paper 1 claim set: interpretability and transportability, not superiority over the whole-transcriptome gate; no prognostic claim; no unsupported subtype or cell-of-origin claim.
- Publish the analysis pipeline openly (MIT, GitHub) and archive a citable release (Zenodo DOI), so the Data sharing statement names a real repository rather than "on request".
- Complete Track B, a leakage-audited benchmark of published MN panels, with explicit attrition, external transfer and refit-permutation analyses.
- Success criteria: seven passing verifiers, frozen Paper 1 results, a logically ordered manuscript and figure package, a public archived repository, and submission declarations completed by the user.

## 2. Completed Items (What has been done)

### 2.1 Paper 1 analysis chain — complete and frozen

Verified chain:

**13 accessions → 2,165 samples → 453 MN samples → 178 unique MN patients → 7 non-redundant cohort sets → Phase A 54/54 → Phase B 86/86 → separability gate passed → subtypes negative → 432,378 cells → 236,107-cell atlas, 16 cell types → deconvolution 1,328 samples, 67/316 significant → preserved WGCNA modules → four-gene panel, ERCB AUC 0.816 → recalibrated Brier 0.096.**

Script-generated values that must not be re-quoted from memory:

| Block | Current verified value |
|---|---|
| Design | 13 accessions; 2,165 samples; 453 MN samples; **178 unique patients**; 7 cohort sets |
| QC | 200 flagged samples retained as annotation only |
| Glomerular meta-DEGs | Control 584; FSGS 18; MCD 60; Other 16; non-MN 18 |
| TI meta-DEGs | Control 323; FSGS 2; MCD 4; Other 7; non-MN 1; TI control has *k* = 1 |
| DL fallback | **0**; never resurrect the incorrect "724 genes" statement |
| ERCB replication | 55.3% glomerular / 32.4% TI for MN versus control |
| Gate | Nested-CV AUC 0.911–0.926; external ERCB 0.861–0.922; RNA-seq 0.928–0.939 |
| Gate robustness | Permutation null median 0.500, 95th 0.674, *p* = 0.048; QC-removed 0.939; reverse 0.772 |
| Gate partition stability | *K* = 1000: 0.9003–0.9032; all: 0.9118–0.9281 |
| Subtypes | Negative: array *k* = 2, PAC 0.00077 but technically confounded; RNA-seq minimum PAC 0.211; 1 shared defining gene, Jaccard 0.0011 |
| scRNA atlas | 432,378 retained cells; 236,107 × 2,019 genes; 16 cell types |
| Deconvolution | 1,328 samples; 0 degenerate fits; 67/316 significant comparisons |
| Signature fragility | Filtered and all-gene signatures disagree on dominant cell type in 1,314/1,328 samples |
| Podocyte result | Rich versus other nephropathies, adjusted *p* = 3.3e-4; not versus healthy glomeruli, *p* = 0.29 |
| Modules | 9 glomerular (power 18) + 7 TI (power 10); all real modules preserved; grey and gold fail |
| MN modules | Pink 84 genes, adjusted *p* = 1.0e-5; red 144 genes, adjusted *p* = 6.6e-4 |
| Candidate space | 708 hubs → 515-gene core; 3,279 with provenance |
| Locked panel | **SLAMF8, TRPC6, MYOZ2, PPP1R3C** |
| Candidate accounting | 515 core → 353 glomerular → 66 excluded for <95% ERCB coverage → **287 tested** |
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

The gate outperforms the panel in ERCB (0.922 versus 0.816); the panel's role is compact interpretability and transportability, not superiority.

### 2.2 Manuscript — recomposed this session

The manuscript is **DOCX-only**: `Manuscript/Paper1_manuscript.docx`. `Manuscript/Paper1_manuscript.md` was deleted at the user's request on 2026-11-19; there is **no Pandoc source any more**, so all edits are made directly in the DOCX with python-docx and validated with OfficeCLI.

Current state:

- **19 pages**, 13 figure plates, **65/65 panels** cited, 44 references.
- Results has 12 subsections. Citation/interpretation order was corrected this session to follow the analysis pipeline:
  - A new subsection **"Sample quality control and batch harmonisation"** now separates design (Fig 1A–C) from QC and ComBat (Fig S1A–D, S2A–B, then Fig 1D–E).
  - The **kidney atlas section moved before "Cellular localisation"**, because the atlas was previously used 12 paragraphs before it was introduced.
  - The **stability/hub paragraph moved before global pathway enrichment**, removing a duplicated `S3B` citation and the S4-before-S3 inversion.
  - `S2D/S2E` now come after Fig 6E instead of splitting Fig 6A from 6B.
  - Resulting first-appearance sequence is monotonic per plate: `1A-C → S1A-D → S2A-B → 1D-E → 2A-F → 3A-E → S2C → 4A,S3A,4B-4E → S3B-E → S4A-E → 5A-E → S5A-E → 6A-E → S2D-E → S6A-E → S7A-E`.
- **Prose humanised**: 54 paragraphs rewritten to remove AI-style tells (all 20 em dashes removed, rhetorical fragments and ordinal scaffolding removed, over-claiming verbs downgraded). Every number (1,031 numeric tokens), every citation bracket, and every figure citation was verified byte-identical afterwards; only citation numbers changed where references were renumbered.
- **Author block added** (user-pasted): Mingyue Yang^1,2,#; Lingzhang Meng^3,#; Yibing Wang^1; Kun Ye^2; three affiliations; equal-contribution and correspondence footnotes. Corresponding authors are **not marked with `*`** on the author line.
- **Back matter**: Authors' Contributions (the drafted text, used intact), Funding (Guangxi Natural Science Foundation, `#桂科AB24010219`), Data sharing statement, Declaration of interests, Declaration of generative AI use.
- **Abstract** (renamed from Summary): Background / Methods / Results / Interpretation, **361 words excluding keywords** (374 incl. keywords); `Results` is used instead of `Findings`; **no Funding line**.
- **Removed sections**: `Research in context`, `Ethics`, `Acknowledgements`.
- **Reference list rebuilt** (see 2.4): 44 entries, strict citation order 1–44, all with DOI and PMID. `Manuscript/references.tsv` was restored from git and rewritten in final order.
- **Figure plates** now carry numbered titles: `Figure 1. Cohort design and analytical independence.` through `Figure 6. …`, and `Figure S1. …` through `Figure S7. …`.

Known manuscript defects still open (details in §3.4):

- `officecli validate` reports **80 schema errors**, all in `styles.xml` (`w:uiPriority` written after `w:qFormat`), introduced when the DOCX was saved by a word processor. The document renders correctly and structural checks pass, but the project's "validate = no errors" gate is broken.
- One empty `Heading 1` (left from the deleted `Contributors` heading) and two empty Body Text paragraphs.
- Back-matter headings are Body Text, not heading styles.
- Typos retained in the user's front matter: `contrinbuted`, `should be address to`, `Big Data And Artificial Intelligence`, stray space before `)` in the correspondence line.
- Data sharing statement as currently pasted by the user: `underhe MIT licence`, a double space before `available`, and a trailing space.
- Reference preamble still says the list is generated from `Manuscript/references.tsv`; that is true again now that the file exists.

### 2.3 Figures — regenerated with numbered titles

- `tools/make_figures.py` now carries a `label` per plate (`Figure 1.` … `Figure S7.`) and renders `label + title + "."`.
- All **13/13 plates rebuilt**, **one page each**, fonts embedded (`pdffonts` reports `emb = yes`), PDF + 300-dpi JPEG.
- `output/figures/` and `output/figures/supplementary/` hold the authoritative plates; `Manuscript/Figures/` holds **26 byte-identical copies** (13 PDF + 13 JPEG).
- Panel content, panel letters and layout untouched; only the title line changed.

### 2.4 References — audited, refreshed, and made sequential

Starting point: 44 entries, 44 DOIs, **42 PMIDs** (`Benjamini–Hochberg 1995` and `Boruta 2010` were the two without, both genuinely outside PubMed), only **2** papers from 2023+ (`Hao 2024`, `Lake 2023`), and one ordering defect (`[42]` first appeared before `[40]`/`[41]`).

Work done:

- Recovered `Manuscript/references.tsv` from git (`git show HEAD:Manuscript/references.tsv`), mapped **44/44 to the DOCX entries by DOI** (the TSV row order does not equal the manuscript order, so DOI is the only safe key).
- Verified every DOI against PubMed E-utilities and OpenAlex. 42/44 had PMIDs; `glmnet` (20808728) and `DeLong` (3203132) were confirmed by PMID/title because their DOI forms are not indexed.
- **Replaced 8 entries** (clinical/review only, plus the two non-PubMed method citations), all at their existing positions:

| Ref | Replaced | With | Year · journal |
|---|---|---|---|
| 1 | Couser 2017, CJASN | Sethi & Fervenza, antigens in MN | 2025 · Nat Rev Nephrol |
| 2 | Ronco 2021, Nat Rev Dis Primers | Mayo Clinic consensus report on MN | 2023 · Kidney Int |
| 13 | Beck 2011, JASN | Barbour et al., treatment nonresponse | 2023 · Clin J Am Soc Nephrol |
| 15 | Hoxha 2012, Kidney Int | Cai et al., podocyte BMP-2 pathway in PMN | 2024 · Adv Sci |
| 16 | Ronco 2015, Lancet | Clair et al., spatial glomerular transcriptome | 2024 · JCI Insight |
| 42 | D'Agati 2011, NEJM (FSGS) | Romagnani et al., Podocytopathies | 2025 · Nat Rev Dis Primers |
| 26 | Benjamini–Hochberg 1995 | Glickman et al., FDR control | 2014 · J Clin Epidemiol |
| 36 | Kursa & Rudnicki (Boruta) 2010 | Degenhardt et al., variable selection evaluation | 2019 · Brief Bioinform |

- **Renumbered** old 42→40, 40→41, 41→42 so first appearance is strictly 1→44. Only two sentences changed, in their brackets: `…clinical difficulty [40]` and `…atlas (16 cell types) [33,41,42]`.
- Rebuilt all 44 entries in Vancouver style with `doi:… PMID: …`; 44/44 now carry both identifiers.
- Year distribution is now **8 papers from 2023–2025**.
- Caveats: the BH and Boruta citations are now PubMed-indexed secondary sources (the originals are not indexed); the "leading cause of ESKD" clause of sentence `[1,2]` is not stated explicitly by either new source; JCI Insight (ref 16) is IF ≈ 9 and may not satisfy a strict "high-impact" filter.

### 2.5 Code repository and archive — live

- **Repository:** https://github.com/yangmylab/MN_transcriptomic_signature — **public**, **MIT**, default branch `main`.
- **Archive DOI:** `10.5281/zenodo.24358387` (version); concept DOI `10.5281/zenodo.24358386`; Zenodo record title, MIT licence and all four authors verified via the Zenodo API.
- **382 tracked files**, 22 top-level entries, repository `.git` ≈ 372 MB. Largest file 2 MB (well under GitHub's 100 MB limit).
- **Included:** `R/` (52 scripts), `tools/`, `metadata/`, `Manuscript/` (DOCX, `references.tsv`, 26 figure copies), `Docs/`, `output/tables/` (240 result CSVs), `output/figures/` (13 plates), `environment/` (`renv.lock`, `sessionInfo.txt`), `reference/`, `LICENSE`, `README.md`, `CITATION.cff`, `.zenodo.json`.
- **Excluded by `.gitignore`:** `work/` (23 GB), `.Rlib/` (4.2 GB), `scrna/extracted/`, `logs/`, per-step intermediate figures, five bulk CSVs >3 MB (`signature_matrix.csv`, `deg_all_significant.csv`, `scrna_pseudobulk_de.csv`, `08f_module_enrichment_all.csv`, `scrna_detection_fraction.csv`), and now `__pycache__/` and `*.pyc`.
- **Commits this session:** `a26c71b` (manuscript + references + README/Zenodo metadata), `64e61c3` (merge initial MIT LICENSE), `ca938e3` (numbered figure titles), `69cd48f` (ignore bytecode caches). Remote HEAD is `69cd48f`.
- **Authentication:** a dedicated deploy key for this repository lives at `~/.ssh/mn_github_deploy_yang` (private half stays on this machine; the public half is registered on the repo with write access). `~/.ssh/config` has `Host github-mn` pointing to **`ssh.github.com:443`**, because port 22 to GitHub is blocked on this network. An earlier `github-yangmylab` host entry uses a different key (`colorectal_deploy`) for another repository; do not reuse it.

### 2.6 Verification completed this session

| Verifier | Result |
|---|---:|
| `99_verify_phase_a.R` | 54 passed, 0 failed |
| `99_verify_phase_b.R` | 86 passed, 0 failed |
| `99_verify_downstream.R` | 32 passed, 0 failed |
| `99_verify_scrna.R` | 42 passed, 0 failed |
| `99_verify_phase_d.R` | 117 passed, 0 failed |
| `99_verify_phase_e.R` | 106 passed, 0 failed |
| `99_verify_trackb.R` | 37 passed, 0 failed |
| **Total** | **474 passed, 0 failed**, all exit 0 |

Additional verified state: 191 figure PDFs / 191 rasters; 240 result CSVs (121 at `output/tables/` top level); no zero-byte outputs; 65/65 panels cited; 44/44 references cited.

### 2.7 Track B — computed, scope undecided

- **12 published panels**, 15 reasoned exclusions (`metadata/trackB_exclusions.csv`, `metadata/trackB_retrieval_log.md`).
- Ladder R1→R6 complete: **12 × 7 = 84 rows** in `output/tables/10c_attrition.csv`.
- Original six-panel preregistered hypothesis **falsified** (median R5 0.821). Expanded post-hoc registry median R5 **0.674**; must not be presented as confirmation.
- **0/12 panels beat their refit null** (*p* = 0.055–1.000); our panel *p* = 0.11 in the same context.
- Internal CV peaks near *K* ≈ 12; external AUC flat 0.767–0.858 from *K* = 2 to 20; best two-gene frontier SULF2 + MYOZ2 = 0.858.

### 2.8 MYOZ2 provenance — completed as a separate P1b test

- Bulk MYOZ2 does not track the SMC/perivascular programme (SMC composite rho 0.01–0.45).
- snRNA 1,686/210,879 cells (0.80%); scRNA 57/21,579 (0.26%); depth OR 8.69 per log10 counts (7.47–10.10).
- No cell-type enrichment after depth and cell-type adjustment.
- Verdict: **tissue-level reproducibility with unresolved cellular source**; not vascular/SMC, and not "podocyte-derived".

## 3. Current Blockers (Where we are stuck)

### 3.1 P0 — Protein-level validation not performed

- **Symptom:** no IF/protein validation of the four-gene panel; MYOZ2 cellular origin unresolved.
- **Impact:** limits translational claims and biological interpretation.
- **Already tried:** bulk reproducibility, atlas localisation, donor-level pseudobulk, dedicated MYOZ2 provenance test; none established a defensible cell of origin.
- **Owner:** user wet lab. **Resume:** `Manuscript/IF_validation_plan.md` — **this file was deleted on 2026-11-19 at the user's request and is no longer in the working tree**. Recover from git history if needed (`git show <commit>:Manuscript/IF_validation_plan.md`), or reconstruct from this handoff. Priority MYOZ2 → TRPC6 → PPP1R3C → SLAMF8.

### 3.2 P0b — Clinical outcomes and subject-to-sample mapping missing

- **Symptom:** eGFR, UPCR, remission/relapse and a subject→sample map are unavailable in downloaded public files.
- **Impact:** prevents prognostic modelling, subtype–outcome analysis and definitive participant linkage.
- **Already tried:** local metadata review; PLA2R recovered only for GSE216841 (6 positive, 6 negative).
- **Owner:** user. Submit the long-lead NIDDK Central Repository R4R request (study DOI `10.58020/f2ae-n094`, DACC `NEPTUNE-STUDY@umich.edu`), including the subject→sample map. The prepared request document `Manuscript/DATA_REQUESTS.md` was also **deleted on 2026-11-19**; recover from git history if needed.

### 3.3 P1 — Submission metadata incomplete

- **Symptom:** ORCIDs, CRediT roles, funding wording, three suggested reviewers and declaration sign-off remain 【To be confirmed — user】. Exact EBioMedicine caps also unconfirmed.
- **Status:** authors and affiliations are now set (Mingyue Yang, Lingzhang Meng, Yibing Wang, Kun Ye; co-first and co-corresponding stated). `Contributors` text is drafted and in place.
- **Impact:** package cannot be submitted without these.
- **Owner:** user supplies metadata; agent updates the DOCX.

### 3.4 P1 — DOCX defects to fix before submission

- **Schema:** `officecli validate` = **80 errors**, all `w:uiPriority` ordering in `styles.xml`. Mechanical repair (reorder `name → aliases → basedOn → next → link → autoRedefine → hidden → uiPriority → semiHidden → unhideWhenUsed → qFormat → …`) restores the gate without semantic change.
- **Structure:** one empty `Heading 1`; two empty Body Text paragraphs; back-matter headings not styled as headings.
- **Format/compliance:** `Abstract`/`Results` instead of `Summary`/`Findings`; no Funding line in the abstract; `Research in context`, `Ethics` and `Acknowledgements` removed; abstract 361 words (cap unconfirmed, likely 250).
- **Text:** `contrinbuted`; `should be address to`; `Big Data And`; stray space in `(yk@gxams.org.cn )`; `The author declares` (singular); no terminal period in Funding and AI declarations; Data sharing statement `underhe` and double space.

### 3.5 Track B registry decision unresolved

- **Symptom:** user must decide whether to ship the 12-panel benchmark or fund another retrieval pass.
- **Impact:** changes Track B manuscript scope and timing, not the computed results.
- **Resume:** `metadata/trackB_panel_register.csv`, `metadata/trackB_retrieval_log.md`.

### 3.6 Optional packages and deferred items

- Absent: `decoupleR`, `InstaPrism`, Python `tabpfn`. Installed: BayesPrism, MuSiC, hdWGCNA.
- **B13:** `scDblFinder()` lacks `set.seed()` in `07b`; matters only if the atlas is rerun (~25 min). Atlas frozen.
- **B14:** many scripts lack a `# Parallelism: not needed — <reason>` waiver; documentation only.
- **P5:** GSE73953 PBMC/Agilent remains unnormalised; optional branch in `02b`.
- **P7:** Paper_2/Paper_3 not yet wired to the shared dataset manifest.

## 4. Next Steps (What to do next)

1. **P0 — Start IF validation**
   - Owner: user wet lab. Trigger: laboratory scheduling.
   - Resume: reconstruct the plan from git history (`git log -- Manuscript/IF_validation_plan.md`) or from §2.8; target order MYOZ2 → TRPC6 → PPP1R3C → SLAMF8, blinded per-glomerulus quantification with prespecified lineage anchors.

2. **P0b — Send the NIDDK R4R request**
   - Owner: user. Trigger: immediate (long lead time).
   - Resume: reconstruct `Manuscript/DATA_REQUESTS.md` from git history. On arrival: write `metadata/neptune_clinical.csv`, append to `wgcna_traits.csv`, rerun `08b`, update the trait-count assertion in `R/99_verify_phase_d.R`.

3. **P1 — Fix the DOCX and finish the submission package**
   - Owner: user for metadata; agent for document updates.
   - Resume: `Manuscript/Paper1_manuscript.docx`. Fix §3.4 items; confirm EBioMedicine caps; add ORCIDs/CRediT/reviewers/declarations; re-run `officecli validate`, structure issues, LibreOffice page count (expect 19 pages, no blank page), and the 13-legend/44-reference/65-panel checks.

4. **P1 — Decide Track B scope**
   - Owner: user. Trigger: before Track B manuscript finalisation.
   - Resume: `metadata/trackB_panel_register.csv`; if another retrieval is funded, resume from `metadata/trackB_retrieval_log.md`, then rerun `10b → 10c → 10d → 99_verify_trackb.R`.

5. **P1 — Align the Zenodo/repository metadata**
   - Owner: agent on request.
   - Zenodo's version field currently reads `MN_transcriptomic_signature`; change to `v1.0.0` in the Zenodo record if desired. `README.md` and `CITATION.cff` still say the DOI is "to be issued"; write `10.5281/zenodo.24358387` into both and push.

6. **P1 — Commit the user's Data sharing edit**
   - Owner: agent on request.
   - `Manuscript/Paper1_manuscript.docx` is modified but uncommitted (user pasted a shortened statement with typos). Fix the typos, then commit and push so the repository copy matches the manuscript.

7. **P2 — Reopen Phase F only with new clinical data**
   - Owner: user supplies data; agent reruns on request.
   - Resume: current negative results in `output/tables/subtype_*` and the negative-subtype legend in the manuscript.

8. **P2 — Housekeeping**
   - Seed `scDblFinder` before any `07b` rerun; add parallelism waivers without changing settled computation; install adoption-test packages only if requested.

## 5. Pitfalls & Anti-patterns (Never repeat these)

### 5.1 Claim limits that must not be overstated

- The gate beats the panel in ERCB; never claim panel superiority.
- The strict refit-permutation null is not significant (*p* = 0.11). Report both nulls.
- Array-to-array transfer is the hard direction; never restore cohort-internal z-scoring to recover old AUCs.
- Tubulointerstitial discrimination is chance (0.449/0.539); there is no TI MN module.
- FSGS is the hard comparator (0.632); the multiclass balanced accuracy is 0.44, so this is MN detection, not a diagnostic multiclassifier.
- MYOZ2 has tissue-level reproducibility with unresolved cellular source; not proven podocyte or SMC.
- SLAMF8 detection is only 9.5% in RNA-seq and 27% in LCM; report the three-gene sensitivity model.
- Deconvolution is fragile; filtered signature primary, all-gene sensitivity; exclude pDC/neutrophil.
- Subtypes are a negative result; the array split tracks technical quality.
- Bulk↔pseudobulk concordance is weak (best Pearson 0.29); pseudobulk DE is hypothesis-grade.
- GSE197307/GSE182380 have no participant linkage; never infer overlap from titles.
- Within-MN PLA2R (6 versus 6) is uninterpretable.
- Expanded Track B is post hoc; do not present median R5 0.674 as preregistered confirmation.
- Do not claim a systematic PubMed/GEO search or an unverified EBioMedicine cap.

### 5.2 New technical lessons from this session

- **A DOCX saved by a word processor can fail schema validation** even though it opens and renders. After any manual edit, run `officecli validate` and `officecli view … issues --type structure`; the 80 `uiPriority` errors came from exactly this.
- **Never edit a paragraph by deleting its runs and then reading `p.text`.** The text is already empty at that point; this silently blanked the Abstract Methods paragraph twice. Read the text into a variable first, then rebuild runs.
- **Collapsing runs destroys `Verbatim Char` styling** (GSE identifiers, package calls rendered in red monospace). Rebuild paragraphs with a token-aware splitter and verify the styled-run count is unchanged (49 → 49).
- **The DOCX is now the only manuscript source.** Do not look for a Markdown file; the Pandoc rebuild path no longer exists.
- **Reference TSV row order ≠ manuscript order.** Map by DOI, never by row position.
- **PubMed DOI search can miss valid articles** (JSS, Biometrics JSTOR forms). Verify `glmnet` and `DeLong` by PMID/title before concluding a PMID does not exist.
- **Network reality:** `eutils.ncbi.nlm.nih.gov`, `api.openalex.org`, `doi.org` and `zenodo.org` work; `api.crossref.org`, Semantic Scholar and direct `github.com:22` do not. GitHub SSH must go through `ssh.github.com:443`.
- **A deploy key is repo-scoped and machine-bound.** It needs "Allow write access" to push; the public half goes on the web page, the private half never leaves the pushing machine. A placeholder key with `...` is invalid.
- **Do not invent identifiers.** PMIDs, DOIs and reference metadata must come from a resolver or the committed TSV; never hand-type.
- **GitHub's repository `size` field can read 0 immediately after a push** while it recomputes; verify by `git ls-remote` and the remote tree, not by the API counter.

### 5.3 Inherited project pitfalls (condensed; full P-1…P-74 list is in git history at `c9977e6` and earlier)

- Compute patient overlap across all releases before merging any matrix; never infer overlap from titles.
- Use `oligo` + Brainarray CDFs for HuGene-2.1-ST; never `affy`.
- Never `match()` each gene against the 1.02M-row probe index; build one lookup.
- Chunk large ST cohorts (`CHUNK = 40`); do not process 244 arrays at once.
- Align platform gene unions with NA fill; never naively `rbind()`.
- Never redirect stdout to a file a script writes with `writeLines()`.
- Data acquisition is manual: never download datasets; present a request list.
- `metafor::rma(method = "REML")` may not converge; wrap and count DL fallbacks.
- Use `fgsea::plotEnrichment`, not `enrichplot::gseaplot2()`.
- Scope verifier figure assertions to the generating step's filename prefix.
- Foreground shells cap at two minutes; run long jobs in the background and inspect logs.
- ERCB codes MN as `MGN`; `pROC::roc()` needs a two-level response.
- PAC ≤ 0.1 is the subtype stability gate, not silhouette.
- Stable splits can be technical; test correlation, mean expression and detection rate.
- `file.exists()` accepts zero-byte files; check `file.size() > 0`.
- Call `JoinLayers()` immediately after Seurat 5 `merge()`.
- Prefix barcodes before merging; rebuild the same key downstream.
- NNLS needs joint quantile normalisation, unit sums, τ > 0.5 filtering, τ on CPM.
- exfat has no hardlinks/symlinks and is case-insensitive; never `ln -f`.
- WGCNA 1.74 `modulePreservation` needs numeric reference indices, two-level `multiData`, `multiColor`, `.qual` and unambiguous module-size mapping.
- Split-half tests must reselect soft power in each half; report Jaccard separately from eigengene bicor.
- Require ≥95% ERCB coverage at panel selection and assert it.
- Summarise each permutation null from its own distribution, including `null_max`.
- Avoid `merge()` collisions by naming derived columns distinctly.
- GEO characteristic lines are per-field vectors; match line length to sample count.
- Ghostscript cannot read `/run/media`; stage through `/tmp`, embed all fonts, verify with `pdffonts`.
- Publisher PDFs can be Cloudflare-403; use live open-access JATS XML.
- Never hand-type citations; `Manuscript/references.tsv` is the only source.
- Assert on assembled prose, not just on successful script completion.
- Draw each nested-CV partition once per repeat and prove test folds are disjoint.
- Fixing transductive standardisation can lower external AUC; never revert to cohort-internal scaling.
- `ggsave(device = "jpeg")` fails on R 4.6.1; use `device = grDevices::jpeg`. `tools/make_figures.py` is unaffected.
- Never infer a panel's content from its filename or an older legend.
- Content-driven page sizing can create a second page; assert one page and delete `<plate>-N.jpg` leftovers.
- Re-derive counts from the filesystem whenever artefacts change and state "as of <date>".
- Never type a number into a script when an artefact supplies it; read it and assert reproduction.
- Main figures must defend claims, not merely show routine QC; renumber downstream references in descending order when inserting a plate.

## 6. Key Resources & Environment Info

### 6.1 Repository and paths

- Project root: `/run/media/meng/MY_Data/MN/Paper_1`; data root: `/run/media/meng/MY_Data/MN/Datasets`.
- `R/_lib.R` resolves the root and data paths; overrides `PAPER1_PROJ` and `PAPER1_DATA`.
- The mount can alternate `/run/media/{meng,yang}/MY_Data`; keep scripts path-independent.
- `work/`, `.Rlib/`, `scrna/extracted/`, `logs/` and five bulk CSVs >3 MB are gitignored as regenerable.
- Public repository: https://github.com/yangmylab/MN_transcriptomic_signature (MIT, `main`).
- Git remote: `origin` = `git@github-mn:yangmylab/MN_transcriptomic_signature.git`; SSH host `github-mn` → `ssh.github.com:443`, key `~/.ssh/mn_github_deploy_yang`.
- Zenodo: version DOI `10.5281/zenodo.24358387`; concept DOI `10.5281/zenodo.24358386`.

### 6.2 Languages and packages

- **R 4.6.1**, Bioconductor 3.23.
- **Python 3.14.4**; figure composition uses LaTeX/TikZ, not matplotlib.
- **Pandoc 3.8.3**, **OfficeCLI 1.0.144**, TeX Live, Ghostscript, LibreOffice for production/QA.
- Two R libraries:
  - `<PROJ>/.Rlib`: makecdfenv, five Brainarray CDFs, oligo, WGCNA 1.74, ConsensusClusterPlus, caret, glmnet, randomForest, xgboost, Boruta, metafor 5.2.1, msigdbr, rms, babelgene.
  - `/home/meng/R/x86_64-pc-linux-gnu-library/4.6`: ~1,701 packages, including Seurat 5.6.0, scDblFinder 1.26.9, harmony 2.0.5, muscat, nnls, limma, edgeR, DESeq2, sva, clusterProfiler, fgsea, ComplexHeatmap, survival.
- Installed optional: BayesPrism, MuSiC, hdWGCNA. Absent: ReactomePA, immunedeconv/xCell/quanTIseq, msigdbdf, liana, decoupleR, InstaPrism, tabpfn.

### 6.3 Reproducibility record and restoration

- Lockfile: `environment/renv.lock` (1,669 packages); notes in `environment/README.md`, `environment/sessionInfo.txt`, `environment/packages_*.csv`.
- Restore:

```r
install.packages("renv")
renv::restore(lockfile = "environment/renv.lock")
```

- Rebuild local Brainarray CDFs (no network):

```bash
Rscript R/00d_install_cdf.R
Rscript R/02a_build_cdf.R
for p in work/cdf/*; do R CMD INSTALL -l .Rlib "$p"; done
export PAPER1_DATA=/path/to/MN/Datasets
```

- Cross-machine: include `R/`, `metadata/`, `output/`, `reference/`, `Manuscript/`, `Docs/` and the `.rds` files under `work/`, plus `MN/Datasets`. Install R 4.6.x + Bioconductor 3.23, restore packages, rebuild `.Rlib`, regenerate `scrna/extracted/` with `Rscript R/07a_scrna_manifest.R`, install TeX Live + Ghostscript, then run all seven verifiers.

### 6.4 Compute resources and key paths

- Linux amd64, 20 threads, 61 GB RAM, ~1.8 TB free; local workstation only, no HPC jobs.
- Project ~32 GB: `work/` 23 GB, `.Rlib/` 4.2 GB, `scrna/` 3.3 GB, `output/` ~650 MB; shared `Datasets` 18 GB; exfat.
- `work/cel/` (16 GB) and `work/rma_cache/` (1.1 GB) may be deleted once `work/rma/` is final.

| Path | Contents |
|---|---|
| `$DATA` | Shared GEO collection: `<GSE>_RAW.tar`, `<GSE>/{raw,processed}/`, `manifest.tsv` |
| `work/sets/` | Seven cohort sets: `*_expr.rds`, `*_expr_symbol.rds`, `*_samples.csv` |
| `work/rma/`, `work/rnaseq/`, `work/combat/`, `work/msigdb/` | Normalised cohorts, ComBat TI matrices, MSigDB v2025.1.Hs GMTs |
| `work/scrna/` | QC objects, 236k-cell atlas, donor × gene pseudobulk matrices |
| `work/wgcna/` | Networks, eigengenes, traits, preservation objects; run `08d` from project root |
| `work/bench/10b_harmonised.rds` | Track B primary and comparison matrices |
| `output/tables/` | 240 result CSVs (121 at top level) |
| `output/figures/`, `output/figures/supplementary/` | 13 plates, PDF + 300-dpi JPEG; LaTeX in `work/figtex/` |
| `Manuscript/` | DOCX, `references.tsv`, `Figures/` (26 files) |
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
| `tools/make_cohort_design.py` | Data-derived main-Fig-1 panels |
| `tools/make_figures.py` | Plate compositor; calls the cohort-design builder; numbered titles |
| `tools/moyz2_provenance.R` | Separate MYOZ2 provenance test (~4 min cold) |
| `tools/make_environment_record.R`, `tools/rerun_*.sh` | Environment record and resume drivers |

`*` means run in background (>2 min), especially `04c`, `04d`, `07b`, `07c`, `07d`, `08c`, `08d`, `08f`, `05c`, `10d`.

### 6.6 Rebuild and verification commands

- Panel revisit: `09a → 09b → 09c → 09d → 09e → 09f`, then `python3 tools/make_figures.py`.
- WGCNA revisit: `08a → 08b → 08c → 08d → 08e` (`08f --force`), or `bash tools/rerun_from_08d.sh`.
- Track B: `10b → 10c → 10d → 99_verify_trackb.R`.
- Figure rebuild: `python3 tools/make_figures.py`; writes PDF + JPEG to `output/figures/` and `supplementary/`, one page per plate, fonts embedded through `/tmp`. Then resync `Manuscript/Figures/`:

```bash
rm -f Manuscript/Figures/*.pdf Manuscript/Figures/*.jpeg
cp output/figures/Fig*.pdf output/figures/Fig*.jpeg \
   output/figures/supplementary/FigS*.pdf output/figures/supplementary/FigS*.jpeg \
   Manuscript/Figures/
```

- Manuscript: the DOCX is edited directly (no Markdown source). Validate with `officecli validate`, `officecli view … issues --type structure`, and a LibreOffice render (expect 19 pages, no blank/orphan).
- Repository sync:

```bash
git add -A && git commit -m "<message>" && git push origin main
```

### 6.7 Network constraints (verified 2026-11-19)

- Reachable: `eutils.ncbi.nlm.nih.gov`, `api.openalex.org`, `doi.org`, `zenodo.org`, `mirrors.aliyun.com`.
- Not reachable: `api.crossref.org`, `api.semanticscholar.org`, direct `github.com:22` (use `ssh.github.com:443`), `www.google.com`.
- Consequence: PMID/DOI resolution works through PubMed E-utilities and OpenAlex; Crossref-only checks do not.

## 7. Startup Instruction for the Next Session

> **When starting a new session, the very first message must be: "Read handoff.md first."**

1. Restate the current state from this document before making changes. If a newer project file postdates this checkpoint, first determine whether it represents unrecorded work.
2. Confirm the frozen state with all seven verifiers, expecting **54 / 86 / 32 / 42 / 117 / 106 / 37**, all exit 0. Check for zero-byte outputs and re-derive any count before quoting it.
3. Do not quote pre-audit results, do not add methods to Paper 1, do not hard-code values already stored in result tables, and do not treat `gsm` as a GEO accession for RNA-seq cohorts.
4. The manuscript is **DOCX-only**; there is no Markdown source. Edit it with python-docx, preserve `Verbatim Char` styling, and re-run `officecli validate` afterwards.
5. B15 is resolved at 12 panels. The live decisions are whether Track B ships at 12, and the P1 submission metadata (ORCIDs, CRediT, funding wording, reviewers, declarations, EBioMedicine caps). P0 IF work and P0b clinical-data requests are externally blocked.
6. Commit `handoff.md` with any completed manuscript or figure changes before closing the session. Do not commit generated work without reviewing the diff first.
