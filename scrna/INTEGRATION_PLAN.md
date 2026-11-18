# Paper 1 — scRNA-seq integration plan

**Question this arm answers:** *which cells carry the bulk MN signature, and can single-cell data make the bulk biomarker panel interpretable and defensible?*

**Division of labour with Paper 2.** Paper 2 builds the publishable MN single-cell/single-nucleus/spatial meta-atlas. Paper 1 does **not** repeat that atlas; it consumes it. Everything here is written so the outputs (signature matrix, cell-type attribution, pseudobulk concordance) can be handed to Paper 2 unchanged.

---

## 1. Why single-cell data belongs in Paper 1

| Paper 1 gap | What scRNA-seq adds | Where it lands |
|---|---|---|
| 584 glomerular meta-DEGs are anonymous tissue-level averages | cell-type attribution of every meta-DEG (podocyte / PEC / endothelial / mesangial / monocyte / T / B / plasma) | new Results subsection + Fig |
| Phase C deconvolution is blocked on CIBERSORTx + a licence token | a kidney-derived **signature matrix** built from the atlas, deconvolved with NNLS — no token, fully reproducible | unblocks Phase C |
| Phase E panel genes are bulk hits; reviewers will ask "expressed where?" | per-cell-type specificity, detection rate, and pseudobulk MN-vs-control per cell type | panel justification + supplement |
| MN vs FSGS is the hard contrast (gate AUC 0.80) | pseudobulk DE in matched cell types tells us whether the residual signal is podocyte-intrinsic or compositional | Discussion + rebuttal armoury |

## 2. Datasets and design

See `DATA_REQUEST.md`. Tier 1 = GSE271188 (snRNA, 37 donors, MN/IgAN/MCD/control/DKD), GSE171458 (6 IMN + 2 healthy), GSE241302 (3 IMN + 1 healthy).

**Comparator logic:** GSE271188 is the only cohort containing MN *and* IgAN/MCD/control measured in the same run, so it is the primary reference. The two smaller scRNA cohorts are replication for podocyte/PEC states, never for composition (too few donors).

## 3. Pipeline (scripts `R/07*`)

| Script | Job | Key outputs |
|---|---|---|
| `07a_scrna_manifest.R` | validate `scrna/raw/**`, extract archives to `scrna/extracted/`, emit one row per 10x sample with disease label | `scrna/samples.csv`, `scrna/extraction.log` |
| `07b_scrna_qc.R` | per-sample import (10x mtx/h5), MAD-based QC on nFeature/nCount/percent.mt, `scDblFinder` doublets, log-normalise, HVGs | `work/scrna/qc/<sample>.rds`, `output/tables/scrna_qc.csv`, QC figures |
| `07c_scrna_integrate_annotate.R` | Harmony on PCA (batch = donor), Leiden/Louvain sweep, canonical kidney markers + `SingleR`/celldex, immune subclustering | `work/scrna/atlas.rds`, `output/tables/scrna_celltype_counts.csv`, UMAP + marker figures |
| `07d_scrna_pseudobulk.R` | aggregate counts per donor × cell type; pseudobulk MN vs control (limma-voom); composition (propeller/miloR) | `output/tables/scrna_pseudobulk_<celltype>.csv`, composition tables |
| `07e_bulk_scrna_mapping.R` | (i) cell-type specificity of meta-DEGs and gate features; (ii) cell-type mean-expression **signature matrix**; (iii) concordance of bulk MN-vs-control logFC with pseudobulk logFC | `output/tables/signature_matrix.csv`, `bulk_deg_celltype_attribution.csv`, `bulk_vs_pseudobulk_concordance.csv` |
| `07f_deconvolution.R` | NNLS deconvolution of all 7 bulk cohort sets with the atlas matrix; fractions ↔ bulk signature score; MN vs comparator tests | `output/tables/deconvolution_fractions.csv`, `deconvolution_tests.csv`, figures |
| `99_verify_scrna.R` | assertion suite for the arm | exit 0 / 1 |

All figures follow the standing conventions: PDF **and** JPEG, no grid lines, filename prefixed with the generating step (`07b_…`, `07e_…`).

## 4. Memory strategy (14 GB box)

GSE271188 is ~262 k nuclei; a merged sparse matrix is ≈ 5 GB before any copy, so the naive route will thrash.

- Never hold the merged count matrix and a merged log-normalised copy at once.
- Integrate on the **HVG subset only** (2 000 genes × ~350 k cells is small).
- Build the full atlas with **Seurat v5 on-disk layers (BPCells)**; if BPCells is unavailable, pseudobulk per sample first and keep only the HVG matrix.
- `07d` runs before any full-atlas DE, so the statistical layer never needs the big object.
- Hard rule: `options(future.globals.maxSize = 4 * 1024^3)` and chunked `IntegrateLayers`.

### 3.1 Package status (installed 2026-10-18)

| Installed | Version | Used by |
|---|---|---|
| `scDblFinder` | 1.26.9 | 07b doublets |
| `miloR` | 2.8.1 | differential abundance (second method) |
| `muscat` | 1.26.0 | pseudobulk helpers / sensitivity |
| `decoupleR` | 2.17.0 | TF/pathway activity (optional) |
| `SoupX` | 1.6.2 | ambient RNA (optional) |
| `nnls`, `quadprog` | 1.6 / 1.5.8 | 07f deconvolution |

`liana` is **not published for Bioconductor 3.23** and is therefore unavailable. Cell–cell communication is outside the Paper-1 arm (it belongs to Paper 2), so this costs nothing here; if Paper 2 needs LIANA, install it from its GitHub repository at that point.

## 5. Statistical standards

- Single-cell DE **only** as pseudobulk (donor = replicate), `limma-voom` with donor as the unit; cell-level Wilcoxon is descriptive only.
- Composition: `propeller` (robust to small n) with `miloR` as a second method once ≥ 6 donors per group; report both, promote only concordant neighbourhoods.
- Deconvolution reported as relative proportions; absolute abundance only if total-cell content can be normalised (it cannot here) — state this.
- Every cell-type claim carries a detection-rate and a donor-count column.

## 6. Integration points with Phases C/D/E

| Phase | Before | After |
|---|---|---|
| C (deconvolution) | blocked on CIBERSORTx | `signature_matrix.csv` + `07f` NNLS; CIBERSORTx becomes optional cross-check |
| D (WGCNA) | modules ↔ MN status | modules ↔ cell-type fractions (module eigengene vs monocyte/podocyte proportion) |
| E (panel) | bulk-selected genes | panel filtered by cell-type specificity + pseudobulk support; off-target/ubiquitous genes dropped |

## 7. What would make this arm fail (stated up front)

1. **Only 6 MN donors in the reference cohort.** Pseudobulk DE for MN is therefore hypothesis-grade; it can corroborate the bulk signature but cannot establish it.
2. **snRNA vs scRNA mixture.** Nuclei under-represent cytoplasmic transcripts (complement factors, immunoglobulins). Cell-type attribution must use the snRNA cohort for nuclear genes and the scRNA cohorts for secreted/cytoplasmic genes; disagreements are reported, not averaged away.
3. **Ambient RNA in kidney biopsies.** Without SoupX/CellBender, podocyte genes leak into every droplet; specificity scores are computed only on cells passing a stringent detection filter, and the caveat is written into the methods.
4. **No participant linkage** between bulk and single-cell cohorts (see handoff B7). All bulk↔single-cell comparisons are at the **signature level**, never at the patient level.
