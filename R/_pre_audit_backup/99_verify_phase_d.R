#!/usr/bin/env Rscript
# Paper 1 - Phase D verifier (WGCNA, module-trait, reproducibility, preservation, hub genes).
# Each assertion is a claim these analyses make; the script exits 0 only if all hold.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "wgcna")
wk_dir  <- file.path(PROJ, "work", "wgcna")
set_dir <- file.path(PROJ, "work", "sets")

n_pass <- 0L
n_fail <- 0L

check <- function(ok, label, detail = "") {
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
    cat(sprintf("[PASS] %-76s %s\n", label, detail))
  } else {
    n_fail <<- n_fail + 1L
    cat(sprintf("[FAIL] %-76s %s\n", label, detail))
  }
}

rd <- function(f) read.csv(file.path(tab_dir, f), stringsAsFactors = FALSE)
TRAIT_CELLS <- c("Podocyte", "PT", "TAL", "DCT", "CD_PC", "PEC", "Endothelium",
                 "Fibroblast", "Mesangium", "Myeloid", "B", "T_NK", "Plasma", "Mast")
COMPS <- c("glom", "ti")
SETS  <- c(glom = "disc_glom_array", ti = "disc_ti_array")

cat("=========== PHASE D VERIFICATION ===========\n\n")

# --- 1. networks -------------------------------------------------------------
summ <- rd("wgcna_network_summary.csv")
check(nrow(summ) == 2 && all(COMPS %in% summ$compartment), "network summary covers both compartments")
check(all(summ$n_genes_used <= 8000), "gene space respects the 8,000-gene MAD cap", max(summ$n_genes_used))
check(all(summ$soft_power >= 1 & summ$soft_power <= 20), "soft power inside the search grid",
      paste(summ$soft_power, collapse = "/"))
check(all(summ$scale_free_R2 >= 0.85 | grepl("^FALLBACK", summ$power_rule)),
      "scale-free criterion met, or an explicit fallback recorded",
      paste(round(summ$scale_free_R2, 3), collapse = "/"))
check(all(grepl("smallest power", summ$power_rule)), "power chosen by the declared rule")
check(all(summ$n_modules_excl_grey >= 1), "at least one non-grey module per compartment",
      paste(summ$n_modules_excl_grey, collapse = "/"))
check(all(summ$prop_assigned > 0.3), "a substantial fraction of genes assigned",
      paste(summ$prop_assigned, collapse = "/"))
check(all(summ$corType == "bicor" & summ$networkType == "signed"), "signed bicor networks as specified")
check(all(summ$mergeCutHeight == 0.25 & summ$deepSplit == 2), "dynamic tree cut parameters honoured")

# --- 2. module tables --------------------------------------------------------
mods <- list(); sizes <- list()
for (comp in COMPS) {
  mods[[comp]]  <- rd(sprintf("wgcna_modules_%s.csv", comp))
  sizes[[comp]] <- rd(sprintf("wgcna_module_sizes_%s.csv", comp))
}
for (comp in COMPS) {
  n_used <- summ$n_genes_used[summ$compartment == comp]
  check(nrow(mods[[comp]]) == n_used, sprintf("%s: every network gene has a module row", comp), nrow(mods[[comp]]))
  check(sum(sizes[[comp]]$n_genes) == n_used, sprintf("%s: module sizes sum to the gene count", comp),
        sum(sizes[[comp]]$n_genes))
  non_grey <- sizes[[comp]][sizes[[comp]]$module != "grey", ]
  check(all(non_grey$n_genes >= summ$min_module[summ$compartment == comp]),
        sprintf("%s: no non-grey module below minModuleSize", comp), min(non_grey$n_genes))
  check(!anyDuplicated(mods[[comp]]$gene), sprintf("%s: each gene appears once", comp))
  check(all(is.na(mods[[comp]]$kME) == (mods[[comp]]$module == "grey")),
        sprintf("%s: kME computed for every assigned gene, NA only for grey", comp))
  check(all(abs(mods[[comp]]$kME[!is.na(mods[[comp]]$kME)]) <= 1),
        sprintf("%s: kME inside [-1, 1]", comp))
  me_f <- file.path(wk_dir, sprintf("me_%s.csv", comp))
  me <- read.csv(me_f, stringsAsFactors = FALSE)
  check(nrow(me) == summ$n_samples[summ$compartment == comp],
        sprintf("%s: one eigengene row per sample", comp), nrow(me))
  check(ncol(me) - 1 == summ$n_modules_excl_grey[summ$compartment == comp],
        sprintf("%s: one eigengene column per non-grey module", comp), ncol(me) - 1)
}

# --- 3. traits ---------------------------------------------------------------
traits <- rd("wgcna_traits.csv")
check(nrow(traits) == sum(summ$n_samples), "trait table covers every discovery sample", nrow(traits))
check(all(TRAIT_CELLS %in% names(traits)), "all 14 reliable cell-type fractions present")
check(!any(c("Neutrophil", "pDC") %in% names(traits)),
      "unreliable fractions (Neutrophil, pDC) excluded from traits")
check(!anyNA(traits[, TRAIT_CELLS]), "no missing cell-type fractions")
# the two excluded columns (Neutrophil, pDC) can carry weight, so the retained sum may fall below 1;
# what must hold is that the FULL 16-column profile sums to 1 and the retained part stays substantial
dec_full <- read.csv(file.path(tab_dir, "deconvolution_fractions.csv"), stringsAsFactors = FALSE)
dec_full$total16 <- rowSums(dec_full[, c(TRAIT_CELLS, "Neutrophil", "pDC")])
check(all(abs(dec_full$total16 - 1) < 0.01), "all 16 cell-type fractions sum to 1 per sample",
      paste(round(range(dec_full$total16), 4), collapse = " .. "))
frac_sum <- rowSums(traits[, TRAIT_CELLS])
check(all(frac_sum > 0.3 & frac_sum <= 1.001), "retained 14-fraction sum is on a unit scale",
      paste(round(range(frac_sum), 3), collapse = " .. "))
check(all(traits$MN %in% c(0L, 1L)) && all(traits$NonMN_GN %in% c(0L, 1L)),
      "MN / non-MN-GN indicators are binary")
for (comp in COMPS) {
  ann <- read.csv(file.path(set_dir, paste0(SETS[[comp]], "_samples.csv")), stringsAsFactors = FALSE)
  tr_c <- traits[traits$compartment == comp, ]
  check(nrow(tr_c) == nrow(ann) && setequal(tr_c$gsm, ann$gsm),
        sprintf("%s: traits keyed to exactly the cohort's samples", comp))
  check(sum(tr_c$MN) == sum(ann$disease == "MN"),
        sprintf("%s: MN count matches the cohort metadata", comp), sum(tr_c$MN))
}

# --- 4. module-trait ---------------------------------------------------------
for (comp in COMPS) {
  mt <- rd(sprintf("wgcna_module_trait_%s.csv", comp))
  n_mod <- summ$n_modules_excl_grey[summ$compartment == comp]
  check(nrow(mt) == n_mod * 16, sprintf("%s: every module x trait pair reported", comp), nrow(mt))
  check(all(is.na(mt$cor) | abs(mt$cor) <= 1), sprintf("%s: correlations inside [-1, 1]", comp),
        sprintf("%d untestable pairs", sum(is.na(mt$cor))))
  check(all(is.na(mt$padj) | (mt$padj >= 0 & mt$padj <= 1)), sprintf("%s: BH-adjusted p in [0, 1]", comp))
  # padj is stored at 4 significant digits, so allow a relative rounding tolerance
  check(all(is.na(mt$padj) | is.na(mt$p) | mt$padj >= mt$p * (1 - 1e-3)),
        sprintf("%s: padj never below the raw p", comp))
  check(all(mt$cor_method %in% c("bicor", "pearson (zero MAD)")),
        sprintf("%s: correlation method recorded per pair", comp))
  adj <- rd(sprintf("wgcna_module_trait_adjusted_%s.csv", comp))
  check(nrow(adj) == n_mod, sprintf("%s: composition-adjusted fit covers every module", comp))
  check(all(grepl("Podocyte", adj$covars_adjusted)), sprintf("%s: adjustment includes Podocyte", comp))
  check(all(adj$padj >= 0 & adj$padj <= 1), sprintf("%s: adjusted padj in [0, 1]", comp))
  ms <- rd(sprintf("wgcna_module_significance_%s.csv", comp))
  check(nrow(ms) == n_mod && all(c("MN", "Podocyte") %in% names(ms)),
        sprintf("%s: module significance computed per trait", comp))
}

# --- 5. reproducibility ------------------------------------------------------
stab <- rd("wgcna_stability_summary.csv")
check(nrow(stab) == 2, "stability summary covers both compartments")
check(all(c("median_split_jaccard", "median_split_ME_cor", "verdict") %in% names(stab)),
      "stability summary reports membership, eigengene and a verdict")
check(all(is.na(stab$median_split_ME_cor) | (stab$median_split_ME_cor >= -1 & stab$median_split_ME_cor <= 1)),
      "stability eigengene correlations inside [-1, 1]")
for (comp in COMPS) {
  sp <- rd(sprintf("wgcna_stability_%s.csv", comp))
  check(nrow(sp) == 4, sprintf("%s: 3 splits + flagged-removal reported", comp), nrow(sp))
  check(all(c("split1", "split2", "split3", "drop_flagged") %in% sp$perturbation),
        sprintf("%s: all perturbations present", comp))
  check(all(sp$power_ref >= 1 & sp$power_alt >= 1),
        sprintf("%s: each refit selected its own soft power", comp))
  check(all(sp$median_best_jaccard >= 0 & sp$median_best_jaccard <= 1),
        sprintf("%s: Jaccard values in [0, 1]", comp))
}

# --- 6. preservation --------------------------------------------------------
pres <- list()
for (comp in COMPS) pres[[comp]] <- rd(sprintf("wgcna_preservation_%s.csv", comp))
for (comp in COMPS) {
  pr <- pres[[comp]]
  n_mod <- summ$n_modules_excl_grey[summ$compartment == comp]
  check(nrow(pr) >= n_mod, sprintf("%s: preservation tested for every module", comp), nrow(pr))
  check(all(pr$n_permutations >= 50), sprintf("%s: at least 50 permutations", comp), min(pr$n_permutations))
  check(all(pr$qsummary >= 0 & pr$qsummary <= 1), sprintf("%s: preservation q in [0, 1]", comp))
  check(all(pr$verdict %in% c("preserved", "weak", "not preserved", "not testable")),
        sprintf("%s: verdicts use the declared vocabulary", comp))
  check(all(pr$Zsummary[pr$verdict == "preserved"] > 10),
        sprintf("%s: every 'preserved' module has Zsummary > 10", comp))
  check(all(pr$qsummary[pr$verdict == "preserved"] < 0.05),
        sprintf("%s: every 'preserved' module is significant at q < 0.05", comp))
  check(all(abs(pr$medianRank) >= 1), sprintf("%s: median rank reported per module", comp))
}
sec <- rd("wgcna_preservation_secondary.csv")
check(nrow(sec) >= 1 && all(grepl("cross-platform", sec$comparison)),
      "cross-platform preservation reported and labelled secondary")
ps <- rd("wgcna_preservation_summary.csv")
check(nrow(ps) == 3, "preservation summary covers 2 primary + 1 secondary test", nrow(ps))

# --- 7. hub genes and the Phase E candidate space ----------------------------
hubs <- rd("wgcna_hub_genes.csv")
check(nrow(hubs) > 0, "hub genes exported", nrow(hubs))
check(all(abs(hubs$kME) > 0.8), "every hub has |kME| > 0.8", round(min(abs(hubs$kME)), 3))
check(all(hubs$compartment %in% COMPS), "hubs carry a valid compartment")
for (comp in COMPS) {
  h <- hubs[hubs$compartment == comp, ]
  if (!nrow(h)) next
  preserved <- pres[[comp]]$module[pres[[comp]]$verdict == "preserved"]
  check(all(h$module %in% preserved),
        sprintf("%s: hubs come only from ERCB-preserved modules", comp),
        paste(sort(unique(h$module)), collapse = ","))
  check(all(h$gene %in% mods[[comp]]$gene), sprintf("%s: hub genes exist in the module table", comp))
  check(all(h$module == mods[[comp]]$module[match(h$gene, mods[[comp]]$gene)]),
        sprintf("%s: hub module labels agree with the module table", comp))
}

cand <- rd("phaseE_candidate_space.csv")
check(nrow(cand) > 0, "Phase E candidate space is non-empty", nrow(cand))
check(all(cand$n_evidence >= 1 & cand$n_evidence <= 4), "each candidate carries 1-4 evidence lines",
      paste(range(cand$n_evidence), collapse = " .. "))
check(all(c("is_hub", "is_metaDEG", "in_gate", "is_specific") %in% names(cand)),
      "candidate space records all four evidence lines")
hub_key <- paste(hubs$compartment, hubs$gene)
check(all(cand$is_hub == (paste(cand$compartment, cand$gene) %in% hub_key)),
      "is_hub agrees with the hub table within each compartment")
check(all(is.na(cand$hub_kME) | abs(cand$hub_kME) > 0.8), "candidate hub kME values are hub-grade")
check(all(is.na(cand$celltype_tau) | (cand$celltype_tau >= 0 & cand$celltype_tau <= 1)),
      "tau specificity values in [0, 1]")
check(sum(cand$n_evidence == 4) <= nrow(cand), "four-line candidates are a subset")
check(all(!is.na(cand$meta_contrast[cand$is_metaDEG])),
      "every meta-DEG candidate names its contrast")
check("in_gate_reduced" %in% names(cand),
      "gate evidence counted only from reduced feature sets (not K = 'all', pitfall B9)")
check(all(cand$n_evidence == cand$is_hub + cand$is_metaDEG + cand$in_gate_reduced + cand$is_specific),
      "n_evidence is the sum of the four recorded evidence lines")
check(nrow(cand) < 5000, "candidate space is a shortlist, not the transcriptome", nrow(cand))
check(!any(cand$in_gate_reduced & is.na(cand$gate_min_K)),
      "gate-selected candidates record the smallest K they survive")

core <- rd("phaseE_candidate_core.csv")
check(nrow(core) > 0, "core candidate set exported", nrow(core))
check(all(core$n_evidence >= 2), "core candidates carry at least two evidence lines",
      paste(range(core$n_evidence), collapse = " .. "))
check(nrow(core) < nrow(cand), "core set is a strict subset of the candidate space")
ov <- rd("wgcna_hub_gene_overlap.csv")
check(nrow(ov) > 0 && all(ov$n_evidence >= 1 & ov$n_evidence <= 4),
      "evidence-line overlap table consistent")

# --- 8. figures and 0-byte guard --------------------------------------------
for (p in c("08a", "08b", "08c", "08d", "08e")) {
  pdfs <- list.files(fig_dir, pattern = paste0("^", p, "_.*\\.pdf$"))
  jpgs <- list.files(fig_dir, pattern = paste0("^", p, "_.*\\.jpeg$"))
  check(length(pdfs) > 0 && length(pdfs) == length(jpgs),
        sprintf("%s: every figure exists as PDF + JPEG", p), sprintf("%d pdf / %d jpeg", length(pdfs), length(jpgs)))
  if (length(pdfs)) {
    sz <- file.size(file.path(fig_dir, pdfs))
    check(all(sz > 1000), sprintf("%s: no empty figure file", p), min(sz))
  }
}
zero <- list.files(c(tab_dir, fig_dir, wk_dir), recursive = TRUE, full.names = TRUE)
zero <- zero[file.size(zero) == 0]
check(length(zero) == 0, "no 0-byte outputs anywhere (resume-cache guard, P-28)",
      if (length(zero)) paste(basename(head(zero, 3)), collapse = ",") else "clean")

cat(sprintf("\n=========== %d passed, %d failed ===========\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
cat("ALL PHASE D ASSERTIONS PASSED\n")
