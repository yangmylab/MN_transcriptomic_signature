#!/usr/bin/env Rscript
# Paper 1 - Phase D, step 08c: is the module structure reproducible at all?
#
# 08a showed that dropping the 3 QC-flagged glomerular samples (2% of the cohort) reshuffled the
# gene partition (ARI ~ 0). That is either a real fragility of the glomerular co-expression
# structure or a knock-on effect of those particular samples, and the difference matters for
# Phase E: hub genes taken from a fragile module cannot anchor a biomarker panel.
#
# This script separates the two explanations with four perturbations per compartment, each
# refitted at the SAME soft power as the full network:
#   split1/2/3  : stratified halves (2 x n/2), stratifying by diagnosis so the split cannot
#                 confound biology; the halves are scored against each other. Each half re-selects
#                 its own soft power (a full-cohort power over-sparsifies a half-size network)
#   drop_flagged: the full cohort minus the QC-flagged samples (the 08a sensitivity fit)
#
# Every module is scored two ways, because gene-membership Jaccard is harsh when a module simply
# splits, while eigengene correlation asks the question that matters (does this axis of
# co-expression still exist?):
#   best Jaccard  : |A n B| / |A u B| against the best-matching module of the other fit
#   best ME cor   : bicor of the two module eigengenes, computed on ALL samples of the full cohort
#
# Outputs
#   output/tables/wgcna_stability_<comp>.csv          one row per perturbation
#   output/tables/wgcna_stability_modules_<comp>.csv  one row per full-fit module per perturbation
#   output/tables/wgcna_stability_summary.csv         both compartments, verdict line
#   output/wgcna/08c_<comp>_stability.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(WGCNA) })

tab_dir <- file.path(OUT, "tables")
wk_dir  <- file.path(PROJ, "work", "wgcna")
fig_dir <- file.path(OUT, "wgcna")

n_threads <- max(1L, parallel::detectCores() - 1L)
allowWGCNAThreads(nThreads = n_threads)

SEED <- 20261019L
MIN_MODULE <- 30L
MERGE_CUT  <- 0.25
DEEP_SPLIT <- 2L
MAX_POUT   <- 0.05
N_SPLITS   <- 3L
JACC_STRONG <- 0.5

fig <- function(stem, width, height, draw) {
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw()
  dev.off()
  jpeg(paste0(stem, ".jpeg"), width = width, height = height, units = "in", res = 300,
       quality = 95, bg = "white")
  draw()
  dev.off()
}

write_tab <- function(x, name) {
  f <- file.path(tab_dir, name)
  write.csv(x, f, row.names = FALSE)
  cat("  wrote", f, "(", nrow(x), "rows )\n")
}

refit <- function(dat, power, seed) {
  set.seed(seed)
  blockwiseModules(dat, power = power, networkType = "signed", TOMType = "signed",
                   corType = "bicor", maxPOutliers = MAX_POUT,
                   maxBlockSize = ncol(dat), minModuleSize = MIN_MODULE,
                   deepSplit = DEEP_SPLIT, mergeCutHeight = MERGE_CUT,
                   numericLabels = TRUE, pamRespectsDendro = FALSE,
                   saveTOMs = FALSE, randomSeed = seed, nThreads = n_threads, verbose = 0)
}

# per-module best match between two colourings of the SAME gene universe
compare_modules <- function(colA, colB, dat) {
  modsA <- setdiff(unique(colA), "grey")
  modsB <- setdiff(unique(colB), "grey")
  if (!length(modsA) || !length(modsB)) return(NULL)
  me_all <- function(cols) moduleEigengenes(dat, colors = cols, excludeGrey = TRUE)$eigengenes
  meA <- me_all(colA)
  meB <- me_all(colB)
  cm  <- suppressWarnings(bicor(meA, meB, use = "p"))
  rows <- lapply(modsA, function(m) {
    A <- names(colA)[colA == m]
    js <- vapply(modsB, function(n) {
      B <- names(colB)[colB == n]
      length(intersect(A, B)) / length(union(A, B))
    }, numeric(1))
    best <- modsB[which.max(js)]
    me_cor <- abs(cm[paste0("ME", m), paste0("ME", best)])
    data.frame(module = m, size = length(A), best_match = best,
               jaccard = round(max(js), 4), me_cor = round(me_cor, 4),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

strat_split <- function(group, seed, take = 0.5) {
  set.seed(seed)
  idx <- unlist(lapply(split(seq_along(group), group), function(i) {
    sample(i, max(2, floor(length(i) * take)))
  }), use.names = FALSE)
  sort(idx)
}

qc_flagged <- read.csv(file.path(tab_dir, "qc_flagged_samples_retained.csv"), stringsAsFactors = FALSE)

# the cohort a compartment's discovery matrix came from (the net rds stores traits without `set`)
SET_OF <- c(glom = "disc_glom_array", ti = "disc_ti_array")

R2_CUT <- 0.85
R2_FALLBACK <- 0.80
POWERS <- c(1:10, 12, 14, 16, 18, 20)

# A subsample is HALF the cohort, so the full-cohort soft power over-sparsifies it (the same
# threshold cuts weaker correlations). Each subsample therefore re-selects its own power with the
# identical rule as 08a - otherwise the stability test would penalise WGCNA for a parameter that no
# one would actually keep fixed.
local_power <- function(dat) {
  sft <- pickSoftThreshold(dat, powerVector = POWERS, networkType = "signed",
                           RsquaredCut = R2_CUT, verbose = 0)
  fi <- sft$fitIndices
  r2 <- fi$SFT.R.sq
  ok <- !is.na(r2) & r2 >= R2_CUT
  if (any(ok)) return(fi$Power[which(ok)[1]])
  ok2 <- !is.na(r2) & r2 >= R2_FALLBACK
  if (any(ok2)) return(fi$Power[which(ok2)[1]])
  fi$Power[which.max(r2)]
}

summary_rows <- list()
for (comp in c("glom", "ti")) {
  cat("\n================ ", toupper(comp), " ================\n", sep = "")

  full <- readRDS(file.path(wk_dir, paste0("net_", comp, ".rds")))
  dat  <- full$dat
  tr   <- full$traits
  pw   <- full$power
  stopifnot(identical(rownames(dat), tr$gsm))

  full_colors <- setNames(full$module_table$module[match(colnames(dat), full$module_table$gene)],
                          colnames(dat))
  n_full_mod <- length(setdiff(unique(full_colors), "grey"))
  cat(sprintf("full fit: %d genes, %d samples, power %d, %d modules (excl. grey)\n",
              ncol(dat), nrow(dat), pw, n_full_mod))

  per_mod <- list()
  pert_rows <- list()

  for (k in seq_len(N_SPLITS)) {
    s1 <- strat_split(tr$group, SEED + k)
    s2 <- setdiff(seq_len(nrow(dat)), s1)
    cat(sprintf("split %d: %d / %d samples\n", k, length(s1), length(s2)))
    pw1 <- local_power(dat[s1, ])
    pw2 <- local_power(dat[s2, ])
    cat(sprintf("  subsample soft powers: %d / %d (full-cohort power %d)\n", pw1, pw2, pw))
    c1 <- labels2colors(refit(dat[s1, ], pw1, SEED + k)$colors)
    c2 <- labels2colors(refit(dat[s2, ], pw2, SEED + 100L + k)$colors)
    # colours index GENES (columns of dat), not samples - naming them with rownames() silently
    # compared sample ids and reported meaningless Jaccard values
    g1 <- setNames(c1, colnames(dat))
    g2 <- setNames(c2, colnames(dat))
    cmp <- compare_modules(g1, g2, dat)
    if (is.null(cmp)) next
    per_mod[[paste0("split", k)]] <- cbind(perturbation = paste0("split", k), cmp)
    pert_rows[[paste0("split", k)]] <- data.frame(
      compartment = comp, perturbation = paste0("split", k),
      n_samples_ref = length(s1), n_samples_alt = length(s2),
      power_ref = pw1, power_alt = pw2, power_full = pw,
      n_modules_ref = length(setdiff(unique(g1), "grey")),
      n_modules_alt = length(cmp$best_match),
      median_best_jaccard = round(median(cmp$jaccard), 4),
      max_best_jaccard = round(max(cmp$jaccard), 4),
      median_best_ME_cor = round(median(cmp$me_cor), 4),
      max_best_ME_cor = round(max(cmp$me_cor), 4),
      n_modules_jaccard_gt_0.5 = sum(cmp$jaccard > JACC_STRONG),
      stringsAsFactors = FALSE)
    cat(sprintf("  median best Jaccard %.3f | median best ME cor %.3f\n",
                median(cmp$jaccard), median(cmp$me_cor)))
  }

  # flagged-sample removal, scored against the full fit
  keep <- !(rownames(dat) %in% qc_flagged$gsm[qc_flagged$set == SET_OF[[comp]]])
  if (sum(!keep) > 0) {
    cat(sprintf("drop_flagged: %d flagged samples removed\n", sum(!keep)))
    pw_f <- local_power(dat[keep, ])
    cf <- labels2colors(refit(dat[keep, ], pw_f, SEED + 7L)$colors)
    gf <- setNames(cf, colnames(dat))
    cmp <- compare_modules(full_colors, gf, dat)
    if (!is.null(cmp)) {
      per_mod[["drop_flagged"]] <- cbind(perturbation = "drop_flagged", cmp)
      pert_rows[["drop_flagged"]] <- data.frame(
        compartment = comp, perturbation = "drop_flagged",
        n_samples_ref = nrow(dat), n_samples_alt = sum(keep),
        power_ref = pw, power_alt = pw_f, power_full = pw,
        n_modules_ref = n_full_mod,
        n_modules_alt = length(setdiff(unique(gf), "grey")),
        median_best_jaccard = round(median(cmp$jaccard), 4),
        max_best_jaccard = round(max(cmp$jaccard), 4),
        median_best_ME_cor = round(median(cmp$me_cor), 4),
        max_best_ME_cor = round(max(cmp$me_cor), 4),
        n_modules_jaccard_gt_0.5 = sum(cmp$jaccard > JACC_STRONG),
        stringsAsFactors = FALSE)
      cat(sprintf("  median best Jaccard %.3f | median best ME cor %.3f\n",
                  median(cmp$jaccard), median(cmp$me_cor)))
    }
  }

  pm  <- do.call(rbind, per_mod)
  pr  <- do.call(rbind, pert_rows)
  rownames(pm) <- NULL
  write_tab(pm, paste0("wgcna_stability_modules_", comp, ".csv"))
  write_tab(pr, paste0("wgcna_stability_", comp, ".csv"))

  # veredict per compartment: reproducible only if split halves keep membership AND the axis
  verdict <- if (median(pr$median_best_jaccard[grepl("^split", pr$perturbation)]) >= 0.3 &&
                 median(pr$median_best_ME_cor[grepl("^split", pr$perturbation)]) >= 0.6) {
    "reproducible"
  } else if (median(pr$median_best_ME_cor[grepl("^split", pr$perturbation)]) >= 0.6) {
    "axes reproducible, membership not"
  } else {
    "not reproducible"
  }
  cat("verdict:", verdict, "\n")

  summary_rows[[comp]] <- data.frame(
    compartment = comp, soft_power = pw, n_modules_full = n_full_mod,
    median_split_jaccard = round(median(pr$median_best_jaccard[grepl("^split", pr$perturbation)]), 4),
    median_split_ME_cor = round(median(pr$median_best_ME_cor[grepl("^split", pr$perturbation)]), 4),
    drop_flagged_jaccard = pr$median_best_jaccard[pr$perturbation == "drop_flagged"],
    drop_flagged_ME_cor  = pr$median_best_ME_cor[pr$perturbation == "drop_flagged"],
    verdict = verdict, stringsAsFactors = FALSE)

  fig(file.path(fig_dir, paste0("08c_", comp, "_stability")), 8, 4.6, function() {
    par(mfrow = c(1, 2), mar = c(7, 4.5, 3, 1))
    bp1 <- barplot(pr$median_best_jaccard, names.arg = pr$perturbation, las = 2,
                   col = "steelblue", border = NA, ylim = c(0, 1),
                   ylab = "median best Jaccard",
                   main = paste0("Module membership - ", comp))
    abline(h = 0.3, lty = 2, col = "grey40")
    abline(h = JACC_STRONG, lty = 3, col = "grey60")
    bp2 <- barplot(pr$median_best_ME_cor, names.arg = pr$perturbation, las = 2,
                   col = "firebrick", border = NA, ylim = c(0, 1),
                   ylab = "median best eigengene bicor",
                   main = paste0("Co-expression axes - ", comp))
    abline(h = 0.6, lty = 2, col = "grey40")
  })

  rm(full, dat)
  gc()
}

summ <- do.call(rbind, summary_rows)
write_tab(summ, "wgcna_stability_summary.csv")
cat("\n--- stability summary ---\n")
print(summ[, c("compartment", "median_split_jaccard", "median_split_ME_cor", "verdict")])
cat("\n08c done\n")
